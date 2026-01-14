"""Config doctor utilities for reporting and appending missing TOML sections."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

try:  # Python 3.11+
    import tomllib  # type: ignore[attr-defined]
except ImportError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]

from typing import get_args, get_origin

from . import config_schema


@dataclass
class InvalidChoice:
    section: str
    field: str
    value: Any
    choices: List[Any]


@dataclass
class DoctorReport:
    missing_sections: List[str]
    missing_fields: Dict[str, List[str]]
    invalid_choices: List[InvalidChoice]


@dataclass
class SectionInfo:
    defaults: Dict[str, Any]
    field_meta: Dict[str, config_schema.FieldMeta]
    field_order: List[str]


@dataclass
class AppendPlan:
    blocks: List[str]
    added_sections: List[str]
    added_fields: Dict[str, List[str]]
    skipped_arrays: Dict[str, List[str]]


def load_toml(path: Path) -> Dict[str, Any]:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def collect_section_info() -> Dict[str, SectionInfo]:
    root = config_schema.ConfigRoot()
    sections: Dict[str, SectionInfo] = {}

    def walk(path: str, cls: type[config_schema.ConfigSection], instance: config_schema.ConfigSection) -> None:
        fields = cls.fields()
        field_order = list(fields.keys())
        defaults = {name: config_schema.ConfigSection._serialize_value(getattr(instance, name)) for name in field_order}
        sections[path] = SectionInfo(defaults=defaults, field_meta=fields, field_order=field_order)
        for name, subcls in cls.subsections().items():
            walk(f"{path}.{name}", subcls, getattr(instance, name))

    walk("params", config_schema.ParamsConfig, root.params)
    return sections


def _get_section(config: Dict[str, Any], path: str) -> Optional[Dict[str, Any]]:
    current: Any = config
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    if isinstance(current, dict):
        return current
    return None


def _is_list_of_dict(meta: config_schema.FieldMeta, default: Any) -> bool:
    origin = get_origin(meta.type)
    if origin in (list, List):
        args = get_args(meta.type)
        if args:
            arg = args[0]
            if arg in (dict, Dict):
                return True
            if get_origin(arg) in (dict, Dict):
                return True
    if isinstance(default, list) and any(isinstance(item, dict) for item in default):
        return True
    return False


def inspect_config(config: Dict[str, Any], sections: Dict[str, SectionInfo]) -> DoctorReport:
    missing_sections: List[str] = []
    missing_fields: Dict[str, List[str]] = {}
    invalid_choices: List[InvalidChoice] = []

    for path, info in sections.items():
        section = _get_section(config, path)
        if section is None:
            missing_sections.append(path)
            continue
        missing: List[str] = []
        for field_name in info.field_order:
            meta = info.field_meta[field_name]
            if field_name not in section:
                missing.append(field_name)
                continue
            if meta.choices:
                value = section[field_name]
                if isinstance(value, list):
                    bad = [v for v in value if v not in meta.choices]
                    if bad:
                        invalid_choices.append(InvalidChoice(path, field_name, value, list(meta.choices)))
                else:
                    if value not in meta.choices:
                        invalid_choices.append(InvalidChoice(path, field_name, value, list(meta.choices)))
        if missing:
            missing_fields[path] = missing

    missing_sections.sort()
    return DoctorReport(
        missing_sections=missing_sections,
        missing_fields=missing_fields,
        invalid_choices=invalid_choices,
    )


def _escape_string(value: str) -> str:
    return (
        value.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\n", "\\n")
        .replace("\t", "\\t")
        .replace("\r", "\\r")
    )


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return str(value)
    if isinstance(value, str):
        return f'"{_escape_string(value)}"'
    if isinstance(value, list):
        return "[" + ", ".join(_format_value(item) for item in value) + "]"
    if value is None:
        return '""'
    return f'"{_escape_string(str(value))}"'


def _render_table(path: str, fields: Dict[str, Any], order: Iterable[str]) -> str:
    lines = [f"[{path}]"]
    for key in order:
        if key not in fields:
            continue
        lines.append(f"{key} = {_format_value(fields[key])}")
    return "\n".join(lines)


def _render_array_tables(path: str, entries: List[Dict[str, Any]]) -> List[str]:
    blocks: List[str] = []
    for entry in entries:
        lines = [f"[[{path}]]"]
        for key, value in entry.items():
            lines.append(f"{key} = {_format_value(value)}")
        blocks.append("\n".join(lines))
    return blocks


def build_append_plan(
    report: DoctorReport,
    sections: Dict[str, SectionInfo],
    append_missing_keys: bool = False,
    include_arrays: bool = False,
) -> AppendPlan:
    blocks: List[str] = []
    added_sections: List[str] = []
    added_fields: Dict[str, List[str]] = {}
    skipped_arrays: Dict[str, List[str]] = {}

    for path in report.missing_sections:
        info = sections[path]
        table_fields: Dict[str, Any] = {}
        array_fields: Dict[str, List[Dict[str, Any]]] = {}
        for field_name in info.field_order:
            default = info.defaults[field_name]
            meta = info.field_meta[field_name]
            if _is_list_of_dict(meta, default) and isinstance(default, list) and default:
                array_fields[field_name] = default
            else:
                table_fields[field_name] = default
        if table_fields:
            blocks.append(_render_table(path, table_fields, info.field_order))
        added_sections.append(path)
        if array_fields:
            if include_arrays:
                for field_name, entries in array_fields.items():
                    blocks.extend(_render_array_tables(f"{path}.{field_name}", entries))
            else:
                skipped_arrays[path] = list(array_fields.keys())

    if append_missing_keys:
        for path, fields in report.missing_fields.items():
            if path in report.missing_sections:
                continue
            info = sections[path]
            table_fields: Dict[str, Any] = {}
            array_fields: Dict[str, List[Dict[str, Any]]] = {}
            for field_name in fields:
                default = info.defaults[field_name]
                meta = info.field_meta[field_name]
                if _is_list_of_dict(meta, default) and isinstance(default, list) and default:
                    array_fields[field_name] = default
                else:
                    table_fields[field_name] = default
            if table_fields:
                blocks.append(_render_table(path, table_fields, info.field_order))
                added_fields[path] = list(table_fields.keys())
            if array_fields:
                if include_arrays:
                    for field_name, entries in array_fields.items():
                        blocks.extend(_render_array_tables(f"{path}.{field_name}", entries))
                    added_fields.setdefault(path, [])
                    added_fields[path].extend(array_fields.keys())
                else:
                    skipped_arrays.setdefault(path, []).extend(array_fields.keys())

    return AppendPlan(
        blocks=blocks,
        added_sections=added_sections,
        added_fields=added_fields,
        skipped_arrays=skipped_arrays,
    )


def format_report(
    config_path: Path,
    report: DoctorReport,
    plan: AppendPlan,
    append_missing_keys: bool,
    include_arrays: bool,
) -> str:
    lines: List[str] = []
    lines.append(f"Config doctor report: {config_path}")
    if report.missing_sections:
        lines.append(f"Missing sections ({len(report.missing_sections)}):")
        for section in report.missing_sections:
            lines.append(f"  - {section}")
    else:
        lines.append("Missing sections: none")

    missing_keys_count = sum(len(v) for v in report.missing_fields.values())
    if missing_keys_count:
        lines.append(f"Missing fields ({missing_keys_count}):")
        for section in sorted(report.missing_fields):
            fields = ", ".join(report.missing_fields[section])
            lines.append(f"  - {section}: {fields}")
    else:
        lines.append("Missing fields: none")

    if report.invalid_choices:
        lines.append(f"Invalid choices ({len(report.invalid_choices)}):")
        for item in report.invalid_choices:
            choices = ", ".join(str(choice) for choice in item.choices)
            lines.append(f"  - {item.section}.{item.field} = {item.value} (allowed: {choices})")
    else:
        lines.append("Invalid choices: none")

    if plan.added_sections:
        lines.append(f"Will append sections ({len(plan.added_sections)}):")
        for section in plan.added_sections:
            lines.append(f"  - {section}")
    else:
        lines.append("Will append sections: none")

    if append_missing_keys and plan.added_fields:
        lines.append("Will append missing fields:")
        for section in sorted(plan.added_fields):
            fields = ", ".join(plan.added_fields[section])
            lines.append(f"  - {section}: {fields}")
    elif append_missing_keys:
        lines.append("Will append missing fields: none")
    else:
        lines.append("Will append missing fields: disabled")

    if plan.skipped_arrays:
        lines.append("Skipped array-of-table defaults:")
        for section in sorted(plan.skipped_arrays):
            fields = ", ".join(plan.skipped_arrays[section])
            suffix = "" if include_arrays else " (use --include-arrays to append)"
            lines.append(f"  - {section}: {fields}{suffix}")

    return "\n".join(lines)


def append_blocks(original: str, blocks: List[str]) -> str:
    if not blocks:
        return original
    chunks = [original.rstrip(), "", "# --- Added by config doctor ---", ""]
    chunks.extend(blocks)
    chunks.append("")
    return "\n".join(chunks)
