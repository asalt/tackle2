"""Rendering helpers for tackle2 HTML reports."""

from __future__ import annotations

import shutil
from functools import lru_cache
from importlib import resources
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape


INDEX_TEMPLATE_NAME = "index.html.j2"
COLLECTION_TEMPLATE_NAME = "collection.html.j2"
COMPARISON_TEMPLATE_NAME = "comparison.html.j2"
AI_TEMPLATE_NAME = "ai_summary.html.j2"
TEMPLATE_PACKAGE = __package__


@lru_cache(maxsize=1)
def _environment() -> Environment:
    template_dir = resources.files(TEMPLATE_PACKAGE) / "templates"
    return Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )


def render_report(context: dict) -> str:
    env = _environment()
    template = env.get_template(INDEX_TEMPLATE_NAME)
    return template.render(**context)


def render_collection(context: dict) -> str:
    env = _environment()
    template = env.get_template(COLLECTION_TEMPLATE_NAME)
    return template.render(**context)


def render_comparison(context: dict) -> str:
    env = _environment()
    template = env.get_template(COMPARISON_TEMPLATE_NAME)
    return template.render(**context)


def render_ai_summary(context: dict) -> str:
    env = _environment()
    template = env.get_template(AI_TEMPLATE_NAME)
    return template.render(**context)


def install_static_assets(destination: Path, *, force: bool = False) -> None:
    base = resources.files(TEMPLATE_PACKAGE)
    static_dir = base / "static"
    target = Path(destination) / "static"
    target.mkdir(parents=True, exist_ok=True)

    for entry in static_dir.iterdir():
        dest = target / entry.name
        if entry.is_dir():
            if dest.exists():
                if force:
                    shutil.rmtree(dest)
                    shutil.copytree(entry, dest)
            else:
                shutil.copytree(entry, dest)
        else:
            if force or not dest.exists():
                shutil.copyfile(entry, dest)


__all__ = ["render_report", "render_collection", "render_comparison", "render_ai_summary", "install_static_assets"]
