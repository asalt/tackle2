import sys
import os
import re
import glob
import shutil
import subprocess
import functools
import logging

import pathlib
from pathlib import Path
import click
from collections import defaultdict

APP_NAME = "tackle2"

# from concurrent.futures import ProcessPoolExecutor, as_completed
# from tqdm import tqdm
import pandas as pd

from . import export_packager
from . import config_schema
from . import config_doctor
from .report.generator import ReportGenerationError, generate_report, serve_directory
from .report.llm import OllamaConfig, OllamaSummarizer

# from concurrent.futures import ThreadPoolExecutor, as_completed

# import pyfaidx
# from pyfaidx import Fasta

# import janitor

# from . import log
# from . import io
# from . import io_external
# from . import modisite
# from .utils import data_generator
# from .constants import VALID_MODI_COLS, get_all_columns
# from .runner import run_pipeline
# from . import mapper
# from . import reduce
logger = logging.getLogger(name=__file__)


# CLI_R = os.path.abspath((os.path.join(os.path.dirname(__file__), "cli.R")))
CLI_R = pathlib.Path(__file__).parent.parent / "R" / "cli.R"
if not CLI_R.exists():
    raise FileNotFoundError(f"{CLI_R} does not exist")


@click.group(chain=True)
def main():
    pass


def _resolve_example(path_components: tuple[str, ...]) -> pathlib.Path:
    base = pathlib.Path(__file__).resolve().parent.parent
    candidate = base.joinpath(*path_components)
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"{candidate} does not exist")


def _copy_example(source: pathlib.Path, destination: pathlib.Path) -> bool:
    if destination.exists():
        click.echo(f"File {destination} already exists, not overwriting")
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(source, destination)
    click.echo(f"Writing {destination}")
    return True


@main.command()
@click.option("-n", "--name", default=f"{APP_NAME}.toml", help="Output filename for the example TOML config.")
@click.option(
    "--include-colormap/--skip-colormap",
    default=False,
    show_default=True,
    help="Copy the example colormap JSON alongside the TOML config.",
)
@click.option(
    "--colormap-name",
    default="colormap.example.json",
    help="Output filename for the optional colormap example.",
)
def get_config(name, include_colormap, colormap_name):
    """Copy example configuration files into the current directory."""

    config_file = _resolve_example(("config", "base.toml"))
    if not name.endswith(".toml"):
        name = name + ".toml"
    config_dest = pathlib.Path.cwd() / name

    copied_any = _copy_example(config_file, config_dest)

    if include_colormap:
        colormap_source = _resolve_example(("config", "colormap.example.json"))
        colormap_dest = pathlib.Path.cwd() / colormap_name
        copied = _copy_example(colormap_source, colormap_dest)
        copied_any = copied_any or copied

    if copied_any:
        click.echo("done")


@main.command("config-doctor")
@click.option(
    "-c",
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="TOML configuration file to inspect.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output TOML file (default: add .doctor.toml to the config name).",
)
@click.option("--in-place", is_flag=True, help="Append updates directly to the input file.")
@click.option("--dry-run", is_flag=True, help="Report missing values without writing output.")
@click.option("--check", is_flag=True, help="Alias for --dry-run.")
@click.option(
    "--append-missing-keys/--skip-missing-keys",
    default=True,
    show_default=True,
    help="Append missing keys for existing sections (reopens tables at the end).",
)
@click.option(
    "--include-arrays",
    is_flag=True,
    help="Append defaults for array-of-table fields (e.g., genesets) when missing.",
)
@click.option("--stdout", "write_stdout", is_flag=True, help="Write updated TOML to stdout.")
def config_doctor_cmd(
    config_path: Path,
    output: Path | None,
    in_place: bool,
    dry_run: bool,
    check: bool,
    append_missing_keys: bool,
    include_arrays: bool,
    write_stdout: bool,
):
    """Check a TOML config and append missing sections to a doctor file."""

    if in_place and output is not None:
        raise click.UsageError("Use --in-place or --output, not both.")

    dry_run = dry_run or check
    sections = config_doctor.collect_section_info()
    config_data = config_doctor.load_toml(config_path)
    report = config_doctor.inspect_config(config_data, sections)
    plan = config_doctor.build_append_plan(
        report,
        sections,
        append_missing_keys=append_missing_keys,
        include_arrays=include_arrays,
    )
    click.echo(
        config_doctor.format_report(
            config_path,
            report,
            plan,
            append_missing_keys=append_missing_keys,
            include_arrays=include_arrays,
        )
    )

    if dry_run:
        return

    if write_stdout:
        updated = config_doctor.append_blocks(config_path.read_text(), plan.blocks)
        click.echo(updated)
        return

    if in_place:
        output_path = config_path
    else:
        if output is None:
            if config_path.name.endswith(".toml"):
                output_path = config_path.with_name(f"{config_path.stem}.doctor.toml")
            else:
                output_path = config_path.with_name(f"{config_path.name}.doctor.toml")
        else:
            output_path = output

    if not plan.blocks:
        click.echo("No sections to append; no output written.")
        return

    updated = config_doctor.append_blocks(config_path.read_text(), plan.blocks)
    output_path.write_text(updated)
    click.echo(f"Wrote {output_path}")


@main.command()
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False),
    help=".toml file with additional parameters for report",
)
@click.option(
    "-i",
    "--interactive",
    type=bool,
    is_flag=True,
    default=True,
    show_default=True,
    help="run in interactive session within python rpy2. button has no effect. always on",
)
@click.option("-v", "--verbose", type=bool, is_flag=True, help="verbose output")
@click.option(
    "-l",
    "--log_level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
    help="log level",
)
def run(config, interactive, verbose, log_level):


    pwd = pathlib.Path('.').absolute()
    logging.info(f"{pwd}")
    config = pathlib.Path(config).absolute()
    if not config.exists():
        raise ValueError(f"Config file {str(config)} does not exist")
    import rpy2.robjects as robjects


    robjects.r(
        f"""
        message("sourcing run file")
        setwd('{str(CLI_R.parent)}')
        """
    )

    robjects.r(
        f"""
        source('{str(CLI_R)}')
        source('utils.R')  # for clean_args
        """
    )
    robjects.r.assign("pwd", str(pwd))
    logger.info(f"Assignning pwd to {pwd} in R environ")
    logger.info(f"trying to load {str(config)} with RcppTOML")
    try:
        robjects.r(
            f"""
            message("trying to load {str(config)}")
            params <- RcppTOML::parseTOML("{str(config)}")
            message("cleaning params")
            cleaned_params <- clean_args(params$params, root_dir = pwd)
            message("running")
            run(cleaned_params)
            """
        )
    except Exception as e:
        print(robjects.r("rlang::last_trace()"))



@main.command()
@click.option(
    "-p",
    "--port",
    type=int,
    default=8765,
    help="port to run the server on",
    show_default=True,
)
def launch_assistant(port):
    from .persistent_llm_server import start_server, run_server

    run_server(port=port)


@main.command()
@click.option(
    "-s",
    "--savedir",
    type=click.Path(exists=True, file_okay=False),
    help=f"Directory containing {APP_NAME} pipeline outputs",
)
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False),
    help="TOML configuration file used for the run",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False),
    help="Target ZIP file path for the packaged export",
)
@click.option("--include-cache/--exclude-cache", default=False, show_default=True)
@click.option("--include-ranks/--exclude-ranks", default=False, show_default=True)
@click.option("--include-hashes/--skip-hashes", default=True, show_default=True)
@click.option(
    "--label",
    type=str,
    help="Custom root directory name inside the archive",
)
@click.option(
    "--split-components/--single-archive",
    default=False,
    show_default=True,
    help="Package each top-level result folder into its own zip inside an output directory",
)
def package(
    savedir,
    config,
    output,
    include_cache,
    include_ranks,
    include_hashes,
    label,
    split_components,
):
    """Package analysis outputs into a distributable archive."""

    if not savedir and not config:
        raise click.UsageError("Provide either --savedir or --config")

    try:
        archive_path, manifest = export_packager.package_results(
            savedir=Path(savedir) if savedir else None,
            config_path=Path(config) if config else None,
            output_path=Path(output) if output else None,
            include_cache=include_cache,
            include_ranks=include_ranks,
            include_hashes=include_hashes,
            export_label=label,
            split_components=split_components,
        )
    except export_packager.PackagingError as exc:
        raise click.ClickException(str(exc))

    descriptor = "package directory" if split_components else "archive"
    click.echo(f"Created {descriptor}: {archive_path}")
    click.echo(f"Files packaged: {manifest['file_count']}")
    click.echo(f"Total size: {manifest['total_size_bytes']:,} bytes")


@main.command()
@click.option(
    "-s",
    "--savedir",
    type=click.Path(exists=True, file_okay=False),
    help=f"Directory containing {APP_NAME} pipeline outputs",
)
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False),
    help="TOML configuration file used for the run",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(file_okay=False, dir_okay=True),
    help="Directory where the report should be written (defaults to <savedir>/report)",
)
@click.option("--serve/--no-serve", default=False, show_default=True, help="Launch a static web server after generation")
@click.option("--port", default=8000, show_default=True, help="Port for --serve")
@click.option("--no-browser", is_flag=True, help="Do not auto-open the browser when serving")
@click.option(
    "--log-level",
    type=click.Choice(["debug", "info", "warning", "error", "critical"], case_sensitive=False),
    default="warning",
    show_default=True,
    help="Logging verbosity for report generation.",
)
@click.option("--llm-model", default=None, help="Ollama model for AI-generated summaries (e.g., 'llama3.2:3b').")
@click.option("--llm-host", default=None, help="Override Ollama host URL (default uses OLLAMA_HOST or localhost).")
@click.option("--llm-temperature", default=0.2, type=float, show_default=True, help="Sampling temperature for the LLM.")
@click.option("--llm-max-tokens", default=512, type=int, show_default=True, help="Maximum tokens requested per summary.")
@click.option("--llm-run-limit", default=6, type=int, show_default=True, help="Collections considered for the overall summary (0 = all).")
@click.option("--llm-collections/--no-llm-collections", default=True, show_default=True, help="Generate per-collection summaries.")
@click.option("--llm-collection-limit", default=3, type=int, show_default=True, help="Maximum collections with individual summaries (0 = all).")
@click.option("--llm-comparisons/--no-llm-comparisons", default=False, show_default=True, help="Generate per-comparison summaries (experimental).")
@click.option("--llm-comparison-limit", default=0, type=int, show_default=True, help="Limit per-collection comparison summaries when enabled (0 = all).")
@click.option("--llm-comparisons-per-collection", default=3, type=int, show_default=True, help="Comparisons per collection included in prompts (0 = all).")
@click.option("--llm-pathway-limit", default=3, type=int, show_default=True, help="Top pathways per comparison shared with the LLM (0 = all).")
def report(
    savedir,
    config,
    output,
    serve,
    port,
    no_browser,
    log_level,
    llm_model,
    llm_host,
    llm_temperature,
    llm_max_tokens,
    llm_run_limit,
    llm_collections,
    llm_collection_limit,
    llm_comparisons,
    llm_comparison_limit,
    llm_comparisons_per_collection,
    llm_pathway_limit,
):
    """Generate an HTML report summarising analysis outputs."""

    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.WARNING),
        format="%(levelname)s:%(name)s:%(message)s",
        force=True,
    )

    summarizer = None
    if llm_model:
        try:
            summarizer_config = OllamaConfig(
                model=llm_model,
                host=llm_host,
                temperature=llm_temperature,
                max_tokens=max(0, llm_max_tokens),
                run_collection_limit=max(0, llm_run_limit),
                enable_collection_summaries=llm_collections,
                collection_summary_limit=max(0, llm_collection_limit),
                enable_comparison_summaries=llm_comparisons,
                comparison_summary_limit=max(0, llm_comparison_limit),
                comparisons_per_collection=max(0, llm_comparisons_per_collection),
                pathways_per_comparison=max(0, llm_pathway_limit),
            )
            summarizer = OllamaSummarizer(summarizer_config)
        except Exception as exc:  # pragma: no cover - defensive
            raise click.ClickException(f"Failed to initialise LLM summariser: {exc}") from exc

    try:
        index_path = generate_report(
            savedir=Path(savedir) if savedir else None,
            config=Path(config) if config else None,
            output=Path(output) if output else None,
            summarizer=summarizer,
        )
    except ReportGenerationError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Report written to {index_path}")

    if serve:
        try:
            serve_directory(index_path.parent, port=port, open_browser=not no_browser)
        except ReportGenerationError as exc:
            raise click.ClickException(str(exc)) from exc

def _serialize_default(value):
    if isinstance(value, (dict, list)):
        import json

        return json.dumps(value, indent=2, sort_keys=True)
    return value


@main.command()
@click.argument("section", required=False)
@click.option("--json", "as_json", is_flag=True, help="Emit JSON instead of text tables")
def describe(section, as_json):
    """Describe configuration sections and defaults."""

    try:
        descriptor = config_schema.describe_section(section)
    except KeyError as exc:
        raise click.ClickException(str(exc)) from exc

    if as_json:
        import json

        click.echo(json.dumps(config_schema.section_to_dict(descriptor), indent=2, sort_keys=True))
        return

    header = descriptor.path
    click.echo(f"[{header}]")
    click.echo(descriptor.description or "")
    click.echo("")

    if descriptor.fields:
        key_width = max(len(field.name) for field in descriptor.fields)
        type_width = max(len(field.type) for field in descriptor.fields)
        for field_info in descriptor.fields:
            default_value = _serialize_default(field_info.default)
            click.echo(
                f"{field_info.name.ljust(key_width)}  {field_info.type.ljust(type_width)}  default={default_value}"
            )
            if field_info.description:
                click.echo(f"  {field_info.description}")
            if field_info.choices:
                click.echo(f"  choices: {', '.join(map(str, field_info.choices))}")
            click.echo("")
    else:
        click.echo("(no options)")
        click.echo("")

    if descriptor.subsections:
        click.echo("Subsections:")
        for subsection in descriptor.subsections:
            suffix = " (array)" if subsection.is_array else ""
            click.echo(f"  - {subsection.path}{suffix}: {subsection.description}")
