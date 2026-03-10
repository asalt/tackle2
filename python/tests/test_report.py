import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from python.cli import main
from python.report.generator import generate_report
from python.report import assets, catalog, summary as report_summary


ROOT = Path(__file__).resolve().parents[2]
GENERATE_SCRIPT = ROOT / "scripts" / "generate_demo_savedir.R"


def _require_rscript() -> str:
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript not available; skipping report integration test")
    if not GENERATE_SCRIPT.exists():
        pytest.skip("Demo savedir script missing")
    return rscript


def _build_savedir(tmp_path: Path) -> Path:
    rscript = _require_rscript()
    destination = tmp_path / "demo"
    result = subprocess.run(
        [rscript, str(GENERATE_SCRIPT), str(destination)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        fixture_savedir = ROOT / "R" / "tests" / "integration_tests" / "test_output"
        if fixture_savedir.exists():
            fallback = destination / "fixture_savedir"
            shutil.copytree(fixture_savedir, fallback)
            return fallback
        pytest.skip(f"Failed to generate demo savedir: {result.stderr}")

    candidates = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not candidates:
        raise AssertionError("generate_demo_savedir produced no output path")

    savedir = Path(candidates[-1])
    if not savedir.exists():
        raise AssertionError(f"Reported savedir does not exist: {savedir}")

    return savedir


def _build_minimal_savedir(tmp_path: Path) -> Path:
    savedir = tmp_path / "mini_savedir"
    table_dir = savedir / "gsea_tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    (savedir / "run.log").write_text("example log\n", encoding="utf-8")
    (table_dir / "H__group_A.tsv").write_text(
        "\n".join(
            [
                "pathway\tNES\tpadj\tpeak_rank_pct\tleading_edge_fraction\tleading_edge_span_pct\tfront_loaded_score\tleadingEdge_genesymbol",
                "HALLMARK_APOPTOSIS\t2.0186087352224678\t0.001234567\t0.1254321\t0.3759876\t0.06254321\t16.812345\tTP53/BAX/CASP3",
                "HALLMARK_TGF_BETA_SIGNALING\t1.423456789\t0.020456789\t0.3409876\t0.2504321\t0.08098765\t4.17654321\tSMAD2/SMAD3/TGFBR1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return savedir


def _write_minimal_table(table_path: Path) -> None:
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table_path.write_text(
        "\n".join(
            [
                "pathway\tNES\tpadj\tleadingEdge_genesymbol",
                "HALLMARK_APOPTOSIS\t2.0\t0.001\tTP53/BAX/CASP3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_combined_table(table_path: Path) -> None:
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table_path.write_text(
        "\n".join(
            [
                "pathway\tpval_group_A\tpadj_group_A\tNES_group_A\tmainpathway_group_A",
                "HALLMARK_APOPTOSIS\t0.001\t0.01\t2.0\tTRUE",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


@pytest.mark.integration
def test_generate_report_end_to_end(tmp_path):
    savedir = _build_savedir(tmp_path)
    output_dir = tmp_path / "report"
    index_path = generate_report(savedir=savedir, output=output_dir)

    assert index_path.exists()
    css_path = index_path.parent / "static" / "report.css"
    assert css_path.exists()

    html = index_path.read_text(encoding="utf-8")
    assert "tackle2 Summary" in html
    assert "HALLMARK" in html
    assert "Pages" in html

    collection_dir = index_path.parent / "collections"
    assert collection_dir.exists()
    detail_pages = list(collection_dir.glob("*.html"))
    assert detail_pages, "Expected at least one collection detail page"
    detail_html = detail_pages[0].read_text(encoding="utf-8")
    assert "Back to Summary" in detail_html

    comparison_dirs = list(collection_dir.glob("*/"))
    assert comparison_dirs, "Expected comparison subdirectories"
    first_comparison = sorted((comparison_dirs[0]).glob("*.html"))
    assert first_comparison, "Expected comparison detail pages"
    comparison_html = first_comparison[0].read_text(encoding="utf-8")
    assert "Back to" in comparison_html


@pytest.mark.integration
def test_cli_report_command(tmp_path):
    savedir = _build_savedir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["report", "--savedir", str(savedir), "--output", str(tmp_path / "web")],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "Report written to" in result.output
    index_path = tmp_path / "web" / "index.html"
    assert index_path.exists()
    collection_dir = index_path.parent / "collections"
    assert any(collection_dir.glob("*.html"))
    assert any(collection_dir.glob("*/*.html"))


def test_prepare_plot_previews_handles_errors(tmp_path):
    dummy_pdf = tmp_path / "plot.pdf"
    previews = assets.prepare_plot_previews([dummy_pdf], tmp_path, tmp_path)
    assert previews == {}


def test_prepare_plot_previews_prioritises_pdfs(tmp_path, monkeypatch):
    pdf_path = tmp_path / "plot.pdf"
    png_path = tmp_path / "plot.png"
    pdf_path.touch()
    png_path.touch()

    calls = []

    def fake_run(cmd, check, capture_output):
        from pathlib import Path
        calls.append(cmd)
        prefix = Path(cmd[-1])
        prefix.with_suffix(".png").write_bytes(b"")
        import subprocess

        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(assets.subprocess, "run", fake_run)
    monkeypatch.setattr(assets.shutil, "which", lambda name: "/usr/bin/pdftoppm")

    previews = assets.prepare_plot_previews(
        [png_path, pdf_path],
        savedir=tmp_path,
        output_dir=tmp_path,
        limit=1,
    )

    assert len(previews) == 1
    key = next(iter(previews))
    assert key == pdf_path.relative_to(tmp_path)
    assert (tmp_path / previews[key]).exists()
    assert len(calls) == 1


def test_report_summaries_command_writes_json(tmp_path, monkeypatch):
    savedir = _build_minimal_savedir(tmp_path)
    summary_dir = tmp_path / "ai"
    calls = {"count": 0}

    class _FakeResponse:
        def __init__(self, payload):
            self._payload = payload

        def read(self):
            return json.dumps(self._payload).encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    def fake_urlopen(request, timeout):
        _ = request, timeout
        calls["count"] += 1
        return _FakeResponse(
            {
                "sessionId": "test-session",
                "message": json.dumps(
                    {
                        "summary": "Run summary from agent.",
                        "bullets": ["H collection has significant pathways."],
                    }
                ),
            }
        )

    monkeypatch.setattr("python.report.llm.urlopen", fake_urlopen)
    monkeypatch.setenv("TACKLE_AGENT_API", "http://agent.test")
    monkeypatch.setenv("TACKLE_AGENT_MODEL_LABEL", "agent-test")

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "report-summaries",
            "--savedir",
            str(savedir),
            "--output",
            str(summary_dir),
            "--include-comparisons",
            "--include-leading-edge",
            "--comparison-limit",
            "1",
            "--pathway-limit",
            "8",
            "--leading-edge-limit",
            "1",
            "--leading-edge-gene-limit",
            "4",
            "--prompt-note",
            "Focus on leading-edge genes when they are available.",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "Summary JSON written to" in result.output
    assert calls["count"] >= 4

    manifest = json.loads((summary_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["generated"] >= 2
    assert manifest["entries"][0]["prompt_variant"] == "default"

    run_doc = json.loads((summary_dir / "run.json").read_text(encoding="utf-8"))
    assert run_doc["summary"]["text"] == "Run summary from agent."
    assert "prompt_text" in run_doc["prompt"]
    assert run_doc["prompt"]["id"] == "run_overview"
    assert run_doc["prompt"]["version"] == "v1"
    assert run_doc["prompt"]["variant"] == "default"
    assert run_doc["prompt"]["user_notes"] == ["Focus on leading-edge genes when they are available."]

    collection_docs = list((summary_dir / "collections").glob("*.json"))
    assert collection_docs, "expected collection summary json"
    comparison_docs = list((summary_dir / "comparisons").glob("*/*.json"))
    assert comparison_docs, "expected comparison summary json"
    comparison_doc = json.loads(comparison_docs[0].read_text(encoding="utf-8"))
    top_pathway_payload = comparison_doc["prompt"]["payload"]["comparison"]["top_pathways"][0]
    assert top_pathway_payload["nes"] == pytest.approx(2.0186)
    assert top_pathway_payload["padj"] == pytest.approx(0.0012346)
    assert top_pathway_payload["peak_rank_pct"] == pytest.approx(0.12543)
    assert top_pathway_payload["leading_edge_fraction"] == pytest.approx(0.37599)
    assert top_pathway_payload["leading_edge_span_pct"] == pytest.approx(0.062543)
    assert top_pathway_payload["front_loaded_score"] == pytest.approx(16.812)
    leading_edge_docs = list((summary_dir / "leading_edge").glob("*/*/*.json"))
    assert leading_edge_docs, "expected leading-edge summary json"
    leading_edge_doc = json.loads(leading_edge_docs[0].read_text(encoding="utf-8"))
    assert leading_edge_doc["prompt"]["id"] == "leading_edge_summary"
    assert leading_edge_doc["prompt"]["payload"]["context"] == "leading_edge"
    assert leading_edge_doc["prompt"]["payload"]["geneset"] == "HALLMARK_APOPTOSIS"
    assert leading_edge_doc["prompt"]["payload"]["nes"] == pytest.approx(2.0186)
    assert leading_edge_doc["prompt"]["payload"]["padj"] == pytest.approx(0.0012346)
    assert leading_edge_doc["prompt"]["payload"]["peak_rank_pct"] == pytest.approx(0.12543)
    assert leading_edge_doc["prompt"]["payload"]["leading_edge_fraction"] == pytest.approx(0.37599)
    assert leading_edge_doc["prompt"]["payload"]["leading_edge_span_pct"] == pytest.approx(0.062543)
    assert leading_edge_doc["prompt"]["payload"]["front_loaded_score"] == pytest.approx(16.812)
    assert len(leading_edge_doc["prompt"]["payload"]["items"]) == 3


def test_build_context_derives_collection_names_from_savedir_layout(tmp_path):
    savedir = tmp_path / "layout_savedir"
    (savedir / "C5_All").mkdir(parents=True, exist_ok=True)
    (savedir / "C2_CP_KEGG_MEDICUS").mkdir(parents=True, exist_ok=True)
    _write_minimal_table(savedir / "gsea_tables" / "C5_All_group_A.tsv")
    _write_minimal_table(savedir / "gsea_tables" / "C2_CP_KEGG_MEDICUS_group_B.tsv")

    artefacts = catalog.scan_savedir(savedir)
    context = report_summary.build_context(savedir, artefacts)
    by_collection = {collection.identifier: collection for collection in context["collections"]}

    assert sorted(by_collection) == ["C2_CP_KEGG_MEDICUS", "C5_All"]
    assert by_collection["C5_All"].comparisons[0].identifier == "group_A"
    assert by_collection["C2_CP_KEGG_MEDICUS"].comparisons[0].identifier == "group_B"


def test_build_context_uses_config_collection_names_when_layout_is_ambiguous(tmp_path):
    savedir = tmp_path / "config_savedir"
    _write_minimal_table(savedir / "gsea_tables" / "C2_CP_KEGG_MEDICUS_group_A.tsv")
    config_path = tmp_path / "report.toml"
    config_path.write_text(
        "\n".join(
            [
                "[params]",
                "savedir = 'ignored'",
                "",
                "[[params.genesets]]",
                "category = 'C2'",
                "subcategory = 'CP:KEGG_MEDICUS'",
                "collapse = true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    artefacts = catalog.scan_savedir(savedir)
    context = report_summary.build_context(savedir, artefacts, config_path=config_path)

    assert len(context["collections"]) == 1
    collection = context["collections"][0]
    assert collection.identifier == "C2_CP_KEGG_MEDICUS"
    assert collection.comparisons[0].identifier == "group_A"


def test_build_context_ignores_combined_all_table_when_real_comparisons_exist(tmp_path):
    savedir = tmp_path / "combined_savedir"
    (savedir / "C2_CP_KEGG_MEDICUS").mkdir(parents=True, exist_ok=True)
    _write_minimal_table(savedir / "gsea_tables" / "C2_CP_KEGG_MEDICUS_group_A.tsv")
    _write_combined_table(savedir / "gsea_tables" / "C2_CP_KEGG_MEDICUS_all.tsv")

    artefacts = catalog.scan_savedir(savedir)
    context = report_summary.build_context(savedir, artefacts)

    assert context["stats"]["tables"] == 1
    assert len(context["collections"]) == 1
    collection = context["collections"][0]
    assert collection.identifier == "C2_CP_KEGG_MEDICUS"
    assert [comparison.identifier for comparison in collection.comparisons] == ["group_A"]


def test_report_consumes_stored_summary_json(tmp_path):
    savedir = _build_minimal_savedir(tmp_path)
    summary_dir = savedir / "report" / "ai"
    summary_dir.mkdir(parents=True, exist_ok=True)
    (summary_dir / "run.json").write_text(
        json.dumps(
            {
                "scope": "run",
                "key": "run",
                "target": {"scope": "run"},
                "prompt": {"prompt_sha256": "abc"},
                "summary": {
                    "text": "Precomputed run summary.",
                    "bullets": ["Stored summary bullet."],
                    "model": "agent-test",
                },
            }
        ),
        encoding="utf-8",
    )
    leading_edge_dir = summary_dir / "leading_edge" / "h" / "group-a"
    leading_edge_dir.mkdir(parents=True, exist_ok=True)
    (leading_edge_dir / "hallmark-apoptosis.json").write_text(
        json.dumps(
            {
                "scope": "leading_edge",
                "key": "H:group_A:HALLMARK_APOPTOSIS",
                "target": {
                    "scope": "leading_edge",
                    "collection_id": "H",
                    "comparison_id": "group_A",
                    "pathway_name": "HALLMARK_APOPTOSIS",
                },
                "prompt": {"prompt_sha256": "def"},
                "summary": {
                    "text": "Leading-edge genes point to an apoptosis-focused signal driven by TP53 and BAX.",
                    "bullets": ["TP53, BAX, and CASP3 anchor the leading-edge set."],
                    "model": "agent-test",
                },
            }
        ),
        encoding="utf-8",
    )

    output_dir = tmp_path / "report"
    index_path = generate_report(
        savedir=savedir,
        output=output_dir,
    )
    html = index_path.read_text(encoding="utf-8")
    assert "Precomputed run summary." in html
    assert "Stored summary bullet." in html
    assert 'href="ai.html"' in html
    ai_html = (output_dir / "ai.html").read_text(encoding="utf-8")
    assert "AI Summaries" in ai_html
    assert "Precomputed run summary." in ai_html
    assert "Leading-edge genes point to an apoptosis-focused signal driven by TP53 and BAX." in ai_html
    comparison_html = (output_dir / "collections" / "h" / "group-a.html").read_text(encoding="utf-8")
    assert "Leading-edge genes point to an apoptosis-focused signal driven by TP53 and BAX." in comparison_html
    assert "TP53, BAX, and CASP3 anchor the leading-edge set." in comparison_html


def test_generate_report_force_refreshes_static_assets(tmp_path):
    savedir = _build_minimal_savedir(tmp_path)
    output_dir = tmp_path / "report"
    generate_report(savedir=savedir, output=output_dir)

    css_path = output_dir / "static" / "report.css"
    original_css = css_path.read_text(encoding="utf-8")
    css_path.write_text("body { color: hotpink; }\n", encoding="utf-8")

    generate_report(savedir=savedir, output=output_dir)
    assert css_path.read_text(encoding="utf-8") == "body { color: hotpink; }\n"

    generate_report(savedir=savedir, output=output_dir, force=True)
    assert css_path.read_text(encoding="utf-8") == original_css


def test_report_summaries_preview_only_does_not_require_agent(tmp_path):
    savedir = _build_minimal_savedir(tmp_path)
    summary_dir = tmp_path / "preview_ai"
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "report-summaries",
            "--savedir",
            str(savedir),
            "--output",
            str(summary_dir),
            "--preview-only",
            "--include-leading-edge",
            "--leading-edge-limit",
            "1",
            "--leading-edge-gene-limit",
            "3",
            "--prompt-variant",
            "brief",
            "--print-prompts",
            "--prompt-note",
            "Preview the prompt only.",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "prompt_id=run_overview" in result.output
    assert "prompt_id=leading_edge_summary" in result.output
    assert "variant=brief" in result.output
    assert "HALLMARK_APOPTOSIS" in result.output
    run_doc = json.loads((summary_dir / "run.json").read_text(encoding="utf-8"))
    assert run_doc["status"] == "preview_only"
    assert run_doc["summary"] is None
    assert run_doc["prompt"]["variant"] == "brief"
    assert run_doc["prompt"]["user_notes"] == ["Preview the prompt only."]
    leading_edge_doc = json.loads(
        next((summary_dir / "leading_edge").glob("*/*/*.json")).read_text(encoding="utf-8")
    )
    assert leading_edge_doc["status"] == "preview_only"
    assert len(leading_edge_doc["prompt"]["payload"]["items"]) == 3
    comparison_docs = list((summary_dir / "comparisons").glob("*/*.json"))
    assert comparison_docs, "expected per-comparison summaries when leading-edge summaries are requested"



def test_list_summary_prompts_command():
    runner = CliRunner()
    result = runner.invoke(main, ["list-summary-prompts"], catch_exceptions=False)
    assert result.exit_code == 0, result.output
    assert "run_overview" in result.output
    assert "collection_summary" in result.output
    assert "comparison_summary" in result.output
    assert "leading_edge_summary" in result.output
