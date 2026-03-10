"""HTML report generation utilities for tackle2 results directories."""

from .generator import ReportGenerationError, generate_report  # noqa: F401
from .llm import AgentApiConfig, AgentApiSummarizer, OllamaConfig, OllamaSummarizer, SummaryRequest  # noqa: F401
from .prompts import list_prompt_specs, resolve_prompt  # noqa: F401
from .summary_store import default_summary_dir, generate_and_store_summaries  # noqa: F401

__all__ = [
    "generate_report",
    "ReportGenerationError",
    "AgentApiConfig",
    "AgentApiSummarizer",
    "OllamaSummarizer",
    "OllamaConfig",
    "SummaryRequest",
    "list_prompt_specs",
    "resolve_prompt",
    "generate_and_store_summaries",
    "default_summary_dir",
]
