"""Python accessors for the R geneset utilities.

The helper mirrors the information in ``docs/README.md`` so users can
discover the bridge via ``help()``. Typical usage matches the README example::

    from tackle2 import get_geneset_collections, geneset_membership

    hallmark = get_geneset_collections([
        {"category": "H", "subcategory": "", "collapse": False},
    ])
    gene_lists = geneset_membership(hallmark["H_"])

``get_geneset_collections`` returns pandas ``DataFrame`` objects fetched by the
underlying R ``geneset_utils`` helpers, while ``geneset_membership`` collapses
any one collection into ``{geneset_name: [genes...]}`` mappings.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass
from functools import lru_cache
from typing import Mapping, MutableMapping, Sequence, Union

import pandas as pd

GenesetSelection = Union[Sequence[Mapping[str, object]], pd.DataFrame]

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def _require_rpy2():
    try:
        from rpy2 import robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
    except ImportError as exc:
        raise ImportError(
            "rpy2 is required for the Python geneset bridge. "
            "Install rpy2 or use the R pipeline directly."
        ) from exc
    return robjects, pandas2ri, localconverter


def _quote_path(path: pathlib.Path) -> str:
    """Return a JSON-quoted string suitable for embedding in R code."""

    return json.dumps(str(path))


def _py_df_to_r(df: pd.DataFrame):
    _, pandas2ri, localconverter = _require_rpy2()
    with localconverter(pandas2ri.converter):
        return pandas2ri.py2rpy(df)


def _r_df_to_pandas(r_df):
    _, pandas2ri, localconverter = _require_rpy2()
    with localconverter(pandas2ri.converter):
        return pandas2ri.rpy2py(r_df)


@lru_cache(maxsize=None)
def _load_geneset_environment(repo_root: pathlib.Path | None = None):
    """Load the lazy-loaded geneset_utils environment from R."""

    root = (repo_root or REPO_ROOT).resolve()
    r_dir = root / "R"
    if not r_dir.exists():
        raise FileNotFoundError(f"R directory not found at {r_dir}")

    robjects, _, _ = _require_rpy2()
    r = robjects.r
    r(f"setwd({_quote_path(root)})")
    r(f"source({_quote_path(r_dir / 'lazyloader.R')})")
    r(
        "if (!requireNamespace('jsonlite', quietly = TRUE)) "
        "stop('jsonlite package is required to bridge python geneset APIs')"
    )
    r("geneset_tools <- get_tool_env('geneset_utils')")
    return robjects.globalenv["geneset_tools"]


def _parse_geneset_specs(env, genesets: Sequence[Mapping[str, object]] | None):
    """Use the R helper to coerce python specs into a tibble."""

    specs = list(genesets or [])
    if not isinstance(specs, list):
        raise TypeError("genesets must be a sequence of dictionaries")
    if not specs:
        raise ValueError("genesets must provide at least one entry")
    for idx, entry in enumerate(specs):
        if not isinstance(entry, Mapping):
            raise TypeError(f"geneset entry at position {idx} must be a mapping")

    robjects, _, _ = _require_rpy2()
    robjects.r.assign(".geneset_spec_json", json.dumps(specs))
    parser = env["geneset_array_to_df"]
    return parser(robjects.r("jsonlite::fromJSON(.geneset_spec_json, simplifyVector = FALSE)"))


def _r_frames_to_pandas_dict(r_list) -> MutableMapping[str, pd.DataFrame]:
    result: MutableMapping[str, pd.DataFrame] = {}
    names = list(r_list.names) if r_list.names is not None else []
    for idx in range(len(r_list)):
        key = names[idx] if idx < len(names) and names[idx] else str(idx)
        result[key] = _r_df_to_pandas(r_list[idx])
    return result


def _r_named_list_to_dict(r_list) -> MutableMapping[str, list[str]]:
    result: MutableMapping[str, list[str]] = {}
    names = list(r_list.names) if r_list.names is not None else []
    for idx in range(len(r_list)):
        key = names[idx] if idx < len(names) and names[idx] else str(idx)
        r_values = r_list[idx]
        values = [str(v) for v in r_values if str(v) != "NA"]
        result[key] = values
    return result


@dataclass
class GenesetAPI:
    """Thin python wrapper around the geneset R module.

    The class loads ``R/lazyloader.R`` via :mod:`rpy2` and exposes methods that
    mirror the R ``geneset_utils`` API.  Instantiate it directly when working
    outside the repository root or rely on the module-level helpers.
    """

    repo_root: pathlib.Path | None = None

    def __post_init__(self) -> None:
        self._env = _load_geneset_environment(self.repo_root)

    def _coerce_spec_frame(self, genesets):
        if isinstance(genesets, pd.DataFrame):
            return _py_df_to_r(genesets)
        return _parse_geneset_specs(self._env, genesets)

    def get_collections(
        self,
        genesets: GenesetSelection,
        *,
        species: str = "Homo sapiens",
    ) -> MutableMapping[str, pd.DataFrame]:
        """Return a dict of pandas DataFrames for the requested selections."""
        selections = self._coerce_spec_frame(genesets)
        collector = self._env["get_collections"]
        r_list = collector(selections, species=species)
        return _r_frames_to_pandas_dict(r_list)

    def geneset_membership(self, geneset_df: pd.DataFrame) -> MutableMapping[str, list[str]]:
        """Collapse a collection DataFrame into ``{geneset_name: [genes]}``."""
        if not isinstance(geneset_df, pd.DataFrame):
            raise TypeError("geneset_df must be a pandas DataFrame")
        converter = self._env["genesets_df_to_list"]
        r_df = _py_df_to_r(geneset_df)
        return _r_named_list_to_dict(converter(r_df))


_DEFAULT_API: GenesetAPI | None = None


def _get_default_api() -> GenesetAPI:
    global _DEFAULT_API
    if _DEFAULT_API is None:
        _DEFAULT_API = GenesetAPI()
    return _DEFAULT_API


def get_geneset_collections(
    genesets: GenesetSelection,
    *,
    species: str = "Homo sapiens",
):
    """Fetch one or more msigdbr collections using the existing R helpers.

    Parameters
    ----------
    genesets:
        Either a sequence of dictionaries with ``category``, ``subcategory``,
        and ``collapse`` keys, or a pandas ``DataFrame`` in the same shape as
        the R ``geneset_utils`` inputs.
    species:
        Species string passed through to :mod:`msigdbr` (defaults to ``"Homo
        sapiens"``).

    Returns
    -------
    dict[str, pandas.DataFrame]
        Mapping of ``{collection_name: dataframe}`` where each dataframe is
        sourced from the R ``get_collections`` helper.

    Example
    -------
    >>> from tackle2 import get_geneset_collections, geneset_membership
    >>> hallmark = get_geneset_collections([
    ...     {"category": "H", "subcategory": "", "collapse": False},
    ... ])
    >>> hallmark["H_"][:2].gs_name.tolist()  # pandas view of the R tibble
    >>> geneset_membership(hallmark["H_"])

    This mirrors the README documentation so ``help()`` exposes a usable
    walkthrough without opening the docs.
    """

    return _get_default_api().get_collections(genesets, species=species)


def geneset_membership(geneset_df: pd.DataFrame):
    """Collapse a geneset dataframe into ``{geneset_name: [genes]}`` via R.

    Pass any dataframe returned by :func:`get_geneset_collections` (or produced
    by the R pipeline) and receive a python dictionary keyed by ``gs_name`` with
    the exact gene membership lists. This is the same transformation the R
    pipeline uses before running fgsea, so downstream python code can re-use
    the curated source-of-truth genesets.
    """

    return _get_default_api().geneset_membership(geneset_df)


__all__ = ["GenesetAPI", "get_geneset_collections", "geneset_membership"]
