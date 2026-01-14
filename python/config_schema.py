"""Configuration schema and helpers for tackle2.

Defines configuration structures using straightforward Python classes with
inheritance and shared mixins. This keeps defaults, CLI descriptions, and
future TOML generation in sync without introducing a custom DSL. **Do not**
rename existing keys without coordinating with the R pipeline.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Union, get_args, get_origin

try:
    from typing import Literal
except ImportError:  # pragma: no cover
    from typing_extensions import Literal  # type: ignore


@dataclass
class FieldMeta:
    """Metadata describing a single configuration option."""

    type: Any
    default: Any
    description: str
    choices: Optional[Sequence[Any]] = None

    def clone(self) -> Any:
        if callable(self.default):
            return self.default()
        if isinstance(self.default, (list, dict, set)):
            return copy.deepcopy(self.default)
        return self.default


def enum_field(
    enum_cls: type[Enum],
    *,
    default: Enum | str,
    description: str,
    extra_values: Optional[Sequence[str]] = None,
) -> FieldMeta:
    values = [member.value for member in enum_cls]
    if extra_values:
        values.extend(extra_values)
    if isinstance(default, Enum):
        default_value = default.value
    else:
        default_value = default
    if default_value not in values:
        raise ValueError(f"Default '{default_value}' not among allowed values {values}")
    return FieldMeta(enum_cls, default_value, description, choices=values)


class ConfigSection:
    """Base class for configuration sections."""

    __fields__: Dict[str, FieldMeta] = {}
    __subsections__: Dict[str, "type[ConfigSection]"] = {}

    def __init__(self, **overrides: Any) -> None:
        self._apply_defaults()
        self._apply_overrides(overrides)

    @classmethod
    def fields(cls) -> Dict[str, FieldMeta]:
        merged: Dict[str, FieldMeta] = {}
        for base in reversed(cls.__mro__):
            if issubclass(base, ConfigSection):
                merged.update(getattr(base, "__fields__", {}))
        return merged

    @classmethod
    def subsections(cls) -> Dict[str, "type[ConfigSection]"]:
        merged: Dict[str, "type[ConfigSection]"] = {}
        for base in reversed(cls.__mro__):
            if issubclass(base, ConfigSection):
                merged.update(getattr(base, "__subsections__", {}))
        return merged

    def _apply_defaults(self) -> None:
        for name, meta in self.fields().items():
            setattr(self, name, meta.clone())
        for name, section_cls in self.subsections().items():
            setattr(self, name, section_cls())

    def _apply_overrides(self, overrides: Dict[str, Any]) -> None:
        for key, value in overrides.items():
            if key in self.subsections():
                section_cls = self.subsections()[key]
                if isinstance(value, ConfigSection):
                    setattr(self, key, value)
                elif isinstance(value, dict):
                    setattr(self, key, section_cls(**value))
                else:
                    setattr(self, key, value)
            elif key in self.fields():
                setattr(self, key, value)
            else:
                raise AttributeError(f"Unknown configuration option '{key}' for {self.__class__.__name__}")

    def to_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {}
        seen: set[str] = set()
        for name in self.fields():
            value = getattr(self, name)
            data[name] = self._serialize_value(value)
            seen.add(name)
        for name in self.subsections():
            if name in seen:
                continue
            value = getattr(self, name)
            data[name] = self._serialize_value(value)
        return data

    @staticmethod
    def _serialize_value(value: Any) -> Any:
        if isinstance(value, ConfigSection):
            return value.to_dict()
        if isinstance(value, list):
            return [ConfigSection._serialize_value(item) for item in value]
        return value

    @classmethod
    def describe(cls, path: str, instance: "ConfigSection", is_array: bool = False) -> "SectionDescriptor":
        fields = []
        for name, meta in cls.fields().items():
            fields.append(
                FieldDescriptor(
                    name=name,
                    type=_type_to_str(meta.type),
                    default=ConfigSection._serialize_value(getattr(instance, name)),
                    description=meta.description,
                    choices=list(meta.choices) if meta.choices else None,
                )
            )

        subsections = []
        for name, subcls in cls.subsections().items():
            subsections.append(
                SubsectionDescriptor(
                    path=f"{path}.{name}",
                    description=(subcls.__doc__ or "").strip(),
                    is_array=False,
                )
            )

        return SectionDescriptor(
            path=path,
            description=(cls.__doc__ or "").strip(),
            fields=fields,
            subsections=subsections,
            is_array=is_array,
        )


class TopNSettings(ConfigSection):
    """Common plot controls for selecting top-N pathways."""

    __fields__ = {
        "limit": FieldMeta(List[int], lambda: [12, 20, 32], "Top-N sizes to plot."),
        "do_individual": FieldMeta(bool, True, "Render per-comparison plots."),
        "do_combined": FieldMeta(bool, True, "Render combined plots."),
    }


class BarplotAdvanced(ConfigSection):
    """Advanced barplot styling controls."""

    __fields__ = {
        "stroke_width": FieldMeta(float, 1.0, "Outline width for barplot panels."),
    }


class BarplotConfig(TopNSettings):
    """Settings for barplot visualisations."""

    __subsections__ = {"advanced": BarplotAdvanced}


class BubbleplotAdvanced(ConfigSection):
    """Advanced bubble plot styling controls."""

    __fields__ = {
        "stroke_width": FieldMeta(float, 0.8, "Outline width for bubble markers."),
        "stroke_alpha": FieldMeta(float, 0.55, "Outline transparency (0-1)."),
    }


class BubbleplotConfig(TopNSettings):
    """Settings for bubble plot visualisations."""

    __fields__ = {
        "glyph": FieldMeta(str, "⁕", "Glyph used to mark padj < 0.05."),
    }
    __subsections__ = {"advanced": BubbleplotAdvanced}


class EnplotConfig(ConfigSection):
    """Enrichment plot toggles."""

    __fields__ = {
        "do_individual": FieldMeta(bool, True, "Render per-contrast enrichment plots."),
        "do_combined": FieldMeta(bool, True, "Render combined enrichment plots."),
        "limit": FieldMeta(int, 10, "Top-N pathways to display."),
        "combine_by": FieldMeta(Union[str, bool], False, "Metadata column to combine by (or false)."),
    }


class HeatmapGSEAConfig(ConfigSection):
    """GSEA heatmap controls."""

    __fields__ = {
        "do": FieldMeta(bool, True, "Render GSEA heatmaps."),
        "limit": FieldMeta(List[int], lambda: [10, 20, 40, 80], "Top-N sizes per heatmap."),
        "cut_by": FieldMeta(Union[str, bool], False, "Metadata column for grouping (or false)."),
        "cluster_rows": FieldMeta(bool, True, "Cluster rows in heatmaps."),
        "cluster_columns": FieldMeta(List[bool], lambda: [False, True], "Cluster columns flags [individual, combined]."),
        "legend_include": FieldMeta(List[str], list, "Metadata columns to include in heatmap legends."),
    }


class HeatmapGeneConfig(ConfigSection):
    """Gene-level heatmap controls."""

    __fields__ = {
        "do": FieldMeta(bool, True, "Render gene-level heatmaps."),
        "limit": FieldMeta(int, 10, "Top-N genes per heatmap."),
    }


class WordcloudConfig(ConfigSection):
    """Pathway wordcloud summary controls (per collection)."""

    __fields__ = {
        "do": FieldMeta(bool, False, "Render wordcloud summaries per collection."),
        "padj_cutoff": FieldMeta(float, 0.25, "Adjusted p-value cutoff for including pathways."),
        "top_n_pathways": FieldMeta(int, 50, "Top-N pathways used to build the wordcloud."),
        "max_words": FieldMeta(int, 70, "Maximum number of words to draw."),
    }


class PCAConfig(ConfigSection):
    """Principal component analysis plotting controls."""

    __fields__ = {
        "do": FieldMeta(bool, False, "Render PCA plots."),
        "width": FieldMeta(float, 7.8, "Figure width (inches)."),
        "height": FieldMeta(float, 7.4, "Figure height (inches)."),
        "col_by": FieldMeta(str, "", "Metadata column for colour."),
        "mark_by": FieldMeta(str, "", "Metadata column for shape."),
        "top_pc": FieldMeta(int, 3, "Number of principal components to display."),
        "max_pc": FieldMeta(int, 3, "Maximum principal component index available."),
        "labSize": FieldMeta(float, 1.8, "Label size for loadings."),
        "pointSize": FieldMeta(float, 4.0, "Point size for samples."),
        "sizeLoadingsNames": FieldMeta(float, 1.4, "Font size for loading names."),
    }


class GenePCAConfig(ConfigSection):
    """Gene-level PCA plotting controls."""

    __fields__ = {
        "do": FieldMeta(bool, False, "Run PCA on the expression matrix."),
        "width": FieldMeta(float, 7.8, "Figure width (inches)."),
        "height": FieldMeta(float, 7.4, "Figure height (inches)."),
        "components": FieldMeta(int, 3, "Number of principal components to include in biplots."),
        "metadata_color": FieldMeta(List[str], list, "Metadata columns used to colour samples."),
        "metadata_shape": FieldMeta(str, "", "Metadata column used to shape samples."),
        "top_loadings": FieldMeta(int, 25, "Number of genes to highlight in loadings heatmaps."),
        "heatmap": FieldMeta(bool, True, "Whether to generate top-loading heatmaps."),
        "cluster_rows": FieldMeta(bool, True, "Cluster genes in the loadings heatmap."),
        "cluster_columns": FieldMeta(List[bool], lambda: [False, True], "Column clustering flags for loadings heatmap."),
        "cut_by": FieldMeta(str, "", "Metadata column used to split the loadings heatmap."),
        "labSize": FieldMeta(float, 1.8, "Label size for gene loadings."),
        "pointSize": FieldMeta(float, 4.0, "Point size for samples."),
        "sizeLoadingsNames": FieldMeta(float, 1.4, "Font size for loading names."),
    }


class GeneUMAPConfig(ConfigSection):
    """Gene-level UMAP controls."""

    __fields__ = {
        "do": FieldMeta(bool, False, "Run UMAP on the expression matrix."),
        "width": FieldMeta(float, 7.2, "Figure width (inches)."),
        "height": FieldMeta(float, 6.4, "Figure height (inches)."),
        "n_neighbors": FieldMeta(int, 15, "Number of neighbors used by UMAP."),
        "min_dist": FieldMeta(float, 0.1, "UMAP min_dist parameter."),
        "metric": FieldMeta(str, "euclidean", "Distance metric passed to UMAP."),
        "seed": FieldMeta(int, 42, "Random seed for reproducibility."),
        "scale": FieldMeta(bool, True, "Z-score features before embedding."),
        "metadata_color": FieldMeta(List[str], list, "Metadata columns used to colour samples."),
        "metadata_shape": FieldMeta(str, "", "Metadata column used for point shapes."),
        "variants": FieldMeta(List[Dict[str, Any]], list, "Optional list of parameter overrides."),
        "point_type": FieldMeta(str, "gene", 'Point type for UMAP output ("gene" or "sample").'),
        "rank_name": FieldMeta(str, "", "Rank/comparison used to colour genes when point_type='gene'."),
    }



class ExtraConfig(ConfigSection):
    """Extra ordering hints for plotting."""

    __fields__ = {
        "rankname_order": FieldMeta(List[str], list, "Explicit ordering for rank names."),
        "samplename_order": FieldMeta(List[str], list, "Explicit ordering for sample names."),
        "sample_order": FieldMeta(List[str], list, "Legacy alias for sample ordering."),
    }


def _default_genesets() -> List[Dict[str, Any]]:
    return [
        {"category": "H", "subcategory": "", "collapse": False},
        {"category": "C2", "subcategory": "CP:KEGG", "collapse": True},
        {"category": "C2", "subcategory": "CP:REACTOME", "collapse": True},
        {"category": "C3", "subcategory": "TFT:GTRD", "collapse": True},
        {"category": "C3", "subcategory": "MIR:MIRDB", "collapse": True},
        {"category": "C5", "subcategory": "GO:MF", "collapse": True},
        {"category": "C5", "subcategory": "GO:BP", "collapse": True},
        {"category": "C5", "subcategory": "GO:CC", "collapse": True},
        {"category": "C5", "subcategory": "All", "collapse": True},
        {"category": "C6", "subcategory": "", "collapse": True},
        {"category": "C7", "subcategory": "IMMUNESIGDB", "collapse": True},
    ]


class AdvancedConfig(ConfigSection):
    """Advanced execution flags."""

    __fields__ = {
        "parallel": FieldMeta(bool, True, "Run fgsea in parallel."),
        "cache": FieldMeta(bool, True, "Cache intermediate results."),
        "replace": FieldMeta(bool, False, "Overwrite existing files."),
        "cachedir": FieldMeta(str, "savedir", "Cache directory (or 'savedir')."),
        "logfile": FieldMeta(str, "savedir", "Log file path (or 'savedir')."),
        "loglevel": FieldMeta(str, "INFO", "Log verbosity."),
        "pivot_gsea_results": FieldMeta(bool, False, "Write wide-format GSEA tables (can be large)."),
        "quiet": FieldMeta(bool, False, "Reduce log chatter/voice prompts."),
    }


class DbConfig(ConfigSection):
    """SQLite persistence controls for GSEA outputs."""

    __fields__ = {
        "enable": FieldMeta(bool, False, "Persist GSEA outputs to SQLite."),
        "path": FieldMeta(str, "savedir", "SQLite file path (or 'savedir' to place in savedir)."),
        "write_results": FieldMeta(bool, True, "Store per-comparison GSEA results."),
        "write_ranks": FieldMeta(bool, False, "Store full rank vectors (can be large)."),
        "write_pathways": FieldMeta(bool, True, "Store pathway membership metadata."),
    }


class ModelConfig(ConfigSection):
    """Model definition for ranks_from='model'."""

    __fields__ = {
        "name": FieldMeta(str, "primary", "Human-readable identifier for the model."),
        "type": FieldMeta(str, "limma", "Backend used to generate ranks (currently only 'limma')."),
        "design": FieldMeta(str, "", "R formula describing the model matrix."),
        "contrasts": FieldMeta(List[str], list, "Limma contrast strings; optionally name=expression."),
    }


class ParamsConfig(ConfigSection):
    """Top-level configuration controlling inputs and plotting."""

    class RankSource(str, Enum):
        VOLCANO = "volcano"
        GCT = "gct"
        MODEL = "model"

    __fields__ = {
        "ranks_from": enum_field(
            RankSource,
            default="",
            description="Data source for ranks (empty to defer).",
            extra_values=[""],
        ),
        "savedir": FieldMeta(str, "", "Base directory for outputs (resolved relative to run cwd)."),
        "volcanodir": FieldMeta(str, "", "Directory containing volcano TSV files."),
        "rankfiledir": FieldMeta(str, "savedir", "Where generated rank files are written; 'savedir' creates savedir/ranks."),
        "gct_path": FieldMeta(str, "", "Path to GCT file when ranks_from='gct'."),
        "model_file": FieldMeta(str, "", "External TOML file containing a [model] block (optional)."),
        "models": FieldMeta(List[Dict[str, Any]], list, "Additional model definitions (same keys as params.model)."),
        "species": FieldMeta(str, "Homo sapiens", "Species name passed to msigdbr."),
        "zscore_emat": FieldMeta(bool, True, "Whether to z-score the expression matrix."),
        "zscore_emat_groupby": FieldMeta(Union[str, bool], False, "Metadata column to group by during z-score normalisation (or false)."),
        "cut_by": FieldMeta(str, "group", "Metadata column used to facet plots."),
        "genesets": FieldMeta(List[Dict[str, Any]], _default_genesets, "Array of gene-set selections (category/subcategory/collapse)."),
    }

    __subsections__ = {
        "barplot": BarplotConfig,
        "bubbleplot": BubbleplotConfig,
        "enplot": EnplotConfig,
        "heatmap_gsea": HeatmapGSEAConfig,
        "heatmap_gene": HeatmapGeneConfig,
        "wordcloud": WordcloudConfig,
        "pca": PCAConfig,
        "pca_gene": GenePCAConfig,
        "umap_gene": GeneUMAPConfig,
        "extra": ExtraConfig,
        "advanced": AdvancedConfig,
        "db": DbConfig,
        "model": ModelConfig,
    }


class ConfigRoot(ConfigSection):
    """Root configuration container."""

    __subsections__ = {"params": ParamsConfig}


@dataclass
class FieldDescriptor:
    name: str
    type: str
    default: Any
    description: str
    choices: Optional[List[Any]] = None


@dataclass
class SubsectionDescriptor:
    path: str
    description: str
    is_array: bool = False


@dataclass
class SectionDescriptor:
    path: str
    description: str
    fields: List[FieldDescriptor]
    subsections: List[SubsectionDescriptor]
    is_array: bool = False


def _type_to_str(tp: Any) -> str:
    origin = get_origin(tp)
    if origin is None:
        if hasattr(tp, "__name__"):
            return tp.__name__
        return str(tp)
    args = get_args(tp)
    if origin in (list, List):
        return f"list[{_type_to_str(args[0])}]"
    if str(origin) == "typing.Literal":
        return "Literal[" + ", ".join(repr(arg) for arg in args) + "]"
    if origin in (Union,):
        return " | ".join(_type_to_str(arg) for arg in args)
    return str(tp)


def describe_section(path: Optional[str] = None) -> SectionDescriptor:
    if not path:
        path = "params"

    parts = path.split(".")
    cls: type[ConfigSection] = ConfigRoot
    instance: ConfigSection | Any = ConfigRoot()
    current_path: List[str] = []
    is_array = False

    for part in parts:
        current_path.append(part)
        subsections = cls.subsections()
        fields = cls.fields()

        if part in subsections:
            cls = subsections[part]
            instance = getattr(instance, part)
            is_array = False
            continue

        if part in fields:
            meta = fields[part]
            value = getattr(instance, part)
            return SectionDescriptor(
                path=".".join(current_path),
                description=meta.description,
                fields=[
                    FieldDescriptor(
                        name=part,
                        type=_type_to_str(meta.type),
                        default=ConfigSection._serialize_value(value),
                        description=meta.description,
                        choices=list(meta.choices) if meta.choices else None,
                    )
                ],
                subsections=[],
                is_array=False,
            )

        raise KeyError(f"Unknown section path '{'.'.join(current_path)}'")

    if not isinstance(instance, ConfigSection):
        raise KeyError(f"Section '{path}' is not a structured configuration block")

    return cls.describe(".".join(current_path), instance, is_array=is_array)


def build_default_config() -> Dict[str, Any]:
    """Return the default configuration as a nested dictionary."""

    return ConfigRoot().to_dict()


def section_to_dict(section: SectionDescriptor) -> Dict[str, Any]:
    return {
        "path": section.path,
        "description": section.description,
        "is_array": section.is_array,
        "fields": [
            {
                "name": field.name,
                "type": field.type,
                "default": field.default,
                "description": field.description,
                "choices": field.choices,
            }
            for field in section.fields
        ],
        "subsections": [
            {
                "path": subsection.path,
                "description": subsection.description,
                "is_array": subsection.is_array,
            }
            for subsection in section.subsections
        ],
    }
