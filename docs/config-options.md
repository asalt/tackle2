# Tackle2 Configuration Reference

This document lists all supported keys in the tackle2 configuration TOML files, including default values and notes on behaviour. You can also explore the schema from the CLI, e.g.:

```
$ tackle2 describe
$ tackle2 describe params.bubbleplot --json
```

## Top-level Structure

```
[params]
[params.barplot]
[params.barplot.advanced]
[params.bubbleplot]
[params.bubbleplot.advanced]
[params.enplot]
[params.heatmap_gsea]
[params.heatmap_gene]
[params.pca]
[params.extra]
[params.genesets]
[params.advanced]
```

## params

- `ranks_from` (string, default `""`): Data source, `"volcano"`, `"gct"`, or `"model"`.
- `savedir` (string, default `./plots`): Root output directory (resolved relative to run cwd).
- `volcanodir` (string): Folder containing volcano TSV files (required when `ranks_from = "volcano"`).
- `rankfiledir` (string, default `"savedir"`): Where generated rank files are written; special value `"savedir"` puts them under `savedir/ranks`.
- `gct_path` (string): Path to GCT file (required when `ranks_from = "gct"`).
- `model_file` (string, default `""`): Path to a TOML file containing a `[model]` table (optional helper when `ranks_from = "model"`).
- `species` (string, default `"Homo sapiens"`): Species name passed to msigdbr.
- `zscore_emat` (bool, default `true`): Z-score expression matrix before analysis.
- `zscore_emat_groupby` (string/bool, default `false`): Grouping column for Z-score normalisation.
- `cut_by` (string, default `"group"`): Metadata column controlling chart faceting.
  
### params.model

- `name` (string, default `"model1"`): Identifier used to label outputs under `model/<type>/<name>/`.
- `type` (string, default `"limma"`): Modelling backend used to generate rank files when `ranks_from = "model"`.
- `design` (string): R formula describing the model matrix (for example `~ 0 + treat`).
- `contrasts` (array of strings, default `[]`): Optional limma contrast expressions; entries may be unnamed (`"treatDrug - treatControl"`) or named (`"MvF = ..."`). When empty, all non-intercept coefficients are exported.
- `volcano_padj_cutoff` (float, default `0.05`): Adjusted P-value threshold used to colour points and compute the `(n / total)` footer in limma volcano plots.
- `volcano_top_n` (int, default `35`): Number of genes automatically labelled on the limma volcano plots, ordered by nominal P-value.

### params.models

An optional array of additional model tables. Each entry accepts the same keys as `params.model` (`name`, `type`, `design`, `contrasts`, `model_file`). When both `params.model` and `params.models` are provided, all models are executed and their ranks concatenated.
Generated limma summaries and volcano PDFs are stored under `savedir/model/<type>/<name>/` inside `tables/`, `volcano_plots/`, and (when expression-derived covariates are used) `metadata/` for the annotated GCT files.

## params.barplot

- `limit` (array of ints, default `[12,20,32]`): `top N` sizes to plot individually/combined.
- `do_individual` (bool, default `true`): Render barplots per comparison.
- `do_combined` (bool, default `true`): Render aggregated barplots.

### params.barplot.advanced

- `stroke_width` (float, default `1.0`): Outline width for bar panels.

## params.bubbleplot

- `limit` (array of ints, default `[12,20,32]`): `top N` sizes for bubble plots.
- `do_individual` (bool, default `true`): Render per-comparison bubbles.
- `do_combined` (bool, default `true`): Render combined bubbles per gene set.
- `glyph` (string, default `"⁕"`): Symbol for `padj < 0.05` markers.

### params.bubbleplot.advanced

- `stroke_alpha` (float, default `0.55`): Outer ring transparency (`0`–`1`).
- `stroke_width` (float, default `0.8`): Outline width for bubble points.

## params.enplot

- `do_individual` (bool, default `true`)
- `do_combined` (bool, default `true`)
- `limit` (int, default `10`)
- `combine_by` (string/bool, default `false`)

## params.heatmap_gsea

- `do` (bool, default `true`)
- `limit` (array of ints, default `[10,20,40,80]`)
- `cut_by` (string): Override for grouping.
- `cluster_rows` (bool, default `true`)
- `cluster_columns` (array of bool, default `[false,true]`)
- `legend_include` (array of strings): Metadata columns to add to legend.

## params.heatmap_gene

- `do` (bool, default `true`)
- `limit` (int, default `10`)

## params.pca (GSEA)

- `do` (bool, default `false`)
- `width` (float, default `7.8`)
- `height` (float, default `7.4`)
- `col_by` (string, default `""`)
- `mark_by` (string, default `""`)
- `top_pc` (int, default `3`)
- `max_pc` (int, default `3`)
- `labSize`, `pointSize`, `sizeLoadingsNames` (floats)

## params.pca_gene

- `do` (bool, default `false`): Enable PCA on the expression matrix.
- `width`/`height` (floats): Figure dimensions for PCA plots.
- `components` (int, default `3`): Number of principal components kept for biplots.
- `metadata_color` (array of strings, default `[]`): One or more metadata columns used to colour samples (one biplot per entry).
- `metadata_shape` (string, default `""`): Metadata column used for point shapes.
- `top_loadings` (int, default `25`): Number of genes rendered in the loadings heatmap.
- `heatmap` (bool, default `true`): Emit a heatmap of top loading genes.
- `cluster_rows` (bool, default `true`): Cluster genes in the PCA loadings heatmap.
- `cluster_columns` (array of bool, default `[false,true]`): Emit heatmaps with column clustering disabled/enabled.
- `cut_by` (string, default `""`): Metadata column used to split the loadings heatmap (leave empty to disable).
- `labSize` / `pointSize` / `sizeLoadingsNames` (floats): Text and point sizes reused for PCA biplots.
- Gene PCA runs once using the parsed GCT matrix; it skips automatically if no GCT was provided.

## params.umap_gene

- `do` (bool, default `false`): Enable UMAP on the expression matrix.
- `width`/`height` (floats): Figure dimensions for UMAP scatter plots.
- `n_neighbors` (int, default `15`): UMAP neighborhood size.
- `min_dist` (float, default `0.1`): Minimum distance between embedded points.
- `metric` (string, default `"euclidean"`): Distance metric passed to `uwot::umap`.
- `seed` (int, default `42`): Random seed for reproducibility.
- `scale` (bool, default `true`): Z-score features per sample before embedding.
- `metadata_color` (array of strings, default `[]`): Metadata columns used for colouring (plots generated per entry).
- `metadata_shape` (string, default `""`): Metadata column used for point shapes.
- `variants` (array of tables, default `[]`): Optional parameter set overrides (each entry may include `name`, `n_neighbors`, `min_dist`, `metadata_color`, etc.).
- `point_type` (string, default `"gene"`): Choose `"gene"` to embed genes or `"sample"` to embed samples.
- `rank_name` (string, default `""`): When `point_type` is `"gene"`, colour genes by the specified rank/comparison (requires matching `.rnk` data).

## params.wordcloud

- `do` (bool, default `false`): Enable per-collection pathway wordcloud summaries from GSEA results.
- `padj_cutoff` (float, default `0.25`): Adjusted p-value threshold; only pathways with `padj <= padj_cutoff` contribute tokens.
- `top_n_pathways` (int, default `50`): Number of pathways (per collection) used to build the wordcloud; sorted by `abs(NES) * -log10(pval)`.
- `max_words` (int, default `70`): Maximum number of distinct words drawn, ordered by aggregated weight.

## params.extra

- `rankname_order` (array of strings): Explicit left-to-right order for rank names (comparisons). When present, every plot that displays comparisons (combined barplots, bubbleplots, GSEA heatmaps, PCA heatmaps, etc.) will honour this sequence after intersecting it with the ranks that were actually generated. A single value of `"sample_order"` is treated as an alias for the `sample_order` vector.
- `samplename_order` (array of strings): Reserved for future use. The field is parsed but currently ignored by the R pipeline.
- `sample_order` (array of strings): Legacy alias for `rankname_order`. If only `sample_order` is provided it is copied to `rankname_order`; the strings must match the comparison names exactly to take effect.

### Rank file helpers

- `rankfiledir/names.txt` (optional): Plain text mapping file used when reading existing `.rnk` files. Each line should take the form `NewLabel=old_file_name.rnk`. Mappings are applied in file order, so you can both rename and reorder comparisons without regenerating ranks. Lines starting with `#` are ignored.

## params.genesets (array of tables)

Each entry defines a gene-set collection:

- `category` (string)
- `subcategory` (string)
- `collapse` (bool)

## params.advanced

- `parallel` (bool, default `true`)
- `cache` (bool, default `true`)
- `replace` (bool, default `false`)
- `cachedir` (string, default `"savedir"` → `savedir/cache`)
- `logfile` (string, default `${savedir}/run.log`)
- `loglevel` (string, default `"INFO"`)
- `pivot_gsea_results` (bool, default `false`)
- `quiet` (bool, default `false`)

## R Options (advanced)

These are R session options that influence naming and presentation. Set them with `options(...)` in your R session before running plots.

- `tackle2_enrichplot_title_wrap` (int, default `40`)
  - Controls the word‑wrap width for enrichplot titles after underscores are replaced with spaces.

- `tackle2_name_map_strip_stems` (bool, default `TRUE`)
  - Enables conservative removal of common lowercase token stems (e.g., `cell` in `cellOCI3`, `treat` in `treatControl`) across a set of labels when generating shorter filenames. Disable if you prefer to keep stems.
