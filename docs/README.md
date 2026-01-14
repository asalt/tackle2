# gsea-web


##  Intro

This is is a versatile tool for running (preranked) Gene Set Enrichment Analysis (GSEA) across combinations of ranks and gene sets.
It mainly consists of a series of wrapper functions to perform the analysis and sift through the results.


clone the repository:
```bash
git clone https://github.com/asalt/tackle2
```


## Features
R Analysis: The main analysis is performed using R, with fgsea for fast gene set enrichment.
Python Integration: Python manages the database, hosts APIs, and stores results.
  This component is still under development.
Testing: Comprehensive testing using pytest for Python and testthat for R.
Visualization: Utilizes ggplot2 for high-quality plots, including faceted barplots, heatmaps, PCA plots, and bubble plots summarizing pathway directionality.

## Output Naming & Shortening

Outputs now use length‑safe, informative filenames. Two strategies keep names concise while preserving the important bits:

- Common prefix/suffix stripping: For a set of comparisons within a collection, shared leading or trailing text is removed (on token boundaries like `_`, `-`, `.`, or spaces). This preserves the unique part of each label. When combined with hashing, filenames stay portable across filesystems.
- Optional token‑stem stripping: If every label contains lowercase stems at the start of tokens (e.g., `cellOCI3`, `treatControl`), those stems (like `cell`, `treat`) are dropped, keeping `OCI3`, `Control`, etc. This is conservative and only applies when the pattern is present across the whole set.

Pathway prefixes (e.g., `HALLMARK_`, `KEGG_`, `REACTOME_`, `GOBP_`, `GOMF_`, `GOCC_`, `MEDICUS_`) are removed from directory and filename components when the collection already encodes this context in the folder structure.

Title readability: Enrichment plot titles replace underscores with spaces and are wrapped to a configurable width.

Advanced R options:

- `options(tackle2_enrichplot_title_wrap = 40)` controls the wrap width for enrichplot titles (defaults to `40`).
- `options(tackle2_name_map_strip_stems = TRUE)` toggles the optional token‑stem stripping used in filename shortening (defaults to `TRUE`).

These options are read at runtime by the R plotting helpers; set them at the top of your session or within your analysis script before plotting.

## Ordering Samples & Comparisons

The order in which comparisons appear determines how combined plots, heatmaps, and PCA overlays read. You can control this in two complementary ways:

1. **`params.extra.rankname_order`** – Add a vector of comparison names (`rankname`s) to your TOML configuration. Every plotting helper intersects this list with the rank names present in the run, preserving your requested left‑to‑right order.

    ```toml
    [params.extra]
    rankname_order = ["Treated_vs_Control", "Knockout_vs_Control", "Rescue_vs_Knockout"]
    ```

    - `sample_order` is a legacy alias; if you only supply `sample_order`, it is copied into `rankname_order` during parameter sanitisation.
    - The strings must match the comparison names exactly (case sensitive). Missing entries are dropped silently, so double-check spelling if your plots fall back to alphabetical order.

2. **`names.txt`** – When you point the run at an existing rank directory, you can rename (and implicitly re-order) comparisons without touching the `.rnk` files. Drop a `names.txt` file next to the ranks with one mapping per line:

    ```text
    TreatmentA=Treated_vs_Control.rnk
    TreatmentB=Knockout_vs_Control.rnk
    Rescue=Rescue_vs_Knockout.rnk
    ```

    - The left side becomes the new comparison label; the right side references the filename (with optional `.rnk` suffix).
    - Mappings are applied from top to bottom, so you can curate both naming and ordering in a single pass. Lines starting with `#` are treated as comments.

When a GCT file is supplied, column metadata (`gct@cdesc`) is leveraged for group annotations. If `rankname_order` lines up with the `group` column, those factors are re-leveled to match your configuration; otherwise, the pipeline falls back to the natural order of the rank files. For volcano-only runs, you can still steer ordering with `rankname_order` or `names.txt`, even though no additional metadata is attached.

Script main entry points are driven by python command line parsing and variable dispatch.

```
## tackle2 example help

~/t/tackle2> tackle2 run --help
Usage: tackle2 run [OPTIONS]

Options:
  -c, --config FILE               .toml file with additional parameters for
                                  report
  -i, --interactive               run in interactive session within python
                                  rpy2
  -v, --verbose                   verbose output
  -l, --log_level [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  log level
  --help                          Show this message and exit.

```


command line interface is here:

```
./python/cli.py

```


`R`:
  - data loading
  - analysis
  - results routing


`python` :
  - database management
  - API hosting
  - results storage in db

### Python geneset API

The `tackle2` python package now exposes the existing R `geneset_utils` helpers via `rpy2`. Import `get_geneset_collections` to pull msigdbr data frames directly from python without invoking the CLI:

```
from tackle2 import get_geneset_collections, geneset_membership

hallmark = get_geneset_collections([
    {"category": "H", "subcategory": "", "collapse": False},
])
gene_lists = geneset_membership(hallmark["H_"])
```

Collections are returned as `pandas.DataFrame` objects, and `geneset_membership` collapses any collection into `{geneset_name: [genes...]}` dictionaries using the same R implementation the pipeline relies on.


test everything with test.sh (no .bat file yet)
python testing with `pytest ./python/tests/`
R tests located in `./R/tests/` using `testthat`, mostly
test coverage << 100%

## Modules

There are a variety of modules to facilitate data loading, geneset retrievement, formatting/wrangling, and more.
The following are descriptions of the available modules loosly in order for running a new analysis

### io.R

Currently there are two options to input data.
First is a tsv file containing statistics regarding a comparison - e.g. t test between two groups.

The function that handles this datasource is:

```
create_rnkfiles_from_volcano
```

After the identifier and value are properly loaded, a rankfile is saved locally for future use.
Identifier currently is expected to be Entrez GeneID.
Value currently expected to be a zero centered quantitative value about the comparison.
Good choices of value are signedlogp value or logFoldChange

Second, for the purpose of ssGSEA, takes gct file format as input.
Each sample (cid) is separated into a separate rank file.
An option for scaling by zscore is available if data are not zero centered.
Generally expected for data to be approximately normally distributed.
Note: setting the GSEA paramater `p` to `0` should allow for proper analysis of non-normalized data, though explicit control over this parameter is not yet available.

## geneset_utils.R

Facitiates fetching of publically available genesets using `msigdbr`.
Caches fetched genesets to a local file.
Support for multiple species is available through `msigdbr`.

### fgsea.R

Execution performed by `run_one` wrapper around `fgsea::fgsea`.
Also has an option to "collapse" pathways to reduce redundancy.
This is performed with `fgsea::collapsePathways`.
No filtering is performed at this stage; if `collapse == TRUE` redundant pathways will be calculated and indicated.

### plot.R


### Custom Colormap (Annotations)

You can override default annotation colors with a JSON colormap file (preferred). Two JSON shapes are supported:

- Object with `global` and `by_column` maps:

```
{
  "global": {
    "Responder": "#4DAF4A",
    "NonResponder": "#984EA3"
  },
  "by_column": {
    "group": {
      "Treated": "#E41A1C",
      "Control": "#377EB8"
    },
    "batch": {
      "BatchA": "#1B9E77",
      "BatchB": "#D95F02"
    }
  }
}
```

- Array of entries with optional `column` field:

```
[
  {"column": "group", "name": "Treated", "color": "#E41A1C"},
  {"column": "group", "name": "Control", "color": "#377EB8"},
  {"name": "Responder", "color": "#4DAF4A"},
  {"name": "NonResponder", "color": "#984EA3"}
]
```

Rules and tips:
- Colors should be hex codes like `#1F77B4`; named R colors are also accepted and converted to hex.
- Matching is exact and case sensitive on `column` and `name`.
- Only discrete annotations use explicit mappings. Continuous annotations (numeric/decimal) use gradients automatically.

You can copy the example files the CLI ships with:

```
tackle2 get-config --include-colormap
```

This writes `tackle2.toml` and `colormap.example.json` to the current directory. Enable the mapping in your config (path is relative to run root):

```
[params.extra]
colormap_file = "config/colormap.json"
```

Backward compatibility: If you point to a `.tsv`/`.csv` file, a two‑column (name,color) or three‑column (column,name,color) table is still supported.



### plot_bubble.R

Bubble plots provide a compact view of pathway directionality for each comparison. Key behaviour:

- Fill colour encodes signed enrichment strength using `1 - pval` (reds for positive NES, blues for negative NES).
- Circle outline indicates significance thresholds (`padj < 0.25`); a centered asterisk marks pathways with `padj < 0.05`.
- Subtitles show the associated rank / contrast along with the `top N` limit used in the selection.
- Output files are written alongside barplots with consistent, length-safe filenames so they remain portable across operating systems.

Enable or disable bubble plotting from configuration by adding a `[params.bubbleplot]` section, e.g.

```
[params.bubbleplot]
limit = [12, 20, 32]
do_individual = true
do_combined = true
```

The limit vector controls how many pathways are shown per figure; the `do_*` toggles mirror the barplot configuration.


## Citations
To cite package ‘fgsea’ in publications use:

  G. Korotkevich, V. Sukhov, A. Sergushichev. Fast gene set enrichment analysis. bioRxiv (2019), doi:10.1101/060012

```bibtex
  @Article{,
    author = {Gennady Korotkevich and Vladimir Sukhov and Alexey Sergushichev},
    title = {Fast gene set enrichment analysis},
    year = {2019},
    doi = {10.1101/060012},
    publisher = {Cold Spring Harbor Labs Journals},
    url = {http://biorxiv.org/content/early/2016/06/20/060012},
    journal = {bioRxiv},
  }
```


To cite package ‘msigdbr’ in publications use:

  Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.

```bibtex

  @Manual{,
    title = {msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format},
    author = {Igor Dolgalev},
    year = {2022},
    note = {R package version 7.5.1},
    url = {https://CRAN.R-project.org/package=msigdbr},
  }
```

To cite ggplot2 in publications, please use

  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.


```bibtex
  @Book{,
    author = {Hadley Wickham},
    title = {ggplot2: Elegant Graphics for Data Analysis},
    publisher = {Springer-Verlag New York},
    year = {2016},
    isbn = {978-3-319-24277-4},
    url = {https://ggplot2.tidyverse.org},
  }
```

To cite package ‘purrr’ in publications use:

  Wickham H, Henry L (2023). _purrr: Functional Programming Tools_. R package version 1.0.2, <https://CRAN.R-project.org/package=purrr>.

```bibtex

@Manual{,
    title = {purrr: Functional Programming Tools},
    author = {Hadley Wickham and Lionel Henry},
    year = {2023},
    note = {R package version 1.0.2},
    url = {https://CRAN.R-project.org/package=purrr},
  }

```

To cite package ‘ComplexHeatmap’ in publications use:

```bibtex

  @Article{,
    title = {Complex heatmaps reveal patterns and correlations in multidimensional genomic data},
    author = {Zuguang Gu and Roland Eils and Matthias Schlesner},
    journal = {Bioinformatics},
    doi = {10.1093/bioinformatics/btw313},
    year = {2016},
  }
```

```bibtex
  @Article{,
    title = {Complex Heatmap Visualization},
    author = {Zuguang Gu},
    doi = {10.1002/imt2.43},
    journal = {iMeta},
    year = {2022},
  }


```
