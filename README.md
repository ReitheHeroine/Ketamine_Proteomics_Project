# Ketamine Proteomics Analysis Pipeline

_MS thesis research, Blanco-Suárez Lab, San Diego State University_

End-to-end reproducible workflow for differential protein abundance and pathway enrichment analysis of ketamine-treated cortical astrocytes (mouse, 24 h post-dose, TMT-labeled DDA LC-MS/MS). Built from Proteome Discoverer exports through publication-quality figures, with a validation framework.

## Highlights

- **Reproducible pipeline**: Snakemake DAG with all parameters, paths, and thresholds externalized to `config.yaml`. Overridable via `--config` at the command line.
- **Statistical rigor under constraints**: When per-sample abundances were unavailable (blocking limma / DESeq2 / permutation testing), I designed a three-tier validation framework combining internal consistency checks, tool cross-validation, and literature/database concordance.
- **Handling of small-n statistics**: empirical characterization of n=3 variance behavior (variability vs. fold-change vs. p-value), documenting why the t-test is valid-but-underpowered with df=2 and how that shapes interpretation.
- **Domain-aware data handling**: separate analytical track for presence/absence proteins (routed around Proteome Discoverer's placeholder ratios of 100 / 0.01 to avoid distorting downstream statistics), cell-type fidelity QC against canonical astrocyte markers, and a documented contamination assessment (oligodendrocyte/neuronal) using published cell-type expression data.
- **Publication-grade visualization**: volcano, MA, variability-diagnostic, GO enrichment lollipop, gene-term network, pathway bar/dot, and composite figures, exported as both static PDF and interactive HTML (plotly + kaleido).
- **Carefully done enrichment**: over-representation analysis via g:Profiler against GO:BP and KEGG with a *custom background* of detected proteins (rather than whole genome) and optional REVIGO integration for GO term redundancy reduction.

## Stack

| Layer | Tools |
|---|---|
| Workflow | Snakemake |
| Analysis | Python 3 (pandas, numpy, scipy) |
| Enrichment | g:Profiler (`gprofiler-python`), REVIGO |
| Visualization | plotly, kaleido, matplotlib |
| Config / reproducibility | YAML, project-relative paths, timestamped logs, conda environment |
| Platform | macOS (Apple Silicon) |

## Project overview

Astrocytes are increasingly implicated in synaptic plasticity and post-stroke circuit remodeling. This project characterizes the cortical astrocyte proteomic response to systemic ketamine (10 mg/kg IP) vs. saline control 24 hours post-treatment, asking which biological processes and pathways are reshaped in the ACSA-2-isolated astrocyte proteome. The long-term motivation is to surface candidate therapeutic targets for stroke and related neurological conditions.

## Methods

- **Input**: Proteome Discoverer 3.1 protein-level exports, stratified by MS2 identification coverage (both conditions / ketamine-only / control-only); label-free quantification via the Minora algorithm with match-between-runs.
- **Differential abundance**: Welch's t-test on log2-transformed ketamine/control ratios; technical replicates (LC-MS reinjections) excluded from effective sample size (n=3 biological per group).
- **Presence/absence proteins**: analyzed on a parallel track and combined only at the gene-set-construction stage for enrichment, so Proteome Discoverer's placeholder ratios do not distort quantitative statistics.
- **Multiple testing**: protein-level p-values reported as exploratory; FDR correction applied at the pathway level (g:Profiler g_SCS).
- **Enrichment**: over-representation analysis via g:Profiler against GO:BP and KEGG, using all detected proteins as a custom statistical background. Optional REVIGO pass for GO:BP semantic clustering (Supek et al., 2011, *PLoS ONE*).
- **QC and validation**: canonical astrocyte marker panel check for cell-type fidelity; contamination assessment for myelin and neurofilament proteins (myelin-removal step was skipped during isolation); internal statistical-consistency diagnostics validating expected n=3 behavior.

## Repository layout

```
.
├── Snakefile                         # Workflow definition (DAG)
├── config.yaml                       # Parameters, thresholds, I/O paths
├── scripts/
│   ├── diff_abundance_analysis.py    # Welch's t-test + presence/absence handling
│   ├── pathway_analysis.py           # g:Profiler ORA + REVIGO integration
│   ├── visualize_diff_abundance.py   # Volcano, MA, variability diagnostics
│   ├── plot_abundance_comparison.py  # Condition-wise abundance comparisons
│   ├── plot_cell_type_fidelity.py    # Astrocyte marker QC
│   ├── plot_top_significant_proteins.py
│   └── protein_summary.py
└── results/                          # Outputs (gitignored)
```

## Running it

```bash
conda activate ketamine_project

snakemake -n                                      # dry run: show the DAG
snakemake --cores 1                               # full pipeline
snakemake --cores 1 --config pval=0.1 fc=2.0      # override thresholds
snakemake --dag | dot -Tpng > workflow_dag.png    # render workflow diagram
```

## Data availability

Raw mass spectrometry data and Proteome Discoverer result files are not included. They will be deposited to a public proteomics repository (e.g., PRIDE) upon thesis publication and manuscript acceptance. The `data/` directory is intentionally gitignored; pipeline code, configuration, and workflow definition are fully public for methods review.

## Publication status

**Methods-only until thesis publication.** Per-protein findings, enrichment results, and quantitative summaries are withheld pending thesis defense and peer-reviewed publication. This README describes the analytical pipeline and does not report biological findings.

## Design decisions worth surfacing

- **Presence/absence vs. quantitative proteins are analyzed separately.** Proteome Discoverer assigns placeholder ratios (100, 0.01) when a protein is identified in only one condition. Treating these as real fold changes would contaminate downstream statistics, so they route through a parallel track and merge only at the gene-set stage.
- **Custom detected-proteome background, not whole genome.** Enrichment runs against proteins actually observable in the experiment, which is more conservative and more faithful than a genome-wide background.
- **Technical replicates excluded from effective n.** LC-MS reinjections assess instrument reproducibility but do not increase biological sample size.
- **Configuration-driven, not code-driven.** Every threshold, path, database choice, and visualization parameter lives in `config.yaml` so parameter sweeps do not require editing analysis code.
- **Honest about the small-n ceiling.** With n=3 (df=2), per-protein variance estimates carry ~58× uncertainty; the pipeline documents this, validates it empirically (near-zero correlation between variability and significance), and anchors biological interpretation at the pathway level rather than overclaiming at the protein level.

## Author

Reina Hastings - MS Bioinformatics and Medical Informatics, San Diego State University
reinahastings13@gmail.com
