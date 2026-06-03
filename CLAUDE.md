# Ketamine Proteomics Project

## Project context
Proteomic analysis of ketamine-treated cortical astrocytes and their role in neural plasticity. PI: Dr. Elena Blanco-Suarez (Blanco-Suarez Lab, SDSU). This is thesis work for an MS in Bioinformatics and Medical Informatics.

## Stack
- Primary: R (tidyverse, limma, DEP, clusterProfiler, ggplot2)
- Supporting: Python (data wrangling, automation), Bash (pipeline glue)
- Documents: RMarkdown for analysis notebooks, docx for project documentation
- Pipeline: Nextflow under evaluation for reproducibility

## Data
- Mass spectrometry proteomics (TMT-labeled, DDA acquisition)
- Experimental groups: ketamine-treated vs control cortical astrocytes
- Downstream: differential expression, pathway enrichment, GO analysis

## File organization conventions
- Raw data: never modified, stored separately from processed
- Analysis scripts: one script per major analysis step
- Outputs: figures go in results/figures/, tables in results/tables/

## R-specific conventions
- Use project-relative paths (here::here() or explicit relative paths), never absolute
- RMarkdown YAML header: include title, author, date, output format
- Chunk naming: descriptive names (e.g., ```{r load-proteomics-data})
- Session info: always print sessionInfo() at the end of RMarkdown documents

## When making changes
- Summarize what you changed and why in a comment or commit message
- If modifying an existing analysis, note what changed in the progress log
- Flag any assumptions about data structure or experimental design for my review
