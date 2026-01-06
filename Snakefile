# title: Snakefile
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-09
# last modified: 2025-12-09
#
# purpose:
#   Snakemake workflow for ketamine proteomics differential abundance analysis.
#   Automates the full pipeline from raw Proteome Discoverer output to
#   publication-ready figures and tables.
#
# usage:
#   # Dry run (shows what would be executed without running)
#   snakemake -n
#
#   # Run the full pipeline
#   snakemake --cores 1
#
#   # Run with a specific p-value threshold
#   snakemake --cores 1 --config pval=0.1
#
#   # Generate a workflow diagram
#   snakemake --dag | dot -Tpng > workflow_dag.png
#
#   # Clean all outputs and start fresh
#   snakemake clean
#
# notes:
#   - Requires: Python 3, pandas, numpy, plotly, kaleido
#   - Input files should be in data/ directory
#   - All outputs go to results/ directory
#   - Logs are written to logs/ directory

# =============================================================================
# CONFIGURATION
# =============================================================================

# Default configuration values (can be overridden via --config on command line)
# Example: snakemake --cores 1 --config pval=0.1 fc=2.0
configfile: 'config.yaml'

# Set defaults if not specified in config
PVAL_THRESHOLD = config.get('pval', 0.05)
FC_THRESHOLD = config.get('fc', 1.5)

# -----------------------------------------------------------------------------
# Directory structure
# -----------------------------------------------------------------------------
DATA_DIR = 'data'
RESULTS_DIR = 'results'
SCRIPTS_DIR = 'scripts'
LOGS_DIR = 'logs'
FIGURES_DIR = f'{RESULTS_DIR}/figures'

# -----------------------------------------------------------------------------
# Input files (from Proteome Discoverer)
# -----------------------------------------------------------------------------
INPUT_FILES = [
    f'{DATA_DIR}/proteins_ketamine.csv',
    f'{DATA_DIR}/proteins_control.csv',
    f'{DATA_DIR}/proteins_control_ketamine.csv'
]

# =============================================================================
# TARGET RULE (defines the final outputs we want)
# =============================================================================

rule all:
    '''
    Default target rule - defines all final outputs.
    Snakemake will work backward from these to determine what needs to run.
    '''
    input:
        # Differential abundance results
        f'{RESULTS_DIR}/combined/all_proteins_categorized.csv',
        f'{RESULTS_DIR}/quantitative/significant.csv',
        f'{RESULTS_DIR}/presence_absence/ketamine_specific.csv',
        f'{RESULTS_DIR}/summary_report.txt',
        # Visualization outputs
        f'{FIGURES_DIR}/volcano_plot.pdf',
        f'{FIGURES_DIR}/volcano_plot.html',
        f'{FIGURES_DIR}/ma_plot.pdf',
        f'{FIGURES_DIR}/ma_plot.html',
        f'{FIGURES_DIR}/summary_bar_chart.pdf',
        f'{FIGURES_DIR}/top_proteins_bar_chart.pdf',
        f'{FIGURES_DIR}/top20_upregulated_table.pdf',
        f'{FIGURES_DIR}/top20_downregulated_table.pdf'


# =============================================================================
# ANALYSIS RULES
# =============================================================================

rule differential_abundance:
    '''
    Run differential abundance analysis on proteomics data.
    Categorizes proteins and applies significance thresholds.
    '''
    input:
        proteins_ket = f'{DATA_DIR}/proteins_ketamine.csv',
        proteins_ctrl = f'{DATA_DIR}/proteins_control.csv',
        proteins_both = f'{DATA_DIR}/proteins_control_ketamine.csv'
    output:
        combined = f'{RESULTS_DIR}/combined/all_proteins_categorized.csv',
        significant = f'{RESULTS_DIR}/quantitative/significant.csv',
        sig_up = f'{RESULTS_DIR}/quantitative/significant_up.csv',
        sig_down = f'{RESULTS_DIR}/quantitative/significant_down.csv',
        ket_specific = f'{RESULTS_DIR}/presence_absence/ketamine_specific.csv',
        ctrl_specific = f'{RESULTS_DIR}/presence_absence/control_specific.csv',
        summary = f'{RESULTS_DIR}/summary_report.txt'
    params:
        pval = PVAL_THRESHOLD,
        fc = FC_THRESHOLD
    log:
        f'{LOGS_DIR}/differential_abundance.log'
    shell:
        '''
        python {SCRIPTS_DIR}/diff_abundance_analysis.py \
            --input_dir {DATA_DIR} \
            --output_dir {RESULTS_DIR} \
            --log_dir {LOGS_DIR} \
            --pval_threshold {params.pval} \
            --fc_threshold {params.fc} \
            2>&1 | tee {log}
        '''


rule visualizations:
    '''
    Generate all visualization outputs (volcano plot, MA plot, tables, etc.)
    '''
    input:
        data = f'{RESULTS_DIR}/combined/all_proteins_categorized.csv'
    output:
        volcano_pdf = f'{FIGURES_DIR}/volcano_plot.pdf',
        volcano_html = f'{FIGURES_DIR}/volcano_plot.html',
        ma_pdf = f'{FIGURES_DIR}/ma_plot.pdf',
        ma_html = f'{FIGURES_DIR}/ma_plot.html',
        summary_bar = f'{FIGURES_DIR}/summary_bar_chart.pdf',
        top_bar_pdf = f'{FIGURES_DIR}/top_proteins_bar_chart.pdf',
        top_bar_html = f'{FIGURES_DIR}/top_proteins_bar_chart.html',
        up_table_pdf = f'{FIGURES_DIR}/top20_upregulated_table.pdf',
        up_table_html = f'{FIGURES_DIR}/top20_upregulated_table.html',
        down_table_pdf = f'{FIGURES_DIR}/top20_downregulated_table.pdf',
        down_table_html = f'{FIGURES_DIR}/top20_downregulated_table.html'
    params:
        pval = PVAL_THRESHOLD
    log:
        f'{LOGS_DIR}/visualizations.log'
    shell:
        '''
        python {SCRIPTS_DIR}/visualize_diff_abundance.py \
            --input {input.data} \
            --output_dir {FIGURES_DIR} \
            --log_dir {LOGS_DIR} \
            --pval_threshold {params.pval} \
            2>&1 | tee {log}
        '''


# =============================================================================
# UTILITY RULES
# =============================================================================

rule clean:
    '''
    Remove all generated outputs to start fresh.
    Usage: snakemake clean
    '''
    shell:
        '''
        rm -rf {RESULTS_DIR}/quantitative/*
        rm -rf {RESULTS_DIR}/presence_absence/*
        rm -rf {RESULTS_DIR}/combined/*
        rm -rf {FIGURES_DIR}/*
        rm -f {RESULTS_DIR}/summary_report.txt
        rm -f {LOGS_DIR}/*.log
        echo "Cleaned all outputs!"
        '''


rule create_dirs:
    '''
    Create the directory structure for the project.
    Usage: snakemake create_dirs
    '''
    shell:
        '''
        mkdir -p {DATA_DIR}
        mkdir -p {RESULTS_DIR}/quantitative
        mkdir -p {RESULTS_DIR}/presence_absence
        mkdir -p {RESULTS_DIR}/combined
        mkdir -p {FIGURES_DIR}
        mkdir -p {SCRIPTS_DIR}
        mkdir -p {LOGS_DIR}
        echo "Directory structure created!"
        '''


rule dag:
    '''
    Generate a visual diagram of the workflow.
    Usage: snakemake dag
    Requires: graphviz (install with: conda install graphviz)
    '''
    output:
        'workflow_dag.png'
    shell:
        '''
        snakemake --dag | dot -Tpng > {output}
        echo "Workflow diagram saved to {output}"
        '''