# title: diff_abundance_analysis.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-09
# last modified: 2025-12-09
#
# purpose:
#   Perform differential abundance analysis on proteomics data from Proteome
#   Discoverer 3.1 output. Categorizes proteins into quantitative (suitable for
#   differential abundance analysis) and presence/absence (detected in only one
#   condition) groups. Applies user-specified significance threshold and generates
#   organized output files for downstream pathway analysis.
#
# inputs:
#   - proteins_ketamine.csv: Proteins with high confidence in ketamine only
#   - proteins_control.csv: Proteins with high confidence in control only
#   - proteins_control_ketamine.csv: Proteins with high confidence in both
#
# outputs:
#   ketamine_project/
#   ├── data/
#   │   └── [input CSV files]
#   ├── scripts/
#   │   └── diff_abundance_analysis.py
#   ├── results/
#   │   ├── quantitative/
#   │   │   ├── all_quantitative_proteins.csv
#   │   │   ├── significant.csv
#   │   │   ├── significant_up.csv
#   │   │   └── significant_down.csv
#   │   ├── presence_absence/
#   │   │   ├── all_presence_absence_proteins.csv
#   │   │   ├── ketamine_specific.csv
#   │   │   └── control_specific.csv
#   │   ├── combined/
#   │   │   └── all_proteins_categorized.csv
#   │   └── summary_report.txt
#   └── logs/
#       └── diff_abundance_YYYYMMDD_HHMMSS.log
#
# usage:
#   python diff_abundance_analysis.py \
#       --input_dir ../data \
#       --output_dir ../results \
#       --log_dir ../logs \
#       --pval_threshold 0.05 \
#       --fc_threshold 1.5
#    
#   copy/paste: python diff_abundance_analysis.py --input_dir ../data --output_dir ../results --log_dir ../logs --pval_threshold 0.05 --fc_threshold 1.5
#
# notes:
#   - Proteins with ratio = 100 or 0.01 are classified as presence/absence
#   - User specifies p-value threshold (common values: 0.05 or 0.1)
#   - Fold-change threshold is applied as an informational flag, not a hard filter

import pandas as pd
import numpy as np
import argparse
import os
import sys
from datetime import datetime

# =============================================================================
# CONFIGURATION AND CONSTANTS
# =============================================================================

# Column names from Proteome Discoverer output
RATIO_COL = 'Abundance Ratio: (ketamine) / (control)'
PVAL_COL = 'Abundance Ratio Adj. P-Value: (ketamine) / (control)'
VARIABILITY_COL = 'Abundance Ratio Variability [%]: (ketamine) / (control)'
ACCESSION_COL = 'Accession'
GENE_COL = 'Gene Symbol'
DESCRIPTION_COL = 'Description'

# Presence/absence ratio values assigned by Proteome Discoverer
PRESENCE_RATIO_HIGH = 100.0
PRESENCE_RATIO_LOW = 0.01

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def setup_logging(log_dir):
    '''
    Initialize the log file with header information.
    
    Parameters:
        log_dir (str): Directory where log file will be created
        
    Returns:
        str: Path to the log file
    '''
    # ---------------------------------------------------------------------
    # Create timestamped log filename
    # ---------------------------------------------------------------------
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f'diff_abundance_{timestamp}.log'
    log_path = os.path.join(log_dir, log_filename)
    
    # ---------------------------------------------------------------------
    # Write log header
    # ---------------------------------------------------------------------
    with open(log_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('KETAMINE PROTEOMICS ANALYSIS PROJECT - ANALYSIS LOG\n')
        f.write('=' * 70 + '\n\n')
        f.write(f'Log file: {log_filename}\n')
        f.write(f'Analysis started: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
    
    return log_path


def log_message(log_path, message, print_to_console=True):
    '''
    Write a message to the log file and optionally print to console.
    
    Parameters:
        log_path (str): Path to the log file
        message (str): Message to write
        print_to_console (bool): Whether to also print to stdout
    '''
    timestamp = datetime.now().strftime('%H:%M:%S')
    formatted_msg = f'[{timestamp}] {message}'
    
    with open(log_path, 'a') as f:
        f.write(formatted_msg + '\n')
    
    if print_to_console:
        print(formatted_msg)


def load_and_merge_data(input_dir, log_path):
    '''
    Load the three input CSV files and merge them into a single dataframe.
    
    Parameters:
        input_dir (str): Directory containing input files
        log_path (str): Path to log file
        
    Returns:
        pd.DataFrame: Merged dataframe with source file annotation
    '''
    log_message(log_path, 'Loading input files...')
    
    # ---------------------------------------------------------------------
    # Load each file and add source annotation
    # ---------------------------------------------------------------------
    files = {
        'ketamine': 'proteins_ketamine.csv',
        'control': 'proteins_control.csv',
        'both': 'proteins_control_ketamine.csv'
    }
    
    dataframes = []
    for source, filename in files.items():
        filepath = os.path.join(input_dir, filename)
        
        if not os.path.exists(filepath):
            log_message(log_path, f'  ERROR: File not found: {filepath}')
            sys.exit(1)
        
        df = pd.read_csv(filepath)
        df['source_file'] = source
        dataframes.append(df)
        log_message(log_path, f'  Loaded {filename}: {len(df)} proteins')
    
    # ---------------------------------------------------------------------
    # Merge all dataframes
    # ---------------------------------------------------------------------
    merged = pd.concat(dataframes, ignore_index=True)
    log_message(log_path, f'  Total proteins after merge: {len(merged)}')
    
    return merged


def categorize_proteins(df, log_path):
    '''
    Categorize proteins into quantitative vs presence/absence groups.
    
    Proteins with ratio = 100 or 0.01 are classified as presence/absence.
    All others with valid ratios are classified as quantitative.
    
    Parameters:
        df (pd.DataFrame): Merged protein dataframe
        log_path (str): Path to log file
        
    Returns:
        pd.DataFrame: Dataframe with category annotations
    '''
    log_message(log_path, 'Categorizing proteins...')
    
    # ---------------------------------------------------------------------
    # Create category column based on ratio values
    # ---------------------------------------------------------------------
    def assign_category(ratio):
        if pd.isna(ratio):
            return 'missing_ratio'
        elif ratio == PRESENCE_RATIO_HIGH:
            return 'presence_absence_ketamine_specific'
        elif ratio == PRESENCE_RATIO_LOW:
            return 'presence_absence_control_specific'
        else:
            return 'quantitative'
    
    df['category'] = df[RATIO_COL].apply(assign_category)
    
    # ---------------------------------------------------------------------
    # Log category counts
    # ---------------------------------------------------------------------
    category_counts = df['category'].value_counts()
    log_message(log_path, '  Category breakdown:')
    for cat, count in category_counts.items():
        log_message(log_path, f'    {cat}: {count}')
    
    return df


def calculate_derived_metrics(df, pval_threshold, fc_threshold, log_path):
    '''
    Calculate log2 fold change and apply significance/FC flags.
    
    Parameters:
        df (pd.DataFrame): Protein dataframe with categories
        pval_threshold (float): P-value threshold for significance
        fc_threshold (float): Fold-change threshold for flagging
        log_path (str): Path to log file
        
    Returns:
        pd.DataFrame: Dataframe with derived metrics
    '''
    log_message(log_path, 'Calculating derived metrics...')
    
    # ---------------------------------------------------------------------
    # Calculate log2 fold change
    # ---------------------------------------------------------------------
    df['log2_fold_change'] = np.log2(df[RATIO_COL])
    
    # ---------------------------------------------------------------------
    # Apply significance flag at user-specified threshold
    # ---------------------------------------------------------------------
    df['significant'] = (df[PVAL_COL] <= pval_threshold) & (df['category'] == 'quantitative')
    
    # ---------------------------------------------------------------------
    # Apply fold-change flag
    # ---------------------------------------------------------------------
    log2_fc_threshold = np.log2(fc_threshold)
    df['meets_fc_threshold'] = np.abs(df['log2_fold_change']) >= log2_fc_threshold
    
    # ---------------------------------------------------------------------
    # Determine direction of change
    # ---------------------------------------------------------------------
    def assign_direction(row):
        if row['category'] == 'presence_absence_ketamine_specific':
            return 'ketamine_specific'
        elif row['category'] == 'presence_absence_control_specific':
            return 'control_specific'
        elif row['category'] == 'quantitative':
            if row['log2_fold_change'] > 0:
                return 'up_in_ketamine'
            elif row['log2_fold_change'] < 0:
                return 'down_in_ketamine'
            else:
                return 'no_change'
        else:
            return 'undetermined'
    
    df['direction'] = df.apply(assign_direction, axis=1)
    
    log_message(log_path, f'  P-value threshold: {pval_threshold}')
    log_message(log_path, f'  Log2 FC threshold: {log2_fc_threshold:.3f} (FC = {fc_threshold})')
    
    return df


def generate_output_tables(df, pval_threshold, output_dir, log_path):
    '''
    Generate all output CSV files organized by category.
    
    Parameters:
        df (pd.DataFrame): Fully annotated protein dataframe
        pval_threshold (float): P-value threshold used for naming
        output_dir (str): Base output directory
        log_path (str): Path to log file
        
    Returns:
        dict: Dictionary containing counts for summary report
    '''
    log_message(log_path, 'Generating output tables...')
    
    counts = {}
    
    # ---------------------------------------------------------------------
    # Define output columns for clean export
    # ---------------------------------------------------------------------
    core_cols = [
        ACCESSION_COL, GENE_COL, DESCRIPTION_COL,
        RATIO_COL, 'log2_fold_change', PVAL_COL, VARIABILITY_COL,
        'Abundances (Grouped): control', 'Abundances (Grouped): ketamine',
        'category', 'direction', 'significant', 'meets_fc_threshold', 'source_file'
    ]
    
    # Filter to columns that exist in dataframe
    export_cols = [c for c in core_cols if c in df.columns]
    
    # ---------------------------------------------------------------------
    # Create subdirectories
    # ---------------------------------------------------------------------
    quant_dir = os.path.join(output_dir, 'quantitative')
    pa_dir = os.path.join(output_dir, 'presence_absence')
    combined_dir = os.path.join(output_dir, 'combined')
    
    for d in [quant_dir, pa_dir, combined_dir]:
        os.makedirs(d, exist_ok=True)
    
    # ---------------------------------------------------------------------
    # QUANTITATIVE PROTEINS
    # ---------------------------------------------------------------------
    log_message(log_path, '  Generating quantitative protein tables...')
    
    quant_df = df[df['category'] == 'quantitative'].copy()
    counts['total_quantitative'] = len(quant_df)
    
    # All quantitative proteins
    quant_df[export_cols].to_csv(
        os.path.join(quant_dir, 'all_quantitative_proteins.csv'),
        index=False
    )
    
    # Significant proteins
    sig = quant_df[quant_df['significant']].copy()
    sig_sorted = sig.sort_values(PVAL_COL)
    sig_sorted[export_cols].to_csv(
        os.path.join(quant_dir, 'significant.csv'),
        index=False
    )
    counts['sig_total'] = len(sig)
    
    # Significant UP in ketamine
    sig_up = sig[sig['direction'] == 'up_in_ketamine'].sort_values(PVAL_COL)
    sig_up[export_cols].to_csv(
        os.path.join(quant_dir, 'significant_up.csv'),
        index=False
    )
    counts['sig_up'] = len(sig_up)
    
    # Significant DOWN in ketamine
    sig_down = sig[sig['direction'] == 'down_in_ketamine'].sort_values(PVAL_COL)
    sig_down[export_cols].to_csv(
        os.path.join(quant_dir, 'significant_down.csv'),
        index=False
    )
    counts['sig_down'] = len(sig_down)
    
    # Count proteins meeting FC threshold
    counts['sig_meets_fc'] = sig['meets_fc_threshold'].sum()
    
    # ---------------------------------------------------------------------
    # PRESENCE/ABSENCE PROTEINS
    # ---------------------------------------------------------------------
    log_message(log_path, '  Generating presence/absence protein tables...')
    
    pa_df = df[df['category'].str.startswith('presence_absence')].copy()
    counts['total_presence_absence'] = len(pa_df)
    
    # All presence/absence proteins
    pa_df[export_cols].to_csv(
        os.path.join(pa_dir, 'all_presence_absence_proteins.csv'),
        index=False
    )
    
    # Ketamine-specific
    ket_specific = pa_df[pa_df['category'] == 'presence_absence_ketamine_specific'].copy()
    ket_specific[export_cols].to_csv(
        os.path.join(pa_dir, 'ketamine_specific.csv'),
        index=False
    )
    counts['ketamine_specific'] = len(ket_specific)
    
    # Control-specific
    ctrl_specific = pa_df[pa_df['category'] == 'presence_absence_control_specific'].copy()
    ctrl_specific[export_cols].to_csv(
        os.path.join(pa_dir, 'control_specific.csv'),
        index=False
    )
    counts['control_specific'] = len(ctrl_specific)
    
    # ---------------------------------------------------------------------
    # COMBINED MASTER FILE
    # ---------------------------------------------------------------------
    log_message(log_path, '  Generating combined master file...')
    
    df[export_cols].to_csv(
        os.path.join(combined_dir, 'all_proteins_categorized.csv'),
        index=False
    )
    counts['total_proteins'] = len(df)
    counts['missing_ratio'] = len(df[df['category'] == 'missing_ratio'])
    
    return counts


def generate_summary_report(counts, pval_threshold, fc_threshold, output_dir, log_path):
    '''
    Generate a human-readable summary report of the analysis.
    
    Parameters:
        counts (dict): Dictionary of counts from output generation
        pval_threshold (float): P-value threshold used
        fc_threshold (float): Fold-change threshold used
        output_dir (str): Output directory
        log_path (str): Path to log file
    '''
    log_message(log_path, 'Generating summary report...')
    
    report_path = os.path.join(output_dir, 'summary_report.txt')
    
    with open(report_path, 'w') as f:
        # ---------------------------------------------------------------------
        # Header
        # ---------------------------------------------------------------------
        f.write('=' * 70 + '\n')
        f.write('KETAMINE PROTEOMICS ANALYSIS PROJECT\n')
        f.write('DIFFERENTIAL ABUNDANCE ANALYSIS - SUMMARY REPORT\n')
        f.write('=' * 70 + '\n\n')
        f.write(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        
        # ---------------------------------------------------------------------
        # Analysis Parameters
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('ANALYSIS PARAMETERS\n')
        f.write('-' * 70 + '\n')
        f.write(f'P-value threshold:            {pval_threshold}\n')
        f.write(f'Fold-change threshold:        {fc_threshold} (|log2FC| >= {np.log2(fc_threshold):.3f})\n')
        f.write(f'Presence/absence high ratio:  {PRESENCE_RATIO_HIGH}\n')
        f.write(f'Presence/absence low ratio:   {PRESENCE_RATIO_LOW}\n\n')
        
        # ---------------------------------------------------------------------
        # Dataset Overview
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('DATASET OVERVIEW\n')
        f.write('-' * 70 + '\n')
        f.write(f'Total proteins analyzed:      {counts["total_proteins"]}\n')
        f.write(f'  Quantitative:               {counts["total_quantitative"]}\n')
        f.write(f'  Presence/Absence:           {counts["total_presence_absence"]}\n')
        f.write(f'  Missing ratio (excluded):   {counts["missing_ratio"]}\n\n')
        
        # ---------------------------------------------------------------------
        # Quantitative Results
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('QUANTITATIVE DIFFERENTIAL ABUNDANCE RESULTS\n')
        f.write('-' * 70 + '\n\n')
        
        f.write(f'Significance threshold (adj. p-value <= {pval_threshold}):\n')
        f.write(f'  Total significant:          {counts["sig_total"]}\n')
        f.write(f'    Upregulated in ketamine:  {counts["sig_up"]}\n')
        f.write(f'    Downregulated in ketamine:{counts["sig_down"]}\n')
        f.write(f'  Meeting FC threshold:       {counts["sig_meets_fc"]} / {counts["sig_total"]}\n\n')
        
        # ---------------------------------------------------------------------
        # Presence/Absence Results
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('PRESENCE/ABSENCE RESULTS\n')
        f.write('-' * 70 + '\n')
        f.write(f'Total presence/absence:       {counts["total_presence_absence"]}\n')
        f.write(f'  Ketamine-specific:          {counts["ketamine_specific"]}\n')
        f.write(f'  Control-specific:           {counts["control_specific"]}\n\n')
        
        # ---------------------------------------------------------------------
        # Output Files
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('OUTPUT FILES GENERATED\n')
        f.write('-' * 70 + '\n')
        f.write('quantitative/\n')
        f.write('  all_quantitative_proteins.csv  - All proteins with real ratios\n')
        f.write('  significant.csv                - Significant proteins\n')
        f.write('  significant_up.csv             - Upregulated significant proteins\n')
        f.write('  significant_down.csv           - Downregulated significant proteins\n\n')
        f.write('presence_absence/\n')
        f.write('  all_presence_absence_proteins.csv - All presence/absence proteins\n')
        f.write('  ketamine_specific.csv          - Detected only in ketamine\n')
        f.write('  control_specific.csv           - Detected only in control\n\n')
        f.write('combined/\n')
        f.write('  all_proteins_categorized.csv   - Master file with all annotations\n\n')
        
        # ---------------------------------------------------------------------
        # Notes
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('NOTES\n')
        f.write('-' * 70 + '\n')
        f.write('- Quantitative proteins include both High/High and High/Peak Found\n')
        f.write('  detection patterns from Proteome Discoverer\n')
        f.write('- Presence/absence proteins have ratio = 100 (ketamine-specific) or\n')
        f.write('  0.01 (control-specific) as assigned by Proteome Discoverer\n')
        f.write('- Fold-change threshold is provided as an informational flag;\n')
        f.write('  it is not used as a hard filter for significance\n')
    
    log_message(log_path, f'  Summary report saved to: {report_path}')


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    '''
    Main entry point for the differential abundance analysis.
    '''
    # ---------------------------------------------------------------------
    # Parse command line arguments
    # ---------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description='Differential abundance analysis for ketamine proteomics data'
    )
    parser.add_argument(
        '--input_dir', '-i',
        type=str,
        required=True,
        help='Directory containing input CSV files from Proteome Discoverer'
    )
    parser.add_argument(
        '--output_dir', '-o',
        type=str,
        required=True,
        help='Directory for output result files (e.g., ketamine_project/results)'
    )
    parser.add_argument(
        '--log_dir', '-l',
        type=str,
        required=True,
        help='Directory for log files (e.g., ketamine_project/logs)'
    )
    parser.add_argument(
        '--pval_threshold', '-p',
        type=float,
        default=0.05,
        help='Adjusted p-value threshold for significance (default: 0.05)'
    )
    parser.add_argument(
        '--fc_threshold', '-fc',
        type=float,
        default=1.5,
        help='Fold-change threshold for flagging (default: 1.5)'
    )
    
    args = parser.parse_args()
    
    # ---------------------------------------------------------------------
    # Create output and log directories, initialize logging
    # ---------------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)
    log_path = setup_logging(args.log_dir)
    
    log_message(log_path, '=' * 50)
    log_message(log_path, 'KETAMINE PROTEOMICS ANALYSIS PROJECT')
    log_message(log_path, 'Differential Abundance Analysis')
    log_message(log_path, '=' * 50)
    log_message(log_path, f'Input directory:  {args.input_dir}')
    log_message(log_path, f'Output directory: {args.output_dir}')
    log_message(log_path, f'Log directory:    {args.log_dir}')
    log_message(log_path, f'P-value threshold: {args.pval_threshold}')
    log_message(log_path, f'FC threshold:     {args.fc_threshold}')
    log_message(log_path, '')
    
    # ---------------------------------------------------------------------
    # Run analysis pipeline
    # ---------------------------------------------------------------------
    
    # Step 1: Load and merge data
    df = load_and_merge_data(args.input_dir, log_path)
    
    # Step 2: Categorize proteins
    df = categorize_proteins(df, log_path)
    
    # Step 3: Calculate derived metrics
    df = calculate_derived_metrics(df, args.pval_threshold, args.fc_threshold, log_path)
    
    # Step 4: Generate output tables
    counts = generate_output_tables(df, args.pval_threshold, args.output_dir, log_path)
    
    # Step 5: Generate summary report
    generate_summary_report(counts, args.pval_threshold, args.fc_threshold, args.output_dir, log_path)
    
    # ---------------------------------------------------------------------
    # Finalize
    # ---------------------------------------------------------------------
    log_message(log_path, '')
    log_message(log_path, '=' * 50)
    log_message(log_path, 'Analysis complete!')
    log_message(log_path, '=' * 50)
    log_message(log_path, f'Analysis finished: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')


if __name__ == '__main__':
    main()