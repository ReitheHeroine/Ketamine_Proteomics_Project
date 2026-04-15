# title: protein_summary.py
# project: Ketamine Astrocyte Proteomics
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-18
# last modified: 2026-03-18
#
# purpose:
#   Generate a comprehensive summary for any given protein (by gene symbol)
#   that includes all proteomics data, pathway/GO term memberships, and
#   cell-type marker status. Queries multiple project data sources:
#     - Proteomics dataset (abundances, fold change, p-values, variability)
#     - Pathway analysis report (GO:BP and KEGG term memberships)
#     - Cell-type marker definitions (marker status and specificity notes)
#
# inputs:
#   - data/all_proteins_categorized.csv
#   - data/cell_type_markers.csv
#   - results/pathway_analysis/pathway_analysis_report.txt
#
# outputs:
#   Prints summary to stdout. Optionally saves to file with --output.
#
# usage examples:
#   python scripts/protein_summary.py Snap25
#   python scripts/protein_summary.py Snap25 Vamp2 Syt1
#   python scripts/protein_summary.py Snap25 --output results/protein_summaries/
#   python scripts/protein_summary.py --all-significant

import pandas as pd
import numpy as np
import argparse
import re
import os
import sys

# ======================================================================
# 1. CONFIGURATION
# ======================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTEINS_FILE = os.path.join(BASE_DIR, 'data', 'all_proteins_categorized.csv')
MARKERS_FILE = os.path.join(BASE_DIR, 'data', 'cell_type_markers.csv')
REPORT_FILE = os.path.join(BASE_DIR, 'results', 'pathway_analysis', 'pathway_analysis_report.txt')

# -- column names --
COL_GENE = 'Gene Symbol'
COL_ACCESSION = 'Accession'
COL_DESC = 'Description'
COL_RATIO = 'Abundance Ratio: (ketamine) / (control)'
COL_LOG2FC = 'log2_fold_change'
COL_PVAL = 'Abundance Ratio Adj. P-Value: (ketamine) / (control)'
COL_VARIABILITY = 'Abundance Ratio Variability [%]: (ketamine) / (control)'
COL_CONTROL = 'Abundances (Grouped): control'
COL_KETAMINE = 'Abundances (Grouped): ketamine'
COL_CATEGORY = 'category'
COL_DIRECTION = 'direction'
COL_SIGNIFICANT = 'significant'

# ======================================================================
# 2. PARSE PATHWAY REPORT
# ======================================================================

def parse_pathway_report(report_path):
    '''
    Parse the pathway analysis report to build a gene-to-terms mapping.
    Returns dict: { gene_symbol: [ (source, term_id, term_name, fdr, gene_count), ... ] }
    '''
    gene_to_terms = {}

    if not os.path.exists(report_path):
        print(f'  Warning: pathway report not found at {report_path}')
        return gene_to_terms

    with open(report_path, 'r') as f:
        text = f.read()

    # -- parse GO:BP terms --
    # Pattern: "  N. term_name\n      FDR = X, Gene count = N\n      Genes: gene1, gene2, ..."
    # Also handle KEGG: "  N. term_name (KEGG:XXXXX)\n     FDR = X, Gene count = N\n     Genes: ..."
    go_section = re.search(
        r'UPREGULATED - Top 10 GO:BP terms:\n(.*?)(?=\nDOWNREGULATED|\n---)',
        text, re.DOTALL
    )
    kegg_section = re.search(
        r'UPREGULATED - All KEGG pathways:\n(.*?)(?=\nDOWNREGULATED|\n---)',
        text, re.DOTALL
    )

    def parse_section(section_text, source):
        '''Parse a section of the report into term entries.'''
        entries = []
        # Match numbered entries with their FDR and gene lines
        pattern = re.compile(
            r'\s*\d+\.\s+(.+?)\n'
            r'\s+FDR\s*=\s*([\d.e+-]+),\s*Gene count\s*=\s*(\d+)\n'
            r'\s+Genes:\s*(.+)',
            re.MULTILINE
        )
        for match in pattern.finditer(section_text):
            term_name_raw = match.group(1).strip()
            fdr = float(match.group(2))
            gene_count = int(match.group(3))
            genes_str = match.group(4).strip()
            gene_list = [g.strip() for g in genes_str.split(',')]

            # Extract term_id from name if present (KEGG format)
            term_id_match = re.search(r'\((\w+:\w+)\)', term_name_raw)
            if term_id_match:
                term_id = term_id_match.group(1)
                term_name = re.sub(r'\s*\(\w+:\w+\)', '', term_name_raw).strip()
            else:
                term_id = ''
                term_name = term_name_raw

            for gene in gene_list:
                if gene not in gene_to_terms:
                    gene_to_terms[gene] = []
                gene_to_terms[gene].append({
                    'source': source,
                    'term_id': term_id,
                    'term_name': term_name,
                    'fdr': fdr,
                    'gene_count': gene_count,
                })

    if go_section:
        parse_section(go_section.group(1), 'GO:BP')
    if kegg_section:
        parse_section(kegg_section.group(1), 'KEGG')

    return gene_to_terms

# ======================================================================
# 3. LOAD DATA
# ======================================================================

def load_data():
    '''Load all data sources.'''
    proteins_df = pd.read_csv(PROTEINS_FILE)
    markers_df = pd.read_csv(MARKERS_FILE)
    gene_to_terms = parse_pathway_report(REPORT_FILE)
    return proteins_df, markers_df, gene_to_terms

# ======================================================================
# 4. GENERATE PROTEIN SUMMARY
# ======================================================================

def generate_summary(gene_symbol, proteins_df, markers_df, gene_to_terms):
    '''Generate a comprehensive summary string for a given gene symbol.'''

    lines = []
    sep = '=' * 70

    # -- case-insensitive lookup --
    match = proteins_df[proteins_df[COL_GENE].str.lower() == gene_symbol.lower()]

    if len(match) == 0:
        lines.append(f'{sep}')
        lines.append(f'  {gene_symbol.upper()} -- NOT FOUND IN DATASET')
        lines.append(f'{sep}')
        lines.append(f'  No protein with gene symbol "{gene_symbol}" was detected')
        lines.append(f'  in the proteomics dataset (903 proteins).')
        return '\n'.join(lines)

    row = match.iloc[0]
    gene = row[COL_GENE]

    # -- header --
    lines.append(f'{sep}')
    lines.append(f'  PROTEIN SUMMARY: {gene.upper()}')
    lines.append(f'{sep}')

    # -- basic info --
    lines.append('')
    lines.append('--- Identification ---')
    lines.append(f'  Gene Symbol:    {gene}')
    lines.append(f'  Accession:      {row[COL_ACCESSION]}')
    desc = row[COL_DESC]
    # Clean up description (truncate after OS=)
    desc_short = re.split(r'\s+OS=', str(desc))[0]
    lines.append(f'  Description:    {desc_short}')

    # -- category and significance --
    lines.append('')
    lines.append('--- Classification ---')
    category = row[COL_CATEGORY]
    direction = row.get(COL_DIRECTION, 'N/A')
    significant = row.get(COL_SIGNIFICANT, 'N/A')
    lines.append(f'  Category:       {category}')
    lines.append(f'  Direction:      {direction}')
    lines.append(f'  Significant:    {significant}')

    # -- quantitative data --
    lines.append('')
    lines.append('--- Differential Abundance ---')
    ratio = row[COL_RATIO]
    log2fc = row[COL_LOG2FC]
    pval = row[COL_PVAL]
    variability = row.get(COL_VARIABILITY, np.nan)
    ctrl_abund = row[COL_CONTROL]
    ket_abund = row[COL_KETAMINE]

    if 'presence_absence' in str(category):
        lines.append(f'  Abundance Ratio:    {ratio} (placeholder -- presence/absence protein)')
        lines.append(f'  log2(FC):           {log2fc:.3f} (placeholder)')
        lines.append(f'  Adj. P-Value:       {pval:.2e} (placeholder)')
        lines.append(f'  Variability:        N/A (presence/absence)')
    else:
        lines.append(f'  Abundance Ratio:    {ratio}')
        lines.append(f'  log2(FC):           {log2fc:.3f}')
        lines.append(f'  Adj. P-Value:       {pval:.2e}')
        if pd.notna(variability):
            lines.append(f'  Variability [%]:    {variability:.1f}%')
        else:
            lines.append(f'  Variability [%]:    N/A')

    lines.append('')
    lines.append('--- Grouped Abundances ---')
    ctrl_str = f'{ctrl_abund:.1f}' if pd.notna(ctrl_abund) else 'N/D (not detected)'
    ket_str = f'{ket_abund:.1f}' if pd.notna(ket_abund) else 'N/D (not detected)'
    lines.append(f'  Control:            {ctrl_str}')
    lines.append(f'  Ketamine:           {ket_str}')

    # -- cell-type marker status --
    lines.append('')
    lines.append('--- Cell-Type Marker Status ---')
    marker_match = markers_df[markers_df['gene_symbol'].str.lower() == gene.lower()]
    if len(marker_match) > 0:
        for _, mrow in marker_match.iterrows():
            lines.append(f'  Cell type:          {mrow["cell_type"]}')
            lines.append(f'  Marker category:    {mrow["category"]}')
            lines.append(f'  Specificity notes:  {mrow["specificity_notes"]}')
    else:
        lines.append(f'  Not listed as a canonical cell-type marker.')

    # -- pathway/GO term memberships --
    lines.append('')
    lines.append('--- Pathway & GO Term Memberships ---')

    # find all terms containing this gene (case-insensitive lookup)
    terms = None
    for g, t in gene_to_terms.items():
        if g.lower() == gene.lower():
            terms = t
            break

    if terms:
        # separate GO:BP and KEGG
        go_terms = [t for t in terms if t['source'] == 'GO:BP']
        kegg_terms = [t for t in terms if t['source'] == 'KEGG']

        if go_terms:
            lines.append(f'')
            lines.append(f'  GO:BP terms ({len(go_terms)}):')
            for t in sorted(go_terms, key=lambda x: x['fdr']):
                tid = f' ({t["term_id"]})' if t['term_id'] else ''
                lines.append(f'    - {t["term_name"]}{tid}')
                lines.append(f'      FDR = {t["fdr"]:.2e}, Gene count = {t["gene_count"]}')

        if kegg_terms:
            lines.append(f'')
            lines.append(f'  KEGG pathways ({len(kegg_terms)}):')
            for t in sorted(kegg_terms, key=lambda x: x['fdr']):
                tid = f' ({t["term_id"]})' if t['term_id'] else ''
                lines.append(f'    - {t["term_name"]}{tid}')
                lines.append(f'      FDR = {t["fdr"]:.2e}, Gene count = {t["gene_count"]}')
    else:
        lines.append(f'  No enriched GO:BP or KEGG terms include this protein.')
        lines.append(f'  (Note: only top 10 GO:BP and all KEGG terms from the report are indexed.)')

    lines.append('')
    return '\n'.join(lines)

# ======================================================================
# 5. MAIN
# ======================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive protein summaries from proteomics data.',
        epilog='Examples:\n'
               '  python scripts/protein_summary.py Snap25\n'
               '  python scripts/protein_summary.py Snap25 Vamp2 Syt1\n'
               '  python scripts/protein_summary.py --all-significant\n'
               '  python scripts/protein_summary.py Snap25 --output results/protein_summaries/',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('genes', nargs='*', help='Gene symbol(s) to summarize')
    parser.add_argument('--all-significant', action='store_true',
                        help='Summarize all significant proteins (p <= 0.05)')
    parser.add_argument('--all-upregulated', action='store_true',
                        help='Summarize all upregulated significant proteins')
    parser.add_argument('--all-downregulated', action='store_true',
                        help='Summarize all downregulated significant proteins')
    default_output = os.path.join(BASE_DIR, 'results', 'protein_summaries')
    parser.add_argument('--output', '-o', type=str, default=default_output,
                        help=f'Output directory (default: results/protein_summaries/)')
    parser.add_argument('--no-save', action='store_true',
                        help='Print to stdout only, do not save files')
    args = parser.parse_args()

    if not args.genes and not args.all_significant and not args.all_upregulated and not args.all_downregulated:
        parser.print_help()
        sys.exit(1)

    # -- load data --
    proteins_df, markers_df, gene_to_terms = load_data()

    # -- determine gene list --
    gene_list = list(args.genes) if args.genes else []

    if args.all_significant:
        sig = proteins_df[proteins_df[COL_SIGNIFICANT] == True]
        gene_list = sorted(sig[COL_GENE].tolist())
    elif args.all_upregulated:
        sig = proteins_df[(proteins_df[COL_SIGNIFICANT] == True) & (proteins_df[COL_DIRECTION] == 'up')]
        gene_list = sorted(sig[COL_GENE].tolist())
    elif args.all_downregulated:
        sig = proteins_df[(proteins_df[COL_SIGNIFICANT] == True) & (proteins_df[COL_DIRECTION] == 'down')]
        gene_list = sorted(sig[COL_GENE].tolist())

    # -- also include presence/absence proteins in --all flags --
    if args.all_significant or args.all_upregulated:
        pa_up = proteins_df[proteins_df[COL_CATEGORY].str.contains('ketamine_specific', na=False)]
        gene_list = sorted(set(gene_list + pa_up[COL_GENE].tolist()))
    if args.all_significant or args.all_downregulated:
        pa_down = proteins_df[proteins_df[COL_CATEGORY].str.contains('control_specific', na=False)]
        gene_list = sorted(set(gene_list + pa_down[COL_GENE].tolist()))

    if not gene_list:
        print('No proteins to summarize.')
        sys.exit(0)

    # -- create output directory --
    save_to_file = not args.no_save
    if save_to_file:
        os.makedirs(args.output, exist_ok=True)

    # -- generate summaries --
    for gene in gene_list:
        summary = generate_summary(gene, proteins_df, markers_df, gene_to_terms)
        print(summary)

        if save_to_file:
            filename = f'{gene.lower()}_summary.txt'
            filepath = os.path.join(args.output, filename)
            with open(filepath, 'w') as f:
                f.write(summary)
            print(f'  >> Saved: {filepath}\n')


if __name__ == '__main__':
    main()
