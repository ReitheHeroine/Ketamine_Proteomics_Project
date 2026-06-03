#!/usr/bin/env python3
# title: revigo_fetch.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-29
# last modified: 2026-05-29
#
# purpose:
#   Programmatically fetch a REVIGO simplified-term table (Supek et al., 2011)
#   from the REVIGO web service (http://revigo.irb.hr/) without manual web-tool
#   interaction. Reads a g:Profiler ORA results CSV (the output of
#   pathway_analysis.py), extracts the GO term ID + adjusted p-value pairs,
#   POSTs them to REVIGO's StartJob endpoint, polls QueryJob for completion,
#   and writes the simplified TSV table in the same schema produced by the
#   REVIGO web tool. The output file can be consumed directly by
#   pathway_analysis.py --revigo mode.
#
#   This script automates the 2026-02-02 manual REVIGO workflow. Algorithm
#   and parameter defaults match the REVIGO web tool, so resulting term lists
#   are equivalent to a manual web-tool run with the same parameters.
#
# inputs:
#   - ORA results CSV with columns: 'term_id' (GO:XXXXXXX) and 'fdr_pvalue'
#     (numeric). Default consumer: results/pathway_analysis/upregulated/
#     upregulated_GO_BP.csv.
#
# outputs:
#   - REVIGO simplified TSV at the path given by --output. Schema:
#       TermID, Name, Value, LogSize, Frequency, Uniqueness, Dispensability,
#       Representative
#     (the exact format that pathway_analysis.py::parse_revigo_output expects)
#
# usage:
#   python revigo_fetch.py \
#       --input ../results/pathway_analysis/upregulated/upregulated_GO_BP.csv \
#       --output ../results/pathway_analysis/revigo/Revigo_BP_Table.tsv \
#       --species_taxon 10090 \
#       --cutoff 0.7 \
#       --measure SIMREL
#
#   copy/paste: python revigo_fetch.py --input ../results/pathway_analysis/upregulated/upregulated_GO_BP.csv --output ../results/pathway_analysis/revigo/Revigo_BP_Table.tsv --species_taxon 10090 --cutoff 0.7 --measure SIMREL
#
# dependencies:
#   pip install requests pandas
#
# notes:
#   - Internet required at run time (POSTs to revigo.irb.hr).
#   - REVIGO assigns server-side numeric JobIDs and tracks them via session
#     cookie; this script preserves cookies across the POST and GET in a
#     single requests.Session().
#   - Default cutoff = 0.7 (REVIGO "medium" simplification). Pass a different
#     value via --cutoff: 0.4 (large), 0.5 (medium-large), 0.7 (medium),
#     0.9 (small). Smaller cutoff = more aggressive clustering = fewer
#     representative terms in the output.
#   - Default similarity measure = SIMREL (REVIGO's recommended default).
#     Other options: RESNIK, LIN, JIANG.
#   - Default species taxon = 10090 (Mus musculus). Use 0 for whole UniProt.
#   - namespace = 1 returns GO:BP; 2 returns GO:CC; 3 returns GO:MF.
#     Defaults to 1 because this is the only namespace the GO:BP ORA
#     results file populates.

import argparse
import json
import os
import sys
import time
from typing import List, Tuple

import pandas as pd
import requests


# --- Configuration & Constants ---

REVIGO_BASE_URL = 'http://revigo.irb.hr'
REVIGO_START_ENDPOINT = f'{REVIGO_BASE_URL}/StartJob'
REVIGO_QUERY_ENDPOINT = f'{REVIGO_BASE_URL}/QueryJob'

# Input CSV columns
GO_ID_COL = 'term_id'
PVAL_COL = 'fdr_pvalue'


# --- Argument Parsing ---

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Fetch a REVIGO simplified term table via the REVIGO HTTP API.'
    )
    p.add_argument('--input', required=True,
                   help='Path to ORA results CSV with term_id and fdr_pvalue columns.')
    p.add_argument('--output', required=True,
                   help='Path to write the REVIGO simplified TSV table.')
    p.add_argument('--species_taxon', type=int, default=10090,
                   help='NCBI taxon ID (10090 = mouse, 9606 = human, '
                        '0 = whole UniProt). Default: 10090.')
    p.add_argument('--cutoff', type=float, default=0.7,
                   help='REVIGO similarity cutoff: 0.4 (large), 0.5, 0.7 (medium), '
                        '0.9 (small). Smaller cutoff produces fewer representatives. '
                        'Default: 0.7.')
    p.add_argument('--measure', default='SIMREL',
                   choices=['SIMREL', 'RESNIK', 'LIN', 'JIANG'],
                   help='Semantic similarity measure. Default: SIMREL.')
    p.add_argument('--value_type', default='pvalue',
                   choices=['pvalue', 'higher', 'lower'],
                   help='REVIGO value-type interpretation. Default: pvalue.')
    p.add_argument('--namespace', type=int, default=1, choices=[1, 2, 3],
                   help='1 = GO:BP, 2 = GO:CC, 3 = GO:MF. Default: 1.')
    p.add_argument('--max_wait_seconds', type=int, default=120,
                   help='Maximum seconds to wait for REVIGO job completion. '
                        'Default: 120.')
    p.add_argument('--poll_interval_seconds', type=float, default=2.0,
                   help='Polling interval for REVIGO job status. Default: 2.0.')
    return p.parse_args()


# --- Data Preparation ---

def load_go_terms(input_path: str) -> List[Tuple[str, float]]:
    '''Read the ORA results CSV and return (GO_ID, p-value) tuples.'''
    df = pd.read_csv(input_path)
    missing = [c for c in (GO_ID_COL, PVAL_COL) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Input file {input_path} is missing required column(s): {missing}. "
            f"Available: {list(df.columns)}"
        )
    # REVIGO expects valid GO:NNNNNNN identifiers and positive p-values.
    df = df.dropna(subset=[GO_ID_COL, PVAL_COL])
    df = df[df[GO_ID_COL].str.startswith('GO:')]
    df = df[df[PVAL_COL] > 0]
    if df.empty:
        raise ValueError(f"No valid GO terms with p-values found in {input_path}.")
    return list(zip(df[GO_ID_COL].tolist(), df[PVAL_COL].tolist()))


# --- REVIGO API ---

def submit_revigo_job(session: requests.Session, terms: List[Tuple[str, float]],
                      species_taxon: int, cutoff: float, measure: str,
                      value_type: str) -> int:
    '''POST to REVIGO StartJob; return the assigned JobID.'''
    # REVIGO expects newline-separated "GO:ID p-value" pairs.
    go_list = '\n'.join(f'{gid} {pval}' for gid, pval in terms)
    resp = session.post(
        REVIGO_START_ENDPOINT,
        data={
            'cutoff': cutoff,
            'valueType': value_type,
            'speciesTaxon': species_taxon,
            'measure': measure,
            'goList': go_list,
        },
        timeout=60,
    )
    resp.raise_for_status()
    payload = resp.json()
    if payload.get('jobid', -1) < 0:
        raise RuntimeError(
            f"REVIGO StartJob refused the request: {payload.get('message', payload)}"
        )
    return int(payload['jobid'])


def poll_until_finished(session: requests.Session, job_id: int,
                        max_wait_seconds: int, poll_interval_seconds: float) -> None:
    '''Poll QueryJob until REVIGO reports the job as finished.'''
    deadline = time.time() + max_wait_seconds
    while time.time() < deadline:
        resp = session.get(
            REVIGO_QUERY_ENDPOINT,
            params={'jobid': job_id, 'type': 'jstatus'},
            timeout=30,
        )
        resp.raise_for_status()
        status = resp.json()
        if status.get('running', 1) == 0:
            return
        time.sleep(poll_interval_seconds)
    raise TimeoutError(
        f"REVIGO job {job_id} did not finish within {max_wait_seconds} seconds."
    )


def fetch_revigo_table(session: requests.Session, job_id: int,
                       namespace: int) -> str:
    '''GET QueryJob with type=Table and return the TSV body.'''
    resp = session.get(
        REVIGO_QUERY_ENDPOINT,
        params={'jobid': job_id, 'namespace': namespace, 'type': 'Table'},
        timeout=60,
    )
    resp.raise_for_status()
    body = resp.text
    if not body.startswith('TermID'):
        raise RuntimeError(
            f"Unexpected REVIGO Table response (first 200 chars): {body[:200]!r}"
        )
    return body


# --- Main ---

def main() -> None:
    args = parse_args()

    # --- Load input ---
    print(f'[revigo_fetch] Loading {args.input}', file=sys.stderr)
    terms = load_go_terms(args.input)
    print(f'[revigo_fetch] {len(terms)} GO terms loaded', file=sys.stderr)

    # --- Submit to REVIGO ---
    print(f'[revigo_fetch] Submitting to REVIGO '
          f'(taxon={args.species_taxon}, cutoff={args.cutoff}, '
          f'measure={args.measure})', file=sys.stderr)
    session = requests.Session()
    job_id = submit_revigo_job(
        session, terms,
        species_taxon=args.species_taxon,
        cutoff=args.cutoff,
        measure=args.measure,
        value_type=args.value_type,
    )
    print(f'[revigo_fetch] JobID = {job_id}', file=sys.stderr)

    # --- Wait for completion ---
    print(f'[revigo_fetch] Polling for completion...', file=sys.stderr)
    poll_until_finished(
        session, job_id,
        max_wait_seconds=args.max_wait_seconds,
        poll_interval_seconds=args.poll_interval_seconds,
    )

    # --- Fetch table ---
    tsv = fetch_revigo_table(session, job_id, args.namespace)

    # --- Write output ---
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, 'w') as f:
        f.write(tsv)
    n_rows = tsv.count('\n') - 1   # subtract header
    print(f'[revigo_fetch] Wrote {n_rows} rows to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
