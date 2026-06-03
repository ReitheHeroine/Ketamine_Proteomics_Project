#!/usr/bin/env python3
"""
Build the ketamine proteomics thesis table workbook.

Style conventions follow project_notes/figure_and_table_style_guide.Rmd
(updated 2026-05-27 with Sections 11-12 governing what does and does not
live inside table cells).

Key compliance points:
- Times New Roman 10 pt body, 11 pt bold headers (cells).
- Caption above the table region: Times New Roman 12 pt bold, NO borders,
  NO fill. Sits as a "Word paragraph above the table" surrogate in xlsx so
  the user can copy it to Word and apply the "Table Title" style.
- Notes below the table region: Times New Roman 10 pt, NO borders, NO fill.
  Italic "Note." prefix via rich text. User applies "Table Note" style in Word.
- Triple-rule APA layout (top medium rule above header, thin rule below
  header, bottom medium rule below last data row; no vertical rules; no
  interior horizontal rules between body rows).
- Gene symbols rendered in protein convention (UPPERCASE ROMAN) per Section
  4.3 because the symbols label protein abundance from mass spec.
- Log2 fold change: 2 decimal places, explicit sign (+1.90 / -0.47).
- Adj. p-value: 2 sig figs, scientific notation with Unicode superscripts
  when below 0.001 (e.g., 1.2 x 10^-6).
- Statistical symbols (p, t, n) italicized via rich text wherever they
  appear in column headers, captions, or notes.
- Missing values: em-dash.
"""

from pathlib import Path
import csv
import re

from openpyxl import Workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path("/Users/reina/Library/Mobile Documents/com~apple~CloudDocs/"
                    "Blanco-Suárez Lab/ketamine_project")
SIG_PATH = PROJECT_ROOT / "results" / "quantitative" / "significant.csv"
KET_PA_PATH = PROJECT_ROOT / "results" / "presence_absence" / "ketamine_specific.csv"
CTRL_PA_PATH = PROJECT_ROOT / "results" / "presence_absence" / "control_specific.csv"
CATEGORIZED_PATH = PROJECT_ROOT / "data" / "all_proteins_categorized.csv"
RAW_CTRL_PATH = PROJECT_ROOT / "data" / "proteins_control.csv"
RAW_BOTH_PATH = PROJECT_ROOT / "data" / "proteins_control_ketamine.csv"
OUTPUT_PATH = PROJECT_ROOT / "results" / "tables" / "ketamine_proteomics_tables.xlsx"

# -----------------------------------------------------------------------------
# Style constants (from figure_and_table_style_guide.Rmd)
# -----------------------------------------------------------------------------
FONT_FAMILY = "Times New Roman"
SIZE_BODY = 10
SIZE_HEADER = 11
SIZE_CAPTION = 12  # Section 12.2: 12 pt bold (was 11 in previous build)
SIZE_NOTE = 10     # Section 12.3: 10 pt (was 9 in previous build)

FILL_HEADER = PatternFill("solid", fgColor="F2F2F2")  # very light gray
THIN = Side(style="thin", color="000000")
MEDIUM = Side(style="medium", color="000000")
NO_BORDER = Border()

# -----------------------------------------------------------------------------
# Rich-text helpers
# -----------------------------------------------------------------------------
def _font(size, bold=False, italic=False):
    return InlineFont(rFont=FONT_FAMILY, sz=size, b=bold, i=italic)


def md_to_rich(text, size, bold=False):
    """Parse *italic* markers in `text` and return CellRichText with the
    surrounding text in the requested style and the asterisk-delimited spans
    in italic."""
    parts = re.split(r'(\*[^*]+\*)', text)
    blocks = []
    regular = _font(size, bold=bold, italic=False)
    italic = _font(size, bold=bold, italic=True)
    for part in parts:
        if not part:
            continue
        if part.startswith('*') and part.endswith('*') and len(part) > 2:
            blocks.append(TextBlock(italic, part[1:-1]))
        else:
            blocks.append(TextBlock(regular, part))
    return CellRichText(blocks)


# -----------------------------------------------------------------------------
# Pathway categories (best guesses; user revises in cells)
# -----------------------------------------------------------------------------
PATHWAY_CATEGORY = {
    # SNARE complex / synaptic vesicle cycle (Tier 1)
    "Snap25": "SNARE complex",
    "Stx1a":  "SNARE complex",
    "Stx1b":  "SNARE complex",
    "Vamp2":  "SNARE complex",
    "Syt1":   "Synaptic vesicle / Ca2+ sensor",
    "Stxbp1": "SNARE regulator (Munc18)",
    "Sv2a":   "Synaptic vesicle protein",
    "Syp":    "Synaptic vesicle protein",
    "Syn2":   "Synaptic vesicle protein",
    "Snap91": "Synaptic vesicle endocytosis",
    "Dnm1":   "Synaptic vesicle endocytosis",
    "Prrt2":  "Synaptic vesicle release regulator",
    # V-ATPase / vesicle acidification
    "Atp6v1g2": "V-ATPase",
    "Atp6v1b2": "V-ATPase",
    # Plasticity / signaling (Tier 2)
    "Gap43":  "Plasticity / signaling",
    "Basp1":  "Plasticity / signaling",
    "Pik3r1": "Plasticity / signaling",
    "Map1b":  "Cytoskeletal / plasticity",
    "Crmp1":  "Axon guidance / plasticity",
    "Nrgn":   "Plasticity / signaling (CaM-binding)",
    # Neurofilament / cytoskeletal (Tier 3 contamination risk)
    "Ina":  "Neurofilament (contamination risk)",
    "Nefl": "Neurofilament (contamination risk)",
    "Nefm": "Neurofilament (contamination risk)",
    "Tubb3": "Neuronal cytoskeleton (contamination risk)",
    # Myelin / oligodendrocyte (Tier 3 contamination risk)
    "Plp1": "Myelin (contamination risk)",
    "Mag":  "Myelin (contamination risk)",
    "Mog":  "Myelin (contamination risk)",
    "Mbp":  "Myelin (contamination risk)",
    "Cnp":  "Myelin (contamination risk)",
    # Ion homeostasis
    "Atp1a3": "Na+/K+ ATPase",
    # Other neuronal / synaptic
    "Thy1": "Neuronal surface marker",
    # Other proteins of interest
    "Snca": "Synaptic / neurodegeneration",
    "Penk": "Neuropeptide precursor",
    "Golga5": "Golgi apparatus",
    "Abcb11": "Transporter (likely artifact)",
    "Atp5if1": "Mitochondrial ATPase regulator",
    # Downregulated
    "Aldh7a1": "Aldehyde metabolism",
    "Akap6":   "PKA anchoring",
    "Tlr1":    "Innate immunity",
    "Jcad":    "Cytoskeletal / vascular",
    # Control-specific P/A (functional best-guesses; user to revise)
    "Carhsp1":  "Calcium signaling",
    "Tspan7":   "Membrane / tetraspanin",
    "Cpe":      "Neuropeptide processing",
    "Cnpy2":    "Secretory pathway",
    "Ftl1":     "Iron storage",
    "Adh5":     "Alcohol / formaldehyde metabolism",
    "Cnn3":     "Cytoskeletal (actin-binding)",
    "Sf1":      "RNA splicing",
    "Top2b":    "DNA topoisomerase",
    "Cdc42":    "Small GTPase / cytoskeletal",
    "Tufm":     "Mitochondrial translation",
    "Golgb1":   "Golgi apparatus",
    "Pgm1":     "Glycolysis / glycogen metabolism",
    "Ssr1":     "ER translocon",
    "Zbtb20":   "Transcription factor",
    "Cast":     "Calpain inhibitor",
    "Fubp1":    "Transcription regulator",
    "Omp":      "Olfactory marker (contamination risk)",
    "Pfkp":     "Glycolysis",
    "Rpl3":     "Ribosomal protein",
    "Rpl13":    "Ribosomal protein",
    "Igkv5-37": "Immunoglobulin (contamination risk)",
    "Sap18":    "HDAC / Sin3 complex",
    "Psmc5":    "Proteasome regulatory subunit",
    "Myo6":     "Motor protein (myosin)",
    "Asb15":    "Ubiquitin ligase adapter",
    "Eif3a":    "Translation initiation",
    "G3bp1":    "Stress granule / RNA-binding",
}


# -----------------------------------------------------------------------------
# Data processing
# -----------------------------------------------------------------------------
def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def format_log2fc(v):
    if v is None:
        return "—"
    return f"{v:+.2f}"


def format_pvalue(v):
    if v is None:
        return "—"
    if v < 1e-300:
        return "< 1 × 10⁻³⁰⁰"
    if v >= 0.001:
        return f"{v:.2g}"
    mantissa_str = f"{v:.1e}"
    mantissa, exponent = mantissa_str.split("e")
    exp_int = int(exponent)
    sup_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    exp_sup = str(exp_int).translate(sup_map)
    return f"{mantissa} × 10{exp_sup}"


def gene_to_protein_label(gene_symbol):
    if not gene_symbol:
        return "—"
    return gene_symbol.upper()


def strip_description(desc):
    if not desc:
        return "—"
    return desc.split(" OS=")[0].strip()


# -----------------------------------------------------------------------------
# Worksheet construction
# -----------------------------------------------------------------------------
def write_caption_paragraph(ws, row, caption_text, n_cols):
    """Caption sits ABOVE the table region. No borders, no fill. Surrogate
    for a Word paragraph above the table (apply 'Table Title' style in Word)."""
    cell = ws.cell(row=row, column=1)
    cell.value = md_to_rich(caption_text, SIZE_CAPTION, bold=True)
    cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    cell.border = NO_BORDER
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=n_cols)
    ws.row_dimensions[row].height = 24


def write_header_row(ws, row, headers, alignments):
    """Header row: top medium rule above, thin rule below. Bold 11 pt with
    light gray fill. Headers passed in as CellRichText."""
    for col_idx, (header, align) in enumerate(zip(headers, alignments), start=1):
        cell = ws.cell(row=row, column=col_idx)
        cell.value = header
        cell.alignment = Alignment(horizontal=align, vertical="center", wrap_text=True)
        cell.border = Border(top=MEDIUM, bottom=THIN)
        cell.fill = FILL_HEADER
    ws.row_dimensions[row].height = 30


def write_body_row(ws, row, values, alignments, is_last_row=False):
    """Body row. No internal horizontal rules. Bottom medium rule on final row."""
    for col_idx, (val, align) in enumerate(zip(values, alignments), start=1):
        cell = ws.cell(row=row, column=col_idx, value=val)
        cell.font = Font(name=FONT_FAMILY, size=SIZE_BODY)
        cell.alignment = Alignment(horizontal=align, vertical="center", wrap_text=True)
        if is_last_row:
            cell.border = Border(bottom=MEDIUM)
        else:
            cell.border = NO_BORDER


def write_note_paragraph(ws, row, note_text, n_cols):
    """Note sits BELOW the table region. No borders, no fill. Surrogate for
    a Word paragraph below the table (apply 'Table Note' style in Word).
    Italic 'Note.' prefix and statistical symbols handled via rich text."""
    cell = ws.cell(row=row, column=1)
    cell.value = md_to_rich(note_text, SIZE_NOTE, bold=False)
    cell.alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)
    cell.border = NO_BORDER
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=n_cols)
    ws.row_dimensions[row].height = 110


def set_column_widths(ws, widths):
    for col_idx, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(col_idx)].width = w


# -----------------------------------------------------------------------------
# Column headers (rich text for italic statistical symbols and subscripts)
# -----------------------------------------------------------------------------
def header_cells_quant():
    return [
        md_to_rich("Gene symbol", SIZE_HEADER, bold=True),
        md_to_rich("UniProt accession", SIZE_HEADER, bold=True),
        md_to_rich("Protein name", SIZE_HEADER, bold=True),
        md_to_rich("Log₂ fold change", SIZE_HEADER, bold=True),
        md_to_rich("Adj. *p*-value", SIZE_HEADER, bold=True),
        md_to_rich("Pathway category", SIZE_HEADER, bold=True),
    ]


def header_cells_pa():
    return [
        md_to_rich("Gene symbol", SIZE_HEADER, bold=True),
        md_to_rich("UniProt accession", SIZE_HEADER, bold=True),
        md_to_rich("Protein name", SIZE_HEADER, bold=True),
        md_to_rich("Detection status", SIZE_HEADER, bold=True),
        md_to_rich("Pathway category", SIZE_HEADER, bold=True),
    ]


ALIGN_QUANT = ["left", "left", "left", "right", "right", "left"]
WIDTHS_QUANT = [13, 16, 50, 18, 20, 32]

ALIGN_PA = ["left", "left", "left", "left", "left"]
WIDTHS_PA = [13, 16, 50, 24, 32]

ALIGN_QPCR = ["left", "left", "right", "right", "right"]
WIDTHS_QPCR = [12, 28, 26, 26, 24]

ALIGN_MR = ["left", "left", "left", "left", "right"]
WIDTHS_MR = [13, 16, 50, 32, 18]


# -----------------------------------------------------------------------------
# Table builders
# -----------------------------------------------------------------------------
def build_quant_sheet(ws, caption_text, proteins, note_text):
    n_cols = 6
    set_column_widths(ws, WIDTHS_QUANT)

    # Row 1: caption paragraph (above table region; no borders)
    write_caption_paragraph(ws, 1, caption_text, n_cols)

    # Row 2: header row (top + bottom rules; gray fill)
    write_header_row(ws, 2, header_cells_quant(), ALIGN_QUANT)

    # Rows 3+: body
    body_start = 3
    last_idx = len(proteins) - 1
    for i, p in enumerate(proteins):
        gene = p["Gene Symbol"]
        values = [
            gene_to_protein_label(gene),
            p["Accession"],
            strip_description(p["Description"]),
            format_log2fc(to_float(p["log2_fold_change"])),
            format_pvalue(to_float(p["Abundance Ratio Adj. P-Value: (ketamine) / (control)"])),
            PATHWAY_CATEGORY.get(gene, ""),
        ]
        write_body_row(ws, body_start + i, values, ALIGN_QUANT,
                       is_last_row=(i == last_idx))

    # Note paragraph immediately below the table (no borders)
    note_row = body_start + len(proteins)
    write_note_paragraph(ws, note_row, note_text, n_cols)

    ws.freeze_panes = f"A{body_start}"


def header_cells_qpcr():
    return [
        md_to_rich("Marker", SIZE_HEADER, bold=True),
        md_to_rich("Target cell type", SIZE_HEADER, bold=True),
        md_to_rich("ACSA-2-positive fraction", SIZE_HEADER, bold=True),
        md_to_rich("Flow-through fraction", SIZE_HEADER, bold=True),
        md_to_rich("Fold below *Slc1a3* (ACSA-2+)", SIZE_HEADER, bold=True),
    ]


def _qpcr_num_format(v):
    if v is None:
        return "General"
    if abs(v) >= 0.001:
        return "0.0000"
    return "0.00E+00"


def _qpcr_fold_format(v):
    if v is None or v == 1.0:
        return "0"
    return "#,##0"


def build_qpcr_sheet(ws, caption_text, rows, note_text):
    n_cols = 5
    set_column_widths(ws, WIDTHS_QPCR)

    write_caption_paragraph(ws, 1, caption_text, n_cols)
    write_header_row(ws, 2, header_cells_qpcr(), ALIGN_QPCR)

    body_start = 3
    last_idx = len(rows) - 1
    for i, r in enumerate(rows):
        marker_rt = md_to_rich(f"*{r['marker']}*", SIZE_BODY, bold=False)
        acsa = r["acsa2_pos"]
        flow = r["flow_through"]
        fold_below = 1.0 if acsa == 1.0 else round(1.0 / acsa)
        values = [marker_rt, r["cell_type"], acsa, flow, fold_below]
        write_body_row(ws, body_start + i, values, ALIGN_QPCR,
                       is_last_row=(i == last_idx))
        row_idx = body_start + i
        # Rich-text marker cell keeps left-aligned wrap
        ws.cell(row=row_idx, column=1).alignment = Alignment(
            horizontal="left", vertical="center", wrap_text=True)
        # Number formats per cell so small values stay readable
        ws.cell(row=row_idx, column=3).number_format = _qpcr_num_format(acsa)
        ws.cell(row=row_idx, column=4).number_format = _qpcr_num_format(flow)
        ws.cell(row=row_idx, column=5).number_format = _qpcr_fold_format(acsa)

    note_row = body_start + len(rows)
    write_note_paragraph(ws, note_row, note_text, n_cols)
    ws.freeze_panes = f"A{body_start}"


def header_cells_mr():
    return [
        md_to_rich("Gene symbol", SIZE_HEADER, bold=True),
        md_to_rich("UniProt accession", SIZE_HEADER, bold=True),
        md_to_rich("Protein name", SIZE_HEADER, bold=True),
        md_to_rich("Detection pattern", SIZE_HEADER, bold=True),
        md_to_rich("# Unique peptides", SIZE_HEADER, bold=True),
    ]


def build_mr_sheet(ws, caption_text, proteins, note_text):
    n_cols = 5
    set_column_widths(ws, WIDTHS_MR)

    write_caption_paragraph(ws, 1, caption_text, n_cols)
    write_header_row(ws, 2, header_cells_mr(), ALIGN_MR)

    body_start = 3
    last_idx = len(proteins) - 1
    for i, p in enumerate(proteins):
        values = [
            gene_to_protein_label(p["Gene Symbol"]),
            p["Accession"],
            strip_description(p["Description"]),
            p["detection_pattern"],
            p["unique_peptides"],
        ]
        write_body_row(ws, body_start + i, values, ALIGN_MR,
                       is_last_row=(i == last_idx))

    note_row = body_start + len(proteins)
    write_note_paragraph(ws, note_row, note_text, n_cols)
    ws.freeze_panes = f"A{body_start}"


def build_pa_sheet(ws, caption_text, proteins, condition_label, note_text):
    n_cols = 5
    set_column_widths(ws, WIDTHS_PA)

    write_caption_paragraph(ws, 1, caption_text, n_cols)
    write_header_row(ws, 2, header_cells_pa(), ALIGN_PA)

    body_start = 3
    last_idx = len(proteins) - 1
    for i, p in enumerate(proteins):
        gene = p["Gene Symbol"]
        values = [
            gene_to_protein_label(gene),
            p["Accession"],
            strip_description(p["Description"]),
            f"Detected only in {condition_label}",
            PATHWAY_CATEGORY.get(gene, ""),
        ]
        write_body_row(ws, body_start + i, values, ALIGN_PA,
                       is_last_row=(i == last_idx))

    note_row = body_start + len(proteins)
    write_note_paragraph(ws, note_row, note_text, n_cols)
    ws.freeze_panes = f"A{body_start}"


# -----------------------------------------------------------------------------
# Caption and note text (markdown * = italic per md_to_rich)
# -----------------------------------------------------------------------------
CAP_1A = (
    "Table 1A. Top 20 quantitative proteins significantly upregulated in "
    "cortical astrocytes following ketamine treatment."
)
CAP_1B = (
    "Table 1B. All quantitative proteins significantly downregulated in "
    "cortical astrocytes following ketamine treatment."
)
CAP_1C = (
    "Table 1C. Proteins detected exclusively in ketamine-treated cortical "
    "astrocytes (presence/absence)."
)
CAP_S1 = (
    "Supplementary Table S1. Proteins detected exclusively in control cortical "
    "astrocytes (presence/absence)."
)
CAP_S2 = (
    "Supplementary Table S2. All 48 quantitative proteins significantly altered "
    "(adj. *p* ≤ 0.05) by ketamine treatment."
)

NOTE_QUANT_UP = (
    "*Note.* Top 20 quantitative proteins significantly upregulated in ketamine "
    "relative to control, ranked by absolute log₂ fold change. Adjusted "
    "*p*-values are from Welch's *t*-test on log₂-transformed abundance ratios "
    "as reported by Proteome Discoverer 3.1 (Tukey HSD post-hoc after ANOVA "
    "across pairwise replicate ratios). *n* = 3 biological replicates per "
    "group. Abbreviations: FC, fold change; Adj., adjusted. Pathway categories "
    "are investigator-assigned based on canonical KEGG and GO:BP membership; "
    "see Section [X] of the manuscript for full pathway enrichment results."
)

NOTE_QUANT_DOWN = (
    "*Note.* All four quantitative proteins significantly downregulated in "
    "ketamine relative to control. Adjusted *p*-values are from Welch's "
    "*t*-test on log₂-transformed abundance ratios as reported by Proteome "
    "Discoverer 3.1 (Tukey HSD post-hoc after ANOVA across pairwise replicate "
    "ratios). *n* = 3 biological replicates per group. Abbreviations: FC, "
    "fold change; Adj., adjusted."
)

NOTE_PA_KET = (
    "*Note.* All 10 proteins identified by high-confidence MS2 only in "
    "ketamine samples (presence/absence). These proteins received placeholder "
    "abundance ratios of 100 from Proteome Discoverer and are not amenable to "
    "standard differential abundance testing; absence from control samples "
    "may reflect either true biological absence or stochastic MS2 sampling "
    "under Data-Dependent Acquisition (DDA). *n* = 3 biological replicates "
    "per group. Abbreviations: MS2, tandem mass spectrometry fragmentation "
    "scan; DDA, Data-Dependent Acquisition; P/A, presence/absence."
)

NOTE_PA_CTRL = (
    "*Note.* Supplementary table. All 28 proteins identified by high-confidence "
    "MS2 only in control samples (presence/absence). These proteins received "
    "placeholder abundance ratios of 0.01 from Proteome Discoverer and are "
    "not amenable to standard differential abundance testing; absence from "
    "ketamine samples may reflect either true biological absence or stochastic "
    "MS2 sampling under Data-Dependent Acquisition (DDA). *n* = 3 biological "
    "replicates per group. Abbreviations: MS2, tandem mass spectrometry "
    "fragmentation scan; DDA, Data-Dependent Acquisition; P/A, presence/absence."
)

NOTE_SUPP_ALL = (
    "*Note.* Supplementary table. Complete list of all 48 quantitative proteins "
    "significantly altered (adj. *p* ≤ 0.05) by ketamine treatment, ranked by "
    "absolute log₂ fold change. Adjusted *p*-values from Welch's *t*-test on "
    "log₂-transformed abundance ratios (Proteome Discoverer 3.1; Tukey HSD "
    "post-hoc after ANOVA). *n* = 3 biological replicates per group. "
    "Abbreviations: FC, fold change; Adj., adjusted."
)

CAP_S3 = (
    "Supplementary Table S3. RT-qPCR relative transcript abundance for cortical "
    "cell-type markers in ACSA-2-positive and flow-through fractions."
)

CAP_S4 = (
    "Supplementary Table S4. High-confidence proteins excluded from the "
    "quantitative differential abundance analysis because Proteome Discoverer "
    "did not return an abundance ratio."
)

NOTE_SUPP_MR = (
    "*Note.* Supplementary table. All 19 proteins identified at high MS2 "
    "confidence that were excluded from the quantitative differential "
    "abundance analysis (846 quantitative + 38 presence/absence + 19 "
    "missing-ratio = 903 high-confidence proteins). These proteins fall into "
    "two patterns: (i) detected at high confidence in the control group only "
    "(*n* = 14), where Proteome Discoverer's pairwise ratio computation could "
    "not produce a quantitative estimate because either no peptide features "
    "were matched in the ketamine samples (n/a) or only a sub-confidence "
    "Peak Found signal was present; and (ii) detected at high confidence in "
    "both groups (*n* = 5) but with insufficient unique-peptide coverage for "
    "Proteome Discoverer to compute a stable ratio. The detection pattern "
    "column reports Proteome Discoverer's \"Found in Sample Group\" status "
    "verbatim. The Unique peptides column reflects the most likely driver of "
    "ratio failure: all 19 proteins were identified by a single unique "
    "peptide, which constrains Proteome Discoverer's ability to compute a "
    "ratio with sufficient precision under unique-and-razor-peptide "
    "quantification. These proteins were also excluded from the presence/absence "
    "category because Proteome Discoverer did not assign them the placeholder "
    "ratios (100 or 0.01) that define that category. *n* = 3 biological "
    "replicates per group. Abbreviations: MS2, tandem mass spectrometry "
    "fragmentation scan."
)

NOTE_SUPP_QPCR = (
    "*Note.* Supplementary table. Relative transcript abundance measured by "
    "RT-qPCR on a parallel cohort of cortical samples processed through the "
    "ACSA-2 magnetic isolation protocol. Values within each fraction are "
    "expressed relative to the astrocyte marker *Slc1a3* (set to 1.0). "
    "The final column reports the fold-depletion of each marker relative to "
    "*Slc1a3* in the ACSA-2-positive fraction (i.e., 1 / column 3); larger "
    "values indicate stronger depletion of the corresponding cell type from "
    "the astrocyte-enriched fraction. Source data: data/Adult ctx astro "
    "miltenyi isolation qPCR characterization Betti analysis.prism. "
    "*n* = 1 biological replicate per fraction; samples "
    "are not the same animals used for proteomics because RNA and protein "
    "yields from a single isolation were insufficient to support both assays. "
    "Marker assignments: *Slc1a3* (GLAST; astrocyte); *Olig2* (oligodendrocyte "
    "lineage); *Iba1*/*Aif1* (microglia); *Syt1* (neuron); *Cspg4* (NG2; "
    "oligodendrocyte precursor cells); *Fgfr4* (purpose to be confirmed). "
    "Abbreviations: ACSA-2, astrocyte cell surface antigen-2; RT-qPCR, reverse "
    "transcription quantitative PCR; OPC, oligodendrocyte precursor cell."
)


# qPCR values transcribed from the Prism source file
# (Adult ctx astro miltenyi isolation qPCR characterization Betti analysis.prism).
# Values are normalized within each fraction to Slc1a3 = 1.0.
QPCR_ROWS = [
    {"marker": "Slc1a3", "cell_type": "Astrocyte",
     "acsa2_pos": 1.0,         "flow_through": 1.0},
    {"marker": "Olig2",  "cell_type": "Oligodendrocyte lineage",
     "acsa2_pos": 0.005775623, "flow_through": 0.137537},
    {"marker": "Iba1",   "cell_type": "Microglia",
     "acsa2_pos": 0.000773161, "flow_through": 1.010132},
    {"marker": "Syt1",   "cell_type": "Neuron",
     "acsa2_pos": 0.000535469, "flow_through": 0.117886},
    {"marker": "Cspg4",  "cell_type": "Oligodendrocyte precursor cell",
     "acsa2_pos": 0.003111245, "flow_through": 0.247552},
    {"marker": "Fgfr4",  "cell_type": "Not yet confirmed",
     "acsa2_pos": 0.0000468059, "flow_through": 0.000934},
]


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def _normalize_pd_status(raw):
    if raw is None:
        return "n/a"
    s = raw.strip()
    return s if s else "n/a"


def load_missing_ratio_proteins():
    """Load the 19 missing-ratio proteins from the categorized output and
    enrich each row with PD's Found-in-Sample-Group status and unique-peptide
    count, both pulled from the raw PD CSVs."""
    categorized = load_csv(CATEGORIZED_PATH)
    missing = [r for r in categorized if r["category"] == "missing_ratio"]

    # Build an accession -> (status_ctrl, status_ket, unique_peptides) lookup
    # from both raw PD CSVs.
    lookup = {}
    for path in (RAW_CTRL_PATH, RAW_BOTH_PATH):
        for r in load_csv(path):
            acc = r["Accession"]
            lookup[acc] = {
                "status_ctrl": _normalize_pd_status(
                    r.get("Found in Sample Group: control")),
                "status_ket": _normalize_pd_status(
                    r.get("Found in Sample Group: ketamine")),
                "unique_peptides": r.get("# Unique Peptides", "").strip(),
            }

    enriched = []
    for r in missing:
        info = lookup.get(r["Accession"], {})
        status_ctrl = info.get("status_ctrl", "n/a")
        status_ket = info.get("status_ket", "n/a")
        upep_raw = info.get("unique_peptides", "")
        try:
            upep = int(float(upep_raw))
        except (TypeError, ValueError):
            upep = None
        r = dict(r)
        r["detection_pattern"] = (
            f"Control: {status_ctrl}; Ketamine: {status_ket}")
        r["unique_peptides"] = upep if upep is not None else "—"
        enriched.append(r)

    # Sort by detection pattern (both-high first), then by gene symbol.
    pattern_order = {
        "Control: High; Ketamine: High": 0,
        "Control: High; Ketamine: Peak Found": 1,
        "Control: High; Ketamine: n/a": 2,
    }
    enriched.sort(key=lambda x: (
        pattern_order.get(x["detection_pattern"], 99),
        x["Gene Symbol"].upper()))
    return enriched


def main():
    sig = load_csv(SIG_PATH)
    ket_pa = load_csv(KET_PA_PATH)
    ctrl_pa = load_csv(CTRL_PA_PATH)
    missing_ratio = load_missing_ratio_proteins()

    for p in sig:
        p["_abs"] = abs(to_float(p["log2_fold_change"]) or 0)
    sig_sorted = sorted(sig, key=lambda x: x["_abs"], reverse=True)

    up = [p for p in sig_sorted if p["direction"] == "up_in_ketamine"]
    down = [p for p in sig_sorted if p["direction"] == "down_in_ketamine"]

    top20_up = up[:20]
    all_4_down = down

    wb = Workbook()
    wb.remove(wb.active)

    ws = wb.create_sheet("Table 1A - Top 20 Upregulated")
    build_quant_sheet(ws, CAP_1A, top20_up, NOTE_QUANT_UP)

    ws = wb.create_sheet("Table 1B - Downregulated")
    build_quant_sheet(ws, CAP_1B, all_4_down, NOTE_QUANT_DOWN)

    ws = wb.create_sheet("Table 1C - Ketamine-specific")
    build_pa_sheet(ws, CAP_1C, ket_pa, "ketamine", NOTE_PA_KET)

    ws = wb.create_sheet("Supp S1 - Control-specific")
    build_pa_sheet(ws, CAP_S1, ctrl_pa, "control", NOTE_PA_CTRL)

    ws = wb.create_sheet("Supp S2 - All Sig Quant")
    build_quant_sheet(ws, CAP_S2, sig_sorted, NOTE_SUPP_ALL)

    ws = wb.create_sheet("Supp S3 - qPCR Cell-Type")
    build_qpcr_sheet(ws, CAP_S3, QPCR_ROWS, NOTE_SUPP_QPCR)

    ws = wb.create_sheet("Supp S4 - Missing-Ratio")
    build_mr_sheet(ws, CAP_S4, missing_ratio, NOTE_SUPP_MR)

    OUTPUT_PATH.parent.mkdir(exist_ok=True, parents=True)
    wb.save(OUTPUT_PATH)
    print(f"Wrote {OUTPUT_PATH}")
    print(f"  Sheets: {wb.sheetnames}")


if __name__ == "__main__":
    main()
