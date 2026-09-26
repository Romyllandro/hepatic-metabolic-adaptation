#!/usr/bin/env python3
"""Build the corrected, originally-formatted Supplementary Tables workbook.

Inputs
  ORIG.xlsx : Supplementary_Tables_ready.xlsx (original submission; house style template)
  NEW.xlsx  : Supplementary_Tables_ready_2026-09-24.xlsx (current content)
Output
  Supplementary_Tables_ready_2026-09-25.xlsx

Strategy
  * Load ORIG as the base so every sheet whose content did not change (S1, S2, SA-SJ, S9,
    Notes_RQ3, Notes_SJ) keeps its original cell styles byte-for-byte.
  * Rebuild every sheet whose shape changed (S3-S5, S3_Alt, S4_Alt, SP1-SP5, SK1-SK6, Contents,
    Notes_RQ1/RQ2/SK) from NEW values, applying the ORIG house style (fonts, fills, borders,
    number formats, row heights, semantic row colours) copied from ORIG template cells.
  * Apply the content corrections listed in CORRECTIONS below. No numeric value is changed
    except those explicitly listed.
"""
import copy, math, sys
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter as L

ORIG, NEW, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
wb = load_workbook(ORIG)
nw = load_workbook(NEW)

# ----------------------------------------------------------------------------- style kit
T = wb["S1_Convergent_All3Diets"]
HAIR = copy.copy(T.cell(3, 1).border)          # body cell border (hair)
THIN = copy.copy(T.cell(2, 2).border)          # header cell border (thin)
NOB = Border()


def fill(rgb):
    return PatternFill("solid", fgColor=rgb) if rgb else PatternFill()


def font(sz=8, b=False, i=False, color=None, name="Arial"):
    return Font(name=name, sz=sz, b=b, i=i, color=color)


def style(c, f=None, fl=None, al=None, bd=None, nf=None):
    if f is not None: c.font = f
    if fl is not None: c.fill = fl
    if al is not None: c.alignment = al
    if bd is not None: c.border = bd
    if nf is not None: c.number_format = nf


def lines(text, width_chars):
    if text is None: return 1
    n = 0
    for part in str(text).split("\n"):
        n += max(1, math.ceil(len(part) / max(width_chars, 1)))
    return n


def colw(ws, c):
    d = ws.column_dimensions[L(c)]
    return d.width if d.width else 8.43


def title_row(ws, text, ncol, fillc, h=28.15, row=1, sz=11):
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=ncol)
    c = ws.cell(row, 1, text)
    style(c, font(sz, b=True, color="FFFFFF"), fill(fillc), Alignment(horizontal="center", vertical="center", wrap_text=True))
    for k in range(2, ncol + 1):
        ws.cell(row, k).fill = fill(fillc)
    w = sum(colw(ws, k) for k in range(1, ncol + 1)) * 1.05
    ws.row_dimensions[row].height = max(h, 15.5 * lines(text, w) + 6)


def bar_row(ws, text, ncol, fillc, row):
    """In-sheet section bar (as ORIG SK3 'Part A — ...')."""
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=ncol)
    c = ws.cell(row, 1, text)
    style(c, font(9, b=True, color="FFFFFF"), fill(fillc), Alignment(horizontal="left", vertical="center", wrap_text=True))
    w = sum(colw(ws, k) for k in range(1, ncol + 1)) * 1.15
    ws.row_dimensions[row].height = max(18.0, 12.5 * lines(text, w) + 5)


def note_row(ws, text, ncol, row, sz=8):
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=ncol)
    c = ws.cell(row, 1, text)
    style(c, font(sz, i=True, color="595959"), None, Alignment(horizontal="left", vertical="top", wrap_text=True))
    w = sum(colw(ws, k) for k in range(1, ncol + 1)) * 1.25
    ws.row_dimensions[row].height = max(15.0, 11.5 * lines(text, w) + 4)


def header_row(ws, row, headers, fillc, sz=9, h=36.0):
    for k, v in enumerate(headers, 1):
        c = ws.cell(row, k, v)
        style(c, font(sz, b=True, color="FFFFFF"), fill(fillc), Alignment(horizontal="center", vertical="center", wrap_text=True), THIN)
    need = max(lines(v, colw(ws, k) * 1.1) for k, v in enumerate(headers, 1))
    ws.row_dimensions[row].height = max(h, 12.5 * need + 6)


def body_rows(ws, r0, rows, specs, rowfill, bold=None, color=None):
    """specs: list of (number_format, horizontal) per column."""
    for i, vals in enumerate(rows):
        r = r0 + i
        fc = rowfill(i, vals)
        b = bold(i, vals) if bold else False
        for k, v in enumerate(vals, 1):
            c = ws.cell(r, k, v)
            nf, hz = specs[k - 1]
            style(c, font(8, b=(b and k == 1), color=color), fill(fc), Alignment(horizontal=hz, vertical="center"), HAIR, nf)


def widths(ws, d):
    for k, v in d.items():
        ws.column_dimensions[k].width = v


def fresh(name, index=None, old=None):
    """Replace/insert a sheet called `name` (removing `old` or same-name sheet first)."""
    for n in {name, old} - {None}:
        if n in wb.sheetnames:
            del wb[n]
    return wb.create_sheet(name)


def grid(ws, r0, r1, c1):
    return [[ws.cell(r, c).value for c in range(1, c1 + 1)] for r in range(r0, r1 + 1)]


def yn(v):
    if v is True: return "Yes"
    if v is False: return "No"
    return v


TIER_FILL = lambda n: ("E2EFDA" if n == 9 else "BDD7EE" if n in (7, 8) else "FFF2CC" if n in (4, 5, 6) else "FFFFFF")
BAND = lambda i, v: "EBF3FB" if i % 2 == 0 else "FFFFFF"
NUM, D4, D2, SCI, CEN, LEFT, RIGHT = "General", "0.0000", "0.00", "0.00E+00", "center", "left", "right"
STRAINS = ["129S1SvImJ", "AJ", "C57BL6J", "CASTEiJ", "DBA2J", "NODShiLtJ", "NZOHlLtJ", "PWKPhJ", "WSBEiJ"]
LOG = []

# ============================================================================= RQ2: S3 / S4
L_HDR = {
    "L2_n_strains_tested_of8": "L2: N Strains Tested (of 8)",
    "L2_n_strains_FDR_lt_0.05_of8": "L2: N Strains FDR<0.05 (of 8)",
    "L2_min_FDR_of8": "L2: Min FDR (of 8)",
    "L2_NZOHlLtJ_FDR_excluded_nonestimable": "L2: NZO/HlLtJ FDR (non-estimable; excluded)",
    "L3_target_strain": "L3: Target Strain",
    "L3_agg_pvalue": "L3a: Aggregate p",
    "L3_agg_FDR": "L3a: Aggregate FDR",
    "L3_rep_pvalue": "L3b: Replicate p",
    "L3_rep_FDR": "L3b: Replicate FDR",
    "L3_supported_both_levels_FDR_lt_0.05": "L3: Supported at Both Levels (FDR<0.05)",
}
L_LEGEND = ("Statistical evidence layers. L2 (within-strain): Welch's t-test with BH-FDR within each strain on its own HFD vs. SCD "
            "replicate flux solutions (replicate-preserving analysis); counts are over the 8 estimable strains (NZO/HlLtJ excluded; HFD n = 1). "
            "L3 (between-strain specificity) for the target strain: L3a = one-sample t-test of the target strain's point estimate against the "
            "other eight strains' point estimates (df = 7; low power); L3b = HC3-robust OLS diet × strain-membership interaction on per-replicate "
            "data (recommended test). For strain-unique reactions the target is the single responsive strain; for reactions responsive in >1 "
            "strain it is the responsive strain with the strongest replicate-level evidence (sanity check, not a joint multi-strain test).")


def rq2_headers(hdr):
    out = []
    for h in hdr:
        h = L_HDR.get(h, h)
        if isinstance(h, str):
            h = h.replace("N Strains\nSignificant", "N Strains\nMagnitude-Responsive")
        out.append(h)
    return out


def rq2_specs(hdr):
    sp = []
    for h in hdr:
        s = str(h)
        if s.startswith("Δ"): sp.append((D4, LEFT))
        elif "FDR" in s and "N Strains" not in s and "Both" not in s or s.endswith(" p"): sp.append((SCI, LEFT))
        elif s in ("Reaction ID", "Reaction Name", "Subsystem", "L3: Target Strain"): sp.append((NUM, LEFT))
        else: sp.append((NUM, CEN))
    return sp


src = nw["S3_Conserved_Core_11"]
hdr = rq2_headers([src.cell(3, c).value for c in range(1, 25)])
rows = [[yn(v) for v in r] for r in grid(src, 4, 14, 24)]
ws = fresh("S3_Conserved_Core_11")
widths(ws, {"A": 12, "B": 38, "C": 30, **{L(k): 12 for k in range(4, 13)}, "M": 12, "N": 11,
            **{L(k): 13 for k in range(15, 25)}, "S": 16, "T": 14})
title_row(ws, "Table S3. Universal HFD-responsive core: reactions with |ΔHFD−SCD| ≥ 0.20 mmol/gDW/h in all 9 mouse strains "
              "(n = 11; standard FBA, aggregate design)", 24, "375623")
header_row(ws, 2, hdr, "375623")
body_rows(ws, 3, rows, rq2_specs(hdr), lambda i, v: "E2EFDA")
note_row(ws, "Consistent Direction = same sign of Δ in all 9 strains; EX_o2_e and O2t exceed the threshold in all 9 strains but reverse "
             "sign in C57BL6J. " + L_LEGEND, 24, 3 + len(rows) + 1)
LOG.append("S3: header 'N Strains Significant' -> 'Magnitude-Responsive'; L2/L3 headers made readable; True/False -> Yes/No")

src = nw["S4_All_Strains_Matrix"]
hdr = rq2_headers([src.cell(3, c).value for c in range(1, 25)])
rows = [[yn(v) for v in r] for r in grid(src, 4, 246, 24)]
ws = fresh("S4_All_Strains_Matrix")
widths(ws, {"A": 12, "B": 38, "C": 30, **{L(k): 12 for k in range(4, 13)}, "M": 11, "N": 12,
            **{L(k): 13 for k in range(15, 25)}, "S": 16, "T": 14})
title_row(ws, "Table S4. All HFD magnitude-responsive reactions across 9 mouse strains (|ΔHFD−SCD| ≥ 0.20 in ≥ 1 strain; union n = 243; "
              "standard FBA, aggregate design)", 24, "2E5090")
header_row(ws, 2, hdr, "2E5090")
body_rows(ws, 3, rows, rq2_specs(hdr), lambda i, v: TIER_FILL(v[12]))
note_row(ws, "Δ values are shown for all strains, including those below threshold. Row colour by number of magnitude-responsive strains: "
             "green = 9 (universal core, Table S3); blue = 7–8; yellow = 4–6; white = 1–3. " + L_LEGEND, 24, 3 + len(rows) + 1)
LOG.append("S4: same header/boolean clean-up as S3; conservation-tier row colours restored")

# ----------------------------------------------------------------------------- S5
src = nw["S5_PerStrain_Summary"]
rows = grid(src, 13, 21, 8)                      # 9 strains; stale duplicate row 22 dropped
assert [r[0] for r in rows] == STRAINS
assert src.cell(22, 1).value == "WSBEiJ" and src.cell(22, 7).value == 15   # the stale replicate-design duplicate
SUPPORT = {"129S1SvImJ": 9, "AJ": 4, "C57BL6J": 3, "CASTEiJ": 14, "DBA2J": 0, "NODShiLtJ": 6, "NZOHlLtJ": "NE", "PWKPhJ": 4, "WSBEiJ": 25}
for r in rows:
    old = r[3]
    r[3] = SUPPORT[r[0]]
    if old != r[3]:
        LOG.append(f"S5: {r[0]} N statistically supported {old} -> {r[3]}")
LOG.append("S5: deleted stale duplicate WSBEiJ row (54 | 1 | 25 | ... | 15, replicate-design values); removed 9 blank rows above header")
hdr = ["Strain", "N Magnitude-Responsive\n(|Δ| ≥ 0.20)", "N Strain-Unique", "N Statistically Supported\n(within-strain FDR<0.05)",
       "Min Δ\n(magnitude-responsive)", "Max Δ\n(magnitude-responsive)", "N Universal Core\nReactions (of 11)", "Within-Strain Inference\nEstimable"]
ws = fresh("S5_PerStrain_Summary")
widths(ws, {"A": 20, "B": 16, "C": 12, "D": 20, "E": 14, "F": 14, "G": 14, "H": 16})
title_row(ws, "Table S5. Per-strain summary of HFD responses (RQ2; HFD vs. SCD)", 8, "1F3864")
header_row(ws, 2, hdr, "1F3864", h=48)
body_rows(ws, 3, rows, [(NUM, LEFT), (NUM, CEN), (NUM, CEN), (NUM, CEN), (D4, LEFT), (D4, LEFT), (NUM, CEN), (NUM, CEN)], BAND)
note_row(ws, "Magnitude columns: aggregate-design standard FBA (Tables S3–S4). N Statistically Supported = number of all 3,726 reactions "
             "reaching Welch's t-test BH-FDR < 0.05 within the strain on its own replicate flux solutions (replicate-preserving analysis; "
             "range 0–25 over the 8 estimable strains). NE = non-estimable (NZO/HlLtJ has HFD n = 1).", 8, 13)

# ----------------------------------------------------------------------------- S3_Alt / S4_Alt
src = nw["S3_Alt_ReplicateDesign_Core15"]
hdr = rq2_headers([src.cell(3, c).value for c in range(1, 15)])
rows = grid(src, 4, 18, 14)
ws = fresh("S3_Alt_ReplicateDesign_Core15")
widths(ws, {"A": 12, "B": 38, "C": 30, **{L(k): 12 for k in range(4, 13)}, "M": 12, "N": 11})
title_row(ws, "Table S3_Alt. Complementary replicate-preserving analysis: universal HFD-responsive core across all 9 mouse strains "
              "(|ΔHFD−SCD| ≥ 0.20; n = 15; standard FBA)", 14, "375623")
header_row(ws, 2, hdr, "375623")
body_rows(ws, 3, rows, rq2_specs(hdr), lambda i, v: "E2EFDA")
note_row(ws, "Per-animal expression profiles were modeled independently (replicate-preserving design) and Δ is the difference of replicate "
             "means. All 11 reactions of the aggregate-design core (Table S3) are contained in this 15-reaction core.", 14, 19)

src = nw["S4_Alt_ReplicateDesign_Union187"]
hdr = rq2_headers([src.cell(3, c).value for c in range(1, 15)])
rows = grid(src, 4, 190, 14)
ws = fresh("S4_Alt_ReplicateDesign_Union187")
widths(ws, {"A": 12, "B": 38, "C": 30, **{L(k): 12 for k in range(4, 13)}, "M": 11, "N": 12})
title_row(ws, "Table S4_Alt. Complementary replicate-preserving analysis: all HFD magnitude-responsive reactions across 9 mouse strains "
              "(union n = 187; standard FBA)", 14, "2E5090")
header_row(ws, 2, hdr, "2E5090")
body_rows(ws, 3, rows, rq2_specs(hdr), lambda i, v: TIER_FILL(v[12]))
note_row(ws, "Blank Δ = below the 0.20 threshold in that strain. Row colour by number of magnitude-responsive strains: green = 9 "
             "(Table S3_Alt); blue = 7–8; yellow = 4–6; white = 1–3. 156 of the 243 aggregate-design union reactions (64.2%) are in "
             "this union.", 14, 3 + len(rows) + 1)
LOG.append("S3_Alt/S4_Alt: removed '[SECONDARY EVIDENCE LAYER]', dates, meeting/advisor references, folder names and project-note paths from titles")

# ============================================================================= SP1-SP5 (pFBA sensitivity; were SL1-SL5)
SUB_FILL = "F2F2F2"
src = nw["SL1_pFBA_Convergent27"]
hdr = [src.cell(2, c).value for c in range(1, 24)]
rows = grid(src, 3, 29, 23)
ws = fresh("SP1_pFBA_Convergent27", old="SL1_pFBA_Convergent27")
widths(ws, {"A": 12, "B": 36, "C": 26, **{L(k): 10 for k in range(4, 24)}, "P": 11, "Q": 12})
title_row(ws, "Table SP1. pFBA sensitivity analysis (RQ1): reactions significant in HFD, KD, and WD vs. SCD under replicate-level pFBA "
              "(n = 27)", 23, "1F3864")
sp = []
for h in hdr:
    s = str(h)
    sp.append((D4, LEFT) if s.startswith("Adjusted Δ") else (D2, LEFT) if s.startswith("Cohen") else (SCI, LEFT) if s.startswith("q ")
              else (NUM, LEFT) if s in ("Reaction ID", "Reaction Name", "Subsystem") else (NUM, CEN))
header_row(ws, 2, hdr, "2E5090")
body_rows(ws, 3, rows, sp, lambda i, v: SUB_FILL if v[16] == "Yes" else "E2EFDA")
note_row(ws, "Sensitivity analysis for Table S1 (Results 3.1; SM10.6). Adjusted Δ from the dataset-adjusted model (Methods 2.3, Eq. 5; "
             "HC3 standard errors). Grey rows = R-group bookkeeping/pseudo-reactions (8 of 27). pFVA columns = number of matched cohorts "
             "in which the parsimonious FVA interval keeps the same direction at 1% / 5% tolerance.", 23, 31)

src = nw["SL2_pFBA_AllSignificant228"]
hdr = [src.cell(2, c).value for c in range(1, 26)]
rows = [r for r in grid(src, 3, src.max_row, 25) if r[0] is not None]
assert len(rows) == 228
ws = fresh("SP2_pFBA_AllSignificant228", old="SL2_pFBA_AllSignificant228")
widths(ws, {"A": 12, "B": 36, "C": 28, **{L(k): 9.5 for k in range(4, 26)}, "X": 13, "Y": 13})
title_row(ws, "Table SP2. pFBA sensitivity analysis (RQ1): all reactions significant in ≥ 1 SCD-referenced contrast under replicate-level "
              "pFBA (n = 228)", 25, "2E5090")
sp = []
for h in hdr:
    s = str(h)
    sp.append((D4, LEFT) if s.startswith("Δ") else (D2, LEFT) if s.startswith("d ") else (SCI, LEFT) if s.startswith("q ")
              else (NUM, LEFT) if s in ("Reaction ID", "Reaction Name", "Subsystem") else (NUM, CEN))
header_row(ws, 2, hdr, "2E5090")
body_rows(ws, 3, rows, sp, lambda i, v: {3: "E2EFDA", 2: "FFF2CC"}.get(v[23], "FFFFFF"))
note_row(ws, "Sensitivity analysis for Table S2 (Results 3.1; SM10.6). Per-contrast counts: HFD–SCD 61, KD–SCD 50, WD–SCD 199, HFD–KD 0, "
             "HFD–WD 155. Row colour: green = significant in all 3 SCD-referenced contrasts (n = 27, Table SP1); yellow = 2; white = 1.",
         25, 3 + len(rows) + 1)

src = nw["SL3_pFBA_UniversalCore16"]
hdr = [src.cell(2, c).value for c in range(1, 18)]
rows = [r for r in grid(src, 3, src.max_row, 17) if r[0] is not None]
assert len(rows) == 16
ws = fresh("SP3_pFBA_UniversalCore16", old="SL3_pFBA_UniversalCore16")
widths(ws, {"A": 12, "B": 38, "C": 28, **{L(k): 12 for k in range(4, 13)}, **{L(k): 13 for k in range(13, 18)}})
title_row(ws, "Table SP3. pFBA sensitivity analysis (RQ2): universal HFD magnitude-responsive core across all 9 mouse strains "
              "(|ΔHFD−SCD| ≥ 0.20; n = 16)", 17, "375623")
sp = [(NUM, LEFT)] * 3 + [(D4, LEFT)] * 9 + [(NUM, CEN)] * 5
header_row(ws, 2, hdr, "375623")
body_rows(ws, 3, rows, sp, lambda i, v: "E2EFDA")
note_row(ws, "Sensitivity analysis for Table S3 (Results 3.1; SM10.6). pFVA columns = number of strains whose targeted parsimonious FVA "
             "interval keeps the same direction at 1% / 5% tolerance (14/16 and 7/16 reactions retained in ≥ 1 strain-consistent "
             "interval, as reported in SM10.6).", 17, 3 + len(rows) + 1)

src = nw["SL4_pFBA_AllStrainsMatrix105"]
hdr = [src.cell(2, c).value for c in range(1, 18)]
rows = [r for r in grid(src, 3, src.max_row, 17) if r[0] is not None]
assert len(rows) == 105
ws = fresh("SP4_pFBA_AllStrainsMatrix105", old="SL4_pFBA_AllStrainsMatrix105")
widths(ws, {"A": 12, "B": 38, "C": 28, **{L(k): 12 for k in range(4, 13)}, **{L(k): 13 for k in range(13, 18)}})
title_row(ws, "Table SP4. pFBA sensitivity analysis (RQ2): all HFD magnitude-responsive reactions across 9 mouse strains "
              "(|ΔHFD−SCD| ≥ 0.20 in ≥ 1 strain; union n = 105)", 17, "2E5090")
header_row(ws, 2, hdr, "2E5090")
body_rows(ws, 3, rows, sp, lambda i, v: TIER_FILL(v[12]))
note_row(ws, "Sensitivity analysis for Table S4 (Results 3.1; SM10.6). Row colour by number of magnitude-responsive strains: green = 9; "
             "blue = 7–8; yellow = 4–6; white = 1–3.", 17, 3 + len(rows) + 1)
LOG.append("SL1-SL4 renamed SP1-SP4 (avoids collision with Supplement Tables SL1/SL2/SL3a/SL3b, the benchmark tables); "
           "'[SENSITIVITY ANALYSIS ...]' / 'revised' wording removed; 137 empty padding rows dropped from SL4 and 222 from SL2")

# ----------------------------------------------------------------------------- SP5
src = nw["SL5_Sensitivity_Summary"]
S = {r: [src.cell(r, c).value for c in range(1, 6)] for r in range(1, src.max_row + 1)}
ws = fresh("SP5_Sensitivity_Summary", old="SL5_Sensitivity_Summary")
widths(ws, {"A": 34, "B": 30, "C": 22, "D": 18, "E": 26})
title_row(ws, "Table SP5. Consolidated RQ1/RQ2 sensitivity and robustness analyses (FBA vs. pFBA, E-Flux cap, primary objective, "
              "FVA/pFVA, and RQ2 magnitude-threshold θ sweep). These characterize solution-selection and threshold sensitivity of the "
              "primary standard-FBA results in Tables S1–S5; they do not replace them (Results 3.1, 3.3; Methods 2.4; SM10.6).", 5, "1F3864")
r = 3


def sp5_block(title, header, data, specs):
    global r
    bar_row(ws, title, 5, "2E5090", r); r += 1
    header_row(ws, r, header + [None] * (5 - len(header)), "2E5090", sz=8.5, h=30)
    for k in range(len(header) + 1, 6):
        ws.cell(r, k).fill = PatternFill(); ws.cell(r, k).border = NOB
    r += 1
    body_rows(ws, r, data, specs, BAND); r += len(data) + 1


sp5_block(S[3][0], S[4], [S[5], S[6]], [(NUM, LEFT), ("0.000", LEFT), ("0.0", LEFT), ("0.000", LEFT), ("0.0", LEFT)])
bar_row(ws, S[8][0], 5, "2E5090", r); r += 1
body_rows(ws, r, [[S[9][0], S[9][1], None, None, None]], [(NUM, LEFT)] * 5, BAND)
ws.merge_cells(start_row=r, start_column=2, end_row=r, end_column=5)
ws.cell(r, 1).alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
ws.row_dimensions[r].height = 24; r += 2
sp5_block(S[11][0], S[12][:3], [S[13][:3], S[14][:3]], [(NUM, LEFT), ("0.000", LEFT), ("0.0", LEFT)])
sp5_block(S[16][0], S[17][:4], [S[18][:4], S[19][:4], S[20][:4]], [(NUM, LEFT), (NUM, CEN), (NUM, CEN), (NUM, CEN)])
assert S[22][0].startswith("E. Targeted parsimonious FVA")
sp5_block("E. Targeted parsimonious FVA (1% tolerance) applied to the 15-reaction replicate-preserving universal core (Table S3_Alt; "
          "includes all 11 aggregate-design core reactions of Table S3) — reaction-level alternate-optima check",
          S[23][:4], [S[k][:4] for k in range(24, 39)], [(NUM, LEFT), (NUM, CEN), (NUM, CEN), (NUM, CEN)])
note_row(ws, "Targeted pFVA values were computed in the pFBA sensitivity run with the same targeted-pFVA procedure used for the 16-reaction "
             "pFBA universal core (Table SP3), applied here to the 15-reaction replicate-preserving core. Alternate-optimum uncertainty near "
             "the pFBA solution persists for this reaction set (0–2 of 9 strains FVA-robust), as for the 27-reaction pFBA RQ1 signature and "
             "the 16-reaction pFBA RQ2 core. These values characterize the solution space of these reactions; they do not re-derive the "
             "standard-FBA results, which involve no parsimony minimization.", 5, r); r += 2
assert S[42][0].startswith("F. RQ2 replicate-level statistical support")
Fdata = [S[k][:2] for k in range(44, 53)]
for row in Fdata:
    if row[0] == "NZOHlLtJ":
        row[1] = "NE (HFD n = 1)"
sp5_block("F. RQ2 within-strain replicate-level statistical support (replicate-preserving standard-FBA analysis; Welch's t-test, BH-FDR "
          "within strain, all 3,726 reactions)", ["Strain", "N reactions FDR_BH < 0.05"], Fdata, [(NUM, LEFT), (NUM, CEN)])
sp5_block("G. RQ2 magnitude-threshold (θ) sensitivity sweep: universal-core and union size across θ = 0.10–0.30 (standard-FBA, "
          "aggregate design; Eq. 7; recomputed from the same per-strain, per-reaction |ΔHFD−SCD| flux differences underlying Tables "
          "S3–S4, at θ ≠ 0.20)",
          ["θ", "Union n (≥ 1/9 strains)", "Universal-core n (9/9 strains)"],
          [[0.10, 421, 14], [0.15, 264, 13], [0.20, 243, 11], [0.25, 205, 9], [0.30, 181, 4]],
          [(D2, LEFT), (NUM, CEN), (NUM, CEN)])
note_row(ws, "Section G reports the primary θ = 0.20 result (243/11; identical to Tables S3–S4) alongside four additional θ values "
             "requested during review (Results 3.3; Methods 2.4). The universal core ranges from 4 to 14 reactions across this range, "
             "and its exact membership is threshold-dependent even though a compact conserved core is robustly recovered at every "
             "tested θ; the primary-θ reaction identities are given in Table S3.", 5, r); r += 1
LOG.append("SP5 (was SL5): section E relabelled to the 15-reaction replicate-preserving core (was 'new ... RQ2 primary result'); "
           "note rewritten; section F '(new ...)' removed; NZO/HlLtJ 0 -> 'NE (HFD n = 1)'; section G added "
           "(2026-09-25): RQ2 theta-sweep (0.10/0.15/0.20/0.25/0.30 -> union 421/264/243/205/181, core 14/13/11/9/4), "
           "synchronizing the Supplementary Tables with the Response to Reviewers letter's theta-sweep claim (Results 3.3; Methods 2.4)")

# ============================================================================= SK1-SK6 (RQ4)
RQ4REF = "Methods 2.4–2.5; Results 3.6; SM8"
src = nw["SK1_Species_Composition"]
rows = grid(src, 3, 59, 8)
for rr in rows:
    rr[5] = {"=TRUE()": "Yes", "=FALSE()": "No"}[rr[5]]
    rr[4] = {"matched": "Matched", "unmatched": "Unmatched"}.get(rr[4], rr[4])
LOG.append("SK1: 57 =TRUE()/=FALSE() formulas converted to Yes/No text (workbook now formula-free)")
hdr = ["Species (metatranscriptome)", "ND/SCD\nAbundance", "DD/HFD\nAbundance", "AGORA2 Model\n(if matched)", "AGORA2\nMatch Status",
       "In Primary\n27-Species Map", "Direction", "Fold Change\n(DD/ND)"]
ws = fresh("SK1_Species_Composition")
widths(ws, {"A": 46, "B": 12, "C": 12, "D": 50, "E": 12, "F": 13, "G": 12, "H": 12})
title_row(ws, "Table SK1. Microbial community composition shifts — ND/SCD vs. DD/HFD (all 57 detected species)", 8, "1F3864", h=26.1)
note_row(ws, "Note: SK1 lists all 57 species detected in the metatranscriptomic data (GSE104913). 27 species (bold) were retained in the "
             "strict primary AGORA2 map used for MICOM community modeling (26 exact mappings plus one verified Bacteroides fragilis strain "
             "synonym; 67.1% of ND/SCD and 71.4% of DD/HFD activity-weighted abundance before renormalization); the remaining species "
             "were not part of the modeled community. Fold change is blank where ND/SCD abundance is zero.", 8, 2, sz=9)
header_row(ws, 3, hdr, "1F3864", sz=8.5, h=43.5)
DIRF = {"Decreased": "BDD7EE", "Increased": "FFC7CE", "Unchanged": "F2F2F2"}
body_rows(ws, 4, rows, [(NUM, LEFT), ("0.000000", RIGHT), ("0.000000", RIGHT), (NUM, LEFT), (NUM, CEN), (NUM, CEN), (NUM, LEFT), ("0.000", RIGHT)],
          lambda i, v: DIRF[v[6]], bold=lambda i, v: v[5] == "Yes", color="000000")
from collections import Counter
dc = Counter(rr[6] for rr in rows)
ms = sum(rr[4] == "Matched" for rr in rows)
pm = sum(rr[5] == "Yes" for rr in rows)
r = 4 + len(rows) + 1
c = ws.cell(r, 1, "Summary statistics")
style(c, font(9, b=True, color="FFFFFF"), fill("1F3864"), Alignment(horizontal="left", vertical="center"))
ws.cell(r, 2).fill = fill("1F3864")
for lab, val in [("Total species detected", len(rows)), ("Species decreased under HFD", dc["Decreased"]),
                 ("Species increased under HFD", dc["Increased"]), ("Species unchanged", dc["Unchanged"]),
                 ("Species with an AGORA2 match", ms), ("Species retained in primary 27-species map", pm)]:
    r += 1
    a, b = ws.cell(r, 1, lab), ws.cell(r, 2, val)
    style(a, font(8, b=True), None, Alignment(horizontal="left", vertical="center"), HAIR)
    style(b, font(8), None, Alignment(horizontal="left", vertical="center"), HAIR)
assert (dc["Decreased"], dc["Increased"], dc["Unchanged"], pm) == (31, 19, 7, 27)

src = nw["SK2_Portal_Metabolites"]
rows = grid(src, 3, 7, 12)
hdr = ["Metabolite", "Metabolite ID", "DD/HFD Flux\n(mmol/gDW/h)", "ND/SCD Flux\n(mmol/gDW/h)", "Importance", "Category",
       "Portal Status\n(ND/SCD)", "Direction\n(ND/SCD)", "Portal Status\n(DD/HFD)", "Direction\n(DD/HFD)", "Fold Change\n(magnitude)",
       "Propagated to Hepatic Model\nas Host Delivery?"]
ws = fresh("SK2_Portal_Metabolites")
widths(ws, {"A": 30, "B": 12, "C": 14, "D": 14, "E": 11, "F": 11, "G": 12, "H": 11, "I": 12, "J": 11, "K": 12, "L": 44})
title_row(ws, "Table SK2. Portal metabolite production — microbiome-derived inputs to the hepatic model (primary configuration)", 12, "375623", h=26.1)
header_row(ws, 2, hdr, "375623", sz=8.5)
body_rows(ws, 3, rows, [(NUM, LEFT), (NUM, LEFT), ("0.000000", RIGHT), ("0.000000", RIGHT)] + [(NUM, LEFT)] * 6 + [("0.000", RIGHT), (NUM, LEFT)],
          lambda i, v: "FFC7CE" if v[4] == "high" else "F2F2F2", bold=lambda i, v: v[4] == "high", color="000000")
note_row(ws, "Primary configuration: 27-species map, MICOM cooperative tradeoff 0.50, community pFBA. Positive community exchange flux = net "
             "microbial export (candidate portal delivery); negative = net microbial uptake. Acetate is the only curated metabolite with a "
             "positive export above threshold (0.9652 → 0.1994 mmol/gDW/h, 79.3% decrease under DD/HFD) and is therefore the only metabolite "
             "propagated to the hepatic model; ammonium and the B vitamins reflect net microbial uptake and were not treated as host "
             "delivery. Red rows = high-importance metabolites.", 12, 9)

src = nw["SK3_Flux_Attribution"]
A = grid(src, 4, 9, 6)
DESC = {"Diet": "Diet is primary driver (VE_diet > VE_micro)", "Diet (opposed by microbiome)": "Diet drives change, microbiome counteracts",
        "Microbiome": "Microbiome is primary driver", "Microbiome (opposed by diet)": "Microbiome drives change, diet counteracts",
        "Stable": "No measurable flux change", "TOTAL / Global": None}
A = [row + [DESC[row[0]]] for row in A]
B = grid(src, 13, 319, 15)
hdrB = ["Reaction ID", "Reaction Name", "Subsystem", "Compartment", "SCD Baseline\nFlux", "HFD (no microbiome)\nFlux",
        "HFD (with microbiome)\nFlux", "Δ Diet", "Δ Microbiome", "Δ Total", "VE Diet\n(%)", "VE Micro\n(%)", "Effect Type",
        "Dominant Driver", "Microbiome\nEffect Size"]
ws = fresh("SK3_Flux_Attribution")
widths(ws, {"A": 26, "B": 40, "C": 26, "D": 14, **{L(k): 13 for k in range(5, 13)}, "M": 13, "N": 26, "O": 12})
title_row(ws, "Table SK3. Reaction-level diet–microbiome flux attribution (all 3,726 reactions; primary configuration)", 15, "2E5090", h=26.1)
bar_row(ws, "Part A — Driver classification summary (47 reactions with |Δmicrobiome| ≥ 0.01; 29 microbiome-dominant = 27 + 2 opposed; "
            "3,322 stable)", 15, "2E5090", 2)
header_row(ws, 3, ["Driver Category", "N Reactions", "% of 3,726", "% of Total Abs\nFlux Change", "Mean VE\nDiet (%)", "Mean VE\nMicro (%)",
                   "Description"], "2E5090", sz=8.5, h=43.5)
AF = {"Diet": ("FFF2CC", True), "Diet (opposed by microbiome)": ("FFEB9C", False), "Microbiome": ("BDD7EE", True),
      "Microbiome (opposed by diet)": ("DEEBF7", False), "Stable": ("F2F2F2", False), "TOTAL / Global": (None, True)}
for i, row in enumerate(A):
    rr = 4 + i if row[0] != "TOTAL / Global" else 10
    fc, b = AF[row[0]]
    for k, v in enumerate(row, 1):
        c = ws.cell(rr, k, v)
        style(c, font(8, b=b and k == 1, color="000000"), fill(fc), Alignment(horizontal="left" if k in (1, 7) else "right", vertical="center"),
              HAIR, "0" if k == 2 else "0.00" if 3 <= k <= 6 else NUM)
ws.merge_cells("G4:O4"); ws.merge_cells("G5:O5"); ws.merge_cells("G6:O6"); ws.merge_cells("G7:O7"); ws.merge_cells("G8:O8")
bar_row(ws, "Part B — 307 reactions with measurable microbiome contribution (VE_micro > 0; 8.2% of 3,726), sorted by absolute microbiome "
            "contribution (descending). The full 3,726-reaction attribution table is available with the code and data repository.", 15, "2E5090", 12)
header_row(ws, 13, hdrB, "2E5090", sz=8.5, h=32.65)
BF = {"Diet": "E2EFDA", "Microbiome": "E2EFDA", "Diet (opposed by microbiome)": "FCE4D6", "Microbiome (opposed by diet)": "FCE4D6", "Stable": "F2F2F2"}
body_rows(ws, 14, B, [(NUM, LEFT)] * 4 + [("0.000000", RIGHT)] * 6 + [("0.00", RIGHT)] * 2 + [(NUM, LEFT)] * 3, lambda i, v: BF[v[13]], color="000000")
assert [a[1] for a in A] == [260, 115, 27, 2, 3322, 3726] and len(B) == 307
assert sum(1 for b in B if b[8] is not None and abs(b[8]) >= 0.01) == 47

src = nw["SK4_Pathway_Synergy"]
rows = grid(src, 3, 41, 7)
ws = fresh("SK4_Pathway_Synergy")
widths(ws, {"A": 36, "B": 12, "C": 14, "D": 13, "E": 12, "F": 12, "G": 24})
title_row(ws, "Table SK4. Pathway-level diet–microbiome interaction classification (primary configuration; 39 subsystems with ≥ 2 "
              "reactions with |Δtotal| > 10⁻³)", 7, "7B2C8B", h=26.1)
header_row(ws, 2, ["Pathway (Subsystem)", "N Total\nReactions", "N Synergistic", "N Antagonistic", "Synergistic\n(%)", "Antagonistic\n(%)",
                   "Pathway Class"], "7B2C8B", sz=8.5)
CF = {"Highly Synergistic": "C6EFCE", "Moderately Synergistic": "E2EFDA", "Mixed": "FFEB9C", "Highly Antagonistic": "FFC7CE"}
body_rows(ws, 3, rows, [(NUM, LEFT), ("0", RIGHT), ("0", RIGHT), ("0", RIGHT), ("0.0", RIGHT), ("0.0", RIGHT), (NUM, LEFT)],
          lambda i, v: CF[v[6]], color="000000")
r = 3 + len(rows) + 1
ws.cell(r, 1, "Pathway class key:").font = font(8, b=True)
for k, (lab, fc) in enumerate(CF.items(), 2):
    c = ws.cell(r, k, lab)
    style(c, font(8), fill(fc), Alignment(horizontal="center", vertical="center", wrap_text=True), HAIR)
ws.row_dimensions[r].height = 24
note_row(ws, "Classification: > 70% synergistic = Highly Synergistic; > 70% antagonistic = Highly Antagonistic; > 50% synergistic = "
             "Moderately Synergistic; otherwise Mixed. Labels are descriptive summaries of the sign relationship between diet- and "
             "microbiome-attributable components, not validated mechanisms (SM8.3).", 7, r + 1)

src = nw["SK5_Pathway_Enrichment"]
rows = grid(src, 3, 57, 13)
ws = fresh("SK5_Pathway_Enrichment")
widths(ws, {"A": 34, "B": 10, "C": 12, "D": 10, "E": 12, "F": 12, "G": 14, "H": 16, "I": 13, "J": 13, "K": 13, "L": 12, "M": 13})
title_row(ws, "Table SK5. Pathway enrichment for diet–microbiome interactions (hypergeometric test, BH-adjusted; 55 subsystems tested, "
              "10 significant at q < 0.05)", 13, "375623", h=26.1)
header_row(ws, 2, ["Subsystem", "N Total", "N Significant", "Expected", "Enrichment\nRatio", "p-value", "Diet-dominated\n(%)",
                   "Microbiome-dominated\n(%)", "Synergistic\n(%)", "Mean Diet\nVariance", "Mean Micro\nVariance", "BH q-value",
                   "Significant\n(q<0.05)"], "375623", sz=8.5)
body_rows(ws, 3, rows, [(NUM, LEFT), ("0", RIGHT), ("0", RIGHT), ("0.00", RIGHT), ("0.000", RIGHT), (SCI, RIGHT), ("0.0", RIGHT),
                        ("0.0", RIGHT), ("0.0", RIGHT), ("0.000", RIGHT), ("0.000", RIGHT), (SCI, RIGHT), (NUM, RIGHT)],
          lambda i, v: "C6EFCE" if str(v[12]).startswith("Yes") else "FFFFFF", color="000000")
assert sum(str(v[12]).startswith("Yes") for v in rows) == 10
note_row(ws, "Background = all 3,726 iMM1415 reactions; 'significant' = reactions with |Δtotal| > 10⁻³ (n = 404). ✱ / green = BH q < 0.05.",
         13, 3 + len(rows) + 1)

src = nw["SK6_Compartment_Enrichment"]
rows = grid(src, 3, 7, 12)
ws = fresh("SK6_Compartment_Enrichment")
widths(ws, {"A": 16, "B": 12, "C": 12, "D": 10, "E": 12, "F": 14, "G": 14, "H": 14, "I": 14, "J": 14, "K": 12, "L": 13})
title_row(ws, "Table SK6. Subcellular compartment enrichment for microbiome-associated flux effects (primary configuration; 5 compartments "
              "tested, 1 significant at q < 0.05)", 12, "1F3864", h=26.1)
header_row(ws, 2, ["Compartment", "N Total\nReactions", "N Significant", "Expected", "Enrichment\nRatio", "p-value\n(hypergeometric)",
                   "Mean Micro\nContribution", "Mean Diet\nContribution", "Mean Micro\nVariance (%)", "Large Micro\nEffects (%)", "BH q-value",
                   "Significant\n(q<0.05)"], "1F3864", sz=8.5)
body_rows(ws, 3, rows, [(NUM, LEFT), ("0", RIGHT), ("0", RIGHT), ("0.00", RIGHT), ("0.000", RIGHT), (SCI, RIGHT), ("0.0000", RIGHT),
                        ("0.0000", RIGHT), ("0.00", RIGHT), ("0.00", RIGHT), (SCI, RIGHT), (NUM, RIGHT)],
          lambda i, v: "C6EFCE" if str(v[11]).startswith("Yes") else "F2F2F2", color="000000")
note_row(ws, "Same hypergeometric framework and 'significant' definition as Table SK5, applied to subcellular compartments.", 12, 9)
LOG.append("SK1-SK6: original house style restored (were 1F4E78, borderless, frozen panes); snake_case headers made readable; "
           "'accepted manuscript', 'v3-corrected', 'frozen', 'Section 2.6' wording removed; SK3 Part A descriptions restored")

# ============================================================================= unchanged-data sheets: text edits only
ws = wb["S1_Convergent_All3Diets"]
ws["A1"] = ("Table S1. Reactions significantly altered in all three high-caloric diets vs. study-specific control (HFD, KD, and WD vs. "
            "SCD; convergent metabolic signature; n = 10)")
ws = wb["S2_All_Significant_Reactions"]
ws["A1"] = "Table S2. All significantly altered reactions across any estimable pairwise dietary comparison (RQ1; n = 450)"
LOG.append("S1/S2: 'Table S.' numbered as Table S1/S2 (data and styling identical to the original)")

EDITS = {
    "SA_Contribution_Summary": {"A1": "Table S-A. Cell-type abundance-weighted contribution to diet-induced metabolic flux changes (WD vs. chow; "
                                      "P99/biomass reference) — Results 3.4",
                                "E2": "N Threshold-\nResponsive Reactions"},
    "SB_CrossStrain_Matrix": {"A1": "Table S-B. Cell-type contribution (%) projected across 9 mouse strain contexts — heuristic reweighting of "
                                    "the single-atlas attribution template (not strain-specific single-cell data)"},
    "SC_BulkCellular_Overlap": {"E2": "N Cell Threshold-\nResponsive"},
    "SD_Conservation_Analysis": {"A1": "Table S-D. Cross-strain variability of projected cell-type contributions — coefficient of variation (CV) analysis"},
    "SE_SigReactions_CellType": {"A1": "Table S-E. Threshold-responsive reactions per cell type (WD vs. chow; ≥ 2 of 3 deterministic criteria)",
                                 "D2": "N Threshold-\nResponsive Reactions"},
    "SH_Reaction_Attribution": {"A1": "Table S-H. Reaction-level cell-type attribution of the 205 bulk reference reactions (primary driver and "
                                      "attribution category, RQ3)"},
}
for sh, cells in EDITS.items():
    for ref, val in cells.items():
        old = wb[sh][ref].value
        wb[sh][ref] = val
        LOG.append(f"{sh}!{ref}: {old!r} -> {val!r}")
wb["SE_SigReactions_CellType"].title = "SE_Responsive_Rxns_CellType"
LOG.append("SE sheet renamed SE_SigReactions_CellType -> SE_Responsive_Rxns_CellType (threshold-responsive, not significance-tested)")

# ----------------------------------------------------------------------------- S9 (original sheet; tier fix + note)
ws = wb["S9_CrossScale_Integration"]
assert [ws.cell(r, 6).value for r in (6, 10, 14)] == ["Moderate"] * 3   # ORIG already correct; NEW had High/Low/—
for r in range(5, 20):
    n = ws.cell(r, 5).value
    exp = "—" if n == "—" else "Universal" if n == 9 else "High" if n >= 5 else "Moderate" if n >= 3 else "Low"
    assert ws.cell(r, 6).value == exp, (r, n, ws.cell(r, 6).value)
ws["A1"] = ("Table S9. Cross-scale integration of 15 anchor reactions — RQ1 diet response traced through genetic conservation (RQ2; "
            "aggregate-design standard FBA, union 243 / universal core 11), cell-type attribution (RQ3), and diet–microbiome decomposition (RQ4)")
ws["A2"] = ("Anchor set: 15 reactions from RQ1 HFD vs SCD with |ΔHFD−SCD| ≥ 0.20 mmol·gDW⁻¹·h⁻¹ (Supplementary Methods SM9.2). "
            "Conservation tier: Universal = 9/9 strains; High ≥ 5/9; Moderate 3–4/9; Low 1–2/9; — = 0/9 (not in the 243-reaction "
            "magnitude-responsive union). VE = variance explained (%). '—' elsewhere = reaction not in the respective downstream analysis.")
ws.row_dimensions[1].height = 36.0
ws.row_dimensions[2].height = 42.0
LOG.append("S9: conservation tiers EX_nh4_e High->Moderate, CO2tm Low->Moderate, NAt '—'->Moderate (restored to SM9.3 rule / original); "
           "'Supplementary Methods S9.2' -> 'SM9.2'; tier legend 'Low < 3/9' -> 'Low 1–2/9; — = 0/9'")

# ----------------------------------------------------------------------------- Notes_RQ3 / Notes_SJ (original sheets; text edits)
ws = wb["Notes_RQ3"]
ws["A1"] = "RQ3 Supplementary Tables (SA–SI) — Data Sources and Methodology (Methods 2.4; Results 3.4; SM7)"
ws["B3"] = ("23 hepatic cell types ranked by abundance-weighted contribution to WD vs. chow flux changes (P99/biomass reference). "
            "Per-cell contribution index = contribution(%) / abundance(%). Source: RQ3_phase2_contribution_analysis.csv")
ws["B4"] = ("Cell-type contribution (%) projected across 9 mouse strain contexts by reweighting the single atlas-derived WD-vs-chow attribution "
            "template; evaluates robustness of the template, not strain-specific single-cell evidence. Source: RQ3_contribution_matrix.csv")
ws["B6"] = "Cross-strain CV per cell type for the projected contributions (low CV = template robust to strain reweighting). Source: RQ3_conservation_analysis.csv"
ws["B7"] = ("Number of threshold-responsive reactions per cell type (deterministic ≥ 2-of-3 rule; see Response criterion below). "
            "Source: RQ3_statistical_tests.csv + RQ3_phase2_contribution_analysis.csv")
ws["B10"] = ("205 bulk reference reactions with cell-type attribution category (unique / cooperative / multicellular / non-cellular). "
             "Source: RQ3_phase3_reaction_attribution.csv")
ws["A13"] = "Response criterion"
ws["B13"] = ("Threshold-responsive (deterministic; not a statistical test): at least two of |v_WD / v_chow| > 1.5, |Δv| > 0.1, and relative "
             "contribution to total flux change > 0.5 (SM7.3). One pseudobulk profile per cell type and condition precludes replicate-level FDR testing.")
ws["B14"] = "Western diet (WD) vs. chow"
for r in (3, 4, 6, 7, 10, 13):
    ws[f"B{r}"].alignment = Alignment(wrap_text=True, vertical="center")
    ws.row_dimensions[r].height = 11.5 * lines(ws[f"B{r}"].value, 90 * 1.25) + 4
LOG.append("Notes_RQ3: 'FDR < 0.10' significance threshold replaced by the SM7.3 threshold-responsive rule; 'WD vs SCD' -> 'WD vs chow'; "
           "SB described as a heuristic projection; title section reference updated")

ws = wb["Notes_SJ"]
ws["A1"] = "SJ Supplementary Tables (SJ1–SJ4) — Source Files and Methods (Results 3.5)"
ws["B3"] = "Key statistics cited in Results 3.5. Part A: global r, p, n, directional agreement. Part B: highlighted subsystem concordance."
LOG.append("Notes_SJ: 'Section 2.5' -> 'Results 3.5'")

# ============================================================================= Notes_RQ1 / Notes_RQ2 / Notes_SK (rebuilt)


def notes_sheet(name, title, items, wA=45, wB=80):
    ws = fresh(name)
    widths(ws, {"A": wA, "B": wB})
    title_row(ws, title, 2, "1F3864")
    header_row(ws, 2, ["Item", "Details"], "2E5090")
    r = 3
    for it in items:
        if it is None:
            r += 1; continue
        a, b = it
        bold = a.endswith(":") or a in ("Source files", "Definitions")
        ca, cb = ws.cell(r, 1, a.rstrip(":")), ws.cell(r, 2, b)
        style(ca, font(8, b=bold), None, Alignment(horizontal="left", vertical="center", wrap_text=True), HAIR)
        style(cb, font(8), None, Alignment(horizontal="left", vertical="center", wrap_text=True), HAIR)
        ws.row_dimensions[r].height = max(12.0, 11.0 * max(lines(b, wB * 1.25), lines(a, wA * 1.25)) + 3)
        r += 1
    return ws


notes_sheet("Notes_RQ1", "Data Sources and Methodology Notes — RQ1 (Tables S1–S2, SP1–SP2, SP5)", [
    ("S1_Convergent_All3Diets", "10 reactions significant (BH-FDR q < 0.05) in each of HFD, KD, and WD vs. the study-specific control diet (SCD); "
                                "9 of 10 are directionally concordant (GLNALANaEx switches direction by diet). Methods 2.3; Results 3.2; SR1.1."),
    ("S2_All_Significant_Reactions", "450 reactions significant in ≥ 1 of the five estimable pairwise contrasts (HFD–SCD, KD–SCD, WD–SCD, HFD–KD, "
                                     "HFD–WD); 262 are significant in ≥ 1 SCD-referenced contrast. KD–WD is non-estimable (no cohort contains both "
                                     "diets). Colour: green = significant in all 3 SCD-referenced contrasts; yellow = 2; white = 1; grey = "
                                     "significant only in a between-diet contrast."),
    ("SP1_pFBA_Convergent27 / SP2_pFBA_AllSignificant228", "pFBA sensitivity analysis of the same cohorts (Results 3.1; SM10.6): 27-reaction shared "
                                     "signature (19 biochemical reactions, 8 bookkeeping/pseudo-reactions) and 228 reactions significant in ≥ 1 "
                                     "SCD-referenced contrast. Reported as a reproducible pFBA-associated signature, not a uniquely constrained flux core."),
    ("SP5_Sensitivity_Summary", "Consolidated FBA-vs-pFBA, E-Flux cap, objective, and FVA/pFVA sensitivity results for RQ1 and RQ2 "
                                 "(SM10.6), plus the RQ2 magnitude-threshold θ sweep (section G; Results 3.3; Methods 2.4)."),
    None,
    ("Source files", ""),
    ("RQ1_reaction_stats_HFD_vs_SCD.csv", "Differential flux statistics, HFD vs SCD (n = 38 vs 37 samples; six cohorts)"),
    ("RQ1_reaction_stats_KD_vs_SCD.csv", "Differential flux statistics, KD vs SCD (KD n = 12; two cohorts)"),
    ("RQ1_reaction_stats_WD_vs_SCD.csv", "Differential flux statistics, WD vs SCD (WD n = 12; one cohort)"),
    ("RQ1_reaction_stats_HFD_vs_KD.csv", "Differential flux statistics, HFD vs KD (two cohorts)"),
    ("RQ1_reaction_stats_HFD_vs_WD.csv", "Differential flux statistics, HFD vs WD (one cohort)"),
    None,
    ("Statistical test", "Standard FBA computed per biological replicate; Welch's t-test with Benjamini–Hochberg FDR within each contrast; "
                         "significance q < 0.05. Cohen's d is reported as the effect size."),
    ("Known limitation", "For contrasts spanning more than one GEO dataset, the primary test pools samples without an explicit dataset covariate. "
                         "A dataset-adjusted model (Methods 2.3, Eq. 5; HC3 standard errors) was applied to the pFBA sensitivity fluxes; see "
                         "Discussion 4.7."),
    ("Flux units", "mmol/gDW/h (millimoles per gram dry weight per hour)"),
    ("Model", "iMM1415 mouse genome-scale metabolic model (3,726 reactions)"),
])

notes_sheet("Notes_RQ2", "Data Sources and Methodology Notes — RQ2 (Tables S3–S5, S3_Alt, S4_Alt, SP3–SP4)", [
    ("S3_Conserved_Core_11", "11 reactions with |mean ΔHFD−SCD| ≥ 0.20 in all 9 strains. Standard FBA with expression averaged within each "
                             "strain–diet group before optimization (aggregate design; Methods 2.4, Eqs. 6–7; SM6). Modules: gas exchange (CO2t, "
                             "O2t, EX_co2_e, EX_o2_e), acetoacetate transport (ACACt2, EX_acac_e), tyrosine degradation (34HPPOR, FUMAC, HGNTOR, "
                             "MACACI, TYRTA)."),
    ("S4_All_Strains_Matrix", "243 reactions magnitude-responsive in ≥ 1 of 9 strains. Conservation spectrum (N strains 1–9): 80, 40, 38, 29, 18, "
                              "13, 9, 5, 11. Colour: green = 9; blue = 7–8; yellow = 4–6; white = 1–3. Tiers follow SM6.4."),
    ("L2 / L3 evidence columns (S3, S4)", L_LEGEND),
    ("S5_PerStrain_Summary", "Per-strain magnitude-responsive, strain-unique, and universal-core counts (aggregate design), plus the number of all "
                             "3,726 reactions reaching within-strain Welch BH-FDR < 0.05 on replicate-preserving flux solutions. NZO/HlLtJ: HFD "
                             "n = 1, within-strain inference non-estimable (NE)."),
    ("S3_Alt / S4_Alt", "Complementary replicate-preserving standard-FBA analysis (per-animal profiles modeled independently): 187-reaction "
                        "union and 15-reaction universal core under the same edge-magnitude definition. All 11 aggregate-design core reactions "
                        "are contained in the 15-reaction core; 156/243 (64.2%) of the aggregate-design union overlaps the 187-reaction union."),
    ("SP3_pFBA_UniversalCore16 / SP4_pFBA_AllStrainsMatrix105", "pFBA sensitivity analysis with the same magnitude threshold: 16-reaction universal "
                        "core and 105-reaction union (Results 3.1; SM10.6). Targeted pFVA retained same-direction intervals for 14/16 reactions "
                        "at 1% and 7/16 at 5% tolerance."),
    None,
    ("Source files", ""),
    ("RQ2_<strain>_edges_HFD_vs_SCD.csv", "Per-strain HFD vs SCD flux differences, aggregate design (one row per reaction–metabolite edge; "
                                          "de-duplicated by reaction for the counts above)"),
    ("flux_pairwise_stats.csv (per strain)", "Replicate-level Welch's t-test statistics with within-strain BH-FDR (replicate-preserving design)"),
    None,
    ("Flux units", "Δ = HFD mean flux − SCD mean flux (mmol/gDW/h); threshold θ = 0.20"),
    ("Model", "iMM1415 mouse genome-scale metabolic model; nine inbred founder strains (GSE182668)"),
])

notes_sheet("Notes_SK", "RQ4 Supplementary Tables (SK1–SK6) — Source Files and Methods (Methods 2.4–2.5; Results 3.6; SM8)", [
    ("SK1_Species_Composition", "All 57 species detected in metatranscriptomic data (GSE104913) with ND/SCD and DD/HFD abundances, direction "
                                "(31 decreased, 19 increased, 7 unchanged) and fold change. 27 species retained in the strict primary AGORA2 map "
                                "(bold). Source: species_activity_weighted_composition.csv; species_to_model_mapping.csv."),
    ("SK2_Portal_Metabolites", "Portal metabolite community exchange fluxes under both conditions (primary configuration). Acetate is the only "
                               "curated metabolite with positive export above threshold and the only one propagated to the hepatic model. "
                               "Source: portal_metabolite_production.csv."),
    ("SK3_Flux_Attribution", "Part A: driver categories across all 3,726 reactions (260 diet; 115 diet opposed by microbiome; 27 microbiome; "
                             "2 microbiome opposed by diet; 3,322 stable). Part B: 307 reactions with VE_micro > 0. 47 reactions exceed "
                             "|Δmicrobiome| ≥ 0.01; 29 are microbiome-dominant. Source: flux_attribution_analysis.csv."),
    ("SK4_Pathway_Synergy", "Subsystems with ≥ 2 reactions with |Δtotal| > 10⁻³ classified as Highly Synergistic / Moderately Synergistic / "
                            "Mixed / Highly Antagonistic. Method: rq4_pathway_enrichment_module.py."),
    ("SK5_Pathway_Enrichment", "Hypergeometric enrichment (background = 3,726 reactions; 'significant' = |Δtotal| > 10⁻³, n = 404) for 55 "
                               "annotated subsystems, BH-FDR corrected. Method: rq4_pathway_enrichment_module.py."),
    ("SK6_Compartment_Enrichment", "Same enrichment test by subcellular compartment (5 compartments)."),
    None,
    ("Primary configuration", "27-species map; MICOM cooperative tradeoff 0.50 with community pFBA; host pFBA with the biomass-associated "
                              "objective; portal scale 0.10; forced coupling."),
    ("Variance explained (VE)", "VE_micro = Δmicro² / (Δdiet² + Δmicro²) × 100; VE_diet = 100 − VE_micro (Methods 2.5, Eq. 11). Computed "
                                "per reaction; means reported for subsets."),
    ("Portal sign convention", "Positive community exchange flux = net microbial export (candidate portal delivery); negative = net microbial uptake."),
], wA=30, wB=95)

# ============================================================================= Contents
ws = wb["Contents"]
for mr in list(ws.merged_cells.ranges):
    ws.unmerge_cells(str(mr))
for row in ws.iter_rows():
    for c in row:
        c.value = None; c.style = "Normal"
for r in list(ws.row_dimensions):
    ws.row_dimensions[r].height = None
RQ1S, RQ2S = "RQ1 — Methods 2.3; Results 3.2", "RQ2 — Methods 2.4; Results 3.3; SM6"
SENS = "Sensitivity — Results 3.1; SM10.6"
RQ3S, SJS, RQ4S = "RQ3 — Methods 2.4; Results 3.4; SM7", "RQ3 validation — Results 3.5", "RQ4 — " + RQ4REF
CONTENTS = [
    ("Table S1", "Reactions significant in all three diets vs. SCD (HFD, KD, WD; n = 10)", RQ1S, "S1_Convergent_All3Diets"),
    ("Table S2", "All reactions significant in any estimable pairwise dietary comparison (n = 450)", RQ1S, "S2_All_Significant_Reactions"),
    ("Table S3", "Universal HFD-responsive core across all 9 strains (n = 11), with L2/L3 statistical evidence", RQ2S, "S3_Conserved_Core_11"),
    ("Table S4", "Reaction × strain HFD flux-change matrix (union n = 243), with L2/L3 statistical evidence", RQ2S, "S4_All_Strains_Matrix"),
    ("Table S5", "Per-strain summary of HFD responses and within-strain statistical support", RQ2S, "S5_PerStrain_Summary"),
    ("Table S3_Alt", "Replicate-preserving analysis: universal core across all 9 strains (n = 15)", RQ2S, "S3_Alt_ReplicateDesign_Core15"),
    ("Table S4_Alt", "Replicate-preserving analysis: reaction × strain matrix (union n = 187)", RQ2S, "S4_Alt_ReplicateDesign_Union187"),
    ("Table SP1", "pFBA sensitivity: reactions shared by HFD, KD, and WD vs. SCD (n = 27)", SENS, "SP1_pFBA_Convergent27"),
    ("Table SP2", "pFBA sensitivity: reactions significant in ≥ 1 SCD-referenced contrast (n = 228)", SENS, "SP2_pFBA_AllSignificant228"),
    ("Table SP3", "pFBA sensitivity: universal magnitude-responsive core across 9 strains (n = 16)", SENS, "SP3_pFBA_UniversalCore16"),
    ("Table SP4", "pFBA sensitivity: reaction × strain matrix (union n = 105)", SENS, "SP4_pFBA_AllStrainsMatrix105"),
    ("Table SP5", "Consolidated FBA/pFBA, cap, objective, and FVA/pFVA sensitivity summary, plus RQ2 θ-threshold sweep", SENS, "SP5_Sensitivity_Summary"),
    ("Table S-A", "Cell-type abundance-weighted contribution summary (WD vs. chow)", RQ3S, "SA_Contribution_Summary"),
    ("Table S-B", "Cell-type contribution projected across 9 strain contexts", RQ3S, "SB_CrossStrain_Matrix"),
    ("Table S-C", "Overlap between bulk and cellular flux signals (F1 scores)", RQ3S, "SC_BulkCellular_Overlap"),
    ("Table S-D", "Cross-strain variability (CV) of projected cell-type contributions", RQ3S, "SD_Conservation_Analysis"),
    ("Table S-E", "Threshold-responsive reactions per cell type", RQ3S, "SE_Responsive_Rxns_CellType"),
    ("Table S-F", "Pathway-level attribution by cell type", RQ3S, "SF_Pathway_Attribution"),
    ("Table S-G", "Driver consistency scores across strains and cell types", RQ3S, "SG_Driver_Consistency"),
    ("Table S-H", "Per-reaction attribution by lineage/function/location", RQ3S, "SH_Reaction_Attribution"),
    ("Table S-I", "scRNA-seq data summary used for cell abundance weighting", RQ3S, "SI_scRNAseq_Summary"),
    ("Table SJ-1", "Bulk vs. single-cell composite validation summary metrics", SJS, "SJ1_Validation_Summary"),
    ("Table SJ-2", "Subsystem-level flux correlations (validation)", SJS, "SJ2_Subsystem_Correlations"),
    ("Table SJ-3", "Per-reaction flux comparison between bulk and composite models", SJS, "SJ3_PerReaction_Comparison"),
    ("Table SJ-4", "Composite model construction details", SJS, "SJ4_Composite_Construction"),
    ("Table SK1", "Gut microbial species composition by diet (all 57 detected species)", RQ4S, "SK1_Species_Composition"),
    ("Table SK2", "Portal metabolite production by the microbial community", RQ4S, "SK2_Portal_Metabolites"),
    ("Table SK3", "Reaction-level diet–microbiome flux attribution", RQ4S, "SK3_Flux_Attribution"),
    ("Table SK4", "Pathway-level diet–microbiome interaction classification", RQ4S, "SK4_Pathway_Synergy"),
    ("Table SK5", "Pathway enrichment for diet–microbiome interactions", RQ4S, "SK5_Pathway_Enrichment"),
    ("Table SK6", "Compartment-level enrichment analysis", RQ4S, "SK6_Compartment_Enrichment"),
    ("Table S9", "Cross-scale integration of 15 RQ1 anchor reactions through RQ2, RQ3, and RQ4", "Cross-scale — Methods 2.6; SM9", "S9_CrossScale_Integration"),
    ("Notes (RQ1)", "Data sources and methodology notes for Tables S1–S2, SP1–SP2, SP5", "RQ1", "Notes_RQ1"),
    ("Notes (RQ2)", "Data sources and methodology notes for Tables S3–S5, S3_Alt, S4_Alt, SP3–SP4", "RQ2", "Notes_RQ2"),
    ("Notes (RQ3)", "Data sources and methodology notes for Tables S-A to S-I", "RQ3", "Notes_RQ3"),
    ("Notes (SJ)", "Source files and methods for Tables SJ-1 to SJ-4", "RQ3 validation", "Notes_SJ"),
    ("Notes (RQ4)", "Source files and methods for Tables SK1–SK6", "RQ4", "Notes_SK"),
]
widths(ws, {"A": 16, "B": 86, "C": 38, "D": 34})
ws.merge_cells("A1:D1"); ws.merge_cells("A2:D2")
c = ws["A1"]; c.value = "Supplementary Tables — Combined"
style(c, font(14, b=True, color="FFFFFF"), fill("1F4E79"), Alignment(horizontal="center", vertical="center"))
ws.row_dimensions[1].height = 30.0
c = ws["A2"]; c.value = "Multi-scale constraint-based modeling reveals conserved and context-dependent hepatic metabolic adaptation to diet"
style(c, font(10, i=True, color="404040"), None, Alignment(horizontal="center", vertical="center"))
ws.row_dimensions[2].height = 20.1
for k, h in enumerate(["Table", "Description", "Section", "Sheet Name"], 1):
    style(ws.cell(4, k, h), font(10, b=True, color="FFFFFF"), fill("4472C4"), Alignment(horizontal="center", vertical="center"), THIN)
ws.row_dimensions[4].height = 18.0
for i, (a, b, s, n) in enumerate(CONTENTS):
    r = 5 + i
    fc = "EBF3FB" if i % 2 == 0 else "FFFFFF"
    for k, v in enumerate((a, b, s, n), 1):
        cc = ws.cell(r, k, v)
        style(cc, font(10, color="1F4E79"), fill(fc), Alignment(horizontal="left" if k == 2 else "center", vertical="center"), THIN)
    ws.row_dimensions[r].height = 16.15
    if n in wb.sheetnames:
        ws.cell(r, 4).hyperlink = f"#'{n}'!A1"

# ============================================================================= order + cleanup
ORDER = [row[3] for row in CONTENTS]
ORDER = ["Contents"] + [n for n in ORDER if n not in ("Notes_RQ1", "Notes_RQ2", "Notes_RQ3", "Notes_SJ", "Notes_SK")] + \
        ["Notes_RQ1", "Notes_RQ2", "Notes_RQ3", "Notes_SJ", "Notes_SK"]
missing = set(ORDER) - set(wb.sheetnames); extra = set(wb.sheetnames) - set(ORDER)
assert not missing and not extra, (missing, extra)
wb._sheets = [wb[n] for n in ORDER]
for s in wb.worksheets:
    s.sheet_view.tabSelected = False
    s.freeze_panes = "A5" if s.title == "S9_CrossScale_Integration" else None
wb.active = 0
wb.worksheets[0].sheet_view.tabSelected = True
wb.save(OUT)
print("\n".join(LOG))
print("sheets:", len(wb.sheetnames))
