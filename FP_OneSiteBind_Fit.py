# -*- coding: utf-8 -*-
# Version 2026-11-24a: pre-release fixes — pooled R², safe report cleanup,
# 1-5 replicate support, probe fit on concatenated data, fit validation.
"""
Per-sample concatenated fit of FP_OneSiteBind for OriginPro's embedded Python.
Pure Python (no numpy/pandas needed).

Data layout: the active worksheet contains repeating blocks of
    [X column][Y1][Y2]...[Yn]        (n = 1...5 replicate Y columns)
(one block per sample). X columns should be designated as X in the
worksheet; if they are not, the script falls back to a fixed GROUP_SIZE
layout (X + GROUP_SIZE-1 Y columns repeating).

Model (user-defined fitting function FP_OneSiteBind):
    A = Amin - ((Kd + P + L - sqrt((Kd + P + L)^2 - 4*P*L)) / (2*P)) * (Amin - Amax)
    L : total concentration of the titrated protein (the X data), µM
    P : total concentration of the labeled probe, FIXED (default 0.05 µM)
    A : fluorescence anisotropy signal (the Y data), mAu
    Fitted per sample: Kd, Amax, Amin (shared by all replicates of a sample).

For EACH sample: the replicate Y columns are concatenated (valid X/Y pairs
only) onto a temporary worksheet and fitted as ONE curve -> one shared
Kd, Amax, Amin per sample. The parameter names/order and the index of P
are discovered automatically from a probe fit, so nothing about the FDF
needs to be hardcoded here.

How to run:
    1. Open your Origin project and click the data worksheet (not a
       report sheet)
    2. View > Code Builder -> open this file -> press F5

What it produces:
    - "FitCurves"  sheet : smooth fit curve for every sample (X, Fit pairs);
      named FitCurves, FitCurves2, ... on repeated runs
    - "FP_Summary" sheet : Sample | <param> | SE | ... | R² | Status
      (FP_Summary2, ... on repeated runs)
    - One graph per sample (clone of TEMPLATE_GRAPH) with the replicate
      datasets and the fit curve; Long Name <sample>_graph

Graphs and result sheets are NEVER deleted — each run creates new ones.
Only the fit report sheets created by the CURRENT run are removed
(see KEEP_REPORTS); pre-existing reports from other analyses are untouched.
"""

import math
import os
import re

import originpro as op

# ========================== USER SETTINGS ==========================
FUNC_NAME = 'FP_OneSiteBind'

# FDF parameter names to fit (shared across the 3 Y columns of a sample,
# because they are concatenated into one curve). Exact FDF spelling.
FIT_PARAMS = ['Kd', 'Amax', 'Amin']

# display units and role (for initial-value guessing) of each parameter
PARAM_INFO = {
    'Kd':   {'role': 'kd',   'units': 'µM'},
    'Amax': {'role': 'rmax', 'units': 'mAu'},
    'Amin': {'role': 'r0',   'units': 'mAu'},
}

FIXED_PARAMS = {'P': 0.05}          # <<< EDIT ME: P in µM (same for all samples)
P_VALUE = FIXED_PARAMS['P']         # convenience: P used for model evaluation/labels

GROUP_SIZE = 4            # X + (GROUP_SIZE-1) Y; fallback layout used ONLY
                          # when the X columns are not designated as X
SET_INITIAL_VALUES = True # False = let the FDF auto-initialization decide

X_TITLE = '[PDZ1], µM'
Y_TITLE = 'FA, mAu'
TEMPLATE_GRAPH = 'Graph1'   # graph in the project used as the output template;
                            # set to '' to build a clean default graph instead
DATA_COLORS = ['Black', 'Red', 'Blue', 'Green', 'Orange']
                          # colors of the replicate Y datasets (cycled if a
                          # sample has more replicates than colors)
FIT_COLOR = 'Black'                      # color of the fit curve
FIT_WIDTH_PT = 3                         # fit curve line width, points
GRAPH_SUFFIX = '_graph'     # long name of each output graph: <sample>_graph
AXIS_TITLE_FONT = 31        # axis title font size (pts)
TICK_LABEL_FONT = 22        # tick label font size (pts)
LEGEND_FONT = 22            # legend font size (pts)
SYMBOL_SHAPE = 2            # data symbol shape (LabTalk set -k list; 2 = circle)
LEGEND_LINE_WIDTH = 2       # legend line width (pts) for the data entries
X_FROM = 0.0001             # fixed X axis scale: from ...
X_TO = 100                  # fixed X axis scale: to ...
N_FIT_POINTS = 300          # points in each smooth fit curve
KEEP_REPORTS = False        # True = keep the fit report sheets created by this
                            # run (default: deleted after reading)
CLEANUP_PREVIOUS = True     # True = delete leftover FPTmp temp sheets from
                            # earlier interrupted runs. Graphs, result sheets
                            # and pre-existing fit reports are NEVER deleted
SAVE_PROJECT = False        # True = save the .opju when done
# ===================================================================


def _norm(name):
    """Normalize a parameter name for comparison (ignore case and characters
    like '*', which Origin sometimes appends in report tables)."""
    return re.sub(r'[^A-Za-z0-9]', '', str(name)).lower()


def _ok(v):
    """True if v is a usable finite number (filters NaN, text, Origin's missing)."""
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return False
    return v == v and abs(v) < 1e250


# Origin color indices (the standard LabTalk color list) — used to build
# line-only legend entries
OCOLOR_INDEX = {'black': 1, 'red': 2, 'green': 3, 'blue': 4, 'cyan': 5,
                'magenta': 6, 'yellow': 7, 'dark yellow': 8, 'navy': 9,
                'purple': 10, 'wine': 11, 'olive': 12, 'dark cyan': 13,
                'royal': 14, 'orange': 15, 'violet': 16, 'pink': 17,
                'white': 18, 'lt gray': 19, 'gray': 20, 'lt yellow': 21,
                'lt cyan': 22, 'lt magenta': 23, 'dark gray': 24}


def _color_index(name):
    """Origin color index for a color name (ints pass through)."""
    if isinstance(name, int):
        return name
    return OCOLOR_INDEX.get(str(name).strip().lower(), 1)


def fp_model(L, Kd, Amax, Amin, P):
    """Exact transcription of the FDF formula (numerically stable form):

        A = Amin - ((Kd+P+L - sqrt((Kd+P+L)^2 - 4*P*L)) / (2*P)) * (Amin - Amax)

    (a - s)/(2P) == 2L/(a + s) with a = Kd+P+L, s = sqrt(a^2 - 4*P*L)
    """
    out = []
    for Lv in L:
        a = Kd + P + Lv
        frac = 2.0 * Lv / (a + math.sqrt(a * a - 4.0 * P * Lv))
        out.append(Amin + frac * (Amax - Amin))
    return out


def find_groups(wks):
    """Return a list of sample groups: {'name', 'xcol', 'ycols'} (1-offset cols)."""
    ncols = wks.cols

    # collect all X-designated columns
    xcols = []
    for i in range(1, ncols + 1):
        op.lt_exec(f'range fpCC = wcol({i}); double fpT = fpCC.type;')
        if int(op.po.LT_get_var('fpT')) == 4:
            xcols.append(i)
    op.lt_exec('del -v fpCC; del -v fpT;')

    if len(xcols) == 1 and xcols[0] == 1:
        # Exactly one X-designated column, in first position: either ONE
        # sample with any number of replicate Y columns, or a multi-sample
        # sheet whose X columns were not designated. If the sheet splits
        # evenly into GROUP_SIZE blocks AND the Long Names of the blocks'
        # first Y columns differ, assume the fixed multi-sample layout.
        def _ln(c):
            return (wks.get_label(c - 1, 'L') or '').strip()
        if (ncols % GROUP_SIZE == 0 and ncols > GROUP_SIZE
                and _ln(2) and _ln(1 + GROUP_SIZE)
                and _ln(2) != _ln(1 + GROUP_SIZE)):
            xcols = list(range(1, ncols + 1, GROUP_SIZE))
    elif len(xcols) < 2:
        # cannot infer the block starts from the X designation:
        # fall back to the fixed periodic layout X + (GROUP_SIZE-1) Y
        if ncols % GROUP_SIZE == 0:
            xcols = list(range(1, ncols + 1, GROUP_SIZE))
        else:
            raise RuntimeError(
                f'Could not determine the sample layout: {ncols} columns, '
                f'X-designated columns: {xcols}. Designate the X columns in the '
                f'worksheet (Column: Set as X) or adjust GROUP_SIZE.')

    bounds = xcols + [ncols + 1]
    groups = []
    for k, xc in enumerate(xcols):
        ycols = []
        for i in range(xc + 1, bounds[k + 1]):
            vals = wks.to_list(i - 1)
            if any(_ok(v) for v in vals):
                ycols.append(i)
        if not ycols:
            continue
        name = (wks.get_label(ycols[0] - 1, 'L') or '').strip()
        if not name:
            name = wks.get_label(ycols[0] - 1, 'G') or f'Sample{k + 1}'
        groups.append({'name': name, 'xcol': xc, 'ycols': ycols})
    if not groups:
        raise RuntimeError('No sample groups with data found.')
    return groups


def get_xy(wks, xcol, ycol):
    """Read one X/Y pair, drop rows with missing values."""
    x = wks.to_list(xcol - 1)
    y = wks.to_list(ycol - 1)
    pairs = [(xv, yv) for xv, yv in zip(x, y) if _ok(xv) and _ok(yv)]
    return [p[0] for p in pairs], [p[1] for p in pairs]


def guess_initial_values(datasets):
    """Rough starting values by role: kd, rmax, r0."""
    ally = [v for _, yy in datasets for v in yy]
    r0 = min(ally)
    rmax = max(ally)
    pos = [v for xx, _ in datasets for v in xx if v > 0]
    kd = math.exp(sum(math.log(v) for v in pos) / len(pos)) if pos else 1.0
    if rmax > r0:
        pairs = sorted((xy for ds in datasets for xy in zip(*ds)))
        ymid = 0.5 * (r0 + rmax)
        for i in range(1, len(pairs)):
            x1, y1 = pairs[i - 1]
            x2, y2 = pairs[i]
            if y1 < ymid <= y2:
                kd = x1 + (ymid - y1) * (x2 - x1) / (y2 - y1)
                break
    return {'kd': kd, 'rmax': rmax, 'r0': r0}


def run_global_fit(wks, xcol, ycols, init_idx, fix_idx):
    """Fit a single (X, Y) dataset (already concatenated by the caller).

    init_idx / fix_idx: lists of (parameter index 1-based, value).
    Everything is done by INDEX — never by name — so parameter names with
    special characters (like 'Kd*') cannot break the LabTalk script.

    Returns (rep, created): rep = the report worksheet (or None) and
    created = names of the FitNL... sheets created by THIS fit (so that
    only they, never pre-existing reports, are cleaned up afterwards).
    """
    before = {sh.name.lower() for sh in wks.get_book()}
    cols = ', '.join(str(c) for c in [xcol] + ycols)
    rng = f'{wks.lt_range()}!({cols})'

    lt = [
        f'range fpAA = {rng};',
        f'nlbegin fpAA {FUNC_NAME} fpTT;',
    ]
    for i, v in init_idx:
        lt.append(f'fpTT.p{i} = {v};')         # initial value
    for i, v in fix_idx:
        lt.append(f'fpTT.f{i} = 1;')           # fix parameter
        lt.append(f'fpTT.p{i} = {v};')         # at this value
    lt += ['nlfit;', 'nlend 1;']
    op.lt_exec('\n'.join(lt))

    # sheets created by this fit (reports + residuals) — determined in a
    # SEPARATE LabTalk-free step; reading a missing tree node inside the
    # fit block would abort the remaining statements of that block
    created = [sh.name for sh in wks.get_book()
               if sh.name.lower().startswith('fitnl')
               and sh.name.lower() not in before]

    # locate the report sheet via __REPORT$ (do NOT confuse it with the
    # residuals sheet, which is named FitNLCurve...)
    rep = None
    rep_name = (op.po.LT_get_str('__REPORT$') or '').strip().strip('!')
    m = re.match(r'^\[([^\]]+)\](.+)$', rep_name or '')
    if m:
        wbk = op.find_book('w', m.group(1))
        if wbk:
            for sh in wbk:
                if sh.name == m.group(2):
                    rep = sh
                    break
    return rep, created


def parse_parameters(rep_wks):
    """Read the Parameters table from the fit report via GetReportTable.

    Origin exports the table to a temp tab-delimited file; we parse it with
    plain Python (no pandas needed).

    Returns an ordered list of (name, value, error) — the order is the
    FDF parameter order.
    """
    ff = op.po.LT_get_str(
        f'GetReportTable({rep_wks.lt_range()}, Parameters, 1, 1, 1)$')
    if not ff:
        raise RuntimeError('Could not export the Parameters table from the fit report.')
    try:
        with open(ff, encoding='utf-8-sig') as fh:
            lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    finally:
        try:
            os.remove(ff)
        except OSError:
            pass

    if len(lines) < 2:
        raise RuntimeError('Parameters table in the fit report is empty.')

    header = [h.strip() for h in lines[0].split('\t')]
    valcol = next((i for i, h in enumerate(header)
                   if 'value' in h.lower()), None)
    if valcol is None:
        raise RuntimeError(f'Unexpected report table layout: {header}')
    # the parameter-name column often has no header text; it sits left of Value
    namecol = next((i for i, h in enumerate(header)
                    if 'parameter' in h.lower()), max(valcol - 1, 0))
    errcol = next((i for i, h in enumerate(header)
                   if 'error' in h.lower()), None)

    rows = []
    for line in lines[1:]:
        cells = [c.strip() for c in line.split('\t')]
        if len(cells) <= max(valcol, namecol):
            continue
        name = cells[namecol]
        if not name:
            continue
        try:
            val = float(cells[valcol])
        except ValueError:
            continue
        err = float('nan')
        if errcol is not None and errcol < len(cells) and cells[errcol]:
            try:
                err = float(cells[errcol])
            except ValueError:
                pass
        rows.append((name, val, err))
    return rows


def delete_fit_sheets(wks, names):
    """Delete the report/residual sheets created by the CURRENT fit.

    Only the sheets named in `names` (captured right after the fit) are
    removed — pre-existing FitNL... sheets from other analyses are untouched.
    """
    if KEEP_REPORTS or not names:
        return
    targets = {n.lower() for n in names}
    for sh in list(wks.get_book()):
        if sh.name.lower() in targets:
            sh.destroy()


def _initial_index_list(order, guesses):
    """Build the [(param index, initial value)] list for FIT_PARAMS."""
    lower = [_norm(n) for n in order]
    init_idx = []
    for name in FIT_PARAMS:
        role = PARAM_INFO.get(name, {}).get('role')
        if _norm(name) not in lower:
            continue
        idx = lower.index(_norm(name)) + 1
        if SET_INITIAL_VALUES and role in guesses:
            init_idx.append((idx, guesses[role]))
    return init_idx


def discover_param_layout(wks, grp):
    """Probe fit (P left untouched) to learn the FDF parameter order.

    Uses the same safe concatenated data preparation as the production
    fits. Returns (order, fix_idx): order = list of parameter names as in
    the FDF, fix_idx = [(index, value)] for every entry of FIXED_PARAMS.
    """
    tmp, _ = _concat_dataset(wks, grp)
    try:
        # no initial values here: the parameter order is not known yet, so
        # the probe relies on the FDF's own auto-initialization
        rep, created = run_global_fit(tmp, 1, [2], [], [])
    finally:
        tmp.destroy()
    if rep is None:
        delete_fit_sheets(wks, created)
        raise RuntimeError('Probe fit failed: no report sheet was created. '
                           'Is the function name correct in the NLFit dialog?')
    try:
        rows = parse_parameters(rep)
    finally:
        delete_fit_sheets(wks, created)
    order = [name for name, _, _ in rows]
    if not order:
        raise RuntimeError('Probe fit produced no parameters — the fit did not run.')

    fix_idx = []
    missing = []
    for pname, pval in FIXED_PARAMS.items():
        idx = next((i + 1 for i, (n, _, _) in enumerate(rows)
                    if _norm(n) == _norm(pname)), None)
        if idx is None:
            missing.append(pname)
        else:
            fix_idx.append((idx, pval))
    if missing:
        raise RuntimeError(
            f'Fixed parameter(s) {missing} not found among the fit parameters '
            f'{order}. The FDF was probably changed — update FIXED_PARAMS or '
            f'the FDF so that the script and the model agree on the fixed '
            f'probe concentration.')
    return order, fix_idx


def _concat_dataset(wks, grp):
    """Concatenate all valid X/Y pairs of a sample onto a temp sheet.

    Returns (tmp_sheet, datasets), where datasets is the list of (x, y)
    pairs per replicate column. The caller must destroy the temp sheet.
    """
    datasets = [get_xy(wks, grp['xcol'], c) for c in grp['ycols']]
    xc, yc = [], []
    for x, y in datasets:
        xc += x
        yc += y
    if not xc:
        raise RuntimeError(
            f'[{grp["name"]}] no valid X/Y pairs found — check that the '
            f'sample block contains numeric data.')
    tmp = wks.get_book().add_sheet('FPTmp', active=False)
    tmp.from_list(0, xc, axis='X')
    tmp.from_list(1, yc, axis='Y')
    return tmp, datasets


def fit_sample(wks, grp, order, fix_idx):
    """Concatenate the sample's replicate Y columns into ONE dataset and fit it.

    Returns (values, r2, status, datasets).
    values: dict param name (as in FIT_PARAMS) -> (value, error).
    status: 'OK' or a note about the fit quality (shown in FP_Summary).
    """
    tmp, datasets = _concat_dataset(wks, grp)
    guesses = guess_initial_values(datasets)
    init_idx = _initial_index_list(order, guesses)

    try:
        rep, created = run_global_fit(tmp, 1, [2], init_idx, fix_idx)
    finally:
        tmp.destroy()
    if rep is None:
        delete_fit_sheets(wks, created)
        raise RuntimeError(f'[{grp["name"]}] fit report sheet not found.')
    try:
        rows = parse_parameters(rep)
    finally:
        delete_fit_sheets(wks, created)

    values = {}
    missing_se = []
    for name in FIT_PARAMS:
        match = next(((n, v, e) for n, v, e in rows
                      if _norm(n) == _norm(name)), None)
        if match is None:
            raise RuntimeError(
                f'[{grp["name"]}] parameter "{name}" not found in the fit report. '
                f'Parameters found: {[n for n, _, _ in rows]}.')
        _, val, err = match
        if not _ok(val):
            raise RuntimeError(
                f'[{grp["name"]}] parameter "{name}" has no finite value '
                f'({val}) — the fit did not converge. Inspect the data and '
                f'check that the titration covers the Kd range.')
        if _ok(err):
            values[name] = (val, err)
        else:
            missing_se.append(name)
            values[name] = (val, float('nan'))

    r2 = r_squared(datasets, tuple(values[n][0] for n in FIT_PARAMS))
    status = 'OK'
    if missing_se:
        status = 'SE unavailable: ' + ', '.join(missing_se)
    if r2 != r2:                      # NaN — no variance in the pooled Y
        status += ('; ' if status != 'OK' else '') + 'R² undefined'
    return values, r2, status, datasets


def r_squared(datasets, params):
    """Pooled R² over ALL replicate observations of one sample.

    The fit is a shared-parameter fit to the concatenated observations, so
    R² is computed on the same pooled data: SSE against the shared model,
    SST against the overall mean of all observations.
    """
    ally = [v for _, y in datasets for v in y]
    if not ally:
        return float('nan')
    ymean = sum(ally) / len(ally)
    ss_res = ss_tot = 0.0
    for x, y in datasets:
        f = fp_model(x, *params, P_VALUE)
        ss_res += sum((yv - fv) ** 2 for yv, fv in zip(y, f))
        ss_tot += sum((yv - ymean) ** 2 for yv in y)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def _is_output_sheet(name):
    n = name.lower()
    return (n.startswith('fitnl') or n.startswith('fitcurves')
            or n.startswith('fp_summary') or n.startswith('fptmp'))


def _unique_sheet_name(wb, base):
    """Return a sheet name based on `base` that is not yet used in the book
    (FitCurves -> FitCurves2 -> ...), so old output sheets are never replaced
    and graphs created in earlier runs keep their data references."""
    names = {sh.name.lower() for sh in wb}
    if base.lower() not in names:
        return base
    k = 2
    while f'{base}{k}'.lower() in names:
        k += 1
    return f'{base}{k}'


def cleanup_previous(wb):
    """Remove leftover FPTmp temp sheets from earlier interrupted runs.

    Graphs, result sheets (FitCurves, FP_Summary) and pre-existing fit
    reports are NEVER deleted — FitNL... reports may belong to unrelated
    analyses performed in the same workbook.
    """
    for sh in list(wb):
        if sh.name.lower().startswith('fptmp'):
            sh.destroy()


def _find_template_graph():
    """Return the graph page to use as the template, or None.

    TEMPLATE_GRAPH (by short name) takes priority; otherwise any existing
    output graph (Long Name ending with GRAPH_SUFFIX) is used - those are
    clones of the template anyway, and all plots/styles are re-applied.
    With TEMPLATE_GRAPH = '' no template is used at all (plain new graphs).
    """
    if not TEMPLATE_GRAPH:
        return None
    gp1 = op.find_graph(TEMPLATE_GRAPH)
    if gp1 is not None:
        return gp1
    for pg in op.po.GetPages():
        if pg.GetLongName().endswith(GRAPH_SUFFIX):
            gp = op.find_graph(pg.GetName())
            if gp is not None:
                return gp
    return None


def clone_template_graph(long_name):
    """Clone the template graph into a NEW window; never touches existing windows.

    The new graph keeps an automatic short name (Graph2, Graph3, ...) and
    gets the given Long Name. If TEMPLATE_GRAPH is missing, an existing
    output graph is used as the template instead; with TEMPLATE_GRAPH = ''
    a plain new graph is created.
    """
    gp1 = _find_template_graph()
    if gp1 is None:
        if TEMPLATE_GRAPH:
            raise RuntimeError(
                f"Template graph \"{TEMPLATE_GRAPH}\" not found and no previous "
                f"output graph (*{GRAPH_SUFFIX}) exists to fall back on. "
                f"Re-create it, or set TEMPLATE_GRAPH = '' to build a plain graph.")
        gp = op.new_graph()
        gp.lname = long_name
        return gp
    before = {p.GetName() for p in op.po.GetPages()}
    gp1.activate()
    op.lt_exec('win -d;')                      # duplicate the (active) template
    new = [p.GetName() for p in op.po.GetPages() if p.GetName() not in before]
    if not new:
        raise RuntimeError('Cloning the template graph failed (win -d).')
    gp = op.find_graph(new[0])
    gp.lname = long_name
    return gp


def _legend_lines_only(gl, n_data):
    """Rewrite the legend so every entry shows a plain LINE + name (no symbols).

    Uses the \\l(L color, thickness, style) legend escape: the entry for data
    plot n gets a line in that plot's color; the fit entry gets the fit color.
    n_data is the number of data plots — every entry after them is the fit.
    """
    lbl = gl.label('legend')
    if lbl is None:
        return
    txt = lbl.text
    if not txt:
        return

    def repl(m):
        try:
            p = int(m.group(1).split('.')[-1])   # plot index (ignore layer part)
        except ValueError:
            return m.group(0)
        if p <= n_data:
            color = _color_index(DATA_COLORS[(p - 1) % len(DATA_COLORS)])
            width = LEGEND_LINE_WIDTH
        else:
            color = _color_index(FIT_COLOR)
            width = FIT_WIDTH_PT
        return f'\\l(L {color}, {width}, 0)'

    new = re.sub(r'\\l\((\d+(?:\.\d+)?)[^)]*\)', repl, txt)
    if new != txt:
        lbl.text = new


def make_graph(wks, fwks, grp, fit_xcol, fit_ycol):
    """Put this sample's data + fit into a clone of TEMPLATE_GRAPH.

    Only the FIRST layer of the clone is styled — use a single-layer
    template graph (additional layers would keep their old plots).
    """
    gp = clone_template_graph(grp['name'] + GRAPH_SUFFIX)
    gl = gp[0]

    for p in list(gl.plot_list()):             # drop old plots from the clone
        p.remove()
    n_data = len(grp['ycols'])
    for c in grp['ycols']:
        gl.add_plot(wks, coly=c - 1, colx=grp['xcol'] - 1, type='s')
    pl = gl.add_plot(fwks, coly=fit_ycol - 1, colx=fit_xcol - 1, type='l')
    pl.color = FIT_COLOR                       # fit curve color
    pl.set_cmd(f'-wp {FIT_WIDTH_PT}')          # fit curve line width (points)

    # the replicate datasets in the configured colors (cycled if there are
    # more replicates than colors), all with the same symbol
    for i, pl in enumerate(gl.plot_list()[:n_data]):
        pl.color = DATA_COLORS[i % len(DATA_COLORS)]
        pl.symbol_kind = SYMBOL_SHAPE

    gl.lt_exec('layer.x.type = 2;')            # log10 X axis
    gl.lt_exec(f'layer.x.label.pt = {TICK_LABEL_FONT};')   # tick label size
    gl.lt_exec(f'layer.y.label.pt = {TICK_LABEL_FONT};')
    # axis lines on the top and on the right
    gl.lt_exec('layer.x.opposite = 1; layer.y.opposite = 1;')
    gl.lt_exec('layer.x.showopposite = 1; layer.y.showopposite = 1;')
    gl.rescale()
    # fix the X scale (after rescale, so it is not overridden); Y stays auto
    gl.lt_exec(f'layer.x.from = {X_FROM}; layer.x.to = {X_TO}; '
               f'layer.x.rescale = 1;')       # 1 = manual (fixed) scale
    gl.lt_exec('legendupdate mode:=reconstruct;')
    _legend_lines_only(gl, n_data)
    lbl = gl.label('legend')
    if lbl is not None:
        lbl.set_int('fsize', LEGEND_FONT)

    # axis titles: enlarge font; set text only if the template left it empty
    for obj, txt in (('xb', X_TITLE), ('yl', Y_TITLE)):
        lb = gl.label(obj)
        if lb is not None:
            lb.set_int('fsize', AXIS_TITLE_FONT)
            if not TEMPLATE_GRAPH or not lb.text.strip():
                lb.text = txt
    gp.show = True


def main():
    wks = op.find_sheet('w')
    if wks is None:
        raise RuntimeError('No active worksheet. Open your Origin project and '
                           'click the data sheet.')
    wks.activate()

    # if a fit report / output sheet is active, switch to the data sheet
    if _is_output_sheet(wks.name):
        for sh in wks.get_book():
            if not _is_output_sheet(sh.name):
                wks = sh
                wks.activate()
                break

    wb = wks.get_book()
    if CLEANUP_PREVIOUS:
        cleanup_previous(wb)

    groups = find_groups(wks)

    # probe fit on the first sample: learn parameter order and P's index
    order, fix_idx = discover_param_layout(wks, groups[0])

    # one sheet holding every sample's smooth fit curve as (X, Fit) pairs;
    # unique name per run so earlier sheets (and graphs using them) survive
    fwks = wb.add_sheet(_unique_sheet_name(wb, 'FitCurves'))

    rows = []                                   # summary rows: (name, values, r2, status)
    for k, grp in enumerate(groups):
        values, r2, status, datasets = fit_sample(wks, grp, order, fix_idx)
        rows.append((grp['name'], values, r2, status))

        # smooth fit curve for this sample
        xall = [v for xx, _ in datasets for v in xx]
        xmin, xmax = min(xall), max(xall)
        if xmin > 0:
            xfit = [xmin * (xmax / xmin) ** (i / (N_FIT_POINTS - 1))
                    for i in range(N_FIT_POINTS)]     # even on log axis
        else:
            step = (xmax - xmin) / (N_FIT_POINTS - 1)
            xfit = [xmin + i * step for i in range(N_FIT_POINTS)]
        cx, cy = 2 * k, 2 * k + 1
        fwks.from_list(cx, xfit, lname=f'{grp["name"]} X', axis='X')
        fwks.from_list(cy, fp_model(xfit, *[values[n][0] for n in FIT_PARAMS], P_VALUE),
                       lname=f'{grp["name"]} Fit', axis='Y')

        make_graph(wks, fwks, grp, cx + 1, cy + 1)

    # --- summary sheet ---------------------------------------------------
    swks = wb.add_sheet(_unique_sheet_name(wb, 'FP_Summary'))
    labels = ['Sample']
    for name in FIT_PARAMS:
        labels += [name, name + ' SE']
    labels += ['R²', 'Status']
    cols = [[r[0] for r in rows]]
    for name in FIT_PARAMS:
        cols.append([r[1][name][0] for r in rows])
        cols.append([r[1][name][1] for r in rows])
    cols.append([r[2] for r in rows])
    cols.append([r[3] for r in rows])
    for j, (lab, data) in enumerate(zip(labels, cols)):
        swks.from_list(j, data)
        swks.set_label(j, lab)

    if SAVE_PROJECT:
        op.save()


if __name__ == '__main__':
    main()
