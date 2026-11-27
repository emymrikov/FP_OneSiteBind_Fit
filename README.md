# FP_OneSiteBind_Fit

An OriginPro Python script for batch fitting of **fluorescence anisotropy**
(FA)/**fluorescence polarization** (FP) binding data with the user-defined fitting function `FP_OneSiteBind`.

For each sample, all replicate Y columns are **concatenated and fitted
together** (one shared parameter set per sample), and the script produces a
results summary, fit curves, and one styled graph per sample.
Fitting, result extraction, and graph generation are automated — no manual
interaction with the NLFit dialog is required. Initial setup (installing
the fitting function, preparing the worksheet and template graph) is still
necessary, see [Quick start](#quick-start).

## The model

The user-defined Origin fitting function `FP_OneSiteBind`
([`FP_OneSiteBind.fdf`](FP_OneSiteBind.fdf)):

```
A = Amin - ((Kd + P + L - sqrt((Kd + P + L)^2 - 4*P*L)) / (2*P)) * (Amin - Amax)
```

where

| Symbol | Meaning | Units |
|--------|---------|-------|
| `L`    | **total** concentration of the titrated protein (the X data) | µM |
| `P`    | **total** concentration of the labeled probe, **fixed** (default 0.05) | µM |
| `A`    | anisotropy signal (the Y data) | mAu |
| `Kd`   | dissociation constant (fitted) | µM |
| `Amax` | anisotropy of the fully bound probe (fitted) | mAu |
| `Amin` | anisotropy of the free probe (fitted) | mAu |

The equation is the exact ligand-depletion (Morrison) form for **1:1
equilibrium binding**: the bound fraction of the labeled probe is
`[PL]/P = (Kd+P+L - sqrt((Kd+P+L)^2-4PL)) / (2P)`. The fixed probe
concentration `P` must be positive and expressed in the same concentration
units as `L` and `Kd`.

The model assumes a linear mixture of free and bound probe signals; if
binding significantly changes the probe's fluorescence intensity, an
intensity correction may be needed before fitting.

## Installing the fitting function

The repository includes the Origin fitting function definition file (FDF)
[`FP_OneSiteBind.fdf`](FP_OneSiteBind.fdf). To install it in Origin:

1. Find your **User Files Folder**: `Help > Open Folder > User Files`
2. Copy `FP_OneSiteBind.fdf` into the `fitfunc` subfolder
3. Restart Origin
4. Verify: press **F9** (`Tools > Fitting Function Organizer`) —
   `FP_OneSiteBind` should appear under *User Defined*

The function defines the parameters `Amin, Amax, Kd, P` with `P` held
fixed. For the production fits the script overrides the stored initial
values (they are guessed from the data) and fixes `P` to
`FIXED_PARAMS['P']` (default **0.05 µM**) for every sample. The stored FDF
values are still used **once**: by the initial probe fit that discovers the
parameter order (see [Implementation notes](#implementation-notes)), so
keep them in a reasonable range for your signal scale.

## The template graph

Every output graph is a **clone of one template graph** that already exists
in the project. This keeps all runs visually consistent and lets you keep
custom page settings (size, position, window state) that the script does
not touch.

### Creating the template

Use a **single-layer** scatter graph — the script styles only the first
layer of the clone:

1. In the data worksheet, select the first block's X column and one Y column
2. `Plot > Basic 2D > Scatter`
3. Make sure the new graph window's **short name** is `Graph1`
   (rename it if Origin gave it a different number) — this is the default
   value of `TEMPLATE_GRAPH` at the top of the script

The template's own styling does not need to be prepared: the script
re-applies dataset colors, symbol shapes, fit-curve style, axis scaling,
axis-line positions, font sizes, and the legend on every clone. What is
**not** reset: the number of layers (extra layers keep their old plots),
axis-title *text* if the template already has one, and page-level settings
(size, position) — which is usually what you want.

The template graph itself is **never modified** by the script.

### What the script does with it

For each sample the script:

1. Activates the template and duplicates it (`win -d`), so the template
   stays untouched and the clone gets a fresh short name (Graph2, Graph3, ...)
2. Sets the clone's Long Name to `<sample>_graph`
3. Removes all plots from the clone's first layer and adds this sample's
   replicate datasets plus the fit curve
4. Applies the full styling: dataset colors (black/red/blue/green/orange,
   cycled), circle symbols, 3-pt black fit line, log X axis fixed to
   1e-4–100, top/right axis lines, font sizes, line-only legend

### Fallback behavior

- If `Graph1` does not exist, the script uses any existing output graph
  (Long Name ending in `_graph`) as the template — those are clones of the
  template anyway, and all plots/styles are re-applied on them
- If neither exists, the script stops with an error
- With `TEMPLATE_GRAPH = ''` the script always builds plain new graphs
  (Long Name still `<sample>_graph`) — no template and no fallback is used
- The template search is by **short name**; the Long Name of the template
  can be anything

## Requirements

- **Windows** with a commercial **OriginPro** installation
- **Origin 2024b or later** — the tick-label font feature requires
  `layer.*.label.pt` from Origin 2024b. Tested with **Origin 2025**; older
  versions have not been smoke-tested
- The embedded Python (Code Builder) and the `originpro` package that ships
  with it — **no numpy / pandas needed**, pure Python + `originpro` only
- The fitting function `FP_OneSiteBind` — included in this repository as
  [`FP_OneSiteBind.fdf`](FP_OneSiteBind.fdf), see
  [Installing the fitting function](#installing-the-fitting-function)
- An **English-language Origin interface** with decimal-point number
  format: the fit-report parser looks for the English column headings
  (*Value*, *Parameter*, *Error*) and converts numbers with Python's
  `float()`. Localized interfaces have not been tested

## Data layout

The active worksheet must contain repeating blocks of columns, one block
per sample:

```
| X | Y1 | Y2 | ... | X | Y1 | Y2 | ... |  (next sample)
```

- **X columns must be designated as X** in the worksheet
  (`Column > Set as X`). Each X-designated column starts a new sample
  block; all following columns with data (up to the next X column) are its
  replicates.
- **1 to 5 replicate Y columns per block** are supported (the default
  styling provides 5 colors; more replicates would cycle through them).
  All replicates of a block are concatenated into one fit.
- A **single sample** (one X-designated column in the first position,
  followed by any number of Y columns) is recognized automatically.
- If no X column is designated, the script falls back to a fixed periodic
  layout of `X + (GROUP_SIZE-1) Y` columns repeating (configurable via
  `GROUP_SIZE`, default 4 = X + 3 Y). This fallback cannot represent
  variable replicate counts — designate your X columns instead.
- The **sample name** is taken from the Long Name (or Comments) of the
  block's first Y column.

## Quick start

1. Download and extract the repository — keep `FP_OneSiteBind_Fit.py` and
   `FP_OneSiteBind.fdf` together
2. Install the FDF (see
   [Installing the fitting function](#installing-the-fitting-function)) and
   confirm it appears in the Fitting Function Organizer (F9)
3. **Work on a copy of your Origin project** — the script modifies the
   open project
4. Prepare the worksheet: designate the X columns, set the sample names in
   the Long Names of the first Y columns
5. Edit `FIXED_PARAMS['P']` at the top of the script — set it to the actual
   final labeled-probe concentration in µM
6. Create the single-layer template graph named `Graph1`
   (see [The template graph](#the-template-graph)), then **re-select the
   data worksheet**
7. `View > Code Builder` — open `FP_OneSiteBind_Fit.py` — press **F5**
8. Inspect the fit quality (R², Status column, residuals) and save the
   project manually if `SAVE_PROJECT` remains `False`

## What it produces

| Output | Description |
|--------|-------------|
| `FitCurves` sheet | fit curve for every sample as (X, Fit) column pairs |
| `FP_Summary` sheet | one row per sample: Sample \| Kd \| Kd SE \| Amax \| Amax SE \| Amin \| Amin SE \| R² \| Status |
| one graph per sample | clone of the template graph (see [The template graph](#the-template-graph)), Long Name `<sample>_graph` |

**Cleanup behavior.** Graphs and result sheets are **never deleted** —
repeated runs create new sheets (`FitCurves2`, `FP_Summary2`, ...) and new
graphs, so earlier results and their graph references survive. The fit
report sheets (`FitNL...`) created by the **current** run are deleted after
reading (set `KEEP_REPORTS = True` to keep them); pre-existing `FitNL...`
reports from other analyses are never touched. `CLEANUP_PREVIOUS = True`
additionally removes leftover `FPTmp` temp sheets from earlier interrupted
runs. Even so: work on a copy of your project.

Graph styling (all configurable at the top of the script):
the replicate datasets as circles in black/red/blue/green/orange, a 3-pt
black fit line, log-scale X axis fixed to 1e-4–100, axis lines on the top
and right, enlarged axis/tick/legend fonts, and a line-only legend.

**X-axis display range.** The graph X axis is fixed to `X_FROM`–`X_TO`
(1e-4–100 µM) on a logarithmic scale. This affects the *display only* —
all data points are used in the fit. Points outside this range are hidden
on the graph, and zero/negative X values cannot be shown on a logarithmic
axis at all.

## Configuration

All user settings are in the `USER SETTINGS` block at the top of the script:

| Setting | Default | Meaning |
|---------|---------|---------|
| `FIXED_PARAMS` | `{'P': 0.05}` | parameters held fixed during the fit — **edit `P` to your probe concentration** |
| `FIT_PARAMS` | `['Kd', 'Amax', 'Amin']` | fitted parameters |
| `DATA_COLORS` | `['Black', 'Red', 'Blue', 'Green', 'Orange']` | colors of the replicate Y datasets (cycled) |
| `FIT_COLOR` / `FIT_WIDTH_PT` | `'Black'` / `3` | fit curve style |
| `AXIS_TITLE_FONT` / `TICK_LABEL_FONT` / `LEGEND_FONT` | 31 / 22 / 22 | font sizes (pt) |
| `X_TITLE` / `Y_TITLE` | `'[PDZ1], µM'` / `'FA, mAu'` | axis titles — **replace `[PDZ1]` with your titrated protein**. Note: if the template graph already has axis-title text, the template text is kept and `X_TITLE`/`Y_TITLE` are ignored |
| `X_FROM` / `X_TO` | 0.0001 / 100 | fixed X axis display range |
| `TEMPLATE_GRAPH` | `'Graph1'` | graph used as the style template; `''` = plain new graphs |
| `GROUP_SIZE` | `4` | fallback layout (X + GROUP_SIZE−1 Y) only when no X column is designated |
| `N_FIT_POINTS` | `300` | points in each fit curve |
| `KEEP_REPORTS` | `False` | keep the fit report sheets created by this run |
| `CLEANUP_PREVIOUS` | `True` | remove leftover `FPTmp` temp sheets before running |
| `SAVE_PROJECT` | `False` | save the `.opju` when done |

## Statistics and interpretation

- Each sample is fitted as a **shared-parameter fit to the pooled
  observations**: all valid X/Y pairs of all replicates are concatenated
  and fitted once, unweighted. Replicates with more valid observations
  contribute more points to the fit. Kd values are *not* averaged across
  separate replicate fits.
- `Kd SE` is the **fit-derived parameter standard error** from the
  covariance matrix of the pooled fit — not the standard deviation of
  independently fitted replicate Kd values.
- `R²` is computed on the same pooled data: the residual sum of squares
  against the shared model, over the total sum of squares against the
  **overall mean** of all observations of the sample.
- The `Status` column reports fits where parameter standard errors were
  unavailable or R² could not be computed. High R² alone is not sufficient
  evidence of a reliable Kd — also inspect the residuals, the parameter
  uncertainties, and whether the titration actually approaches saturation.

## Example

[`examples/FA_example_data.csv`](examples/FA_example_data.csv) contains a
synthetic single-sample titration (13 ligand concentrations × 3 replicates,
Gaussian noise) generated from the model with **Kd = 0.5 µM, Amax = 250
mAu, Amin = 50 mAu, P = 0.05 µM**. Import it into a worksheet with the
first row as Long Names, designate column 1 as X, and run the script — you
should obtain approximately **Kd ≈ 0.49 µM, Amax ≈ 249 mAu, Amin ≈ 50 mAu,
R² ≈ 0.999**.

## Implementation notes

- Fit reports are parsed via LabTalk `GetReportTable(...)` into a temp file
  and read with plain Python; fit parameters are addressed **by index**,
  never by name, so special characters in FDF parameter names cannot break
  the script.
- The parameter order and the index of the fixed parameter are discovered
  by a **probe fit** on the first sample. The probe uses the same
  concatenated data preparation as the production fits and relies on the
  FDF's auto-initialization (the parameter order is not known yet, so
  initial values cannot be passed by name).
- Replicate Y columns are explicitly concatenated on a temporary worksheet
  before fitting — passing `(x, y1, y2, y3)` as a single range would make
  Origin interpret columns 3/4 as error columns.
- Fits are validated after parsing: fitted parameter values must be finite
  (otherwise the script stops with an error naming the sample), and missing
  standard errors are reported in the `Status` column instead of being
  silently ignored.
- Tick-label font size uses `layer.*.label.pt`, available since
  **Origin 2024b**.

## Troubleshooting

| Problem | Cause / fix |
|---------|-------------|
| "Probe fit failed: no report sheet was created" | `FP_OneSiteBind` is not installed — see [Installing the fitting function](#installing-the-fitting-function) |
| "No active worksheet" | Click the **data** worksheet before running (not a graph or report sheet) |
| "Template graph ... not found" | Create the `Graph1` template (single-layer scatter) or set `TEMPLATE_GRAPH = ''` |
| "Could not determine the sample layout" | Designate the X columns (`Column > Set as X`) or adjust `GROUP_SIZE` |
| "Fixed parameter(s) ... not found" | The FDF was changed — update `FIXED_PARAMS` or the FDF so both agree on the probe concentration |
| "parameter ... has no finite value" | The fit did not converge — check the data of the named sample and that the titration covers the Kd range |
| Wrong probe concentration in results | `FIXED_PARAMS['P']` was not edited — set it to your actual probe concentration |

## Citing

If you use this script, please cite it (see [`CITATION.cff`](CITATION.cff)
for the machine-readable record). If you use the `FP_OneSiteBind` fitting
model in published research, please also cite the publication in which the
model was used:

> Mymrikov E.V., Wirth C., Heinicke J.I., Goll J., Kern B.A., Steck C.,
> Iaroslavtceva A.K., Mühlethaler T., Köttgen A., Hunte C. (2026).
> **Molecular Determinants of Selective and High-affinity Binding of the
> Scaffold Protein PDZK1 to the Urate Transporter URAT1.**
> *Journal of Molecular Biology* 438(5), 169615.
> DOI: [10.1016/j.jmb.2025.169615](https://doi.org/10.1016/j.jmb.2025.169615)
> (PMID: [41453723](https://pubmed.ncbi.nlm.nih.gov/41453723/))

## License

[MIT](LICENSE)
