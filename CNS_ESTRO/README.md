# CNS_ESTRO: CNS test-patient planning for the ESTRO 2027 abstract

Plans a CNS case in matRad with every **machine × plan variant** of the study, and
evaluates **target DVH** (physical, RBE 1.1, McNamara) and **OAR dose–LET metrics**.
The same code runs on (1) public or synthetic data on any computer and (2) UHN
patients on the hospital machine only.

Requires the `let_opt` branch (Remo's optimisation quantities: `LETd` = LET × dose).
Dirty-dose plans also need `matRad_DirtyDose.m` (letAgent); without it they are skipped.

```matlab
% from the matRadDWA root
matRad_rc; addpath('CNS_ESTRO');
run('CNS_ESTRO/run_CNS_ESTRO.m')          % synthetic smoke test, then edit the case list
```

## Start here: training on matRad's HEAD_AND_NECK (`matRad_example5_DWA_HN_training.m`)

Rehearses the UHN workflow on public data, cell by cell:
1. a 9-field photon IMRT plan on HEAD_AND_NECK, exported to DICOM (CT + RTSTRUCT +
   RTDOSE), as the stand-in for the clinical RP1/RD1;
2. import of that folder with `cns_importDicom` (the step you run on the UHN machine);
3. canonical structures (`cases/HN_matRad_training.json`: CTV = GTV, PTV = PTV70,
   70 Gy(RBE) / 35 fx);
4. a DWA dose-only plan (example 5 flow) and **GUI checkpoint**: total RBE 1.1 dose,
   McNamara dose, LET and the photon dose are added to `resultGUI` for the GUI dropdown;
5. the same with a LET × dose objective, **GUI checkpoint** with both plans and their difference;
6. a 2 × 3 slice figure for slides (`cns_plotSlice`: CT + wash + contours);
7. the full study on this case (both machines, all plans, photon reference) and the cohort table.

`S.fast = true` uses 5 mm spots and a 4 mm grid. matRad blocks RTPLAN export on purpose, so
the test folder has no RP file; for the UHN folder `cns_importDicom` imports RP1 if it can
and falls back to CT + RS + RD otherwise (only the dose is needed).

## Clinical photon plan (RD1)

`cns_importDicom` returns the RTDOSE as `refDose` (on the CT grid) and saves it with
`ct` and `cst`. `cns_runCase` then evaluates it with the same metrics as machine
`photon`, plan `clinical` (RBE = 1, no LET), overlays it in black on the proton DVHs, and the
cohort table gets a "photon clinical" column. A clinical RD1 is usually a PLAN (total) dose;
matRad divides it by the RTPLAN fractions when it imports RP1, and `refDose.toTotal` undoes
this. The run prints the median PTV dose of the photon plan: check it against the
prescription.

## UHN CNS ependymoma case: three scripts (v1)

Run on the UHN PC only, in this order, each with your hospital-storage copy of
`cases/UHN_CNS_ependymoma_TEMPLATE_v1.json` and `study_UHN_CNS_v1.json`:

| Script | What it does | Output (under `<workDir>/<alias>/`) |
|---|---|---|
| `UHN_CNS_01_photonPlan_DVH_v1.m` | imports CT + RS1 + RP1 + RD1 once; maps structures; evaluates the delivered photon plan (RBE 1); constraint check; DVH figure; opens matRadGUI | `photon_v1/` |
| `UHN_CNS_02_protonPB_doseOnly_3machines_v1.m` | beam-geometry figure for the candidate beam sets; optional beam screen (Generic, coarse); dose-only plans for DWA_proton_v2, Generic, IBA_cyclotron (60 mm range shifter), normalised to the photon plan's PTV D95; comparison table + DVH/LVH | `dij_*.mat`, `plans/`, `doseOnly_v1/` |
| `UHN_CNS_03_protonPB_LETxD_weightSearch_3machines_v1.m` | per machine: gradient descent on log10(w) of the LET x dose weight from w = 0.1 (`cns_letWeightSearch`), best plan vs dose-only vs photon; search plot, DVH/LVH, table | `LETxD_search_v1/` |

Set `setenv('MATRAD_ROOT', 'X:\matRadDWA_LET')` once per session, then run each script.

**Beams.** matRad (LPS): gantry 0 = from anterior, 90 = from the patient's LEFT (right side of the
axial image), couch 270 + gantry 90 = vertex. "2 o'clock" on the axial image = gantry 60.
Candidates in `study.beamCandidates`: A 2 + 4 o'clock (60/120), B 2 o'clock + vertex, C left lateral
+ vertex, D 2 + 10 o'clock (crosses the right orbit; for contrast).

**Weight search** (`study.letSearch`):
J(w) = mean over optic nerves, chiasm and brainstem of [D x LETd D2](w) / [D x LETd D2](dose-only)
+ lambda x (squared coverage violations in %). Coverage: CTV D98 >= 95 % Rx, D98 within 1 % of
dose-only, D2 at most 2 % above dose-only. Steepest descent with central finite differences in
log10(w), Armijo backtracking, w in [0.001, 3], at most 14 plans per machine; the best of all
evaluated weights is reported. The LET objective uses a nonzero reference: the 5th percentile of
LET x dose in the CTV of the dose-only plan (total course; matRad divides every objective dose
level by numOfFractions, `matRad_initOptimization` line 79, and `matRad_LETd` adds no transform).
J is not guaranteed smooth or convex: the search finds a local optimum at the resolution of the
step, and every evaluated weight is kept in `search_<machine>_v1.csv`.

**Statistics.** One patient: descriptive values only, no tests (Wilcoxon only for n >= 5, in
`cns_cohortSummary`).

**Figures.** `cns_journalFigure` / `cns_saveJournalFigure`: Arial 7-8 pt, 17.5 cm wide,
Okabe-Ito colours (photon black, DWA blue, Generic vermillion, IBA bluish green), vector PDF +
300 dpi PNG + csv of every plotted value.

## What one run does (`cns_runCase`)

1. **Privacy check** (`cns_checkPrivacy`): alias must be a study code; a patient case
   refuses to run unless the case file, data and `workDir` are outside the repository
   and every folder in `ai_connected_folders.txt`, and `reb.protocol` is set.
2. **Structures** (`cns_prepareStructures`): canonical names from regexes or
   `structures.map`; PTV = CTV + 3 mm if absent; Brain from HU if absent (GLIS-RT);
   derived `Brain_healthy` (brain − CTV), `BrainMinusPTV`, `BrainRing_0_15`
   (LET-objective shell), `Ring_0_5/5_10/10_20` (displacement check).
3. **Per machine** (`DWA_proton_v2`, `Generic`): stf + dose and LET × dose dij
   (HongPB, double Gaussian, cached as `dij_<machine>.mat`).
4. **Per plan** (`study.plans`): `doseOnly`, `LETxD_*` (SquaredOverdosing(p, 0) on
   LET × dose in brainstem, chiasm, optic nerves and the 0–15 mm brain shell),
   `dirty_*`. Target penalties as in the water phantom (p_cov = 1000).
5. **Normalisation**: PTV D95 (RBE 1.1) = 95 % Rx, as in the water-phantom runs.
6. **Evaluation** (`cns_evaluatePlan`), all total doses:

| Structure | Quantity | Metrics |
|---|---|---|
| CTV, PTV, GTV | phys, RBE11, McN (α/β 10) | D98, D95, D50, D2, Dmean, V95 %, V107 %; DVH; dose-weighted LETd |
| CTV (if scenarios) | RBE11 | worst-case D95, D98 |
| each OAR + helper | RBE11, McN (α/β 2) | Dmean, D2, D0.03cc, Dmax, V20Gy, V40Gy; DVH |
| | LETd (dose ≥ 1 Gy(RBE)) | dose-weighted mean, P50, P98; LVH |
| | DxLETd | RBE11·LETd smoothed with a 6 mm Gaussian, voxels ≥ 5 Gy(RBE): mean, D2 (Heinzelmann 2026) |
| | hotSpot | V(D ≥ X & LETd ≥ Y) in cc and %: brainstem 50 Gy & 2.8, healthy brain 43 Gy & 4.6, all OARs 20 & 4 and 30 & 2.5 |
| | dirty | dirty-dose Dmean, D2 (if the dij has it) |

Outputs in `<workDir>/<alias>/`: `metrics.csv` (long table), `dvh.csv`, `lvh.csv`,
`run_info.json`, `plans/<machine>__<plan>_w.mat`, `figures/` (target and OAR DVH, OAR
LVH, hot-spot bars; dose/LET slices for non-patient cases only).

## Cohort table (`cns_cohortSummary`)

Picks per case and machine **LETselected**: the LET plan with the lowest brainstem
D×LETd D2 among those keeping CTV D98 ≥ 95 % and D2 ≤ 107 % (`study.selection`). Then
median [range] per machine × plan, and paired differences vs `doseOnly` (same machine)
and vs `Generic` (same plan), with a signed-rank p for n ≥ 5. `summary_abstract.md` is
the table for the abstract.

## Data

**Synthetic** (`cases/CNS_synthetic.json`): head with skull, CTV 3 mm from a
brainstem-like cylinder; vary `gapMm`.

**Public (TCIA)**, import once with `cns_importDicom(dicomDir, outMat)`:
- GLIS-RT, doi:10.7937/TCIA.T905-ZQ20: 230 glioma cases, CT + RTSTRUCT with GTV, CTV,
  brainstem, chiasm, optic nerves, eyes, cochleae, lacrimal glands; no brain contour
  (the HU proxy is used). Registration needed.
- Burdenko-GBM-Progression, doi:10.7937/E1QP-D183: 180 GBM cases with CTV, PTV, brain,
  brainstem, chiasm, optics, plus RTPLAN/RTDOSE. Needs a TCIA restricted licence
  agreement, so check its terms before sharing results or figures.

**UHN patients** (`cases/CNS_UHN_TEMPLATE.json`): copy the template to hospital storage,
one file per patient; import with `cns_importDicom` on the UHN machine; set
`reb.protocol`, `structures.map`, beams and prescription. Everything (case files,
`.mat`, `workDir`, figures) stays on hospital storage. No AI tool on that machine.

## Tested

Octave 8.4 on the pushed `let_opt` branch (3e604042), synthetic case, Generic and
DWA_proton_v2, dose-only and LETxD_0p3, at a 4 mm grid with 6 mm spots and 120 IPOPT
iterations (smoke-test settings, not the study settings). Everything ran end to end:
metrics, DVH/LVH tables, figures, cohort table, LETselected, privacy guards. Dirty-dose
plans were skipped because `matRad_DirtyDose` is not on GitHub. Not yet run in MATLAB.

**Memory.** The first attempt at the study settings (3 mm grid, DWA 3 mm spots, dose + LET
dij) ran out of 7 GB. For a real head CT, use the UHN workstation and check RAM, or do a
first pass with `doseCalc.doseGridMm = 4`.

Also tested in Octave on HEAD_AND_NECK: photon plan, import of CT + RTSTRUCT + a
PLAN-summed RTDOSE (made with pydicom; the dose came back within 0.001 Gy), structure
mapping, Generic dose-only and LET × dose plans, the photon reference rows, overlaid DVHs,
slice figures and the cohort table. Octave cannot *write* RTSTRUCT/RTDOSE and reads
multi-frame RTDOSE differently, so sections 1–2 of the training script need MATLAB. The
DWA dij for this phantom did not fit in 7 GB of RAM; use the Mac or PC.

## Relation to letAgent

Same privacy rules and alias check as letAgent, and the same optimisation quantities, but
standalone: letAgent is not on the GitHub branch. This is the fixed-design study for the
abstract; letAgent is the multicriteria search on top.

## Things to check before the abstract

- **LET source**: DWA/IBA LET tables are the Generic-range-matched stopgap; say so in
  Methods. Pencil-beam LET; recompute the selected plans with MCsquare.
- **Hot-spot thresholds** (50 Gy & 2.8, 43 Gy & 4.6) are evaluated on RBE 1.1 total dose.
  Confirm this is how the source papers define them.
- **6 mm Gaussian**: implemented as σ = 6 mm. Check Heinzelmann 2026 for σ vs FWHM.
- **LET penalties** (0.03–1) were set from the water phantom (0.1–0.3 usable, ≥ 10 not).
  Check the LETselected plan is not always the smallest or largest penalty; if it is,
  widen the sweep.
- **OAR objectives** in `study_ESTRO2027.json` are generic; replace with the UHN CNS protocol.
- Beam angles are per case (`beams` in the case file).

## IPOPT on Mac ("IPOPT mex interface not available for MATLAB")
All optimisations now go through `cns_optimize.m`. Before each plan it tests IPOPT itself; if the
test fails it prints the real reason, unloads the mex (`clear ipopt`) and tests again.
If it still fails, the run stops with the diagnosis (MATLAB architecture, mex extension, which
ipopt file was found). Common causes and fixes:
- **Apple Silicon MATLAB (`computer('arch')` = `maca64`)**: matRad ships only `ipopt.mexmaci64`
  (Intel). Either run the Intel MATLAB build (Rosetta), compile IPOPT for maca64, or run on the PC.
- **Intermittent failure on Intel/Rosetta MATLAB** (the mex left in a bad state, e.g. after Ctrl-C
  or a crashed run): the automatic `clear ipopt` usually fixes it; otherwise restart MATLAB.
- **Last resort**: add `"fallbackOptimizer": "fmincon"` under `"optimisation"` in the study JSON
  (needs the Optimization Toolbox). Every plan that used fmincon is flagged in the log, in
  `run_info.json` (`optimizer`) and in the weight-search trace. Do not mix optimisers in the final
  machine comparison; rerun those plans with IPOPT.

Figure saving on Mac: hidden figures sometimes report a zero-size canvas ("The output width and
height must be greater than 0"). The savers now resize the figure, try `exportgraphics`, and finally
redraw it visibly before giving up. The plan results are not affected by this message.

## UHN_CNS_02b_beamScreen_slices_v2.m (beam screen, TPS-style slices)
Run after the beam screen (UHN_CNS_02 with `S.runBeamScreen = true`). Nothing is re-optimised: it
reads each `<alias>_scr<k>` folder, normalises the plan exactly as the screen table did, and draws
the photon plan + sets A-D on the same slice (axial, sagittal, coronal), centred on the CTV and on
the chiasm. Isodoses 50/95/107/120 % (120 % shows the hot spots). Output:
`beamScreen_v1/figures_slices_v2/fig_screen_<plane>_<centre>_v2.pdf/.png` (+ `_slices.mat`).

## Eye-sparing beam screen (study_UHN_CNS_beamScreen_v2.json)
Same as study_UHN_CNS_v1 but with beam sets that never enter through the orbit (source directions
checked with matRad_getRotationMatrix): C left lateral + vertex; E lateral 10 deg posterior +
superior-left oblique (couch 300); F superior-left 45 deg (couch 315) + vertex; G lateral 10 deg
posterior + vertex. Results go to `beamScreen_v2/` (study field `beamScreen.folder`), so the first
screen is kept. Run UHN_CNS_02 (runBeamScreen = true) and UHN_CNS_02b with S.studyFile set to this file.
**`beams.setId` is only a label**: copy the set's gantryAngles and couchAngles into the case file.
cns_buildPln now stops if setId names a set of the study but the angles differ.
