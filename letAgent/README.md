# letAgent: LET-aware multi-criteria plan search on matRad (MVP)

For one case, letAgent explores the trade-off between target dose, OAR high-LET
burden and robustness automatically, and returns a short list of
Pareto-optimal plans with an interactive report. It runs entirely on the
machine that holds the case data.

```
case JSON ──► letAgent_prepareCase.m   dose, LET×dose and dirty-dose influence matrices,
             (MATLAB/Octave, once)    all robustness scenarios, cached per beam set
                    │
                    ▼
python -m letagent search  ──►  propose recipes ──► letAgent_runBatch.m ──► metrics JSON
   (search + rules)        ◄──  Pareto front, rules  ◄───────────────────────┘
                    │
                    ▼
report.html · pareto.csv · runs.csv · agent_log.md · plans/<runId>_w.mat
```

* **Engine (MATLAB):** matRad with Remo's optimisation-quantity framework. It adds
  a `dirtyDose` quantity (`matRad/optimization/quantities/matRad_DirtyDose.m`,
  Seckler et al. 2026 Eq. 4) beside the existing `LETd` quantity. Note that
  `LETd` optimises LET×dose (`dij.mLETDose * w`); in the agent it is called
  `LETxDose`.
* **Agent (Python 3.8+, standard library only):** it never sees images or
  identifiers, only the metrics JSON.

## Privacy protocol (UHN: AI must never see patient data)

The work is split into two zones. Code moves one way only: from the development
zone to the hospital zone, through git.

| | Development zone (Mac, GitHub, AI assistants) | Hospital zone (UHN machine, no AI tools) |
|---|---|---|
| Holds | This code; matRad phantoms; synthetic cases; public data (TCIA) | DICOMs (even anonymised), imported `ct`/`cst` `.mat` files, patient case JSONs, `letAgent_runs`, reports, logs |
| Who works there | Developers and AI assistants | The REB-approved researcher only |
| Crosses over | Code (`git pull` on the UHN machine) | Nothing. A failure is reproduced on a synthetic case with the same settings, or described in your own words |

Anonymised head CTs still count as patient data here: a face can be rebuilt from a
head CT, and DICOM files can keep identifiers in private tags, UIDs and dates.

Guards in the code:
* **Data class.** Every case file sets `dataClass`: `patient`, `synthetic` or
  `public`. A case loaded from a `.mat` file is a patient case unless
  `dataClassNote` says where the data comes from.
* **Data location.** For patient cases, letAgent (Python and MATLAB) refuses to
  run unless:
  * the case file, `source.path` and `workDir` are all outside this repository
    and outside every folder in `config/ai_connected_folders.txt` (or in the
    `LETAGENT_AI_FOLDERS` environment variable);
  * `workDir` is set explicitly;
  * `reb.protocol` holds your REB number.
* **Commit guard.** Run `python3 letAgent/tools/git_guard.py --install` once per
  clone. It blocks committing:
  * DICOM files;
  * `letAgent_runs*` outputs;
  * `.mat`, `.raw`, `.mhd` and `.nii` files outside the machine and phantom folders;
  * patient case files.

  It is also in `.pre-commit-config.yaml`. On Windows, git hooks run in Git Bash;
  if `python3` is not found there, edit `.git/hooks/pre-commit` to call `python`.
* **Aliases.** The alias must be a study code (`CNS_03`). Aliases that look like
  an MRN or a date are refused.
* **Scrubbing.** The agent reads metrics only. File paths, long numbers and dates
  are scrubbed from its logs.
* **Banner.** Patient reports carry a banner with the REB protocol number.

Things the guards cannot catch, so they need care:
* matRad's DICOM import and MCsquare/FRED write into `userdata/` inside this
  repository: `userdata/patients`, and `userdata/MCsquare/MC2patientCT.*`.
  * On the UHN machine that is acceptable, because no AI tool is connected there.
    Still, point the DICOM import's save dialog to hospital storage.
  * On any machine linked to an AI assistant, never import or recalculate patient
    data.
* Do not connect a hospital-storage folder to an AI tool. When you connect a new
  folder to one, add it to `config/ai_connected_folders.txt`.

## Quick start (phantom, about 15–30 min)

In MATLAB, from the repository root:

```matlab
matRad_rc; addpath('letAgent/matlab');
letAgent_prepareCase('letAgent/cases/TG119_DWA.json');   % optional; search does it too
```

In a terminal:

```bash
cd letAgent/agent
python -m letagent search --case ../cases/TG119_DWA.json --engine matlab --budget 30 --workers 2
```

Then open `letAgent_runs/TG119_DWA/report.html`.

On Windows, if `matlab` is not on the PATH, add
`--engine-cmd "C:/Program Files/MATLAB/R2024b/bin/matlab.exe"`. GNU Octave
also works (`--engine octave`); it is slower, but it is how this code was
tested.

## Running a UHN patient

1. Import the patient into matRad on the UHN machine. Use `matRad_importDicom`
   or the GUI, then `save('D:/letAgent_cases/CNS_01.mat','ct','cst')`. Use
   hospital storage outside this repository, not `userdata/patients`.
2. Copy `cases/UHN_CNS_TEMPLATE.json` to hospital storage outside this repository. Set:
   * `reb.protocol` to your REB number;
   * `alias`, `source.path` and `workDir`;
   * structure names exactly as they appear in `cst(:,2)`;
   * the prescription;
   * OAR limits (`hardLimit: true` makes a limit a pass/fail check);
   * which OARs get LET objectives (`let: true`);
   * the beam sets to compare.
3. `python -m letagent search --case D:/.../CNS_01.json --engine matlab --workers 2`
4. Review `report.html`. Move the sliders to your clinical limits, and the
   page names the Pareto plan with the lowest OAR LET burden that meets them.
   That plan's spot weights are in `plans/<runId>_w.mat`. Recalculate it with
   MCsquare before drawing conclusions.

Sizing: one process holds one beam set's influence matrices in RAM, for all
robustness scenarios. With `wcScen` (9 scenarios) at 3 mm, start with
`--workers 1` and check the memory before adding workers.

## Case file reference

| Field | Meaning |
|---|---|
| `alias` | Study code: the only case label that appears in outputs |
| `source` | `matRadPhantom` (`name`), `matFile` (`path` to a .mat with `ct`, `cst`), or `synthetic` (`cnsStress`, with `gapMm` etc.) |
| `machine` | matRad machine, e.g. `DWA_proton_v2`, `IBA_cyclotron`, `Generic` |
| `prescription` | `totalDoseGy`, `numOfFractions` |
| `structures.targets[]` | `name`, `doseGy`, `penalty` (squared deviation) |
| `structures.oars[]` | `name`, `maxDoseGy` + `penalty` (squared overdosing), `meanDosePenalty`, `let` (gets the LET objective), `hardLimit` |
| `beamSets[]` | `id`, `gantryAngles`, `couchAngles` |
| `stf`, `doseGridMm` | Spot spacing and dose grid |
| `robustness` | `nomScen`, or `wcScen` with `shiftSDmm`, `rangeRelSDpct`, `rangeAbsSDmm` |
| `metrics` | `dirtyDoseLETThreshold` (keV/µm), `letdMinDoseGy`, `letVolume.minDoseGy` / `minLETd` (for V(D & LETd)) |
| `engine.maxIter` | IPOPT iterations per plan |

## What the agent does

Objectives, all minimised:
* **coverage loss:** 100·(1 − D95/Rx), nominal;
* **OAR LET burden:** by default dirty-dose D2% summed over the `let` OARs;
  `--let-metric VDLpct` or `LETdP98` switches it;
* **robust coverage loss:** the same coverage loss in the worst-case scenario,
  when the case has scenarios.

Each search parameter is a recipe:
* beam set;
* LET formulation (none / dirty dose / LET×dose);
* LET penalty;
* OAR-dose penalty scale;
* robustness (none / COWC);
* coverage as an objective or as a hard DVH constraint.

Rounds:
1. **Anchor round:** a dose-only plan, plus one plan per LET formulation, for
   each beam set.
2. **Later rounds:** half refinements of Pareto plans (perturbing them, and
   filling the largest gaps along the front), half exploration.
3. **Stop** when the budget is used or the hypervolume stops improving.

Between rounds the rules (`agent/letagent/rules.py`) read the results, change
the search space and explain why in the log:

| Rule | Trigger | Action |
|---|---|---|
| coverage-erosion | LET plan loses > 2 points of coverage vs dose-only | Cap the LET penalty; add plans with coverage as a hard DVH constraint |
| let-unneeded | Dose-only plan already has negligible OAR LET burden | Stop LET exploration for that beam set |
| let-inactive | 10× penalty range changes the LET burden < 10% | Widen the penalty range once, then drop the combination |
| formulation | One LET formulation never reaches the front | Drop it and focus on the other |
| robustness | Front plans lose much more coverage in the worst case | Add robust (COWC) plans |
| beam-set | A beam set never reaches the front | Stop spending runs on it |

## Limitations of this MVP

* Pencil-beam LET only; verify chosen plans with MCsquare.
* Robustness is handled with matRad's scenario-based worst case, which is slow.
  The scenario-free robust MCO of Cristoforetti et al. 2026 is the planned
  replacement.
* The rules are a first version; adjust the thresholds in
  `agent/letagent/agent.py` (`DEFAULT_SETTINGS`).
* MATLAB is still needed (or Octave). A pure-Python engine (pyRadPlan) is the
  next step toward a commercial product.

## Public and synthetic test data (no patient data needed)

* matRad phantoms: TG119, PROSTATE, HEAD_AND_NECK, LIVER (`matRadPhantom`).
* `CNS_stress_synthetic.json`: a target abutting a brainstem-like box.
  Vary `gapMm` to make the trade-off easier or harder.
* TCIA GLIS-RT: 230 glioma patients with planning CT and targets. Seckler et
  al. 2026 used its anaplastic astrocytoma case with a brainstem. Download it
  and import it with `matRad_importDicom`; it can then be shared freely
  under the TCIA terms.
