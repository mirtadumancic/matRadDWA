**DWA machine, final version**

Part 3 — Double-Gaussian lateral model and interpolated energy layers

*Mirta Dumančić · working notes · 30 September 2026 · for the group meeting of 1 October*

**Scope.** Part 1 built the first DWA machine for matRad from the TOPAS/RayStation exports; Part 2 compared it with matRad’s Generic proton machine in a water phantom. This part records the version of the machine we intend to use from now on, DWA_proton_doubleGauss_interp, and the three changes that define it: a double-Gaussian lateral model, interpolated energy layers, and a unit fix in the energy spectrum. It also corrects three statements in Part 2.

Summary for the meeting

1.  **Lateral model: double Gaussian.** The DWA spot is a narrow core (σ ≈ 1.2–2.3 mm at isocentre) on a broad halo (σ ≈ 19–23 mm) that carries 38–52% of the fluence above 50 MeV. The old single Gaussian was fitted to the RMS of the whole profile, which oversized the spot about 4× and erased the beam’s focusing.

2.  **Energy layers: 58 → 201.** The native TOPAS layers are 1.4–2.8 Bragg-peak widths apart, so no weighting can make a flat SOBP; that is the ripple in Remo’s SOBP. With 143 interpolated layers the spacing is ≤ 0.5 peak widths (Generic: 0.44) and the ripple of an ideal 1D SOBP falls from 10–29% to below 0.7%.

3.  **The interpolation reproduces held-out TOPAS layers** to within 0.2 mm in range (leave-one-out, with twice the real gap).

4.  **Energy-spread unit fix.** matRad stores the spread in percent; our script stored MeV. It only affects Monte Carlo recalculation.

5.  **Validated in matRad (30 Sept).** Single pencil beams at 20, 122 and 181 MeV reproduce the TOPAS depth curves (Bragg peak within 1 mm) and deliver the halo exactly as modelled (matRad/model = 1.00), after fixing a bug that made matRad silently ignore the halo (Section 7).

6.  **Still missing:** absolute dose normalisation, LET data, and the in-water nuclear halo — all need new TOPAS output from Alaina.

1\. Versions and file names

Every change to the machine gets a new commissioning script and a new machine name, so older machines stay on disk and any result can be traced to, and re-run with, the machine that produced it.

| **Machine name (pln.machine)** | **Commissioning script**                         | **What it is**                                                                                                                       |
|--------------------------------|--------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------|
| DWA_proton                     | matRad_commissionDWAmachine.m                    | Single Gaussian (RMS of the air profile), 58 layers. Used in Parts 1–2 and in Remo’s SOBP.                                           |
| DWA_proton_doubleGauss         | matRad_commissionDWAmachine_doubleGauss.m        | Double Gaussian, 58 layers.                                                                                                          |
| DWA_proton_doubleGauss_interp  | matRad_commissionDWAmachine_doubleGauss_interp.m | Double Gaussian, 201 layers (58 native + 143 interpolated), energy spread in percent, no sigma field (Section 7.2). Current version. |

Each machine is saved as protons\_\<name\>.mat. The four analysis scripts now take the machine from one line, dwaMachine = '…', and put that name into every output folder and file:

| **Script**                        | **Output (for the new machine)**                                                                                       |
|-----------------------------------|------------------------------------------------------------------------------------------------------------------------|
| …\_singleBeam.m                   | singleBeam_DWA_proton_doubleGauss_interp/ DWA_proton_doubleGauss_interp_PDD, …\_beamseye_slice, …\_depth_slice         |
| matRad_DWAvsGeneric_BOXPHANTOM.m  | DWAvsGeneric_150MeV_DWA_proton_doubleGauss_interp/ …\_vsGeneric_PDD, \_distalFalloff, \_lateralPeak, \_lateralPeak_log |
| matRad_DWAvsGeneric_energySweep.m | DWAvsGeneric_sweep_DWA_proton_doubleGauss_interp/ per-energy panels, 5 summaries, metrics CSV                          |
| matRad_DWAvsGeneric_fineFalloff.m | DWAvsGeneric_fineFalloff_DWA_proton_doubleGauss_interp/ falloff figure, edge figure, metrics CSV                       |

The two sweep scripts have a new switch, dwaLayers = 'native' (default), which runs only the 58 TOPAS layers so that results line up one-to-one with the older machines; 'all' includes the interpolated ones. The outputs from Parts 1–2 were written before this scheme, to untagged folders such as DWAvsGeneric_sweep/. They were produced with DWA_proton and should be renamed by hand, e.g. to DWAvsGeneric_sweep_DWA_proton/.

2\. New information this week

From Alaina (TOPAS / beam design)

- x_iso_1rms_radius_mm in the beam summary is the 1σ spot radius in vacuum — a design target of 1.000 mm at every energy, not a measured quantity.

- The beam crosses 55 cm of air to isocentre (15 cm from the exit window to the snout, 40 cm from the snout to isocentre). The only material in the path is the IC64-6 monitor chamber; the snout bore clears the beam.

- Spot profiles were scored in air; Bragg peaks were scored in water with the surface at isocentre.

- She leans towards the double Gaussian, for consistency with Generic, pending Remo’s view.

From Remo (DKFZ)

- Generic is not a real machine: its data are TOPAS-generated, with optics loosely modelled on a synchrotron. (This corrects Part 2 — see Section 7.)

- The machine file “looks good”. Two issues: the DWA energy spacing is coarse and produces SOBP artifacts, and the depth-dose normalisation (Z) is still open.

- Three routes to LET-based planning: a colleague’s published implementation (diverged from current matRad), his own custom-quantity implementation (dose + LET simultaneously, not yet in official matRad), or pyRadPlan (Python, native custom quantities, same machine format).

3\. Lateral model: double Gaussian

**Why a single Gaussian fails.** Each measured air profile is a narrow core on a broad, low-amplitude halo. The previous script summarised it by its second moment (RMS), which is dominated by the halo: at 149 MeV the RMS is 4.99 mm, the core 1.27 mm. Worse, across the five scoring planes the core varies by 85% (a clear waist near isocentre) while the RMS varies by 6%, so matRad was handed an almost parallel beam.

**Why the halo cannot be dropped either.** On a line profile the halo looks small, but a Gaussian’s 2D integral scales with a·σ², so a halo about 17× wider than the core carries far more fluence than the line profile suggests. The correct weight is the 2D share a₂σ₂² / (a₁σ₁² + a₂σ₂²); the 1D share a·σ understates it about 8×. A model-free radial integral of the measured profile (∫2πr·f(r) dr) agrees with the fitted weights to 0.02. The halo is physical, not a scoring floor: between r = 20 and 60 mm the fluence falls by a factor 0.033, as a σ ≈ 20 mm Gaussian would, where a flat background would give 1.0.

**Fit.** Stage 1 fits a parabola to log(fluence) over the points above 10% of the maximum (a Gaussian is a parabola in log space), giving the core. Stage 2 subtracts the core and fits the residual beyond 3 core-σ as a second Gaussian on the same axis. All 580 profiles fit without failure.

| **Energy** | **σ core at iso** | **Highland prediction** | **σ halo at iso** | **Halo weight (2D)** |
|------------|-------------------|-------------------------|-------------------|----------------------|
| 50 MeV     | 2.27 mm           | 2.08 mm                 | 19.3 mm           | 0.38                 |
| 92 MeV     | 1.50 mm           | 1.42 mm                 | 21.9 mm           | 0.51                 |
| 149 MeV    | 1.27 mm           | 1.19 mm                 | 22.6 mm           | 0.50                 |
| 219 MeV    | 1.17 mm           | 1.10 mm                 | 23.3 mm           | 0.46                 |

**Independent check on the core.** Starting from Alaina’s 1.0 mm vacuum spot and adding multiple Coulomb scattering (Highland) in the 55 cm of air and the IC64-6 foils (124 µm polyimide, 0.8 µm Al, 0.2 µm Au, 50.6 cm lever arm) predicts the fitted core to 6–9%. The core is therefore explained by known scattering; the halo is the part that is not.

**Mapping onto matRad.** matRad’s double-Gaussian engine adds the same initial beam width to both components, σ_narrow² = σ₁² + σ_ini² and σ_broad² = σ₂² + σ_ini², with σ_ini taken from initFocus. So:

initFocus.sigma = core sigma vs distance (5 planes) % real source size and focusing

sigma1 = in-water MCS (matRad's Highland form)

sigma2 = sqrt(sigmaHalo_iso^2 - sigmaCore_iso^2 + sigma1.^2)

weight = halo 2D weight at isocentre (constant with depth)

(no 'sigma' field - its presence makes matRad use a single Gaussian, Section 7.2)

4\. Interpolated energy layers

4.1 The problem: layers further apart than a peak is wide

A flat SOBP needs neighbouring Bragg peaks to overlap. The DWA’s tiny energy spread (0.044–0.060%) makes its peaks exceptionally narrow — the sharp distal falloff found in Part 2 — so it needs closely spaced layers. The native TOPAS layers are 3 MeV apart below 98 MeV and 4–5 MeV above, which is 1.4–2.8 times the 80%–80% peak width everywhere (Figure 1). Generic sits at 0.44.

| **150–210 mm (Remo’s SOBP)** | **layers** | **spacing** | **80–80 peak width** | **spacing / width** |
|------------------------------|------------|-------------|----------------------|---------------------|
| DWA, native                  | 7          | 7–10 mm     | 4.5–5.6 mm           | 1.5–2.1             |
| Generic                      | 21         | 3 mm        | 6.8 mm               | 0.44                |



**Figure 1.** *Range step between consecutive layers divided by the local 80–80 Bragg-peak width. Black: the 58 native TOPAS layers. Red: after interpolation. Below R80 ≈ 40 mm the red curve rises because of the 1 mm minimum step (Section 4.2); the peaks there are only one or two 1 mm bins wide in the source data.*

4.2 The rule

New layers are inserted between each pair of native layers so that the range step is at most interpSpacingFraction = 0.5 times the local 80–80 width, but never finer than interpMinStep_mm = 1 mm, the sampling of the TOPAS depth curves. The width comes from a smooth fit, W80 = 0.0396 · R80^0.935 (mm), made where the peaks are resolved (R80 \> 60 mm), because below that the measured widths are only one or two bins and jump around. The result is 201 layers: 58 native and 143 interpolated.

The fraction 0.5 was chosen with a 1D test: optimise non-negative layer weights for a flat SOBP using the depth curves directly. Spacing 1.0 still left 1–2% ripple; 0.5 gave under 0.7% at every depth tested, at a spacing close to Generic’s. Refining the floor from 1 mm to 0.5 mm changed nothing measurable.

4.3 How a layer is interpolated

**Depth dose — range-scaled blending.** For a target range R between native ranges Rₐ \< R \< R_b, each neighbour is stretched in depth so that its range lands on R, and the two are blended:

t = (R - Ra) / (Rb - Ra)

D(z) = (1-t) \* Da(z \* Ra/R) + t \* Db(z \* Rb/R)

Stretching aligns the two peaks before they are averaged, so the blend is a single peak rather than a double hump, and the entrance region scales consistently. Plain linear blending without alignment would smear the peak.

**Energy.** From a local power law R = a·E^p through the two neighbours, rounded to 0.01 MeV like the native list (matRad matches energies by equality).

**Lateral parameters and energy spread.** Core σ at each plane, halo width and weight, and dE/E are blended linearly in t; the in-water MCS term is recomputed for the new energy exactly as for native layers. These parameters vary smoothly with energy (Figure 4), so linear blending is adequate.

**Peak position.** peakPos is now taken to sub-millimetre precision (a parabola through the maximum and its two neighbours) for every layer. On the raw 1 mm grid, 13 interpolated layers would have shared a peakPos with their neighbour, and matRad selects layers for a target by peakPos. Native peak positions move by at most 0.5 mm, within one sampling bin.

**Provenance.** Every layer carries machine.data(i).interpolated (true/false), and machine.meta.interpolation records the method, the settings and the native energy list.

4.4 Validation

**Leave-one-out.** Each native layer (except the two at the ends) was deleted and rebuilt from its neighbours — which are two steps apart, twice the gap the interpolation actually bridges, so this is a pessimistic test (Figure 2). Across all 56: R80 within 0.19 mm (median 0.03 mm); 80–80 width +0.14 mm median, i.e. slightly broadened, as expected from averaging two peaks; distal 80–20 falloff +0.06 mm median. Above 100 MeV the entrance plateau agrees within 1.7% and the width within 0.23 mm. Below 100 MeV agreement is looser (plateau differences up to 7%, widths up to 0.7 mm), because those peaks are only 1–2 mm wide and sampled at 1 mm, so the source curves themselves are coarse there.



**Figure 2.** *Leave-one-out test. Left and centre: a native TOPAS curve (black), the two neighbours used to rebuild it (grey) and the rebuilt curve (red dashed). Right: rebuilt minus TOPAS for R80 and the 80–80 width, for all 56 held-out layers.*

**SOBP.** Figure 3 repeats the 1D SOBP optimisation with both layer sets. With the native layers the flattest achievable 150–210 mm SOBP has a 15% ripple with a period of about 9–10 mm — the pattern in Remo’s figure. Its distal edge also falls short (R80 209.3 mm, against 213.7 mm with interpolated layers), because the deepest native layer that fits stops short and the next one overshoots; this matches the ~6 mm shallower edge Remo saw. With 201 layers the ripple is 0.1%, and below 0.7% at every depth tested.



**Figure 3.** *1D SOBP optimisation on the laterally integrated depth curves. Left: 150–210 mm target (grey) with native (black) and interpolated (red) layers. Right: ripple (max − min over the target, excluding 1 mm at each edge) for five target depths.*

**What this test does and does not show.** It uses laterally integrated depth doses and an unconstrained optimiser, i.e. an ideal broad field. It shows that the layer set can in principle produce a flat SOBP; it is not a matRad plan. The real check is Remo’s SOBP re-run with the new machine.



**Figure 4.** *Lateral parameters for all 201 layers: core σ at isocentre (left) and 2D halo weight (right). Black: fitted to the native profiles. Red: interpolated.*

5\. Energy-spread unit fix

matRad_MCemittanceBaseData treats energySpectrum.sigma as a relative spread in percent: it copies the field straight into EnergySpread, and when it has to estimate the spread itself it computes spreadInMeV / MeanEnergy \* 100. Our scripts stored dEoverE_pct/100 \* E, i.e. MeV, which overstated the spread by a factor E/100 (≈1.5× at 150 MeV). The new script stores dEoverE_pct directly. The analytical pencil-beam dose does not read this field; MCsquare and TOPAS recalculation do.

6\. How to build and use the machine

1.  Run matRad_commissionDWAmachine_doubleGauss_interp.m (the three file paths at the top point to Alaina’s CSVs). The console should report: 58 energy blocks; 580 profiles with 0 failed fits; a width model close to W80 = 0.0396·R80^0.935; “58 native + 143 interpolated = 201 layers”; and a spacing ratio of at most ≈0.52 above 60 mm. These numbers come from a line-by-line Python mirror of the script.

2.  Check the two QA figures: core vs vacuum, halo width and weight; and range vs energy, layer spacing and lateral parameters for all layers.

3.  Save: obj = matRad_MCemittanceBaseData(machine); obj.saveMatradMachine('DWA_proton_doubleGauss_interp'); This fits beam optics per layer, so it takes about 3.5× longer than for 58 layers. The earlier “negative distance” air-correction warning may reappear.

4.  The file is written to userdata/machines/, where matRad finds it. For git, also copy it to matRad/basedata/ (userdata/\* is in .gitignore). matRad reads basedata first, so after any rebuild copy it again.

5.  Plan with pln.machine = 'DWA_proton_doubleGauss_interp';. The console must say “Using a double Gaussian pencil-beam kernel model”. The analysis scripts also set pln.propDoseCalc.lateralModel = 'double' as a safeguard, and a geometric lateral cutoff of 100 mm (about four halo sigmas).

**To check in matRad before trusting SOBPs.** How matRad’s spot-placement code picks energy layers for a target — in particular whether it thins layers to a fixed longitudinal spacing (look for longitudinalSpotSpacing in the particle stf generator). If it does, the default must not undo the finer grid.

7\. Single-beam validation in matRad (30 September)

The machine was built, saved and fired as single pencil beams into BOXPHANTOM with matRad_example5_protons_DWA_BOXPHANTOM_singleBeam.m (150 MeV) and matRad_DWA_singleBeam_MinMidMax.m (lowest, middle and highest native layer whose peak fits in the 240 mm water block). The first runs looked wrong in three ways; each had a clear cause.

7.1 Compare summed dose, not central-axis dose

**The central-axis curve of this beam has almost no Bragg peak, and that is correct.** The TOPAS Bragg curves are dose summed over the whole lateral plane at each depth. A single voxel on the axis receives that dose divided by the area the beam covers, and this beam widens a lot: at 181 MeV the core σ grows from 1.2 mm at the surface to 4.8 mm at the peak (Gottschalk’s rule of thumb, 2.3% of the range, gives 5.0 mm), a 16× larger area against a 4.8× higher summed dose. On the axis the entrance is therefore about three times hotter than the Bragg peak. With the old single-Gaussian machine (σ ≈ 5 mm everywhere) the two curves looked alike, which is why this never showed before.

The single-beam script now computes the laterally summed dose on a finer grid (1 mm along the beam, 2 mm across) and overlays the machine’s own input curve. On the 3 mm grid the summed curve was distorted: a 3 mm step cannot land on a 4.5 mm-wide peak, and summing 3 mm voxels over a 1.2 mm core inflates the entrance by 6–13%. The water surface is now read from the CT instead of being assumed 120 mm upstream of isocentre (it is 121.5 mm).



**Figure 5.** *150.16 MeV (an interpolated layer): matRad’s laterally summed dose on the 1 mm grid (blue) lies on the machine input curve (black). The grey dashed curve is the central-axis dose, normalised to its own maximum.*

7.2 matRad was silently dropping the halo

**Cause.** matRad does not choose single or double Gaussian from meta.dataType. Its default lateralModel = 'fast' uses a single Gaussian whenever machine.data has a sigma field (matRad_ParticlePencilBeamEngineAbstract.m, lines 130–133). The commissioning script kept sigma = in-water MCS “as a fallback”, so matRad computed the core only and discarded the halo, about half the fluence. Generic has no sigma field and was never affected.

**Evidence.** On a log scale the lateral profile followed the core exactly and dropped to zero where the halo should dominate (Figure 6); the depth view widened in steps, the signature of the dosimetric cutoff computed for a narrow single Gaussian. The depth-dose shape could not reveal it: the halo fraction is constant with depth, so dropping it scales every depth equally.



**Figure 6.** *Before the fix: lateral profiles at the Bragg peak (log scale). Circles: matRad. Black: the machine’s double Gaussian; blue dotted: its core; orange dashed: its halo. matRad followed the core only (halo ratio 0.00 and 0.04).*

**Fix.** The sigma field was removed from both saved copies of the machine (no rebuild needed, and no change to the physics) and from the commissioning script. The file had not left this computer, so it kept its name. The analysis scripts additionally force lateralModel = 'double'. After the fix the halo ratio (matRad ÷ model, 15–30 mm off axis) is 1.00 at all three energies.

7.3 Results at three energies

| **Layer** | **Energy** | **Bragg peak, matRad** | **Bragg peak, TOPAS** | **Halo ratio (matRad/model)** |
|-----------|------------|------------------------|-----------------------|-------------------------------|
| MIN       | 20.01 MeV  | 2.5 mm                 | 3.5 mm                | 1.00                          |
| MID       | 122.54 MeV | 108.5 mm               | 109.5 mm              | 1.00                          |
| MAX       | 181.01 MeV | 216.5 mm               | 216.5 mm              | 1.00                          |

The Bragg peak agrees within the 1 mm sampling at MID and MAX. MIN cannot be judged in this phantom: its whole 3.5 mm range fits in the first 3 mm CT voxel, and the first 1 mm dose voxel straddles the surface. (The first version of the summary CSV reported entrance dose 0 and a summed/input ratio of 0 for the same reason; the script now starts these comparisons at 1 mm depth.)

7.4 Single spot versus a small field

To show that the machine produces the familiar Bragg peak once spots overlap, the script adds a 5 × 5 field of spots 3 mm apart at the same energy, built by summing the single spot shifted by whole voxels (exact here because the phantom is uniform water across the beam; not valid in a patient). In the field, neighbouring spots replace the dose each one scatters off its own axis, and the Bragg peak becomes the hottest region again.

| **On-axis dose, Bragg peak ÷ entrance** | **single spot**        | **5 × 5 field, 3 mm spacing** |
|-----------------------------------------|------------------------|-------------------------------|
| MID 122.54 MeV                          | 1.16                   | 4.52                          |
| MAX 181.01 MeV                          | 0.22 (entrance hotter) | 2.59 (Bragg peak hotter)      |



**Figure 7.** *181.01 MeV, 5 × 5 spots 3 mm apart: dose in the plane through the beam axis, depth left to right. Solid white line: water surface; dashed: Bragg peak of the summed dose. The faint wider region towards the peak is the halo.*

This behaviour is specific to very narrow beams and matches what is known from proton minibeams. It matters for how DWA plans should be read: single-spot central-axis dose is not a meaningful check of this machine; summed dose, or dose in a field, is.

8\. Corrections to Part 2

- **Generic is not a cyclotron with a degrader.** Part 2 §5 explained Generic’s low-energy spot growth (σ ≈ 9.6 mm at 32 MeV) as the signature of degrading a fixed high-energy beam. Remo confirms Generic is a synthetic, TOPAS-generated machine with synchrotron-inspired optics. The lateral difference is between two beam-model design choices, not between two accelerator technologies, and the comparison is against a synthetic reference, not a clinical machine.

- **The halo is not 6% of the fluence.** Part 2 said the core carries about 94% of the fluence. That was the 1D share; the 2D share, which is what matters for dose, puts 38–52% in the halo above 50 MeV (Section 3).

- **The lateral and entrance-dose results in Part 2 used the single-Gaussian machine.** The FWHM ratios, the core-substitution table and the entrance-dose comparison all change with the double-Gaussian model and must be re-run. The range and distal-falloff results do not depend on the lateral model and stand.

9\. Open items

Needs new TOPAS output (Alaina)

1.  **Absolute depth-dose normalisation.** The curves are scaled to 100 at the peak. matRad’s Z should be dose per primary (MeV·cm²/g). Needed: per-primary integrated depth dose, or the AbsoluteDosimetry curve if one was made for RayStation.

2.  **LETd vs depth in water, per energy,** in the same geometry as the Bragg peaks. Every LET route needs machine.data(i).LET; the DWA machine has none.

3.  **Origin of the halo.** Does the TOPAS source at the exit window already contain it? (The nozzle file also lists the exit window material as Air, under a comment describing 100 µm Mylar.)

4.  **In-water lateral data,** to model the nuclear halo that grows with depth (Generic carries one; ours holds only the in-air halo).

5.  **IDD truncation.** Bragg curves were scored within r \< 60 mm, missing 0.3% (50 MeV) to 2.3% (219 MeV) of the dose.

For us / Remo

1.  Re-run Remo’s SOBP and the lateral profile at mid-SOBP with the new machine, with the double-Gaussian kernel active (check the console line).

2.  Re-run the DWA vs Generic sweeps (Part 2) with the new machine; the lateral and entrance results there all change.

3.  Choose the LET route (colleague’s implementation, Remo’s custom quantities, or pyRadPlan). The machine file works for all three once it has LET data.

4.  Resolve BAMStoIsoDist: 400 mm (snout to isocentre, used now) versus 550 mm (exit window to isocentre).

5.  Steps 2–4 of Remo’s validation plan: SOBP, TG119, patient.
