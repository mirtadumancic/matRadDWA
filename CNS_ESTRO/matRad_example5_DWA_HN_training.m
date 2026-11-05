%% matRad_example5_DWA_HN_training.m
% Training / rehearsal for the UHN CNS case on public data (matRad HEAD_AND_NECK).
% Example 5 (protons, constant RBE) with the DWA machine, run on DICOM files
% that have the same layout as the UHN patient (CT + RS + RP + RD):
%
%   1  photon IMRT plan on HEAD_AND_NECK, exported to DICOM  (stand-in for RP1/RD1)
%   2  import the DICOM folder with cns_importDicom          (the UHN step)
%   3  canonical structure names (CTV = GTV, PTV = PTV70)
%   4  DWA proton plan, dose-only, example-5 style  ->  GUI checkpoint
%   5  the same with a LET x dose objective         ->  GUI checkpoint
%   6  slice figures for slides (dose, McNamara, LETd, photon)
%   7  the full ESTRO pipeline on this case (both machines, all plans, photon reference)
%
% Run cell by cell (Ctrl+Enter) from the matRadDWA root on the let_opt branch.
% Everything here is public phantom data, so the Mac is fine.

%% 0  Setup ---------------------------------------------------------------
matRad_cfg = matRad_rc;
here = fileparts(which('cns_runCase'));
if isempty(here), addpath(fullfile(pwd, 'CNS_ESTRO')); here = fileparts(which('cns_runCase')); end

S = struct();
S.fast = true;                 % true: 5 mm spots, 4 mm grid, 150 iterations (minutes)
                               % false: study settings (3 mm, 500 iterations) for final pictures
                               % RAM: the DWA dose + LET dij for this phantom needs > 7 GB even
                               % in fast mode; at 3 mm plan on 32 GB or more
S.machine = 'DWA_proton_v2';   % or 'Generic' / 'IBA_cyclotron'
S.dicomDir = fullfile(matRad_cfg.primaryUserFolder, 'patients', 'HN_training_dicom');  % userdata: git-ignored
S.matFile  = cns_expandPath('~/Work/CNS_ESTRO_public/HN_training.mat');
S.caseFile = fullfile(here, 'cases', 'HN_matRad_training.json');
S.studyFile = fullfile(here, 'study_ESTRO2027.json');

caseCfg = cns_readJson(S.caseFile);
study = cns_readJson(S.studyFile);
if S.fast
    study.stf.bixelWidthMm.(S.machine) = 5;
    study.stf.longitudinalSpotSpacingMm = 4;
    study.doseCalc.doseGridMm = 4;
    study.optimisation.maxIter = 150;
end
matRad_cfg.defaults.propOpt.maxIter = study.optimisation.maxIter;
nFx = caseCfg.prescription.numOfFractions;
if ~exist(fileparts(S.matFile), 'dir'), mkdir(fileparts(S.matFile)); end

%% 1  Photon "clinical" plan on HEAD_AND_NECK and DICOM export -------------
% The phantom ships with objectives (PTV70 70 Gy, PTV63 63 Gy, parotids,
% skin). A 9-field IMRT plan stands in for the delivered photon plan (RP1/RD1).
load('HEAD_AND_NECK.mat');                       % ct, cst

plnPh = struct();
plnPh.radiationMode = 'photons';
plnPh.machine = 'Generic';
plnPh.bioModel = 'none';
plnPh.numOfFractions = nFx;
plnPh.propStf.gantryAngles = 0:40:320;
plnPh.propStf.couchAngles = zeros(size(plnPh.propStf.gantryAngles));
plnPh.propStf.bixelWidth = 5;
plnPh.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
gPh = 3 + 2 * S.fast;
plnPh.propDoseCalc.doseGrid.resolution = struct('x', gPh, 'y', gPh, 'z', gPh);
plnPh.propOpt.quantityOpt = 'physicalDose';
plnPh.propOpt.runDAO = false;
plnPh.propSeq.runSequencing = false;
plnPh.multScen = 'nomScen';

stfPh = matRad_generateStf(ct, cst, plnPh);
dijPh = matRad_calcDoseInfluence(ct, cst, stfPh, plnPh);
resPh = matRad_fluenceOptimization(dijPh, cst, plnPh);
clear dijPh

% Export CT + RTSTRUCT + RTDOSE. Only the physical dose is written.
% matRad blocks RTPLAN export on purpose (exported plans must never be used
% clinically), so the test folder has no RP file; the UHN folder will, and
% cns_importDicom falls back to CT + RS + RD if the clinical RTPLAN won't import.
% Note: pass resultGUI as a property; the exporter's 5th constructor
% argument is stored in the wrong field in this matRad version.
if exist(S.dicomDir, 'dir'), delete(fullfile(S.dicomDir, '*.dcm')); else, mkdir(S.dicomDir); end
dcmExp = matRad_DicomExporter(ct, cst, plnPh, stfPh);
dcmExp.resultGUI = struct('physicalDose', resPh.physicalDose);
dcmExp.dicomDir = S.dicomDir;
dcmExp.matRad_exportDicom();
fprintf('Exported %d DICOM files to %s\n', numel(dir(fullfile(S.dicomDir, '*.dcm'))), S.dicomDir);

%% 2  Import the DICOM folder (exactly what you will do with the UHN patient)
[ct, cst, refDose] = cns_importDicom(S.dicomDir, S.matFile);
% refDose.cube is the photon dose on the CT grid; refDose.toTotal turns it into
% the total dose (NaN here: matRad exports per-fraction dose and there is no
% RTPLAN, so the case file's 35 fractions are used). Check: median PTV70 dose
% should be close to 70 Gy. A clinical RD1 is usually a PLAN (total) dose.
clear plnPh stfPh resPh dcmExp          % free memory before the proton dij

%% 3  Canonical structures, PTV/brain/rings -------------------------------
[cst, structInfo] = cns_prepareStructures(ct, cst, caseCfg, study);
disp(structInfo.map);
% Brain: no contour in this phantom, so a HU-based proxy is built (flagged).

%% 4  DWA proton plan, dose-only (example 5 flow) --------------------------
pln = cns_buildPln(ct, cst, caseCfg, study, S.machine);   % HongPB, double Gaussian, LET on, RBE 1.1
disp(pln.propStf); disp(pln.propDoseCalc);

stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);          % dose + LET x dose influence matrices

planDose = struct('id', 'doseOnly', 'let', struct('mode', 'none'));
cstOpt = cns_buildCst(cst, caseCfg, study, planDose, pln);
resultGUI = matRad_fluenceOptimization(dij, cstOpt, pln);

% Total-dose cubes for the GUI (matRad shows dose per fraction by default)
resultGUI = cns_addTotalsForGUI(resultGUI, ct, cst, nFx, study);
resultGUI.photon_dose_total = refDose.cube * refDose.toTotal;
if isnan(refDose.toTotal), resultGUI.photon_dose_total = refDose.cube * nFx; end

% ---- GUI checkpoint 1 ----------------------------------------------------
% In the GUI: pick the cube in the "Plot" dropdown (RBExDose_total, McN_total,
% LET, photon_dose_total), move through slices, tick structures, open the DVH
% (Show DVH / Plan analysis). Then run the next cell.
matRadGUI;

%% 5  Same machine with a LET x dose objective -----------------------------
planLET = struct('id', 'LETxD_0p3', 'let', struct('mode', 'LETxDose', 'penalty', 0.3));
cstLET = cns_buildCst(cst, caseCfg, study, planLET, pln);
resLET = matRad_fluenceOptimization(dij, cstLET, pln);
resLET = cns_addTotalsForGUI(resLET, ct, cst, nFx, study);
% keep both plans in one resultGUI so the GUI can switch between them
for f = {'RBExDose_total', 'McN_total', 'LET'}
    resultGUI.([f{1} '_LETopt']) = resLET.(f{1});
end
resultGUI.RBExDose_diff_LETopt_minus_doseOnly = resLET.RBExDose_total - resultGUI.RBExDose_total;

% ---- GUI checkpoint 2: compare *_LETopt with the dose-only cubes ---------
matRadGUI;

%% 6  Slice figures for slides ---------------------------------------------
% Axial slice through the CTV centre; same window for all dose panels.
iC = find(strcmp(cst(:,2), 'CTV'), 1);
[~, ~, zz] = ind2sub(ct.cubeDim, cst{iC,4}{1});
z = round(median(zz));
rx = caseCfg.prescription.totalDoseGyRBE;
panels = {resultGUI.photon_dose_total, 'Photon (clinical stand-in) [Gy]', [0 1.1*rx]; ...
          resultGUI.RBExDose_total,    [S.machine ' dose-only [Gy(RBE)]'], [0 1.1*rx]; ...
          resultGUI.RBExDose_total_LETopt, [S.machine ' LET\timesD-optimised [Gy(RBE)]'], [0 1.1*rx]; ...
          resultGUI.McN_total,         'McNamara, dose-only [Gy(RBE)]', [0 1.1*rx]; ...
          resultGUI.LET,               'LETd, dose-only [keV/\mum]', [0 8]; ...
          resultGUI.LET_LETopt,        'LETd, LET\timesD-optimised [keV/\mum]', [0 8]};
plotNames = {'CTV', 'PTV', 'Brainstem', 'Chiasm', 'OpticNerve_L', 'OpticNerve_R', 'SpinalCord'};
% (cns_plotSlice blends the wash into the CT itself: it prints the same in
% MATLAB and Octave. The GUI and matRad_plotSlice are fine for browsing.)
fig = figure('Color', 'w', 'Position', [50 50 1500 900]);
for k = 1:size(panels, 1)
    ax = subplot(2, 3, k);
    D = panels{k,1};
    if ~isempty(strfind(panels{k,2}, 'LETd')), D(resultGUI.RBExDose_total < 1) = 0; end   % LET where dose >= 1 Gy(RBE)
    thr = 0.1; if ~isempty(strfind(panels{k,2}, 'LETd')), thr = 0.02; end
    cns_plotSlice(ax, ct, cst, D, z, panels{k,3}, jet(64), plotNames, '', k == 1, thr);
    title(ax, panels{k,2}, 'Interpreter', 'tex');
end
figDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias, 'training_figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
print(fig, fullfile(figDir, ['slices_' S.machine '.png']), '-dpng', '-r200');
fprintf('Saved %s\n', fullfile(figDir, ['slices_' S.machine '.png']));

%% 7  The full ESTRO pipeline on this case ----------------------------------
% Both machines, every plan of study_ESTRO2027.json, normalisation, all metrics,
% plus the photon reference ("photon / clinical" rows). This repeats the dij
% per machine, so it takes longer than sections 4-5.
studyRun = S.studyFile;
if S.fast                                 % write the fast settings to a temporary study file
    studyRun = [tempname '.json'];
    cns_writeJson(studyRun, study);
end
out = cns_runCase(S.caseFile, studyRun);
cns_cohortSummary({out.caseDir}, studyRun, fullfile(fileparts(out.caseDir), 'cohort_training'));
type(fullfile(fileparts(out.caseDir), 'cohort_training', 'summary_abstract.md'));
