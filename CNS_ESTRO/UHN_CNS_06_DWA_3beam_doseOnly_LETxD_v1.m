%% UHN_CNS_06_DWA_3beam_doseOnly_LETxD_v1.m
% DWA only, custom 3-beam arrangement, three plans compared:
%   Photon (delivered)  |  DWA dose-only  |  DWA LET x D (fixed weight, default 0.1)
%
% One run does everything:
%   1  dose + LET influence matrix for DWA with these beams, dose-only plan (cns_runCase)
%   2  LET x D plan at the fixed weight (warm start from dose-only, same normalisation)
%   3  metrics table (photon / dose-only / LET x D) + constraint check
%   4  Figure: DVH + LVH + dose-LET bars (journal style)
%   5  TPS-style slices (axial, sagittal, coronal): dose wash + isodoses, LETd, Delta dose
%
% Uses its own folder <workDir>/<alias>_DWA3b/ so the cached 2-beam (G) data are not mixed in.
% RUN ON THE UHN PC ONLY (patient data). Time: DWA dij ~10-30 min (3 beams), 2 optimisations
% ~5-15 min, figures ~2 min.

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';
S.studyFile = '';                                 % '' = CNS_ESTRO/study_UHN_CNS_v1.json
S.machine   = 'DWA_proton_v2';
S.gantryAngles = [270 100 75];
S.couchAngles  = [0   0   300];
S.setLabel  = '3beam_G270_G100_G75C300';
S.letWeight = 0.1;                                % LET x D penalty weight (fixed, no search)
S.planes    = {'axial', 'sagittal', 'coronal'};
S.crop      = {'Brain', 'Eye_L', 'Eye_R'};       % slices: brain + eyes + margin
S.cropMarginMm = 15;
S.background = 'dark';

matRadRoot = getenv('MATRAD_ROOT');
if isempty(matRadRoot)
    here = mfilename('fullpath');
    if isempty(here), here = which('cns_runCase'); end
    if isempty(here), here = fullfile(pwd, 'x'); end
    matRadRoot = fileparts(fileparts(here));
end
if ~exist(fullfile(matRadRoot, 'matRad_rc.m'), 'file')
    error(['matRad root not found (%s). Set it once per session with ' ...
        'setenv(''MATRAD_ROOT'', ''<path to matRadDWA>'').'], matRadRoot);
end
cd(matRadRoot); matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;
addpath(fullfile(matRadRoot, 'CNS_ESTRO'));
if isempty(S.studyFile), S.studyFile = fullfile(matRadRoot, 'CNS_ESTRO', 'study_UHN_CNS_v1.json'); end
if ~exist(S.caseFile, 'file')
    [f, p] = uigetfile('*.json', 'Case file of the patient (made with cns_newCaseFile_v1)');
    if isequal(f, 0), error('Case file not found: %s', S.caseFile); end
    S.caseFile = fullfile(p, f);
end
caseCfg0 = cns_readJson(S.caseFile);
study    = cns_readJson(S.studyFile);
cns_checkPrivacy(caseCfg0, S.caseFile);
mc = S.machine; tag = strrep(mc, '-', '_');
rx  = caseCfg0.prescription.totalDoseGyRBE;
nFx = caseCfg0.prescription.numOfFractions;

% case + study variants for this run (written next to the results, outside the repo)
caseCfg = caseCfg0;
caseCfg.alias = [caseCfg0.alias '_DWA3b'];
caseCfg.beams = struct('setId', S.setLabel, 'gantryAngles', S.gantryAngles, 'couchAngles', S.couchAngles);
caseDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias);
if ~exist(caseDir, 'dir'), mkdir(caseDir); end
study.machines = {mc};
study.plans = cns_toCell(study.plans);
study.plans = study.plans(cellfun(@(p) strcmp(cns_getfield(p, 'let.mode', 'none'), 'none'), study.plans));
caseFile = fullfile(caseDir, [caseCfg.alias '.json']);
studyFile = fullfile(caseDir, 'study_DWA3b_used.json');
cns_writeJson(caseFile, caseCfg);
cns_writeJson(studyFile, study);
matRad_cfg.defaults.propOpt.maxIter = cns_getfield(study, 'optimisation.maxIter', 500);
outDir = fullfile(caseDir, 'results_v1');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('[UHN_CNS_06] %s, beams gantry [%s] couch [%s], LET x D weight %g\n', mc, ...
    num2str(S.gantryAngles), num2str(S.couchAngles), S.letWeight);

%% 1  Dose-only plan (+ photon evaluation) ---------------------------------------
out = cns_runCase(caseFile, studyFile);
rows = out.rows;
ok = strcmp(rows(:,2), mc) & strcmp(rows(:,3), 'doseOnly');
if ~any(ok), error('The DWA dose-only plan failed: see the log above.'); end

%% 2  LET x D plan at the fixed weight ---------------------------------------------
P = load(fullfile(caseDir, 'prepared.mat'));
refDose = []; if isfield(P, 'refDose'), refDose = P.refDose; end
D = load(fullfile(caseDir, ['dij_' mc '.mat']), 'pln', 'dij');
pln = D.pln;
tmp = cns_buildPln(P.ct, P.cst, caseCfg, study, mc);
pln.multScen = tmp.multScen;
W = load(fullfile(caseDir, 'plans', [mc '__doseOnly_w.mat']), 'w');
cstD = matRad_resizeCstToGrid(P.cst, D.dij.ctGrid.x, D.dij.ctGrid.y, D.dij.ctGrid.z, ...
    D.dij.doseGrid.x, D.dij.doseGrid.y, D.dij.doseGrid.z);
wDoseN = cns_normaliseWeights(W.w, D.dij, cstD, caseCfg, study, refDose, P.cst);

% LET x D reference: 5th percentile of LET x dose in the CTV of the dose-only plan (as UHN_CNS_03)
C = cns_getfield(study, 'letSearch', struct());
covS = cns_getfield(cns_getfield(C, 'coverage', struct()), 'structure', 'CTV');
iC = find(strcmp(cstD(:,2), covS), 1);
LxD = full(D.dij.mLETDose{1} * wDoseN) * nFx;
x = sort(LxD(cstD{iC,4}{1}));
ref = cns_getfield(C, 'refFactor', 1) * x(max(1, ceil(0.05 * numel(x))));
clear LxD x
letS = cns_toCell(cns_getfield(C, 'letStructures', cns_getfield(study, 'letStructures', {})));
fprintf('[UHN_CNS_06] LET x D reference %.1f keV/um*Gy, structures: %s\n', ref, strjoin(letS, ', '));

Pdef = struct('id', 'LETxD_best', 'let', struct('mode', 'LETxDose', 'penalty', S.letWeight, ...
    'ref', ref, 'structures', {letS}));
cstRun = cns_buildCst(P.cst, caseCfg, study, Pdef, pln);
[res, ~, optUsed] = cns_optimize(D.dij, cstRun, pln, study, W.w);
wLET = cns_normaliseWeights(res.wUnsequenced, D.dij, cstD, caseCfg, study, refDose, P.cst);
wDose = wDoseN; %#ok<NASGU>
save(fullfile(caseDir, 'plans', sprintf('%s__LETxD_w%g_v1.mat', mc, S.letWeight)), 'wLET', 'wDose', 'ref', 'optUsed');

E0 = cns_evaluatePlan(D.dij, P.cst, wDoseN, caseCfg, study);
E1 = cns_evaluatePlan(D.dij, P.cst, wLET, caseCfg, study);
rows = rows(~(strcmp(rows(:,2), mc) & strcmp(rows(:,3), 'LETxD_best')), :);
rows = [rows; repmat({caseCfg.alias, mc, 'LETxD_best'}, size(E1.rows, 1), 1), E1.rows];
cns_writeCsvTable(fullfile(outDir, 'metrics_DWA3b_v1.csv'), ...
    {'alias', 'machine', 'plan', 'structure', 'quantity', 'metric', 'value', 'unit'}, rows);

% cubes for the slices (total course, RBE 1.1)
R0 = matRad_calcCubes(wDoseN, D.dij);  R1 = matRad_calcCubes(wLET, D.dij);
Dp0 = 1.1 * nFx * R0.physicalDose;  L0 = R0.LET;  L0(~isfinite(L0)) = 0;
Dp1 = 1.1 * nFx * R1.physicalDose;  L1 = R1.LET;  L1(~isfinite(L1)) = 0;
clear D R0 R1

%% 3  Table ----------------------------------------------------------------------
wl = sprintf('LETxD w=%g', S.letWeight);
cols = struct('machine', {'photon', mc, mc}, 'plan', {'clinical', 'doseOnly', 'LETxD_best'}, ...
    'label', {'Photon', 'DWA dose-only', ['DWA ' wl]});
if isempty(refDose), cols = cols(2:3); end
cns_compareTable(rows, study.summary.keyMetrics, cols, fullfile(outDir, 'table_DWA3b_v1'), ...
    cns_getfield(study, 'constraintsCheck.list', {}));
type(fullfile(outDir, 'table_DWA3b_v1.md'));

%% 4  DVH / LVH / bars figure -------------------------------------------------------
set(0, 'DefaultFigureVisible', 'off');
dvh = []; lvh = [];
for E = {E0, E1; 'doseOnly', 'LETxD_best'}
    e = E{1};
    for i = 1:numel(e.dvh), e.dvh(i).machine = mc; e.dvh(i).plan = E{2}; end
    for i = 1:numel(e.lvh), e.lvh(i).machine = mc; e.lvh(i).plan = E{2}; end
    dvh = [dvh, orderfields(e.dvh)]; lvh = [lvh, orderfields(e.lvh)]; %#ok<AGROW>
end
try
    ph = out.dvh(strcmp({out.dvh.machine}, 'photon'));
    if ~isempty(ph)
        ph = rmfield(ph, setdiff(fieldnames(ph), fieldnames(dvh)));
        dvh = [dvh, orderfields(ph, dvh(1))];
    end
catch err
    fprintf(2, '[UHN_CNS_06] photon DVH not added (%s)\n', err.message);
end
barM = struct( ...
    'structure', {'OpticNerve_L', 'Chiasm', 'OpticNerve_L', 'Chiasm', 'BrainRing_0_15'}, ...
    'quantity',  {'DxLETd', 'DxLETd', 'McN', 'McN', 'hotSpot'}, ...
    'metric',    {'D2', 'D2', 'D0_03cc', 'D0_03cc', 'V_D20_L4_cc'}, ...
    'label',     {'Optic nerve L', 'Chiasm', 'Optic nerve L', 'Chiasm', 'Brain ring 0-15 mm'}, ...
    'sublabel',  {'D\timesLET_d D_{2%} (Gy\cdotkeV/\mum)', 'D\timesLET_d D_{2%} (Gy\cdotkeV/\mum)', ...
                  'McNamara D_{0.03cm^3} (Gy(RBE))', 'McNamara D_{0.03cm^3} (Gy(RBE))', ...
                  'V(D\geq20 Gy, LET_d\geq4) (cm^3)'});
o = struct('letdMinDoseGy', study.evaluation.letdMinDoseGy, 'lvh', lvh, ...
    'lvhStructs', {{'OpticNerve_L', 'Chiasm', 'Brainstem', 'BrainRing_0_15'}});
try
    cns_figDVHandBars(dvh, rows, {mc}, {'CTV', 'OpticNerve_L', 'Chiasm', 'Brainstem', 'Brain_healthy'}, ...
        barM, rx, fullfile(outDir, 'fig_DVH_doseLET_DWA3b_v1'), o);
catch err
    fprintf(2, '[UHN_CNS_06] DVH figure failed (%s); the table and slices are still made\n', err.message);
end

%% 5  TPS-style slices -------------------------------------------------------------
pan = struct('cube', {}, 'type', {}, 'title', {}, 'row', {}, 'col', {}, 'doseForMask', {});
if ~isempty(refDose)
    t = refDose.toTotal; if isnan(t), t = nFx; end
    pan(end+1) = struct('cube', double(refDose.cube) * t, 'type', 'dose', 'title', 'Photon (delivered)', 'row', 1, 'col', 1, 'doseForMask', []);
end
wLab = sprintf('LET\\timesD (w = %g)', S.letWeight);
pan(end+1) = struct('cube', Dp0, 'type', 'dose', 'title', 'DWA: dose-only', 'row', 1, 'col', 2, 'doseForMask', []);
pan(end+1) = struct('cube', Dp1, 'type', 'dose', 'title', ['DWA: ' wLab], 'row', 1, 'col', 3, 'doseForMask', []);
pan(end+1) = struct('cube', L0, 'type', 'let', 'title', 'LET_d, dose-only', 'row', 2, 'col', 2, 'doseForMask', Dp0);
pan(end+1) = struct('cube', L1, 'type', 'let', 'title', ['LET_d, ' wLab], 'row', 2, 'col', 3, 'doseForMask', Dp1);
pan(end+1) = struct('cube', Dp1 - Dp0, 'type', 'diff', 'title', '\DeltaDose (LET\timesD - dose-only)', ...
    'row', 3, 'col', 3, 'doseForMask', max(Dp0, Dp1));
contours = struct( ...
    'name',  {'GTV', 'CTV', 'PTV', 'OpticNerve_L', 'OpticNerve_R', 'Chiasm', 'Brainstem', 'Eye_L', 'Eye_R', 'Body'}, ...
    'color', {[0.15 0.55 0.15], [0.85 0.10 0.85], [0.30 0.60 1.00], [1.00 0.85 0.10], [1.00 0.85 0.10], ...
              [0.20 0.95 0.95], [0.95 0.60 0.20], [0.60 0.80 1.00], [0.60 0.80 1.00], [0.20 0.75 0.20]}, ...
    'label', {'', 'CTV', '', 'ON_L', '', 'Chiasm', 'Brainstem', '', '', ''}, ...
    'labelOffsetMm', {[0 0], [-22 -18], [0 0], [25 -10], [0 0], [18 -20], [20 -15], [0 0], [0 0], [0 0]});
for i = 1:numel(S.planes)
    oS = struct('plane', S.planes{i}, 'centreOn', 'CTV', 'rxGy', rx, 'washMinPct', 50, ...
        'isoPct', [50 95 107], 'letWindow', [1 8], 'letMinDoseGy', 0.1 * rx, 'diffWindow', [-5 5], ...
        'diffUnit', '\DeltaDose (Gy(RBE))', 'contours', contours, 'crop', {S.crop}, ...
        'cropMarginMm', S.cropMarginMm, 'background', S.background);
    base = fullfile(outDir, sprintf('fig_slices_%s_DWA3b_v1', S.planes{i}));
    try
        cns_figSliceRS(P.ct, P.cst, pan, oS, base);
        fprintf('[UHN_CNS_06] %s.png\n', base);
    catch err
        fprintf(2, '[UHN_CNS_06] slice figure %s failed: %s\n', S.planes{i}, err.message);
    end
end
set(0, 'DefaultFigureVisible', 'on');
fprintf('[UHN_CNS_06] done (optimiser %s). Results: %s\n', optUsed, outDir);
