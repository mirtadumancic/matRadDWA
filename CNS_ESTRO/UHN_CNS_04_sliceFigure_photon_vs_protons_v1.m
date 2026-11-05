%% UHN_CNS_04_sliceFigure_photon_vs_protons_v1.m
% CODE 4: the "money plot". One representative slice, TPS (RayStation-like) style:
%
%   row 1  dose wash + isodoses:  Photon (delivered) | Proton dose-only | Proton LET x dose
%   row 2  LETd:                  (photon: none)      | LETd dose-only   | LETd LET x dose
%   row 3  difference:                                                       Delta dose (LET x D - dose-only)
%
% Needs UHN_CNS_01, 02 and 03 (it reads prepared.mat, dij_<machine>.mat and
% LETxD_search_v1/plan_<machine>_LETxDbest_v1.mat). Run anywhere the patient data may be.
% Output: <workDir>/<alias>/figures_money_v1/  (PDF + 600 dpi PNG + slice .mat)
% PATIENT IMAGE: the default crop to the brain hides the face profile on a sagittal
% slice. Check consent and journal rules before any patient image leaves hospital storage.

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';
S.studyFile = '';                     % '' = CNS_ESTRO/study_UHN_CNS_v1.json
S.machine   = 'DWA_proton_v2';        % proton machine to show
S.planes    = {'sagittal', 'axial'};  % one figure per plane
S.slice     = struct('sagittal', [], 'axial', [], 'coronal', []);   % [] = through the CTV centre
S.background = 'dark';                % 'dark' (TPS look) or 'light' (print)
S.rows      = 2;                      % 1 = dose only, 2 = + LETd, 3 = + Delta dose (LET x D - dose-only)
S.crop      = 'Brain';                % 'Brain' hides the face on sagittal slices; 'Body' shows all
S.washMinPct = 50;                    % lowest dose shown in colour (% of Rx)
% structures drawn (RGB colours as in the TPS screenshot; label '' = no arrow label)
S.contours = struct( ...
    'name',  {'GTV', 'CTV', 'PTV', 'OpticNerve_L', 'Chiasm', 'Brainstem', 'Eye_L', 'Body'}, ...
    'color', {[0.15 0.55 0.15], [0.85 0.10 0.85], [0.30 0.60 1.00], [1.00 0.85 0.10], ...
              [0.20 0.95 0.95], [0.95 0.60 0.20], [0.60 0.80 1.00], [0.20 0.75 0.20]}, ...
    'label', {'GTV', 'CTV', '', 'ON_L', 'Chiasm', 'Brainstem', '', ''}, ...
    'labelOffsetMm', {[20 -18], [-22 -18], [0 0], [25 -10], [18 -20], [20 -15], [0 0], [0 0]});

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
caseCfg = cns_readJson(S.caseFile);
study   = cns_readJson(S.studyFile);
cns_checkPrivacy(caseCfg, S.caseFile);
caseDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias);
outDir = fullfile(caseDir, 'figures_money_v1');
if ~exist(outDir, 'dir'), mkdir(outDir); end
rx = caseCfg.prescription.totalDoseGyRBE;
nFx = caseCfg.prescription.numOfFractions;

%% 1  Dose and LET cubes on the CT grid -----------------------------------------
P = load(fullfile(caseDir, 'prepared.mat'));
ct = P.ct; cst = P.cst;
t = P.refDose.toTotal; if isnan(t), t = nFx; end
Dphoton = double(P.refDose.cube) * t;

D = load(fullfile(caseDir, ['dij_' S.machine '.mat']), 'dij');
tag = strrep(S.machine, '-', '_');
best = fullfile(caseDir, 'LETxD_search_v1', ['plan_' tag '_LETxDbest_v1.mat']);
if exist(best, 'file')
    W = load(best);
    wDose = W.wDoseNorm; wLET = W.wBest; wLab = sprintf('LET\\timesD (w = %.2g)', W.bestW);
else
    warning('No LET x dose result for %s yet (run UHN_CNS_03): showing dose-only only.', S.machine);
    Wd = load(fullfile(caseDir, 'plans', [S.machine '__doseOnly_w.mat']));
    cstD = matRad_resizeCstToGrid(cst, D.dij.ctGrid.x, D.dij.ctGrid.y, D.dij.ctGrid.z, ...
        D.dij.doseGrid.x, D.dij.doseGrid.y, D.dij.doseGrid.z);
    wDose = cns_normaliseWeights(Wd.w, D.dij, cstD, caseCfg, study, P.refDose, cst);
    wLET = []; S.rows = min(S.rows, 2);
end
R0 = matRad_calcCubes(wDose, D.dij);
Dp0 = 1.1 * nFx * R0.physicalDose;  L0 = R0.LET;  L0(~isfinite(L0)) = 0;
if ~isempty(wLET)
    R1 = matRad_calcCubes(wLET, D.dij);
    Dp1 = 1.1 * nFx * R1.physicalDose;  L1 = R1.LET;  L1(~isfinite(L1)) = 0;
end
clear D
mName = strrep(S.machine, '_', ' ');

%% 2  Panels -------------------------------------------------------------------
pan = struct('cube', {}, 'type', {}, 'title', {}, 'row', {}, 'col', {}, 'doseForMask', {});
pan(end+1) = struct('cube', Dphoton, 'type', 'dose', 'title', 'Photon (delivered)', 'row', 1, 'col', 1, 'doseForMask', []);
pan(end+1) = struct('cube', Dp0, 'type', 'dose', 'title', [mName ': dose-only'], 'row', 1, 'col', 2, 'doseForMask', []);
if ~isempty(wLET)
    pan(end+1) = struct('cube', Dp1, 'type', 'dose', 'title', [mName ': ' wLab], 'row', 1, 'col', 3, 'doseForMask', []);
end
if S.rows >= 2
    pan(end+1) = struct('cube', L0, 'type', 'let', 'title', 'LET_d, dose-only', 'row', 2, 'col', 2, 'doseForMask', Dp0);
    if ~isempty(wLET)
        pan(end+1) = struct('cube', L1, 'type', 'let', 'title', ['LET_d, ' wLab], 'row', 2, 'col', 3, 'doseForMask', Dp1);
    end
end
if S.rows >= 3 && ~isempty(wLET)
    pan(end+1) = struct('cube', Dp1 - Dp0, 'type', 'diff', 'title', '\DeltaDose (LET\timesD - dose-only)', ...
        'row', 3, 'col', 3, 'doseForMask', max(Dp0, Dp1));
end

%% 3  Figures, one per plane -------------------------------------------------
set(0, 'DefaultFigureVisible', 'off');
for i = 1:numel(S.planes)
    o = struct('plane', S.planes{i}, 'slice', S.slice.(S.planes{i}), 'rxGy', rx, ...
        'washMinPct', S.washMinPct, 'isoPct', [50 95 107], 'letWindow', [1 8], ...
        'letMinDoseGy', 0.1 * rx, 'diffWindow', [-5 5], 'diffUnit', '\DeltaDose (Gy(RBE))', ...
        'contours', S.contours, 'crop', S.crop, 'background', S.background);
    base = fullfile(outDir, sprintf('fig_money_%s_%s_v1', S.planes{i}, tag));
    cns_figSliceRS(ct, cst, pan, o, base);
    fprintf('[UHN_CNS_04] %s.png\n', base);
end
set(0, 'DefaultFigureVisible', 'on');
