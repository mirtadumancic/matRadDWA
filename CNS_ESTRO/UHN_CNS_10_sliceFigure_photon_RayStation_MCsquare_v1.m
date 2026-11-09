%% UHN_CNS_10_sliceFigure_photon_RayStation_MCsquare_v1.m
% Three-panel slice figure for one patient, dose in % of the prescription:
%   Photon (delivered) | DWA, RayStation MC | DWA, MCsquare recalculated
% Nothing is recalculated: the doses are read from the files written by UHN_CNS_09 v2
% (rs_doses_<alias>_v2_PATIENT_DATA_keep_on_PC.mat) and from prepared.mat. Takes about a minute.
%   clear all; run('CNS_ESTRO/run_UHN_CNS_10_CNS02_v1.m')
% Output: <case>_DWAv4rob/<S.doseFolder>/fig_slices_<plane>_<alias>_v2.(pdf|png) and _slices.mat
% PATIENT IMAGES: the output stays on hospital storage; check consent / congress rules before use.
S = struct();
S.caseFile   = 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/CNS_02.json';
S.doseFolder = 'results_RayStation_v2_LVH';   % folder with the rs_doses_*.mat (or 'results_RayStation_v2')
S.planes     = {'axial'};                     % any of 'axial', 'sagittal', 'coronal'
S.slice      = [];                            % [] = slice through the centre of the CTV; or a slice index
S.crop       = {'Body'};                      % whole patient outline on the slice; {'Brain','Eye_L','Eye_R'} = brain only
S.cropMarginMm = 8;
S.background = 'dark';                        % 'dark' (planning-system look) | 'light'
S.widthCm    = 17.5;
if exist('S_override', 'var') && isstruct(S_override)
    fo = fieldnames(S_override); for k = 1:numel(fo), S.(fo{k}) = S_override.(fo{k}); end
end

here = fileparts(mfilename('fullpath')); addpath(here);
caseCfg0 = cns_readJson(S.caseFile); cns_checkPrivacy(caseCfg0, S.caseFile);
alias = caseCfg0.alias; rx = caseCfg0.prescription.totalDoseGyRBE; nFx = caseCfg0.prescription.numOfFractions;
caseDir = fullfile(cns_expandPath(caseCfg0.workDir), [alias '_DWAv4rob']);
outDir = fullfile(caseDir, S.doseFolder);
fDose = fullfile(outDir, sprintf('rs_doses_%s_v2_PATIENT_DATA_keep_on_PC.mat', alias));
if ~exist(fDose, 'file'), error('Not found: %s (run UHN_CNS_09 v2 first, or set S.doseFolder).', fDose); end

P = load(fullfile(caseDir, 'prepared.mat'), 'ct', 'cst', 'refDose');
R = load(fDose, 'Drs', 'Drsm_mcfit');
if isempty(R.Drsm_mcfit), error('%s has no MCsquare recalculation of the RayStation plan (Drsm_mcfit is empty).', fDose); end
if ~isequal(size(R.Drs), P.ct.cubeDim) || ~isequal(size(R.Drsm_mcfit), P.ct.cubeDim)
    error('The dose cubes are not on the CT grid of prepared.mat.');
end

pan = struct('cube', {}, 'type', {}, 'title', {}, 'row', {}, 'col', {}, 'doseForMask', {});
if ~isempty(P.refDose)
    t = P.refDose.toTotal; if isnan(t), t = nFx; end
    pan(end+1) = struct('cube', single(double(P.refDose.cube) * t), 'type', 'dose', 'title', 'Photon (delivered)', 'row', 1, 'col', 1, 'doseForMask', []);
else
    fprintf(2, '[UHN_CNS_10] no photon dose in prepared.mat: two panels only.\n');
end
pan(end+1) = struct('cube', single(R.Drs), 'type', 'dose', 'title', 'DWA, RayStation MC', 'row', 1, 'col', numel(pan) + 1, 'doseForMask', []);
pan(end+1) = struct('cube', single(R.Drsm_mcfit), 'type', 'dose', 'title', 'DWA, MCsquare recalculated', 'row', 1, 'col', numel(pan) + 1, 'doseForMask', []);

contours = struct( ...
    'name',  {'GTV', 'CTV', 'PTV', 'OpticNerve_L', 'OpticNerve_R', 'Chiasm', 'Brainstem', 'Eye_L', 'Eye_R', 'Body'}, ...
    'color', {[0.15 0.55 0.15], [0.85 0.10 0.85], [0.30 0.60 1.00], [1.00 0.85 0.10], [1.00 0.85 0.10], ...
              [0.20 0.95 0.95], [0.95 0.60 0.20], [0.60 0.80 1.00], [0.60 0.80 1.00], [0.20 0.75 0.20]}, ...
    'label', {'', 'CTV', '', 'ON_L', 'ON_R', 'Chiasm', 'Brainstem', '', '', ''}, ...
    'labelOffsetMm', {[0 0], [-22 -18], [0 0], [25 -10], [-25 -10], [18 -20], [20 -15], [0 0], [0 0], [0 0]});
vis0 = get(0, 'DefaultFigureVisible'); set(0, 'DefaultFigureVisible', 'off');
for ip = 1:numel(S.planes)
    oS = struct('plane', S.planes{ip}, 'centreOn', 'CTV', 'rxGy', rx, 'washMinPct', 50, 'isoPct', [50 95 107], ...
        'contours', contours, 'crop', {S.crop}, 'cropMarginMm', S.cropMarginMm, 'background', S.background, 'widthCm', S.widthCm);
    if ~isempty(S.slice), oS.slice = S.slice; end
    base = fullfile(outDir, sprintf('fig_slices_%s_%s_v2', S.planes{ip}, alias));
    try
        cns_figSliceRS_v2(P.ct, P.cst, pan, oS, base);
        fprintf('[UHN_CNS_10] %s.png\n', base);
    catch err
        fprintf(2, '[UHN_CNS_10] slice figure %s failed: %s\n', S.planes{ip}, getReport(err));
    end
end
set(0, 'DefaultFigureVisible', vis0);
