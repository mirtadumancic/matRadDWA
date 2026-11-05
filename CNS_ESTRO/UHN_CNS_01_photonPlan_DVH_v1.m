%% UHN_CNS_01_photonPlan_DVH_v1.m
% CODE 1 of 3: the delivered photon plan in matRad. Imports, visualises and evaluates it.
%
%   1  import the patient folder once (CT + RS1 + RP1 + RD1; MR and REG are ignored)
%   2  structure mapping (CTVp_5400 -> CTV, PTVp_5400 -> PTV, ...), body outline from the CT
%   3  DVH metrics of the photon dose (RBE = 1) with the same metric set as the
%      proton plans: targets D98/D95/D50/D2/Dmean/V95/V107, OARs Dmean/D2/D0.03cc/Dmax/V20/V40
%   4  constraint check (study.constraintsCheck: replace with the UHN protocol)
%   5  journal-style figures (vector PDF + 300 dpi PNG + csv)
%   6  matRadGUI with the photon dose
%
% RUN ON THE UHN PC ONLY (patient data). From the matRad root:
%   setenv('MATRAD_ROOT', 'X:\matRadDWA_LET');      % once per session
%   run(fullfile(getenv('MATRAD_ROOT'), 'CNS_ESTRO', 'UHN_CNS_01_photonPlan_DVH_v1.m'))
% Time: import 1-3 min, evaluation < 1 min.
% Output: <workDir>/<alias>/photon_v1/  (metrics, DVH csv, tables, figures)

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';   % your copy of cases/UHN_CNS_ependymoma_TEMPLATE_v1.json
S.studyFile = '';            % '' = CNS_ESTRO/study_UHN_CNS_v1.json
S.openGUI   = true;
S.reimport  = false;         % true: import the DICOM folder again even if the .mat exists

matRadRoot = getenv('MATRAD_ROOT');
if isempty(matRadRoot)       % default: the matRadDWA folder that contains CNS_ESTRO/
    here = mfilename('fullpath');
    if isempty(here), here = which('cns_runCase'); end      % running cell by cell
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
    if isequal(f, 0)
        error(['Case file not found: %s\nMake one with caseFile = cns_newCaseFile_v1() ' ...
            'and set S.caseFile to it.'], S.caseFile);
    end
    S.caseFile = fullfile(p, f);
end
caseCfg = cns_readJson(S.caseFile);
study   = cns_readJson(S.studyFile);
dataClass = cns_checkPrivacy(caseCfg, S.caseFile);       % refuses AI-reachable folders for patient data
outDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias, 'photon_v1');
figDir = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
fprintf('[UHN_CNS_01] case %s (%s) -> %s\n', caseCfg.alias, dataClass, outDir);

%% 1  Import (once) -----------------------------------------------------------
matFile = cns_expandPath(caseCfg.source.path);
if ~exist(matFile, 'file') || S.reimport
    dcmDir = cns_expandPath(cns_getfield(caseCfg, 'source.dicomDir', ''));
    if isempty(dcmDir), error('Set source.dicomDir in the case file.'); end
    [ct, cst, refDose] = cns_importDicom(dcmDir, matFile);
else
    L = load(matFile); ct = L.ct; cst = L.cst; refDose = L.refDose; clear L
end
if isempty(refDose), error('No RTDOSE was imported: check that RD1 is in source.dicomDir.'); end
fprintf('[UHN_CNS_01] RTDOSE: %s, %s, factor to total %g; RTPLAN fractions: %g; case file: %d\n', ...
    refDose.doseType, refDose.summationType, refDose.toTotal, refDose.numOfFractions, ...
    caseCfg.prescription.numOfFractions);

%% 2  Structures ---------------------------------------------------------------
[cst, structInfo] = cns_prepareStructures(ct, cst, caseCfg, study);
fid = fopen(fullfile(outDir, 'structures_v1.csv'), 'w');
fprintf(fid, 'name_in_matRad,name_in_RTSTRUCT,volume_cc\n');
fn = fieldnames(structInfo.map);
for i = 1:size(cst, 1)
    orig = '';
    for j = 1:numel(fn), if strcmp(fn{j}, cst{i,2}), orig = structInfo.map.(fn{j}); end, end
    fprintf(fid, '%s,%s,%.3f\n', cst{i,2}, strtrim(orig), numel(cst{i,4}{1}) * ...
        ct.resolution.x * ct.resolution.y * ct.resolution.z / 1000);
end
fclose(fid);

%% 3  Photon plan metrics (RBE = 1, so RBE11 and McN rows equal the physical dose; LET = 0)
E = cns_evaluateReferenceDose(ct, cst, refDose, caseCfg, study);
rows = [repmat({caseCfg.alias, 'photon', 'clinical'}, size(E.rows, 1), 1), E.rows];
fid = fopen(fullfile(outDir, 'metrics_photon_v1.csv'), 'w');
fprintf(fid, 'alias,machine,plan,structure,quantity,metric,value,unit\n');
for i = 1:size(rows, 1)
    if any(strcmp(rows{i,5}, {'McN', 'LETd', 'DxLETd', 'hotSpot'})), continue; end   % meaningless for photons
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%.6g,%s\n', rows{i,1:6}, rows{i,7}, rows{i,8});
end
fclose(fid);
dvh = E.dvh; for i = 1:numel(dvh), dvh(i).machine = 'photon'; dvh(i).plan = 'clinical'; end
save(fullfile(outDir, 'photon_eval_v1.mat'), 'rows', 'dvh');

%% 4  Table + constraint check -----------------------------------------------
cols = struct('machine', {'photon'}, 'plan', {'clinical'}, 'label', {'Photon (delivered)'});
K = cns_toCell(study.summary.keyMetrics);
K = K(cellfun(@(k) any(strcmp(k.quantity, {'RBE11'})), K));       % physical-dose metrics only
cns_compareTable(rows, K, cols, fullfile(outDir, 'table_photon_v1'), ...
    cns_getfield(study, 'constraintsCheck.list', {}));
type(fullfile(outDir, 'table_photon_v1.md'));

%% 5  Figures ---------------------------------------------------------------
rx = caseCfg.prescription.totalDoseGyRBE;
ser = struct('machine', 'photon', 'plan', 'clinical', 'quantity', 'RBE11', 'label', 'Photon (delivered)', ...
    'color', cns_okabeIto('photon'), 'style', '-');
dvhS = cns_toCell(study.figures.dvhStructures);
dvhS = dvhS(cellfun(@(s) any(strcmp({dvh.structure}, s)), dvhS));
set(0, 'DefaultFigureVisible', 'off');
cns_figCurves(dvh, 'dose', 'volPct', dvhS, ser, 'Dose (Gy)', 'Volume (%)', [0 1.15 * rx], ...
    fullfile(figDir, 'fig_photon_DVH_v1'), '');
set(0, 'DefaultFigureVisible', 'on');
fprintf('[UHN_CNS_01] figures in %s\n', figDir);

%% 6  GUI ----------------------------------------------------------------------
% resultGUI holds the photon dose PER FRACTION (matRad convention) and the total.
if S.openGUI
    clear pln stf dij resultGUI                       % stale objects crash matRadGUI
    t = refDose.toTotal; if isnan(t), t = caseCfg.prescription.numOfFractions; end
    resultGUI = struct();
    resultGUI.physicalDose = refDose.cube * t / caseCfg.prescription.numOfFractions;
    resultGUI.photon_dose_total = refDose.cube * t;
    matRad_cfg.disableGUI = false;
    matRadGUI;
end
