%% UHN_CNS_02_protonPB_doseOnly_3machines_v1.m
% CODE 2 of 3: physical-dose-based proton planning (constant RBE 1.1, pencil beam)
% with three machines: DWA_proton_v2 (3 mm spots), Generic (~HIT synchrotron, 5 mm)
% and IBA_cyclotron (5 mm, 60 mm range shifter for shallow spots).
%
%   A  beam-geometry figure for the candidate beam sets (check angles before any dose)
%   B  optional quick beam screen (Generic, coarse grid, dose-only) -> table to choose a set
%   C  dose-only plans for the three machines with the beams of the case file,
%      normalised to the photon plan's PTV D95 (study.normalisation), all metrics,
%      photon plan evaluated alongside
%   D  journal-style comparison: DVH, LET-volume histograms, table with constraint check
%
% Needs UHN_CNS_01 first (it creates the imported .mat).
% RUN ON THE UHN PC ONLY. From the matRad root:
%   setenv('MATRAD_ROOT', 'X:\matRadDWA_LET');
%   run(fullfile(getenv('MATRAD_ROOT'), 'CNS_ESTRO', 'UHN_CNS_02_protonPB_doseOnly_3machines_v1.m'))
% Time (3 mm grid): dose + LET influence matrices ~5-20 min per machine (DWA longest),
% optimisation ~2-10 min per plan. Memory: the DWA dij is the largest (several GB).
% Output: <workDir>/<alias>/  (cns_runCase: dij_<machine>.mat, plans/, metrics.csv, dvh.csv, lvh.csv)
%         <workDir>/<alias>/doseOnly_v1/  (tables and journal figures)
% The dij caches and the dose-only weights are reused by UHN_CNS_03.

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';
S.studyFile = '';             % '' = CNS_ESTRO/study_UHN_CNS_v1.json
S.runBeamScreen = true;      % true: section B (about 5-10 min per candidate set)
S.machines = {};              % {} = study.machines; e.g. {'Generic'} for a first test

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
dataClass = cns_checkPrivacy(caseCfg, S.caseFile);
caseDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias);
outDir = fullfile(caseDir, 'doseOnly_v1');
figDir = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~isempty(S.machines), study.machines = S.machines; end
set(0, 'DefaultFigureVisible', 'off');

%% A  Beam geometry -------------------------------------------------------------
L = load(cns_expandPath(caseCfg.source.path));
[cstP, ~] = cns_prepareStructures(L.ct, L.cst, caseCfg, study);
sets = cns_toCell(study.beamCandidates.sets);
inCase = struct('id', ['case_' cns_getfield(caseCfg, 'beams.setId', 'beams')], ...
    'gantryAngles', caseCfg.beams.gantryAngles, 'couchAngles', caseCfg.beams.couchAngles);
cns_figBeamGeometry(L.ct, cstP, [{inCase}, sets], fullfile(figDir, 'fig_beamGeometry_v1'));
fprintf('[UHN_CNS_02] check %s.png before going on: arrows point from the source to the target.\n', ...
    fullfile(figDir, 'fig_beamGeometry_v1'));

%% B  Optional beam screen (Generic, coarse, dose-only) -----------------------
if S.runBeamScreen
    B = study.beamScreen;
    scr = study;
    scr.machines = {B.machine};
    scr.doseCalc.doseGridMm = B.doseGridMm;
    scr.stf.bixelWidthMm.(B.machine) = B.bixelWidthMm;
    scr.stf.longitudinalSpotSpacingMm = B.longitudinalSpotSpacingMm;
    scr.optimisation.maxIter = B.maxIter;
    scrDir = fullfile(caseDir, 'beamScreen_v1');
    if ~exist(scrDir, 'dir'), mkdir(scrDir); end
    scrStudy = fullfile(scrDir, 'study_beamScreen.json');
    cns_writeJson(scrStudy, scr);
    allRows = {};
    cols = struct('machine', {}, 'plan', {}, 'label', {});
    for k = 1:numel(sets)
        c = caseCfg;
        c.alias = sprintf('%s_scr%d', caseCfg.alias, k);   % alias max 32 characters
        c.workDir = scrDir;
        c.beams = struct('gantryAngles', sets{k}.gantryAngles, 'couchAngles', sets{k}.couchAngles);
        cf = fullfile(scrDir, [c.alias '.json']);
        cns_writeJson(cf, c);
        o = cns_runCase(cf, scrStudy);
        r = o.rows;
        r(:,2) = {sets{k}.id};                     % column = beam set
        allRows = [allRows; r]; %#ok<AGROW>
        cols(end+1) = struct('machine', sets{k}.id, 'plan', 'doseOnly', 'label', sets{k}.id); %#ok<AGROW>
    end
    K = cns_toCell(study.summary.keyMetrics);
    K = K(cellfun(@(k) any(strcmp(k.quantity, {'RBE11', 'DxLETd'})), K));
    cns_compareTable(allRows, K, cols, fullfile(scrDir, 'table_beamScreen_v1'), ...
        cns_getfield(study, 'constraintsCheck.list', {}));
    type(fullfile(scrDir, 'table_beamScreen_v1.md'));
    fprintf(['[UHN_CNS_02] Beam screen done. Pick a set, copy its gantry/couch angles into the case file ' ...
        '(beams), then run this script again with S.runBeamScreen = false.\n']);
    return;
end

%% C  Dose-only plans, three machines -------------------------------------------
st = study;
st.plans = {struct('id', 'doseOnly', 'let', struct('mode', 'none'))};
runStudy = fullfile(caseDir, 'study_doseOnly_v1_used.json');     % the exact settings used
if ~exist(caseDir, 'dir'), mkdir(caseDir); end
cns_writeJson(runStudy, st);
out = cns_runCase(S.caseFile, runStudy);

%% D  Tables and journal figures ---------------------------------------------
rows = out.rows;
machines = cns_toCell(st.machines);
cols = struct('machine', {'photon'}, 'plan', {'clinical'}, 'label', {'Photon'});
for m = 1:numel(machines)
    cols(end+1) = struct('machine', machines{m}, 'plan', 'doseOnly', 'label', machines{m}); %#ok<AGROW>
end
cns_compareTable(rows, study.summary.keyMetrics, cols, fullfile(outDir, 'table_doseOnly_v1'), ...
    cns_getfield(study, 'constraintsCheck.list', {}));
type(fullfile(outDir, 'table_doseOnly_v1.md'));

rx = caseCfg.prescription.totalDoseGyRBE;
ser = struct('machine', {}, 'plan', {}, 'quantity', {}, 'label', {}, 'color', {}, 'style', {});
ser(end+1) = struct('machine', 'photon', 'plan', 'clinical', 'quantity', 'RBE11', ...
    'label', 'Photon (delivered)', 'color', cns_okabeIto('photon'), 'style', '-');
for m = 1:numel(machines)
    ser(end+1) = struct('machine', machines{m}, 'plan', 'doseOnly', 'quantity', 'RBE11', ...
        'label', [strrep(machines{m}, '_', ' ') ' (RBE 1.1)'], 'color', cns_okabeIto(machines{m}), 'style', '-'); %#ok<SAGROW>
    ser(end+1) = struct('machine', machines{m}, 'plan', 'doseOnly', 'quantity', 'McN', ...
        'label', [strrep(machines{m}, '_', ' ') ' (McNamara)'], 'color', cns_okabeIto(machines{m}), 'style', '--'); %#ok<SAGROW>
end
dvhS = cns_toCell(study.figures.dvhStructures);
cns_figCurves(out.dvh, 'dose', 'volPct', dvhS, ser, 'Dose (Gy(RBE))', 'Volume (%)', [0 1.15 * rx], ...
    fullfile(figDir, 'fig_doseOnly_DVH_v1'), '');
serL = ser(~strcmp({ser.machine}, 'photon') & strcmp({ser.quantity}, 'RBE11'));
for j = 1:numel(serL), serL(j).quantity = ''; serL(j).label = strrep(machines{j}, '_', ' '); end
cns_figCurves(out.lvh, 'letd', 'volCc', cns_toCell(study.figures.lvhStructures), serL, ...
    'LET_d (keV/\mum)', ...
    'Volume (cm^3)', [0 10], fullfile(figDir, 'fig_doseOnly_LVH_v1'), ...
    sprintf('LET_d in voxels with RBE 1.1 dose >= %g Gy(RBE)', study.evaluation.letdMinDoseGy));
set(0, 'DefaultFigureVisible', 'on');
fprintf('[UHN_CNS_02] tables and figures in %s\n', outDir);

%% E  Optional: look at one plan in matRadGUI -------------------------------
% m = 'DWA_proton_v2';
% Lp = load(fullfile(caseDir, 'prepared.mat')); Ld = load(fullfile(caseDir, ['dij_' m '.mat']));
% Lw = load(fullfile(caseDir, 'plans', [m '__doseOnly_w.mat']));
% ct = Lp.ct; cst = Lp.cst; pln = Ld.pln; stf = Ld.stf;
% resultGUI = matRad_calcCubes(Lw.w, Ld.dij);
% resultGUI = cns_addTotalsForGUI(resultGUI, ct, cst, pln.numOfFractions, study);
% matRad_cfg.disableGUI = false; matRadGUI;
