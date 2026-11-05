%% UHN_CNS_02b_beamScreen_slices_v2.m
% BEAM SCREEN, SLICES: see the screened plans as a planner would, in the TPS (RayStation-like)
% style: colour wash + isodose lines + structure outlines, on axial, sagittal and coronal slices.
%
% One figure per plane (and per centre structure), same slice for every panel:
%     Photon (delivered) | A | B
%     C                  | D |
% Isodoses 50 / 95 / 107 / 120 % of Rx: the 120 % line shows where the hot spots sit
% (e.g. inside the chiasm / optic nerve where they overlap the PTV).
%
% Nothing is re-optimised: it reads what UHN_CNS_02 (S.runBeamScreen = true) saved in
%   <workDir>/<alias>/beamScreen_v1/<alias>_scr<k>/  (prepared.mat, dij_Generic.mat, plans/)
% and normalises each plan exactly as the screen table did.
% Output: <workDir>/<alias>/beamScreen_v1/figures_slices_v2/fig_screen_<plane>_<centre>_v2.pdf/.png
% Time: about 1 min per beam set (loads one dij at a time) + a few s per figure.
% PATIENT IMAGES: keep the output on hospital storage (the default crop to the brain hides the face).

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';
S.studyFile = '';                            % '' = CNS_ESTRO/study_UHN_CNS_v1.json (beam set names)
S.planes    = {'axial', 'sagittal', 'coronal'};
S.centreOn  = {'CTV', 'Chiasm'};             % one set of figures through the centre of each
S.isoPct    = [50 95 107 120];               % isodose lines (% of Rx)
S.washMinPct = 50;                           % lowest dose shown in colour
S.background = 'dark';                       % 'dark' (TPS look) or 'light' (print)
S.crop      = 'Brain';                       % 'Brain' hides the face; 'Body' shows all
S.showPhoton = true;                         % first panel = delivered photon plan
S.contours = struct( ...
    'name',  {'GTV', 'CTV', 'PTV', 'OpticNerve_L', 'OpticNerve_R', 'Chiasm', 'Brainstem', 'Eye_L', 'Eye_R', 'Body'}, ...
    'color', {[0.15 0.55 0.15], [0.85 0.10 0.85], [0.30 0.60 1.00], [1.00 0.85 0.10], [1.00 0.85 0.10], ...
              [0.20 0.95 0.95], [0.95 0.60 0.20], [0.60 0.80 1.00], [0.60 0.80 1.00], [0.20 0.75 0.20]}, ...
    'label', {'', 'CTV', '', 'ON_L', '', 'Chiasm', 'Brainstem', '', '', ''}, ...
    'labelOffsetMm', {[0 0], [-22 -18], [0 0], [25 -10], [0 0], [18 -20], [20 -15], [0 0], [0 0], [0 0]});

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
scrDir  = fullfile(caseDir, 'beamScreen_v1');
outDir  = fullfile(scrDir, 'figures_slices_v2');
if ~exist(outDir, 'dir'), mkdir(outDir); end
scrStudyFile = fullfile(scrDir, 'study_beamScreen.json');
if ~exist(scrStudyFile, 'file')
    error('No beam screen found in %s. Run UHN_CNS_02 with S.runBeamScreen = true first.', scrDir);
end
scrStudy = cns_readJson(scrStudyFile);           % the settings the screen really used
machine = scrStudy.machines; if iscell(machine), machine = machine{1}; end
rx  = caseCfg.prescription.totalDoseGyRBE;
nFx = caseCfg.prescription.numOfFractions;
sets = cns_toCell(study.beamCandidates.sets);

%% 1  Dose cubes of every screened set (one dij in memory at a time) --------------
pan = struct('cube', {}, 'type', {}, 'title', {}, 'row', {}, 'col', {}, 'doseForMask', {});
ct = []; cst = [];
slot = 0;
if S.showPhoton, slot = 1; end
for k = 1:numel(sets)
    d = fullfile(scrDir, sprintf('%s_scr%d', caseCfg.alias, k));
    fW = fullfile(d, 'plans', [machine '__doseOnly_w.mat']);
    fD = fullfile(d, ['dij_' machine '.mat']);
    fP = fullfile(d, 'prepared.mat');
    if ~exist(fW, 'file') || ~exist(fD, 'file') || ~exist(fP, 'file')
        fprintf(2, '[UHN_CNS_02b] set %s: no screen result in %s, skipped\n', sets{k}.id, d);
        continue;
    end
    P = load(fP);
    if isempty(ct)
        ct = P.ct; cst = P.cst;
        if S.showPhoton
            t = P.refDose.toTotal; if isnan(t), t = nFx; end
            pan(end+1) = struct('cube', double(P.refDose.cube) * t, 'type', 'dose', ...
                'title', 'Photon (delivered)', 'row', 1, 'col', 1, 'doseForMask', []); %#ok<SAGROW>
        end
    end
    W = load(fW);
    D = load(fD, 'dij');
    cstD = matRad_resizeCstToGrid(P.cst, D.dij.ctGrid.x, D.dij.ctGrid.y, D.dij.ctGrid.z, ...
        D.dij.doseGrid.x, D.dij.doseGrid.y, D.dij.doseGrid.z);
    c = caseCfg; c.alias = sprintf('%s_scr%d', caseCfg.alias, k);
    w = cns_normaliseWeights(W.w, D.dij, cstD, c, scrStudy, P.refDose, P.cst);
    R = matRad_calcCubes(w, D.dij);
    clear D
    slot = slot + 1;
    lab = cns_getfield(sets{k}, 'label', sets{k}.id);
    lab = strtrim(regexprep(lab, '[\(,].*$', ''));            % short title: text before ( or ,
    lab = sprintf('%s: %s', sets{k}.id(1), lab);
    pan(end+1) = struct('cube', 1.1 * nFx * R.physicalDose, 'type', 'dose', 'title', lab, ...
        'row', ceil(slot / 3), 'col', mod(slot - 1, 3) + 1, 'doseForMask', []); %#ok<SAGROW>
    fprintf('[UHN_CNS_02b] set %s: max dose %.1f Gy(RBE) (%.0f %% of Rx)\n', sets{k}.id, ...
        max(pan(end).cube(:)), 100 * max(pan(end).cube(:)) / rx);
end
if numel(pan) < 1 + S.showPhoton
    error('No screened plan found in %s.', scrDir);
end

%% 2  Figures: one per plane and centre structure --------------------------------
set(0, 'DefaultFigureVisible', 'off');
for c = 1:numel(S.centreOn)
    if ~any(strcmp(cst(:,2), S.centreOn{c})) || isempty(cst{find(strcmp(cst(:,2), S.centreOn{c}), 1), 4}{1})
        fprintf(2, '[UHN_CNS_02b] %s not found or empty: no figures centred on it\n', S.centreOn{c});
        continue;
    end
    for i = 1:numel(S.planes)
        o = struct('plane', S.planes{i}, 'centreOn', S.centreOn{c}, 'rxGy', rx, ...
            'washMinPct', S.washMinPct, 'isoPct', S.isoPct, 'contours', S.contours, ...
            'crop', S.crop, 'background', S.background);
        base = fullfile(outDir, sprintf('fig_screen_%s_%s_v2', S.planes{i}, S.centreOn{c}));
        cns_figSliceRS(ct, cst, pan, o, base);
        fprintf('[UHN_CNS_02b] %s.png\n', base);
    end
end
set(0, 'DefaultFigureVisible', 'on');
fprintf('[UHN_CNS_02b] done: %s\n', outDir);
