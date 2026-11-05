%% UHN_CNS_05_figure2_DVH_doseLET_v1.m
% FIGURE 2: (a) DVHs of the CTV and OARs, (b) LET-volume histograms (optional), (c) dose-LET bars, for
%   photon (delivered) + DWA / Generic / IBA, each dose-only and LET x D-optimised.
%
% Reads (nothing is recomputed):
%   photon_v1/photon_eval_v1.mat                    (UHN_CNS_01)
%   LETxD_search_v1/metrics_<machine>_v1.csv        (UHN_CNS_03, one per machine)
%   LETxD_search_v1/dvh_lvh_<machine>_v1.mat        (UHN_CNS_03, one per machine)
% so UHN_CNS_03 can be run machine by machine (e.g. DWA first) and this figure
% uses whatever machines are finished.
% Output: <workDir>/<alias>/figures_money_v1/fig2_DVH_doseLET_v1.pdf/.png + csv,
%         table_fig2_v1.md/.csv (all plotted numbers and the constraint check)

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';
S.studyFile = '';                     % '' = CNS_ESTRO/study_UHN_CNS_v1.json
S.machines  = {'DWA_proton_v2', 'Generic', 'IBA_cyclotron'};   % order of the colours / bars
S.dvhStructures = {'CTV', 'OpticNerve_L', 'Chiasm', 'Brainstem', 'Brain_healthy'};
S.showLVH = true;             % row (b): LET-volume histograms (drop it for the one-figure ESTRO version)
S.lvhStructures = {'OpticNerve_L', 'Chiasm', 'Brainstem', 'BrainRing_0_15'};
S.barMetrics = struct( ...
    'structure', {'OpticNerve_L', 'Chiasm', 'OpticNerve_L', 'Chiasm', 'BrainRing_0_15'}, ...
    'quantity',  {'DxLETd', 'DxLETd', 'McN', 'McN', 'hotSpot'}, ...
    'metric',    {'D2', 'D2', 'D0_03cc', 'D0_03cc', 'V_D20_L4_cc'}, ...
    'label',     {'Optic nerve L', 'Chiasm', 'Optic nerve L', 'Chiasm', 'Brain ring 0-15 mm'}, ...
    'sublabel',  {'D\timesLET_d D_{2%} (Gy\cdotkeV/\mum)', 'D\timesLET_d D_{2%} (Gy\cdotkeV/\mum)', ...
                  'McNamara D_{0.03cm^3} (Gy(RBE))', 'McNamara D_{0.03cm^3} (Gy(RBE))', ...
                  'V(D\geq20 Gy, LET_d\geq4) (cm^3)'});

matRadRoot = getenv('MATRAD_ROOT');
if isempty(matRadRoot)
    here = mfilename('fullpath');
    if isempty(here), here = which('cns_runCase'); end
    if isempty(here), here = fullfile(pwd, 'x'); end
    matRadRoot = fileparts(fileparts(here));
end
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
srcDir = fullfile(caseDir, 'LETxD_search_v1');
outDir = fullfile(caseDir, 'figures_money_v1');
if ~exist(outDir, 'dir'), mkdir(outDir); end
rx = caseCfg.prescription.totalDoseGyRBE;

%% 1  Collect the data -----------------------------------------------------------
Ph = load(fullfile(caseDir, 'photon_v1', 'photon_eval_v1.mat'));     % rows, dvh
rows = Ph.rows;
dvh = orderfields(Ph.dvh);
lvh = [];
have = {};
for m = 1:numel(S.machines)
    tag = strrep(S.machines{m}, '-', '_');
    fM = fullfile(srcDir, ['metrics_' tag '_v1.csv']);
    fD = fullfile(srcDir, ['dvh_lvh_' tag '_v1.mat']);
    if ~exist(fM, 'file') || ~exist(fD, 'file')
        fprintf('[UHN_CNS_05] %s: no UHN_CNS_03 result yet, left out\n', S.machines{m});
        continue;
    end
    fid = fopen(fM); fgetl(fid);
    C = textscan(fid, '%s %s %s %s %s %s %f %s', 'Delimiter', ',');
    fclose(fid);
    rows = [rows; [C{1:6}, num2cell(C{7}), C{8}]]; %#ok<AGROW>
    L = load(fD);
    dvh = [dvh, orderfields(L.dvhM)]; %#ok<AGROW>
    lvh = [lvh, orderfields(L.lvhM)]; %#ok<AGROW>
    have{end+1} = S.machines{m}; %#ok<AGROW>
end
if isempty(have), error('No UHN_CNS_03 results found in %s.', srcDir); end

%% 2  Figure + table -------------------------------------------------------------
set(0, 'DefaultFigureVisible', 'off');
base = fullfile(outDir, 'fig2_DVH_doseLET_v1');
o = struct('letdMinDoseGy', study.evaluation.letdMinDoseGy);
if S.showLVH, o.lvh = lvh; o.lvhStructs = S.lvhStructures; end
cns_figDVHandBars(dvh, rows, have, S.dvhStructures, S.barMetrics, rx, base, o);
set(0, 'DefaultFigureVisible', 'on');

cols = struct('machine', {'photon'}, 'plan', {'clinical'}, 'label', {'Photon'});
for m = 1:numel(have)
    cols(end+1) = struct('machine', have{m}, 'plan', 'doseOnly', 'label', [have{m} ' dose']); %#ok<SAGROW>
    cols(end+1) = struct('machine', have{m}, 'plan', 'LETxD_best', 'label', [have{m} ' LETxD']); %#ok<SAGROW>
end
cns_compareTable(rows, study.summary.keyMetrics, cols, fullfile(outDir, 'table_fig2_v1'), ...
    cns_getfield(study, 'constraintsCheck.list', {}));
fprintf('[UHN_CNS_05] %s.png  (machines: %s)\n', base, strjoin(have, ', '));
