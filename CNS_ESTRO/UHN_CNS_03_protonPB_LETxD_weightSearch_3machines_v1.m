%% UHN_CNS_03_protonPB_LETxD_weightSearch_3machines_v1.m
% CODE 3 of 3: LET x dose optimised proton planning, three machines, with a gradient-descent
% search for the weight of the LET x dose term (start w = 0.1).
%
% For each machine (DWA_proton_v2, Generic, IBA_cyclotron):
%   - the dose + LET influence matrix and the dose-only plan from UHN_CNS_02 are reused;
%   - LET objective: SquaredOverdosing(w, ref) on matRad's 'LETd' quantity, which is
%     LET x dose (dij.mLETDose * w), in the optic nerves, chiasm, brainstem and the
%     0-15 mm brain ring outside the PTV (study.letSearch.letStructures). ref = the
%     5th percentile of LET x dose in the CTV of the dose-only plan (total course),
%     so only OAR voxels hotter than the target's low end are penalised;
%   - cns_letWeightSearch: steepest descent on log10(w) of
%       J(w) = mean_s [D x LETd D2]_s(w) / [D x LETd D2]_s(dose-only) + lambda * coverage penalty
%     (central finite differences, Armijo backtracking, w in [0.001, 3], <= 14 plans);
%   - the best plan is evaluated with every metric next to the dose-only and photon plans.
%
% RUN ON THE UHN PC ONLY, after UHN_CNS_01 and UHN_CNS_02, with the same study file.
%   setenv('MATRAD_ROOT', 'X:\matRadDWA_LET');
%   run(fullfile(getenv('MATRAD_ROOT'), 'CNS_ESTRO', 'UHN_CNS_03_protonPB_LETxD_weightSearch_3machines_v1.m'))
% Time: up to maxEvals optimisations per machine (each ~2-10 min at 3 mm): plan on an
% overnight run for three machines. Each machine's result is saved as soon as it finishes.
% Output: <workDir>/<alias>/LETxD_search_v1/ : search_<machine>_v1.csv (every weight tried),
%         _trace.txt (the descent), metrics.csv, table_LETxD_v1.md/.csv (+ constraint check),
%         plan_<machine>_LETxDbest_v1.mat (spot weights for MCsquare recalculation), figures/

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';
S.studyFile = '';             % '' = CNS_ESTRO/study_UHN_CNS_v1.json (must match UHN_CNS_02)
S.machines  = {};             % {} = study.machines
S.w0        = 0.1;            % starting weight of the LET x dose term

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
matRad_cfg.logLevel = 1;
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
study.letSearch.w0 = S.w0;
machines = cns_toCell(study.machines);
if ~isempty(S.machines), machines = S.machines; end
caseDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias);
outDir = fullfile(caseDir, 'LETxD_search_v1');
figDir = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
cns_writeJson(fullfile(outDir, 'study_LETxD_v1_used.json'), study);
matRad_cfg.defaults.propOpt.maxIter = cns_getfield(study, 'optimisation.maxIter', 500);

P = load(fullfile(caseDir, 'prepared.mat'));        % made by UHN_CNS_02 (cns_runCase)
refDose = []; if isfield(P, 'refDose'), refDose = P.refDose; end

%% 1  Weight search per machine ---------------------------------------------
allRows = {};
search = struct();
for m = 1:numel(machines)
    mc = machines{m};
    dijFile = fullfile(caseDir, ['dij_' mc '.mat']);
    wFile = fullfile(caseDir, 'plans', [mc '__doseOnly_w.mat']);
    if ~exist(dijFile, 'file') || ~exist(wFile, 'file')
        error('Run UHN_CNS_02 first: %s or %s is missing.', dijFile, wFile);
    end
    fprintf('\n[UHN_CNS_03] ===== %s: LET x dose weight search from w = %g\n', mc, S.w0);
    D = load(dijFile, 'pln', 'dij');
    pln = D.pln;
    tmp = cns_buildPln(P.ct, P.cst, caseCfg, study, mc);
    pln.multScen = tmp.multScen;                       % classdef objects are not restored from .mat in Octave
    W = load(wFile, 'w');
    R = cns_letWeightSearch(D.dij, P.cst, pln, caseCfg, study, W.w, refDose);

    % save this machine immediately (an overnight run may be interrupted)
    tag = strrep(mc, '-', '_');
    cns_writeCsvTable(fullfile(outDir, ['search_' tag '_v1.csv']), R.header, R.table);
    fid = fopen(fullfile(outDir, ['search_' tag '_v1_trace.txt']), 'w');
    fprintf(fid, 'LET x dose reference: %.3f keV/um*Gy (total course)\n', R.ref);
    fprintf(fid, '%s\n', R.trace{:});
    fprintf(fid, 'best w = %.5g, J = %.5f (dose-only J = %.5f), %d evaluations. %s\n', ...
        R.best.w, R.best.J, R.base.J, R.nEvaluations, R.note);
    fclose(fid);
    wBest = R.best.weights; wDoseNorm = R.base.w; ref = R.ref; bestW = R.best.w; %#ok<NASGU>
    save(fullfile(outDir, ['plan_' tag '_LETxDbest_v1.mat']), 'wBest', 'wDoseNorm', 'ref', 'bestW');

    n0 = size(R.base.rows, 1); n1 = size(R.best.rows, 1);
    rowsM = [repmat({caseCfg.alias, mc, 'doseOnly'}, n0, 1), R.base.rows; ...
        repmat({caseCfg.alias, mc, 'LETxD_best'}, n1, 1), R.best.rows];
    cns_writeCsvTable(fullfile(outDir, ['metrics_' tag '_v1.csv']), ...
        {'alias', 'machine', 'plan', 'structure', 'quantity', 'metric', 'value', 'unit'}, rowsM);   % per machine: UHN_CNS_05 merges runs
    allRows = [allRows; rowsM]; %#ok<AGROW>
    search.(tag) = R;
    clear D
end

%% 2  Photon reference rows, combined metrics ---------------------------------
if ~isempty(refDose)
    E = cns_evaluateReferenceDose(P.ct, P.cst, refDose, caseCfg, study);
    allRows = [allRows; repmat({caseCfg.alias, 'photon', 'clinical'}, size(E.rows, 1), 1), E.rows];
end
cns_writeCsvTable(fullfile(outDir, 'metrics.csv'), {'alias', 'machine', 'plan', 'structure', 'quantity', 'metric', 'value', 'unit'}, allRows);

cols = struct('machine', {}, 'plan', {}, 'label', {});
if ~isempty(refDose), cols(end+1) = struct('machine', 'photon', 'plan', 'clinical', 'label', 'Photon'); end
for m = 1:numel(machines)
    tag = strrep(machines{m}, '-', '_');
    cols(end+1) = struct('machine', machines{m}, 'plan', 'doseOnly', 'label', [machines{m} ' dose']); %#ok<SAGROW>
    cols(end+1) = struct('machine', machines{m}, 'plan', 'LETxD_best', ...
        'label', sprintf('%s LETxD w=%.3g', machines{m}, search.(tag).best.w)); %#ok<SAGROW>
end
cns_compareTable(allRows, study.summary.keyMetrics, cols, fullfile(outDir, 'table_LETxD_v1'), ...
    cns_getfield(study, 'constraintsCheck.list', {}));
type(fullfile(outDir, 'table_LETxD_v1.md'));

%% 3  Journal figures ---------------------------------------------------------
set(0, 'DefaultFigureVisible', 'off');
% (a) the search: J, burden ratio and CTV D98 against w, one colour per machine
f = cns_journalFigure(17.5, 5.5);
labels = {'Objective J', 'D \times LET_d burden (rel. to dose-only)', 'CTV D_{98} (Gy(RBE))'};
colsIdx = [4 5 7];
csv = {};
for p = 1:3
    ax = subplot(1, 3, p); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'XScale', 'log', 'GridAlpha', 0.15);
    for m = 1:numel(machines)
        tag = strrep(machines{m}, '-', '_');
        Tb = search.(tag).table;
        w = cell2mat(Tb(:, 2)); y = cell2mat(Tb(:, colsIdx(p)));
        sel = w > 0; [ws, o] = sort(w(sel)); ys = y(sel); ys = ys(o);
        plot(ax, ws, ys, 'o-', 'Color', cns_okabeIto(machines{m}), 'MarkerSize', 3, ...
            'MarkerFaceColor', cns_okabeIto(machines{m}), 'DisplayName', strrep(machines{m}, '_', ' '));
        y0 = y(~sel);
        if ~isempty(y0)
            plot(ax, xlim(ax), [y0(1) y0(1)], ':', 'Color', cns_okabeIto(machines{m}), 'HandleVisibility', 'off');
        end
        ib = find(abs(w - search.(tag).best.w) < 1e-12, 1);
        if ~isempty(ib)
            plot(ax, w(ib), y(ib), 'p', 'Color', cns_okabeIto(machines{m}), 'MarkerSize', 8, ...
                'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
        end
        if p == 1
            for i = 1:numel(w), csv(end+1, :) = [{machines{m}}, Tb(i, 1:9)]; end %#ok<SAGROW>
        end
    end
    xlabel(ax, 'LET \times dose weight w'); ylabel(ax, labels{p});
    if p == 3
        plot(ax, xlim(ax), 0.95 * caseCfg.prescription.totalDoseGyRBE * [1 1], 'k--', 'HandleVisibility', 'off');
    end
    if p == 1, legend(ax, 'Location', 'best'); end
end
cns_saveJournalFigure(f, fullfile(figDir, 'fig_LETxD_weightSearch_v1'), ...
    {'machine', 'evaluation', 'w', 'log10w', 'J', 'burdenRatio', 'Pcov', 'CTV_D98', 'CTV_D2', 'CTV_D95'}, csv);

% (b) DVH and LVH: dose-only (solid) vs LET x dose best (dashed), photon in black
dvhAll = []; lvhAll = [];
for m = 1:numel(machines)
    tag = strrep(machines{m}, '-', '_');
    D = load(fullfile(caseDir, ['dij_' machines{m} '.mat']), 'dij');
    dvhM = []; lvhM = [];
    for pl = {'doseOnly', 'LETxD_best'}
        if strcmp(pl{1}, 'doseOnly'), w = search.(tag).base.w; else, w = search.(tag).best.weights; end
        E2 = cns_evaluatePlan(D.dij, P.cst, w, caseCfg, study);
        for i = 1:numel(E2.dvh), E2.dvh(i).machine = machines{m}; E2.dvh(i).plan = pl{1}; end
        for i = 1:numel(E2.lvh), E2.lvh(i).machine = machines{m}; E2.lvh(i).plan = pl{1}; end
        dvhM = [dvhM, orderfields(E2.dvh)]; lvhM = [lvhM, orderfields(E2.lvh)]; %#ok<AGROW>
    end
    save(fullfile(outDir, ['dvh_lvh_' tag '_v1.mat']), 'dvhM', 'lvhM');   % per machine: UHN_CNS_05 merges runs
    dvhAll = [dvhAll, dvhM]; lvhAll = [lvhAll, lvhM]; %#ok<AGROW>
    clear D
end
if ~isempty(refDose)
    for i = 1:numel(E.dvh), E.dvh(i).machine = 'photon'; E.dvh(i).plan = 'clinical'; end
    dvhAll = [dvhAll, orderfields(E.dvh)];
end
ser = struct('machine', {}, 'plan', {}, 'quantity', {}, 'label', {}, 'color', {}, 'style', {});
if ~isempty(refDose)
    ser(end+1) = struct('machine', 'photon', 'plan', 'clinical', 'quantity', 'RBE11', 'label', 'Photon', ...
        'color', cns_okabeIto('photon'), 'style', '-');
end
serL = ser([]);
for m = 1:numel(machines)
    nm = strrep(machines{m}, '_', ' ');
    ser(end+1) = struct('machine', machines{m}, 'plan', 'doseOnly', 'quantity', 'RBE11', ...
        'label', [nm ' dose-only'], 'color', cns_okabeIto(machines{m}), 'style', '-'); %#ok<SAGROW>
    ser(end+1) = struct('machine', machines{m}, 'plan', 'LETxD_best', 'quantity', 'RBE11', ...
        'label', [nm ' LET\timesD'], 'color', cns_okabeIto(machines{m}), 'style', '--'); %#ok<SAGROW>
    serL(end+1) = struct('machine', machines{m}, 'plan', 'doseOnly', 'quantity', '', ...
        'label', [nm ' dose-only'], 'color', cns_okabeIto(machines{m}), 'style', '-'); %#ok<SAGROW>
    serL(end+1) = struct('machine', machines{m}, 'plan', 'LETxD_best', 'quantity', '', ...
        'label', [nm ' LET\timesD'], 'color', cns_okabeIto(machines{m}), 'style', '--'); %#ok<SAGROW>
end
rx = caseCfg.prescription.totalDoseGyRBE;
cns_figCurves(dvhAll, 'dose', 'volPct', cns_toCell(study.figures.dvhStructures), ser, ...
    'Dose (Gy(RBE))', 'Volume (%)', [0 1.15 * rx], fullfile(figDir, 'fig_LETxD_DVH_v1'), '');
cns_figCurves(lvhAll, 'letd', 'volCc', cns_toCell(study.figures.lvhStructures), serL, ...
    'LET_d (keV/\mum)', ...
    'Volume (cm^3)', [0 10], fullfile(figDir, 'fig_LETxD_LVH_v1'), ...
    sprintf('LET_d in voxels with RBE 1.1 dose >= %g Gy(RBE)', study.evaluation.letdMinDoseGy));
set(0, 'DefaultFigureVisible', 'on');
fprintf('\n[UHN_CNS_03] done: %s\n', outDir);
