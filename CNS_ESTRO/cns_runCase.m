function out = cns_runCase(caseFile, studyFile, opts)
% cns_runCase  Plan one CNS case with every machine x plan variant of the study and evaluate it.
%
%   out = cns_runCase('cases/CNS_synthetic.json', 'study_ESTRO2027.json')
%   out = cns_runCase(caseFile, studyFile, struct('forceDij', true, 'reusePlans', false))
%
% For each machine in study.machines:
%   1) steering file + dose / LET x dose influence matrices (cached per machine)
%   2) every plan in study.plans (dose-only reference + LET-modified plans)
%   3) normalisation (study.normalisation, e.g. PTV D95 = 95 % Rx)
%   4) cns_evaluatePlan: CTV/PTV DVH metrics in physical, RBE 1.1 and
%      McNamara dose; OAR dose, LETd, smoothed D x LETd and dose-LET hot spots
%
% Writes to <workDir>/<alias>/ :
%   metrics.csv   long table: alias, machine, plan, structure, quantity, metric, value, unit
%   dvh.csv, lvh.csv, run_info.json, plans/<machine>__<plan>_w.mat, figures/*.png
% and caches prepared.mat and dij_<machine>.mat (large; never commit them).
%
% PRIVACY. cns_checkPrivacy runs first. For a UHN case (dataClass "patient")
% it refuses to run unless the case file, the data and workDir are on
% hospital storage, outside this repository and any AI-connected folder.

if nargin < 3, opts = struct(); end
here = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(studyFile), studyFile = fullfile(here, 'study_ESTRO2027.json'); end

matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;
isOctave = exist('OCTAVE_VERSION', 'builtin') > 0;

caseCfg = cns_readJson(caseFile);
study = cns_readJson(studyFile);
dataClass = cns_checkPrivacy(caseCfg, caseFile);
alias = caseCfg.alias;
caseDir = fullfile(cns_expandPath(cns_getfield(caseCfg, 'workDir', fullfile(pwd, 'CNS_ESTRO_runs'))), alias);
planDir = fullfile(caseDir, 'plans');
if ~exist(planDir, 'dir'), mkdir(planDir); end
fprintf('\n[cns] ===== case %s (%s) -> %s\n', alias, dataClass, caseDir);

forceDij = logical(cns_getfield(opts, 'forceDij', cns_getfield(study, 'forceDij', false)));
reusePlans = logical(cns_getfield(opts, 'reusePlans', cns_getfield(study, 'reusePlans', true)));
maxIter = cns_getfield(study, 'optimisation.maxIter', 500);
matRad_cfg.defaults.propOpt.maxIter = maxIter;

% ---- case geometry
prepFile = fullfile(caseDir, 'prepared.mat');
if exist(prepFile, 'file') && ~forceDij
    S = load(prepFile);
    ct = S.ct; cst = S.cst; structInfo = S.structInfo;
    refDose = []; if isfield(S, 'refDose'), refDose = S.refDose; end
else
    [ct, cst, refDose] = cns_loadCase(caseCfg);
    [cst, structInfo] = cns_prepareStructures(ct, cst, caseCfg, study);
    saveMat(prepFile, isOctave, 'ct', ct, 'cst', cst, 'structInfo', structInfo, 'refDose', refDose);
end

machines = cns_toCell(cns_getfield(opts, 'machines', study.machines));
plans = cns_toCell(study.plans);
hasDirty = exist('matRad_DirtyDose', 'class') == 8;
needDirty = any(cellfun(@(p) strcmp(cns_getfield(p, 'let.mode', 'none'), 'dirtyDose'), plans));
dirtyThr = cns_getfield(study, 'evaluation.dirtyDoseThresholdKeVum', 2);

rows = {};
dvhAll = {}; lvhAll = {};
slices = struct('machine', {}, 'plan', {}, 'cubes', {}, 'masks', {});
info = struct('alias', alias, 'dataClass', dataClass, 'study', cns_getfield(study, 'studyId', ''), ...
    'started', datestr(now, 31), 'structureNotes', {structInfo.notes}, ...
    'brainIsProxy', structInfo.brainIsProxy, 'runs', {{}});

for m = 1:numel(machines)
    machine = machines{m};
    dijFile = fullfile(caseDir, ['dij_' machine '.mat']);
    % (Nov 8) spots only on study.stf.spotTargets (e.g. CTV): the robust margin is added by
    % matRad; spots on the PTV + robust margin double the spot count and run out of memory
    spotTargets = cns_toCell(cns_getfield(study, 'stf.spotTargets', {}));
    spotKey = strjoin(sort(spotTargets), ',');
    useCache = exist(dijFile, 'file') && ~forceDij;
    if useCache
        Sk = whos('-file', dijFile);
        cachedKey = '';
        if any(strcmp({Sk.name}, 'spotKey')), Sc = load(dijFile, 'spotKey'); cachedKey = Sc.spotKey; end
        if ~strcmp(cachedKey, spotKey)
            fprintf('[cns] %s: cached dij has spot targets {%s}, study asks {%s}: recomputing\n', ...
                machine, cachedKey, spotKey);
            useCache = false;
        end
    end
    dijFresh = ~useCache;
    if useCache
        fprintf('[cns] %s: loading cached dij\n', machine);
        S = load(dijFile, 'pln', 'stf', 'dij');
        pln = S.pln; dij = S.dij; clear S
        tmp = cns_buildPln(ct, cst, caseCfg, study, machine);
        pln.multScen = tmp.multScen;     % classdef objects come back as structs in Octave
    else
        pln = cns_buildPln(ct, cst, caseCfg, study, machine);
        fprintf('[cns] %s: steering + dose/LET influence matrices...\n', machine);
        t0 = tic;
        cstStf = cst;
        if ~isempty(spotTargets)
            for iS = 1:size(cstStf, 1)
                if strcmp(cstStf{iS,3}, 'TARGET') && ~any(strcmp(cstStf{iS,2}, spotTargets))
                    cstStf{iS,3} = 'OAR';
                end
            end
            fprintf('[cns] spots placed on: %s (+ matRad robustness margin)\n', spotKey);
        end
        stf = matRad_generateStf(ct, cstStf, pln);
        clear cstStf
        dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
        if needDirty && hasDirty
            dij = matRad_DirtyDose.addDirtyDoseDij(dij, dirtyThr);
        end
        fprintf('[cns] %s: %d spots, dij in %.0f s\n', machine, dij.totalNumOfBixels, toc(t0));
        saveMat(dijFile, isOctave, 'pln', pln, 'stf', stf, 'dij', dij, 'spotKey', spotKey);
    end
    cstD = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
        dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
    masks = struct();
    for sName = {'CTV', 'PTV', 'Brainstem'}
        ix = find(strcmp(cstD(:,2), sName{1}), 1);
        if ~isempty(ix)
            mk = false(dij.doseGrid.dimensions); mk(cstD{ix,4}{1}) = true; masks.(sName{1}) = mk;
        end
    end

    for p = 1:numel(plans)
        P = plans{p};
        mode = cns_getfield(P, 'let.mode', 'none');
        onlyFor = cns_toCell(cns_getfield(P, 'machines', {}));
        if ~isempty(onlyFor) && ~any(strcmp(onlyFor, machine)), continue; end
        if strcmp(mode, 'dirtyDose') && ~(hasDirty && isfield(dij, 'mDirtyDose'))
            fprintf('[cns] %s / %s: skipped (matRad_DirtyDose not on the path)\n', machine, P.id);
            continue;
        end
        tag = [machine '__' P.id];
        wFile = fullfile(planDir, [tag '_w.mat']);
        run = struct('machine', machine, 'plan', P.id, 'status', 'ok', 'notes', {{}}, 'seconds', 0);
        t0 = tic;
        try
            wOk = false;
            if exist(wFile, 'file') && reusePlans && ~dijFresh       % never reuse weights on a new dij
                S = load(wFile, 'w'); w = S.w; wOk = numel(w) == dij.totalNumOfBixels;
            end
            if wOk
                run.notes{end+1} = 'weights reused';
            else
                [cstRun, notes] = cns_buildCst(cst, caseCfg, study, P, pln);
                run.notes = [run.notes, notes];
                fprintf('[cns] %s / %s: optimising\n', machine, P.id);
                [resultGUI, ~, optUsed] = cns_optimize(dij, cstRun, pln, study);
                run.optimizer = optUsed;
                if ~strcmpi(optUsed, 'IPOPT'), run.notes{end+1} = ['optimiser ' optUsed]; end
                w = resultGUI.wUnsequenced;
                run.optimizerStatus = optStatus(resultGUI);
                clear resultGUI                       % 10 dose-grid result sets (9 scenarios)
                saveMat(wFile, isOctave, 'w', w);
            end
            [w, run.normalisationFactor, run.normalisationGoalGy] = ...
                cns_normaliseWeights(w, dij, cstD, caseCfg, study, refDose, cst);
            E = cns_evaluatePlan(dij, cst, w, caseCfg, study);
            n = size(E.rows, 1);
            rows = [rows; repmat({alias, machine, P.id}, n, 1), E.rows]; %#ok<AGROW>
            dvhAll{end+1} = tagRows(E.dvh, machine, P.id); %#ok<AGROW>
            lvhAll{end+1} = tagRows(E.lvh, machine, P.id); %#ok<AGROW>
            slices(end+1) = struct('machine', machine, 'plan', P.id, 'cubes', E.cubes, 'masks', masks); %#ok<AGROW>
        catch err
            run.status = 'failed';
            run.error = err.message;
            fprintf(2, '[cns] %s / %s FAILED: %s\n', machine, P.id, err.message);
            for kS = 1:min(4, numel(err.stack))      % where it failed (function, line)
                fprintf(2, '        at %s line %d\n', err.stack(kS).name, err.stack(kS).line);
            end
        end
        run.seconds = toc(t0);
        info.runs{end+1} = run;
    end
    clear dij
end

% ---- clinical reference plan (e.g. the delivered photon plan), same metrics
if ~isempty(refDose)
    try
        refMachine = cns_getfield(caseCfg, 'reference.label', 'photon');
        E = cns_evaluateReferenceDose(ct, cst, refDose, caseCfg, study);
        n = size(E.rows, 1);
        rows = [rows; repmat({alias, refMachine, 'clinical'}, n, 1), E.rows];
        dvhAll{end+1} = tagRows(E.dvh, refMachine, 'clinical');
        lvhAll{end+1} = tagRows(E.lvh, refMachine, 'clinical');
        mk = struct();
        for sName = {'CTV', 'PTV', 'Brainstem'}
            ix = find(strcmp(cst(:,2), sName{1}), 1);
            if ~isempty(ix), m0 = false(ct.cubeDim); m0(cst{ix,4}{1}) = true; mk.(sName{1}) = m0; end
        end
        slices(end+1) = struct('machine', refMachine, 'plan', 'clinical', 'cubes', E.cubes, 'masks', mk);
        info.runs{end+1} = struct('machine', refMachine, 'plan', 'clinical', 'status', 'ok', ...
            'notes', {{sprintf('RTDOSE %s, %s, factor to total %g', refDose.doseType, ...
            refDose.summationType, refDose.toTotal)}}, 'seconds', 0);
    catch err
        fprintf(2, '[cns] reference dose FAILED: %s\n', err.message);
    end
end

dvh = [dvhAll{:}];
lvh = [lvhAll{:}];
writeMetrics(fullfile(caseDir, 'metrics.csv'), rows);
writeCurves(fullfile(caseDir, 'dvh.csv'), alias, dvh, 'dose', 'doseGyRBE');
writeCurves(fullfile(caseDir, 'lvh.csv'), alias, lvh, 'letd', 'letdKeVum');
info.finished = datestr(now, 31);
cns_writeJson(fullfile(caseDir, 'run_info.json'), info);

try
    cns_plotCase(caseDir, alias, rows, dvh, lvh, slices, caseCfg, study, dataClass);
catch err
    fprintf(2, '[cns] figures failed: %s\n', err.message);
end
fprintf('[cns] case %s done: %d metric rows in %s\n', alias, size(rows, 1), caseDir);
out = struct('caseDir', caseDir, 'rows', {rows}, 'dvh', dvh, 'lvh', lvh, 'info', info);
end

% =========================================================================
function [w, f] = normalise(w, dij, cstD, caseCfg, study)
% Scale the spot weights so that <structure> <metric> (RBE 1.1) = fraction x Rx.
f = 1;
N = cns_getfield(study, 'normalisation', struct());
if isempty(fieldnames(N)) || strcmp(cns_getfield(N, 'structure', 'none'), 'none'), return; end
ix = find(strcmp(cstD(:,2), N.structure), 1);
if isempty(ix), return; end
rx = caseCfg.prescription.totalDoseGyRBE;
d = full(dij.physicalDose{1} * w) * caseCfg.prescription.numOfFractions * 1.1;
x = sort(d(cstD{ix,4}{1}), 'descend');
pct = str2double(regexprep(N.metric, '^D', ''));
cur = x(max(1, min(numel(x), ceil(pct / 100 * numel(x)))));
if cur > 0
    f = N.fractionOfRx * rx / cur;
    w = w * f;
end
end

function s = optStatus(resultGUI)
s = '';
try
    if isfield(resultGUI.info, 'status'), s = num2str(resultGUI.info.status); end
catch
end
end

function S = tagRows(S, machine, plan)
for i = 1:numel(S)
    S(i).machine = machine;
    S(i).plan = plan;
end
end

function saveMat(f, isOctave, varargin)
S = struct();
for i = 1:2:numel(varargin), S.(varargin{i}) = varargin{i+1}; end
if isOctave
    save('-v7', f, '-struct', 'S');
else
    save(f, '-struct', 'S', '-v7.3');
end
end

function writeMetrics(f, rows)
fid = fopen(f, 'w');
fprintf(fid, 'alias,machine,plan,structure,quantity,metric,value,unit\n');
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%.6g,%s\n', rows{i,1:6}, rows{i,7}, rows{i,8});
end
fclose(fid);
end

function writeCurves(f, alias, S, xField, xName)
fid = fopen(f, 'w');
hasCc = ~isempty(S) && isfield(S, 'volCc');
q = ~isempty(S) && isfield(S, 'quantity');
fprintf(fid, 'alias,machine,plan,structure,quantity,%s,volPct%s\n', xName, ifelse(hasCc, ',volCc', ''));
for i = 1:numel(S)
    if q, qq = S(i).quantity; else, qq = 'LETd'; end
    x = S(i).(xField); v = S(i).volPct;
    last = find(v > 0, 1, 'last');
    if isempty(last), last = 1; end
    for k = 1:min(numel(x), last + 1)
        if hasCc
            fprintf(fid, '%s,%s,%s,%s,%s,%.4g,%.5g,%.5g\n', alias, S(i).machine, S(i).plan, ...
                S(i).structure, qq, x(k), v(k), S(i).volCc(k));
        else
            fprintf(fid, '%s,%s,%s,%s,%s,%.4g,%.5g\n', alias, S(i).machine, S(i).plan, ...
                S(i).structure, qq, x(k), v(k));
        end
    end
end
fclose(fid);
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end
