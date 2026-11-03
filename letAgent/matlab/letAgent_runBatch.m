function letAgent_runBatch(batchFile, outFile)
% letAgent_runBatch  Optimise a list of plan recipes and write their metrics.
%
%   letAgent_runBatch('round_03_w1.json', 'round_03_w1_results.json')
%
% batch.caseFile  path to the case JSON
% batch.runs      list of {runId, beamSet, recipe}
%
% Each plan's spot weights are saved next to the cache (plans/<runId>_w.mat)
% so the chosen plan can be recalculated, exported or checked with Monte
% Carlo later. Only metrics go into outFile.

matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;
matRad_cfg.logLevel = 1;

batch = letAgent_readJson(batchFile);
caseCfg = letAgent_readJson(batch.caseFile);
letAgent_checkAlias(caseCfg.alias);
letAgent_checkDataLocation(caseCfg, batch.caseFile);
caseDir = letAgent_caseDir(caseCfg);
planDir = fullfile(caseDir, 'plans');
if ~exist(planDir, 'dir'), mkdir(planDir); end

maxIter = letAgent_getfield(caseCfg, 'engine.maxIter', matRad_cfg.defaults.propOpt.maxIter);
matRad_cfg.defaults.propOpt.maxIter = maxIter;

runs = letAgent_toCell(batch.runs);
results = cell(1, numel(runs));
beamSetIds = cellfun(@(r) r.beamSet, runs, 'UniformOutput', false);
[uIds, ~, grp] = unique(beamSetIds);

for g = 1:numel(uIds)
    S = load(fullfile(caseDir, ['cache_' uIds{g} '.mat']), 'ct', 'cst', 'pln', 'dij');
    S.pln.multScen = letAgent_buildMultScen(S.ct, caseCfg);
    for k = find(grp(:)' == g)
        r = runs{k};
        res = struct('runId', r.runId, 'beamSet', r.beamSet, 'recipe', r.recipe);
        tStart = tic;
        try
            [cstRun, notes] = letAgent_buildCst(S.cst, caseCfg, r.recipe, S.pln);
            plnRun = S.pln;
            plnRun.propOpt.optimizer = 'IPOPT';
            resultGUI = matRad_fluenceOptimization(S.dij, cstRun, plnRun);
            w = resultGUI.wUnsequenced;
            res.status = 'ok';
            res.notes = notes;
            res.optimizer = optimizerInfo(resultGUI);
            res.metrics = letAgent_metrics(S.dij, S.cst, w, caseCfg, S.pln);
            save(fullfile(planDir, [r.runId '_w.mat']), 'w');
        catch err
            res.status = 'failed';
            res.error = err.message;
        end
        res.seconds = toc(tStart);
        results{k} = res;
        fprintf('[letAgent] %s %s (%.0f s)\n', r.runId, res.status, res.seconds);
        letAgent_writeJson(outFile, struct('runs', {results(~cellfun(@isempty, results))}));
    end
    clear S
end
end

function s = optimizerInfo(resultGUI)
s = struct('iterations', NaN, 'status', '');
try
    info = resultGUI.info;
    if isfield(info, 'iter'), s.iterations = info.iter; end
    if isfield(info, 'status'), s.status = num2str(info.status); end
catch
end
end
