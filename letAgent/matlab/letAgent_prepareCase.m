function summary = letAgent_prepareCase(caseFile)
% letAgent_prepareCase  Precompute everything the search needs for one case.
%
%   summary = letAgent_prepareCase('cases/TG119.json')
%
% For every beam set in the case file this builds the steering geometry and
% the dose, LET x dose and dirty dose influence matrices (for all robustness
% scenarios), and caches them in <workDir>/<alias>/cache_<beamSetId>.mat.
% This is the expensive step; every plan in the search reuses the cache.
%
% Everything stays on this computer. The only file meant to be read by the
% agent is case_summary.json, which holds the alias, beam-set ids, bixel and
% scenario counts and the structure names from the case file. No image data,
% no patient identifiers.

matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;

caseCfg = letAgent_readJson(caseFile);
letAgent_checkAlias(caseCfg.alias);
letAgent_checkDataLocation(caseCfg, caseFile);
outDir = letAgent_caseDir(caseCfg);
if ~exist(outDir, 'dir'), mkdir(outDir); end

[ct, cst] = letAgent_loadCase(caseCfg);
cst = letAgent_checkStructures(cst, caseCfg);

beamSets = letAgent_toCell(caseCfg.beamSets);
summary = struct('alias', caseCfg.alias, 'beamSets', {{}}, 'numOfScenarios', [], ...
    'structures', {letAgent_structureNames(caseCfg)}, 'dirtyDoseLETThreshold', ...
    letAgent_getfield(caseCfg, 'metrics.dirtyDoseLETThreshold', 2));

for b = 1:numel(beamSets)
    bs = beamSets{b};
    cacheFile = fullfile(outDir, ['cache_' bs.id '.mat']);
    pln = letAgent_buildPln(ct, cst, caseCfg, bs);

    if exist(cacheFile, 'file') && ~letAgent_getfield(caseCfg, 'forceRecompute', false)
        fprintf('[letAgent] beam set %s: cache exists, skipping (set forceRecompute to redo)\n', bs.id);
        S = load(cacheFile, 'info');
        info = S.info;
    else
        fprintf('[letAgent] beam set %s: dose + LET influence matrices...\n', bs.id);
        tStart = tic;
        stf = matRad_generateStf(ct, cst, pln);
        engine = DoseEngines.matRad_ParticleHongPencilBeamEngine(pln);
        engine.calcLET = true;
        dij = engine.calcDoseInfluence(ct, cst, stf);
        dij = matRad_DirtyDose.addDirtyDoseDij(dij, summary.dirtyDoseLETThreshold);

        info = struct('beamSetId', bs.id, 'numOfBixels', dij.totalNumOfBixels, ...
            'numOfScenarios', pln.multScen.totNumScen, 'dijSeconds', toc(tStart));
        if matRad_cfg.isOctave
            save('-v7', cacheFile, 'ct', 'cst', 'pln', 'stf', 'dij', 'info');
        else
            save(cacheFile, 'ct', 'cst', 'pln', 'stf', 'dij', 'info', '-v7.3');
        end
        fprintf('[letAgent] beam set %s: %d bixels, %d scenarios, %.0f s\n', bs.id, ...
            info.numOfBixels, info.numOfScenarios, info.dijSeconds);
    end
    summary.beamSets{end+1} = struct('id', bs.id, 'gantryAngles', bs.gantryAngles, ...
        'numOfBixels', info.numOfBixels, 'numOfScenarios', info.numOfScenarios);
    summary.numOfScenarios = info.numOfScenarios;
end

letAgent_writeJson(fullfile(outDir, 'case_summary.json'), summary);
fprintf('[letAgent] case %s ready in %s\n', caseCfg.alias, outDir);
end

% -------------------------------------------------------------------------
function pln = letAgent_buildPln(ct, cst, caseCfg, bs)
pln = struct();
pln.radiationMode = letAgent_getfield(caseCfg, 'radiationMode', 'protons');
pln.machine = caseCfg.machine;
pln.numOfFractions = caseCfg.prescription.numOfFractions;
pln.propStf.gantryAngles = bs.gantryAngles(:)';
if isfield(bs, 'couchAngles')
    pln.propStf.couchAngles = bs.couchAngles(:)';
else
    pln.propStf.couchAngles = zeros(size(pln.propStf.gantryAngles));
end
pln.propStf.bixelWidth = letAgent_getfield(caseCfg, 'stf.bixelWidth', 5);
pln.propStf.longitudinalSpotSpacing = letAgent_getfield(caseCfg, 'stf.longitudinalSpotSpacing', 3);
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
g = letAgent_getfield(caseCfg, 'doseGridMm', 3);
pln.propDoseCalc.doseGrid.resolution = struct('x', g, 'y', g, 'z', g);
pln.propOpt.runDAO = 0;
pln.propSeq.runSequencing = 0;
pln.propOpt.quantityOpt = 'physicalDose';
pln.bioModel = 'none';

pln.multScen = letAgent_buildMultScen(ct, caseCfg);
end
