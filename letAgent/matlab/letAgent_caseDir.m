function d = letAgent_caseDir(caseCfg)
% letAgent_caseDir  <workDir>/<alias>; workDir defaults to ./letAgent_runs
workDir = letAgent_getfield(caseCfg, 'workDir', fullfile(pwd, 'letAgent_runs'));
d = fullfile(workDir, caseCfg.alias);
end
