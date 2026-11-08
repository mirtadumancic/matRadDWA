%% run_UHN_CNS_07_smoke_v1.m
% Quick end-to-end test of UHN_CNS_07 on ONE patient before the overnight run
% (about 5-10 min): robust 9-scenario dij, IPOPT, clinical goals (PB nominal and
% worst), photon reference, MCsquare, tables and slice figures, all with coarse
% settings (study spot spacing, 80 iterations, 1e5 MC protons). Written to a separate
% <alias>_DWAv4rob_SMOKE folder, so the caches of the real run are not touched.
% The numbers are NOT results: this only shows that every step runs.
%
% RUN ON THE UHN PC:  clear all; cd X:\matRadDWA_LET; run('CNS_ESTRO/run_UHN_CNS_07_smoke_v1.m')

dataRoot = 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data';
pat = struct('caseFile', fullfile(dataRoot, 'Patient1', 'CNS_01.json'), ...
    'gantryAngles', [270 100 75], 'couchAngles', [0 0 300], 'label', 'CNS_01');   % Patient 1

here = fileparts(mfilename('fullpath'));
script = fullfile(here, 'UHN_CNS_07_DWAv4_robustDoseOnly_clinicalGoals_v1.m');
S_override = struct('caseFile', pat.caseFile, 'gantryAngles', pat.gantryAngles, ...
    'couchAngles', pat.couchAngles, 'setLabel', ['RS_approved_' pat.label], ...
    'geometryOnly', false, 'runMC', true, 'smoke', true, 'forceDij', true);
tSmoke = tic;
try
    run(script);
    ok = true; msg = '';
catch err
    ok = false;
    msg = sprintf('%s (line %d of %s)', err.message, err.stack(1).line, err.stack(1).name);
end
fprintf('\n================ SMOKE TEST %s: %s after %.1f min ================\n', ...
    pat.label, ifelse(ok, 'PASSED', 'FAILED'), toc(tSmoke) / 60);
if ~ok, fprintf(2, '%s\n', msg); end
if exist('outDir', 'var') && exist(outDir, 'dir')
    f = dir(fullfile(outDir, '*.*')); f = f(~[f.isdir]);
    fprintf('Files written in %s:\n', outDir); fprintf('  %s\n', f.name);
end
fprintf('Check: robust plan, PB worst-scenario table, MCsquare table and slice figures should all be listed.\n');

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end
