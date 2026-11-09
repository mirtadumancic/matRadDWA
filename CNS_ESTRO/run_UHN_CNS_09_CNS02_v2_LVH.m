%% run_UHN_CNS_09_CNS02_v2_LVH.m
% LETd-volume histograms for CNS_02 with MCsquare (no dose threshold): (a) the matRad robust plan and
% (b) the RayStation plan recalculated in matRad with the MU calibration already fitted in
% results_RayStation_v2 (reused, so the 24 pencil-beam runs are skipped). About 15-20 min.
%   clear all; matRad_rc; run('CNS_ESTRO/run_UHN_CNS_09_CNS02_v2_LVH.m')
% Send back from results_RayStation_v2_LVH: LVH_MCsquare_CNS_02_v2.csv, LETd_metrics_MCsquare_CNS_02_v2.csv
S_override = struct( ...
    'caseFile', 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/CNS_02.json', ...
    'rsDir',    'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/RayStation_DWA', ...  % <-- same folder as before
    'runMC',    true, 'nHistMC', 2e7, 'nHistMC_RSplan', 5e7, ...
    'exactEnergies', true, 'fitMU', true, 'reuseFit', true, 'skipUnfitMC', true, ...
    'calcLVH', true, 'outTag', '_LVH');
if ~exist(S_override.rsDir, 'dir'), error('Put the RayStation RP/RD files in %s (or edit rsDir).', S_override.rsDir); end
here = fileparts(mfilename('fullpath'));
run(fullfile(here, 'UHN_CNS_09_RayStation_vs_matRad_v2.m'));
