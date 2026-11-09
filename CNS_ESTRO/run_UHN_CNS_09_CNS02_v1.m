%% run_UHN_CNS_09_CNS02_v1.m
% RayStation DWA plan vs matRad for CNS_02 (patient 2). Runs UHN_CNS_09 with the settings below.
% 1) Copy Alaina's RayStation folder (RP*.dcm + 4 RD*.dcm) to rsDir on the PC (hospital storage).
% 2) clear all; matRad_rc; run('CNS_ESTRO/run_UHN_CNS_09_CNS02_v1.m')
% Takes ~1-2 h (two MCsquare runs; the RayStation plan has ~110 000 spots).
S_override = struct( ...
    'caseFile', 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/CNS_02.json', ...
    'rsDir',    'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/RayStation_DWA', ...  % <-- check
    'runMC',    true, ...
    'nHistMC',  2e7, ...           % matRad robust plan
    'nHistMC_RSplan', 5e7);        % RayStation plan recalculated in matRad
if ~exist(S_override.rsDir, 'dir'), error('Put the RayStation RP/RD files in %s (or edit rsDir).', S_override.rsDir); end
here = fileparts(mfilename('fullpath'));
run(fullfile(here, 'UHN_CNS_09_RayStation_vs_matRad_v1.m'));
