%% run_UHN_CNS_09_CNS02_v2_airOutside.m
% As run_UHN_CNS_09_CNS02_v2.m, but with everything outside the Body contour set to air in all
% matRad calculations (test of whether RayStation ignored the mask / head rest / couch).
% Results go to results_RayStation_v2_airOutside, next to results_RayStation_v2 (material included).
%   clear all; matRad_rc; run('CNS_ESTRO/run_UHN_CNS_09_CNS02_v2_airOutside.m')
S_override = struct( ...
    'caseFile', 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/CNS_02.json', ...
    'rsDir',    'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/RayStation_DWA', ...  % <-- same folder as before
    'runMC',    true, 'nHistMC', 2e7, 'nHistMC_RSplan', 5e7, ...
    'exactEnergies', true, 'fitMU', true, 'nBins', 8, ...
    'outsideBodyAir', true, 'outTag', '_airOutside');
if ~exist(S_override.rsDir, 'dir'), error('Put the RayStation RP/RD files in %s (or edit rsDir).', S_override.rsDir); end
here = fileparts(mfilename('fullpath'));
run(fullfile(here, 'UHN_CNS_09_RayStation_vs_matRad_v2.m'));
