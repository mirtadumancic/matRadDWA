%% run_UHN_CNS_09_CNS02_v2.m
% RayStation DWA plan vs matRad for CNS_02 (patient 2), with UHN_CNS_09 v2 (exact RayStation
% energies, fitted MU calibration per energy, comparison at the same CTV D95).
%   clear all; matRad_rc; run('CNS_ESTRO/run_UHN_CNS_09_CNS02_v2.m')
% About 45-60 min: 24 pencil-beam runs (3 beams x 8 energy bins) + 3 MCsquare runs.
S_override = struct( ...
    'caseFile', 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/CNS_02.json', ...
    'rsDir',    'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient2/RayStation_DWA', ...  % <-- same folder as in the v1 driver
    'runMC',    true, ...
    'nHistMC',  2e7, ...           % matRad robust plan
    'nHistMC_RSplan', 5e7, ...     % RayStation plan recalculated in matRad
    'exactEnergies', true, 'fitMU', true, 'nBins', 8);
if ~exist(S_override.rsDir, 'dir'), error('Put the RayStation RP/RD files in %s (or edit rsDir).', S_override.rsDir); end
here = fileparts(mfilename('fullpath'));
run(fullfile(here, 'UHN_CNS_09_RayStation_vs_matRad_v2.m'));
