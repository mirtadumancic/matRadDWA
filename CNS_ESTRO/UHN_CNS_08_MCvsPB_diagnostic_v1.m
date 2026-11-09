%% UHN_CNS_08_MCvsPB_diagnostic_v1.m
% Why is MCsquare ~7-10 % below the pencil beam (PB) in the CTV for CNS_01 (nominal smoke
% plan, DWA_proton_v4, 3 mm spots), when the two agree within 1 % in the water phantom?
% Same plan (stf + weights of the SMOKEnom run), forward dose with both engines in 4 set-ups:
%   A  patient CT, PB ignores densities outside the contours (matRad default) | MC sees all
%   B  patient CT, PB includes everything outside the contours (mask, head rest, couch)
%   C  patient CT with everything outside the Body set to air in BOTH engines
%   D  water patient: HU 0 inside the Body, air outside, in BOTH engines
% Reading: if MC/PB -> 1 in D but not A-C, it is the CT/heterogeneity handling;
%          if MC/PB stays ~0.93 even in D, it is the beam model / spot handling in MCsquare.
% Output (numbers only, no images): <SMOKEnom case>/diagnostic_MCvsPB_v1.csv
% RUN ON THE UHN PC after the nominal smoke test:
%   clear all; cd X:\matRadDWA_LET; matRad_rc; run('CNS_ESTRO/UHN_CNS_08_MCvsPB_diagnostic_v1.m')

caseFile = 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient1/CNS_01.json';
nHist = 2e6;
mcScanner = 'matRad_hlutConsistent';

matRadRoot = fileparts(fileparts(mfilename('fullpath')));
cd(matRadRoot); matRad_rc; addpath(fullfile(matRadRoot, 'CNS_ESTRO'));
matRad_cfg = MatRad_Config.instance(); matRad_cfg.disableGUI = true;
caseCfg0 = cns_readJson(caseFile); cns_checkPrivacy(caseCfg0, caseFile);
caseDir = fullfile(cns_expandPath(caseCfg0.workDir), [caseCfg0.alias '_DWAv4rob_SMOKEnom']);
P = load(fullfile(caseDir, 'prepared.mat'), 'ct', 'cst');
D = load(fullfile(caseDir, 'dij_DWA_proton_v4.mat'), 'pln', 'stf');
W = load(fullfile(caseDir, 'plans', 'DWA_proton_v4__robustDoseOnly_w.mat'), 'w');
ct0 = P.ct; cst = P.cst; stf = D.stf; w = W.w;
for i = 1:numel(stf), stf(i).machine = 'DWA_proton_v4'; end
mask = @(n) cst{find(strcmp(cst(:,2), n), 1), 4}{1};
body = false(ct0.cubeDim); body(mask('Body')) = true;

% set-up CTs
ctAir = ct0; ctAir.cubeHU{1}(~body) = -1000;
ctWat = ct0; ctWat.cubeHU{1}(body) = 0; ctWat.cubeHU{1}(~body) = -1000;
% (matRad rebuilds the stopping-power cube from cubeHU in every engine: matRad_calcWaterEqD)
setups = {'A_patient_PBignoresOutside', ct0, true; 'B_patient_PBincludesOutside', ct0, false; ...
          'C_outsideBodyAir', ctAir, true; 'D_waterPatient', ctWat, true};

rows = {};
for k = 1:size(setups, 1)
    [name, ct, ign] = setups{k, :};
    plnPB = D.pln; plnPB.multScen = 'nomScen';
    plnPB.propDoseCalc.calcLET = 0; plnPB.propDoseCalc.ignoreOutsideDensities = ign;
    rPB = matRad_calcDoseForward(ct, cst, stf, plnPB, w);
    plnMC = D.pln; plnMC.multScen = 'nomScen';
    plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', 0, 'numHistoriesDirect', nHist, ...
        'doseGrid', D.pln.propDoseCalc.doseGrid);
    matRad_cfg.defaults.MCsquareScanner = mcScanner;
    rMC = matRad_calcDoseForward(ct, cst, stf, plnMC, w);
    matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'MCsquareScanner');
    cd(matRadRoot);
    for sName = {'CTV', 'GTV', 'PTV', 'Body'}
        ix = find(strcmp(cst(:,2), sName{1}), 1); if isempty(ix), continue; end
        v = cst{ix,4}{1};
        pb = mean(rPB.physicalDose(v)); m = mean(rMC.physicalDose(v));
        rows(end+1, :) = {name, sName{1}, pb, m, m / pb}; %#ok<SAGROW>
        fprintf('%-28s %-5s  PB %.4g  MC %.4g  MC/PB %.3f\n', name, sName{1}, pb, m, m / pb);
    end
end
T = cell2table(rows, 'VariableNames', {'setup', 'structure', 'meanPhysDose_PB', 'meanPhysDose_MC', 'MC_over_PB'});
writetable(T, fullfile(caseDir, 'diagnostic_MCvsPB_v1.csv'));
disp(T); fprintf('Saved %s\n', fullfile(caseDir, 'diagnostic_MCvsPB_v1.csv'));
