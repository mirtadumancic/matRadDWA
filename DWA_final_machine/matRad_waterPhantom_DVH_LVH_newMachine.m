%% matRad_waterPhantom_DVH_LVH_newMachine.m
% DVH and LVH of the robust DWA water-phantom plan (centre target) with the NEW
% machine DWA_proton_v4 (v3 dose tables + MCsquare LET table), for the 4-panel
% comparison with RayStation ("DWA 1 mm matRad NEW").
%
% Plan: robust dose-only plan from matRad_waterPhantom_LETopt_v5_robust.m
% (CTV-based COWC, 9 scenarios, worst-scenario CTV D95 = 95 % Rx, 5 mm spots).
% v3 and v4 have identical dose tables, so the dose is that plan's dose; the
% LETd is recalculated with the v4 LET table, pencil beam, nominal scenario,
% LETd in every voxel with dose (no 1 % cut-off, as RayStation).
%
% Output (no cubes, histograms only): DWA_final_machine/WaterPhantom_newMachine/
%   DVH_LVH_DWA_proton_v4_center.csv  (region, quantity, x, volume_pct)
%   metrics_DWA_proton_v4_center.csv
% Run on the PC (where the v5 plan is), from X:\matRadDWA_LET:
%   matRad_rc; run('DWA_final_machine/matRad_waterPhantom_DVH_LVH_newMachine.m')

clear; clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir); matRad_rc;
matRad_cfg = MatRad_Config.instance(); matRad_cfg.disableGUI = true;

planFile = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v5_robust', 'plan_DWA_proton_v3_center_LETw0_bw5.mat');
newMachine = 'DWA_proton_v4';
nFx = 30; res_mm = 3;
outDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_newMachine');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% phantom (identical to the planning scripts)
semiLR = 250; semiAP = 150; halfSI = 150; airMargin = 30;
Ny = ceil((2*semiAP + 2*airMargin)/res_mm); Nx = ceil((2*semiLR + 2*airMargin)/res_mm);
Nz = ceil((2*halfSI + 2*airMargin)/res_mm);
ct = struct(); ct.resolution = struct('x', res_mm, 'y', res_mm, 'z', res_mm);
ct.cubeDim = [Ny Nx Nz]; ct.numOfCtScen = 1; ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
c0 = [mean(ct.x) mean(ct.y) mean(ct.z)];
[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
body = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 <= 1 & abs(Zg - c0(3)) <= halfSI;
clear Yg Xg Zg
rED = 0.0012*ones(ct.cubeDim); rED(body) = 1; HU = -1000*ones(ct.cubeDim); HU(body) = 0;
ct.cube = {rED}; ct.cubeHU = {HU}; ct.hlut = [-1024 0.00; 0 1.00; 3000 2.5];

%% plan
S = load(planFile);
tMask = false(ct.cubeDim); tMask(double(S.masks.tMask)) = true;
outside = false(ct.cubeDim); outside(double(S.masks.outsideEll)) = true;
cst = cell(2, 6);
cst(1,:) = {0, 'TARGET', 'TARGET', {find(tMask)},         struct('Priority', 1, 'alphaX', 0.1, 'betaX', 0.01, 'Visible', 1), {}};
cst(2,:) = {1, 'BODY',   'OAR',    {find(body & ~tMask)}, struct('Priority', 2, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1), {}};
pln = S.pln; pln.machine = newMachine; pln.multScen = 'nomScen';
pln.propDoseCalc.calcLET = 1;
stf = S.stf; for i = 1:numel(stf), stf(i).machine = newMachine; end

matRad_cfg.defaults.doseWeightingThreshold = 0;        % LETd wherever dose > 0
res = matRad_calcDoseForward(ct, cst, stf, pln, S.w(:));
matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'doseWeightingThreshold');
cd(rootDir);
if isfield(res, 'RBExD'), D = res.RBExD; elseif isfield(res, 'RBExDose'), D = res.RBExDose; else, D = 1.1*res.physicalDose; end
D = nFx * S.sc * D;
L = res.LET; L(~isfinite(L) | D <= 0) = 0;
fprintf('dose check vs saved plan dose: max |diff| in CTV %.3f Gy\n', max(abs(D(tMask) - double(S.Dtot(tMask)))));

%% histograms (cumulative, % of ROI volume)
gD = (0:0.02:65)'; gL = (0:0.01:35)';
H = table();
regs = {'CTV', tMask; 'PhantomMinusCTV', outside};
for r = 1:2
    d = D(regs{r,2}); l = L(regs{r,2}); n = numel(d);
    vD = 100 * arrayfun(@(x) nnz(d >= x), gD) / n;
    vL = 100 * arrayfun(@(x) nnz(l >= x), gL) / n;
    H = [H; table(repmat(regs(r,1), numel(gD), 1), repmat({'dose_Gy'}, numel(gD), 1), gD, vD, ...
                  'VariableNames', {'region', 'quantity', 'x', 'volume_pct'})]; %#ok<AGROW>
    H = [H; table(repmat(regs(r,1), numel(gL), 1), repmat({'LETd_keVum'}, numel(gL), 1), gL, vL, ...
                  'VariableNames', {'region', 'quantity', 'x', 'volume_pct'})]; %#ok<AGROW>
end
writetable(H, fullfile(outDir, 'DVH_LVH_DWA_proton_v4_center.csv'));
v = sort(D(tMask)); q = @(p) v(max(1, ceil(p*numel(v))));
M = table(q(0.02), q(0.05), q(0.5), q(0.98), mean(v), mean(L(tMask)), ...
    mean(D(tMask)/1.1)/mean(D(outside)/1.1), mean(L(tMask))/mean(L(outside)), 100*mean(L(outside) > 0), ...
    'VariableNames', {'CTV_D98', 'CTV_D95', 'CTV_D50', 'CTV_D2', 'CTV_Dmean', 'CTV_LETd_mean', ...
                      'DinDout', 'LETratio', 'pctOutside_LETpos'});
writetable(M, fullfile(outDir, 'metrics_DWA_proton_v4_center.csv')); disp(M);
fprintf('Saved to %s\n', outDir);
