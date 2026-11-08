%% matRad_waterPhantom_MCsquareLET_RSstyle.m
% Recalculates a saved water-phantom plan with MCsquare and scores LETd the
% way RayStation's Monte Carlo does, so that the matRad and RayStation
% LET-volume histograms can be compared like with like.
%
% RayStation MC LETd (Wagenaar et al. 2020, PMB 65 025006, RayStation v6R;
% to be confirmed for RayStation 14 with Alaina / the reference manual):
%   - Monte Carlo transport; no LET input data (Alaina: only spot profiles,
%     Bragg peaks and absolute dose are entered for commissioning)
%   - dose-averaged, unrestricted LET in water (water at the local density)
%   - primary + secondary protons, deuterons and alphas; heavier fragments
%     deposited locally (no LET contribution)
%   - stopping power sampled along each step (stopping-power method)
%   - LETd reported in every voxel that receives dose (no dose cut-off;
%     ~23 % of Phantom - CTV has LETd in Alaina's export)
%
% matRad settings used here to match:
%   - MCsquare transport (protons, deuterons, alphas transported; heavier
%     fragments local) - same particle set as RayStation
%   - LET_Calculation_Method = 'StopPow' (stopping-power method, RayStation-
%     like); 'DepositedEnergy' is also run as a sensitivity check
%   - doseWeightingThreshold = 0: LETd in every voxel with dose > 0
%     (matRad default is 1 % of the maximum dose; needs the small local change
%     in matRad/util/matRad_calcCubes.m from 8 Nov 2026)
% Four LETd versions of the same plan are compared, to separate the two steps:
%   PB_table_cut1pct    pencil beam, machine LET table, matRad default 1 % cut (as planned)
%   PB_table_noCut      step 2 only: same, no dose cut-off
%   MC_StopPow          steps 1 + 2: MCsquare, RayStation-like scoring, no cut-off
%   MC_DepositedEnergy  sensitivity to the LET scoring method
%
% Run on the Windows PC (MCsquare) from X:\matRadDWA_LET:
%   matRad_rc; run('DWA_final_machine/matRad_waterPhantom_MCsquareLET_RSstyle.m')
% Output: DWA_final_machine/WaterPhantom_MCLET_RSstyle/
%   <plan>_MCLET.mat  (Dtot, LET per method, LET_PB, masks; single precision)
%   <plan>_LVH.csv    (cumulative LVHs, CTV and Phantom - CTV, several dose cut-offs)
%   <plan>_summary.csv

clear; clc;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;
if any(rootDir == ' '), error('MCsquare fails on paths with spaces: run from X:\\matRadDWA_LET (subst).'); end

%% ---------------- settings ------------------------------------------------
planFile   = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v2', ...
                      'plan_DWA_proton_v2_center_LETw0_cov.mat');   % dose-only plan, centre target
nHist      = 1e8;            % LETd in low-dose voxels needs many histories; test with 1e7 first
letMethods = {'StopPow', 'DepositedEnergy'};   % first = RayStation-like
doseCut    = [0 0.001 0.01 0.05 0.10];         % LVH dose cut-offs, fraction of Rx
Rx_total   = 54; nFx = 30;
targetDiam_mm = 50; marginLat_mm = 3; rangeUnc = 0.035;   % as in the planning scripts
res_mm     = 3;
letGrid    = (0:0.05:40)';                     % keV/um
% The plan was optimised with DWA_proton_v2, whose depth-dose tables are relative
% (max = 100), so its spot weights are not proton numbers. MCsquare reads weights as
% protons, so the energy layers get the wrong relative intensities (about 17 % tilt
% across this target). Rescaling each weight by Edep_v2 / Edep_v3 of its energy layer
% gives the proton numbers that reproduce the same pencil-beam dose with the
% calibrated DWA_proton_v3.
rescaleToCalibrated = true;
calibratedMachine   = 'DWA_proton_v3';

outDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_MCLET_RSstyle');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, planName] = fileparts(planFile);
if rescaleToCalibrated, outTag = '_wCalibrated'; else, outTag = '_wRaw'; end

%% ---------------- phantom (identical to the planning scripts) -------------
semiLR = 250; semiAP = 150; halfSI = 150; airMargin = 30;
Ny = ceil((2*semiAP + 2*airMargin)/res_mm); Nx = ceil((2*semiLR + 2*airMargin)/res_mm);
Nz = ceil((2*halfSI + 2*airMargin)/res_mm);
ct = struct();
ct.resolution.x = res_mm; ct.resolution.y = res_mm; ct.resolution.z = res_mm;
ct.cubeDim = [Ny Nx Nz]; ct.numOfCtScen = 1; ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
c0 = [mean(ct.x) mean(ct.y) mean(ct.z)];
[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
body = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 <= 1 & abs(Zg - c0(3)) <= halfSI;
ell  = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 + ((Zg - c0(3))/halfSI).^2 <= 1;
rED = 0.0012 * ones(ct.cubeDim); rED(body) = 1;
HU  = -1000  * ones(ct.cubeDim); HU(body)  = 0;
ct.cube = {rED}; ct.cubeHU = {HU};
ct.hlut = [-1024 0.00; 0 1.00; 3000 2.5];

%% ---------------- plan ----------------------------------------------------
S = load(planFile);
w = S.w(:); pln = S.pln; tc = S.tc(:)'; d = S.d(:)'; d = d / norm(d);
sc = 1; if isfield(S, 'sc'), sc = S.sc; end
if isfield(S, 'LETc'), LET_PB = double(S.LETc); elseif isfield(S, 'LET'), LET_PB = double(S.LET); else, LET_PB = []; end
D_PB = double(S.Dtot);

% structures: target used for spot placement = planning target, as in planning
R  = targetDiam_mm/2;
rr = sqrt((Xg - tc(1)).^2 + (Yg - tc(2)).^2 + (Zg - tc(3)).^2);
tMask = rr <= R & body;
ax = (Xg - tc(1))*d(1) + (Yg - tc(2))*d(2) + (Zg - tc(3))*d(3);
pp = sqrt(max(rr.^2 - ax.^2, 0));
mAx = marginLat_mm + rangeUnc * norm(tc - entryPoint(ct, body, tc, d));
pMask = ((pp/(R + marginLat_mm)).^2 + (ax/(R + mAx)).^2 <= 1 & body) | tMask;
if isfield(S, 'masks') && isstruct(S.masks) && isfield(S.masks, 'pMask')   % use the saved planning target if present
    pm = S.masks.pMask; pMask = false(ct.cubeDim);
    if islogical(pm) && isequal(size(pm), ct.cubeDim), pMask = pm; else, pMask(pm(:)) = true; end
end
outside = ell & ~tMask;                         % = RayStation "Ext-CTVcen5"
clear Yg Xg Zg rr ax pp
cst = cell(2, 6);
cst(1,:) = {0, 'TARGET', 'TARGET', {find(pMask)},          struct('Priority', 1, 'alphaX', 0.1, 'betaX', 0.01, 'Visible', 1), {}};
cst(2,:) = {1, 'BODY',   'OAR',    {find(body & ~pMask)},  struct('Priority', 2, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1), {}};
stf = matRad_generateStf(ct, cst, pln);
for i = 1:numel(stf), stf(i).machine = pln.machine; end
if sum([stf.totalNumOfBixels]) ~= numel(w)
    error('Spot grid does not match the plan (%d spots vs %d weights).', sum([stf.totalNumOfBixels]), numel(w));
end
fprintf('%s: %d spots, %s\n', planName, numel(w), pln.machine);

% spot weights -> proton numbers (see settings)
wMC = w;
if rescaleToCalibrated
    mOld = matRad_loadMachine(struct('radiationMode', 'protons', 'machine', pln.machine));
    mNew = matRad_loadMachine(struct('radiationMode', 'protons', 'machine', calibratedMachine));
    Eold = [mOld.data.energy]; Enew = [mNew.data.energy];
    if numel(Eold) ~= numel(Enew) || max(abs(Eold - Enew)) > 1e-6, error('Energy layers of %s and %s differ.', pln.machine, calibratedMachine); end
    cE = arrayfun(@(i) trapz(mNew.data(i).depths, mNew.data(i).Z) / trapz(mOld.data(i).depths, mOld.data(i).Z), 1:numel(Eold));
    eSpot = zeros(numel(w), 1); n = 0;
    for i = 1:numel(stf), for j = 1:stf(i).numOfRays, for k = 1:stf(i).numOfBixelsPerRay(j)
        n = n + 1; eSpot(n) = stf(i).ray(j).energy(k);
    end, end, end
    [~, ixE] = min(abs(eSpot - Eold), [], 2);
    wMC = w ./ cE(ixE)';
    fprintf('  weights rescaled to %s: calibration factor %.3f (shallowest layer) to %.3f (deepest layer)\n', ...
        calibratedMachine, max(cE(unique(ixE))), min(cE(unique(ixE))));
end

%% ---------------- MCsquare, RayStation-like LET ---------------------------
matRad_cfg.defaults.doseWeightingThreshold = 0;      % LETd wherever dose > 0
plnMC = pln;
plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', 1, 'numHistoriesDirect', nHist, ...
                            'doseGrid', pln.propDoseCalc.doseGrid);
LETmc = struct(); Dmc = [];
for m = 1:numel(letMethods)
    engine = DoseEngines.matRad_ParticleMCsquareEngine(plnMC);
    engine.config.LET_Calculation_Method = letMethods{m};
    t0 = tic;
    res = engine.calcDoseForward(ct, cst, stf, wMC);
    cd(rootDir);
    if isfield(res, 'RBExD'), D = res.RBExD; elseif isfield(res, 'RBExDose'), D = res.RBExDose; else, D = 1.1 * res.physicalDose; end
    D = nFx * sc * D;
    L = res.LET; L(~isfinite(L) | D <= 0) = 0;
    if isempty(Dmc), Dmc = D; end
    LETmc.(letMethods{m}) = L;
    fprintf('  %-16s %.0f s: dose>0 in %.1f %% of grid, LET>0 in %.1f %%; Phantom-CTV: dose>0 %.1f %%, LET>0 %.1f %%\n', ...
        letMethods{m}, toc(t0), 100*mean(D(:) > 0), 100*mean(L(:) > 0), 100*mean(D(outside) > 0), 100*mean(L(outside) > 0));
end
% pencil beam (machine LET table) recalculated with the same zero threshold
plnPB = pln;
plnPB.propDoseCalc = struct('engine', 'HongPB', 'lateralModel', 'double', 'calcLET', 1, ...
                            'doseGrid', pln.propDoseCalc.doseGrid);
res = matRad_calcDoseForward(ct, cst, stf, plnPB, w);
cd(rootDir);
if isfield(res, 'RBExD'), D = res.RBExD; elseif isfield(res, 'RBExDose'), D = res.RBExDose; else, D = 1.1 * res.physicalDose; end
D_PB0 = nFx * sc * D;
LET_PB0 = res.LET; LET_PB0(~isfinite(LET_PB0) | D_PB0 <= 0) = 0;
fprintf('  %-16s Phantom-CTV: dose>0 %.1f %%, LET>0 %.1f %%\n', 'HongPB, no cut', ...
    100*mean(D_PB0(outside) > 0), 100*mean(LET_PB0(outside) > 0));
matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'doseWeightingThreshold');   % back to matRad default

%% ---------------- LVHs and summary ----------------------------------------
regs = {'CTV', tMask; 'PhantomMinusCTV', outside};
srcs = cell(0, 3);
if ~isempty(LET_PB), srcs(end+1,:) = {'PB_table_cut1pct', LET_PB, D_PB}; end   % as planned: matRad default 1 % cut
srcs(end+1,:) = {'PB_table_noCut', LET_PB0, D_PB0};                              % step 2 only
for m = 1:numel(letMethods), srcs(end+1,:) = {['MC_' letMethods{m}], LETmc.(letMethods{m}), Dmc}; end %#ok<SAGROW>
H = table(); Ssum = table();
for s = 1:size(srcs, 1)
    for r = 1:size(regs, 1)
        for c = doseCut
            m = regs{r,2} & srcs{s,3} > c * Rx_total;
            vRef = nnz(regs{r,2});                     % % of the whole ROI, as in RayStation
            Lm = srcs{s,2}(m);
            vol = 100 * arrayfun(@(x) nnz(Lm >= x), letGrid) / vRef;
            H = [H; table(repmat(srcs(s,1), numel(letGrid), 1), repmat(regs(r,1), numel(letGrid), 1), ...
                repmat(c, numel(letGrid), 1), letGrid, vol, ...
                'VariableNames', {'source', 'region', 'doseCut_fracRx', 'LETd_keVum', 'volume_pctROI'})]; %#ok<AGROW>
            Ssum = [Ssum; table(srcs(s,1), regs(r,1), c, 100*nnz(m)/vRef, 100*nnz(Lm > 0)/vRef, ...
                mean(Lm(Lm > 0)), prctile(Lm(Lm > 0), 99), max([Lm; 0]), ...
                'VariableNames', {'source', 'region', 'doseCut_fracRx', 'pctROI_aboveCut', 'pctROI_LETpos', ...
                                  'LETd_mean_LETpos', 'LETd_p99', 'LETd_max'})]; %#ok<AGROW>
        end
    end
end
writetable(H, fullfile(outDir, [planName outTag '_LVH.csv']));
writetable(Ssum, fullfile(outDir, [planName outTag '_summary.csv']));
disp(Ssum(Ssum.doseCut_fracRx == 0 | Ssum.doseCut_fracRx == 0.01, :));

Dtot = single(Dmc); LET_MC = structfun(@single, LETmc, 'UniformOutput', false);
LET_PB = single(LET_PB); D_PB = single(D_PB); LET_PB0 = single(LET_PB0); D_PB0 = single(D_PB0); %#ok<NASGU>
masks = struct('tMask', find(tMask), 'outside', find(outside)); %#ok<NASGU>
save(fullfile(outDir, [planName outTag '_MCLET.mat']), 'Dtot', 'LET_MC', 'LET_PB', 'D_PB', 'LET_PB0', 'D_PB0', 'masks', 'tc', 'd', 'nHist', 'letMethods', 'wMC', '-v7.3');
fprintf('\nSaved to %s\n', outDir);

%% ---------------- helpers -------------------------------------------------
function e = entryPoint(ct, body, tc, d)
    F = griddedInterpolant({ct.y, ct.x, ct.z}, double(body), 'nearest', 'nearest');
    s = 0;
    while F(tc(2) - s*d(2), tc(1) - s*d(1), tc(3) - s*d(3)) > 0.5 && s < 1000
        s = s + 0.5;
    end
    e = tc - s * d;
end
