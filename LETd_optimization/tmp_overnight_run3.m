%% matRad_waterPhantom_LETopt_forMCSquare.m
% = matRad_waterPhantom_LETopt_v2.m plus a second MCsquare result row,
% 'MCrecalc_norm' (MC dose rescaled to D95 = 51.3 Gy), added Oct 28 after the
% first MCsquare test on the PC. Results still go to WaterPhantom_LETopt_v2
% so they sit next to the earlier _fine / _cov runs.
%
% Water-phantom validation of LETd-based planning (Remo's let_opt branch) and
% machine comparison. Oct 28, 2026.
%
% VERSION 2 (Oct 28, after the first DWA center run of v1):
%  - new setting targetUnderdoseWeight (v1 fixed it at 10): weight of the
%    SquaredUnderdosing(0.95 Rx) objective on the planning target. Raising it
%    (e.g. 1000) stops the LET objective from eroding the distal margin.
%  - new metrics PTV_V95 and PTV_D98 (planning target = 5 cm target + margin).
%  - each run writes metrics_<machine>_<loc><runTag>.csv; runTag also goes
%    into the plan file names, so runs never overwrite each other.
%  - output folder WaterPhantom_LETopt_v2 (v1 results stay in WaterPhantom_LETopt).
%
% PHANTOM AND PLANNING: same as matRad_waterPhantom_Bui2026.m (oval water
% phantom 50 x 30 x 30 cm, one 5 cm spherical target per plan, one field,
% 54 Gy(RBE) / 30 fx, RBE 1.1, planning margin 3 mm + 3.5 % of depth as a
% stand-in for robust optimisation, every plan scaled to D95 = 51.3 Gy).
%
% WHAT IS NEW
%  1. LETd optimisation (let_opt branch). Remo's 'LETd' optimisation quantity
%     is LET x dose (dij.mLETDose * w, keV/um x Gy). For every LET weight in
%     letWeights, the dose objectives stay the same and a
%     SquaredOverdosing(weight, 0) objective on LET x dose is added to the
%     tissue DISTAL to the planning target (MARGIN_DIST and SHELL_DIST, the
%     downstream half of the 0-15 mm shell). Weight 0 = the reference plan.
%     The dose-influence matrix is calculated once per machine x target and
%     reused for all weights.
%  2. Optional MCsquare recalculation of the final dose (mcRecalcFor), with
%     the same spot weights. Two rows: 'MCrecalc' (same scaling as the PB plan;
%     absolute dose is NOT comparable because the commissioned machines have
%     relative depth-dose data) and 'MCrecalc_norm' (MC dose rescaled to
%     D95 = 51.3 Gy: the shape comparison to use).
%  3. Optional MCsquare optimisation (mcOpt): dij from MCsquare instead of
%     HongPB (slow; run overnight).
%  4. More metrics (see header of evalPlan at the end): target D95/D98/D2,
%     HI, LETd in the target, D_in/D_out, LET_in/LET_out, lateral penumbra,
%     distal fall-off, integral dose, V20, dose-LET hotspot volumes, LETd and
%     McNamara RBE-weighted dose in the distal 15 mm.
%
% RUN ORDER SUGGESTION (edit the settings block):
%   a) DWA, center only, letWeights [0 1 10 100], mcRecalcFor = []  (first check)
%   b) DWA, posterior and oblique
%   c) IBA_cyclotron and Generic, all locations, letWeights 0 (machine comparison)
%   d) mcRecalcFor = [0 10] for DWA center/posterior (PB vs MC validation)
%
% Works on the let_opt branch. On the old branch it runs LET weight 0 only.
%
% OUTPUT: DWA_final_machine/WaterPhantom_LETopt_v2/
%   metrics_<machine>_<loc><runTag>.csv (one row per plan and dose engine;
%   rewritten after every plan, so partial runs are usable; separate runs
%   never overwrite each other), plan_<machine>_<loc>_LETw<w>.mat
%   (dose and LETd cubes, PB and MC, plus the evaluation masks).

clear; clc; close all;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;
set(0, 'DefaultFigureVisible', 'off');

%% ---------------- settings --------------------------------------------
machines   = {'DWA_proton_v2'};                 % add 'IBA_cyclotron', 'Generic'
locations  = {'anterior'};                        % 'anterior','center','posterior','oblique'
letWeights = struct('DWA_proton_v2', [0 0.1 0.3 1], ...   % 0 = dose-only reference
                    'IBA_cyclotron', 0, 'Generic', 0);
mcRecalcFor = [];        % LET weights whose final dose is recalculated with MCsquare, e.g. [0 10]
nHistFinal  = 1e7;       % MCsquare histories for the recalculation
mcOpt       = false;     % true: optimise on an MCsquare dij instead of HongPB (slow)
nHistBeamlet = 1e4;      % MCsquare histories per spot when mcOpt = true

bixelWidth = struct('DWA_proton_v2', 3, 'IBA_cyclotron', 5, 'Generic', 5);   % mm
longSpotSpacing_mm = 3;
res_mm     = 3;
targetDiam_mm = 50;
isoDepth_mm = 150;
marginLat_mm = 3;  rangeUnc = 0.035;     % planning margin
distal_mm  = 15;         % depth of the distal evaluation / LET-objective shell
Rx_total   = 54;  nFx = 30;
abTarget   = 10; abNormal = 2;           % alpha/beta for McNamara evaluation [Gy]
targetUnderdoseWeight = 1000;   % weight of SquaredUnderdosing(0.95 Rx) on the planning target;
                              % raise (e.g. 1000) to stop the LET objective eroding the distal margin

runTag  = '_cov';           % appended to the csv name, e.g. '_MC' for an MC-recalculation run,
                        % so a later run never overwrites an earlier one
outDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v2');
if ~exist(outDir, 'dir'), mkdir(outDir); end
newAPI = exist('matRad_OptimizationQuantity', 'class') == 8;   % let_opt branch?
fprintf('Optimisation-quantity API (let_opt branch): %d\n', newAPI);

%% ---------------- MCsquare needs text metadata in the machine files -------
if ~isempty(mcRecalcFor) || mcOpt
    for m = 1:numel(machines)
        mf = fullfile(rootDir, 'matRad', 'basedata', ['protons_' machines{m} '.mat']);
        if ~exist(mf, 'file'), continue; end
        S = load(mf, 'machine');
        need = {'description', machines{m}; 'created_by', 'Mirta Dumancic'; 'created_on', datestr(now, 'yyyy-mm-dd')};
        changed = false;
        for k = 1:size(need, 1)
            if ~isfield(S.machine.meta, need{k,1}) || isempty(S.machine.meta.(need{k,1}))
                S.machine.meta.(need{k,1}) = need{k,2}; changed = true;
            end
        end
        if changed, machine = S.machine; save(mf, 'machine'); end
    end
    clear S machine
end

%% ---------------- phantom ----------------------------------------------
semiLR = 250; semiAP = 150; halfSI = 150; airMargin = 30;   % mm
Ny = ceil((2*semiAP + 2*airMargin)/res_mm); Nx = ceil((2*semiLR + 2*airMargin)/res_mm);
Nz = ceil((2*halfSI + 2*airMargin)/res_mm);
ct = struct();
ct.resolution.x = res_mm; ct.resolution.y = res_mm; ct.resolution.z = res_mm;
ct.cubeDim = [Ny Nx Nz];
ct.numOfCtScen = 1;
ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
c0 = [mean(ct.x) mean(ct.y) mean(ct.z)];
[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
body = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 <= 1 & abs(Zg - c0(3)) <= halfSI;
ell  = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 + ((Zg - c0(3))/halfSI).^2 <= 1;   % oval phantom
rED = 0.0012 * ones(ct.cubeDim);  rED(body) = 1;
HU  = -1000  * ones(ct.cubeDim);  HU(body)  = 0;
ct.cube = {rED};  ct.cubeHU = {HU};
ct.hlut = [-1024 0.00; 0 1.00; 3000 2.5];
voxVol_cc = (res_mm/10)^3;

dirOf = @(g) beamDirection(ct, body, c0, g);
tgt.anterior  = c0 + (50  - semiAP) * dirOf(0);
tgt.center    = c0;
tgt.posterior = c0 + (250 - semiAP) * dirOf(0);
tgt.oblique   = c0 + [150 0 0];
gantry = struct('anterior', 0, 'center', 0, 'posterior', 0, 'oblique', 330);

%% ---------------- loop ----------------------------------------------------
T = table();
for m = 1:numel(machines)
    mName = machines{m};
    mach = matRad_loadMachine(struct('radiationMode', 'protons', 'machine', mName));
    minRange = min([mach.data.peakPos]);
    if isfield(letWeights, mName), wList = letWeights.(mName); else, wList = 0; end
    if ~newAPI && any(wList > 0)
        warning('Old branch: LET optimisation not available, running weight 0 only.'); wList = 0;
    end

    for L = 1:numel(locations)
        loc = locations{L}; tc = tgt.(loc); g = gantry.(loc);
        d = dirOf(g);
        entry = entryPoint(ct, body, tc, d);
        depthToCentre = norm(tc - entry);
        if minRange > depthToCentre - targetDiam_mm/2 + 2
            fprintf('\n[%s / %s] skipped: needs a range shifter\n', mName, loc); continue;
        end

        % ---- structures (disjoint)
        R  = targetDiam_mm/2;
        rr = sqrt((Xg - tc(1)).^2 + (Yg - tc(2)).^2 + (Zg - tc(3)).^2);
        tMask = rr <= R & body;
        dist  = rr - R;
        ax = (Xg - tc(1))*d(1) + (Yg - tc(2))*d(2) + (Zg - tc(3))*d(3);   % along the beam, + downstream
        pp = sqrt(max(rr.^2 - ax.^2, 0));
        mAx = marginLat_mm + rangeUnc*depthToCentre;
        pMask = ((pp/(R + marginLat_mm)).^2 + (ax/(R + mAx)).^2 <= 1 & body) | tMask;
        mOff = marginLat_mm;
        down = ax > 0;
        margin = body & ~pMask & dist <= mOff + 2;
        shell  = body & ~pMask & dist > mOff + 2 & dist <= mOff + 12;
        far    = body & ~pMask & dist > mOff + 12;
        % evaluation regions
        outsideEll = ell & ~tMask;
        distalEval = body & ~tMask & down & dist <= distal_mm;       % 0-15 mm behind the 5 cm target

        cst = cell(6, 6);
        cst(1,:) = {0, 'TARGET',      'TARGET', {find(pMask)},            tissue(1, abTarget), {}};
        cst(2,:) = {1, 'MARGIN_DIST', 'OAR',    {find(margin & down)},    tissue(2, abNormal), {}};
        cst(3,:) = {2, 'MARGIN_LAT',  'OAR',    {find(margin & ~down)},   tissue(2, abNormal), {}};
        cst(4,:) = {3, 'SHELL_DIST',  'OAR',    {find(shell & down)},     tissue(3, abNormal), {}};
        cst(5,:) = {4, 'SHELL_LAT',   'OAR',    {find(shell & ~down)},    tissue(3, abNormal), {}};
        cst(6,:) = {5, 'FAR',         'OAR',    {find(far)},              tissue(4, abNormal), {}};

        % ---- plan
        pln = struct();
        pln.radiationMode  = 'protons';
        pln.machine        = mName;
        pln.bioModel       = 'constRBE';
        pln.numOfFractions = nFx;
        pln.propStf.gantryAngles = g;
        pln.propStf.couchAngles  = 0;
        pln.propStf.bixelWidth   = bixelWidth.(mName);
        pln.propStf.longitudinalSpotSpacing = longSpotSpacing_mm;
        pln.propStf.isoCenter    = entry + isoDepth_mm * d;
        if mcOpt
            pln.propDoseCalc.engine = 'MCsquare';
            pln.propDoseCalc.numHistoriesPerBeamlet = nHistBeamlet;
        else
            pln.propDoseCalc.engine       = 'HongPB';
            pln.propDoseCalc.lateralModel = 'double';
        end
        pln.propDoseCalc.calcLET = 1;
        pln.propDoseCalc.doseGrid.resolution = struct('x', res_mm, 'y', res_mm, 'z', res_mm);
        if newAPI, pln.propOpt.quantityOpt = 'RBExDose'; else, pln.propOpt.quantityOpt = 'RBExD'; end
        pln.propOpt.runDAO = 0;  pln.propSeq.runSequencing = 0;
        pln.multScen = 'nomScen';

        fprintf('\n===== %s / %s (gantry %d, depth %.0f mm): dij =====\n', mName, loc, g, depthToCentre);
        stf = matRad_generateStf(ct, cst, pln);
        tDij = tic;
        dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
        tDij = toc(tDij);

        Tloc = table();                       % rows of this machine x location only
        for wLET = wList
            cst = setObjectives(cst, Rx_total, wLET, newAPI, targetUnderdoseWeight);
            fprintf('\n----- %s / %s: LET weight %g -----\n', mName, loc, wLET);
            tOpt = tic;
            resultGUI = matRad_fluenceOptimization(dij, cst, pln);
            tOpt = toc(tOpt);
            w = resultGUI.w;

            % ---- pencil-beam (or MC-dij) result, scaled to D95 = 51.3 Gy
            [Dtot, LETc] = totalDose(resultGUI, nFx, ct.cubeDim);
            v = sort(Dtot(tMask));  sc = 0.95*Rx_total / v(max(1, ceil(0.05*numel(v))));
            Dtot = sc * Dtot;
            engineName = 'PB'; if mcOpt, engineName = 'MCopt'; end
            geo = struct('tc', tc, 'd', d, 'ct', ct, 'tMask', tMask, 'outsideEll', outsideEll, ...
                         'distal', distalEval, 'pMask', pMask, 'Rx', Rx_total, 'nFx', nFx, 'voxVol', voxVol_cc, ...
                         'abT', abTarget, 'abN', abNormal);
            r = evalPlan(Dtot, LETc, geo);
            Tloc = addRow(Tloc, mName, loc, g, depthToCentre, wLET, engineName, r, numel(w), sum(w)*sc, tOpt, tDij);
            fprintf('  %s: D_in/D_out %.1f  LET_in/LET_out %.1f  distal LETd %.2f  V(D>=20 & LET>=4) %.1f cc\n', ...
                engineName, r.DinDout, r.LETratio, r.LETd_distal_mean, r.V_D20_LET4_out_cc);
            DtotMC = []; LETMC = [];

            % ---- MCsquare recalculation of the same plan
            if ~mcOpt && any(abs(mcRecalcFor - wLET) < 1e-12)
                plnMC = pln;
                plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', 1, ...
                    'numHistoriesDirect', nHistFinal, 'doseGrid', pln.propDoseCalc.doseGrid);
                fprintf('  MCsquare recalculation, %.0e histories ...\n', nHistFinal);
                tMC = tic;
                rMC = matRad_calcDoseForward(ct, cst, stf, plnMC, w);
                tMC = toc(tMC);
                cd(rootDir);
                [DtotMC, LETMC] = totalDose(rMC, nFx, ct.cubeDim);
                DtotMC = sc * DtotMC;                              % same plan, same scaling
                r2 = evalPlan(DtotMC, LETMC, geo);
                Tloc = addRow(Tloc, mName, loc, g, depthToCentre, wLET, 'MCrecalc', r2, numel(w), sum(w)*sc, tMC, NaN);
                fprintf('  MC: D_in/D_out %.1f  D95 %.1f Gy  distal LETd %.2f\n', r2.DinDout, r2.D95, r2.LETd_distal_mean);
                % The commissioned machines carry relative (not absolute) depth-dose
                % data, so PB and MC differ by a constant factor in absolute dose.
                % Row 'MCrecalc_norm': the MC dose rescaled to D95 = 51.3 Gy like the PB
                % plan (shape comparison). The factor is sumW_scaled(norm)/sumW_scaled(MCrecalc).
                v = sort(DtotMC(tMask));  scMC = 0.95*Rx_total / v(max(1, ceil(0.05*numel(v))));
                r3 = evalPlan(scMC * DtotMC, LETMC, geo);
                Tloc = addRow(Tloc, mName, loc, g, depthToCentre, wLET, 'MCrecalc_norm', r3, numel(w), sum(w)*sc*scMC, tMC, NaN);
                fprintf('  MC (rescaled to D95, factor %.2f): D_in/D_out %.1f  PTV D98 %.1f Gy  distal LETd %.2f\n', ...
                    scMC, r3.DinDout, r3.PTV_D98, r3.LETd_distal_mean);
            end

            % ---- save
            planFile = fullfile(outDir, sprintf('plan_%s_%s_LETw%g%s.mat', mName, loc, wLET, runTag));
            masks = struct('tMask', uint32(find(tMask)), 'distal', uint32(find(distalEval)), ...
                           'outsideEll', uint32(find(outsideEll)), 'pMask', uint32(find(pMask)));
            Dtot = single(Dtot); LETc = single(LETc); DtotMC = single(DtotMC); LETMC = single(LETMC); %#ok<NASGU>
            save(planFile, 'Dtot', 'LETc', 'DtotMC', 'LETMC', 'w', 'sc', 'pln', 'tc', 'd', 'c0', 'masks', 'wLET', '-v7.3');
            writetable(Tloc, fullfile(outDir, sprintf('metrics_%s_%s%s.csv', mName, loc, runTag)));
        end
        T = [T; Tloc];
        clear dij
    end
end
disp(T);
set(0, 'DefaultFigureVisible', 'on');
fprintf('\nSaved to %s\n', outDir);

%% ---------------- local functions ---------------------------------------
function cst = setObjectives(cst, Rx, wLET, newAPI, wUnder)
% dose objectives as in matRad_waterPhantom_Bui2026.m, plus LET x dose on the distal tissue
    oT = {DoseObjectives.matRad_MinDVH(100, Rx, 95), DoseObjectives.matRad_SquaredUnderdosing(wUnder, 0.95*Rx), ...
          DoseObjectives.matRad_MinDVH(5, 52.5, 99), DoseObjectives.matRad_SquaredDeviation(50, Rx), ...
          DoseObjectives.matRad_SquaredOverdosing(100, 57)};
    oM = {DoseObjectives.matRad_SquaredOverdosing(100, 57)};
    oS = {DoseObjectives.matRad_SquaredOverdosing(10, 53), DoseObjectives.matRad_SquaredOverdosing(1, 20)};
    oF = {DoseObjectives.matRad_SquaredOverdosing(10, 27), DoseObjectives.matRad_SquaredOverdosing(1, 20)};
    cst{1,6} = oT;  cst{2,6} = oM;  cst{3,6} = oM;  cst{4,6} = oS;  cst{5,6} = oS;  cst{6,6} = oF;
    if newAPI && wLET > 0
        for row = [2 4]                       % MARGIN_DIST, SHELL_DIST
            o = DoseObjectives.matRad_SquaredOverdosing(wLET, 0);
            o.quantity = 'LETd';              % LET x dose in Remo's branch
            cst{row,6}{end+1} = o;
        end
    end
end

function [Dtot, LET] = totalDose(res, nFx, dims)
    if isfield(res, 'RBExD'), D = res.RBExD;
    elseif isfield(res, 'RBExDose'), D = res.RBExDose;
    else, D = 1.1 * res.physicalDose; end
    if ~isequal(size(D), dims), error('Result not on the CT grid'); end
    Dtot = nFx * D;
    LET = zeros(dims);
    if isfield(res, 'LET'), LET = res.LET; end
    LET(~isfinite(LET) | Dtot <= 0) = 0;
end

function r = evalPlan(Dtot, LET, g)
% Dtot: total dose, Gy(RBE 1.1); LET: LETd, keV/um (0 where no dose)
%  target: D95 D98 D2 Dmean (Gy), V95 (%), HI = (D2-D98)/D50, LETd mean and 2nd percentile
%  D_in/D_out and LET_in/LET_out over the oval phantom (voxels without LET = 0)
%  penumbra 80-20 lateral (cm), distal fall-off 80-20 along the beam (mm)
%  integral dose outside (Gy L), V20 outside (cc), V(D>=20 & LETd>=4) and V(D>=5.4 & LETd>=4) outside (cc)
%  distal 15 mm: mean dose, mean LETd (dose > 1 % Rx), LETd 98th percentile, mean LET x dose,
%  McNamara RBE-weighted mean dose (alpha/beta 2); target McNamara D95 and mean (alpha/beta 10)
    Dphys = Dtot / 1.1;  t = g.tMask;  o = g.outsideEll;  ds = g.distal;
    v = sort(Dtot(t));  n = numel(v);  q = @(p) v(max(1, min(n, ceil(p*n))));
    r.D95 = q(0.05); r.D98 = q(0.02); r.D2 = q(0.98); r.D50 = q(0.5); r.Dmean = mean(v);
    r.V95 = 100*mean(v >= 0.95*g.Rx);  r.HI = (r.D2 - r.D98)/r.D50;
    pv = sort(Dtot(g.pMask)); r.PTV_V95 = 100*mean(pv >= 0.95*g.Rx); r.PTV_D98 = pv(max(1, ceil(0.02*numel(pv))));
    lt = sort(LET(t)); r.LETd_target_mean = mean(lt); r.LETd_target_min2 = lt(max(1, ceil(0.02*numel(lt))));
    r.DinDout  = mean(Dphys(t)) / mean(Dphys(o));
    r.LETratio = mean(LET(t)) / mean(LET(o));
    ct = g.ct;
    F = griddedInterpolant({ct.y, ct.x, ct.z}, Dtot, 'linear', 'none');
    lat = cross(g.d, [0 0 1]); lat = lat / norm(lat);
    s = (-120:0.1:120)';  P = g.tc + s * lat;
    prof = F(P(:,2), P(:,1), P(:,3)); prof(isnan(prof)) = 0; prof = 100*prof/mean(prof(abs(s) < 5));
    r.penumbra_cm = penumbra8020(s, prof) / 10;
    s2 = (0:0.1:150)';  P = g.tc + s2 * g.d;
    pd = F(P(:,2), P(:,1), P(:,3)); pd(isnan(pd)) = 0; pd = 100*pd/mean(pd(s2 < 5));
    r.distalFalloff_mm = cross1(s2, pd, 20) - cross1(s2, pd, 80);
    r.integralDose_out_GyL = sum(Dphys(o)) * g.voxVol / 1000;
    r.V20_out_cc = nnz(Dtot(o) >= 20) * g.voxVol;
    r.V_D20_LET4_out_cc  = nnz(Dtot(o) >= 20  & LET(o) >= 4) * g.voxVol;
    r.V_D5p4_LET4_out_cc = nnz(Dtot(o) >= 5.4 & LET(o) >= 4) * g.voxVol;
    r.D_distal_mean = mean(Dtot(ds));
    hasD = ds & Dtot > 0.01*g.Rx;
    ld = sort(LET(hasD));
    if isempty(ld), ld = 0; end
    r.LETd_distal_mean = mean(ld);  r.LETd_distal_98 = ld(max(1, ceil(0.98*numel(ld))));
    r.LETxD_distal_mean = mean(Dphys(ds) .* LET(ds));
    Dmc = mcNamara(Dphys, LET, g.nFx, g.abN);  r.DMcN_distal_mean = mean(Dmc(ds));  r.DMcN_out_mean = mean(Dmc(o));
    Dmt = mcNamara(Dphys, LET, g.nFx, g.abT);  v = sort(Dmt(t));
    r.DMcN_target_D95 = v(max(1, ceil(0.05*numel(v))));  r.DMcN_target_mean = mean(v);
end

function D = mcNamara(Dphys, LET, nFx, ab)
% McNamara et al. 2015 variable RBE, total RBE-weighted dose
    d = Dphys / nFx;
    RBEmax = 0.99064 + 0.35605 ./ ab .* LET;
    RBEmin = 1.1012 - 0.0038703 .* sqrt(ab) .* LET;
    RBE = (sqrt(ab^2 + 4*d*ab.*RBEmax + 4*d.^2.*RBEmin.^2) - ab) ./ (2*d);
    RBE(d <= 0) = 0;
    D = nFx * d .* RBE;
end

function T = addRow(T, mName, loc, g, depth, wLET, eng, r, nSpots, sumW, tOpt, tDij)
    row = struct2table(r);
    meta = table({mName}, {loc}, g, depth, wLET, {eng}, nSpots, sumW, tOpt, tDij, ...
        'VariableNames', {'machine','location','gantry','centreDepth_mm','LETweight','engine', ...
                          'nSpots','sumW_scaled','time_opt_s','time_dij_s'});
    T = [T; [meta row]];
end

function s = tissue(priority, alphaBeta)
    s = struct('TissueClass', 1, 'alphaX', 0.1, 'betaX', 0.1/alphaBeta, ...
               'Priority', priority, 'Visible', 1, 'visibleColor', [1 0 0]);
end

function d = beamDirection(ct, body, c0, g)
    [Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
    t = (Xg - c0(1)).^2 + (Yg - c0(2)).^2 + (Zg - c0(3)).^2 <= 10^2 & body;
    cstD = cell(1, 6);
    cstD(1,:) = {0, 'T', 'TARGET', {find(t)}, struct('TissueClass',1,'alphaX',0.1,'betaX',0.01, ...
        'Priority',1,'Visible',1,'visibleColor',[1 0 0]), {}};
    p = struct('radiationMode','protons','machine','Generic','numOfFractions',1, ...
               'bioModel','constRBE','multScen','nomScen');
    p.propStf.gantryAngles = g; p.propStf.couchAngles = 0; p.propStf.bixelWidth = 10;
    p.propStf.isoCenter = c0;
    st = matRad_generateStf(ct, cstD, p);
    d = st(1).isoCenter - st(1).sourcePoint; d = d / norm(d);
end

function e = entryPoint(ct, body, tc, d)
    F = griddedInterpolant({ct.y, ct.x, ct.z}, double(body), 'nearest', 'nearest');
    s = 0;
    while F(tc(2) - s*d(2), tc(1) - s*d(1), tc(3) - s*d(3)) > 0.5 && s < 1000
        s = s + 0.5;
    end
    e = tc - s * d;
end

function p = penumbra8020(s, prof)
    c = find(s >= 0, 1);
    sr = s(c:end);         qr = prof(c:end);
    sl = -flipud(s(1:c));  ql = flipud(prof(1:c));
    p = mean([cross1(sr, qr, 20) - cross1(sr, qr, 80), cross1(sl, ql, 20) - cross1(sl, ql, 80)]);
end

function xc = cross1(xs, qs, lev)
% first crossing below lev, searching outward from the start (skips the plateau)
    k = find(qs < lev, 1);
    if isempty(k) || k == 1, xc = NaN; return; end
    xc = interp1(qs([k-1 k]), xs([k-1 k]), lev);
end
