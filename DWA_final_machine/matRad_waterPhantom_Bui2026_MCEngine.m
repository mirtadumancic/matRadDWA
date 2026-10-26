%% matRad_waterPhantom_Bui2026_MCEngine.m
% MONTE CARLO VERSION of matRad_waterPhantom_Bui2026.m (Oct 24, 2026).
% Same phantom, targets, objectives, normalisation and metrics; the dose is
% computed with matRad's MCsquare Monte Carlo engine instead of HongPB.
% MCsquare uses only the in-air beam model of each machine (initFocus /
% emittance + energy spread) and simulates the in-water scatter and nuclear
% halo itself (sigma1/sigma2/weight are not used). That matches how the
% machines were set up in RayStation: single Gaussian in air, no in-water
% lateral profiles, the dose engine does the rest.
%   mcMode = 'dij'    : optimise on an MCsquare dose-influence matrix, then
%                       recalculate the final dose with MCsquare (forward,
%                       nHistFinal histories, no sparse cut-off)
%   mcMode = 'recalc' : optimise on HongPB (as the PB script), then
%                       recalculate the final dose with MCsquare (faster)
% All metrics use the final MCsquare dose; DinDout_ellipsoid_opt is the same
% metric on the dose the optimiser saw, for comparison.
% Outputs: DWA_final_machine/WaterPhantom_Bui2026_MCEngine/
%
% ---- original header ----
% matRad_waterPhantom_Bui2026.m
% Water-phantom case study after Bui et al. 2026 (Phys. Med. Biol. 71 175030),
% section 2.3 and appendix A.4, for comparison with her RayStation Figs 2-4:
%   D_in/D_out, LETd_in/LETd_out and lateral penumbra (80-20 %, cm)
%   for the 5 cm target at 4 locations: anterior, center, posterior, oblique.
%
% PHANTOM (paper 2.3): oval water phantom, 30 cm AP x 50 cm LR x 30 cm SI
% (elliptic cylinder along SI), in air. One spherical target per plan:
%   anterior   center 5 cm deep along the gantry-0 beam (beam from anterior)
%   center     15 cm deep
%   posterior  25 cm deep
%   oblique    15 cm to the patient's right of the phantom centre, gantry 330
%
% PLANNING (paper A.4, approximated in matRad - see objectives below):
%   54 Gy(RBE) in 30 fractions, constant RBE 1.1, physical-dose based.
%   Single field per target. Robust optimisation (3 mm / 3.5 %) is an option
%   (robust = true) but OFF by default: it multiplies memory and run time.
%   After optimisation each plan is scaled so that D95 of the target = 51.3 Gy
%   (95 % of the target receives 95 % of the prescription), as in the paper.
%
% METRICS (paper 2.3):
%   D_in/D_out   mean dose inside the target / mean dose in the rest of the
%                phantom (body minus target)
%   LET ratio    mean LETd inside / mean LETd outside, where voxels without
%                LET count as 0 keV/um (RayStation convention, paper 2.3/4)
%   penumbra     80-20 % distance on the lateral dose profile through the
%                target centre in the axial plane, perpendicular to the beam,
%                profile normalised to the dose at the target centre, mean of
%                both sides, in cm
%
% GEOMETRY NOTE: matRad's initFocus tables cover 200 mm either side of
% isocentre. Each plan therefore puts the isocentre on the beam axis through
% the target centre, isoDepth_mm below the entry point (default 150 mm), so the
% phantom surface always sits at Z = -150 mm. Spot sizes at the surface are
% the same for every target. RayStation placed the snout >= 10 cm from the
% surface; the exact air gap is not reproduced.
%
% MACHINES: DWA_proton_v2 (1 mm) and IBA_cyclotron (paper's reference line);
% add 'Generic' if wanted. IBA's lowest energy (99.6 MeV, R80 ~77 mm) cannot
% reach the anterior target without a range shifter; that plan is skipped.
% The paper's 1 mm machine started at 68 MeV, so RayStation has no 1 mm
% anterior point either.
%
% Outputs: DWA_final_machine/WaterPhantom_Bui2026/
%   metrics.csv, metrics_bars.png/.fig, doseMaps.png/.fig, plan_<machine>_<loc>.mat

clear; clc; close all;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.disableGUI = true;          % no optimiser progress window: runs unattended
set(0, 'DefaultFigureVisible', 'off'); % figures are saved to disk, not shown

%% ---------------- settings --------------------------------------------
machines   = {'DWA_proton_v2', 'IBA_cyclotron'};          % add 'Generic' if wanted
labels     = {'DWA 1 mm (matRad MC)', 'IBA cyclotron (matRad MC)'};
bixelWidth = struct('DWA_proton_v2', 3, 'IBA_cyclotron', 5, 'Generic', 5);   % mm
longSpotSpacing_mm = 3;
res_mm     = 3;               % CT and dose grid
targetDiam_mm = 50;
locations  = {'anterior', 'center', 'posterior', 'oblique'};
isoDepth_mm = 150;            % isocentre this far below the entry point
robust     = false;           % true: 3 mm / 3.5 % worst case (slow, memory heavy)
planMargin = true;            % cheap stand-in for robust optimisation: plan on the target
                              % + 3 mm laterally and + (3 mm + 3.5 % of depth) along the beam;
                              % coverage (D95) and all metrics stay on the 5 cm target
marginLat_mm = 3;  rangeUnc = 0.035;
saveResultGUI = true;
mcMode       = 'dij';         % 'dij' or 'recalc' (see header)
nHistBeamlet = 1e4;           % MCsquare histories per spot for the dij
nHistFinal   = 1e7;           % MCsquare histories for the final forward dose

Rx_total  = 54;  nFx = 30;    % Gy(RBE), fractions
outDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_Bui2026_MCEngine');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- MCsquare needs text metadata in the machine files -------
% writeMCsquareData prints meta.description / created_by / created_on into the
% BDL header. Our commissioned machines lack some of these; add them once
% (metadata only, the beam data are untouched).
for m = 1:numel(machines)
    mf = fullfile(rootDir, 'matRad', 'basedata', ['protons_' machines{m} '.mat']);
    S = load(mf, 'machine');
    need = {'description', machines{m}; 'created_by', 'Mirta Dumancic'; 'created_on', datestr(now, 'yyyy-mm-dd')};
    changed = false;
    for k = 1:size(need, 1)
        if ~isfield(S.machine.meta, need{k,1}) || isempty(S.machine.meta.(need{k,1}))
            S.machine.meta.(need{k,1}) = need{k,2}; changed = true;
        end
    end
    if changed
        machine = S.machine; save(mf, 'machine');
        fprintf('Added BDL header metadata to %s\n', mf);
    end
end
clear S machine

%% ---------------- phantom ----------------------------------------------
semiLR = 250; semiAP = 150; halfSI = 150;      % mm
margin = 30;                                   % air around the phantom
Ny = ceil((2*semiAP + 2*margin)/res_mm); Nx = ceil((2*semiLR + 2*margin)/res_mm);
Nz = ceil((2*halfSI + 2*margin)/res_mm);
ct = struct();
ct.resolution.x = res_mm; ct.resolution.y = res_mm; ct.resolution.z = res_mm;
ct.cubeDim = [Ny Nx Nz];
ct.numOfCtScen = 1;
ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
c0 = [mean(ct.x) mean(ct.y) mean(ct.z)];       % phantom centre, world [x y z]
[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);        % cube order [y x z]
body = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 <= 1 & abs(Zg - c0(3)) <= halfSI;
rED = 0.0012 * ones(ct.cubeDim);  rED(body) = 1;
HU  = -1000  * ones(ct.cubeDim);  HU(body)  = 0;
ct.cube = {rED};  ct.cubeHU = {HU};
ct.hlut = [-1024 0.00; 0 1.00; 3000 2.5];
fprintf('Phantom: %d x %d x %d voxels (%d mm), body %d voxels\n', Ny, Nx, Nz, res_mm, nnz(body));

%% ---------------- beam directions (gantry 0 and 330) --------------------
dirOf = @(g) beamDirection(ct, body, c0, g);
dir0   = dirOf(0);
dir330 = dirOf(330);
fprintf('gantry 0 beam direction [x y z] = %s\n', mat2str(round(dir0, 3)));

% target centres (world, mm)
tgt.anterior  = c0 + (50  - semiAP) * dir0;
tgt.center    = c0;
tgt.posterior = c0 + (250 - semiAP) * dir0;
tgt.oblique   = c0 + [150 0 0];                % 15 cm lateral; with gantry 330 the beam crosses
                                               % from the middle of the phantom towards the target
                                               % (path ~167 mm), as drawn in Fig. 1A of the paper
ell = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 + ((Zg - c0(3))/halfSI).^2 <= 1;   % ellipsoid variant
gantry = struct('anterior', 0, 'center', 0, 'posterior', 0, 'oblique', 330);

%% ---------------- plan every machine x location -------------------------
rows = {};
mapStore = struct();
for m = 1:numel(machines)
    mName = machines{m};
    pTmp = struct('radiationMode', 'protons', 'machine', mName);
    mach = matRad_loadMachine(pTmp);
    minRange = min([mach.data.peakPos]);

    for L = 1:numel(locations)
        loc = locations{L};
        tc  = tgt.(loc);
        g   = gantry.(loc);

        % ---- beam geometry first (needed for the planning margin)
        d = dirOf(g);
        entry = entryPoint(ct, body, tc, d);
        depthToCentre = norm(tc - entry);

        % ---- structures (disjoint, so matRad's overlap priorities do not matter)
        tMask  = (Xg - tc(1)).^2 + (Yg - tc(2)).^2 + (Zg - tc(3)).^2 <= (targetDiam_mm/2)^2 & body;
        dist   = sqrt((Xg - tc(1)).^2 + (Yg - tc(2)).^2 + (Zg - tc(3)).^2) - targetDiam_mm/2;
        if planMargin
            R  = targetDiam_mm/2;
            ax = (Xg - tc(1))*d(1) + (Yg - tc(2))*d(2) + (Zg - tc(3))*d(3);          % along the beam
            pp = sqrt(max((Xg - tc(1)).^2 + (Yg - tc(2)).^2 + (Zg - tc(3)).^2 - ax.^2, 0));
            mAx = marginLat_mm + rangeUnc*depthToCentre;
            pMask = ((pp/(R + marginLat_mm)).^2 + (ax/(R + mAx)).^2 <= 1 & body) | tMask;
            mOff = marginLat_mm;
        else
            pMask = tMask; mOff = 0;
        end
        margin = body & ~pMask & dist <= mOff + 2;               % 0-2 mm around the planning target
        shell  = body & ~pMask & dist > mOff + 2 & dist <= mOff + 12;   % the 1 cm fall-off zone
        far    = body & ~pMask & dist > mOff + 12;               % beyond the fall-off zone
        outside    = body & ~tMask;                              % "outside the target" for the METRICS
        outsideEll = ell & ~tMask;                               % same, ellipsoidal phantom

        cst = cell(4, 6);
        cst(1,:) = {0, 'TARGET',  'TARGET', {find(pMask)},  tissue(1, 10), {}};
        cst(2,:) = {1, 'MARGIN',  'OAR',    {find(margin)}, tissue(2, 2),  {}};
        cst(3,:) = {2, 'SHELL',   'OAR',    {find(shell)},  tissue(2, 2),  {}};
        cst(4,:) = {3, 'FAR',     'OAR',    {find(far)},    tissue(2, 2),  {}};

        % ---- objectives (approximating paper appendix A.4; weights are matRad's scale)
        % target: min DVH 54 Gy @ 95 % (w 100, robust), min dose 51.3 (w 10),
        %         min DVH 52.5 @ 99 % (w 5), plus a soft max of 57 Gy (the paper's constraint)
        oT1 = DoseObjectives.matRad_MinDVH(100, Rx_total, 95);
        oT2 = DoseObjectives.matRad_SquaredUnderdosing(10, 0.95*Rx_total);
        oT3 = DoseObjectives.matRad_MinDVH(5, 52.5, 99);
        oT4 = DoseObjectives.matRad_SquaredDeviation(50, Rx_total);
        oT5 = DoseObjectives.matRad_SquaredOverdosing(100, 57);
        % outside (beyond 2 mm): max 20 Gy (w 1, robust), max 53 Gy, fall-off 54 -> 27 Gy in 1 cm
        oS1 = DoseObjectives.matRad_SquaredOverdosing(10, 53);
        oS2 = DoseObjectives.matRad_SquaredOverdosing(1, 20);
        oF1 = DoseObjectives.matRad_SquaredOverdosing(10, 27);
        oF2 = DoseObjectives.matRad_SquaredOverdosing(1, 20);
        oM1 = DoseObjectives.matRad_SquaredOverdosing(100, 57);
        if robust, oT1.robustness = 'COWC'; oS2.robustness = 'COWC'; oF2.robustness = 'COWC'; end
        cst{1,6} = {oT1, oT2, oT3, oT4, oT5};
        cst{2,6} = {oM1};
        cst{3,6} = {oS1, oS2};
        cst{4,6} = {oF1, oF2};

        % ---- isocentre on the beam axis through the target, isoDepth below entry
        iso = entry + isoDepth_mm * d;
        if minRange > depthToCentre - targetDiam_mm/2 + 2
            fprintf('\n[%s / %s] skipped: proximal edge at %.0f mm, lowest range %.0f mm (needs a range shifter)\n', ...
                mName, loc, depthToCentre - targetDiam_mm/2, minRange);
            rows(end+1, :) = {mName, loc, g, depthToCentre, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<SAGROW>
            continue;
        end

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
        pln.propStf.isoCenter    = iso;
        if strcmp(mcMode, 'dij')
            pln.propDoseCalc.engine = 'MCsquare';
            pln.propDoseCalc.numHistoriesPerBeamlet = nHistBeamlet;
        else
            pln.propDoseCalc.engine       = 'HongPB';
            pln.propDoseCalc.lateralModel = 'double';
        end
        pln.propDoseCalc.calcLET      = 1;
        pln.propDoseCalc.doseGrid.resolution = struct('x', res_mm, 'y', res_mm, 'z', res_mm);
        pln.propOpt.quantityOpt  = 'RBExD';
        pln.propOpt.runDAO       = 0;
        pln.propSeq.runSequencing = 0;
        if robust
            pln.multScen = matRad_multScen(ct, 'wcScen');
            pln.multScen.shiftSD    = [3 3 3];
            pln.multScen.rangeRelSD = 3.5;
            pln.multScen.rangeAbsSD = 0;
        else
            pln.multScen = 'nomScen';
        end

        fprintf('\n===== %s / %s (gantry %d, centre depth %.0f mm) =====\n', mName, loc, g, depthToCentre);
        stf = matRad_generateStf(ct, cst, pln);
        dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
        resultGUI = matRad_fluenceOptimization(dij, cst, pln);
        clear dij

        % ---- dose the optimiser saw (for comparison only)
        if isfield(resultGUI, 'RBExD'), Dopt = resultGUI.RBExD; else, Dopt = 1.1 * resultGUI.physicalDose; end

        % ---- final dose: MCsquare forward calculation with the optimised weights
        plnMC = pln;
        plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', 1, ...
            'numHistoriesDirect', nHistFinal, 'doseGrid', pln.propDoseCalc.doseGrid);
        wOpt = resultGUI.w;
        fprintf('  MCsquare forward dose, %.0e histories ...\n', nHistFinal);
        resultGUI = matRad_calcDoseForward(ct, cst, stf, plnMC, wOpt);
        resultGUI.w = wOpt;

        % ---- normalise: D95(target) = 95 % of Rx
        if isfield(resultGUI, 'RBExD')
            Dfx = resultGUI.RBExD;
        else
            Dfx = 1.1 * resultGUI.physicalDose;          % constant RBE 1.1
        end
        if ~isequal(size(Dfx), ct.cubeDim)
            error('Result is not on the CT grid (%s vs %s).', mat2str(size(Dfx)), mat2str(ct.cubeDim));
        end
        Dtot = nFx * Dfx;
        v = sort(Dtot(tMask));  D95 = v(max(1, ceil(0.05*numel(v))));
        sc = 0.95*Rx_total / D95;
        Dtot = sc * Dtot;  Dphys = Dtot / 1.1;
        DoptT = nFx * Dopt;  v = sort(DoptT(tMask));  DoptT = DoptT * 0.95*Rx_total / v(max(1, ceil(0.05*numel(v))));
        DinDout_ell_opt = mean(DoptT(tMask)) / mean(DoptT(outsideEll));
        LET = zeros(ct.cubeDim);
        if isfield(resultGUI, 'LET'), LET = resultGUI.LET; end
        LET(~isfinite(LET) | Dphys <= 0) = 0;          % no LET -> 0 keV/um (RayStation)

        % ---- metrics
        DinDout = mean(Dphys(tMask)) / mean(Dphys(outside));
        LETratio = mean(LET(tMask)) / mean(LET(outside));
        DinDout_ell  = mean(Dphys(tMask)) / mean(Dphys(outsideEll));
        LETratio_ell = mean(LET(tMask)) / mean(LET(outsideEll));
        lat = cross(d, [0 0 1]); lat = lat / norm(lat);          % in the axial plane, perpendicular to the beam
        s = (-120:0.1:120)';
        F = griddedInterpolant({ct.y, ct.x, ct.z}, Dtot, 'linear', 'none');
        P = tc + s * lat;
        prof = F(P(:,2), P(:,1), P(:,3)); prof(isnan(prof)) = 0;
        prof = 100 * prof / mean(prof(abs(s) < 5));
        pen_cm = penumbra8020(s, prof) / 10;
        V95 = 100 * mean(Dtot(tMask) >= 0.95*Rx_total);
        Dmean_out = mean(Dtot(outside));
        rows(end+1, :) = {mName, loc, g, depthToCentre, DinDout, LETratio, DinDout_ell, LETratio_ell, ...
                          pen_cm, V95, Dmean_out, mean(LET(tMask)), numel(resultGUI.w), DinDout_ell_opt}; %#ok<SAGROW>
        fprintf('  [MCsquare final dose; optimiser dose gave ellipsoid D_in/D_out = %.1f]\n', DinDout_ell_opt);
        fprintf('  D_in/D_out = %.1f (ellipsoid %.1f)   LET_in/LET_out = %.1f (ellipsoid %.1f)   penumbra = %.2f cm   V95 = %.1f %%\n', ...
            DinDout, DinDout_ell, LETratio, LETratio_ell, pen_cm, V95);

        % ---- keep an axial slice for the dose-map figure
        [~, kz] = min(abs(ct.z - tc(3)));
        mapStore.(mName).(loc) = struct('D', Dtot(:,:,kz), 'tc', tc, 'g', g);

        if saveResultGUI
            planFile = fullfile(outDir, sprintf('plan_%s_%s.mat', mName, loc));
            w = resultGUI.w; %#ok<NASGU>
            save(planFile, 'Dtot', 'LET', 'w', 'pln', 'iso', 'tc', '-v7.3');
        end
        clear resultGUI
    end
end

%% ---------------- table -------------------------------------------------
T = cell2table(rows, 'VariableNames', {'machine','location','gantry','centreDepth_mm', ...
    'DinDout','LETinLETout','DinDout_ellipsoid','LETinLETout_ellipsoid', ...
    'penumbra8020_cm','V95_pct','Dmean_outside_Gy','LETmean_target','nSpots','DinDout_ellipsoid_opt'});
disp(T);
writetable(T, fullfile(outDir, 'metrics.csv'));

%% ---------------- figure: metrics per location ---------------------------
f1 = figure('Color', 'w', 'Position', [40 40 1800 430]);
metricsList = {'DinDout', 'DinDout_ellipsoid', 'LETinLETout', 'penumbra8020_cm'};
yLab = {'D_{in}/D_{out} (elliptic cylinder)', 'D_{in}/D_{out} (ellipsoid)', 'LET_{in}/LET_{out}', 'penumbra 80-20 % [cm]'};
% RayStation values for the 5 cm target, DWA 1 mm, read off Fig. 2 of Bui et al. 2026
paperDinDout    = struct('anterior', NaN, 'center', 64,   'posterior', 30.8, 'oblique', 49.7);  % DWA 1 mm
paperDinDoutIBA = struct('anterior', 116, 'center', 53,   'posterior', 26.5, 'oblique', 40.3);  % IBA dashed line
cmap = [0.75 0.14 0.23; 0.93 0.55 0.13; 0.13 0.36 0.67];
for k = 1:numel(metricsList)
    subplot(1, numel(metricsList), k); hold on; grid on; box on;
    Y = nan(numel(locations), numel(machines));
    for L = 1:numel(locations)
        for m = 1:numel(machines)
            r = strcmp(T.machine, machines{m}) & strcmp(T.location, locations{L});
            if any(r), Y(L, m) = T.(metricsList{k})(r); end
        end
    end
    b = bar(1:numel(locations), Y);
    for m = 1:numel(machines), b(m).FaceColor = cmap(m,:); end
    if startsWith(metricsList{k}, 'DinDout')
        pv = cellfun(@(l) paperDinDout.(l), locations);
        plot((1:numel(locations)) - 0.14, pv, 'kd', 'MarkerFaceColor', cmap(1,:), 'MarkerSize', 8);
        pv = cellfun(@(l) paperDinDoutIBA.(l), locations);
        plot((1:numel(locations)) + 0.14, pv, 'kd', 'MarkerFaceColor', cmap(2,:), 'MarkerSize', 8);
    end
    set(gca, 'XTick', 1:numel(locations), 'XTickLabel', locations);
    ylabel(yLab{k}); title(yLab{k});
    if k == 1, legend([labels, {'RayStation DWA 1 mm (Fig. 2)', 'RayStation IBA (Fig. 2, dashed)'}], 'Location', 'northeast', 'FontSize', 8); end
end
sgtitle('Water phantom after Bui et al. 2026, 5 cm target (matRad, MCsquare); diamonds: RayStation (Fig. 2), DWA 1 mm and IBA');
saveFig(f1, fullfile(outDir, 'metrics_bars')); close(f1);

%% ---------------- figure: axial dose maps --------------------------------
f2 = figure('Color', 'w', 'Position', [40 40 1500 300*numel(machines)]);
for m = 1:numel(machines)
    for L = 1:numel(locations)
        subplot(numel(machines), numel(locations), (m-1)*numel(locations) + L);
        if isfield(mapStore, machines{m}) && isfield(mapStore.(machines{m}), locations{L})
            M = mapStore.(machines{m}).(locations{L});
            imagesc(ct.x - c0(1), ct.y - c0(2), M.D); axis image; hold on;
            contour(ct.x - c0(1), ct.y - c0(2), M.D, [0.2 0.5 0.95]*Rx_total, 'k', 'LineWidth', 0.7);
            th = linspace(0, 2*pi, 100);
            plot(M.tc(1) - c0(1) + targetDiam_mm/2*cos(th), M.tc(2) - c0(2) + targetDiam_mm/2*sin(th), 'w--');
            caxis([0 1.05*Rx_total]); colormap(gca, turbo);
        else
            text(0.5, 0.5, 'not planned', 'HorizontalAlignment', 'center', 'Units', 'normalized'); axis off;
        end
        title(sprintf('%s, %s', strrep(machines{m}, '_', '\_'), locations{L}), 'FontSize', 9);
    end
end
sgtitle('Axial dose (total, Gy(RBE)); isodoses 20 / 50 / 95 % of 54 Gy');
saveFig(f2, fullfile(outDir, 'doseMaps')); close(f2);
set(0, 'DefaultFigureVisible', 'on');
fprintf('\nSaved to %s\n', outDir);

%% ---------------- local functions ---------------------------------------
function s = tissue(priority, alphaBeta)
    s = struct('TissueClass', 1, 'alphaX', 0.1, 'betaX', 0.1/alphaBeta, ...
               'Priority', priority, 'Visible', 1, 'visibleColor', [1 0 0]);
end

function d = beamDirection(ct, body, c0, g)
% unit vector from source to isocentre for gantry angle g (couch 0), taken
% from matRad's own steering geometry with a small dummy target at c0
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
% walk from the target centre back towards the source until leaving the body
    F = griddedInterpolant({ct.y, ct.x, ct.z}, double(body), 'nearest', 'nearest');
    s = 0;
    while F(tc(2) - s*d(2), tc(1) - s*d(1), tc(3) - s*d(3)) > 0.5 && s < 1000
        s = s + 0.5;
    end
    e = tc - s * d;
end

function p = penumbra8020(s, prof)
% 80-20 % distance on both sides of the centre, averaged [mm]
    c = find(s >= 0, 1);
    sr = s(c:end);         qr = prof(c:end);
    sl = -flipud(s(1:c));  ql = flipud(prof(1:c));
    p = mean([cross1(sr, qr, 20) - cross1(sr, qr, 80), cross1(sl, ql, 20) - cross1(sl, ql, 80)]);
end

function xc = cross1(xs, qs, lev)
    k = find(qs < lev, 1);
    if isempty(k) || k == 1, xc = NaN; return; end
    xc = interp1(qs([k-1 k]), xs([k-1 k]), lev);
end

function saveFig(f, base)
    exportgraphics(f, [base '.png'], 'Resolution', 150);
    savefig(f, [base '.fig']);             % open later with openfig(file, 'visible')
end
