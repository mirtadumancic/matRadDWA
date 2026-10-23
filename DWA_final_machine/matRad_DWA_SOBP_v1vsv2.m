%% matRad_DWA_SOBP_v1vsv2.m
% SOBP check of the DWA beam model: v1 vs v2 vs Generic, BOXPHANTOM.
%
% Question: does v2 (depth-dependent in-water halo) fix the shallow lateral
% penumbra and the ~20% shoulder Remo saw with v1 (constant in-air halo,
% weight ~0.5)?
%
% Design: the SAME uniform-fluence SOBP for every machine, no optimiser.
%  - matRad places spots on the 63 x 63 x 63 mm BOXPHANTOM target (3 mm grid).
%  - Each energy layer gets ONE weight, the same for every spot in the layer,
%    so the lateral fluence is flat and the field edge is set by the beam
%    model alone. An optimiser would sharpen the edges differently per
%    machine and hide exactly what we want to see.
%  - Layer weights: non-negative least squares on the machine's own IDDs
%    (data.Z) for a flat SOBP over the target depth (90-153 mm).
%  - Forward dose (no dij), so memory stays small even with 100 mm cut-offs.
%
% Outputs (folder SOBP_v1vsv2 next to this script):
%  - central-axis depth dose, lateral profiles at mid-SOBP (linear and log),
%    LETd (machines that have LET), dose maps with isodoses
%  - SOBP_v1vsv2_metrics.csv: flatness, distal R80/R20, lateral 80-20
%    penumbra, 50% field width, dose 10/20/30 mm outside the field edge
%
% Caveat: dose and sampling grid are 3 mm (the BOXPHANTOM CT); profiles are
% interpolated. Penumbra values carry a ~1 mm grid blur, the same for all
% machines, so the COMPARISON is sound but absolute penumbrae are slightly
% overestimated. Set doseRes = 2 for a sharper (slower) run.

clear; clc; close all;
rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
matRad_rc;

%% ---------------- settings --------------------------------------------
machines = {'DWA_proton_v2', 'IBA_cyclotron', 'Generic'};
labels   = {'DWA v2 (depth-dependent halo)', 'IBA cyclotron', 'Generic'};
cols     = {[0.60 0.60 0.60], [0.75 0.14 0.23], [0.13 0.36 0.67]};

bixelWidth   = 3;        % mm, lateral spot spacing (DWA core sigma ~1.3-2 mm)
doseRes      = 3;        % mm, dose grid
isoDepth_mm  = 120;      % entrance surface -> isocentre (BOXPHANTOM)
tgtProx_mm   = 90;       % target proximal edge (depth)
tgtDist_mm   = 153;      % target distal edge (depth)
tgtHalf_mm   = 31.5;     % target half-width laterally
plateauMargin_mm = 4;    % flatness evaluated this far inside the target

depths_mm  = (0:0.5:200)';
lateral_mm = (-90:0.25:90)';

outDir = fullfile(rootDir, 'DWA_final_machine', 'SOBP_v1vsv2');
if ~exist(outDir, 'dir'), mkdir(outDir); end

load('BOXPHANTOM.mat');                       % ct, cst
ct = matRad_getWorldAxes(ct);
ixTarget = find(strcmp(cst(:,3), 'TARGET'), 1);

res = struct([]);

%% ---------------- run each machine --------------------------------------
for m = 1:numel(machines)
    pln = struct();
    pln.radiationMode = 'protons';
    pln.machine       = machines{m};
    pln.bioModel      = 'constRBE';
    pln.multScen      = 'nomScen';
    pln.numOfFractions = 1;

    pln.propDoseCalc.engine                 = 'HongPB';
    pln.propDoseCalc.lateralModel           = 'double';
    pln.propDoseCalc.geometricLateralCutOff = 100;   % mm, keeps the halo
    pln.propDoseCalc.calcLET                = 1;     % skipped with a warning if the machine has no LET
    pln.propDoseCalc.doseGrid.resolution    = struct('x', doseRes, 'y', doseRes, 'z', doseRes);

    pln.propStf.gantryAngles = 0;
    pln.propStf.couchAngles  = 0;
    pln.propStf.bixelWidth   = bixelWidth;
    pln.propStf.isoCenter    = matRad_getIsoCenter(cst, ct, 0);
    pln.propOpt.runDAO        = 0;
    pln.propSeq.runSequencing = 0;

    stf     = matRad_generateStf(ct, cst, pln);
    machine = matRad_loadMachine(pln);

    % ---- layer weights: NNLS on the machine's own IDDs for a flat SOBP
    E      = unique(cell2mat(cellfun(@(v) v(:)', {stf(1).ray.energy}, 'UniformOutput', false)));
    zFit   = (tgtProx_mm + 2 : 1 : tgtDist_mm + 1)';   % fit window reaches past the distal edge so no layer peak is left unconstrained
    A      = zeros(numel(zFit), numel(E));
    allE   = [machine.data.energy];
    for k = 1:numel(E)
        [~, ix] = min(abs(allE - E(k)));
        dk = machine.data(ix);
        zk = double(dk.depths(:)) + double(dk.offset);
        A(:, k) = interp1(zk, double(dk.Z(:)), zFit, 'linear', 0);
    end
    A = A ./ max(A(:));
    layerW = lsqnonneg(A, ones(numel(zFit), 1));

    % ---- bixel weight vector in matRad order (ray by ray, energy by energy)
    w = [];
    for j = 1:numel(stf(1).ray)
        for e = reshape(stf(1).ray(j).energy, 1, [])
            w(end+1, 1) = layerW(abs(E - e) < 1e-6); %#ok<AGROW>
        end
    end
    fprintf('\n[%s] %d layers (%d with weight > 0), %d rays, %d spots\n', machines{m}, ...
        numel(E), nnz(layerW), numel(stf(1).ray), numel(w));

    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);
    D = resultGUI.physicalDose;
    if ~isequal(size(D), ct.cubeDim)
        error('Dose cube is not on the CT grid (%s vs %s).', mat2str(size(D)), mat2str(ct.cubeDim));
    end
    hasLET = isfield(resultGUI, 'LET');

    % ---- sampling geometry (cubes are [y x z])
    beamDir = stf(1).isoCenter - stf(1).sourcePoint; beamDir = beamDir / norm(beamDir);
    entry   = stf(1).isoCenter - isoDepth_mm * beamDir;
    FD = griddedInterpolant({ct.y, ct.x, ct.z}, D, 'linear', 'none');
    sampleD = @(P) max(FD(P(:,2), P(:,1), P(:,3)), 0);
    pdd = sampleD(entry + depths_mm * beamDir);

    plateau = depths_mm >= tgtProx_mm + plateauMargin_mm & depths_mm <= tgtDist_mm - plateauMargin_mm;
    Dref = mean(pdd(plateau));
    pdd  = 100 * pdd / Dref;

    midDepth = 0.5 * (tgtProx_mm + tgtDist_mm);
    latZ = 100 * sampleD(entry + midDepth * beamDir + lateral_mm * [0 0 1]) / Dref;
    latX = 100 * sampleD(entry + midDepth * beamDir + lateral_mm * [1 0 0]) / Dref;

    if hasLET
        L = resultGUI.LET; L(~isfinite(L)) = 0;
        FL = griddedInterpolant({ct.y, ct.x, ct.z}, L, 'linear', 'none');
        sampleL = @(P) FL(P(:,2), P(:,1), P(:,3));
        letCAX = sampleL(entry + depths_mm * beamDir);
        letLat = sampleL(entry + midDepth * beamDir + lateral_mm * [0 0 1]);
    else
        letCAX = nan(size(depths_mm)); letLat = nan(size(lateral_mm));
    end

    % ---- central slice for the dose map (plane through iso, along z)
    [~, ixX] = min(abs(ct.x - stf(1).isoCenter(1)));
    slice = squeeze(D(:, ixX, :)) / Dref * 100;

    % ---- metrics
    lz = lateralMetrics(lateral_mm, latZ, tgtHalf_mm);
    lx = lateralMetrics(lateral_mm, latX, tgtHalf_mm);
    res(m).machine     = machines{m};
    res(m).label       = labels{m};
    res(m).nLayers     = nnz(layerW);
    res(m).flatness    = max(pdd(plateau)) - min(pdd(plateau));          % % of plateau
    res(m).entrance    = pdd(find(depths_mm >= 2, 1));
    res(m).R80         = distalCrossing(depths_mm, pdd, 80);
    res(m).R20         = distalCrossing(depths_mm, pdd, 20);
    res(m).distal8020  = res(m).R20 - res(m).R80;
    res(m).pen8020     = mean([lz.pen8020 lx.pen8020]);
    res(m).width50     = mean([lz.width50 lx.width50]);
    res(m).out10       = mean([lz.out10 lx.out10]);
    res(m).out20       = mean([lz.out20 lx.out20]);
    res(m).out30       = mean([lz.out30 lx.out30]);
    res(m).LETmid      = ifelse(hasLET, interp1(depths_mm, letCAX, midDepth), NaN);
    res(m).LETatR80    = ifelse(hasLET, interp1(depths_mm, letCAX, res(m).R80), NaN);
    res(m).pdd = pdd; res(m).latZ = latZ; res(m).latX = latX;
    res(m).letCAX = letCAX; res(m).letLat = letLat; res(m).slice = slice;
end

%% ---------------- table -------------------------------------------------
T = struct2table(rmfield(res, {'pdd','latZ','latX','letCAX','letLat','slice','label'}));
disp(T);
writetable(T, fullfile(outDir, 'SOBP_v1vsv2_metrics.csv'));

%% ---------------- figure 1: depth dose + lateral profiles ---------------
f1 = figure('Color', 'w', 'Position', [60 60 1500 460]);
subplot(1,3,1); hold on; grid on; box on;
patch([tgtProx_mm tgtDist_mm tgtDist_mm tgtProx_mm], [0 0 120 120], [0.93 0.93 0.93], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
for m = 1:numel(res)
    plot(depths_mm, res(m).pdd, 'LineWidth', 1.8, 'Color', cols{m}, 'DisplayName', res(m).label);
end
xlabel('depth in water [mm]'); ylabel('dose [% of SOBP plateau]'); ylim([0 115]); xlim([0 180]);
title('Central-axis depth dose'); legend('Location', 'southwest');

subplot(1,3,2); hold on; grid on; box on;
xline([-tgtHalf_mm tgtHalf_mm], ':', 'HandleVisibility', 'off');
for m = 1:numel(res)
    plot(lateral_mm, res(m).latZ, 'LineWidth', 1.8, 'Color', cols{m}, 'DisplayName', res(m).label);
end
xlabel('lateral position [mm]'); ylabel('dose [% of plateau]'); xlim([-70 70]); ylim([0 115]);
title(sprintf('Lateral profile at mid-SOBP (%.0f mm)', 0.5*(tgtProx_mm+tgtDist_mm)));

subplot(1,3,3); hold on; grid on; box on;
xline([-tgtHalf_mm tgtHalf_mm], ':', 'HandleVisibility', 'off');
for m = 1:numel(res)
    semilogy(lateral_mm, max(res(m).latZ, 1e-2), 'LineWidth', 1.8, 'Color', cols{m});
end
set(gca, 'YScale', 'log'); ylim([0.1 150]); xlim([-90 90]);
xlabel('lateral position [mm]'); ylabel('dose [% of plateau], log');
title('Same, log scale: shoulder and low-dose bath');
sgtitle('Uniform-fluence SOBP in water, 63 mm cube target, no optimiser');
saveFig(f1, fullfile(outDir, 'SOBP_v1vsv2_profiles'));

%% ---------------- figure 2: dose maps -----------------------------------
f2 = figure('Color', 'w', 'Position', [60 60 1500 430]);
depthAx = (ct.y(:)' - entry(2)) * sign(beamDir(2));   % depth along the beam (gantry 0: y axis)
for m = 1:numel(res)
    subplot(1, numel(res), m);
    imagesc(depthAx, ct.z - stf(1).isoCenter(3), res(m).slice'); axis image; hold on;
    contour(depthAx, ct.z - stf(1).isoCenter(3), res(m).slice', [20 50 80 95], 'k', 'LineWidth', 0.8);
    rectangle('Position', [tgtProx_mm -tgtHalf_mm tgtDist_mm-tgtProx_mm 2*tgtHalf_mm], ...
        'EdgeColor', 'w', 'LineStyle', '--');
    caxis([0 110]); colormap(gca, turbo); xlim([0 180]); ylim([-80 80]);
    xlabel('depth [mm]'); ylabel('lateral [mm]'); title(res(m).label);
    if m == numel(res), cb = colorbar; cb.Label.String = '% of plateau'; end
end
sgtitle('Dose maps, isodoses 20 / 50 / 80 / 95 %');
saveFig(f2, fullfile(outDir, 'SOBP_v1vsv2_doseMaps'));

%% ---------------- figure 3: LETd (machines that have LET) --------------
hasAnyLET = arrayfun(@(r) any(isfinite(r.letCAX) & r.letCAX > 0), res);
if any(hasAnyLET)
    f3 = figure('Color', 'w', 'Position', [60 60 1100 420]);
    subplot(1,2,1); hold on; grid on; box on;
    for m = find(hasAnyLET)
        plot(depths_mm, res(m).letCAX, 'LineWidth', 1.8, 'Color', cols{m}, 'DisplayName', res(m).label);
    end
    xline([tgtProx_mm tgtDist_mm], ':', 'HandleVisibility', 'off');
    xlabel('depth in water [mm]'); ylabel('LET_d [keV/\mum]'); xlim([0 180]);
    title('Central-axis LET_d'); legend('Location', 'northwest');
    subplot(1,2,2); hold on; grid on; box on;
    for m = find(hasAnyLET)
        plot(lateral_mm, res(m).letLat, 'LineWidth', 1.8, 'Color', cols{m});
    end
    xline([-tgtHalf_mm tgtHalf_mm], ':', 'HandleVisibility', 'off');
    xlabel('lateral position [mm]'); ylabel('LET_d [keV/\mum]'); xlim([-70 70]);
    title('LET_d across the field at mid-SOBP');
    sgtitle('LET_d (DWA v2 uses the range-matched Generic LET stopgap)');
    saveFig(f3, fullfile(outDir, 'SOBP_v1vsv2_LETd'));
end

fprintf('\nSaved figures and metrics to %s\n', outDir);

%% ---------------- local functions ---------------------------------------
function r = distalCrossing(z, d, level)
    k = find(d >= level, 1, 'last');
    if isempty(k) || k == numel(d), r = NaN; return; end
    r = interp1(d([k k+1]), z([k k+1]), level);
end

function s = lateralMetrics(x, p, edge)
    % p in % of plateau; normalise to the central plateau (|x| < 0.5 edge)
    p0 = mean(p(abs(x) < 0.5*edge)); q = 100 * p / p0;
    c  = find(x >= 0, 1);
    xr = x(c:end); qr = q(c:end);                 % right side, centre outward
    xl = -flipud(x(1:c)); ql = flipud(q(1:c));    % left side, mirrored
    cross = @(xs, qs, lev) interp1(qs([find(qs < lev, 1)-1, find(qs < lev, 1)]), ...
                                   xs([find(qs < lev, 1)-1, find(qs < lev, 1)]), lev);
    r80 = cross(xr, qr, 80); r50 = cross(xr, qr, 50); r20 = cross(xr, qr, 20);
    l80 = cross(xl, ql, 80); l50 = cross(xl, ql, 50); l20 = cross(xl, ql, 20);
    s.pen8020 = mean([r20 - r80, l20 - l80]);
    s.width50 = r50 + l50;
    at = @(d) mean([interp1(x, p, edge + d), interp1(x, p, -edge - d)]);
    s.out10 = at(10); s.out20 = at(20); s.out30 = at(30);
end

function v = ifelse(c, a, b)
    if c, v = a; else, v = b; end
end

function saveFig(f, base)
    exportgraphics(f, [base '.png'], 'Resolution', 150);
    savefig(f, [base '.fig']);
end
