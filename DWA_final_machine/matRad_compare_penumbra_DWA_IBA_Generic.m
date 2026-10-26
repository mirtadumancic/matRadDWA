%% matRad_compare_penumbra_DWA_IBA_Generic.m
% Lateral profiles and 80-20% penumbra of the SAME field with three machines:
% DWA_proton_v2 (1 mm spot), IBA_cyclotron (4.6-11 mm spot) and Generic.
%
% Field: uniform-fluence SOBP on the BOXPHANTOM target (63 mm cube, 90-153 mm
% deep), no optimiser, so the field edge is set by the beam model alone.
% Layer weights: NNLS on each machine's own IDDs for a flat SOBP; every spot
% in a layer gets the same weight; forward dose (no dij).
%
% Outputs (DWA_final_machine/Penumbra_DWA_IBA_Generic/):
%   penumbra_profiles.png/.fig   lateral profiles at 5 depths, linear and log
%   penumbra_vs_depth.png/.fig   80-20% penumbra and 50% width vs depth
%   penumbra_vs_depth.csv        the numbers behind the second figure
%   penumbra_summary.csv         one row per machine at the reference depths
%
% Penumbra definition: profile normalised to its own central dose at that
% depth (mean over |x| < 10 mm); distance between the 80% and 20% crossings,
% averaged over left/right and over the two lateral axes (x and z).
%
% Dose grid: doseRes (default 2 mm). Profiles are interpolated to 0.1 mm;
% the grid still blurs every machine by roughly the same ~0.5-1 mm.

clear; clc; close all;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_rc;

%% ---------------- settings --------------------------------------------
machines = {'DWA_proton_v2', 'IBA_cyclotron', 'Generic'};
labels   = {'DWA (1 mm)', 'IBA cyclotron', 'Generic'};
cols     = {[0.75 0.14 0.23], [0.93 0.55 0.13], [0.13 0.36 0.67]};

bixelWidth  = 3;      % mm, same for all machines (DWA needs <= ~3 mm)
doseRes     = 2;      % mm, dose grid
isoDepth_mm = 120;    % entrance surface -> isocentre (BOXPHANTOM)
tgtProx_mm  = 90;     % target depth range
tgtDist_mm  = 153;
tgtHalf_mm  = 31.5;   % target half-width

refDepths_mm  = [20 60 95 122 148];        % profiles shown at these depths
sweepDepths_mm = (5:5:150)';               % penumbra vs depth
lateral_mm    = (-90:0.1:90)';

outDir = fullfile(rootDir, 'DWA_final_machine', 'Penumbra_DWA_IBA_Generic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

load('BOXPHANTOM.mat');                    % ct, cst
ct = matRad_getWorldAxes(ct);

res = struct([]);

%% ---------------- run each machine --------------------------------------
for m = 1:numel(machines)
    pln = struct();
    pln.radiationMode  = 'protons';
    pln.machine        = machines{m};
    pln.bioModel       = 'constRBE';
    pln.multScen       = 'nomScen';
    pln.numOfFractions = 1;
    pln.propDoseCalc.engine                 = 'HongPB';
    pln.propDoseCalc.lateralModel           = 'double';
    pln.propDoseCalc.geometricLateralCutOff = 100;
    pln.propDoseCalc.calcLET                = 0;
    pln.propDoseCalc.doseGrid.resolution    = struct('x', doseRes, 'y', doseRes, 'z', doseRes);
    pln.propStf.gantryAngles = 0;
    pln.propStf.couchAngles  = 0;
    pln.propStf.bixelWidth   = bixelWidth;
    pln.propStf.isoCenter    = matRad_getIsoCenter(cst, ct, 0);
    pln.propOpt.runDAO        = 0;
    pln.propSeq.runSequencing = 0;

    stf     = matRad_generateStf(ct, cst, pln);
    machine = matRad_loadMachine(pln);

    % ---- layer weights (flat SOBP from the machine's own IDDs)
    E    = unique(cell2mat(cellfun(@(v) v(:)', {stf(1).ray.energy}, 'UniformOutput', false)));
    zFit = (tgtProx_mm + 2 : 1 : tgtDist_mm + 1)';
    A    = zeros(numel(zFit), numel(E));
    allE = [machine.data.energy];
    for k = 1:numel(E)
        [~, ix] = min(abs(allE - E(k)));
        dk = machine.data(ix);
        A(:, k) = interp1(double(dk.depths(:)) + double(dk.offset), double(dk.Z(:)), zFit, 'linear', 0);
    end
    layerW = lsqnonneg(A ./ max(A(:)), ones(numel(zFit), 1));
    w = [];
    for j = 1:numel(stf(1).ray)
        for e = reshape(stf(1).ray(j).energy, 1, [])
            w(end+1, 1) = layerW(abs(E - e) < 1e-6); %#ok<AGROW>
        end
    end
    fprintf('\n[%s] %d layers, %d spots\n', machines{m}, nnz(layerW), numel(w));

    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);
    D = resultGUI.physicalDose;
    if ~isequal(size(D), ct.cubeDim)
        error('Dose cube is not on the CT grid (%s vs %s).', mat2str(size(D)), mat2str(ct.cubeDim));
    end

    beamDir = stf(1).isoCenter - stf(1).sourcePoint; beamDir = beamDir / norm(beamDir);
    entry   = stf(1).isoCenter - isoDepth_mm * beamDir;
    FD = griddedInterpolant({ct.y, ct.x, ct.z}, D, 'linear', 'none');
    sampleD = @(P) max(FD(P(:,2), P(:,1), P(:,3)), 0);
    profAt  = @(d, ax) sampleD(entry + d * beamDir + lateral_mm * ax);

    % ---- profiles at the reference depths (z axis), own-centre normalised
    prof = zeros(numel(lateral_mm), numel(refDepths_mm));
    for q = 1:numel(refDepths_mm)
        p = profAt(refDepths_mm(q), [0 0 1]);
        prof(:, q) = 100 * p / mean(p(abs(lateral_mm) < 10));
    end

    % ---- penumbra and widths vs depth (x and z averaged)
    nD = numel(sweepDepths_mm);
    pen = nan(nD,1); w50 = nan(nD,1); w95 = nan(nD,1); w20 = nan(nD,1);
    for q = 1:nD
        sz = edgeMetrics(lateral_mm, profAt(sweepDepths_mm(q), [0 0 1]));
        sx = edgeMetrics(lateral_mm, profAt(sweepDepths_mm(q), [1 0 0]));
        pen(q) = mean([sz.pen8020 sx.pen8020]);
        w50(q) = mean([sz.w50 sx.w50]);
        w95(q) = mean([sz.w95 sx.w95]);
        w20(q) = mean([sz.w20 sx.w20]);
    end

    res(m).machine = machines{m}; res(m).label = labels{m};
    res(m).prof = prof; res(m).pen = pen; res(m).w50 = w50; res(m).w95 = w95; res(m).w20 = w20;
end

%% ---------------- tables ------------------------------------------------
Tsweep = table(sweepDepths_mm, 'VariableNames', {'depth_mm'});
for m = 1:numel(res)
    tag = matlab.lang.makeValidName(res(m).machine);
    Tsweep.(['pen8020_' tag]) = res(m).pen;
    Tsweep.(['width50_' tag]) = res(m).w50;
    Tsweep.(['width95_' tag]) = res(m).w95;
end
writetable(Tsweep, fullfile(outDir, 'penumbra_vs_depth.csv'));

Tsum = table();
for m = 1:numel(res)
    row = table(string(res(m).machine), 'VariableNames', {'machine'});
    for q = 1:numel(refDepths_mm)
        row.(sprintf('pen8020_at_%dmm', refDepths_mm(q))) = interp1(sweepDepths_mm, res(m).pen, refDepths_mm(q));
    end
    row.pen8020_mean_inTarget = mean(res(m).pen(sweepDepths_mm >= tgtProx_mm & sweepDepths_mm <= tgtDist_mm), 'omitnan');
    Tsum = [Tsum; row]; %#ok<AGROW>
end
disp(Tsum);
writetable(Tsum, fullfile(outDir, 'penumbra_summary.csv'));

%% ---------------- figure 1: lateral profiles ----------------------------
f1 = figure('Color', 'w', 'Position', [40 40 1650 640]);
nQ = numel(refDepths_mm);
for q = 1:nQ
    for row = 1:2
        ax = subplot(2, nQ, (row-1)*nQ + q); hold on; grid on; box on;
        xline([-tgtHalf_mm tgtHalf_mm], ':', 'HandleVisibility', 'off');
        for m = 1:numel(res)
            y = res(m).prof(:, q);
            if row == 2, y = max(y, 1e-2); end
            plot(lateral_mm, y, 'LineWidth', 1.6, 'Color', cols{m}, 'DisplayName', res(m).label);
        end
        if row == 1
            ylim([0 110]); xlim([-60 60]);
            title(sprintf('depth %d mm%s', refDepths_mm(q), ...
                ternary(refDepths_mm(q) >= tgtProx_mm, ' (in target)', ' (entrance)')));
            yline([20 80], '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        else
            set(ax, 'YScale', 'log'); ylim([0.1 150]); xlim([-90 90]);
        end
        if q == 1, ylabel(ternary(row == 1, 'dose [% of centre]', 'dose [% of centre], log')); end
        if row == 2, xlabel('lateral position [mm]'); end
    end
end
legend(subplot(2, nQ, 1), 'Location', 'south');
sgtitle('Lateral profiles of the same uniform-fluence SOBP (63 mm target, 90-153 mm deep)');
saveFig(f1, fullfile(outDir, 'penumbra_profiles'));

%% ---------------- figure 2: penumbra vs depth ---------------------------
f2 = figure('Color', 'w', 'Position', [40 40 1250 450]);
subplot(1,2,1); hold on; grid on; box on;
patch([tgtProx_mm tgtDist_mm tgtDist_mm tgtProx_mm], [0 0 40 40], [0.93 0.93 0.93], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
for m = 1:numel(res)
    plot(sweepDepths_mm, res(m).pen, 'o-', 'LineWidth', 1.6, 'MarkerSize', 4, ...
        'Color', cols{m}, 'DisplayName', res(m).label);
end
ylim([0 max(cellfun(@max, {res.pen}))*1.15]);
xlabel('depth in water [mm]'); ylabel('lateral penumbra 80-20% [mm]');
title('Penumbra vs depth (grey: target)'); legend('Location', 'northwest');
subplot(1,2,2); hold on; grid on; box on;
for m = 1:numel(res)
    plot(sweepDepths_mm, res(m).w95, '-',  'LineWidth', 1.6, 'Color', cols{m}, 'DisplayName', [res(m).label ', 95%']);
    plot(sweepDepths_mm, res(m).w50, '--', 'LineWidth', 1.2, 'Color', cols{m}, 'DisplayName', [res(m).label ', 50%']);
end
yline(2*tgtHalf_mm, ':', 'target width', 'HandleVisibility', 'off');
xlabel('depth in water [mm]'); ylabel('field width [mm]');
title('Field width at 95% (solid) and 50% (dashed)'); legend('Location', 'southeast', 'FontSize', 8);
sgtitle('Same field, three beam models');
saveFig(f2, fullfile(outDir, 'penumbra_vs_depth'));

fprintf('\nSaved to %s\n', outDir);

%% ---------------- local functions ---------------------------------------
function s = edgeMetrics(x, p)
% 80-20 penumbra and 95/50/20% widths of a profile normalised to its centre
    p0 = mean(p(abs(x) < 10));
    q  = 100 * p(:) / p0;  x = x(:);
    c  = find(x >= 0, 1);
    xr = x(c:end);            qr = q(c:end);
    xl = -flipud(x(1:c));     ql = flipud(q(1:c));
    r = @(lev) crossing(xr, qr, lev);  l = @(lev) crossing(xl, ql, lev);
    s.pen8020 = mean([r(20) - r(80), l(20) - l(80)]);
    s.w95 = r(95) + l(95);
    s.w50 = r(50) + l(50);
    s.w20 = r(20) + l(20);
end

function xc = crossing(xs, qs, lev)
% first point, going outward from the centre, where the profile drops below lev
    k = find(qs < lev, 1);
    if isempty(k) || k == 1, xc = NaN; return; end
    xc = interp1(qs([k-1 k]), xs([k-1 k]), lev);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function saveFig(f, base)
    exportgraphics(f, [base '.png'], 'Resolution', 150);
    savefig(f, [base '.fig']);
end
