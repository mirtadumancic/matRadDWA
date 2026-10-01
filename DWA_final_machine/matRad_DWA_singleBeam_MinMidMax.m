%% matRad_DWA_singleBeam_MinMidMax.m
%
% Single-pencil-beam check of a DWA machine at THREE energies - lowest,
% middle and highest usable - in matRad's BOXPHANTOM. For each energy it
% produces the four plots of matRad_example5_protons_DWA_BOXPHANTOM_singleBeam.m:
%
%   1. depth dose: matRad laterally integrated (1 mm grid) vs the machine's
%      own input curve, plus the central-axis curve for reference
%   2. beam's-eye view at the Bragg-peak depth
%   3. depth view through the beam axis
%   4. lateral profile at the Bragg-peak depth, log scale, with the
%      machine's analytic double Gaussian (core, halo, total) overlaid
%   5. depth view of a small 5 x 5 spot field (3 mm spacing) at the same
%      energy. For a single ~1 mm spot the hottest VOXEL is near the
%      entrance, because the beam spreads out with depth (panel 3); in a
%      field, neighbouring spots replace the dose scattered off each axis and
%      the Bragg peak is the hottest region again.
%
% Every file name carries the machine, the label (MIN/MID/MAX) and the
% energy actually used. A 3 x 5 overview (rows = energies, columns = plots)
% is saved as well.
%
% LATERAL CUTOFF - WHY THIS SCRIPT CHANGES IT
% With matRad's default cutoffs the halo was being thrown away: in the
% single-beam run at 150 MeV the lateral dose dropped to exactly zero beyond
% ~12 mm, where the machine's halo (weight ~0.5, sigma ~23 mm) should still
% give a few percent of the central dose out to 40+ mm. The depth view
% showed it too - the dose region widened in steps with depth, the signature
% of matRad's depth-dependent dosimetric cutoff. The depth-dose SHAPE cannot
% reveal this (the halo fraction is constant with depth, so cutting it
% scales every depth by the same factor). Here the dosimetric cutoff is
% switched off (1 = keep everything) and the geometric cutoff widened to
% 100 mm, about 4 halo sigmas. Panel 4 then shows whether matRad delivers
% the halo the machine describes, and the console prints the ratio.
%
% Output: <folder of this script>/singleBeam_MinMidMax_<machine>/

%% ------------------------- configuration -------------------------------
scriptDir = fileparts(mfilename('fullpath'));
cd('/Users/mirtadumancic/Work/MATLAB_Projects/matRadDWA');
run('matRad_rc');
load('BOXPHANTOM.mat');

dwaMachine = 'DWA_proton_doubleGauss_interp';

% 'auto' = lowest usable, middle and highest usable NATIVE (TOPAS) layer, so
% each is compared against a real TOPAS curve rather than an interpolated
% one. Or give three energies in MeV, e.g. [50 110 180].
energiesRequested = 'auto';
labels            = {'MIN','MID','MAX'};

waterThick_mm  = 240;   % BOXPHANTOM water block along the beam
peakMargin_mm  = 15;    % a usable layer must peak at least this far inside it

geometricLateralCutOff_mm = 100;  % ~4 halo sigmas
dosimetricLateralCutOff   = 1;    % 1 = no dosimetric cutoff (keep the halo).
                                  % If this matRad version rejects exactly 1,
                                  % use 0.9999.

fineRes = struct('x', 2, 'y', 1, 'z', 2);   % depth-curve grid; beam along y

% Panel 5: a small scanned field, fieldN x fieldN spots one CT voxel (3 mm)
% apart, same energy. Built by SUPERPOSITION of the single spot shifted by
% whole voxels - exact here because BOXPHANTOM is uniform water across the
% beam, so a laterally displaced spot deposits the same dose, displaced
% (beam divergence at SAD = 10 m changes a 6 mm offset by <0.2 mm over 25 cm).
% Not valid in a heterogeneous patient: there each spot must be computed.
fieldN = 5;

outDir = fullfile(scriptDir, ['singleBeam_MinMidMax_' dwaMachine]);
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------- plan skeleton -------------------------------
pln.radiationMode = 'protons';
pln.machine       = dwaMachine;
pln.bioModel      = 'constRBE';
pln.multScen      = 'nomScen';
pln.numOfFractions = 1;
pln.propDoseCalc.calcLET = 0;
pln.propDoseCalc.engine  = 'HongPB';
% FORCE the double-Gaussian kernel. matRad's default lateralModel ('fast')
% picks a SINGLE Gaussian whenever machine.data has a 'sigma' field - and
% the DWA double-Gaussian machines carry one (the in-water MCS term, kept as
% a fallback). Without this line matRad silently ignores the halo and
% computes the core only, whatever meta.dataType says.
pln.propDoseCalc.lateralModel = 'double';
pln.propDoseCalc.geometricLateralCutOff = geometricLateralCutOff_mm;
pln.propDoseCalc.dosimetricLateralCutOff = dosimetricLateralCutOff;
pln.propDoseCalc.doseGrid.resolution = struct('x', ct.resolution.x, ...
    'y', ct.resolution.y, 'z', ct.resolution.z);
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles  = 0;
pln.propStf.bixelWidth   = 5;
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);
pln.propStf.generator    = 'ParticleSingleSpot';
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;
pln.propOpt.quantityOpt   = 'physicalDose';

plnFine = pln;
plnFine.propDoseCalc.doseGrid.resolution = fineRes;

%% ------------------------- choose the energies -------------------------
machine = matRad_loadMachine(pln);
allE  = [machine.data.energy];
allPk = [machine.data.peakPos];
isNat = true(size(allE));
if isfield(machine.data,'interpolated')
    isNat = ~arrayfun(@(d) isequal(d.interpolated, true), machine.data);
end
% machine.data is a column, so arrayfun returns a column while [..] gives
% rows; force all three to rows, or '&' expands them into an N x N matrix
allE = allE(:)';  allPk = allPk(:)';  isNat = isNat(:)';
if ischar(energiesRequested)
    usable = isNat & allPk <= waterThick_mm - peakMargin_mm;
    Eu = allE(usable);  Pu = allPk(usable);
    [~, iMid] = min(abs(Pu - (min(Pu) + max(Pu))/2));   % halfway in RANGE
    energiesRequested = [min(Eu), Eu(iMid), max(Eu)];
end
fprintf('\nEnergies: MIN %.2f, MID %.2f, MAX %.2f MeV\n', energiesRequested);

%% ------------------------- composite figure ----------------------------
hComp = figure('Color','w','Position',[30 50 2300 1300]);
tl = tiledlayout(hComp, 3, 5, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('%s - single pencil beam in BOXPHANTOM at three energies', ...
    strrep(dwaMachine,'_','\_')), 'FontWeight','bold');

surfaceY = [];
summary = table('Size',[3 10], 'VariableTypes', [{'string'} repmat({'double'},1,9)], ...
    'VariableNames', {'label','E_MeV','braggDepth_matRad_mm','braggDepth_input_mm', ...
    'IDDoverInput_min','IDDoverInput_max','haloRatio_matRadOverModel','entrance_pct', ...
    'axisPeakOverEntrance_spot','axisPeakOverEntrance_field'});

for k = 1:3
    %% ---- dose calculation
    pln.propStf.energy     = energiesRequested(k);
    plnFine.propStf.energy = energiesRequested(k);
    stf  = matRad_generateStf(ct, cst, pln);
    eUse = stf(1).ray.energy;
    tag  = sprintf('%s_%s_E%06.2fMeV', dwaMachine, labels{k}, eUse);
    fprintf('\n=== %s: %.2f MeV ===\n', labels{k}, eUse);

    dij      = matRad_calcDoseInfluence(ct, cst, stf, pln);
    doseCube = reshape(full(dij.physicalDose{1}), dij.doseGrid.dimensions);
    dijFine  = matRad_calcDoseInfluence(ct, cst, stf, plnFine);

    %% ---- geometry (once): beam direction and water surface from the CT
    iso     = stf(1).isoCenter;
    beamDir = -stf(1).sourcePoint / norm(stf(1).sourcePoint);   % sourcePoint is iso-relative
    if abs(beamDir(2) - 1) > 1e-3
        error('This script assumes gantry 0 / couch 0 (beam along +y).');
    end
    [~, ixAx] = min(abs(dij.doseGrid.x - iso(1)));
    [~, izAx] = min(abs(dij.doseGrid.z - iso(3)));
    if isempty(surfaceY)
        if isfield(ct,'cube') && ~isempty(ct.cube)
            isWater = ct.cube{1}(:, ixAx, izAx) > 0.5;
        else
            isWater = ct.cubeHU{1}(:, ixAx, izAx) > -500;
        end
        iyWater  = find(isWater, 1, 'first');
        surfaceY = dij.doseGrid.y(iyWater) - ct.resolution.y/2;
        isoDepth = iso(2) - surfaceY;
        fprintf('Water surface at y = %.1f mm; isocentre %.1f mm deep\n', surfaceY, isoDepth);
    end

    %% ---- depth curves
    dimsF = dijFine.doseGrid.dimensions;
    [vIdx, ~, vDose] = find(dijFine.physicalDose{1}(:,1));
    [iyF, ~, ~] = ind2sub(dimsF, vIdx);
    idd   = accumarray(iyF, vDose, [dimsF(1) 1]) * fineRes.x * fineRes.z;
    iddZ  = dijFine.doseGrid.y(:) - surfaceY;
    keep  = iddZ >= 0 & iddZ <= 280;
    iddZ  = iddZ(keep);  idd = idd(keep);
    [iddMax, iMax] = max(idd);
    braggDepth = iddZ(iMax);

    cax   = squeeze(doseCube(:, ixAx, izAx));
    caxZ  = dij.doseGrid.y(:) - surfaceY;
    okC   = caxZ >= 0 & caxZ <= 280;
    caxZ  = caxZ(okC);  cax = cax(okC);

    [~, eIx] = min(abs(allE - eUse));
    d  = machine.data(eIx);
    inZ = d.depths(:);  inD = d.Z(:);
    [inMax, inMaxI] = max(inD);
    isInterp = ~isNat(eIx);

    inOnGrid = interp1(inZ, inD/inMax, iddZ, 'linear', NaN);
    % Compare from 1 mm depth (the first 1 mm voxel straddles the surface and
    % matRad leaves it empty) to 5 mm before the peak (on the steep proximal
    % edge a 1 mm depth offset alone moves the ratio by several percent).
    okCmp    = iddZ >= 1 & iddZ < inZ(inMaxI) - 5 & isfinite(inOnGrid) & inOnGrid > 0.05;
    ratio    = (idd(okCmp)/iddMax) ./ inOnGrid(okCmp);
    if isempty(ratio), ratio = NaN; end
    fprintf('Bragg peak: matRad %.1f mm, input %.1f mm; IDD/input %.3f - %.3f\n', ...
        braggDepth, inZ(inMaxI), min(ratio), max(ratio));

    %% ---- slices on the coarse (CT) grid
    [~, iyPk] = min(abs(dij.doseGrid.y - (surfaceY + braggDepth)));
    zPkVox    = dij.doseGrid.y(iyPk) - surfaceY;         % depth of that voxel centre
    sliceBEV  = squeeze(doseCube(iyPk,:,:));
    bevMax    = 1.001 * max(sliceBEV(:));
    cubeMax   = 1.001 * max(doseCube(:));

    %% ---- small scanned field by superposition (panel 5)
    halfN = (fieldN - 1)/2;
    fieldCube = zeros(size(doseCube));
    for di = -halfN:halfN
        for dk = -halfN:halfN
            fieldCube = fieldCube + circshift(doseCube, [0 di dk]);   % [y x z]: shift x and z
        end
    end
    axSpot  = squeeze(doseCube(:, ixAx, izAx));
    axField = squeeze(fieldCube(:, ixAx, izAx));
    [~, iyEnt] = min(abs(dij.doseGrid.y - (surfaceY + ct.resolution.y/2)));   % first water voxel
    fprintf('On-axis dose, Bragg peak / entrance: single spot %.2f, %dx%d field %.2f\n', ...
        axSpot(iyPk)/axSpot(iyEnt), fieldN, fieldN, axField(iyPk)/axField(iyEnt));

    %% ---- analytic lateral model at the Bragg-peak voxel depth
    SSD      = norm(stf(1).sourcePoint) - isoDepth;     % source-to-surface distance
    sigIni   = interp1(d.initFocus.dist(1,:), d.initFocus.sigma(1,:), SSD, 'linear', 'extrap');
    s1       = interp1(d.depths, d.sigma1, zPkVox);
    s2       = interp1(d.depths, d.sigma2, zPkVox);
    w        = interp1(d.depths, d.weight, zPkVox);
    sN2 = s1^2 + sigIni^2;  sB2 = s2^2 + sigIni^2;
    xLat   = dij.doseGrid.x(:) - iso(1);
    latM   = squeeze(doseCube(iyPk, :, izAx));  latM = latM(:);
    xFine  = linspace(-60, 60, 1201)';
    core   = (1-w)/(2*pi*sN2) * exp(-xFine.^2/(2*sN2));
    halo   =   w  /(2*pi*sB2) * exp(-xFine.^2/(2*sB2));
    model0 = (1-w)/(2*pi*sN2) + w/(2*pi*sB2);
    latM   = latM / latM(abs(xLat) == min(abs(xLat)));     % normalise at the axis
    haloZone = abs(xLat) >= 15 & abs(xLat) <= 30;
    modelAtX = ((1-w)/(2*pi*sN2)*exp(-xLat.^2/(2*sN2)) + w/(2*pi*sB2)*exp(-xLat.^2/(2*sB2))) / model0;
    haloRatio = mean(latM(haloZone) ./ modelAtX(haloZone));
    fprintf('Halo check (|x| = 15-30 mm): matRad / model = %.2f  (1 = halo delivered, 0 = cut off)\n', haloRatio);
    fprintf('  model at Bragg peak: sigma core %.2f mm, sigma halo %.1f mm, halo weight %.2f\n', ...
        sqrt(sN2), sqrt(sB2), w);

    summary(k,:) = {string(labels{k}), eUse, braggDepth, inZ(inMaxI), min(ratio), max(ratio), ...
        haloRatio, 100*idd(find(iddZ >= 1, 1))/iddMax, axSpot(iyPk)/axSpot(iyEnt), axField(iyPk)/axField(iyEnt)};

    %% ---- pack what the plotting functions need
    S = struct('label',labels{k}, 'E',eUse, 'machine',dwaMachine, 'isInterp',isInterp, ...
        'iddZ',iddZ, 'idd',idd/iddMax, 'inZ',inZ, 'inD',inD/inMax, 'caxZ',caxZ, 'cax',cax/max(cax), ...
        'braggDepth',braggDepth, 'isoDepth',isoDepth, ...
        'iyPk',iyPk, 'ixAx',ixAx, 'izAx',izAx, 'iyWater',iyWater, 'zPkVox',zPkVox, ...
        'bevMax',bevMax, 'cubeMax',cubeMax, 'surfaceY',surfaceY, 'iso',iso, ...
        'xLat',xLat, 'latM',latM, 'xFine',xFine, 'core',core/model0, 'halo',halo/model0, ...
        'haloRatio',haloRatio);

    %% ---- individual figures, then the same panels in the overview row
    f = figure('Color','w','Position',[100 100 1000 420]);
    drawDepthCurves(gca, S, true);
    saveFigure(f, fullfile(outDir, [tag '_1_depthDose']));  close(f);

    f = figure('Color','w');
    drawBEV(gca, ct, cst, doseCube, dij, S);
    saveFigure(f, fullfile(outDir, [tag '_2_beamseye']));  close(f);

    f = figure('Color','w');
    drawDepthView(gca, ct, cst, doseCube, dij, S);
    saveFigure(f, fullfile(outDir, [tag '_3_depthView']));  close(f);

    f = figure('Color','w','Position',[100 100 700 420]);
    drawLateral(gca, S, true);
    saveFigure(f, fullfile(outDir, [tag '_4_lateralLog']));  close(f);

    S5 = S;
    S5.cubeMax   = 1.001 * max(fieldCube(:));
    S5.viewTitle = sprintf('%s %.2f MeV: %dx%d spots, 3 mm apart', S.label, S.E, fieldN, fieldN);
    f = figure('Color','w');
    drawDepthView(gca, ct, cst, fieldCube, dij, S5);
    saveFigure(f, fullfile(outDir, [tag '_5_fieldDepthView']));  close(f);

    drawDepthCurves(nexttile(tl), S, k == 1);
    drawBEV(nexttile(tl), ct, cst, doseCube, dij, S);
    drawDepthView(nexttile(tl), ct, cst, doseCube, dij, S);
    drawLateral(nexttile(tl), S, k == 1);
    drawDepthView(nexttile(tl), ct, cst, fieldCube, dij, S5);
end

saveFigure(hComp, fullfile(outDir, [dwaMachine '_MinMidMax_overview']));
writetable(summary, fullfile(outDir, [dwaMachine '_MinMidMax_summary.csv']));
disp(summary);
fprintf('\nAll output in %s\n', outDir);

%% ========================================================================
%% Plotting functions - each draws into the axes it is given
%% ========================================================================
function drawDepthCurves(ax, S, showLegend)
    axes(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    h1 = plot(ax, S.iddZ, 100*S.idd, '.-', 'LineWidth', 1.3, 'MarkerSize', 7, ...
        'DisplayName', 'matRad, laterally integrated (1 mm grid)');
    h2 = plot(ax, S.inZ, 100*S.inD, 'k-', 'LineWidth', 1, 'DisplayName', ...
        ternary(S.isInterp, 'machine input (interpolated layer)', 'machine input (TOPAS)'));
    h3 = plot(ax, S.caxZ, 100*S.cax, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, ...
        'DisplayName', 'matRad, central axis (own max = 100)');
    xMax = min(280, S.braggDepth + max(15, 0.3*S.braggDepth));
    xlim(ax, [0 xMax]); ylim(ax, [0 105]);
    xlabel(ax, 'Depth in water [mm]'); ylabel(ax, 'Relative dose [%]');
    title(ax, sprintf('%s %.2f MeV: depth dose (peak %.1f mm)', S.label, S.E, S.braggDepth));
    if showLegend, legend(ax, [h1 h2 h3], 'Location', 'east', 'FontSize', 8); end
    set(ax, 'TickDir','out', 'LineWidth',0.9, 'FontSize',10);
end

function drawBEV(ax, ct, cst, doseCube, dij, S)
    m = 6;   % +-18 mm
    matRad_plotSlice(ct, 'axesHandle', ax, 'cst', cst, 'cubeIdx', 1, 'dose', doseCube, ...
        'plane', 1, 'slice', S.iyPk, 'alpha', 0.75, 'doseWindow', [0 S.bevMax], ...
        'voiSelection', zeros(size(cst,1),1), 'colorBarLabel', 'Dose [Gy]');
    axis(ax, 'on');
    xlim(ax, S.izAx + [-m m]);  ylim(ax, S.ixAx + [-m m]);
    xlabel(ax, 'z - z_{iso} [mm]', 'Color','k');  ylabel(ax, 'x - x_{iso} [mm]', 'Color','k');
    setNiceMmTicks(ax, 'x', dij.doseGrid.z, S.izAx, m, 5, S.iso(3));
    setNiceMmTicks(ax, 'y', dij.doseGrid.x, S.ixAx, m, 5, S.iso(1));
    styleSliceAxes(ax);
    title(ax, sprintf('%s %.2f MeV: beam''s-eye at %.0f mm', S.label, S.E, S.zPkVox), 'Color','k');
end

function drawDepthView(ax, ct, cst, doseCube, dij, S) %#ok<INUSL>
% Dose in the plane containing the beam axis, drawn with DEPTH ON THE
% HORIZONTAL AXIS (0 = water surface, increasing to the right), so it reads
% the same way as the depth-dose plot. matRad_plotSlice draws slices as
% images - vertical axis increasing DOWNWARD - which made the earlier
% version look inverted; this one is drawn directly with imagesc.
% Colour = dose per voxel. On the beam axis that is highest near the
% ENTRANCE for this narrow beam (see the grey central-axis curve); the
% Bragg peak dominates only in the laterally summed dose.
    depth = dij.doseGrid.y(:) - S.surfaceY;        % voxel-centre depth [mm]
    lat   = dij.doseGrid.z(:) - S.iso(3);          % lateral offset [mm]
    sl    = squeeze(doseCube(:, S.ixAx, :));       % [depth, lateral]
    halfW = min(60, max(20, 0.25*(S.braggDepth + 25)));
    inD   = depth >= -3 & depth <= S.braggDepth + 25;
    inL   = abs(lat) <= halfW;
    imagesc(ax, depth(inD), lat(inL), sl(inD, inL).');
    set(ax, 'YDir', 'normal');
    axis(ax, 'image');
    colormap(ax, jet(256));  caxis(ax, [0 S.cubeMax]);
    cb = colorbar(ax);  cb.Label.String = 'Dose [Gy]';
    hold(ax, 'on');
    xline(ax, 0, 'w-', 'LineWidth', 1);                              % water surface
    xline(ax, S.braggDepth, 'w--', 'LineWidth', 0.8);               % Bragg peak (summed dose)
    xlabel(ax, 'Depth in water [mm]');  ylabel(ax, 'z - z_{iso} [mm]');
    set(ax, 'TickDir','out', 'Layer','top', 'Box','on', 'FontSize',10);
    if isfield(S, 'viewTitle')
        title(ax, S.viewTitle);
    else
        title(ax, sprintf('%s %.2f MeV: single spot, depth view', S.label, S.E));
    end
end

function drawLateral(ax, S, showLegend)
    axes(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    h1 = semilogy(ax, S.xLat, max(S.latM, 1e-7), 'o', 'MarkerSize', 4, 'DisplayName', 'matRad (3 mm grid)');
    h2 = semilogy(ax, S.xFine, S.core + S.halo, 'k-', 'LineWidth', 1, 'DisplayName', 'machine model, total');
    h3 = semilogy(ax, S.xFine, S.core, ':', 'Color', [0 0.45 0.74], 'LineWidth', 1.2, 'DisplayName', 'core');
    h4 = semilogy(ax, S.xFine, S.halo, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, 'DisplayName', 'halo');
    set(ax, 'YScale', 'log');
    xlim(ax, [-60 60]); ylim(ax, [1e-4 1.5]);
    xlabel(ax, 'x - x_{iso} [mm]'); ylabel(ax, 'relative dose');
    title(ax, sprintf('%s %.2f MeV: lateral at %.0f mm (halo ratio %.2f)', ...
        S.label, S.E, S.zPkVox, S.haloRatio));
    if showLegend, legend(ax, [h1 h2 h3 h4], 'Location', 'south', 'FontSize', 8); end
    set(ax, 'TickDir','out', 'LineWidth',0.9, 'FontSize',10);
end

%% ========================================================================
%% Helpers (same as in the single-beam script)
%% ========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end


%% Local function: publication-quality figure export
function saveFigure(figHandle, baseName)
% saveas(...,'png') rasterises at screen resolution (~96 dpi), which is why
% the exported PNGs look soft. exportgraphics re-renders at a chosen
% resolution and trims the surrounding whitespace. A vector PDF is written
% alongside it - that one stays sharp at any zoom and is what to use in a
% document or a slide. The .fig is kept so the figure can be reopened and
% edited in MATLAB later.
    set(figHandle, 'Color', 'w', 'InvertHardcopy', 'off');
    if ~isempty(which('exportgraphics'))   % R2020a and later
        exportgraphics(figHandle, [baseName '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
        exportgraphics(figHandle, [baseName '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    else
        print(figHandle, [baseName '.png'], '-dpng', '-r300');
        print(figHandle, [baseName '.pdf'], '-dpdf', '-painters');
    end
    savefig(figHandle, [baseName '.fig']);
    fprintf('Saved %s .png (300 dpi) / .pdf (vector) / .fig\n', baseName);
end

%% Local function: make axis ticks visible on top of an image
function styleSliceAxes(ax)
% matRad_plotSlice draws the CT/dose as an image filling the axes. Default
% MATLAB ticks point INWARD, so they are painted over by that image and
% disappear - which is why the slice figures showed tick labels but no tick
% marks. 'Layer','top' draws the axis rulers above the image and TickDir
% 'out' moves the ticks outside the image entirely.
    set(ax, 'Layer','top', 'TickDir','out', 'TickLength',[0.012 0.012], ...
            'LineWidth',0.9, 'Box','on', 'FontSize',11, 'XColor','k', 'YColor','k');
end

%% Local function: clean, evenly-spaced axis ticks over a zoomed window
function setNiceMmTicks(axHandle, whichAxis, gridVecMm, centerIx, marginVox, nTicks, originMm)
% Places nTicks evenly-spaced ticks (in voxel-index units, since
% matRad_plotSlice draws in index coordinates) across [centerIx-marginVox,
% centerIx+marginVox], snapped so the LABELS are round numbers of mm.
%
% originMm (optional) shifts the labels to a chosen reference point:
% label = gridVecMm - originMm. This matters because matRad's world
% origin is NOT the centre of the phantom. matRad_getWorldAxes builds the
% axes as firstVox = -(dim/2)*res, i.e. it puts the first voxel CENTRE at
% exactly -N/2*res, so for BOXPHANTOM x runs -240:3:237 - an axis whose
% midpoint is -1.5 mm, not 0. On top of that the target spans an odd
% number of voxels (70..90), so its centroid - which is what
% matRad_getIsoCenter returns - lands at -3 mm in every axis. Labelling
% against the raw world coordinate therefore puts "0" 3 mm off the beam
% axis. Passing originMm = stf(1).isoCenter(k) makes 0 mean the beam axis;
% passing the entrance-surface coordinate makes the axis read as depth.
    if nargin < 7 || isempty(originMm)
        originMm = 0;
    end

    nVox = numel(gridVecMm);
    lo = max(1, centerIx - marginVox);
    hi = min(nVox, centerIx + marginVox);

    labelVals = gridVecMm(:).' - originMm;     % what the ticks will SAY
    spanLo = min(labelVals(lo), labelVals(hi));
    spanHi = max(labelVals(lo), labelVals(hi));

    % Choose a round step that gives about nTicks ticks across the visible
    % span, then place ticks only at round labels that actually fall INSIDE
    % that span - otherwise a rounded label can land just outside the zoom
    % window and the tick silently disappears.
    niceSteps = [1 2 5 10 15 20 25 30 40 50 60 75 100 150 200 250 500];
    rawStep   = (spanHi - spanLo) / max(nTicks - 1, 1);
    step      = niceSteps(find(niceSteps >= rawStep, 1));
    if isempty(step), step = rawStep; end

    wantLabels = ceil(spanLo/step)*step : step : floor(spanHi/step)*step;
    if numel(wantLabels) < 2
        wantLabels = [spanLo spanHi];
    end

    % snap each desired label to the nearest voxel inside the window
    tickIx = zeros(size(wantLabels));
    for k = 1:numel(wantLabels)
        [~, rel] = min(abs(labelVals(lo:hi) - wantLabels(k)));
        tickIx(k) = lo + rel - 1;
    end
    [tickIx, keep] = unique(tickIx, 'stable');
    tickLabels = arrayfun(@(v) sprintf('%g', v), wantLabels(keep), 'UniformOutput', false);

    if strcmpi(whichAxis, 'x')
        set(axHandle, 'XTick', tickIx, 'XTickLabel', tickLabels);
    else
        set(axHandle, 'YTick', tickIx, 'YTickLabel', tickLabels);
    end
end
