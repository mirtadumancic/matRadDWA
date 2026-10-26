%% Water Phantom Validation: single beam, single energy, DWA_proton vs BOXPHANTOM
%
% Adapted from matRad_example5_protons_DWA_BOXPHANTOM.m for Step 1 of the
% DWA machine validation plan (per Remo Cristoforetti, DKFZ): fire ONE
% pencil beam at ONE energy into matRad's BOXPHANTOM using the DWA_proton
% machine, and plot the resulting percentage/central-axis depth-dose
% (PDD) curve.
%
% BOXPHANTOM.mat geometry (checked directly against the .mat file):
%   - 160^3 voxel cube, 3 mm isotropic resolution -> 480x480x480 mm cube
%   - Water block ('BODY' OAR): 118.5-358.5 mm in x/y/z (240 mm thick),
%     surrounded by air.
%   - 'OuterTarget' (TARGET): 208.5-268.5 mm in x/y/z, a 60 mm cube
%     centered in the water block.
%   - matRad_getIsoCenter places the isocenter at the target's bounding-box
%     center, i.e. (238.5, 238.5, 238.5) mm.
%   => entrance-surface-to-isocenter depth ~120 mm; entrance-to-distal-
%      target-edge depth ~150 mm; full water thickness 240 mm.
% testEnergy_MeV = 150 is kept below because 150 MeV protons range out to
% ~158 mm in water, i.e. right around the distal target depth - a sensible
% physical pick for this phantom, not an arbitrary placeholder.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set matRad runtime configuration
cd '/Users/mirtadumancic/Work/MATLAB_Projects/matRadDWA';
run('matRad_rc');

%% Load BOXPHANTOM instead of a patient dataset
load('BOXPHANTOM.mat');

%% Treatment Plan
% radiationMode/machine: use the commissioned DWA machine
pln.radiationMode = 'protons';
pln.machine       = 'DWA_proton';
pln.bioModel      = 'constRBE';
pln.multScen      = 'nomScen';

pln.propDoseCalc.calcLET = 0;
pln.propDoseCalc.engine  = 'HongPB';

% single fraction is enough for a physics/PDD check
pln.numOfFractions = 1;

% ONE beam, straight into the phantom (gantry = 0, couch = 0)
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles  = 0;
pln.propStf.bixelWidth   = 5;
pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO       = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings - keep the CT resolution for a clean profile
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y;
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;

pln.propOpt.quantityOpt = 'physicalDose';   % raw physics PDD, no RBE weighting

%% Generate Beam Geometry STF: ONE ray, ONE energy, by construction
% This matRad build auto-falls back to the OOP 'ParticleSingleSpot' stf
% generator (matRad_StfGeneratorParticleSingleBeamlet) whenever the
% default generator isn't available - and that generator already does
% exactly what Step 1 needs: a single ray at BEV [0 0 0], single energy.
% We select it explicitly (rather than relying on the fallback) and give
% it the energy directly via pln.propStf.energy. This also sidesteps a
% bug in that generator's automatic energy-selection path (a ray-tracing
% calculation through a synthetic single-voxel target) that otherwise
% throws "Arrays have incompatible sizes for this operation" in
% setBeamletEnergies.
%
% 150 MeV is a sensible pick for BOXPHANTOM: its water range (~158 mm)
% lands just past the target's distal edge (~150 mm deep), so the Bragg
% peak should sit inside/just beyond the OuterTarget volume rather than
% stopping short in the entrance region or overshooting out the back of
% the water block (which ends at 240 mm depth).
testEnergy_MeV = 150;

pln.propStf.generator = 'ParticleSingleSpot';
pln.propStf.energy    = testEnergy_MeV;

stf = matRad_generateStf(ct,cst,pln);

% the generator snaps testEnergy_MeV to the closest energy layer actually
% available in the DWA_proton machine file - report what was actually used
fprintf('Using energy %.2f MeV (closest available to requested %.1f MeV)\n', ...
    stf(1).ray.energy, testEnergy_MeV);

%% Dose Calculation
% Single bixel -> unit weight gives the dose of this one pencil beam directly,
% no optimization needed.
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
doseCube = reshape(full(dij.physicalDose{1}), dij.doseGrid.dimensions);

%% Extract and plot the PDD (central-axis depth-dose curve)
% Sample along the true beam axis (source -> isocenter), so this works
% regardless of gantry/couch angle, rather than assuming a fixed cube axis.
sourcePoint = stf(1).sourcePoint;
targetPoint = stf(1).isoCenter;
beamDir     = (targetPoint - sourcePoint) / norm(targetPoint - sourcePoint);

% IMPORTANT: sourcePoint sits SAD=10000mm upstream of isocenter (a fixed
% numerical convention for the pencil-beam formalism, not a real physical
% distance - see the machine's meta.description). Sampling depths_mm
% directly from sourcePoint would stay in air the whole time and never
% reach the phantom. Instead, start from the water ENTRANCE surface,
% which sits ~120mm before isocenter along the beam (the same value used
% for the "Isocenter (120mm)" reference line below).
isoDepth_mm = 120;   % entrance-surface-to-isocenter depth, from the BOXPHANTOM geometry
entryPoint  = targetPoint - isoDepth_mm * beamDir;

% water block is 240 mm thick along the beam axis; sample a bit past it
% (280 mm) so the falloff into the exit-side air is visible too
depths_mm = (0:dij.doseGrid.resolution.x/2:280)';
coords    = entryPoint + depths_mm * beamDir;   % Nx3, world [x y z] mm, depths_mm=0 at the entrance surface

% Use matRad's OWN world-coordinate vectors for the dose grid rather than
% reconstructing them by hand: matRad's internal coordinate convention
% isn't simply "distance from the cube corner" (it's offset/centered
% differently), and stf(1).sourcePoint/isoCenter are given in that same
% internal frame. Reading dij.doseGrid.x/.y/.z directly guarantees the
% interpolation grid lines up exactly with where sourcePoint/isoCenter
% - and therefore doseCube - actually live, instead of silently sampling
% the wrong region the way a hand-built (i-0.5)*res grid can.
% IMPORTANT: matRad stores cubes as [y, x, z] (rows = y, columns = x),
% see ct.cubeDim = [Ny Nx Nz]. The grid vectors must therefore be passed
% in that order, and the query points swapped to (y, x, z) to match.
% Passing {x, y, z} silently transposes the cube: in this phantom x and y
% grids are identical, so nothing errors - but the 'PDD' line then runs
% laterally through the spot and returns a lateral profile instead.
F = griddedInterpolant({dij.doseGrid.y, dij.doseGrid.x, dij.doseGrid.z}, doseCube, 'linear', 'none');
pdd = F(coords(:,2), coords(:,1), coords(:,3));
pdd(isnan(pdd)) = 0;   % outside the dose grid -> treat as zero dose

[peakVal, peakIx] = max(pdd);
fprintf('Bragg peak depth: %.1f mm (%.2f%% of max reached at entrance)\n', ...
    depths_mm(peakIx), 100*pdd(1)/peakVal);

% Sanity check against the known BOXPHANTOM geometry: water block is
% 0-240 mm deep from the entrance surface, target sits at 90-150 mm deep.
if depths_mm(peakIx) > 240
    warning('Bragg peak (%.1f mm) falls beyond the water block (240 mm) - overshooting into air. Pick a lower testEnergy_MeV.', depths_mm(peakIx));
elseif depths_mm(peakIx) < 90
    warning('Bragg peak (%.1f mm) falls short of the target (90-150 mm deep) - pick a higher testEnergy_MeV.', depths_mm(peakIx));
end

figure('Color','w');
plot(depths_mm, 100*pdd/peakVal, 'LineWidth', 2);
grid on; box on; hold on;
xline(120, '--', 'Isocenter (120 mm)', 'LabelVerticalAlignment','bottom');
xline([90 150], ':', {'Target proximal','Target distal'});
xlabel('Depth in water [mm]');
ylabel('Relative dose [%]');
title(sprintf('DWA\\_proton PDD @ %.1f MeV, single pencil beam, BOXPHANTOM', stf(1).ray.energy));
xlim([0 depths_mm(end)]);
set(gca, 'TickDir','out', 'TickLength',[0.012 0.012], 'LineWidth',0.9, 'FontSize',11);
saveFigure(gcf, 'DWA_proton_PDD');

%% 2D colorwash dose slices, anchored on the TRUE dose maximum voxel
% matRad cubes are indexed [y, x, z]. matRad_plotSlice: plane 1 fixes the
% 1st index (y), plane 2 fixes the 2nd (x), plane 3 fixes the 3rd (z).
% The beam here travels along +y (gantry 0), so plane 1 = beam's-eye view
% and plane 2 = depth-dose view along the beam axis.
% Slicing at the nominal isocenter voxel can miss the true 3D dose peak
% by a voxel or two (the beam isn't perfectly axis-aligned), which shows
% up as the colorwash never reaching the top of its own color scale. So
% find the actual peak voxel first and anchor both slices on it.
[peakDoseVal, peakLinIx] = max(doseCube(:));
[iyPeak, ixPeak, izPeak] = ind2sub(size(doseCube), peakLinIx);
fprintf('Peak dose voxel (array indices): [%d %d %d], value=%.4g\n', iyPeak, ixPeak, izPeak, peakDoseVal);

% Zoom window: matRad_plotSlice draws in VOXEL-INDEX axis coordinates
% (the mm tick labels are just relabeled text on top of that), so the
% zoom has to be applied as an index range, not a mm range. +-20 voxels
% (=60mm at 3mm resolution) around the peak comfortably frames the 60mm
% target box plus a little margin, cropping out most of the surrounding
% air/water background.
marginVox = 20;

% Give OuterTarget a distinct, fixed contour color instead of relying on
% matRad_plotSlice's default contourColorMap sampling. matRad_plotVoiContourSlice
% falls back to each structure's own cst{i,5}.visibleColor when no
% 'contourColorMap' is passed, so set that directly, looked up by name
% (not a hardcoded row index) so this still works if cst's row order changes.
%
% BODY is deliberately left OUT of the contour entirely (via voiSelection
% below), not just given a color: BODY is the full water block, and the
% diagnostic check confirmed its boundary sits far outside this +-60 voxel
% zoom window (the whole zoomed view is inside BODY, nowhere near its
% edge) - so a BODY contour line has nothing to trace here and would
% either be invisible (correctly) or misleadingly implied a boundary
% exists nearby. Only OuterTarget's boundary is actually within view.
bodyIx   = find(strcmpi(cst(:,2), 'BODY'), 1);
targetIx = find(strcmpi(cst(:,2), 'OuterTarget'), 1);
if ~isempty(targetIx), cst{targetIx,5}.visibleColor = [0.10 0.30 0.85]; end   % blue

% NOTE: matRad_plotSlice documents 'voiSelection' as "logicals", but its
% own input validator is @(x) isnumeric(x) || isempty(x) - and in MATLAB a
% logical array is NOT isnumeric, so passing true(N,1) is rejected with
% "The value of 'voiSelection' is invalid". Use a double 0/1 vector: it
% passes the validator, and matRad_plotVoiContourSlice only ever uses it
% in a logical context (if ... && selection(s)), where 0/1 works fine.
voiSelection = ones(size(cst,1), 1);
if ~isempty(bodyIx), voiSelection(bodyIx) = 0; end

% matRad_plotVoiContourSlice only appends a cell to its output for
% non-'IGNORED', selected cst rows, in row order - so the position of a
% structure's handle in hContour is its rank among those rows, not its
% raw cst row index. Compute that mapping once so the legend below points
% at the right handle.
eligibleRows = find(~strcmpi(cst(:,3), 'IGNORED') & logical(voiSelection));
targetPos    = find(eligibleRows == targetIx, 1);

% Beam's-eye view (plane=1, fixes y=iyPeak, the Bragg peak depth):
% squeeze(doseCube(iyPeak,:,:)) -> rows = x, columns = z. Shows the
% round spot cross-section at its hottest depth.
figure('Color','w');
[~,~,~,hContourCor,~] = matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, 'dose', doseCube, ...
    'plane', 1, 'slice', iyPeak, 'alpha', 0.75, 'doseWindow', [0 peakDoseVal], ...
    'voiSelection', voiSelection, 'colorBarLabel', 'Dose [Gy] (unit pencil-beam weight)');
title(sprintf('DWA\\_proton dose, beam''s-eye view at Bragg peak depth @ %.1f MeV', stf(1).ray.energy));
xlim(gca, izPeak + [-marginVox, marginVox]);
ylim(gca, ixPeak + [-marginVox, marginVox]);
% matRad_plotSlice sets axis label text via its own GUI theme color,
% which can end up invisible against our white figure background -
% re-apply the labels explicitly in plain black.
% matRad_plotSlice internally calls axis(axesHandle,'off') to hide the CT
% slice's box, and never turns it back on - which hides the axis box,
% ticks AND label text no matter what we set them to. Switch it back on
% before (re)applying labels.
axis(gca, 'on');
xlabel(gca, 'z [mm]', 'Color', 'k');
ylabel(gca, 'x [mm]', 'Color', 'k');   % rows of plane-1 image = x (lateral)
set(gca, 'XColor', 'k', 'YColor', 'k');
% Replace matRad's default 10 ticks spanning the whole (unzoomed) grid
% with a handful of clean, evenly-spaced ticks confined to the zoomed
% window actually shown, labeled with their real mm positions.
setNiceMmTicks(gca, 'x', dij.doseGrid.z, izPeak, marginVox, 5);
setNiceMmTicks(gca, 'y', dij.doseGrid.x, ixPeak, marginVox, 5);
styleSliceAxes(gca);
% A structure's contour cell can come back empty ({}) if it happens not
% to intersect this particular slice - guard against that before indexing.
if ~isempty(targetPos) && ~isempty(hContourCor{targetPos})
    legend(hContourCor{targetPos}(1), {'OuterTarget'}, 'TextColor', 'k', 'Location', 'northeast');
end
saveFigure(gcf, 'DWA_proton_beamseye_slice');

% Depth view (plane=2, fixes x=ixPeak, through the beam axis):
% squeeze(doseCube(:,ixPeak,:)) -> rows = y (depth), columns = z.
% Shows the entrance plateau, the Bragg peak and the lateral spread.
figure('Color','w');
[~,~,~,hContourSag,~] = matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, 'dose', doseCube, ...
    'plane', 2, 'slice', ixPeak, 'alpha', 0.75, 'doseWindow', [0 peakDoseVal], ...
    'voiSelection', voiSelection, 'colorBarLabel', 'Dose [Gy] (unit pencil-beam weight)');
title(sprintf('DWA\\_proton dose, depth view through beam axis @ %.1f MeV', stf(1).ray.energy));
xlim(gca, izPeak + [-marginVox, marginVox]);
ylim(gca, iyPeak + [-marginVox, marginVox]);
axis(gca, 'on');
xlabel(gca, 'z [mm]', 'Color', 'k');
ylabel(gca, 'y [mm] (beam direction)', 'Color', 'k');   % rows of plane-2 image = y (depth)
set(gca, 'XColor', 'k', 'YColor', 'k');
setNiceMmTicks(gca, 'x', dij.doseGrid.z, izPeak, marginVox, 5);
setNiceMmTicks(gca, 'y', dij.doseGrid.y, iyPeak, marginVox, 5);
styleSliceAxes(gca);
if ~isempty(targetPos) && ~isempty(hContourSag{targetPos})
    legend(hContourSag{targetPos}(1), {'OuterTarget'}, 'TextColor', 'k', 'Location', 'northeast');
end
saveFigure(gcf, 'DWA_proton_depth_slice');

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
function setNiceMmTicks(axHandle, whichAxis, gridVecMm, centerIx, marginVox, nTicks)
% Places nTicks evenly-spaced ticks (in voxel-index units, since
% matRad_plotSlice draws in index coordinates) across [centerIx-marginVox,
% centerIx+marginVox], but snaps each one to the nearest voxel whose real
% mm coordinate (from gridVecMm, e.g. dij.doseGrid.x/.y/.z) is a round
% number of mm - so the labels read like "90, 105, 120, 135, 150" instead
% of whatever odd values the raw grid resolution would otherwise produce.
    nVox = numel(gridVecMm);
    idxWindow = [max(1, centerIx - marginVox), min(nVox, centerIx + marginVox)];
    targetIx  = round(linspace(idxWindow(1), idxWindow(2), nTicks));

    mmAtTarget = gridVecMm(targetIx);
    roundTo    = 10; % mm
    mmNice     = round(mmAtTarget / roundTo) * roundTo;

    tickIx = zeros(size(mmNice));
    for k = 1:numel(mmNice)
        [~, tickIx(k)] = min(abs(gridVecMm - mmNice(k)));
    end
    tickIx = unique(tickIx, 'stable');
    tickLabels = arrayfun(@(v) sprintf('%g', gridVecMm(v)), tickIx, 'UniformOutput', false);

    if strcmpi(whichAxis, 'x')
        set(axHandle, 'XTick', tickIx, 'XTickLabel', tickLabels);
    else
        set(axHandle, 'YTick', tickIx, 'YTickLabel', tickLabels);
    end
end