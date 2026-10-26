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
%% Which DWA machine to plan with.
% 'DWA_proton_doubleGauss_interp' - NEWEST: double Gaussian + 201 energy layers
%                                   (58 native + 143 interpolated), from
%                                   matRad_commissionDWAmachine_doubleGauss_interp.m
% 'DWA_proton_doubleGauss'        - double Gaussian, 58 native layers
% 'DWA_proton'                    - OLDEST: single Gaussian (halo discarded,
%                                   spot oversized ~4x), 58 native layers
% All three sit in basedata. Every output folder and file name below carries
% this machine name, so runs with different machines never overwrite each other.
dwaMachine = 'DWA_proton_doubleGauss_interp';

% Output folder and file-name prefix, both tagged with the machine name
outDir = ['singleBeam_' dwaMachine];
if ~exist(outDir,'dir'), mkdir(outDir); end

pln.radiationMode = 'protons';
pln.machine       = dwaMachine;
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
%
% (a) On the CT grid (3 mm). Used for the colour-wash slices, because
%     matRad_plotSlice draws a dose cube on the CT grid.
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
doseCube = reshape(full(dij.physicalDose{1}), dij.doseGrid.dimensions);

% (b) On a finer grid, for the depth curves only. The 3 mm grid is too
%     coarse for this beam in both directions:
%       - along the beam, a 3 mm step cannot land on a Bragg peak that is
%         only ~4.5 mm wide (80-80), so the sampled maximum is lower than the
%         true one and every other point, normalised to it, reads too high;
%       - across the beam, the ~1.2 mm core is narrower than a voxel, so
%         summing 3 mm voxels over-counts the entrance dose by 6-13%.
%     1 mm along the beam and 2 mm across fixes both (the across-beam error
%     drops below 0.1%). Only the non-zero voxels of the sparse dij are read,
%     so the finer grid does not need a dense cube in memory.
fineRes = struct('x', 2, 'y', 1, 'z', 2);     % beam travels along y (gantry 0)
plnFine = pln;
plnFine.propDoseCalc.doseGrid.resolution = fineRes;
dijFine = matRad_calcDoseInfluence(ct,cst,stf,plnFine);

%% Geometry: beam direction and the WATER SURFACE, found from the CT
sourcePoint = stf(1).sourcePoint;
targetPoint = stf(1).isoCenter;
% stf.sourcePoint is stored RELATIVE to the isocentre, so the beam direction
% is simply minus it (isoCenter - sourcePoint mixes the two frames; it only
% worked because SAD = 10 m dwarfs the 3 mm isocentre offset)
beamDir     = -sourcePoint / norm(sourcePoint);
if abs(beamDir(2) - 1) > 1e-3
    error('This script assumes gantry 0 / couch 0 (beam along +y).');
end

% The surface used to be ASSUMED 120 mm upstream of isocentre. It is now read
% from the CT: the first voxel along the beam axis that contains water. On
% the 3 mm CT grid the assumption could be a voxel off, which shifts every
% depth by up to 3 mm. (The coarse dose grid is the CT grid here.)
[~, ixIsoCT] = min(abs(dij.doseGrid.x - targetPoint(1)));
[~, izIsoCT] = min(abs(dij.doseGrid.z - targetPoint(3)));
if isfield(ct,'cube') && ~isempty(ct.cube)
    isWater = ct.cube{1}(:, ixIsoCT, izIsoCT) > 0.5;       % relative stopping power
else
    isWater = ct.cubeHU{1}(:, ixIsoCT, izIsoCT) > -500;    % Hounsfield units
end
iyWater     = find(isWater, 1, 'first');
surfaceY    = dij.doseGrid.y(iyWater) - ct.resolution.y/2;   % upstream face of that voxel
isoDepth_mm = targetPoint(2) - surfaceY;
entryPoint  = [targetPoint(1) surfaceY targetPoint(3)];
fprintf('Water surface at y = %.1f mm: isocentre is %.1f mm deep\n', surfaceY, isoDepth_mm);

%% Depth curves
% --- Laterally INTEGRATED depth dose (IDD), fine grid. This, not the
% central-axis curve, is what the TOPAS Bragg curves are: dose summed over
% the whole lateral plane at each depth. For a wide spot the two look alike.
% For the DWA's ~1.2 mm core they do not: the pencil beam widens with depth
% (about 1.2 -> 4 mm by the peak), so the on-axis dose falls roughly as
% 1/sigma^2 and the Bragg peak all but disappears from the central-axis
% curve. That is the physics of a narrow beam, not a machine error.
dimsF = dijFine.doseGrid.dimensions;                     % [Ny Nx Nz]
[vIdx, ~, vDose] = find(dijFine.physicalDose{1}(:,1));
[iyF, ~, ~] = ind2sub(dimsF, vIdx);
voxAreaF    = dijFine.doseGrid.resolution.x * dijFine.doseGrid.resolution.z;
iddFine     = accumarray(iyF, vDose, [dimsF(1) 1]) * voxAreaF;   % Gy*mm^2 per unit weight
iddDepth_mm = dijFine.doseGrid.y(:) - surfaceY;                   % voxel centre depth
inWater     = iddDepth_mm >= 0 & iddDepth_mm <= 280;
iddDepth_mm = iddDepth_mm(inWater);
iddFine     = iddFine(inWater);

% --- Central axis, coarse grid. The beam axis sits exactly on a 3 mm voxel
% centre but between 2 mm voxels, so the coarse grid is the right one here.
% matRad stores cubes as [y, x, z]; the interpolant takes {y, x, z}.
depths_mm = (0:0.5:280)';
coords    = entryPoint + depths_mm * beamDir;
F   = griddedInterpolant({dij.doseGrid.y, dij.doseGrid.x, dij.doseGrid.z}, doseCube, 'linear', 'none');
pdd = F(coords(:,2), coords(:,1), coords(:,3));
pdd(isnan(pdd)) = 0;

% --- the machine's own depth curve for this energy (TOPAS, or interpolated)
machine   = matRad_loadMachine(pln);
[~, eIx]  = min(abs([machine.data.energy] - stf(1).ray.energy));
inDepth   = machine.data(eIx).depths(:);
inZ       = machine.data(eIx).Z(:);
isInterp  = isfield(machine.data,'interpolated') && isequal(machine.data(eIx).interpolated, true);

[iddPeakVal, iddPeakIx] = max(iddFine);
[inPeakVal,  inPeakIx]  = max(inZ);
braggDepth_mm = iddDepth_mm(iddPeakIx);
fprintf('matRad IDD (fine grid): Bragg peak at %.1f mm, entrance %.1f%% of peak\n', ...
    braggDepth_mm, 100*iddFine(1)/iddPeakVal);
fprintf('Machine input curve:    Bragg peak at %.1f mm, entrance %.1f%% of peak%s\n', ...
    inDepth(inPeakIx), 100*inZ(1)/inPeakVal, ternary(isInterp,' (interpolated layer)',''));

% Shape check: matRad IDD vs input curve, both normalised to their peak.
% Should be ~1.00 everywhere. A ratio that drifts with depth would mean the
% lateral cutoff is clipping part of the halo.
inOnGrid = interp1(inDepth, inZ/inPeakVal, iddDepth_mm, 'linear', NaN);
okCmp    = iddDepth_mm < inDepth(inPeakIx) - 5 & isfinite(inOnGrid);
ratio    = iddFine(okCmp)/iddPeakVal ./ inOnGrid(okCmp);
fprintf('matRad IDD / input curve, entrance to 5 mm before the peak: %.3f - %.3f\n', min(ratio), max(ratio));

if braggDepth_mm > 240
    warning('Bragg peak (%.1f mm) falls beyond the water block (240 mm). Pick a lower testEnergy_MeV.', braggDepth_mm);
elseif braggDepth_mm < isoDepth_mm - 30
    warning('Bragg peak (%.1f mm) falls short of the target. Pick a higher testEnergy_MeV.', braggDepth_mm);
end

figure('Color','w','Position',[100 100 1000 420]);
hIDD = plot(iddDepth_mm, 100*iddFine/iddPeakVal, '.-', 'LineWidth', 1.3, 'MarkerSize', 7, ...
    'DisplayName', 'matRad, laterally integrated (IDD, 1 mm grid)');
grid on; box on; hold on;
hIn  = plot(inDepth, 100*inZ/inPeakVal, 'k-', 'LineWidth', 1, ...
    'DisplayName', ternary(isInterp, 'machine input curve (interpolated layer)', 'machine input curve (TOPAS)'));
hCax = plot(depths_mm, 100*pdd/max(pdd), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, ...
    'DisplayName', 'matRad, central axis only (own max = 100)');
xline(isoDepth_mm, '--', sprintf('Isocenter (%.1f mm)', isoDepth_mm), 'LabelVerticalAlignment','bottom');
xline(isoDepth_mm + [-30 30], ':', {'Target proximal','Target distal'});
xlabel('Depth in water [mm]');
ylabel('Relative dose [%]');
title(sprintf('%s depth dose @ %.2f MeV, single pencil beam, BOXPHANTOM', ...
    strrep(dwaMachine,'_','\_'), stf(1).ray.energy));
xlim([0 depths_mm(end)]); ylim([0 105]);
% the region beyond the Bragg peak is empty, so the legend goes there
legend([hIDD hIn hCax], 'Location', 'east');
set(gca, 'TickDir','out', 'TickLength',[0.012 0.012], 'LineWidth',0.9, 'FontSize',11);
saveFigure(gcf, fullfile(outDir, [dwaMachine '_PDD']));

%% 2D colorwash dose slices, anchored on the TRUE dose maximum voxel
% matRad cubes are indexed [y, x, z]. matRad_plotSlice: plane 1 fixes the
% 1st index (y), plane 2 fixes the 2nd (x), plane 3 fixes the 3rd (z).
% The beam here travels along +y (gantry 0), so plane 1 = beam's-eye view
% and plane 2 = depth-dose view along the beam axis.
% Slicing at the nominal isocenter voxel can miss the true 3D dose peak
% by a voxel or two (the beam isn't perfectly axis-aligned), which shows
% up as the colorwash never reaching the top of its own color scale. So
% find the actual peak voxel first and anchor both slices on it.
% For a narrow beam the hottest voxel sits near the ENTRANCE (see the PDD
% note above), so the beam's-eye slice is taken at the IDD Bragg-peak depth
% instead, and the lateral position from the hottest voxel in that slice.
[~, iyPeak] = min(abs(dij.doseGrid.y - (surfaceY + braggDepth_mm)));   % coarse-grid voxel at the Bragg peak
peakSlice   = squeeze(doseCube(iyPeak,:,:));
[~, lin2]   = max(peakSlice(:));
[ixPeak, izPeak] = ind2sub(size(peakSlice), lin2);
% Colour windows, one per figure, ending slightly ABOVE that figure's own
% maximum (a window ending exactly at the maximum can drop the hottest voxel
% out of the colormap - the white voxel in the earlier figures).
%   beam's-eye view: this slice's maximum, so the spot uses the full scale
%   depth view: the whole cube's maximum, which for this narrow beam is on
%   the axis near the entrance - that is where the red is
bevDoseMax   = 1.001 * max(peakSlice(:));
depthDoseMax = 1.001 * max(doseCube(:));
fprintf('Bragg-peak slice: voxel [%d %d %d], max %.4g Gy; cube max %.4g Gy\n', ...
    iyPeak, ixPeak, izPeak, bevDoseMax/1.001, depthDoseMax/1.001);

% Zoom window: matRad_plotSlice draws in VOXEL-INDEX axis coordinates
% (the mm tick labels are just relabeled text on top of that), so the
% zoom has to be applied as an index range, not a mm range. +-20 voxels
% (=60mm at 3mm resolution) around the peak comfortably frames the 60mm
% target box plus a little margin, cropping out most of the surrounding
% air/water background.
marginVox    = 20;   % depth view, across the beam (the target box is 60 mm wide)
marginBEV    = 6;    % beam's-eye view: +-18 mm, the spot is ~4 mm sigma at the peak

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
    'plane', 1, 'slice', iyPeak, 'alpha', 0.75, 'doseWindow', [0 bevDoseMax], ...
    'voiSelection', voiSelection, 'colorBarLabel', 'Dose [Gy] (unit pencil-beam weight)');
title(sprintf('%s, beam''s-eye view at Bragg-peak depth (%.0f mm) @ %.2f MeV', ...
    strrep(dwaMachine,'_','\_'), braggDepth_mm, stf(1).ray.energy));
xlim(gca, izPeak + [-marginBEV, marginBEV]);
ylim(gca, ixPeak + [-marginBEV, marginBEV]);
% matRad_plotSlice sets axis label text via its own GUI theme color,
% which can end up invisible against our white figure background -
% re-apply the labels explicitly in plain black.
% matRad_plotSlice internally calls axis(axesHandle,'off') to hide the CT
% slice's box, and never turns it back on - which hides the axis box,
% ticks AND label text no matter what we set them to. Switch it back on
% before (re)applying labels.
axis(gca, 'on');
xlabel(gca, 'z - z_{iso} [mm]', 'Color', 'k');
ylabel(gca, 'x - x_{iso} [mm]', 'Color', 'k');   % rows of plane-1 image = x (lateral)
set(gca, 'XColor', 'k', 'YColor', 'k');
% Replace matRad's default 10 ticks spanning the whole (unzoomed) grid
% with a handful of clean, evenly-spaced ticks confined to the zoomed
% window actually shown. Labels are referenced to the ISOCENTER, not to
% matRad's world origin - see the note at setNiceMmTicks below - so that
% 0 marks the beam axis.
setNiceMmTicks(gca, 'x', dij.doseGrid.z, izPeak, marginBEV, 5, stf(1).isoCenter(3));
setNiceMmTicks(gca, 'y', dij.doseGrid.x, ixPeak, marginBEV, 5, stf(1).isoCenter(1));
styleSliceAxes(gca);
% A structure's contour cell can come back empty ({}) if it happens not
% to intersect this particular slice - guard against that before indexing.
if ~isempty(targetPos) && ~isempty(hContourCor{targetPos})
    legend(hContourCor{targetPos}(1), {'OuterTarget'}, 'TextColor', 'k', ...
        'Location', 'northeastoutside', 'AutoUpdate', 'off');
end
saveFigure(gcf, fullfile(outDir, [dwaMachine '_beamseye_slice']));

% Depth view (plane=2, fixes x=ixPeak, through the beam axis):
% squeeze(doseCube(:,ixPeak,:)) -> rows = y (depth), columns = z.
% Shows the entrance plateau, the Bragg peak and the lateral spread.
figure('Color','w');
[~,~,~,hContourSag,~] = matRad_plotSlice(ct, 'axesHandle', gca, 'cst', cst, 'cubeIdx', 1, 'dose', doseCube, ...
    'plane', 2, 'slice', ixPeak, 'alpha', 0.75, 'doseWindow', [0 depthDoseMax], ...
    'voiSelection', voiSelection, 'colorBarLabel', 'Dose [Gy] (unit pencil-beam weight)');
title(sprintf('%s, depth view through beam axis @ %.2f MeV', ...
    strrep(dwaMachine,'_','\_'), stf(1).ray.energy));
% show the whole path: from just before the water surface to ~20 mm past
% the Bragg peak
iyTop = max(1, iyWater - 2);
iyBot = min(size(doseCube,1), iyPeak + 7);
xlim(gca, izPeak + [-marginVox, marginVox]);
ylim(gca, [iyTop iyBot]);
axis(gca, 'on');
xlabel(gca, 'z - z_{iso} [mm]', 'Color', 'k');
% The vertical axis here is the beam direction, so label it as DEPTH IN
% WATER measured from the phantom entrance surface - the same quantity the
% PDD figure above uses. That makes the two figures directly comparable:
% the Bragg peak sits at the same number on both.
ylabel(gca, 'Depth in water [mm]', 'Color', 'k');   % rows of plane-2 image = y (beam direction)
set(gca, 'XColor', 'k', 'YColor', 'k');
setNiceMmTicks(gca, 'x', dij.doseGrid.z, izPeak, marginVox, 5, stf(1).isoCenter(3));
setNiceMmTicks(gca, 'y', dij.doseGrid.y, round((iyTop+iyBot)/2), ceil((iyBot-iyTop)/2), 6, entryPoint(2));
styleSliceAxes(gca);
if ~isempty(targetPos) && ~isempty(hContourSag{targetPos})
    % outside the plot box, top right: matRad_plotSlice's axis settings make
    % an inside legend land unpredictably (it sat on the beam before)
    legend(hContourSag{targetPos}(1), {'OuterTarget'}, 'TextColor', 'k', ...
        'Location', 'northeastoutside', 'AutoUpdate', 'off');
end
saveFigure(gcf, fullfile(outDir, [dwaMachine '_depth_slice']));

%% Lateral profile at the Bragg-peak depth: core and halo
% On a linear colour scale the halo is invisible: it carries ~half the
% fluence but spreads it over sigma ~20 mm, so each voxel gets ~1% of the
% central value. A log scale shows both components.
latX_mm = dij.doseGrid.x(:) - stf(1).isoCenter(1);
latProf = squeeze(doseCube(iyPeak, :, izPeak));
latProf = latProf(:) / max(latProf);
figure('Color','w','Position',[100 100 900 360]);
subplot(1,2,1); plot(latX_mm, latProf, 'o-', 'MarkerSize', 3.5); grid on; box on;
xlim([-60 60]); xlabel('x - x_{iso} [mm]'); ylabel('relative dose'); title('linear');
subplot(1,2,2); semilogy(latX_mm, max(latProf,1e-6), 'o-', 'MarkerSize', 3.5); grid on; box on;
xlim([-60 60]); ylim([1e-4 1.2]); xlabel('x - x_{iso} [mm]'); ylabel('relative dose'); title('log - halo visible');
sgtitle(sprintf('%s lateral profile at %.0f mm depth, %.2f MeV (3 mm dose grid)', ...
    strrep(dwaMachine,'_','\_'), braggDepth_mm, stf(1).ray.energy));
saveFigure(gcf, fullfile(outDir, [dwaMachine '_lateral_profile_at_peak']));

%% Local helper
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
