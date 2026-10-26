%% Cross-validation: DWA_proton vs matRad's Generic proton machine, BOXPHANTOM
%
% Step 1 (Water Phantom Validation) of the DWA commissioning plan, second
% half: fire the SAME single pencil beam into the SAME water phantom with
% two different machine files and compare the beam models directly.
%
% Why a single pencil beam in a water box is the right place for this:
% it isolates the machine file. There is no optimisation, no second beam,
% no heterogeneity - so any difference between the two curves comes from
% the commissioned beam data and nothing else. (TG119 later deliberately
% adds the optimiser and realistic anatomy; if the beam models disagree
% there, you would not be able to tell which of the three caused it.)
%
% Comparator: 'Generic' is matRad's reference proton machine (DKFZ,
% double-Gaussian kernels fitted to MC). Note two structural differences
% that are worth keeping in mind when reading the output:
%   - machine.meta.dataType: Generic = 'doubleGauss', DWA = 'singleGauss'.
%     Generic models a narrow core PLUS a broad halo (at 148.72 MeV, at the
%     Bragg peak: sigma1 = 3.39 mm, sigma2 = 13.95 mm, halo weight 0.085).
%     The DWA air profiles show the same core+halo structure, but the DWA
%     machine file currently collapses it into one Gaussian.
%   - machine.meta.BAMStoIsoDist: Generic = 1000 mm, DWA = 400 mm.
%
% matRad has no Heidelberg-specific proton machine - and HIT is a
% synchrotron, not a cyclotron. The other protons_*.mat files in basedata
% are for the Monte Carlo engines (MCsquare / TOPAS), not for HongPB.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set matRad runtime configuration
cd '/Users/mirtadumancic/Work/MATLAB_Projects/matRadDWA';
run('matRad_rc');

load('BOXPHANTOM.mat');

%% What to compare
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
outDir = ['DWAvsGeneric_150MeV_' dwaMachine];
if ~exist(outDir,'dir'), mkdir(outDir); end

machines            = {dwaMachine, 'Generic'};
machineLabels       = {strrep(dwaMachine,'_','\_'), 'Generic'};
requestedEnergy_MeV = 150;    % each machine snaps to its own nearest layer
isoDepth_mm         = 120;    % entrance surface -> isocenter, BOXPHANTOM geometry

% Sampling for the extracted curves. The dose grid is 3 mm, so these are
% interpolated between voxels - fine for comparing two curves computed on
% the SAME grid, but do not read the penumbra numbers below as if they had
% sub-millimetre accuracy in an absolute sense.
depths_mm  = (0:0.5:280)';
lateral_mm = (-60:0.25:60)';

res = struct('machine',{},'energy',{},'pdd',{},'latPeak',{},'latIso',{}, ...
             'peakDepth',{},'r80',{},'r20',{},'entrancePct',{});

%% Run the same beam through each machine
for m = 1:numel(machines)

    pln = struct();
    pln.radiationMode = 'protons';
    pln.machine       = machines{m};
    pln.bioModel      = 'constRBE';
    pln.multScen      = 'nomScen';

    pln.propDoseCalc.calcLET = 0;
    pln.propDoseCalc.engine  = 'HongPB';
    % FORCE the double-Gaussian kernel. matRad's default lateralModel ('fast')
    % picks a SINGLE Gaussian whenever machine.data has a 'sigma' field - and
    % the DWA double-Gaussian machines carry one (the in-water MCS term, kept as
    % a fallback). Without this line matRad silently ignores the halo and
    % computes the core only, whatever meta.dataType says.
    pln.propDoseCalc.lateralModel = 'double';
    pln.propDoseCalc.geometricLateralCutOff = 100;   % mm; the halo sigma is ~20-23 mm
    pln.numOfFractions       = 1;

    pln.propStf.gantryAngles  = 0;
    pln.propStf.couchAngles   = 0;
    pln.propStf.bixelWidth    = 5;
    pln.propStf.isoCenter     = matRad_getIsoCenter(cst,ct,0);
    pln.propOpt.runDAO        = 0;
    pln.propSeq.runSequencing = 0;

    pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;
    pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y;
    pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;
    pln.propOpt.quantityOpt = 'physicalDose';

    % single ray, single energy, by construction
    pln.propStf.generator = 'ParticleSingleSpot';
    pln.propStf.energy    = requestedEnergy_MeV;

    stf = matRad_generateStf(ct,cst,pln);
    fprintf('\n[%s] using %.2f MeV (nearest layer to %.1f MeV requested)\n', ...
        machines{m}, stf(1).ray.energy, requestedEnergy_MeV);

    dij      = matRad_calcDoseInfluence(ct,cst,stf,pln);
    doseCube = reshape(full(dij.physicalDose{1}), dij.doseGrid.dimensions);

    % Geometry of this beam. matRad cubes are indexed [y, x, z]
    % (dij.doseGrid.dimensions = [numel(y) numel(x) numel(z)]), so the
    % interpolant takes the grid vectors in that order and the query
    % points must be swapped to (y, x, z) to match.
    beamDir    = (stf(1).isoCenter - stf(1).sourcePoint);
    beamDir    = beamDir / norm(beamDir);
    entryPoint = stf(1).isoCenter - isoDepth_mm * beamDir;

    F = griddedInterpolant({dij.doseGrid.y, dij.doseGrid.x, dij.doseGrid.z}, ...
                            doseCube, 'linear', 'none');
    % max(...,0) does double duty: MATLAB's two-argument max ignores NaN,
    % so points sampled outside the dose grid (NaN) become 0, and any tiny
    % negative interpolation artefact is clamped to 0 as well.
    sampleAt = @(P) max(F(P(:,2), P(:,1), P(:,3)), 0);

    % --- central-axis depth dose
    pdd = sampleAt(entryPoint + depths_mm * beamDir);

    % --- lateral profiles: along z, at the Bragg peak depth and at isocenter
    [~, kPeak]   = max(pdd);
    peakDepth_mm = depths_mm(kPeak);
    latAtDepth = @(D) sampleAt(entryPoint + D*beamDir + lateral_mm * [0 0 1]);
    latPeak = latAtDepth(peakDepth_mm);
    latIso  = latAtDepth(isoDepth_mm);

    % --- distal metrics off the depth-dose curve
    r80 = distalCrossing(depths_mm, pdd, 0.80);
    r20 = distalCrossing(depths_mm, pdd, 0.20);

    res(m).machine     = machines{m};
    res(m).energy      = stf(1).ray.energy;
    res(m).pdd         = pdd;
    res(m).latPeak     = latPeak;
    res(m).latIso      = latIso;
    res(m).peakDepth   = peakDepth_mm;
    res(m).r80         = r80;
    res(m).r20         = r20;
    res(m).entrancePct = 100 * pdd(1) / max(pdd);
end

%% Metrics table
fprintf('\n==================== DWA_proton vs Generic ====================\n');
fprintf('%-26s %14s %14s\n', '', machines{1}, machines{2});
rowf = @(name,a,b,fmt) fprintf(['%-26s ' fmt ' ' fmt '\n'], name, a, b);
rowf('energy used [MeV]',        res(1).energy,      res(2).energy,      '%14.2f');
rowf('Bragg peak depth [mm]',    res(1).peakDepth,   res(2).peakDepth,   '%14.1f');
rowf('distal R80 [mm]',          res(1).r80,         res(2).r80,         '%14.1f');
rowf('distal R20 [mm]',          res(1).r20,         res(2).r20,         '%14.1f');
rowf('distal falloff R80-R20',   res(1).r20-res(1).r80, res(2).r20-res(2).r80, '%14.2f');
rowf('entrance dose [% of peak]',res(1).entrancePct, res(2).entrancePct, '%14.1f');

for tag = {'latPeak','latIso'}
    where = ternary(strcmp(tag{1},'latPeak'),'at Bragg peak','at isocenter');
    f1 = lateralMetrics(lateral_mm, res(1).(tag{1}));
    f2 = lateralMetrics(lateral_mm, res(2).(tag{1}));
    fprintf('  -- lateral, %s --\n', where);
    rowf('   FWHM [mm]',             f1.fwhm,     f2.fwhm,     '%14.2f');
    rowf('   penumbra 80-20 [mm]',   f1.pen8020,  f2.pen8020,  '%14.2f');
    rowf('   sigma-equivalent [mm]', f1.sigmaEq,  f2.sigmaEq,  '%14.2f');
end
fprintf('==============================================================\n');
fprintf(['NOTE: dose is per unit pencil-beam weight and the Z kernels are not\n' ...
         'normalised to a common absolute scale, so only SHAPES are comparable\n' ...
         'here - every curve below is normalised to its own maximum.\n']);

%% Figure 1: depth-dose overlay
col = {[0.16 0.47 0.84], [0.92 0.41 0.20]};   % blue, orange
figure('Color','w'); hold on; grid on; box on;
for m = 1:numel(res)
    plot(depths_mm, 100*res(m).pdd/max(res(m).pdd), 'LineWidth', 2, 'Color', col{m}, ...
        'DisplayName', sprintf('%s @ %.2f MeV', machineLabels{m}, res(m).energy));
end
xline(isoDepth_mm, '--', 'Isocenter', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
xline([90 150], ':', {'Target proximal','Target distal'}, 'HandleVisibility','off');
xlabel('Depth in water [mm]'); ylabel('Relative dose [% of own maximum]');
title('Central-axis depth dose, single pencil beam in BOXPHANTOM');
legend('Location','northwest'); xlim([0 200]);
set(gca,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',11);
saveFigure(gcf, fullfile(outDir,[dwaMachine '_vsGeneric_PDD']));

%% Figure 2: distal edge, zoomed
figure('Color','w'); hold on; grid on; box on;
for m = 1:numel(res)
    plot(depths_mm, 100*res(m).pdd/max(res(m).pdd), 'LineWidth', 2, 'Color', col{m}, ...
        'DisplayName', sprintf('%s', machineLabels{m}));
end
yline([80 20], ':', {'80%','20%'}, 'HandleVisibility','off');
xlabel('Depth in water [mm]'); ylabel('Relative dose [% of own maximum]');
title('Distal falloff - sensitive to the modelled energy spread');
legend('Location','northeast');
if all(isfinite([res.r80])) && all(isfinite([res.r20]))
    xlim([min([res.r80])-15, max([res.r20])+10]);
end
ylim([0 105]);
set(gca,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',11);
saveFigure(gcf, fullfile(outDir,[dwaMachine '_vsGeneric_distalFalloff']));

%% Figure 3: lateral profiles at the Bragg peak depth
figure('Color','w'); hold on; grid on; box on;
for m = 1:numel(res)
    plot(lateral_mm, 100*res(m).latPeak/max(res(m).latPeak), 'LineWidth', 2, 'Color', col{m}, ...
        'DisplayName', sprintf('%s', machineLabels{m}));
end
xlabel('Lateral offset from beam axis, z [mm]'); ylabel('Relative dose [% of own maximum]');
title('Lateral profile at the Bragg peak depth');
legend('Location','northeast'); xlim([-40 40]);
set(gca,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',11);
saveFigure(gcf, fullfile(outDir,[dwaMachine '_vsGeneric_lateralPeak']));

%% Figure 4: same lateral profiles on a log axis - shows the halo, or its absence
figure('Color','w'); hold on; grid on; box on;
for m = 1:numel(res)
    semilogy(lateral_mm, res(m).latPeak/max(res(m).latPeak), 'LineWidth', 2, 'Color', col{m}, ...
        'DisplayName', sprintf('%s', machineLabels{m}));
end
set(gca,'YScale','log');
xlabel('Lateral offset from beam axis, z [mm]'); ylabel('Relative dose [fraction of own maximum]');
title('Lateral profile, log scale - Generic is doubleGauss, DWA is singleGauss');
legend('Location','northeast'); xlim([-50 50]); ylim([1e-4 1.5]);
set(gca,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',11);
saveFigure(gcf, fullfile(outDir,[dwaMachine '_vsGeneric_lateralPeak_log']));

%% ---------------------------------------------------------------- helpers
function d = distalCrossing(depths, curve, frac)
% Depth at which the curve falls through frac*max on the DISTAL side of the
% Bragg peak, linearly interpolated between the two bracketing samples.
    c = curve / max(curve);
    [~, k] = max(c);
    j = k - 1 + find(c(k:end) < frac, 1);
    if isempty(j) || j <= 1
        d = NaN; return;
    end
    d = interp1(c([j j-1]), depths([j j-1]), frac);
end

function s = lateralMetrics(x, prof)
% FWHM, 80-20 penumbra (mean of the two flanks) and the sigma a Gaussian of
% that FWHM would have. Measured off the profile directly rather than by
% fitting, so a core+halo profile is not forced into a single Gaussian.
    p = prof / max(prof);
    s.fwhm    = widthAt(x, p, 0.5);
    w80       = widthAt(x, p, 0.8);
    w20       = widthAt(x, p, 0.2);
    s.pen8020 = (w20 - w80) / 2;          % per flank
    s.sigmaEq = s.fwhm / 2.355;
end

function w = widthAt(x, p, frac)
% full width of the profile at a given fraction of its maximum
    [~, k] = max(p);
    iR = k - 1 + find(p(k:end) < frac, 1);
    iL = find(p(1:k) < frac, 1, 'last');
    if isempty(iR) || isempty(iL), w = NaN; return; end
    xR = interp1(p([iR iR-1]), x([iR iR-1]), frac);
    xL = interp1(p([iL iL+1]), x([iL iL+1]), frac);
    w  = xR - xL;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function saveFigure(figHandle, baseName)
% 300 dpi PNG + vector PDF + the .fig, rather than saveas() which
% rasterises at screen resolution.
    set(figHandle, 'Color', 'w', 'InvertHardcopy', 'off');
    if ~isempty(which('exportgraphics'))   % R2020a and later
        exportgraphics(figHandle, [baseName '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
        exportgraphics(figHandle, [baseName '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    else
        print(figHandle, [baseName '.png'], '-dpng', '-r300');
        print(figHandle, [baseName '.pdf'], '-dpdf', '-painters');
    end
    savefig(figHandle, [baseName '.fig']);
    fprintf('Saved %s .png / .pdf / .fig\n', baseName);
end
