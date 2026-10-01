%% DWA_proton vs Generic: energy sweep over all DWA layers, BOXPHANTOM
%
% Step 1 (Water Phantom Validation), full version. For EVERY energy layer
% in the DWA machine file:
%   1. fire a single pencil beam into BOXPHANTOM with DWA_proton,
%   2. read back the energy it actually used,
%   3. ask the Generic machine for THAT energy (matched-achieved-energy
%      approach - Generic snaps to its own nearest layer),
%   4. compare the two dose distributions.
%
% Matching on the ACHIEVED energy matters. If both machines are simply
% asked for the same nominal energy they can snap to different layers: at
% a request of 150 MeV, DWA picks 149.16 and Generic picks 150.35, a
% 1.19 MeV mismatch worth ~0.8 mm of range. Matching on what DWA actually
% delivered cuts the mismatch to a median of ~0.44 MeV across the sweep,
% so a residual range difference can be attributed to the range-energy
% relation rather than to the two energy grids not lining up.
%
% TWO CLASSES OF ENERGY CANNOT BE COMPARED IN THIS PHANTOM, and the script
% excludes them by default rather than producing meaningless numbers:
%   - DWA layers below Generic's lowest layer (31.73 MeV). DWA's 20.01 to
%     29.02 MeV layers would snap to 31.73, a mismatch of up to 11.7 MeV.
%   - DWA layers whose Bragg peak falls beyond the 240 mm water block
%     (roughly 195 MeV and above, peak > 246 mm). The peak is in air, so
%     range and falloff metrics are meaningless.
% Set energyMode = 'all' to run them anyway; the status column in the CSV
% records why each was excluded.
%
% Expect roughly 45 usable layers out of 58, spanning ~32 to ~190 MeV.
% Runtime is dominated by 2 dose calculations per layer.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd '/Users/mirtadumancic/Work/MATLAB_Projects/matRadDWA';
run('matRad_rc');
load('BOXPHANTOM.mat');

%% ------------------------------- configuration
outDir               = 'DWAvsGeneric_sweep';  % all output lands here
energyMode           = 'valid';               % 'valid' | 'all'
savePerEnergyFigures = true;                  % one 2x2 panel per energy
perEnergyFormats     = {'png'};               % {'png'} or {'png','pdf'} - 45 energies adds up
summaryFormats       = {'png','pdf','fig'};

isoDepth_mm     = 120;   % entrance surface -> isocenter (BOXPHANTOM geometry)
waterThickness  = 240;   % mm of water along the beam
peakDepthMargin = 15;    % require the peak this far inside the block, so the
                         % distal R20 is still measured in water

% Sample ON the dose-grid voxel centres rather than on a finer arbitrary
% grid. The dose is only known every 3 mm; anything finer is linear
% interpolation dressed up as data, and it hides where the calculation
% actually happened. Because griddedInterpolant is LINEAR here, sampling
% coarsely costs nothing in the metrics: R80/R20/FWHM are found by linear
% interpolation between the two bracketing samples, which for grid-aligned
% samples is exactly the same number the fine sampling produced. The plots
% below draw markers at these points so the sampling is visible.
%
% This lands exactly on voxel centres: the y grid runs -240:3:237, the
% entrance surface sits at y = -123 (a grid point), and the beam axis sits
% at z = -3 (also a grid point), so steps of one voxel stay on the grid.
depths_mm  = (0:ct.resolution.y:280)';
lateral_mm = (-60:ct.resolution.z:60)';

if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------------- which energies to run
plnD = basePln(ct, cst, 'DWA_proton');
plnG = basePln(ct, cst, 'Generic');
machD = matRad_loadMachine(plnD);
machG = matRad_loadMachine(plnG);

E_DWA  = [machD.data.energy];
E_Gen  = [machG.data.energy];

% peakPos is what tells us whether a layer's Bragg peak still lands inside
% the water block. It was added to the DWA machine file during
% commissioning; if this errors, the machine file on the path predates that.
if ~isfield(machD.data,'peakPos')
    error(['machine.data has no peakPos field. Re-run matRad_commissionDWAmachine.m ' ...
           '(which now sets data(i).peakPos) and copy the result into matRad/basedata.']);
end
pk_DWA = [machD.data.peakPos];

belowGeneric = E_DWA < min(E_Gen);
overshoot    = pk_DWA > (waterThickness - peakDepthMargin);

status = repmat({'ok'}, 1, numel(E_DWA));
status(belowGeneric) = {'below Generic energy range'};
status(overshoot)    = {'Bragg peak beyond water block'};
status(belowGeneric & overshoot) = {'below Generic range; peak beyond water'};

switch energyMode
    case 'valid', runMask = ~belowGeneric & ~overshoot;
    case 'all',   runMask = true(size(E_DWA));
    otherwise,    error('energyMode must be ''valid'' or ''all''');
end

fprintf('\nDWA layers: %d total, %d excluded, %d to run (%.2f - %.2f MeV)\n', ...
    numel(E_DWA), sum(~runMask), sum(runMask), min(E_DWA(runMask)), max(E_DWA(runMask)));
for k = find(~runMask)
    fprintf('   skipping %7.2f MeV : %s\n', E_DWA(k), status{k});
end

%% ------------------------------- the sweep
n = numel(E_DWA);
R = struct('E_request',num2cell(nan(1,n)),'E_DWA',num2cell(nan(1,n)), ...
           'E_Gen',num2cell(nan(1,n)),'dE',num2cell(nan(1,n)),'status',status, ...
           'peak_DWA',num2cell(nan(1,n)),'peak_Gen',num2cell(nan(1,n)), ...
           'r80_DWA',num2cell(nan(1,n)),'r80_Gen',num2cell(nan(1,n)), ...
           'fall_DWA',num2cell(nan(1,n)),'fall_Gen',num2cell(nan(1,n)), ...
           'ent_DWA',num2cell(nan(1,n)),'ent_Gen',num2cell(nan(1,n)), ...
           'fwhm_DWA',num2cell(nan(1,n)),'fwhm_Gen',num2cell(nan(1,n)), ...
           'pen_DWA',num2cell(nan(1,n)),'pen_Gen',num2cell(nan(1,n)));

tSweep = tic;
for k = find(runMask)

    fprintf('\n===== layer %d/%d : DWA request %.2f MeV =====\n', ...
        find(find(runMask)==k), sum(runMask), E_DWA(k));

    % --- DWA at this layer
    [dwa, eDWA] = runOneBeam(ct, cst, 'DWA_proton', E_DWA(k), ...
                             isoDepth_mm, depths_mm, lateral_mm);

    % --- Generic, asked for the energy DWA actually delivered
    [gen, eGen] = runOneBeam(ct, cst, 'Generic', eDWA, ...
                             isoDepth_mm, depths_mm, lateral_mm);

    fprintf('   DWA %.2f MeV  |  Generic %.2f MeV  |  mismatch %.2f MeV\n', ...
        eDWA, eGen, abs(eGen-eDWA));

    R(k).E_request = E_DWA(k);
    R(k).E_DWA = eDWA;  R(k).E_Gen = eGen;  R(k).dE = eGen - eDWA;
    R(k).peak_DWA = dwa.peakDepth;   R(k).peak_Gen = gen.peakDepth;
    R(k).r80_DWA  = dwa.r80;         R(k).r80_Gen  = gen.r80;
    R(k).fall_DWA = dwa.r20-dwa.r80; R(k).fall_Gen = gen.r20-gen.r80;
    R(k).ent_DWA  = dwa.entrancePct; R(k).ent_Gen  = gen.entrancePct;
    R(k).fwhm_DWA = dwa.lat.fwhm;    R(k).fwhm_Gen = gen.lat.fwhm;
    R(k).pen_DWA  = dwa.lat.pen8020; R(k).pen_Gen  = gen.lat.pen8020;

    if savePerEnergyFigures
        makePerEnergyFigure(dwa, gen, eDWA, eGen, depths_mm, lateral_mm, ...
                            isoDepth_mm, outDir, perEnergyFormats);
    end
end
fprintf('\nSweep finished in %.1f min\n', toc(tSweep)/60);

%% ------------------------------- results table
T = struct2table(R);
T = movevars(T, 'status', 'After', 'dE');
T.dPeak   = T.peak_Gen - T.peak_DWA;      % Generic minus DWA
T.dFall   = T.fall_Gen - T.fall_DWA;
T.ratioFWHM = T.fwhm_DWA ./ T.fwhm_Gen;
writetable(T, fullfile(outDir,'DWAvsGeneric_metrics.csv'));
fprintf('Wrote %s\n', fullfile(outDir,'DWAvsGeneric_metrics.csv'));

ok = ~isnan([R.E_DWA]);
fprintf('\n--- across %d compared layers ---\n', sum(ok));
fprintf('energy mismatch |dE|      : median %.2f, max %.2f MeV\n', median(abs(T.dE(ok))), max(abs(T.dE(ok))));
fprintf('peak depth, Generic - DWA : median %+.2f, range %+.2f .. %+.2f mm\n', ...
    median(T.dPeak(ok)), min(T.dPeak(ok)), max(T.dPeak(ok)));
fprintf('distal falloff R80-R20    : DWA median %.2f mm, Generic median %.2f mm\n', ...
    median(T.fall_DWA(ok)), median(T.fall_Gen(ok)));
fprintf('lateral FWHM at peak      : DWA median %.2f mm, Generic median %.2f mm (ratio %.2f)\n', ...
    median(T.fwhm_DWA(ok)), median(T.fwhm_Gen(ok)), median(T.ratioFWHM(ok)));

%% ------------------------------- summary figures across all energies
blue = [0.16 0.47 0.84]; orange = [0.92 0.41 0.20]; ink = [0.04 0.04 0.04];
E = T.E_DWA(ok);

% 1. range agreement
figure('Color','w'); hold on; grid on; box on;
plot(E, T.peak_DWA(ok), '-o','Color',blue,  'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','DWA\_proton');
plot(E, T.peak_Gen(ok), '-s','Color',orange,'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Bragg peak depth in water [mm]');
title('Range vs energy, matched achieved energy'); legend('Location','northwest');
tidyAxes(gca); saveFigure(gcf, fullfile(outDir,'summary_range_vs_energy'), summaryFormats);

% 2. range difference, with the residual energy mismatch underneath it
figure('Color','w');
subplot(2,1,1); hold on; grid on; box on;
plot(E, T.dPeak(ok), '-o','Color',ink,'LineWidth',1.6,'MarkerSize',4,'MarkerFaceColor','w');
yline(0,'-','Color',[.6 .6 .6],'HandleVisibility','off');
ylabel('\Deltapeak, Generic - DWA [mm]'); title('Range difference and the energy mismatch behind it');
tidyAxes(gca);
subplot(2,1,2); hold on; grid on; box on;
stem(E, T.dE(ok), 'Color',[.45 .45 .45],'MarkerSize',3,'MarkerFaceColor','w');
yline(0,'-','Color',[.6 .6 .6]);
xlabel('Energy [MeV]'); ylabel('\DeltaE, Generic - DWA [MeV]');
tidyAxes(gca); saveFigure(gcf, fullfile(outDir,'summary_range_difference'), summaryFormats);

% 3. distal falloff - the energy-spread signature
figure('Color','w'); hold on; grid on; box on;
plot(E, T.fall_DWA(ok), '-o','Color',blue,  'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','DWA\_proton');
plot(E, T.fall_Gen(ok), '-s','Color',orange,'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Distal falloff R80-R20 [mm]');
title('Distal falloff - sensitive to the modelled energy spread'); legend('Location','northwest');
tidyAxes(gca); saveFigure(gcf, fullfile(outDir,'summary_distal_falloff'), summaryFormats);

% 4. lateral width at the Bragg peak
figure('Color','w'); hold on; grid on; box on;
plot(E, T.fwhm_DWA(ok), '-o','Color',blue,  'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','DWA\_proton');
plot(E, T.fwhm_Gen(ok), '-s','Color',orange,'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Lateral FWHM at Bragg peak [mm]');
title('Lateral spot width at the Bragg peak'); legend('Location','northeast');
tidyAxes(gca); saveFigure(gcf, fullfile(outDir,'summary_lateral_fwhm'), summaryFormats);

% 5. entrance dose
figure('Color','w'); hold on; grid on; box on;
plot(E, T.ent_DWA(ok), '-o','Color',blue,  'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','DWA\_proton');
plot(E, T.ent_Gen(ok), '-s','Color',orange,'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w','DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Central-axis entrance dose [% of peak]');
title('Entrance dose on the central axis'); legend('Location','northwest');
tidyAxes(gca); saveFigure(gcf, fullfile(outDir,'summary_entrance_dose'), summaryFormats);

fprintf('\nAll output written to ./%s\n', outDir);

%% ================================================================ helpers

function pln = basePln(ct, cst, machineName)
    pln = struct();
    pln.radiationMode = 'protons';
    pln.machine       = machineName;
    pln.bioModel      = 'constRBE';
    pln.multScen      = 'nomScen';
    pln.propDoseCalc.calcLET = 0;
    pln.propDoseCalc.engine  = 'HongPB';
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
    pln.propStf.generator   = 'ParticleSingleSpot';
end

function [out, energyUsed] = runOneBeam(ct, cst, machineName, requestE, ...
                                        isoDepth_mm, depths_mm, lateral_mm)
% One single-spot dose calculation, with the curves and metrics extracted.
    pln = basePln(ct, cst, machineName);
    pln.propStf.energy = requestE;

    stf = matRad_generateStf(ct,cst,pln);
    energyUsed = stf(1).ray.energy;

    dij      = matRad_calcDoseInfluence(ct,cst,stf,pln);
    doseCube = reshape(full(dij.physicalDose{1}), dij.doseGrid.dimensions);

    beamDir    = stf(1).isoCenter - stf(1).sourcePoint;
    beamDir    = beamDir / norm(beamDir);
    entryPoint = stf(1).isoCenter - isoDepth_mm * beamDir;

    % matRad cubes are indexed [y, x, z] - grid vectors in that order, and
    % the query points swapped to (y, x, z) to match.
    F = griddedInterpolant({dij.doseGrid.y, dij.doseGrid.x, dij.doseGrid.z}, ...
                            doseCube, 'linear', 'none');
    % two-argument max ignores NaN, so off-grid samples become 0
    sampleAt = @(P) max(F(P(:,2), P(:,1), P(:,3)), 0);

    out.pdd = sampleAt(entryPoint + depths_mm * beamDir);
    [~, kPk] = max(out.pdd);
    out.peakDepth   = depths_mm(kPk);
    out.r80         = distalCrossing(depths_mm, out.pdd, 0.80);
    out.r20         = distalCrossing(depths_mm, out.pdd, 0.20);
    out.entrancePct = 100 * out.pdd(1) / max(out.pdd);
    out.latPeak     = sampleAt(entryPoint + out.peakDepth*beamDir + lateral_mm*[0 0 1]);
    out.lat         = lateralMetrics(lateral_mm, out.latPeak);
end

function makePerEnergyFigure(dwa, gen, eDWA, eGen, depths_mm, lateral_mm, ...
                             isoDepth_mm, outDir, formats)
% One 2x2 panel per energy: depth dose, distal edge, lateral, lateral on a
% log axis. The energy goes in the filename so the folder sorts sensibly.
    blue = [0.16 0.47 0.84]; orange = [0.92 0.41 0.20];
    % markers sit on the voxel centres where the dose was actually computed;
    % the straight segments between them are all the data supports
    pointStyle = {'LineWidth',1.4,'MarkerSize',3.5,'MarkerFaceColor','w'};
    nD = dwa.pdd/max(dwa.pdd); nG = gen.pdd/max(gen.pdd);
    lD = dwa.latPeak/max(dwa.latPeak); lG = gen.latPeak/max(gen.latPeak);
    labD = sprintf('DWA %.2f MeV', eDWA);
    labG = sprintf('Generic %.2f MeV', eGen);

    figure('Color','w','Position',[100 100 1100 800]);

    subplot(2,2,1); hold on; grid on; box on;
    plot(depths_mm, 100*nD, '-o','Color',blue,  pointStyle{:},'DisplayName',labD);
    plot(depths_mm, 100*nG, '-s','Color',orange,pointStyle{:},'DisplayName',labG);
    xline(isoDepth_mm,'--','Isocenter','HandleVisibility','off');
    xlabel('Depth in water [mm]'); ylabel('Dose [% of own max]');
    title('Central-axis depth dose'); legend('Location','northwest');
    xlim([0 min(280, max(dwa.peakDepth,gen.peakDepth)+40)]); tidyAxes(gca);

    subplot(2,2,2); hold on; grid on; box on;
    plot(depths_mm, 100*nD, '-o','Color',blue,  pointStyle{:},'DisplayName','DWA');
    plot(depths_mm, 100*nG, '-s','Color',orange,pointStyle{:},'DisplayName','Generic');
    yline([80 20],':',{'80%','20%'},'HandleVisibility','off');
    lo = min([dwa.r80 gen.r80]); hi = max([dwa.r20 gen.r20]);
    if isfinite(lo) && isfinite(hi), xlim([lo-12 hi+8]); end
    ylim([0 105]);
    xlabel('Depth in water [mm]'); ylabel('Dose [% of own max]');
    title(sprintf('Distal edge (falloff %.2f vs %.2f mm)', dwa.r20-dwa.r80, gen.r20-gen.r80));
    tidyAxes(gca);

    subplot(2,2,3); hold on; grid on; box on;
    plot(lateral_mm, 100*lD, '-o','Color',blue,  pointStyle{:},'DisplayName','DWA');
    plot(lateral_mm, 100*lG, '-s','Color',orange,pointStyle{:},'DisplayName','Generic');
    xlabel('Lateral offset z [mm]'); ylabel('Dose [% of own max]');
    title(sprintf('Lateral at Bragg peak (FWHM %.1f vs %.1f mm)', dwa.lat.fwhm, gen.lat.fwhm));
    xlim([-40 40]); tidyAxes(gca);

    subplot(2,2,4); hold on; grid on; box on;
    plot(lateral_mm, lD, '-o','Color',blue,  pointStyle{:},'DisplayName','DWA');
    plot(lateral_mm, lG, '-s','Color',orange,pointStyle{:},'DisplayName','Generic');
    set(gca,'YScale','log'); ylim([1e-4 1.5]); xlim([-50 50]);
    xlabel('Lateral offset z [mm]'); ylabel('Dose [fraction of own max]');
    title('Lateral, log scale (halo shows here)'); tidyAxes(gca);

    sgtitle(sprintf('DWA %.2f MeV  vs  Generic %.2f MeV   (\\DeltaE = %+.2f MeV)', ...
        eDWA, eGen, eGen-eDWA), 'FontWeight','bold');

    base = fullfile(outDir, sprintf('DWAvsGeneric_E%06.2fMeV', eDWA));
    saveFigure(gcf, base, formats);
    close(gcf);
end

function d = distalCrossing(depths, curve, frac)
% Depth at which the curve falls through frac*max on the DISTAL side.
    c = curve / max(curve);
    [~, k] = max(c);
    j = k - 1 + find(c(k:end) < frac, 1);
    if isempty(j) || j <= 1, d = NaN; return; end
    d = interp1(c([j j-1]), depths([j j-1]), frac);
end

function s = lateralMetrics(x, prof)
% Measured off the profile rather than fitted, so a core+halo shape is not
% forced into a single Gaussian.
    p = prof / max(prof);
    s.fwhm    = widthAt(x, p, 0.5);
    s.pen8020 = (widthAt(x, p, 0.2) - widthAt(x, p, 0.8)) / 2;   % per flank
    s.sigmaEq = s.fwhm / 2.355;
end

function w = widthAt(x, p, frac)
    [~, k] = max(p);
    iR = k - 1 + find(p(k:end) < frac, 1);
    iL = find(p(1:k) < frac, 1, 'last');
    if isempty(iR) || isempty(iL) || iR < 2 || iL >= numel(p), w = NaN; return; end
    xR = interp1(p([iR iR-1]), x([iR iR-1]), frac);
    xL = interp1(p([iL iL+1]), x([iL iL+1]), frac);
    w  = xR - xL;
end

function tidyAxes(ax)
    set(ax,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',10);
end

function saveFigure(figHandle, baseName, formats)
% 300 dpi PNG and/or vector PDF, rather than saveas() at screen resolution.
    set(figHandle, 'Color','w', 'InvertHardcopy','off');
    haveEG = ~isempty(which('exportgraphics'));
    for f = 1:numel(formats)
        switch formats{f}
            case 'png'
                if haveEG, exportgraphics(figHandle,[baseName '.png'],'Resolution',300,'BackgroundColor','white');
                else,      print(figHandle,[baseName '.png'],'-dpng','-r300'); end
            case 'pdf'
                if haveEG, exportgraphics(figHandle,[baseName '.pdf'],'ContentType','vector','BackgroundColor','white');
                else,      print(figHandle,[baseName '.pdf'],'-dpdf','-painters'); end
            case 'fig'
                savefig(figHandle,[baseName '.fig']);
        end
    end
end