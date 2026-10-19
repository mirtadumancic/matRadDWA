%% matRad_DWAvsGeneric_v2.m
%
% DWA vs Generic, single pencil beams in BOXPHANTOM, one per DWA energy
% layer - version 2, for the double-Gaussian DWA machine.
%
% Replaces matRad_DWAvsGeneric_energySweep.m AND
% matRad_DWAvsGeneric_fineFalloff.m (both kept unchanged for the record).
% Produces the same family of plots, in a NEW folder:
%     <folder of this script>/DWAvsGeneric_v2_<dwaMachine>/
%
% WHAT CHANGED AND WHY
% 1. Depth metrics come from the LATERALLY INTEGRATED depth dose (IDD), not
%    the central axis. The v1 scripts took the Bragg-peak depth as the
%    maximum of the central-axis curve. That worked for the old
%    single-Gaussian DWA (sigma ~5 mm everywhere) but not for the
%    double-Gaussian machine: its ~1.2 mm core spreads to ~5 mm by the peak,
%    so at high energy the central-axis maximum sits near the ENTRANCE, and
%    peak depth, R80, R20 and "entrance % of peak" would all be wrong. The
%    IDD is what the TOPAS Bragg curves are, so it compares like with like
%    and treats both machines identically.
% 2. One fine dose grid: 1 mm along the beam (the sampling of the TOPAS
%    curves) and 2 mm across it (resolves the 1.2 mm core; on 3 mm voxels
%    the lateral sum over-counts it by 6-13%). Dose is read from the sparse
%    dij directly, so no 27 M-voxel cube is ever built.
% 3. The water surface is read from the CT (isocentre is 121.5 mm deep, not
%    the 120 mm the v1 scripts assumed).
% 4. Double-Gaussian kernel forced (lateralModel = 'double'), dosimetric
%    lateral cutoff off and geometric cutoff 100 mm, for BOTH machines, so
%    neither halo is clipped.
% 5. New plot: share of the dose outside r = 10 mm at the Bragg peak - the
%    halo. Generic carries a nuclear halo, the DWA an in-air halo; FWHM alone
%    sees neither.
% 6. The "entrance dose" plot is now the Bragg curve's entrance-to-peak
%    ratio (IDD at 2 mm / IDD at the peak). The v1 plot used central-axis
%    dose, which for a narrow beam mostly measures how fast the spot widens.
%
% Energy matching is unchanged: DWA is fired first, then Generic is asked
% for the energy DWA actually used (Generic snaps to its nearest layer).

%% ------------------------------- configuration
scriptDir = fileparts(mfilename('fullpath'));
cd('/Users/mirtadumancic/Work/MATLAB_Projects/matRadDWA');
run('matRad_rc');
load('BOXPHANTOM.mat');

dwaMachine = 'DWA_proton_doubleGauss_interp';
dwaLayers  = 'native';    % 'native' = the 58 TOPAS layers (comparable with v1);
                          % 'all' = also the 143 interpolated layers
fineRes    = struct('x', 2, 'y', 1, 'z', 2);   % mm; beam travels along y

waterThick_mm   = 240;
peakMargin_mm   = 15;     % layer is used only if its peak is this far inside the water
haloRadius_mm   = 10;     % "outside" radius for the halo share
edgeEnergy_MeV  = 95.06;  % energy for the distal-edge close-up figure
resolvedFall_mm = 1.2;    % DWA falloffs below this are drawn hollow: they are at
                          % the 1 mm sampling of the TOPAS input curves themselves

savePerEnergyFigures = true;
perEnergyFormats     = {'png'};
summaryFormats       = {'png','pdf','fig'};

outDir = fullfile(scriptDir, ['DWAvsGeneric_v2_' dwaMachine]);
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------------- geometry: water surface from the CT
ct  = matRad_getWorldAxes(ct);
iso = matRad_getIsoCenter(cst, ct, 0);
[~, ixC] = min(abs(ct.x - iso(1)));
[~, izC] = min(abs(ct.z - iso(3)));
if isfield(ct,'cube') && ~isempty(ct.cube)
    colWater = ct.cube{1}(:, ixC, izC) > 0.5;
else
    colWater = ct.cubeHU{1}(:, ixC, izC) > -500;
end
surfaceY = ct.y(find(colWater, 1, 'first')) - ct.resolution.y/2;
fprintf('Water surface at y = %.1f mm; isocentre %.1f mm deep\n', surfaceY, iso(2) - surfaceY);

%% ------------------------------- which DWA layers
machD = matRad_loadMachine(basePln(ct, cst, dwaMachine, fineRes));
machG = matRad_loadMachine(basePln(ct, cst, 'Generic',  fineRes));
E_DWA  = arrayfun(@(d) d.energy,  machD.data);  E_DWA  = E_DWA(:)';
pk_DWA = arrayfun(@(d) d.peakPos, machD.data);  pk_DWA = pk_DWA(:)';
E_Gen  = arrayfun(@(d) d.energy,  machG.data);
isNat  = true(size(E_DWA));
if isfield(machD.data,'interpolated')
    isNat = ~arrayfun(@(d) isequal(d.interpolated, true), machD.data);
    isNat = isNat(:)';
end
use = E_DWA >= min(E_Gen) & pk_DWA <= waterThick_mm - peakMargin_mm;
if strcmp(dwaLayers,'native'), use = use & isNat; end
E_run = E_DWA(use);
fprintf('Running %d DWA layers, %.2f - %.2f MeV\n', numel(E_run), min(E_run), max(E_run));

%% ------------------------------- the sweep
n = numel(E_run);
R = repmat(struct('E_DWA',NaN,'E_Gen',NaN,'dE',NaN, ...
    'peak_DWA',NaN,'peak_Gen',NaN,'r80_DWA',NaN,'r80_Gen',NaN, ...
    'fall_DWA',NaN,'fall_Gen',NaN,'entr_DWA',NaN,'entr_Gen',NaN, ...
    'fwhm_DWA',NaN,'fwhm_Gen',NaN,'halo_DWA',NaN,'halo_Gen',NaN), n, 1);
C = cell(n, 2);   % curves per energy, for the per-energy and edge figures

tAll = tic;
for k = 1:n
    fprintf('\n=== %d/%d: DWA %.2f MeV ===\n', k, n, E_run(k));
    [d, eD] = runOneBeam(ct, cst, dwaMachine, E_run(k), fineRes, surfaceY, iso, haloRadius_mm);
    [g, eG] = runOneBeam(ct, cst, 'Generic',  eD,       fineRes, surfaceY, iso, haloRadius_mm);
    fprintf('   DWA %.2f | Generic %.2f MeV | peak %.1f vs %.1f mm | falloff %.2f vs %.2f mm | halo %.0f%% vs %.0f%%\n', ...
        eD, eG, d.peak, g.peak, d.fall, g.fall, 100*d.halo, 100*g.halo);
    R(k) = struct('E_DWA',eD,'E_Gen',eG,'dE',eG-eD, ...
        'peak_DWA',d.peak,'peak_Gen',g.peak,'r80_DWA',d.r80,'r80_Gen',g.r80, ...
        'fall_DWA',d.fall,'fall_Gen',g.fall,'entr_DWA',d.entr,'entr_Gen',g.entr, ...
        'fwhm_DWA',d.fwhm,'fwhm_Gen',g.fwhm,'halo_DWA',d.halo,'halo_Gen',g.halo);
    C{k,1} = d;  C{k,2} = g;
    if savePerEnergyFigures
        perEnergyFigure(d, g, eD, eG, dwaMachine, outDir, perEnergyFormats);
    end
    fprintf('   %.1f min elapsed, ~%.1f min left\n', toc(tAll)/60, toc(tAll)/60*(n-k)/k);
end

T = struct2table(R);
T.dPeak = T.peak_Gen - T.peak_DWA;
T.dR80  = T.r80_Gen  - T.r80_DWA;
writetable(T, fullfile(outDir, [dwaMachine '_vsGeneric_v2_metrics.csv']));

fprintf('\n--- %d layers ---\n', n);
fprintf('|dE| median %.2f MeV; R80 Generic-DWA median %+.2f mm (%+.2f .. %+.2f)\n', ...
    median(abs(T.dE)), median(T.dR80), min(T.dR80), max(T.dR80));
fprintf('falloff median: DWA %.2f, Generic %.2f mm\n', median(T.fall_DWA), median(T.fall_Gen));
fprintf('FWHM median: DWA %.2f, Generic %.2f mm; halo share median: DWA %.0f%%, Generic %.0f%%\n', ...
    median(T.fwhm_DWA), median(T.fwhm_Gen), 100*median(T.halo_DWA), 100*median(T.halo_Gen));

%% ------------------------------- summary figures
blue = [0.16 0.47 0.84]; orange = [0.92 0.41 0.20]; ink = [0.1 0.1 0.1];
E = T.E_DWA;
labD = strrep(dwaMachine,'_','\_');
base = @(name) fullfile(outDir, [dwaMachine '_v2_' name]);
sty  = {'LineWidth',1.6,'MarkerSize',5};

% 1. range vs energy
figure('Color','w'); hold on; grid on; box on;
plot(E, T.r80_DWA, '-o','Color',blue,  'MarkerFaceColor','w', sty{:}, 'DisplayName',labD);
plot(E, T.r80_Gen, '-s','Color',orange,'MarkerFaceColor','w', sty{:}, 'DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Range R80 in water [mm]');
title('Range vs energy (laterally integrated dose, 1 mm grid)'); legend('Location','northwest');
tidyAxes(gca); saveFigure(gcf, base('summary_range_vs_energy'), summaryFormats);

% 2. range difference + the energy mismatch behind it
figure('Color','w','Position',[100 100 800 600]);
subplot(2,1,1); hold on; grid on; box on;
plot(E, T.dR80,  '-o','Color',ink,'MarkerFaceColor','w', sty{:}, 'DisplayName','R80');
plot(E, T.dPeak, ':s','Color',[.5 .5 .5],'MarkerSize',4, 'DisplayName','Bragg peak');
yline(0,'-','Color',[.6 .6 .6],'HandleVisibility','off');
ylabel('Generic - DWA [mm]'); legend('Location','best');
title('Range difference and the energy mismatch behind it'); tidyAxes(gca);
subplot(2,1,2); hold on; grid on; box on;
stem(E, T.dE, 'Color',[.45 .45 .45],'MarkerSize',3,'MarkerFaceColor','w');
yline(0,'-','Color',[.6 .6 .6]);
xlabel('Energy [MeV]'); ylabel('\DeltaE, Generic - DWA [MeV]'); tidyAxes(gca);
saveFigure(gcf, base('summary_range_difference'), summaryFormats);

% 3. distal falloff
unres = T.fall_DWA < resolvedFall_mm;
figure('Color','w'); hold on; grid on; box on;
plot(E, T.fall_Gen, '-s','Color',orange,'MarkerFaceColor',orange, sty{:}, 'DisplayName','Generic');
plot(E, T.fall_DWA, '-','Color',blue, 'LineWidth',1.6, 'HandleVisibility','off');
plot(E(~unres), T.fall_DWA(~unres), 'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',6,'DisplayName',labD);
if any(unres)
    plot(E(unres), T.fall_DWA(unres), 'o','Color',blue,'MarkerFaceColor','w','MarkerSize',6, ...
        'DisplayName',sprintf('DWA, falloff < %.1f mm (at the input sampling)', resolvedFall_mm));
end
ylim([0 max([T.fall_Gen; T.fall_DWA])*1.1]);
xlabel('Energy [MeV]'); ylabel('Distal falloff R80 - R20 [mm]');
title('Distal falloff (laterally integrated dose, 1 mm grid)'); legend('Location','northeast');
tidyAxes(gca); saveFigure(gcf, base('summary_distal_falloff'), summaryFormats);

% 4. distal edge close-up at one energy
[~, kE] = min(abs(E - edgeEnergy_MeV));
d = C{kE,1}; g = C{kE,2};
figure('Color','w'); hold on; grid on; box on;
plot(d.z, 100*d.idd/d.pkVal, '-o','Color',blue,  'MarkerFaceColor','w','LineWidth',1.4,'MarkerSize',4, ...
    'DisplayName',sprintf('DWA %.2f MeV', d.E));
plot(g.z, 100*g.idd/g.pkVal, '-s','Color',orange,'MarkerFaceColor','w','LineWidth',1.4,'MarkerSize',4, ...
    'DisplayName',sprintf('Generic %.2f MeV', g.E));
yline([80 20], ':', {'80%','20%'}, 'HandleVisibility','off');
xlim([min(d.r80,g.r80)-20, max(d.r80,g.r80)+12]); ylim([0 105]);
xlabel('Depth in water [mm]'); ylabel('Laterally integrated dose [% of peak]');
title(sprintf('Distal edge at %.2f MeV, sampled every 1 mm', d.E)); legend('Location','southwest');
tidyAxes(gca); saveFigure(gcf, base('distal_edge_example'), summaryFormats);

% 5. Bragg curve entrance-to-peak
figure('Color','w'); hold on; grid on; box on;
plot(E, 100*T.entr_DWA, '-o','Color',blue,  'MarkerFaceColor','w', sty{:}, 'DisplayName',labD);
plot(E, 100*T.entr_Gen, '-s','Color',orange,'MarkerFaceColor','w', sty{:}, 'DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Entrance dose (2 mm) [% of Bragg peak]');
title('Bragg curve entrance-to-peak (laterally integrated)'); legend('Location','northwest');
tidyAxes(gca); saveFigure(gcf, base('summary_entrance_to_peak'), summaryFormats);

% 6. lateral FWHM at the Bragg peak
figure('Color','w'); hold on; grid on; box on;
plot(E, T.fwhm_DWA, '-o','Color',blue,  'MarkerFaceColor','w', sty{:}, 'DisplayName',labD);
plot(E, T.fwhm_Gen, '-s','Color',orange,'MarkerFaceColor','w', sty{:}, 'DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel('Lateral FWHM at Bragg peak [mm]');
title('Lateral spot width at the Bragg peak (2 mm grid)'); legend('Location','northeast');
tidyAxes(gca); saveFigure(gcf, base('summary_lateral_fwhm'), summaryFormats);

% 7. halo: share of the dose outside r = haloRadius at the Bragg peak
figure('Color','w'); hold on; grid on; box on;
plot(E, 100*T.halo_DWA, '-o','Color',blue,  'MarkerFaceColor','w', sty{:}, 'DisplayName',labD);
plot(E, 100*T.halo_Gen, '-s','Color',orange,'MarkerFaceColor','w', sty{:}, 'DisplayName','Generic');
xlabel('Energy [MeV]'); ylabel(sprintf('Dose outside r = %d mm [%% of plane]', haloRadius_mm));
title('Halo: dose outside the core at the Bragg-peak depth'); legend('Location','best');
tidyAxes(gca); saveFigure(gcf, base('summary_halo_share'), summaryFormats);

fprintf('\nAll output in %s  (%.1f min)\n', outDir, toc(tAll)/60);

%% ================================================================ helpers
function pln = basePln(ct, cst, machineName, res)
    pln = struct();
    pln.radiationMode = 'protons';
    pln.machine       = machineName;
    pln.bioModel      = 'constRBE';
    pln.multScen      = 'nomScen';
    pln.numOfFractions = 1;
    pln.propDoseCalc.calcLET  = 0;
    pln.propDoseCalc.engine   = 'HongPB';
    pln.propDoseCalc.lateralModel = 'double';          % both machines are double Gaussian
    pln.propDoseCalc.geometricLateralCutOff  = 100;    % mm, ~4 halo sigmas
    pln.propDoseCalc.dosimetricLateralCutOff = 1;      % keep the whole halo
    pln.propDoseCalc.doseGrid.resolution = res;
    pln.propStf.gantryAngles  = 0;
    pln.propStf.couchAngles   = 0;
    pln.propStf.bixelWidth    = 5;
    pln.propStf.isoCenter     = matRad_getIsoCenter(cst, ct, 0);
    pln.propStf.generator     = 'ParticleSingleSpot';
    pln.propOpt.runDAO        = 0;
    pln.propSeq.runSequencing = 0;
    pln.propOpt.quantityOpt   = 'physicalDose';
end

function [o, eUsed] = runOneBeam(ct, cst, machineName, requestE, res, surfaceY, iso, haloR)
% One single-spot calculation on the fine grid; all metrics read from the
% sparse dij (non-zero voxels only).
    pln = basePln(ct, cst, machineName, res);
    pln.propStf.energy = requestE;
    stf   = matRad_generateStf(ct, cst, pln);
    eUsed = stf(1).ray.energy;
    dij   = matRad_calcDoseInfluence(ct, cst, stf, pln);

    G = dij.doseGrid;   dims = G.dimensions;            % [Ny Nx Nz]
    [vIdx, ~, v] = find(dij.physicalDose{1}(:,1));
    [iy, ix, iz] = ind2sub(dims, vIdx);

    % laterally integrated depth dose, depth from the water surface
    idd = accumarray(iy, v, [dims(1) 1]) * res.x * res.z;
    z   = G.y(:) - surfaceY;
    in  = z >= 0 & z <= 280;
    z = z(in);  idd = idd(in);
    o.E = eUsed;  o.z = z;  o.idd = idd;

    % peak: parabola through the maximum and its neighbours (1 mm sampling)
    [~, m] = max(idd);
    if m > 1 && m < numel(idd)
        den = idd(m-1) - 2*idd(m) + idd(m+1);
        dz  = 0.5*(idd(m-1) - idd(m+1))/den;
        o.peak  = z(m) + dz*(z(2)-z(1));
        o.pkVal = idd(m) - 0.25*(idd(m-1) - idd(m+1))*dz;
    else
        o.peak = z(m);  o.pkVal = idd(m);
    end
    o.r80  = distalCrossing(z, idd/o.pkVal, m, 0.80);
    o.r20  = distalCrossing(z, idd/o.pkVal, m, 0.20);
    o.fall = o.r20 - o.r80;
    [~, i2] = min(abs(z - 2));
    o.entr = idd(i2) / o.pkVal;

    % lateral, in the Bragg-peak plane
    [~, yPk] = min(abs(G.y - (surfaceY + o.peak)));
    inPl = iy == yPk;
    xs = G.x(ix(inPl)) - iso(1);   zs = G.z(iz(inPl)) - iso(3);   vs = v(inPl);
    xs = xs(:); zs = zs(:); vs = vs(:);
    o.halo = sum(vs(xs.^2 + zs.^2 > haloR^2)) / sum(vs);
    % x profile: mean of the two z columns either side of the axis (the 2 mm
    % grid puts the axis between voxels)
    zc = unique(zs);  [~, srt] = sort(abs(zc));  zNear = zc(srt(1:min(2,end)));
    xGrid = (G.x(:) - iso(1));
    prof  = zeros(numel(xGrid),1);
    for zz = zNear(:)'
        sel  = abs(zs - zz) < 1e-6;
        xSel = xs(sel);  vSel = vs(sel);
        [tf, loc] = ismember(round(xSel*1000), round(xGrid*1000));
        prof = prof + accumarray(loc(tf), vSel(tf), [numel(xGrid) 1]);
    end
    o.latX = xGrid;  o.lat = prof / max(prof);
    o.fwhm = widthAt(o.latX, o.lat, 0.5);
end

function d = distalCrossing(z, c, m, frac)
    j = m - 1 + find(c(m:end) < frac, 1);
    if isempty(j) || j <= 1, d = NaN; return; end
    d = interp1(c([j j-1]), z([j j-1]), frac);
end

function w = widthAt(x, p, frac)
    [~, k] = max(p);
    iR = k - 1 + find(p(k:end) < frac, 1);
    iL = find(p(1:k) < frac, 1, 'last');
    if isempty(iR) || isempty(iL) || iR < 2 || iL >= numel(p), w = NaN; return; end
    w = interp1(p([iR iR-1]), x([iR iR-1]), frac) - interp1(p([iL iL+1]), x([iL iL+1]), frac);
end

function perEnergyFigure(d, g, eD, eG, dwaMachine, outDir, formats)
    blue = [0.16 0.47 0.84]; orange = [0.92 0.41 0.20];
    st = {'LineWidth',1.3,'MarkerSize',3,'MarkerFaceColor','w'};
    figure('Color','w','Position',[100 100 1100 800]);
    subplot(2,2,1); hold on; grid on; box on;
    plot(d.z, 100*d.idd/d.pkVal, '-o','Color',blue,  st{:}, 'DisplayName',sprintf('DWA %.2f MeV',eD));
    plot(g.z, 100*g.idd/g.pkVal, '-s','Color',orange,st{:}, 'DisplayName',sprintf('Generic %.2f MeV',eG));
    xlim([0 min(280, max(d.r20,g.r20)+30)]); ylim([0 105]);
    xlabel('Depth in water [mm]'); ylabel('Integrated dose [% of peak]');
    title('Laterally integrated depth dose'); legend('Location','northwest'); tidyAxes(gca);
    subplot(2,2,2); hold on; grid on; box on;
    plot(d.z, 100*d.idd/d.pkVal, '-o','Color',blue,  st{:});
    plot(g.z, 100*g.idd/g.pkVal, '-s','Color',orange,st{:});
    yline([80 20],':',{'80%','20%'});
    xlim([min(d.r80,g.r80)-15, max(d.r20,g.r20)+8]); ylim([0 105]);
    xlabel('Depth in water [mm]'); ylabel('Integrated dose [% of peak]');
    title(sprintf('Distal edge (R80-R20 %.2f vs %.2f mm)', d.fall, g.fall)); tidyAxes(gca);
    subplot(2,2,3); hold on; grid on; box on;
    plot(d.latX, 100*d.lat, '-o','Color',blue,  st{:});
    plot(g.latX, 100*g.lat, '-s','Color',orange,st{:});
    xlim([-30 30]); xlabel('x - x_{iso} [mm]'); ylabel('Dose [% of own max]');
    title(sprintf('Lateral at Bragg peak (FWHM %.1f vs %.1f mm)', d.fwhm, g.fwhm)); tidyAxes(gca);
    subplot(2,2,4); hold on; grid on; box on;
    semilogy(d.latX, max(d.lat,1e-6), '-o','Color',blue,  st{:});
    semilogy(g.latX, max(g.lat,1e-6), '-s','Color',orange,st{:});
    set(gca,'YScale','log'); ylim([1e-4 1.5]); xlim([-60 60]);
    xlabel('x - x_{iso} [mm]'); ylabel('Dose [fraction of own max]');
    title(sprintf('Lateral, log scale (outside 10 mm: %.0f%% vs %.0f%%)', 100*d.halo, 100*g.halo)); tidyAxes(gca);
    sgtitle(sprintf('%s %.2f MeV vs Generic %.2f MeV (\\DeltaE = %+.2f MeV)', ...
        strrep(dwaMachine,'_','\_'), eD, eG, eG-eD), 'FontWeight','bold');
    saveFigure(gcf, fullfile(outDir, sprintf('%s_v2_E%06.2fMeV', dwaMachine, eD)), formats);
    close(gcf);
end

function tidyAxes(ax)
    set(ax,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',10);
end

function saveFigure(figHandle, baseName, formats)
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
