%% Distal falloff on a fine dose grid: DWA_proton vs Generic, BOXPHANTOM
%
% WHY THIS SCRIPT EXISTS
% The energy sweep measured the distal falloff (R80-R20) on the default
% 3 mm dose grid and produced an obviously unphysical result: 18 of 43 DWA
% layers reported a falloff of EXACTLY 1.80 mm, with spikes in between.
%
% 1.80 mm is not a beam property, it is the voxel size. If the whole
% distal edge collapses into one 3 mm voxel interval, the dose falls
% linearly from 100% at the peak voxel to ~0 at the next, so R80 lands at
% 0.2*3 = 0.6 mm past the peak, R20 at 0.8*3 = 2.4 mm, and the difference
% is 0.6*3 = 1.80 mm for every such layer. The CSV confirmed it: those
% layers all had R80 - peakDepth = 0.6000 mm. The "spikes" are simply the
% layers where the peak sat differently relative to a voxel boundary, so
% the two crossings straddled two intervals instead of one.
%
% The DWA's true falloff is ~2.4 mm (from the raw TOPAS curves), i.e.
% SMALLER than one voxel - unresolvable on that grid. Generic's ~4 mm
% spans more than a voxel, which is why its curve came out smooth. The
% comparison was therefore not DWA vs Generic but unresolved vs resolved.
%
% WHAT THIS SCRIPT DOES DIFFERENTLY
%   - Refines the dose grid ALONG THE BEAM ONLY. The falloff is a depth
%     quantity; refining laterally would cost 27x the memory for nothing.
%     1 mm isotropic on this phantom is 480^3 = 110.6 M voxels (885 MB as
%     a full cube). 1 mm along y with x,z left at 3 mm is 12.3 M voxels.
%   - Reads the dose by INDEXING THE SPARSE dij directly at voxel centres,
%     instead of reshape(full(...)). No large cube is ever materialised,
%     and the samples are the computed voxels themselves rather than
%     interpolations between them.
%   - Runs each energy on BOTH grids, so the improvement is visible
%     directly rather than by comparison with an earlier run.
%   - Flags any falloff still spanning fewer than 2 voxel intervals as
%     unresolved, rather than reporting it as a measurement.
%
% WHY 1 mm AND NOT FINER
% The TOPAS depth-dose curves in DWA_PristineBraggPeaks.csv are themselves
% sampled every 1 mm, so data(i).Z is a 1 mm kernel. Refining the dose grid
% below 1 mm only interpolates that same curve more finely - it adds no
% information. 1 mm along the beam matches the input data and is the
% sensible floor.
%
% A SEPARATE, HARDER LIMIT AT LOW ENERGY
% The DWA's true falloff runs from ~0.6 mm at 20 MeV to ~2.7 mm at 163 MeV.
% Below roughly 60 MeV it is at or under the 1 mm sampling of the source
% data itself, so those layers cannot be validated from this dataset at any
% dose-grid resolution. That is a data limitation, not a grid one, and the
% fix would be finer-sampled TOPAS curves - not a finer dose grid here.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd '/Users/mirtadumancic/Work/MATLAB_Projects/matRadDWA';
run('matRad_rc');
load('BOXPHANTOM.mat');

%% ------------------------------- configuration
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

outDir        = ['DWAvsGeneric_fineFalloff_' dwaMachine];
dwaLayers     = 'native';   % 'native' = 58 TOPAS layers only (1:1 with older machines); 'all'
resAlongBeam  = 1.0;    % mm along y. Matches the 1 mm sampling of the TOPAS
                        % source curves - finer buys no real information.
energyList    = 'all';  % 'all' = every usable DWA layer; or a vector e.g. [50 100 150]
runCoarse     = false;  % also run the 3 mm grid for comparison. The artifact is
                        % already understood and the coarse numbers are in the
                        % sweep CSV, so this defaults off - it halves the runtime.
makeEdgeFigure = true;  % second figure: the distal edge at one energy
isoDepth_mm   = 120;
waterThick    = 240;
peakMargin    = 15;

if ~exist(outDir,'dir'), mkdir(outDir); end

machD = matRad_loadMachine(basePln(ct,cst,dwaMachine,3));
machG = matRad_loadMachine(basePln(ct,cst,'Generic',3));
E_DWA = [machD.data.energy];  pk_DWA = [machD.data.peakPos];
usable = E_DWA >= min([machG.data.energy]) & pk_DWA <= (waterThick - peakMargin);
if strcmp(dwaLayers,'native') && isfield(machD.data,'interpolated')
    usable = usable & ~[machD.data.interpolated];
end

if ischar(energyList) || isstring(energyList)
    energyList = E_DWA(usable);          % every layer we can legitimately compare
end
nRuns = numel(energyList) * (2 + 2*runCoarse);
fprintf('\nFine-falloff run: %d energies, %.2f mm along the beam\n', numel(energyList), resAlongBeam);
fprintf('Grid along beam: %d voxels (vs %d at 3 mm)\n', ...
    round(480/resAlongBeam), round(480/ct.resolution.y));
fprintf('%d dose calculations queued%s\n', nRuns, ...
    repmat(' (coarse comparison off)', 1, ~runCoarse));

%% ------------------------------- the sweep
nE = numel(energyList);
Rows = struct('E_DWA',[],'E_Gen',[], ...
              'fall_DWA_coarse',[],'fall_Gen_coarse',[], ...
              'fall_DWA_fine',[],  'fall_Gen_fine',[], ...
              'span_DWA_fine',[],  'span_Gen_fine',[], ...
              'gridLimited',[]);
Rows = repmat(Rows, 1, nE);
keep = struct('depth',{},'dwa',{},'gen',{},'E',{});

tAll = tic;
for k = 1:nE

    fprintf('\n=== %d/%d : requesting %.2f MeV ===\n', k, nE, energyList(k));

    tK = tic;

    % fine grid along the beam. DWA first, then Generic asked for the energy
    % DWA actually delivered (matched achieved energy).
    [fD, eD] = runBeam(ct, cst, dwaMachine, energyList(k), resAlongBeam, isoDepth_mm);
    [fG, eG] = runBeam(ct, cst, 'Generic',    eD,            resAlongBeam, isoDepth_mm);

    if runCoarse
        cD = runBeam(ct, cst, dwaMachine, eD, ct.resolution.y, isoDepth_mm);
        cG = runBeam(ct, cst, 'Generic',    eD, ct.resolution.y, isoDepth_mm);
    else
        cD.fall = NaN; cG.fall = NaN;
    end

    Rows(k).E_DWA = eD;  Rows(k).E_Gen = eG;
    Rows(k).fall_DWA_coarse = cD.fall;  Rows(k).fall_Gen_coarse = cG.fall;
    Rows(k).fall_DWA_fine   = fD.fall;  Rows(k).fall_Gen_fine   = fG.fall;
    Rows(k).span_DWA_fine   = fD.fall / resAlongBeam;   % voxel intervals spanned
    Rows(k).span_Gen_fine   = fG.fall / resAlongBeam;
    Rows(k).gridLimited     = Rows(k).span_DWA_fine < 2;

    fprintf('   DWA %.2f MeV vs Generic %.2f MeV | falloff %5.2f vs %5.2f mm (%.1f voxels)\n', ...
        eD, eG, fD.fall, fG.fall, Rows(k).span_DWA_fine);
    if k == 1
        fprintf('   [%.0f s for this energy -> roughly %.0f min for all %d]\n', ...
            toc(tK), toc(tK)*numel(energyList)/60, numel(energyList));
    end
    if Rows(k).gridLimited
        fprintf(['   *** falloff spans only %.1f voxels - NOT resolved. Note the TOPAS\n' ...
                 '       curves are 1 mm sampled, so a finer dose grid will not fix this. ***\n'], ...
                 Rows(k).span_DWA_fine);
    end

    keep(end+1) = struct('depth',fD.depth,'dwa',fD.pdd,'gen',fG.pdd,'E',eD); %#ok<SAGROW>
end
fprintf('\nDone in %.1f min\n', toc(tAll)/60);

T = struct2table(Rows);
writetable(T, fullfile(outDir,[dwaMachine '_fineFalloff_metrics.csv']));
fprintf('Wrote %s\n', fullfile(outDir,[dwaMachine '_fineFalloff_metrics.csv']));

res = ~logical([Rows.gridLimited]);
fprintf('\n--- resolved layers only (%d of %d) ---\n', sum(res), numel(Rows));
fprintf('DWA     falloff: %.2f .. %.2f mm (median %.2f)\n', ...
    min([Rows(res).fall_DWA_fine]), max([Rows(res).fall_DWA_fine]), median([Rows(res).fall_DWA_fine]));
fprintf('Generic falloff: %.2f .. %.2f mm (median %.2f)\n', ...
    min([Rows(res).fall_Gen_fine]), max([Rows(res).fall_Gen_fine]), median([Rows(res).fall_Gen_fine]));
fprintf('DWA is sharper by a median factor of %.2f\n', ...
    median([Rows(res).fall_Gen_fine] ./ [Rows(res).fall_DWA_fine]));

%% ------------------------------- the figure
% One panel, two series. DWA points whose falloff spans fewer than 2 voxel
% intervals are drawn HOLLOW: there the number is read off roughly one
% sampling interval, and the DWA falloff is comparable to the 1 mm spacing
% of the TOPAS curves the machine was built from, so it is shown but not
% presented as a measurement.
%
% Deliberately NOT shaded as a vertical band. A band spans both series, and
% Generic is fully resolved at those energies (its falloff there is >5 mm,
% i.e. five samples) - shading the strip would wrongly imply both curves
% are suspect. The limitation belongs to the DWA series alone, so it is
% marked on those markers alone.
blue=[0.16 0.47 0.84]; orange=[0.92 0.41 0.20]; ink=[0.10 0.10 0.10];
E      = [Rows.E_DWA];
fDWA   = [Rows.fall_DWA_fine];
fGen   = [Rows.fall_Gen_fine];
unres  = logical([Rows.gridLimited]);

figure('Color','w','Position',[100 100 900 620]); hold on; grid on; box on;
set(gca,'Layer','top');

plot(E, fGen, '-', 'Color',orange,'LineWidth',2,'HandleVisibility','off');
plot(E, fDWA, '-', 'Color',blue,  'LineWidth',2,'HandleVisibility','off');

% Generic: filled throughout - resolved at every energy shown
plot(E, fGen, 's','Color',orange,'MarkerFaceColor',orange, ...
     'MarkerSize',6,'LineWidth',1.2,'DisplayName','Generic');

% DWA: filled where resolved, hollow where not
plot(E(~unres), fDWA(~unres), 'o','Color',blue,'MarkerFaceColor',blue, ...
     'MarkerSize',6,'LineWidth',1.2,'DisplayName','DWA\_proton');
if any(unres)
    plot(E(unres), fDWA(unres), 'o','Color',blue,'MarkerFaceColor','w', ...
         'MarkerSize',6,'LineWidth',1.2, ...
         'DisplayName','DWA, falloff < 2 mm (not resolved)');
end

xlabel('Proton energy [MeV]');
ylabel('Distal falloff, R80 - R20 [mm]');
title(sprintf('Distal dose falloff in water, %.1f mm dose grid along the beam', resAlongBeam));
if ~isempty(which('subtitle'))
    subtitle(['Hollow markers: the DWA falloff there is comparable to the 1 mm sampling ' ...
              'of the TOPAS input curves'], 'FontSize',10,'Color',[0.40 0.40 0.40]);
end
legend('Location','northeast');
xlim([min(E)-5 max(E)+5]); ylim([0 max([fDWA fGen])*1.12]);
set(gca,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',12, ...
        'XColor',ink,'YColor',ink);
saveFig(gcf, fullfile(outDir,[dwaMachine '_distal_falloff_vs_energy']));

%% optional second figure: the distal edge itself at one energy
if makeEdgeFigure && ~isempty(keep)
    m = ceil(numel(keep)/2);
    figure('Color','w'); hold on; grid on; box on;
    nD = keep(m).dwa/max(keep(m).dwa); nG = keep(m).gen/max(keep(m).gen);
    plot(keep(m).depth, 100*nD, '-o','Color',blue,  'LineWidth',1.6,'MarkerSize',3.5,'MarkerFaceColor','w','DisplayName','DWA');
    plot(keep(m).depth, 100*nG, '-s','Color',orange,'LineWidth',1.6,'MarkerSize',3.5,'MarkerFaceColor','w','DisplayName','Generic');
    yline([80 20],':',{'80%','20%'},'HandleVisibility','off');
    [~,kk]=max(nD); xlim(keep(m).depth(kk)+[-20 15]); ylim([0 105]);
    xlabel('Depth in water [mm]'); ylabel('Dose [% of own max]');
    title(sprintf('Distal edge at %.2f MeV, sampled every %.1f mm', keep(m).E, resAlongBeam));
    legend('Location','southwest'); tidyAxes(gca);
    saveFig(gcf, fullfile(outDir,[dwaMachine '_distal_edge_fine']));
end

fprintf('\nOutput in ./%s\n', outDir);
if any(unres)
    fprintf(['\nNOTE: %d of %d energies have a falloff spanning under 2 voxels and are\n' ...
             'drawn hollow. Below ~%.0f MeV the DWA falloff is at or under the 1 mm\n' ...
             'sampling of the TOPAS source curves, so it cannot be validated from this\n' ...
             'dataset at ANY dose-grid resolution. Finer TOPAS output would be needed.\n'], ...
        sum(unres), numel(Rows), max(E(unres)));
end

%% ================================================================ helpers

function pln = basePln(ct, cst, machineName, resY)
    pln = struct();
    pln.radiationMode = 'protons';
    pln.machine       = machineName;
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
    % Anisotropic on purpose: fine along the beam (y), unchanged laterally.
    pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x;
    pln.propDoseCalc.doseGrid.resolution.y = resY;
    pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z;
    pln.propOpt.quantityOpt = 'physicalDose';
    pln.propStf.generator   = 'ParticleSingleSpot';
end

function [out, energyUsed] = runBeam(ct, cst, machineName, requestE, resY, isoDepth_mm)
    pln = basePln(ct, cst, machineName, resY);
    pln.propStf.energy = requestE;

    stf = matRad_generateStf(ct,cst,pln);
    energyUsed = stf(1).ray.energy;
    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

    beamDir    = stf(1).isoCenter - stf(1).sourcePoint;
    beamDir    = beamDir / norm(beamDir);
    entryPoint = stf(1).isoCenter - isoDepth_mm * beamDir;

    % Sample at the dose-grid voxel centres along the beam. The y grid runs
    % from ctGrid.y(1) in steps of resY and the entrance surface sits on a
    % grid point, so stepping by resY stays exactly on voxel centres.
    depth = (0 : resY : 280)';
    P     = entryPoint + depth * beamDir;

    out.depth = depth;
    out.pdd   = sampleVoxels(dij, dij.physicalDose{1}, P);
    out.r80   = distalCrossing(depth, out.pdd, 0.80);
    out.r20   = distalCrossing(depth, out.pdd, 0.20);
    out.fall  = out.r20 - out.r80;
end

function v = sampleVoxels(dij, doseSparse, P)
% Read the sparse dij directly at the voxel containing each point, instead
% of reshape(full(...)) - which on a 12 M voxel grid would allocate ~100 MB
% (and ~885 MB if the grid were refined isotropically).
    gx = dij.doseGrid.x; gy = dij.doseGrid.y; gz = dij.doseGrid.z;
    ix = round((P(:,1)-gx(1)) / (gx(2)-gx(1))) + 1;
    iy = round((P(:,2)-gy(1)) / (gy(2)-gy(1))) + 1;
    iz = round((P(:,3)-gz(1)) / (gz(2)-gz(1))) + 1;
    inside = ix>=1 & ix<=numel(gx) & iy>=1 & iy<=numel(gy) & iz>=1 & iz<=numel(gz);
    v = zeros(size(P,1),1);
    lin = sub2ind(dij.doseGrid.dimensions, iy(inside), ix(inside), iz(inside));
    v(inside) = full(doseSparse(lin));
end

function d = distalCrossing(depths, curve, frac)
    if max(curve) <= 0, d = NaN; return; end
    c = curve / max(curve);
    [~, k] = max(c);
    j = k - 1 + find(c(k:end) < frac, 1);
    if isempty(j) || j <= 1, d = NaN; return; end
    d = interp1(c([j j-1]), depths([j j-1]), frac);
end

function tidyAxes(ax)
    set(ax,'TickDir','out','TickLength',[0.012 0.012],'LineWidth',0.9,'FontSize',11);
end

function saveFig(figHandle, baseName)
    set(figHandle,'Color','w','InvertHardcopy','off');
    if ~isempty(which('exportgraphics'))
        exportgraphics(figHandle,[baseName '.png'],'Resolution',300,'BackgroundColor','white');
        exportgraphics(figHandle,[baseName '.pdf'],'ContentType','vector','BackgroundColor','white');
    else
        print(figHandle,[baseName '.png'],'-dpng','-r300');
    end
    savefig(figHandle,[baseName '.fig']);
end
