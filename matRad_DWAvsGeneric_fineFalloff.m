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
outDir        = 'DWAvsGeneric_fineFalloff';
resAlongBeam  = 1.0;    % mm along y. Matches the 1 mm sampling of the TOPAS
                        % source curves - finer buys no real information.
energyList    = [];     % [] = a spread across the usable range; or e.g. [50 100 150]
isoDepth_mm   = 120;
waterThick    = 240;
peakMargin    = 15;

if ~exist(outDir,'dir'), mkdir(outDir); end

machD = matRad_loadMachine(basePln(ct,cst,'DWA_proton',3));
machG = matRad_loadMachine(basePln(ct,cst,'Generic',3));
E_DWA = [machD.data.energy];  pk_DWA = [machD.data.peakPos];
usable = E_DWA >= min([machG.data.energy]) & pk_DWA <= (waterThick - peakMargin);

if isempty(energyList)
    cand = E_DWA(usable);
    energyList = cand(round(linspace(1, numel(cand), 8)));   % 8 across the range
end
fprintf('\nFine-falloff run: %d energies, %.2f mm along the beam\n', numel(energyList), resAlongBeam);
fprintf('Grid along beam: %d voxels (vs %d at 3 mm)\n', ...
    round(480/resAlongBeam), round(480/ct.resolution.y));

%% ------------------------------- run both grids, both machines
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

    % coarse grid, as the sweep ran it
    [cD, eD] = runBeam(ct, cst, 'DWA_proton', energyList(k), ct.resolution.y, isoDepth_mm);
    [cG, eG] = runBeam(ct, cst, 'Generic',    eD,            ct.resolution.y, isoDepth_mm);

    % fine grid along the beam, same two machines, same matched energy
    fD = runBeam(ct, cst, 'DWA_proton', eD, resAlongBeam, isoDepth_mm);
    fG = runBeam(ct, cst, 'Generic',    eD, resAlongBeam, isoDepth_mm);

    Rows(k).E_DWA = eD;  Rows(k).E_Gen = eG;
    Rows(k).fall_DWA_coarse = cD.fall;  Rows(k).fall_Gen_coarse = cG.fall;
    Rows(k).fall_DWA_fine   = fD.fall;  Rows(k).fall_Gen_fine   = fG.fall;
    Rows(k).span_DWA_fine   = fD.fall / resAlongBeam;   % voxel intervals spanned
    Rows(k).span_Gen_fine   = fG.fall / resAlongBeam;
    Rows(k).gridLimited     = Rows(k).span_DWA_fine < 2;

    fprintf('   falloff DWA  : %5.2f mm (3 mm grid)  ->  %5.2f mm (%.1f mm grid, %.1f voxels)\n', ...
        cD.fall, fD.fall, resAlongBeam, Rows(k).span_DWA_fine);
    fprintf('   falloff Gen  : %5.2f mm (3 mm grid)  ->  %5.2f mm\n', cG.fall, fG.fall);
    if Rows(k).gridLimited
        fprintf(['   *** falloff spans only %.1f voxels - NOT resolved. Note the TOPAS\n' ...
                 '       curves are 1 mm sampled, so a finer dose grid will not fix this. ***\n'], ...
                 Rows(k).span_DWA_fine);
    end

    keep(end+1) = struct('depth',fD.depth,'dwa',fD.pdd,'gen',fG.pdd,'E',eD); %#ok<SAGROW>
end
fprintf('\nDone in %.1f min\n', toc(tAll)/60);

T = struct2table(Rows);
writetable(T, fullfile(outDir,'fineFalloff_metrics.csv'));
disp(T);

%% ------------------------------- figure 1: the artifact, before and after
blue=[0.16 0.47 0.84]; orange=[0.92 0.41 0.20];
figure('Color','w'); hold on; grid on; box on;
plot([Rows.E_DWA],[Rows.fall_DWA_coarse],'--o','Color',blue,'LineWidth',1.3,'MarkerSize',5, ...
     'MarkerFaceColor','w','DisplayName','DWA, 3 mm grid (artifact)');
plot([Rows.E_DWA],[Rows.fall_DWA_fine], '-o','Color',blue,'LineWidth',2,'MarkerSize',5, ...
     'MarkerFaceColor',blue,'DisplayName',sprintf('DWA, %.1f mm grid',resAlongBeam));
plot([Rows.E_DWA],[Rows.fall_Gen_coarse],'--s','Color',orange,'LineWidth',1.3,'MarkerSize',5, ...
     'MarkerFaceColor','w','DisplayName','Generic, 3 mm grid');
plot([Rows.E_DWA],[Rows.fall_Gen_fine], '-s','Color',orange,'LineWidth',2,'MarkerSize',5, ...
     'MarkerFaceColor',orange,'DisplayName',sprintf('Generic, %.1f mm grid',resAlongBeam));
yline(1.8,':','1.80 mm = 0.6 x 3 mm voxel','Color',[.45 .45 .45],'HandleVisibility','off');
xlabel('Energy [MeV]'); ylabel('Distal falloff R80-R20 [mm]');
title('Distal falloff: the 3 mm grid quantises it, the fine grid resolves it');
legend('Location','northwest'); tidyAxes(gca);
saveFig(gcf, fullfile(outDir,'falloff_coarse_vs_fine'));

% figure 2: distal edge at the middle energy, showing the actual samples
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
saveFig(gcf, fullfile(outDir,'distal_edge_fine'));

fprintf('\nOutput in ./%s\n', outDir);
if any([Rows.gridLimited])
    fprintf(['\nNOTE: %d of %d energies still have a falloff spanning under 2 voxels.\n' ...
             'Expect this below ~60 MeV: the DWA falloff there (<1.2 mm) is at or under the\n' ...
             '1 mm sampling of the TOPAS source curves themselves, so it cannot be validated\n' ...
             'from this dataset at ANY dose-grid resolution. Finer TOPAS output would be\n' ...
             'needed, not a finer grid here. Treat those layers as unmeasured.\n'], ...
        sum([Rows.gridLimited]), numel(Rows));
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
