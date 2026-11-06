%% matRad_BraggPeaks_vsTOPAS_RayStation.m
% Integrated depth-dose curves (IDD) of single DWA spots in water, computed by
% matRad, for comparison with the TOPAS "ground truth" (Alaina's commissioning
% data, DWA_PristineBraggPeaks.csv) and with RayStation (pencil beam and MC).
%
% Geometry = commissioning geometry of the TOPAS / RayStation files:
%   single spot, gantry 0, water surface AT the isocentre (iso-to-surface 0 mm),
%   isocentre-to-snout 400 mm, dose integrated over a circular detector of
%   120 mm diameter (r <= 60 mm), 1 mm depth bins.
%
% Engines (one run per engine; set `engine` below):
%   'HongPB'   - matRad pencil beam (double Gaussian). Run on the Mac.
%   'MCsquare' - matRad Monte Carlo. Run on the Windows PC only (the Mac
%                MCsquare binary is x86 and fails on Apple Silicon), from a path
%                without spaces (X:\matRadDWA_LET, via subst X: "M:\My Documents").
%
% Output: DWA_final_machine/BraggPeaks_vsTOPAS_RS/IDD_<machine>_<engine>.csv
%   columns: energy_MeV, depth_mm (from the water surface), idd_r60, idd_all
%   (idd_all = whole plane, for reference). Curves are not normalised here.
%
% Run from the matRadDWA folder:
%   matRad_rc; run('DWA_final_machine/matRad_BraggPeaks_vsTOPAS_RayStation.m')

clear; clc;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_cfg = matRad_rc;
matRad_cfg.disableGUI = true;

%% ---------------- settings ------------------------------------------------
machineName = 'DWA_proton_v2';
engine      = 'HongPB';          % 'HongPB' (Mac) or 'MCsquare' (PC)
energies    = 'all';             % 'all' (58 TOPAS energies) or a vector in MeV, e.g. [70 100 150 200 228.84]
nHistMC     = 1e6;               % MCsquare primaries per energy
detRadius_mm = 60;               % commissioning detector: 120 mm diameter
waterDepth_mm = 430;             % > range of 229 MeV (~326 mm)
latHalf_mm  = 80;                % half-width of the water box laterally (> detector radius)
resLat_mm   = 2;  
resDepth_mm = 1;
airGap_mm   = 20;                % air in front of the water surface inside the CT

outDir = fullfile(rootDir, 'DWA_final_machine', 'BraggPeaks_vsTOPAS_RS');
if ~exist(outDir, 'dir'), mkdir(outDir); end
outFile = fullfile(outDir, sprintf('IDD_%s_%s_noAirOffset.csv', machineName, engine));

%% ---------------- MCsquare needs text metadata in the machine file -------
if strcmp(engine, 'MCsquare')
    mf = fullfile(rootDir, 'matRad', 'basedata', ['protons_' machineName '.mat']);
    S = load(mf, 'machine');
    need = {'description', machineName; 'created_by', 'Mirta Dumancic'; 'created_on', datestr(now, 'yyyy-mm-dd')};
    changed = false;
    for k = 1:size(need, 1)
        if ~isfield(S.machine.meta, need{k,1}) || isempty(S.machine.meta.(need{k,1}))
            S.machine.meta.(need{k,1}) = need{k,2}; changed = true;
        end
    end
    if changed, machine = S.machine; save(mf, 'machine'); clear machine; end
    clear S
    if any(rootDir == ' '), warning('MCsquare fails on paths with spaces: run from X:\\ (subst).'); end
end

%% ---------------- water box: beam along +y, surface at the isocentre ------
Nx = round(2*latHalf_mm/resLat_mm); Nz = Nx;
Ny = round((airGap_mm + waterDepth_mm)/resDepth_mm);
ct = struct();
ct.resolution = struct('x', resLat_mm, 'y', resDepth_mm, 'z', resLat_mm);
ct.cubeDim = [Ny Nx Nz];
ct.numOfCtScen = 1;
ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
iyW = round(airGap_mm/resDepth_mm) + 1;                 % first water voxel (depth index)
water = false(ct.cubeDim); water(iyW:end, :, :) = true;
% Air in front of the surface must have ~zero water-equivalent thickness: the dose engines
% recompute density from HU with ct.hlut, and HU -1000 with this HLUT gives rED 0.023
% (20 mm air = 0.47 mm WET, which shifted every Bragg peak 0.45 mm shallower).
rED = zeros(ct.cubeDim); rED(water) = 1;
HU  = -1024*ones(ct.cubeDim); HU(water) = 0;
ct.cube = {rED}; ct.cubeHU = {HU};
ct.hlut = [-1024 0.00; 0 1.00; 3000 2.5];
surfaceY = ct.y(iyW) - resDepth_mm/2;
iso = [mean(ct.x), surfaceY, mean(ct.z)];               % water surface at the isocentre

[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
tgt = water & abs(Yg - (surfaceY + 100)) <= 2 & abs(Xg - iso(1)) <= 2 & abs(Zg - iso(3)) <= 2;
cst = cell(2, 6);
cst(1,:) = {0, 'TARGET', 'TARGET', {find(tgt)},   struct('Priority', 1, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1), {}};
cst(2,:) = {1, 'WATER',  'OAR',    {find(water)}, struct('Priority', 2, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1), {}};
r2 = (Xg - iso(1)).^2 + (Zg - iso(3)).^2;
inDet = squeeze(r2(1, :, :)) <= detRadius_mm^2;          % [Nx Nz] detector disc
voxArea = resLat_mm^2;                                    % mm^2 per voxel column
depth = ct.y(:) - surfaceY;                               % mm from the surface (voxel centres)

%% ---------------- plan --------------------------------------------------
pln = struct();
pln.radiationMode = 'protons';  pln.machine = machineName;
pln.bioModel = 'constRBE';      pln.multScen = 'nomScen';  pln.numOfFractions = 1;
pln.propStf.gantryAngles = 0;   pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 5;     pln.propStf.isoCenter = iso;
pln.propStf.generator = 'ParticleSingleSpot';
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
pln.propDoseCalc.calcLET = 0;
if strcmp(engine, 'HongPB')
    pln.propDoseCalc.engine = 'HongPB';
    pln.propDoseCalc.lateralModel = 'double';
    pln.propDoseCalc.geometricLateralCutOff = 100;   % mm, keeps the halo
    pln.propDoseCalc.dosimetricLateralCutOff = 1;    % no dosimetric cut-off
    % matRad adds the WEPL of the nozzle-to-skin air minus machine.meta.fitAirOffset.
    % DWA_proton_v2 has fitAirOffset = 0 although the TOPAS IDDs already include the
    % 400 mm snout-to-surface air, so every peak came out 0.0011*400 = 0.44 mm shallower.
    % Here (surface at the isocentre = commissioning geometry) switching the correction off
    % is exact. The permanent fix is fitAirOffset = 400 in the next machine (DWA_proton_v3).
    pln.propDoseCalc.airOffsetCorrection = false;
else
    pln.propDoseCalc.engine = 'MCsquare';
    pln.propDoseCalc.numHistoriesDirect = nHistMC;
end
pln.propOpt.runDAO = 0; pln.propSeq.runSequencing = 0;

machine = matRad_loadMachine(pln);
allE = [machine.data.energy];
if ischar(energies)
    % the 58 energies of the TOPAS commissioning data (nearest machine layer)
    Etopas = [20.01 23.01 26.02 29.02 32.02 35.02 38.02 41.03 44.03 47.03 50.03 53.03 ...
              56.03 59.04 62.04 65.04 68.04 71.04 74.05 77.05 80.05 83.05 86.05 89.06 ...
              92.06 95.06 98.06 99.65 102.42 106.64 110.11 115.16 118.76 122.54 127.46 131.59 ...
              136.15 140.72 144.90 149.16 153.14 158.19 163.12 167.89 172.00 176.32 181.01 185.86 ...
              190.30 195.22 199.75 204.67 209.29 214.25 219.20 224.05 226.49 228.84];
    energies = Etopas;
end

%% ---------------- loop over energies ------------------------------------
T = table();
for k = 1:numel(energies)
    [dE, ix] = min(abs(allE - energies(k)));
    if dE > 0.05, fprintf('  %.2f MeV: nearest machine layer %.2f MeV\n', energies(k), allE(ix)); end
    pln.propStf.energy = allE(ix);
    stf = matRad_generateStf(ct, cst, pln);
    t0 = tic;
    res = matRad_calcDoseForward(ct, cst, stf, pln, 1);
    cd(rootDir);
    D = res.physicalDose;
    if ~isequal(size(D), ct.cubeDim), error('Dose not on the CT grid.'); end
    idd60  = zeros(Ny, 1); iddAll = squeeze(sum(sum(D, 2), 3)) * voxArea;
    for iy = 1:Ny
        sl = squeeze(D(iy, :, :));
        idd60(iy) = sum(sl(inDet)) * voxArea;
    end
    keep = depth >= 0;
    n = nnz(keep);
    T = [T; table(repmat(energies(k), n, 1), depth(keep), idd60(keep), iddAll(keep), ...
        'VariableNames', {'energy_MeV', 'depth_mm', 'idd_r60', 'idd_all'})]; %#ok<AGROW>
    [~, ipk] = max(idd60(keep)); dk = depth(keep);
    fprintf('%3d/%d  %7.2f MeV  peak at %6.1f mm  (%.1f s)\n', k, numel(energies), energies(k), dk(ipk), toc(t0));
    writetable(T, outFile);                 % save as we go
end
fprintf('\nSaved %s\n', outFile);
