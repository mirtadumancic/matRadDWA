%% matRad_LETd_vsDepth_check.m
% Dose-averaged LET (LETd) versus depth for single DWA spots in water:
% matRad pencil beam (LET from the machine table) vs matRad MCsquare (LET
% scored in transport). Purpose: find why the LET-volume histogram outside
% the target differs so much between matRad and RayStation (water phantom,
% centre 5 cm target), while the CTV LVH agrees.
%
% What the plan files already show (DWA_proton_v2, centre target):
%   - matRad assigns LETd to 1.4 % of the grid but dose to 6.9 %;
%     RayStation assigns LETd to ~23 % of Phantom - CTV.
%   - matRad LETd never exceeds ~12-14 keV/um: the machine LET table (the
%     Generic stopgap) peaks at ~14 keV/um just beyond the Bragg peak and is
%     held flat beyond it. RayStation's LVH tail goes past 25 keV/um.
%   - The pencil beam applies the same LET at all lateral positions of a spot
%     (LET depends on depth only); Monte Carlo LETd rises in the penumbra.
% This script measures each of these effects per depth.
%
% Geometry: as the Bragg-peak check (single spot, gantry 0, water surface at
% the isocentre), but 1 mm voxels in all directions to resolve the penumbra.
%
% Output: DWA_final_machine/LETd_vsDepth/LETd_<machine>_<engine>.csv, columns
%   energy_MeV, depth_mm,
%   D_axis, LET_axis           central 1 mm column (beam axis)
%   D_x5, LET_x5 ... D_x20     columns 5, 10, 20 mm off axis (penumbra / halo)
%   LETd_r60                   sum(D*LET)/sum(D) over r <= 60 mm (IDD-weighted)
%   nDose, nLET                voxels in the r <= 60 mm slice with dose > 0 / LET > 0
%   nDose1pct                  voxels with dose > 1 % of the slice maximum of the spot
% Run HongPB on the Mac, MCsquare on the PC (from X:\, path without spaces):
%   matRad_rc; run('DWA_final_machine/matRad_LETd_vsDepth_check.m')

clear; clc;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_cfg = matRad_rc;
matRad_cfg.disableGUI = true;

%% ---------------- settings ------------------------------------------------
machineName  = 'DWA_proton_v3';
engine       = 'HongPB';          % 'HongPB' (Mac) or 'MCsquare' (PC)
energies     = [70.0 110.0 150.0 200.0];   % nearest machine layers are used
nHistMC      = 5e6;               % LET in low-dose regions needs more histories than dose
offsets_mm   = [5 10 20];         % lateral sample positions
detRadius_mm = 60;
waterDepth_mm = 330; latHalf_mm = 40;      % 200 MeV range ~ 260 mm
res_mm       = 1;  airGap_mm = 10;

outDir = fullfile(rootDir, 'DWA_final_machine', 'LETd_vsDepth');
if ~exist(outDir, 'dir'), mkdir(outDir); end
outFile = fullfile(outDir, sprintf('LETd_%s_%s.csv', machineName, engine));
if strcmp(engine, 'MCsquare') && any(rootDir == ' ')
    error('MCsquare fails on paths with spaces: run from X:\\matRadDWA_LET (subst).');
end

%% ---------------- water box: beam along +y, surface at the isocentre ------
Nx = round(2*latHalf_mm/res_mm); Nz = Nx;
Ny = round((airGap_mm + waterDepth_mm)/res_mm);
ct = struct();
ct.resolution = struct('x', res_mm, 'y', res_mm, 'z', res_mm);
ct.cubeDim = [Ny Nx Nz];
ct.numOfCtScen = 1;
ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
iyW = round(airGap_mm/res_mm) + 1;
water = false(ct.cubeDim); water(iyW:end, :, :) = true;
rED = 0.0012*ones(ct.cubeDim); rED(water) = 1;
HU  = -1000*ones(ct.cubeDim);  HU(water) = 0;
ct.cube = {rED}; ct.cubeHU = {HU};
ct.hlut = [-1024 0.00; 0 1.00; 3000 2.5];
surfaceY = ct.y(iyW) - res_mm/2;
% isocentre on a voxel edge would put the beam axis between columns: put it on a voxel centre
iso = [ct.x(round(Nx/2)), surfaceY, ct.z(round(Nz/2))];

[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
tgt = water & abs(Yg - (surfaceY + 100)) <= 2 & abs(Xg - iso(1)) <= 2 & abs(Zg - iso(3)) <= 2;
cst = cell(2, 6);
cst(1,:) = {0, 'TARGET', 'TARGET', {find(tgt)},   struct('Priority', 1, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1), {}};
cst(2,:) = {1, 'WATER',  'OAR',    {find(water)}, struct('Priority', 2, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1), {}};
clear Yg Xg Zg tgt
ix0 = round(Nx/2); iz0 = round(Nz/2);                      % beam-axis column
[XX, ZZ] = ndgrid(ct.x - iso(1), ct.z - iso(3));
inDet = XX.^2 + ZZ.^2 <= detRadius_mm^2;
depth = ct.y(:) - surfaceY;

%% ---------------- plan ----------------------------------------------------
pln = struct();
pln.radiationMode = 'protons';  pln.machine = machineName;
pln.bioModel = 'constRBE';      pln.multScen = 'nomScen';  pln.numOfFractions = 1;
pln.propStf.gantryAngles = 0;   pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = 5;     pln.propStf.isoCenter = iso;
pln.propStf.generator = 'ParticleSingleSpot';
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
pln.propDoseCalc.calcLET = 1;
if strcmp(engine, 'HongPB')
    pln.propDoseCalc.engine = 'HongPB';
    pln.propDoseCalc.lateralModel = 'double';
    pln.propDoseCalc.geometricLateralCutOff = 100;
    pln.propDoseCalc.dosimetricLateralCutOff = 1;
else
    pln.propDoseCalc.engine = 'MCsquare';
    pln.propDoseCalc.numHistoriesDirect = nHistMC;
end
pln.propOpt.runDAO = 0; pln.propSeq.runSequencing = 0;
machine = matRad_loadMachine(pln);
allE = [machine.data.energy];

%% ---------------- loop ----------------------------------------------------
T = table();
for k = 1:numel(energies)
    [~, ie] = min(abs(allE - energies(k)));
    pln.propStf.energy = allE(ie);
    stf = matRad_generateStf(ct, cst, pln);
    t0 = tic;
    res = matRad_calcDoseForward(ct, cst, stf, pln, 1);
    cd(rootDir);
    if ~isfield(res, 'LET')
        error('No LET in the result. Fields: %s', strjoin(fieldnames(res), ', '));
    end
    D = res.physicalDose; L = res.LET; L(~isfinite(L)) = 0;

    n = Ny; z = depth;
    col = @(A, dx) squeeze(A(:, ix0 + dx, iz0));
    row = table(repmat(allE(ie), n, 1), z, col(D, 0), col(L, 0), ...
        'VariableNames', {'energy_MeV', 'depth_mm', 'D_axis', 'LET_axis'});
    for o = offsets_mm
        row.(sprintf('D_x%d', o))   = col(D, round(o/res_mm));
        row.(sprintf('LET_x%d', o)) = col(L, round(o/res_mm));
    end
    LETd = zeros(n, 1); nD = zeros(n, 1); nL = zeros(n, 1); nD1 = zeros(n, 1);
    for iy = 1:n
        d = squeeze(D(iy, :, :)); l = squeeze(L(iy, :, :));
        d = d(inDet); l = l(inDet);
        if sum(d) > 0, LETd(iy) = sum(d .* l) / sum(d); end
        nD(iy) = nnz(d > 0); nL(iy) = nnz(l > 0); nD1(iy) = nnz(d > 0.01 * max(d));
    end
    row.LETd_r60 = LETd; row.nDose = nD; row.nLET = nL; row.nDose1pct = nD1;
    row = row(z >= 0, :);
    T = [T; row]; %#ok<AGROW>
    [~, ipk] = max(row.D_axis);
    fprintf('%7.2f MeV: peak %.0f mm, LETd(r60) at peak %.2f, max LETd(r60) %.2f, max LET on axis %.2f keV/um (%.0f s)\n', ...
        allE(ie), row.depth_mm(ipk), row.LETd_r60(ipk), max(row.LETd_r60), max(row.LET_axis), toc(t0));
    writetable(T, outFile);
end
fprintf('\nSaved %s\n', outFile);
