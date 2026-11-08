%% matRad_scoreLETtable_MCsquare.m
% Scores LETd versus depth with MCsquare for single DWA_proton_v3 spots in water,
% to build a machine LET table from Monte Carlo (planned DWA_proton_v4).
% Same geometry as matRad_LETd_vsDepth_check.m (1 mm voxels, surface at the
% isocentre, LETd in every voxel with dose). For the table the laterally
% integrated, dose-weighted LETd is used: LETd(z) = sum(D*LET)/sum(D) over the
% r <= 60 mm disc (column LETd_r60), which is what a depth-only pencil-beam LET
% table represents. The IDD (sum of dose over the disc) is saved too.
% MCsquare LET: stopping-power method (matRad default); it includes the
% secondary particles MCsquare transports (protons, deuterons, alphas).
%
% Energies: ~20 of the 58 TOPAS energies across the spectrum (the build script
% interpolates between them in relative depth z/R80 for all 200 machine layers).
% Run on the PC from X:\matRadDWA_LET (MCsquare; about 3-5 min per energy):
%   matRad_rc; run('DWA_final_machine/matRad_scoreLETtable_MCsquare.m')
% Output: DWA_final_machine/LETtable_MCsquare/LETscore_DWA_proton_v3_MCsquare.csv
%   (saved after every energy, so a stopped run keeps what it has)

clear; clc;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_cfg = matRad_rc;
matRad_cfg.disableGUI = true;

%% ---------------- settings ------------------------------------------------
machineName  = 'DWA_proton_v3';
engines      = {'MCsquare'};
energies     = [20.01 29.02 38.02 47.03 56.03 65.04 74.05 83.05 92.06 102.42 115.16 127.46 140.72 153.14 167.89 181.01 195.22 209.29 219.20 228.84];   % 20 of the 58 TOPAS energies
nHistMC      = 5e6;               % LET in low-dose regions needs more histories than dose
offsets_mm   = [5 10 20];         % lateral sample positions
detRadius_mm = 60;
waterDepth_mm = 380; latHalf_mm = 40;      % 220 MeV range ~ 305 mm, plus the distal tail
res_mm       = 1;  airGap_mm = 10;

outDir = fullfile(rootDir, 'DWA_final_machine', 'LETtable_MCsquare');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if any(strcmp(engines, 'MCsquare')) && any(rootDir == ' ')
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

% LETd in every voxel with dose (RayStation-like), not only above 1 % of the maximum
matRad_cfg.defaults.doseWeightingThreshold = 0;

for iEng = 1:numel(engines)
engine  = engines{iEng};
outFile = fullfile(outDir, sprintf('LETscore_%s_%s.csv', machineName, engine));
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
    pln.propDoseCalc.numHistoriesDirect = nHistMC;   % MCsquare LET: stopping-power method (matRad default, RayStation-like)
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
    LETd = zeros(n, 1); IDD = zeros(n, 1); nD = zeros(n, 1); nL = zeros(n, 1); nD1 = zeros(n, 1);
    for iy = 1:n
        d = squeeze(D(iy, :, :)); l = squeeze(L(iy, :, :));
        d = d(inDet); l = l(inDet);
        IDD(iy) = sum(d) * res_mm^2;
        if sum(d) > 0, LETd(iy) = sum(d .* l) / sum(d); end
        nD(iy) = nnz(d > 0); nL(iy) = nnz(l > 0); nD1(iy) = nnz(d > 0.01 * max(d));
    end
    row.IDD_r60 = IDD; row.LETd_r60 = LETd; row.nDose = nD; row.nLET = nL; row.nDose1pct = nD1;
    row = row(z >= 0, :);
    T = [T; row]; %#ok<AGROW>
    [~, ipk] = max(row.D_axis);
    fprintf('%7.2f MeV: peak %.0f mm, LETd(r60) at peak %.2f, max LETd(r60) %.2f, max LET on axis %.2f keV/um (%.0f s)\n', ...
        allE(ie), row.depth_mm(ipk), row.LETd_r60(ipk), max(row.LETd_r60), max(row.LET_axis), toc(t0));
    writetable(T, outFile);
end
fprintf('\nSaved %s\n', outFile);
end   % engines
matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'doseWeightingThreshold');
