%% matRad_makeMachine_DWA_v4_MCLET.m
% Builds DWA_proton_v4 = DWA_proton_v3 with the LET table replaced by LETd
% scored with MCsquare (matRad_scoreLETtable_MCsquare.m, 20 energies, 8 Nov 2026).
% DWA_proton_v3 is not modified (new name for the new version).
%
% WHAT CHANGES: machine.data(i).LET for all 200 layers. Nothing else (depth-dose
% tables, calibration, spot sizes, energy spectrum, fitAirOffset are those of v3).
%
% HOW
%  - Scored quantity: laterally integrated, dose-weighted LETd over the r <= 60 mm
%    disc, LETd(z) = sum(D*LET)/sum(D) (column LETd_r60), 1 mm depth bins, LETd in
%    every voxel with dose. This is what a depth-only pencil-beam LET table means.
%  - MCsquare definition: stopping-power method, includes the secondary particles
%    MCsquare transports (protons, deuterons, alphas). The entrance plateau is
%    therefore higher than a protons-only LETd (1.5-2.0 keV/um at high energy).
%  - Each scored curve is put on relative depth u = z / R80 (R80 from its own IDD).
%    Plateau noise (u < 0.85) is smoothed with a 5-point moving mean; the Bragg-peak
%    region is left as scored.
%  - Distal to the peak the LETd is forced to be non-decreasing (cumulative maximum), which
%    removes dips caused by the 1 mm bins at low energy (fix after the first build, 8 Nov).
%  - Beyond the end of the dose (IDD < 1 % of its maximum distal to the peak) the
%    LETd is held at its value at that point (MC LET there is noise from a handful of
%    protons; the pencil beam weights LET with dose, so this only matters for LVHs
%    without a dose cut-off).
%  - For each machine layer: R80 from its own depth-dose table, u = depths / R80,
%    linear interpolation in energy between the two neighbouring scored energies
%    (clamped outside 20.01-228.84 MeV).
% LIMITATIONS: 1 mm scoring bins resolve the LET rise poorly below ~40 MeV (range
% < 15 mm); depth-only LET (no lateral variation); MCsquare's particle set.
%
% Run on the Mac or the PC from the matRadDWA folder:
%   matRad_rc; run('DWA_final_machine/matRad_makeMachine_DWA_v4_MCLET.m')
% Output: protons_DWA_proton_v4.mat next to protons_DWA_proton_v3.mat;
%   DWA_final_machine/LETtable_MCsquare/QA_v4/ (csv + png: old vs new LET tables)

clear; clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir); matRad_rc;

parentName = 'DWA_proton_v3';  newName = 'DWA_proton_v4';
scoreFile  = fullfile(rootDir, 'DWA_final_machine', 'LETtable_MCsquare', 'LETscore_DWA_proton_v3_MCsquare.csv');
qaDir      = fullfile(rootDir, 'DWA_final_machine', 'LETtable_MCsquare', 'QA_v4');
if ~exist(qaDir, 'dir'), mkdir(qaDir); end
plateauU = 0.85;  smoothN = 5;  tailFrac = 0.01;

%% ---------------- scored curves on relative depth -------------------------
S = readtable(scoreFile);
Es = unique(S.energy_MeV);
sc = struct('E', {}, 'u', {}, 'L', {}, 'R80', {});
for k = 1:numel(Es)
    g = S(abs(S.energy_MeV - Es(k)) < 1e-6, :);
    z = g.depth_mm; I = g.IDD_r60; L = g.LETd_r60;
    R80 = r80of(z, I);
    [~, ip] = max(I);
    iEnd = ip - 1 + find(I(ip:end) < tailFrac * max(I), 1);          % end of dose
    if isempty(iEnd), iEnd = numel(z); end
    L(iEnd:end) = L(iEnd);                                             % hold beyond the dose
    L(ip:end) = cummax(L(ip:end));   % distal to the peak LETd only rises (no dip from 1 mm binning / noise)
    pl = z / R80 < plateauU;
    Ls = movmean(L, smoothN, 'Endpoints', 'shrink');
    L(pl) = Ls(pl);
    sc(k).E = Es(k); sc(k).u = z / R80; sc(k).L = L; sc(k).R80 = R80;
end
fprintf('%d scored energies, %.2f-%.2f MeV\n', numel(sc), sc(1).E, sc(end).E);

%% ---------------- new table for every layer --------------------------------
mf = fullfile(rootDir, 'matRad', 'basedata', ['protons_' parentName '.mat']);
if ~exist(mf, 'file'), mf = fullfile(rootDir, 'userdata', 'machines', ['protons_' parentName '.mat']); end
P = load(mf, 'machine'); machine = P.machine;
oldLET = {machine.data.LET};
Eall = [sc.E];
for i = 1:numel(machine.data)
    d = machine.data(i); E = double(d.energy);
    z = double(d.depths(:)); R80 = r80of(z, double(d.Z(:)));
    u = z / R80;
    kHi = find(Eall >= E, 1); if isempty(kHi), kHi = numel(Eall); end
    kLo = max(1, kHi - 1); if Eall(kHi) == E || kHi == 1, kLo = kHi; end
    Llo = interp1(sc(kLo).u, sc(kLo).L, u, 'linear', 'extrap');
    Lhi = interp1(sc(kHi).u, sc(kHi).L, u, 'linear', 'extrap');
    if kLo == kHi, t = 0; else, t = (E - Eall(kLo)) / (Eall(kHi) - Eall(kLo)); end
    t = min(max(t, 0), 1);
    Lnew = (1 - t) * Llo + t * Lhi;
    % no extrapolation beyond the scored range of u: hold the end values
    uLo = max(sc(kLo).u(1), sc(kHi).u(1)); uHi = min(sc(kLo).u(end), sc(kHi).u(end));
    Lnew(u < uLo) = Lnew(find(u >= uLo, 1));
    Lnew(u > uHi) = Lnew(find(u <= uHi, 1, 'last'));
    Lnew = max(Lnew, 0);
    machine.data(i).LET = reshape(Lnew, size(d.LET));
end

machine.meta.machine = newName;
machine.meta.parent  = parentName;
machine.meta.created_on = datestr(now, 'dd-mmm-yyyy');
machine.meta.description = [newName ', built from ' parentName ' by matRad_makeMachine_DWA_v4_MCLET.m. ' ...
    'CHANGE: LET table = MCsquare LETd (stopping-power method, protons + deuterons + alphas, ' ...
    'dose-weighted over r <= 60 mm, no dose threshold) scored at 20 energies and interpolated in z/R80. ' ...
    'Everything else as ' parentName '. Parent description: ' P.machine.meta.description];
outFile = fullfile(fileparts(mf), ['protons_' newName '.mat']);
if exist(outFile, 'file'), error('%s exists already. Delete or rename it first.', outFile); end
save(outFile, 'machine');
fprintf('Saved %s\n', outFile);

%% ---------------- QA ---------------------------------------------------------
showE = [30 70 100 150 200 228];
fig = figure('Color', 'w', 'Position', [60 60 1200 420]);
Q = table();
for q = 1:numel(showE)
    [~, i] = min(abs([machine.data.energy] - showE(q)));
    d = machine.data(i); z = d.depths(:); R80 = r80of(z, double(d.Z(:)));
    subplot(2, 3, q); hold on; box on;
    plot(z, oldLET{i}(:), 'r-', 'LineWidth', 1.1);
    plot(z, d.LET(:), 'b-', 'LineWidth', 1.1);
    xline(R80, ':'); xlim([0 R80 + 10]); ylim([0 min(20, 1.2*max(d.LET(z < R80 + 3)))]);
    title(sprintf('%.1f MeV', d.energy)); xlabel('Depth (mm)'); ylabel('LET_d (keV/\mum)');
    if q == 1, legend({parentName, newName}, 'Location', 'northwest', 'Interpreter', 'none'); end
    for f = [0.2 0.5 0.8 0.95 1.0]
        Q = [Q; table(d.energy, f, interp1(z, oldLET{i}(:), f*R80), interp1(z, d.LET(:), f*R80), ...
            'VariableNames', {'energy_MeV', 'z_over_R80', 'LET_v3', 'LET_v4'})]; %#ok<AGROW>
    end
end
writetable(Q, fullfile(qaDir, [newName '_LET_v3_vs_v4.csv']));
exportgraphics(fig, fullfile(qaDir, [newName '_LET_v3_vs_v4.png']), 'Resolution', 200);
disp(Q);

function R = r80of(z, I)
    [m, ip] = max(I);
    k = ip - 1 + find(I(ip:end) < 0.8 * m, 1);
    if isempty(k) || k <= ip, R = z(ip); return; end
    R = interp1(I([k-1 k]), z([k-1 k]), 0.8 * m);
end
