%% matRad_makeMachines_v3.m
% Builds DWA_proton_v3 from DWA_proton_v2 and IBA_cyclotron_v2 from
% IBA_cyclotron. Parents are not modified (new names for new versions).
% Run it from anywhere: it switches to the matRadDWA root itself.
%
% WHAT CHANGES (and only this)
% (1) meta.fitAirOffset = 400 mm (= meta.BAMStoIsoDist).
%     The TOPAS depth-dose curves were scored with 400 mm of air between the
%     nozzle and the water surface (surface at the isocentre). matRad's pencil
%     beam adds 0.0011 * (nozzleToSkin - fitAirOffset) mm of WEPL for air.
%     With fitAirOffset = 0 it counted that air twice: ranges were 0.44 mm short
%     (Bragg-peak check, Nov 6 2026). With 400, the commissioning geometry gets
%     zero correction and other SSDs get only the difference.
%     MCsquare does not use this (energySpectrum is given, see below).
% (2) Absolute dose calibration of Z (MeV cm^2 / g per primary, as Generic).
%     v2 / v1 Z tables are relative (max = 100). Each layer is scaled so the
%     energy it deposits, integral(Z dz), equals f(E) * E, where f(E) is the
%     fraction of the initial energy deposited in water by matRad's Generic
%     proton machine at the same energy (0.91 at 236 MeV ... 0.99 at 50 MeV;
%     the rest leaves as neutrons and gammas). Below Generic's lowest energy
%     (31.7 MeV) f is held at its first value, capped at 1.
%     Cross-check printed below: the calibrated entrance dose at 2 mm is
%     compared with Generic's at the same energy (expected within ~1 %).
%     Pencil-beam plan dose distributions do not change (only dose per unit
%     weight). MCsquare dose is absolute (protons), so PB and MC can now be
%     compared in Gy.
%
% UNCHANGED: depths, Z shape, peakPos, offset, sigma1/sigma2/weight,
% initFocus, energySpectrum, LET (stopgap), energies, interpolated layers.
% energySpectrum is already Alaina's source (DWA: achieved energy within
% 0.005 MeV, dE/E_nozzle exact); the script only checks that it is present.
%
% Output: protons_<newName>.mat next to the parent machine file, plus
% DWA_final_machine/QA_machines_v3/<newName>_calibration.csv and .png

clear; clc;
rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
matRad_rc;

%% ---------------- options -------------------------------------------------
build = { ...  parent            new name
    'DWA_proton_v2',  'DWA_proton_v3'; ...
    'IBA_cyclotron',  'IBA_cyclotron_v2'};
fitAirOffset_mm = 400;
qaDir = fullfile(rootDir, 'DWA_final_machine', 'QA_machines_v3');
if ~exist(qaDir, 'dir'), mkdir(qaDir); end

%% ---------------- Generic: deposited-energy fraction f(E) -----------------
G  = loadMachineFile(rootDir, 'Generic');
gE = arrayfun(@(g) double(g.energy), G.data(:));
gF = arrayfun(@(g) trapz(double(g.depths(:)), double(g.Z(:))) / 10, G.data(:)) ./ gE;  % MeV / MeV
gEnt = arrayfun(@(g) interp1(double(g.depths(:)), double(g.Z(:)), 2), G.data(:));    % Z at 2 mm
[gE, o] = sort(gE); gF = gF(o); gEnt = gEnt(o);
fOfE = @(E) min(interp1(gE, gF, min(max(E, gE(1)), gE(end))), 1);

%% ---------------- build ----------------------------------------------------
for b = 1:size(build, 1)
    [m, parentFile] = loadMachineFile(rootDir, build{b,1});
    newName = build{b,2};
    nL = numel(m.data);
    qa = zeros(nL, 7);   % E, peakPos, scale, f(E), Zmax, entrance 2 mm, Generic entrance 2 mm
    for i = 1:nL
        d = m.data(i);
        E = double(d.energy);
        z = double(d.depths(:)); Z = double(d.Z(:));
        if ~isfield(d, 'energySpectrum') || isempty(d.energySpectrum)
            warning('%s layer %.2f MeV has no energySpectrum: MCsquare will fit one.', build{b,1}, E);
        end
        Edep = trapz(z, Z) / 10;                 % depths in mm -> cm
        c = fOfE(E) * E / Edep;
        m.data(i).Z = reshape(Z * c, size(d.Z));
        gEi = NaN; if E >= gE(1), gEi = interp1(gE, gEnt, E); end
        qa(i,:) = [E double(d.peakPos) c fOfE(E) max(Z)*c interp1(z, Z*c, 2) gEi];
    end

    % meta
    m.meta.machine      = newName;
    m.meta.parent       = build{b,1};
    m.meta.fitAirOffset = fitAirOffset_mm;
    m.meta.created_on   = datestr(now, 'dd-mmm-yyyy');
    if ~isfield(m.meta, 'created_by') || isempty(m.meta.created_by), m.meta.created_by = 'Mirta Dumancic'; end
    old = ''; if isfield(m.meta, 'description'), old = m.meta.description; end
    m.meta.description = [newName ', built from ' build{b,1} ' by matRad_makeMachines_v3.m. ' ...
        'CHANGES: (1) meta.fitAirOffset = 400 mm (TOPAS data include the 400 mm nozzle-to-surface air; ' ...
        'fitAirOffset = 0 added 0.44 mm of WEPL twice in the pencil beam). ' ...
        '(2) Z calibrated to MeV cm^2/g per primary: integral(Z dz) = f(E)*E, with f(E) the deposited-energy ' ...
        'fraction of matRad Generic at the same energy. UNCHANGED: everything else. ' ...
        'Parent description: ' old];

    % save next to the parent
    machine = m; %#ok<NASGU>
    outFile = fullfile(fileparts(parentFile), ['protons_' newName '.mat']);
    if exist(outFile, 'file'), error('%s exists already. Delete or rename it first.', outFile); end
    save(outFile, 'machine');
    fprintf('\nSaved %s\n', outFile);

    % QA
    T = array2table(qa, 'VariableNames', {'energy_MeV','peakPos_mm','scale','fDeposited', ...
        'Zmax_MeVcm2g','Z2mm_MeVcm2g','Z2mm_Generic_sameE'});
    T.entranceRatio_vsGeneric = T.Z2mm_MeVcm2g ./ T.Z2mm_Generic_sameE;
    writetable(T, fullfile(qaDir, [newName '_calibration.csv']));
    ok = isfinite(T.entranceRatio_vsGeneric) & T.energy_MeV >= 50;
    fprintf('%s: entrance dose at 2 mm / Generic (same E, >= 50 MeV): median %.3f, range %.3f-%.3f\n', ...
        newName, median(T.entranceRatio_vsGeneric(ok)), min(T.entranceRatio_vsGeneric(ok)), max(T.entranceRatio_vsGeneric(ok)));

    fig = figure('Color', 'w', 'Position', [60 60 1100 420]);
    subplot(1,2,1); hold on; box on;
    showE = [70 100 150 200 225];
    cols = lines(numel(showE));
    for q = 1:numel(showE)
        [~, i] = min(abs(qa(:,1) - showE(q)));
        [~, k] = min(abs(gE - qa(i,1)));
        plot(m.data(i).depths, m.data(i).Z, '-', 'Color', cols(q,:), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s %.0f MeV', strrep(newName,'_','\_'), qa(i,1)));
        plot(G.data(o(k)).depths, G.data(o(k)).Z, '--', 'Color', cols(q,:), 'LineWidth', 0.9, ...
            'DisplayName', sprintf('Generic %.0f MeV', gE(k)));
    end
    xlabel('Depth in water (mm)'); ylabel('Z (MeV cm^2 g^{-1} per proton)'); legend('Location', 'northwest', 'FontSize', 7);
    subplot(1,2,2); hold on; box on;
    plot(T.energy_MeV, T.entranceRatio_vsGeneric, 'o-', 'MarkerSize', 3);
    yline(1, ':'); ylim([0.9 1.1]);
    xlabel('Energy (MeV)'); ylabel('Entrance dose at 2 mm / Generic (same energy)');
    sgtitle(sprintf('%s: absolute calibration check', newName), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(qaDir, [newName '_calibration.png']), 'Resolution', 200);
end
fprintf('\nDone. Parents are untouched. Next: rerun matRad_BraggPeaks_vsTOPAS_RayStation.m with\n');
fprintf('machineName = ''DWA_proton_v3'' and airOffsetCorr = true (the default): expect dR80 = 0.\n');

%% ---------------- helpers ---------------------------------------------------
function [machine, f] = loadMachineFile(rootDir, name)
    cand = {fullfile(rootDir, 'matRad', 'basedata', ['protons_' name '.mat']), ...
            fullfile(rootDir, 'userdata', 'machines', ['protons_' name '.mat'])};
    ix = find(cellfun(@(c) exist(c, 'file') == 2, cand), 1);
    if isempty(ix), error('protons_%s.mat not found in matRad/basedata or userdata/machines.', name); end
    f = cand{ix};
    S = load(f, 'machine'); machine = S.machine;
end
