%% matRad_makeDWAmachine_v2.m
% Builds the machine DWA_proton_v2 from DWA_proton_doubleGauss_interp (v1).
% Run it from anywhere: it switches to the matRadDWA root itself.
%
% WHAT CHANGES FROM v1 (and only this)
% (1) initFocus fixed. v1 stored initFocus.dist = SAD - Z, which mirrors the
%     in-air spot about isocentre (matRad's dist increases downstream, and so
%     does Alaina's Z). v2 uses dist = SAD + Z. The beam now converges
%     downstream at high energy, as the TRANSOPTR optics say. The emittance
%     for MC engines is refitted from the corrected initFocus on save.
% (2) In-water halo (sigma2, weight) now changes with depth. v1 carried the
%     in-air halo from the air profiles at a constant weight of 0.33-0.52 at
%     every depth, which gave Remo's shallow penumbra. That weight is not
%     physically plausible for a beam that passes only a monitor chamber and
%     air, and its meaning (central cut vs projection) is unresolved. v2 drops
%     it and uses the in-water nuclear halo from matRad's Generic proton
%     kernel at the SAME RANGE and the same fraction of range (z / R). The
%     Generic table is smoothed along depth and across energy first, because
%     its own fits switch between two halo families around 120-150 MeV.
%     The nuclear halo is set by water physics, not by the machine, so this
%     is a reasonable first approximation. Both machines then share one
%     in-water halo model, and a DWA vs Generic comparison isolates the
%     beam itself (spot size, convergence, depth dose).
% (3) Z corrected for halo dose outside the TOPAS scoring radius (r = 60 mm):
%     Z = Z / (1 - w * exp(-r^2 / (2 sigma2tot^2))). At most ~1% (220 MeV).
% (4) OPTIONAL (addLET = true): data(i).LET, the dose-averaged LET from
%     Generic at the same range, mapped by z / R. A STOPGAP so LET can be
%     computed and optimised before TOPAS or MCsquare LET exists. DWA has a
%     narrower energy spread than Generic, so the distal LET will differ.
%
% UNCHANGED: depths, Z shape, peakPos, sigma1 (in-water core MCS, Highland),
% initFocus sigma values, energySpectrum, energies, interpolated layers.
%
% KNOWN GAPS (still): no DWA-specific in-air halo; halo and LET borrowed from
% Generic; Z relative, not absolute. To be replaced by in-water TOPAS scoring
% or Remo's MCsquare machine.

clear; clc;
% Work from the matRadDWA root (parent of DWA_final_machine), wherever MATLAB is
rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
matRad_rc;                                   % matRad on the path

% ---------------- options -------------------------------------------------
v1Name   = 'DWA_proton_doubleGauss_interp';
v2Name   = 'DWA_proton_v2';
addLET   = true;      % stopgap LET from Generic (see header, point 4)
rDet_mm  = 60;        % TOPAS IDD scoring radius
tMax     = 1.02;      % halo table stops at 1.02 x peakPos (Generic is noisy beyond)
dT       = 0.005;     % normalised-depth step of the halo table
nSmDepth = 5;         % moving median along depth (Generic layers)
nSmE     = 9;         % moving median across Generic energies ...
nSmE2    = 5;         % ... then moving mean across energies
qaE      = [40 70 110 150 220];   % energies shown in the QA figure

% ---------------- load ----------------------------------------------------
S  = load(fullfile(rootDir,'userdata','machines',['protons_' v1Name '.mat']));
m1 = S.machine;
S  = load(fullfile(rootDir,'matRad','basedata','protons_Generic.mat'));
G  = S.machine;
SAD = double(m1.meta.SAD);

% ---------------- Generic halo table on (range, z/R) ----------------------
gR = arrayfun(@(g) double(g.peakPos), G.data);
[gR, gOrd] = sort(gR(:));
GD = G.data(gOrd);
T  = 0:dT:tMax;
W  = zeros(numel(GD), numel(T));
S2 = zeros(numel(GD), numel(T));
for k = 1:numel(GD)
    dk = double(GD(k).depths(:));
    wk = movmedian(double(GD(k).weight(:)), nSmDepth);
    sk = movmedian(double(GD(k).sigma2(:)), nSmDepth);
    W(k,:)  = interp1(dk, wk, T*gR(k), 'linear', wk(end));
    S2(k,:) = interp1(dk, sk, T*gR(k), 'linear', sk(end));
end
W  = movmean(movmedian(W,  nSmE, 1), nSmE2, 1);
S2 = movmean(movmedian(S2, nSmE, 1), nSmE2, 1);

% ---------------- build v2 layer by layer ---------------------------------
nL = numel(m1.data);
dataCell = cell(1, nL);
qa = zeros(nL, 6);   % energy, R, w(R/2), w(R), max Z corr %, sigmaIni(iso)
for i = 1:nL
    d = m1.data(i);
    z = double(d.depths(:));
    R = double(d.peakPos);

    % (1) initFocus: v1 dist = SAD - Z  ->  Z = SAD - dist  ->  v2 dist = SAD + Z
    Zair  = SAD - double(d.initFocus.dist(:));        % double(): dist is uint16 in v1
    dist2 = SAD + Zair;
    [dist2, o] = sort(dist2);
    sig = double(d.initFocus.sigma(:));
    d.initFocus.dist  = dist2(:)';
    d.initFocus.sigma = sig(o)';
    if isfield(d.initFocus, 'emittance')         % v1 emittance was fitted to the
        d.initFocus = rmfield(d.initFocus, 'emittance');   % mirrored spot: refit on save
    end
    sIni = interp1(d.initFocus.dist, d.initFocus.sigma, SAD);
    d.initFocus.SisFWHMAtIso = sIni * 2*sqrt(2*log(2));

    % (2) halo from Generic at the same range and z/R
    rq = min(max(R, gR(1)), gR(end));
    tq = min(max(z / R, 0), tMax);
    wRow = interp1(gR, W,  rq);
    sRow = interp1(gR, S2, rq);
    w  = interp1(T, wRow, tq);
    s2 = interp1(T, sRow, tq);
    s1 = double(d.sigma1(:));
    s2 = max(s2, s1 + 1);                 % halo always wider than core
    w  = min(max(w, 0), 1);

    % (3) Z: add back halo dose outside the scoring radius
    fOut = w .* exp(-rDet_mm^2 ./ (2*(s2.^2 + sIni^2)));
    Zv   = double(d.Z(:)) ./ (1 - fOut);

    d.Z      = reshape(Zv, size(d.Z));
    d.sigma2 = reshape(s2, size(d.depths));
    d.weight = reshape(w,  size(d.depths));

    % (4) LET stopgap
    if addLET
        if R <= gR(1),       a = 1;          b = 1; t = 0;
        elseif R >= gR(end), a = numel(gR);  b = a; t = 0;
        else
            a = find(gR <= R, 1, 'last'); b = a + 1;
            t = (R - gR(a)) / (gR(b) - gR(a));
        end
        La = letAt(GD(a), gR(a), R, z);
        Lb = letAt(GD(b), gR(b), R, z);
        d.LET = reshape((1-t)*La + t*Lb, size(d.depths));
    end

    dataCell{i} = d;
    [~, iH] = min(abs(z - 0.5*R)); [~, iR] = min(abs(z - R));
    qa(i,:) = [double(d.energy) R w(iH) w(iR) 100*max(fOut(z <= 1.05*R)) sIni];
end
m2      = m1;
m2.data = [dataCell{:}];

% ---------------- meta -----------------------------------------------------
m2.meta.machine     = v2Name;
m2.meta.created_on  = datestr(now, 'dd-mmm-yyyy');
m2.meta.parent      = v1Name;
m2.meta.description = [ ...
 'DWA_proton_v2, built from ' v1Name ' by matRad_makeDWAmachine_v2.m. ' ...
 'CHANGES: (1) initFocus fixed to dist = SAD + Z (v1 mirrored it about isocentre). ' ...
 '(2) In-water halo sigma2/weight now depth dependent, taken from matRad Generic at ' ...
 'the same range and z/R (Generic table smoothed in depth and energy); the in-air ' ...
 'halo of v1 (constant weight 0.33-0.52) is dropped as physically implausible and ' ...
 'unresolved. (3) Z corrected for halo dose outside the r = 60 mm TOPAS scoring ' ...
 'radius (max ~1%). (4) ' iif(addLET, ['LET (keV/um) added as a STOPGAP from Generic ' ...
 'at the same range and z/R; DWA energy spread is narrower, so distal LET is approximate. '], ...
 'No LET. ') ...
 'UNCHANGED: depths, peakPos, sigma1 (Highland core MCS), initFocus sigma values, ' ...
 'energySpectrum, 200 layers (58 native + 142 interpolated). ' ...
 'KNOWN GAPS: no DWA-specific in-air halo; halo and LET borrowed from Generic; Z relative. ' ...
 'Replace with in-water TOPAS scoring or the MCsquare-based machine when available. ' ...
 'Parent description: ' m1.meta.description];

% ---------------- QA figure ------------------------------------------------
fig = figure('Color','w','Position',[50 50 1400 720]);
tl  = {'sigma1: core MCS in water [mm] (unchanged)', ...
       'sigma2: halo width [mm] (v1 dashed, v2 solid)', ...
       'weight: halo fraction (v1 dashed, v2 solid)', ...
       'Z correction: halo dose outside r = 60 mm [%]', ...
       'LET_d (stopgap, from Generic) [keV/\mum]', ...
       'initFocus \sigma vs Z from iso (+ downstream) [mm] (v1 dashed, v2 solid)'};
ax = gobjects(1,6);
for p = 1:6, ax(p) = subplot(2,3,p); hold on; grid on; title(tl{p},'FontSize',9); end
cols = lines(numel(qaE));
for q = 1:numel(qaE)
    [~, i] = min(abs(arrayfun(@(x) double(x.energy), m1.data) - qaE(q)));
    a1 = m1.data(i); a2 = m2.data(i);
    z = double(a1.depths(:)); R = double(a1.peakPos); m = z <= 1.08*R; c = cols(q,:);
    nm = sprintf('%.0f MeV', a1.energy);
    plot(ax(1), z(m), a1.sigma1(m), '-', 'Color', c, 'DisplayName', nm);
    plot(ax(2), z(m), a1.sigma2(m), '--', 'Color', c); plot(ax(2), z(m), a2.sigma2(m), '-', 'Color', c);
    plot(ax(3), z(m), a1.weight(m), '--', 'Color', c); plot(ax(3), z(m), a2.weight(m), '-', 'Color', c);
    plot(ax(4), z(m), 100*(double(a2.Z(m))./double(a1.Z(m)) - 1), '-', 'Color', c);
    if addLET, plot(ax(5), z(m), a2.LET(m), '-', 'Color', c); end
    plot(ax(6), double(a1.initFocus.dist) - SAD, a1.initFocus.sigma, '--o', 'Color', c, 'MarkerSize', 3);
    plot(ax(6), a2.initFocus.dist - SAD, a2.initFocus.sigma, '-o', 'Color', c, 'MarkerSize', 3);
end
for p = 1:5, xlabel(ax(p), 'depth in water [mm]'); end
xlabel(ax(6), 'Z from isocentre [mm]'); legend(ax(1), 'Location', 'northwest');
sgtitle(sprintf('%s vs %s', v2Name, v1Name), 'Interpreter', 'none');
qaDir = fullfile(rootDir, 'DWA_final_machine', ['QA_' v2Name]);
if ~exist(qaDir, 'dir'), mkdir(qaDir); end
exportgraphics(fig, fullfile(qaDir, [v2Name '_QA.png']), 'Resolution', 150);
writetable(array2table(qa, 'VariableNames', ...
    {'energy_MeV','peakPos_mm','weight_halfRange','weight_atPeak','maxZcorr_pct','sigmaIni_iso_mm'}), ...
    fullfile(qaDir, [v2Name '_QA.csv']));
fprintf('weight at R/2: %.3f (%.0f MeV) ... %.3f (%.0f MeV)\n', qa(1,3), qa(1,1), qa(end,3), qa(end,1));
fprintf('max Z correction: %.2f %%\n', max(qa(:,5)));

% ---------------- save (refits MC emittance from the fixed initFocus) -----
machine = m2;
obj = matRad_MCemittanceBaseData(machine);
obj.saveMatradMachine(v2Name);
fprintf('\nSaved protons_%s.mat. Use pln.machine = ''%s''.\n', v2Name, v2Name);
fprintf('v1 (%s) is untouched for comparison.\n', v1Name);

% ---------------- local functions -----------------------------------------
function L = letAt(g, Rg, R, z)
    dg = double(g.depths(:)); Lg = double(g.LET(:));
    L  = interp1(dg, Lg, z * Rg / R, 'linear', Lg(end));
end
function out = iif(c, a, b)
    if c, out = a; else, out = b; end
end
