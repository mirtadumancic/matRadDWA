%% matRad_commissionIBAmachine.m
% Commissions the IBA cyclotron (RayStation machine GTR102_AddedRS) into
% matRad as the machine 'IBA_cyclotron', with the same recipe as
% DWA_proton_v2, so the two can be compared on equal terms.
%
% INPUT (TOPAS commissioning data in RayStation export format, as for the DWA),
% in matRadDWA/IBA_data/:
%   Measured_GTR102_AddedRS_PristineBraggPeaks.csv  30 IDDs in water,
%       99.59-226.35 MeV, detector diameter 81.6 mm, surface at isocentre
%   Measured_GTR102_AddedRS_SpotProfiles.csv        300 air profiles, X and Y
%       at Z = -200, -100, 0, +100, +200 mm (+Z downstream), snout-iso 400 mm
%   (Measured_GTR102_AddedRS_AbsoluteDosimetry.csv  dose/MU; not used by
%       matRad's pencil-beam engine, kept for reference)
%
% RECIPE (identical to DWA_proton_v2):
%  - Z: the measured IDDs, prepended with depth 0 (linear extrapolation of
%    the first two points), then corrected for halo dose outside the
%    detector radius (40.8 mm).
%  - initFocus: core sigma of the air profiles (X and Y averaged) vs
%    distance, dist = SAD + Z (correct orientation, no mirroring).
%  - sigma1: in-water MCS core broadening (matRad's Highland form).
%  - sigma2, weight: in-water nuclear halo from matRad Generic at the same
%    range and z/R (Generic table smoothed in depth and energy). The air
%    profiles here only extend ~±18 mm, too narrow to fit an air halo.
%  - LET: range-matched from Generic (stopgap, as for DWA_proton_v2).
%  - Energy layers: 30 native + interpolated layers so the range step is at
%    most 0.5 x the local 80-80 peak width (same rule as the DWA).
%  - energySpectrum: dE/E linear from 0.67 % (100 MeV) to 0.27 % (226 MeV),
%    the values quoted for this beam library in Bui et al 2026. Used only by
%    MC engines; the pencil-beam dose does not read it.
%
% Run from anywhere: the script switches to the matRadDWA root itself.

clear; clc;
rootDir = fileparts(mfilename('fullpath'));                 % works from the root or a subfolder
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_rc;

%% ----------------------- settings ---------------------------------------
inDir     = fullfile(rootDir, 'IBA_data');
braggFile = fullfile(inDir, 'Measured_GTR102_AddedRS_PristineBraggPeaks.csv');
spotFile  = fullfile(inDir, 'Measured_GTR102_AddedRS_SpotProfiles.csv');

machineName   = 'IBA_cyclotron';
isoToSnout_mm = 400;
SAD_mm        = 10000;      % numerical convention, as for the DWA machines
rDet_mm       = 81.6/2;     % IDD detector radius

energyMatchTol_MeV = 0.05;
coreThreshold  = 0.10;
haloCoreSigmas = 3.0;

doInterpolate         = true;
interpSpacingFraction = 0.5;
interpMinStep_mm      = 1.0;

addLET   = true;
tMax     = 1.02;  dT = 0.005;  nSmDepth = 5;  nSmE = 9;  nSmE2 = 5;

spreadE_MeV = [100 226.35];
spreadPct   = [0.67 0.27];
spreadModel = @(E) min(max(interp1(spreadE_MeV, spreadPct, E, 'linear', 'extrap'), 0.2), 0.8);

%% ----------------------- STEP 1: depth-dose curves ----------------------
fprintf('Parsing %s ...\n', braggFile);
braggBlocks = parseRayStationBlocks(braggFile);
nE = numel(braggBlocks);
braggEnergy = nan(nE,1); braggDepths = cell(nE,1); braggZ = cell(nE,1);
for k = 1:nE
    b = braggBlocks(k);
    braggEnergy(k) = str2double(b.meta('Nominal beam energy [MeV/A]:'));
    d = vertcat(b.data{:});
    d = d(~any(isnan(d),2), :);
    [z, o] = sort(d(:,1)); Zk = d(o,2);
    [z, iu] = unique(z); Zk = Zk(iu);
    if z(1) > 0                                   % add depth 0
        Z0 = Zk(1) - (Zk(2)-Zk(1)) / (z(2)-z(1)) * z(1);
        z = [0; z]; Zk = [max(Z0, 0); Zk];
    end
    braggDepths{k} = z; braggZ{k} = Zk;
end
[braggEnergy, o] = sort(braggEnergy);
braggDepths = braggDepths(o); braggZ = braggZ(o);
fprintf('  %d energies, %.2f - %.2f MeV\n', nE, braggEnergy(1), braggEnergy(end));

%% ----------------------- STEP 2: air profiles -> core sigma -------------
fprintf('Parsing %s ...\n', spotFile);
spotBlocks = parseRayStationBlocks(spotFile);
nSpot = numel(spotBlocks);
spotEnergy = nan(nSpot,1); spotZpos = nan(nSpot,1); spotCore = nan(nSpot,1);
for k = 1:nSpot
    b = spotBlocks(k);
    spotEnergy(k) = str2double(b.meta('Nominal beam energy [MeV/A]:'));
    ax            = string(b.meta('Curve type [Depth/X/Y/AbsoluteDosimetry]:'));
    spotZpos(k)   = str2double(b.meta('Air profile Z pos [mm]:'));
    d = vertcat(b.data{:}); d = d(~any(isnan(d),2), :);
    if ax == "X", pos = d(:,1); else, pos = d(:,2); end
    spotCore(k) = fitCoreAndHalo(pos, d(:,3), coreThreshold, haloCoreSigmas);
end
fprintf('  %d profiles, %d core fits failed\n', nSpot, sum(~isfinite(spotCore)));
uniqueZ = unique(spotZpos);
isoIdx  = find(uniqueZ == 0, 1);
if isempty(isoIdx), error('No air profile at Z = 0.'); end

%% ----------------------- STEP 3: native layers --------------------------
data = struct([]);
for i = 1:nE
    E = braggEnergy(i);
    z = braggDepths{i}; Zi = braggZ{i};
    [~, iM] = max(Zi);
    coreVsZ = nan(numel(uniqueZ),1);
    for j = 1:numel(uniqueZ)
        m = abs(spotEnergy - E) < energyMatchTol_MeV & spotZpos == uniqueZ(j);
        coreVsZ(j) = mean(spotCore(m), 'omitnan');
    end
    s = struct();
    s.energy  = E;
    s.depths  = z;
    s.Z       = Zi;
    s.offset  = 0;
    s.peakPos = z(iM);
    s.initFocus.dist  = (SAD_mm + uniqueZ(:))';      % +Z downstream -> larger dist
    s.initFocus.sigma = coreVsZ(:)';
    s.initFocus.unit  = 'mm';
    s.initFocus.SisFWHMAtIso = coreVsZ(isoIdx) * 2*sqrt(2*log(2));
    s.sigma1 = matRad_calcSigmaLatMCS_water(z, E);
    s.sigma2 = s.sigma1 + 1;                          % placeholder, set in STEP 5
    s.weight = zeros(size(z));                        % placeholder, set in STEP 5
    s.energySpectrum.type  = 'gaussian';
    s.energySpectrum.mean  = E;
    s.energySpectrum.sigma = spreadModel(E);          % percent
    s.interpolated = false;
    if isempty(data), data = s; else, data(end+1) = s; end %#ok<AGROW>
end

%% ----------------------- STEP 4: interpolated layers --------------------
R80native = nan(nE,1); W80native = nan(nE,1);
for i = 1:nE, [R80native(i), W80native(i)] = braggMetrics(data(i).depths, data(i).Z); end
okW = isfinite(W80native);
pW  = polyfit(log(R80native(okW)), log(W80native(okW)), 1);
widthModel = @(R) exp(pW(2)) .* R.^pW(1);
fprintf('  80-80 width model: W80 = %.4f * R80^%.3f mm\n', exp(pW(2)), pW(1));
if doInterpolate
    newData = data(1);
    for i = 1:nE-1
        Ra = R80native(i); Rb = R80native(i+1);
        step = max(interpSpacingFraction * widthModel(0.5*(Ra+Rb)), interpMinStep_mm);
        nIns = ceil((Rb - Ra) / step) - 1;
        for j = 1:nIns
            newData(end+1) = interpolateLayer(data(i), data(i+1), Ra, Rb, Ra + (Rb-Ra)*j/(nIns+1)); %#ok<SAGROW>
        end
        newData(end+1) = data(i+1); %#ok<SAGROW>
    end
    allE = [newData.energy];
    if numel(unique(allE)) ~= numel(allE) || any(diff(allE) <= 0)
        error('Interpolated energies are not unique/ascending.');
    end
    data = newData;
end
fprintf('  %d layers (%d native + %d interpolated)\n', numel(data), nE, numel(data) - nE);

%% ----------------------- STEP 5: halo, Z correction, LET from Generic ---
S = load(fullfile(rootDir, 'matRad', 'basedata', 'protons_Generic.mat'));  G = S.machine;
gR = arrayfun(@(g) double(g.peakPos), G.data);
[gR, gOrd] = sort(gR(:)); GD = G.data(gOrd);
T  = 0:dT:tMax;
W  = zeros(numel(GD), numel(T)); S2 = zeros(numel(GD), numel(T));
for k = 1:numel(GD)
    dk = double(GD(k).depths(:));
    wk = movmedian(double(GD(k).weight(:)), nSmDepth);
    sk = movmedian(double(GD(k).sigma2(:)), nSmDepth);
    W(k,:)  = interp1(dk, wk, T*gR(k), 'linear', wk(end));
    S2(k,:) = interp1(dk, sk, T*gR(k), 'linear', sk(end));
end
W  = movmean(movmedian(W,  nSmE, 1), nSmE2, 1);
S2 = movmean(movmedian(S2, nSmE, 1), nSmE2, 1);

nL = numel(data); dataCell = cell(1, nL);
qa = zeros(nL, 7);   % E, R80, peakPos, sigmaIni_iso, w(R/2), w(R), maxZcorr%
for i = 1:nL
    d = data(i); z = double(d.depths(:)); R = double(d.peakPos);
    sIni = interp1(d.initFocus.dist, d.initFocus.sigma, SAD_mm);
    rq = min(max(R, gR(1)), gR(end));  tq = min(max(z / R, 0), tMax);
    w  = interp1(T, interp1(gR, W,  rq), tq);
    s2 = interp1(T, interp1(gR, S2, rq), tq);
    s2 = max(s2, d.sigma1(:) + 1);  w = min(max(w, 0), 1);
    fOut = w .* exp(-rDet_mm^2 ./ (2*(s2.^2 + sIni^2)));
    d.Z      = double(d.Z(:)) ./ (1 - fOut);
    d.sigma2 = s2;
    d.weight = w;
    if addLET
        if R <= gR(1),       a = 1;         b = 1; t = 0;
        elseif R >= gR(end), a = numel(gR); b = a; t = 0;
        else, a = find(gR <= R, 1, 'last'); b = a + 1; t = (R - gR(a)) / (gR(b) - gR(a));
        end
        d.LET = (1-t) * letAt(GD(a), gR(a), R, z) + t * letAt(GD(b), gR(b), R, z);
    end
    dataCell{i} = d;
    [~, iH] = min(abs(z - 0.5*R)); [~, iR] = min(abs(z - R));
    qa(i,:) = [d.energy braggMetrics(z, d.Z) R sIni w(iH) w(iR) 100*max(fOut(z <= 1.05*R))];
end
data = [dataCell{:}];

%% ----------------------- STEP 6: meta -----------------------------------
meta.radiationMode = 'protons';
meta.machine       = machineName;
meta.SAD           = SAD_mm;
meta.BAMStoIsoDist = isoToSnout_mm;
meta.MCcode        = 'TOPAS';
meta.dataType      = 'doubleGauss';
meta.fitAirOffset  = 0;
meta.LUTspotSize   = [0 100; 0 0];
meta.created_on    = datestr(now, 'dd-mmm-yyyy');
meta.interpolation.nNative        = nE;
meta.interpolation.nTotal         = numel(data);
meta.interpolation.nativeEnergies = braggEnergy(:)';
meta.description = [ ...
 'IBA cyclotron proton machine (TOPAS commissioning data in RayStation format, GTR102_AddedRS: 30 IDDs ' ...
 '99.59-226.35 MeV, detector 81.6 mm; air profiles at Z = -200..+200 mm, snout-iso 400 mm), ' ...
 'built by matRad_commissionIBAmachine.m with the same recipe as DWA_proton_v2: initFocus = ' ...
 'core of the air profiles (X/Y averaged), dist = SAD + Z; sigma1 = Highland MCS in water; ' ...
 'sigma2/weight = in-water halo from matRad Generic at the same range and z/R; Z corrected for ' ...
 'halo dose outside r = 40.8 mm; LET range-matched from Generic (stopgap); interpolated layers ' ...
 'to a range step <= 0.5 x the 80-80 width. energySpectrum.sigma (percent) linear 0.67 % at ' ...
 '100 MeV to 0.27 % at 226 MeV (Bui et al 2026), approximate, MC engines only. Absolute ' ...
 'dose/MU data not used. KNOWN GAPS: halo and LET borrowed from Generic; Z relative.'];

machine.meta = meta;
machine.data = data;

%% ----------------------- STEP 7: QA -------------------------------------
qaDir = fullfile(rootDir, 'DWA_final_machine', ['QA_' machineName]);
if ~exist(qaDir, 'dir'), mkdir(qaDir); end
isInt = [data.interpolated];
fig = figure('Color', 'w', 'Position', [50 50 1400 720]);
subplot(2,3,1); hold on; grid on;
plot(qa(~isInt,1), qa(~isInt,2), 'ko', 'MarkerFaceColor', 'k'); plot(qa(isInt,1), qa(isInt,2), 'r.');
xlabel('energy [MeV]'); ylabel('R80 [mm]'); title('Range: native (black), interpolated (red)');
subplot(2,3,2); hold on; grid on;
plot(qa(:,1), qa(:,4), '.-'); xlabel('energy [MeV]'); ylabel('\sigma in air at iso [mm]');
title('Spot size in air at isocentre (core)');
subplot(2,3,3); hold on; grid on;
qE = [100 130 160 190 226]; cols = lines(numel(qE));
for q = 1:numel(qE)
    [~, i] = min(abs([data.energy] - qE(q)));
    plot(data(i).initFocus.dist - SAD_mm, data(i).initFocus.sigma, '-o', 'Color', cols(q,:), ...
        'DisplayName', sprintf('%.0f MeV', data(i).energy));
end
xlabel('Z from isocentre (+ downstream) [mm]'); ylabel('\sigma [mm]'); title('initFocus'); legend('Location', 'best');
subplot(2,3,4); hold on; grid on;
for q = 1:numel(qE)
    [~, i] = min(abs([data.energy] - qE(q))); m = data(i).depths <= 1.08*data(i).peakPos;
    plot(data(i).depths(m), data(i).Z(m) / max(data(i).Z), 'Color', cols(q,:));
end
xlabel('depth [mm]'); ylabel('Z (normalised)'); title('Depth dose');
subplot(2,3,5); hold on; grid on;
for q = 1:numel(qE)
    [~, i] = min(abs([data.energy] - qE(q))); m = data(i).depths <= 1.08*data(i).peakPos;
    plot(data(i).depths(m), data(i).weight(m), 'Color', cols(q,:));
end
xlabel('depth [mm]'); ylabel('halo weight'); title('Halo weight (from Generic)');
subplot(2,3,6); hold on; grid on;
if addLET
    for q = 1:numel(qE)
        [~, i] = min(abs([data.energy] - qE(q))); m = data(i).depths <= 1.08*data(i).peakPos;
        plot(data(i).depths(m), data(i).LET(m), 'Color', cols(q,:));
    end
end
xlabel('depth [mm]'); ylabel('LET_d [keV/\mum]'); title('LET (stopgap, from Generic)');
sgtitle(machineName, 'Interpreter', 'none');
exportgraphics(fig, fullfile(qaDir, [machineName '_QA.png']), 'Resolution', 150);
writetable(array2table(qa, 'VariableNames', {'energy_MeV','R80_mm','peakPos_mm', ...
    'sigmaIni_iso_mm','weight_halfRange','weight_atPeak','maxZcorr_pct'}), ...
    fullfile(qaDir, [machineName '_QA.csv']));
fprintf('spot sigma in air at iso: %.2f mm (%.0f MeV) ... %.2f mm (%.0f MeV)\n', ...
    qa(1,4), qa(1,1), qa(end,4), qa(end,1));
fprintf('max Z correction: %.2f %%\n', max(qa(:,7)));

%% ----------------------- STEP 8: save (fits MC emittance) ---------------
obj = matRad_MCemittanceBaseData(machine);
obj.saveMatradMachine(machineName);
copyfile(fullfile(rootDir, 'userdata', 'machines', ['protons_' machineName '.mat']), ...
         fullfile(rootDir, 'matRad', 'basedata', ['protons_' machineName '.mat']));
fprintf('\nSaved protons_%s.mat (userdata/machines and matRad/basedata). Use pln.machine = ''%s''.\n', ...
    machineName, machineName);

%% ----------------------- local functions --------------------------------
function L = letAt(g, Rg, R, z)
    dg = double(g.depths(:)); Lg = double(g.LET(:));
    L  = interp1(dg, Lg, z * Rg / R, 'linear', Lg(end));
end

function s = interpolateLayer(a, b, Ra, Rb, Rt)
% New layer with R80 = Rt between native layers a and b (same rule as the
% DWA machine): depth curves stretched so both ranges land on Rt, then
% blended; energy from a local power law; initFocus and energy spread blended
% linearly; in-water MCS recomputed. Halo and LET are set later (STEP 5).
    t = (Rt - Ra) / (Rb - Ra);
    s = a;
    p = log(Rb/Ra) / log(b.energy/a.energy);
    E = round(a.energy * (Rt/Ra)^(1/p) * 100) / 100;
    s.energy = E;
    z  = (0:0.5:max(b.depths))';            % uniform grid covering both neighbours
    fa = stretchCurve(a.depths(:), a.Z(:), z * Ra/Rt);
    fb = stretchCurve(b.depths(:), b.Z(:), z * Rb/Rt);
    Z  = (1-t)*fa + t*fb;
    Z  = Z / max(Z) * ((1-t)*max(a.Z) + t*max(b.Z));
    s.depths = z; s.Z = Z; s.offset = 0;
    [~, iM] = max(Z); s.peakPos = z(iM);
    s.sigma1 = matRad_calcSigmaLatMCS_water(z, E);
    s.sigma2 = s.sigma1 + 1;
    s.weight = zeros(size(z));
    s.initFocus.sigma        = (1-t)*a.initFocus.sigma + t*b.initFocus.sigma;
    s.initFocus.SisFWHMAtIso = (1-t)*a.initFocus.SisFWHMAtIso + t*b.initFocus.SisFWHMAtIso;
    s.energySpectrum.mean  = E;
    s.energySpectrum.sigma = (1-t)*a.energySpectrum.sigma + t*b.energySpectrum.sigma;
    s.interpolated = true;
end


function [sigmaCore, sigmaHalo, haloWeight2D] = fitCoreAndHalo(pos, fluence, coreThr, nSig)
% Fits a measured lateral profile as a narrow core plus a broad halo,
% without an optimiser or any toolbox, in two deterministic stages.
%
% Stage 1 - core. A Gaussian is a parabola in log space, so a quadratic is
% fitted to log(fluence) over the points above coreThr of the peak. The
% 10% default was chosen by testing all 580 profiles in the dataset: it
% succeeds on every one, while a 20% threshold fails on the seven narrowest
% spots (>=209 MeV at Z = +100/+200 mm) where too few samples clear the cut.
%
% Stage 2 - halo. The fitted core is subtracted and the residual beyond
% nSig core sigmas is fitted as a Gaussian centred on the same axis. With
% the centre fixed, log(residual) is LINEAR in u = (x-x0)^2, so this is a
% straight-line fit rather than another curved one.
%
% haloWeight2D is what matRad's 'weight' field expects: the halo's share of
% the 2D fluence, a2*s2^2 / (a1*s1^2 + a2*s2^2). Note the sigma SQUARED -
% the profile is a line cut through a 2D spot, and using the 1D share
% (a*sigma) understates the halo by roughly an order of magnitude. Checked
% against a model-free radial integral of the measured profile,
% 2*pi*r*f(r): the two agree to within 0.02 in weight.
%
% Returns NaN for sigmaCore if the core cannot be fitted; returns a valid
% core with sigmaHalo = NaN and weight 0 if only the halo fit fails, so a
% single-Gaussian fallback stays possible.

    w = max(fluence(:), 0);
    x = pos(:);
    sigmaCore = NaN; sigmaHalo = NaN; haloWeight2D = 0;

    if sum(w) <= 0 || max(w) <= 0
        return;
    end

    % ---- stage 1: core
    core = w >= coreThr * max(w);
    if nnz(core) < 4
        return;
    end
    c = polyfit(x(core), log(w(core)), 2);
    if c(1) >= 0            % not a downward parabola -> no Gaussian core
        return;
    end
    s1 = sqrt(-1/(2*c(1)));
    x0 = -c(2) / (2*c(1));
    a1 = exp(c(3) - c(2)^2/(4*c(1)));

    if ~isfinite(s1) || s1 < 0.2 || s1 > 15
        return;
    end
    sigmaCore = s1;

    % ---- stage 2: halo, from the residual well outside the core
    resid = w - a1 * exp(-(x-x0).^2 ./ (2*s1^2));
    outer = abs(x-x0) > nSig*s1 & resid > 0;
    if nnz(outer) < 5
        return;
    end

    u = (x(outer)-x0).^2;
    q = polyfit(u, log(resid(outer)), 1);
    if q(1) >= 0
        return;
    end
    s2 = sqrt(-1/(2*q(1)));
    a2 = exp(q(2));

    if ~isfinite(s2) || s2 <= s1 || s2 > 100
        return;
    end

    W = a2*s2^2 / (a1*s1^2 + a2*s2^2);
    if ~isfinite(W) || W <= 0 || W >= 0.95
        return;
    end

    sigmaHalo    = s2;
    haloWeight2D = W;
end
function f = stretchCurve(depths, Z, q)
% linear interpolation of a depth curve at query depths q; before the first
% sample the entrance value is held, beyond the last the last value is held
    f = interp1(depths, Z, q, 'linear');
    f(q < depths(1))   = Z(1);
    f(q > depths(end)) = Z(end);
end
function [R80, W80] = braggMetrics(depths, Z)
% distal 80% range and proximal-to-distal 80% width, on a 0.01 mm grid
    zf = (depths(1):0.01:depths(end))';
    Df = interp1(depths(:), Z(:), zf, 'linear');
    [M, iMax] = max(Df);
    iD  = find(Df(iMax:end) < 0.8*M, 1, 'first');
    R80 = zf(iMax + iD - 1);
    iP  = find(Df(1:iMax) < 0.8*M, 1, 'last');
    if isempty(iP)
        W80 = NaN;       % entrance above 80% of peak: width undefined
    else
        W80 = R80 - zf(iP);
    end
end
function blocks = parseRayStationBlocks(filepath)
    fid = fopen(filepath, 'r');
    if fid == -1, error('Could not open file: %s', filepath); end
    rawLines = {};
    tline = fgetl(fid);
    while ischar(tline)
        rawLines{end+1} = tline; %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    blocks = struct('meta', {}, 'data', {});
    curMeta = containers.Map(); curData = {}; inData = false;

    for i = 1:numel(rawLines)
        line = rawLines{i};
        if startsWith(strtrim(line), '#'), continue; end   % RayStation export header lines
        if startsWith(line, 'Origin of data')
            if ~isempty(curMeta) || ~isempty(curData)
                blocks(end+1) = struct('meta', curMeta, 'data', {curData}); %#ok<AGROW>
            end
            curMeta = containers.Map(); curData = {}; inData = false;
        end

        parts = strsplit(line, ';');
        key = strtrim(parts{1});

        if contains(key, 'curve') && (contains(key,'Measured') || contains(key,'curve ('))
            inData = true; continue;
        end

        if inData
            tokens = strtrim(parts);
            tokens = tokens(~cellfun(@isempty, tokens));
            if isempty(tokens), continue; end
            nums = str2double(tokens);
            if any(isnan(nums)), inData = false; continue; end
            curData{end+1} = nums; %#ok<AGROW>
        else
            val = '';
            if numel(parts) > 1, val = strtrim(parts{2}); end
            curMeta(key) = val;
        end
    end
    if ~isempty(curMeta) || ~isempty(curData)
        blocks(end+1) = struct('meta', curMeta, 'data', {curData});
    end
    % keep only real curve blocks (drops anything before the first 'Origin of data')
    keep = arrayfun(@(b) isKey(b.meta, 'Nominal beam energy [MeV/A]:') && ~isempty(b.data), blocks);
    blocks = blocks(keep);
end
function sigmaMCS = matRad_calcSigmaLatMCS_water(depthZ_mm, primaryEnergy_MeV)
% matRad's own analytical Highland/Fermi-Eyges form (Gottschalk 1992),
% water radiation length 36.3 cm, Bragg-Kleeman alpha = 0.0022, p = 1.77.
alpha     = 2.2e-3;
p         = 1.77;
radLength = 36.3;

origSize   = size(depthZ_mm);
depthZ_row = depthZ_mm(:)' ./ 10;
range      = alpha * primaryEnergy_MeV^p;

sigma1  = @(z) 14.1^2 / radLength * (1 + 1/9 * log10(z ./ radLength)).^2;
sigma21 = @(z) 1 ./ (1 - 2/p) .* (range.^(1 - 2/p) .* (range - z).^2 - (range - z).^(3 - 2/p));
sigma22 = @(z) -2 * (range - z) ./ (2 - 2/p) .* (range.^(2 - 2/p) - (range - z).^(2 - 2/p));
sigma23 = @(z) 1 ./ (3 - 2/p) .* (range.^(3 - 2/p) - (range - z).^(3 - 2/p));
sigmaTot = @(z) alpha^(1/p) / 2 * sqrt(sigma1(z) .* (sigma21(z) + sigma22(z) + sigma23(z)));

sigmaBeyond = sigmaTot(range);
isBelowR    = depthZ_row <= range;
isBeyondR   = depthZ_row > range;

sigmaMCS = 10 .* (sigmaTot(depthZ_row) .* isBelowR + sigmaBeyond .* isBeyondR);
sigmaMCS(depthZ_row == 0) = 0;
sigmaMCS = real(sigmaMCS);
sigmaMCS = reshape(sigmaMCS, origSize);
end
