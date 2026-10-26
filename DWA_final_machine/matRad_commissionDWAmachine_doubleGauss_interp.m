%% matRad_commissionDWAmachine_doubleGauss_interp.m
%
% Builds the DWA proton machine for matRad with
%   (1) a DOUBLE-GAUSSIAN lateral model (core + halo), and
%   (2) INTERPOLATED ENERGY LAYERS between the 58 native TOPAS energies.
%
% OUTPUT:  protons_DWA_proton_doubleGauss_interp.mat
% USE IT:  pln.machine = 'DWA_proton_doubleGauss_interp';
%
% VERSION HISTORY (each version has its own script AND its own machine name,
% so old and new results can always be told apart and re-run):
%   matRad_commissionDWAmachine.m                -> DWA_proton
%       single Gaussian (second moment of the air profile), 58 layers
%   matRad_commissionDWAmachine_doubleGauss.m    -> DWA_proton_doubleGauss
%       double Gaussian, 58 layers
%   matRad_commissionDWAmachine_doubleGauss_interp.m (THIS FILE)
%                                                -> DWA_proton_doubleGauss_interp
%       double Gaussian, 201 layers (58 native + 143 interpolated),
%       energy-spread unit fix
%
% CHANGES RELATIVE TO matRad_commissionDWAmachine_doubleGauss.m
%
% (a) Interpolated energy layers - STEP 4b.
%     The native TOPAS layers are 3 MeV apart below 98 MeV and ~4-5 MeV
%     above, i.e. 1-12 mm apart in range. The DWA Bragg peak is very narrow
%     (energy spread 0.04-0.06%), so the layers sit 1.4-2.8 x further apart
%     than a single peak is wide (80%-80% width). No set of spot weights can
%     fill such gaps: in a 1D test the flattest achievable SOBP had a
%     10-29% ripple at every depth, which is the ripple Remo saw.
%     Generic, for comparison, has spacing/width = 0.44.
%     New layers are inserted so that the range step is at most
%     interpSpacingFraction (0.5) x the local 80-80 peak width, with a floor
%     of interpMinStep_mm (1 mm, the sampling of the TOPAS depth curves).
%     Resulting ripple in the same 1D test: < 0.7% at every depth tested.
%     Method: range-scaled interpolation between the two neighbouring native
%     curves (see interpolateLayer below). Validated leave-one-out: each
%     native layer was rebuilt from its neighbours two steps apart (twice the
%     real gap, so a pessimistic test): R80 within 0.2 mm, distal falloff
%     within 0.06 mm (median), 80-80 width +0.14 mm (median).
%
% (b) energySpectrum.sigma is now stored in PERCENT (relative 1-sigma), which
%     is what matRad_MCemittanceBaseData expects (it copies the field
%     straight into EnergySpread, and computes that field itself as
%     spreadInMeV / MeanEnergy * 100). The previous script stored MeV, which
%     overstated the spread by a factor E/100. Only Monte Carlo engines read
%     this field; the analytical pencil-beam dose is unaffected.
%
% (c) No data(i).sigma field (removed 30 Sept 2026, same machine name - the
%     file had not left this computer). With it present, matRad's default
%     lateralModel = 'fast' silently used a single Gaussian (the core only)
%     and dropped the halo. Without it, matRad picks the double Gaussian by
%     default, as it does for Generic.
%
% ---------------------------------------------------------------------------
% BACKGROUND FROM THE DOUBLE-GAUSSIAN VERSION (unchanged)
%
% WHY THIS VERSION EXISTS
% The DWA spot is a narrow core sitting on a broad scattered halo. The
% previous script characterised each profile by its second moment (RMS),
% which was wrong in two compounding ways:
%
%   1. It oversized the spot. At 149 MeV the RMS gives 4.99 mm while the
%      core is 1.27 mm.
%   2. It destroyed the beam optics. Across the 400 mm of measured planes
%      the core varies by 85% (2.19 -> 0.98 mm at 149 MeV, a clear waist
%      near isocenter) while the RMS varies by 6%. matRad was therefore
%      handed an essentially parallel, non-focusing beam.
%
% And collapsing the profile to ONE Gaussian - at any width - discards the
% halo, which is not a small correction. Because a Gaussian's 2D integral
% scales as a*sigma^2 rather than a*sigma, a halo ~17x wider than the core
% carries far more fluence than a line profile suggests: 33% at 50 MeV
% rising to ~50% at 92-219 MeV. Verified model-free by integrating the
% measured profile radially as 2*pi*r*f(r); fitted and direct weights agree
% to within 0.02.
%
% THE HALO IS PHYSICAL, NOT A SCORING ARTIFACT
% From r = 20 to 60 mm the measured fluence falls by a factor 0.033. A flat
% scoring background would give 1.000; a Gaussian of sigma 20-23 mm gives
% 0.018-0.049. So it is a genuine broad scattered component. Its source is
% not yet identified - 124 um of polyimide in the IC64-6 monitor chamber
% cannot produce a component this large, and Alaina confirmed the snout
% bore clears the beam. Open question (does the TOPAS source at the exit
% window already contain it?), but it does not block commissioning: the
% halo is measured, so it can be modelled whatever produces it.
%
% CROSS-CHECK ON THE CORE
% Starting from the 1.0 mm vacuum spot in the beam summary, a Highland
% calculation of the 55 cm air path plus the IC64-6 (124 um polyimide,
% 0.8 um Al, 0.2 um Au, 50.6 cm lever arm) predicts 2.08 / 1.42 / 1.19 /
% 1.10 mm at 50 / 92 / 149 / 219 MeV. The core fit here returns 2.27 /
% 1.50 / 1.27 / 1.17 mm - agreement to 6-9% by a completely independent
% route, which is what justifies treating the core as the source size.
%
% WHAT IS STILL MISSING
%   - The in-water NUCLEAR halo. Generic carries this in its own sigma2
%     with a weight growing 0.001 -> 0.073 from 32 to 150 MeV. Here sigma2
%     carries the IN-AIR halo only, roughly constant with depth, because
%     that is all the data supports. The lateral model is improved, not
%     complete.
%   - IDD truncation. The Bragg curves were scored in a 120 mm diameter
%     detector (r < 60 mm), so 0.3% (50 MeV) to 2.3% (219 MeV) of the
%     fluence falls outside and is missing from Z.
%   - Z normalisation to MeV*cm^2/(g*primary), still open with Remo. The
%     TOPAS curves are scaled to 100 at the peak; absolute per-primary IDDs
%     are needed from Alaina.
%   - LET. No machine.data.LET field, so matRad cannot plan on or report
%     LET with this machine. Needs LETd vs depth in water from TOPAS.
%
% Copyright 2026 - DWA machine commissioning, built from RayStation/TOPAS
% raw exports per Bui et al 2026 Phys. Med. Biol. 71 175030.

%% ----------------------- USER CONFIGURATION ---------------------------

braggFile   = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD_local/DWA_Project/AlainaDWAFiles/TOPAS Output/DWA_PristineBraggPeaks.csv';
spotFile    = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD_local/DWA_Project/AlainaDWAFiles/TOPAS Output/DWA_SpotProfiles.csv';
twissFile   = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD_local/DWA_Project/AlainaDWAFiles/TOPAS Input/AlainasBeamSummaryFinalized.csv';

% Machine name. matRad stores machines as <radiationMode>_<machineName>.mat
% and finds them by that filename when pln.machine is set, so this decides
% both the output file and how plans refer to it:
%
%     machineName = 'DWA_proton_doubleGauss_interp'
%       -> protons_DWA_proton_doubleGauss_interp.mat
%       -> pln.machine = 'DWA_proton_doubleGauss_interp'
%
% Deliberately a NEW name. DWA_proton (single Gaussian) and
% DWA_proton_doubleGauss (58 layers) stay on disk untouched, so every
% version can be run against the others in the same session.
machineName = 'DWA_proton_doubleGauss_interp';

isoToSnout_mm = 400;   % isocenter-to-snout, from the CSV headers
SAD_mm        = 10000; % fixed numerical convention (matches protons_generic_TOPAS.mat),
                       % NOT a physical distance - keeps the beam effectively
                       % parallel for the pencil-beam formalism. Per Remo, Sept 2026.
sourceToIso_mm = 550;  % real source-to-isocenter distance, per Remo. NOTE: this is
                       % NOT currently what goes into meta.BAMStoIsoDist - see there.

energyMatchTol_MeV = 0.05;

coreThreshold   = 0.10;  % core fitted to points above this fraction of the peak
haloCoreSigmas  = 3.0;   % halo fitted beyond this many core sigmas

% --- energy-layer interpolation (STEP 4b)
doInterpolate         = true;  % false reproduces the 58 native layers only
interpSpacingFraction = 0.5;   % max range step as a fraction of the local 80-80
                               % Bragg-peak width. Generic sits at 0.44. In a 1D
                               % SOBP test: 1.0 -> ~1-2% ripple, 0.5 -> <0.7%.
interpMinStep_mm      = 1.0;   % never step finer than the 1 mm sampling of the
                               % TOPAS depth curves (sub-mm layers added nothing
                               % measurable to SOBP flatness in the 1D test)

%% ========================================================================
%% STEP 1: Bragg peak (depth-dose) file
%% ========================================================================
fprintf('Parsing %s ...\n', braggFile);
braggBlocks = parseRayStationBlocks(braggFile);
nE = numel(braggBlocks);
fprintf('  found %d energy blocks\n', nE);

braggEnergy = nan(nE,1); braggDepths = cell(nE,1); braggZ = cell(nE,1);
for k = 1:nE
    b = braggBlocks(k);
    braggEnergy(k) = str2double(b.meta('Nominal beam energy [MeV/A]:'));
    d = vertcat(b.data{:});
    valid = ~any(isnan(d),2);
    braggDepths{k} = d(valid,1);
    braggZ{k}      = d(valid,2);
end
[braggEnergy, sortIdx] = sort(braggEnergy);
braggDepths = braggDepths(sortIdx);
braggZ      = braggZ(sortIdx);

%% ========================================================================
%% STEP 2: Spot profiles - fit a CORE and a HALO to each
%% ========================================================================
fprintf('Parsing %s ...\n', spotFile);
spotBlocks = parseRayStationBlocks(spotFile);
nSpot = numel(spotBlocks);
fprintf('  found %d profile blocks\n', nSpot);

spotEnergy = nan(nSpot,1);
spotAxis   = strings(nSpot,1);
spotZpos   = nan(nSpot,1);
spotCore   = nan(nSpot,1);   % sigma of the narrow component [mm]
spotHalo   = nan(nSpot,1);   % sigma of the broad component  [mm]
spotWeight = nan(nSpot,1);   % 2D fluence fraction in the broad component

for k = 1:nSpot
    b = spotBlocks(k);
    spotEnergy(k) = str2double(b.meta('Nominal beam energy [MeV/A]:'));
    spotAxis(k)   = string(b.meta('Curve type [Depth/X/Y/AbsoluteDosimetry]:'));
    spotZpos(k)   = str2double(b.meta('Air profile Z pos [mm]:'));

    d = vertcat(b.data{:});
    valid = ~any(isnan(d),2);
    d = d(valid,:);
    if spotAxis(k) == "X", pos = d(:,1); else, pos = d(:,2); end
    fluence = d(:,3);

    [spotCore(k), spotHalo(k), spotWeight(k)] = ...
        fitCoreAndHalo(pos, fluence, coreThreshold, haloCoreSigmas);
end

nFailed = sum(~isfinite(spotCore));
fprintf('  core+halo fit: %d of %d profiles failed\n', nFailed, nSpot);
if nFailed > 0
    warning('%d profiles failed the core/halo fit - inspect before trusting the machine.', nFailed);
end

uniqueZ = unique(spotZpos);
fprintf('  physical distances from nozzle (mm, reference only): %s\n', ...
    mat2str((isoToSnout_mm + uniqueZ)'));

%% ========================================================================
%% STEP 3: Twiss / beam summary
%% ========================================================================
fprintf('Parsing %s ...\n', twissFile);
twiss = readtable(twissFile, 'NumHeaderLines', 2, 'ReadVariableNames', false);
twissHeader = readcell(twissFile, 'Range', '1:1');
twiss.Properties.VariableNames = matlab.lang.makeValidName(twissHeader);
fprintf('  found %d rows\n', height(twiss));

%% ========================================================================
%% STEP 4: Assemble machine.data
%% ========================================================================
fprintf('Assembling machine.data (double Gaussian) ...\n');

data = repmat(struct(), nE, 1);
qaTable = table('Size',[nE 7], ...
    'VariableTypes', repmat({'double'},1,7), ...
    'VariableNames', {'energy','sigmaCore_iso_mm','sigmaHalo_iso_mm', ...
                      'haloWeight2D','sigma_vacuum_beamSummary_mm', ...
                      'dEoverE_pct','energySpread_MeV'});

isoIdx = find(uniqueZ == 0, 1);
if isempty(isoIdx)
    error('No profile at Z = 0 (isocenter) - the halo reference plane is undefined.');
end

for i = 1:nE
    E = braggEnergy(i);
    data(i).energy  = E;
    data(i).depths  = braggDepths{i};
    data(i).Z       = braggZ{i};
    data(i).offset  = 0;

    % sub-mm peak position (parabola through the maximum and its two
    % neighbours). On the raw 1 mm grid, 13 interpolated layers would share
    % a peakPos with the next layer, and matRad selects layers by peakPos.
    data(i).peakPos = peakPosition(data(i).depths, data(i).Z);

    % --- in-water MCS broadening, matRad's own Highland/Fermi-Eyges form
    sigmaMCS = matRad_calcSigmaLatMCS_water(data(i).depths, E);

    % --- core / halo vs distance, X and Y averaged
    coreVsZ = nan(numel(uniqueZ),1);
    haloVsZ = nan(numel(uniqueZ),1);
    wVsZ    = nan(numel(uniqueZ),1);
    for j = 1:numel(uniqueZ)
        m = abs(spotEnergy-E) < energyMatchTol_MeV & spotZpos == uniqueZ(j);
        coreVsZ(j) = mean(spotCore(m),  'omitnan');
        haloVsZ(j) = mean(spotHalo(m),  'omitnan');
        wVsZ(j)    = mean(spotWeight(m),'omitnan');
    end

    % --- initFocus carries the CORE, which is the real source size
    data(i).initFocus.dist  = (SAD_mm - uniqueZ(:))';
    data(i).initFocus.sigma = coreVsZ(:)';
    data(i).initFocus.unit  = 'mm';
    [data(i).initFocus.dist, sIdx] = sort(data(i).initFocus.dist);   % MUST be ascending:
    data(i).initFocus.sigma = data(i).initFocus.sigma(sIdx);          % matRad_calcSigmaIni
                                                                      % returns NaN otherwise
    data(i).initFocus.SisFWHMAtIso = coreVsZ(isoIdx) * 2.3548;

    % --- double-Gaussian kernel.
    % matRad forms:  sigmaNarrow^2 = sigma1^2 + sigmaIni^2
    %                sigmaBroad^2  = sigma2^2 + sigmaIni^2
    % with the SAME sigmaIni (from initFocus) on both. With initFocus
    % carrying the core, the narrow component needs only the in-water MCS,
    % and the broad one needs the extra air-side width of the halo on top:
    sigCoreIso = coreVsZ(isoIdx);
    sigHaloIso = haloVsZ(isoIdx);
    extraHaloSq = max(sigHaloIso^2 - sigCoreIso^2, 0);

    data(i).sigma1 = sigmaMCS(:);
    data(i).sigma2 = sqrt(extraHaloSq + sigmaMCS(:).^2);
    data(i).weight = repmat(wVsZ(isoIdx), numel(data(i).depths), 1);
    % weight is depth-independent because the air profiles give no depth
    % information. The real in-water nuclear halo does grow with depth
    % (Generic: 0.001 -> 0.073 over 32-150 MeV) and is NOT included.

    % NO data(i).sigma field. An earlier draft kept sigma = sigmaMCS as a
    % single-Gaussian fallback, but matRad's default lateral model ('fast')
    % uses a single Gaussian whenever that field exists, whatever
    % meta.dataType says - so the halo was silently dropped (found 30 Sept
    % 2026 in the single-beam check: lateral profile = core only).

    % --- energy spectrum from the beam summary.
    % sigma is RELATIVE, IN PERCENT - the convention matRad_MCemittanceBaseData
    % uses (EnergySpread = spreadInMeV / MeanEnergy * 100). The previous
    % script stored MeV here (dEoverE_pct/100*E), overstating it by E/100.
    [~, tRow] = min(abs(twiss.achieved_final_kinetic_energy - E));
    dEoverE_pct = twiss.dE_E_nozzle(tRow);
    data(i).energySpectrum.type  = 'gaussian';
    data(i).energySpectrum.mean  = E;
    data(i).energySpectrum.sigma = dEoverE_pct;

    % provenance flag: native TOPAS layer
    data(i).interpolated = false;

    qaTable.energy(i)                      = E;
    qaTable.sigmaCore_iso_mm(i)            = sigCoreIso;
    qaTable.sigmaHalo_iso_mm(i)            = sigHaloIso;
    qaTable.haloWeight2D(i)                = wVsZ(isoIdx);
    qaTable.sigma_vacuum_beamSummary_mm(i) = twiss.x_iso_1rms_radius_mm(tRow);
    qaTable.dEoverE_pct(i)                 = dEoverE_pct;
    qaTable.energySpread_MeV(i)            = dEoverE_pct/100 * E;
end

%% ========================================================================
%% STEP 4b: Interpolated energy layers
%% ========================================================================
% Native ranges and 80-80 widths
R80native = nan(nE,1); W80native = nan(nE,1);
for i = 1:nE
    [R80native(i), W80native(i)] = braggMetrics(data(i).depths, data(i).Z);
end

% Smooth width model W80 = c * R80^q, fitted where the peak is resolved by
% the 1 mm sampling (R80 > 60 mm). Below that the measured widths are only
% 1-2 bins and jump around, so the model - not the raw width - sets the step.
okW = R80native > 60 & isfinite(W80native);
pW  = polyfit(log(R80native(okW)), log(W80native(okW)), 1);
widthModel = @(R) exp(pW(2)) .* R.^pW(1);
fprintf('  80-80 width model: W80 = %.4f * R80^%.3f  (mm)\n', exp(pW(2)), pW(1));

nativeRatio = diff(R80native) ./ widthModel(0.5*(R80native(1:end-1)+R80native(2:end)));
fprintf('  native layers: range step / peak width = %.2f - %.2f (median %.2f)\n', ...
    min(nativeRatio), max(nativeRatio), median(nativeRatio));

if doInterpolate
    newData = data([]);          % empty struct array with the same fields
    isNative = false(0,1);
    for i = 1:nE
        newData(end+1) = data(i); %#ok<SAGROW>
        isNative(end+1,1) = true; %#ok<SAGROW>
        if i == nE, break; end

        Ra = R80native(i); Rb = R80native(i+1);
        step = max(interpSpacingFraction * widthModel(0.5*(Ra+Rb)), interpMinStep_mm);
        nIns = ceil((Rb - Ra) / step) - 1;
        for j = 1:nIns
            Rt = Ra + (Rb - Ra) * j / (nIns + 1);
            newData(end+1) = interpolateLayer(data(i), data(i+1), Ra, Rb, Rt); %#ok<SAGROW>
            isNative(end+1,1) = false; %#ok<SAGROW>
        end
    end

    % matRad matches energies by equality (intersect/ismember), so they must
    % be unique after rounding and strictly ascending
    allE = [newData.energy];
    if numel(unique(allE)) ~= numel(allE) || any(diff(allE) <= 0)
        error('Interpolated energies are not unique/ascending - reduce interpMinStep_mm resolution or check the rounding.');
    end

    data = newData(:);
    if any(diff([data.peakPos]) <= 0)
        warning('peakPos is not strictly increasing across layers - check the interpolation.');
    end
    fprintf('  interpolation: %d native + %d interpolated = %d layers\n', ...
        nnz(isNative), nnz(~isNative), numel(data));

    R80all = arrayfun(@(d) braggMetrics(d.depths, d.Z), data);
    finalRatio = diff(R80all) ./ widthModel(0.5*(R80all(1:end-1)+R80all(2:end)));
    fprintf('  after interpolation: range step / peak width <= %.2f where R80 > 60 mm\n', ...
        max(finalRatio(R80all(2:end) > 60)));
end

%% ========================================================================
%% STEP 5: machine.meta
%% ========================================================================
meta.radiationMode = 'protons';
meta.machine       = machineName;   % must match the name passed to saveMatradMachine
meta.SAD           = SAD_mm;
meta.BAMStoIsoDist = isoToSnout_mm;   % 400 mm, snout-to-iso, matching the geometry
                                      % the spot profiles were scored in. The 550 mm
                                      % source-to-iso figure is NOT used here - that
                                      % discrepancy is still open with Remo.
meta.MCcode        = 'TOPAS';
meta.dataType      = 'doubleGauss';
meta.fitAirOffset  = 0;
meta.LUTspotSize   = [0 100; 0 0];

meta.interpolation.enabled         = doInterpolate;
meta.interpolation.method          = 'range-scaled blend of neighbouring native depth curves';
meta.interpolation.spacingFraction = interpSpacingFraction;
meta.interpolation.minStep_mm      = interpMinStep_mm;
meta.interpolation.nNative         = nE;
meta.interpolation.nTotal          = numel(data);
meta.interpolation.nativeEnergies  = braggEnergy(:)';

meta.description = [ ...
 'DWA (dielectric wall accelerator) proton machine, built from raw TOPAS/RayStation ' ...
 'commissioning exports (Bui et al 2026 Phys. Med. Biol. 71 175030) by ' ...
 'matRad_commissionDWAmachine_doubleGauss_interp.m. ' ...
 'ENERGY LAYERS: ' num2str(nE) ' native TOPAS energies plus ' num2str(numel(data)-nE) ...
 ' interpolated (data(i).interpolated = true), inserted so the range step is at most ' ...
 num2str(interpSpacingFraction) ' x the local 80-80 Bragg-peak width (floor ' ...
 num2str(interpMinStep_mm) ' mm). The native spacing was 1.4-2.8 x the peak width and ' ...
 'produced a 10-29% SOBP ripple. Depth curves are range-scaled blends of the two ' ...
 'neighbouring native curves; lateral parameters and energy spread are linear blends. ' ...
 'energySpectrum.sigma is RELATIVE IN PERCENT. ' ...
 'LATERAL MODEL: double Gaussian. Each measured air profile is fitted as a narrow ' ...
 'core plus a broad halo. initFocus.sigma carries the CORE (the real source size); ' ...
 'sigma1 is the in-water MCS term; sigma2 additionally carries the air-side halo ' ...
 'width referenced to isocenter; weight is the 2D fluence fraction in the halo ' ...
 '(0.33 at 50 MeV rising to ~0.50 at 92-219 MeV). Using the second moment of the ' ...
 'whole profile instead - as the previous single-Gaussian script did - oversized ' ...
 'the spot ~4x, flattened the beam convergence, and discarded roughly half the ' ...
 'fluence. CORE VALIDATED independently: Highland scattering through the 55 cm air ' ...
 'path plus the IC64-6 monitor chamber, starting from the 1.0 mm vacuum spot, ' ...
 'predicts the fitted core widths to 6-9%. ' ...
 'KNOWN GAPS: (1) the in-water NUCLEAR halo is not modelled - sigma2 carries only ' ...
 'the in-air halo, depth-independent, since the exports contain no in-water lateral ' ...
 'data; Generic by contrast carries a nuclear halo whose weight grows with depth. ' ...
 '(2) the Bragg curves were scored in a 120 mm diameter detector, so 0.3% (50 MeV) ' ...
 'to 2.3% (219 MeV) of the fluence lies outside r = 60 mm and is missing from Z. ' ...
 '(3) Z is NOT normalised to MeV*cm^2/(g*primary) - open with Remo; harmless for ' ...
 'relative dose in matRad, not for absolute dose or a TOPAS recomputation. ' ...
 '(4) the physical origin of the halo is unidentified - too large for the monitor ' ...
 'chamber foils, and the snout bore clears the beam - open with Alaina. ' ...
 '(5) no LET data: machine.data has no LET field, so LET cannot be computed. ' ...
 'SAD = ' num2str(SAD_mm) ' mm is a FIXED NUMERICAL CONVENTION, not a physical ' ...
 'distance. BAMStoIsoDist = ' num2str(isoToSnout_mm) ' mm (snout-to-iso); the 550 mm ' ...
 'source-to-iso alternative is unresolved.'];

machine.meta = meta;
machine.data = data;

%% ========================================================================
%% QA
%% ========================================================================
fprintf('\n--- QA: core / halo fit at isocenter ---\n');
disp(qaTable(1:4:end,:));

figure('Color','w','Position',[100 100 1150 420]);

subplot(1,3,1); hold on; box on; grid on;
plot(qaTable.energy, qaTable.sigmaCore_iso_mm, 'o-','LineWidth',1.5,'DisplayName','core (used as initFocus)');
plot(qaTable.energy, qaTable.sigma_vacuum_beamSummary_mm, 'k--','LineWidth',1.5,'DisplayName','vacuum, beam summary');
xlabel('Energy [MeV]'); ylabel('\sigma at isocenter [mm]');
legend('Location','best'); title('Core vs vacuum spot');
set(gca,'TickDir','out');

subplot(1,3,2); hold on; box on; grid on;
plot(qaTable.energy, qaTable.sigmaHalo_iso_mm, 's-','LineWidth',1.5);
xlabel('Energy [MeV]'); ylabel('halo \sigma at isocenter [mm]');
title('Halo width'); set(gca,'TickDir','out');

subplot(1,3,3); hold on; box on; grid on;
plot(qaTable.energy, 100*qaTable.haloWeight2D, 'd-','LineWidth',1.5);
xlabel('Energy [MeV]'); ylabel('halo share of 2D fluence [%]');
ylim([0 60]); title('Halo weight'); set(gca,'TickDir','out');

% --- QA of the interpolation: native vs interpolated layers
isInt  = [data.interpolated];
allE   = [data.energy];
allR80 = arrayfun(@(d) braggMetrics(d.depths, d.Z), data)';
allW   = arrayfun(@(d) d.weight(1), data)';
allCoreIso = arrayfun(@(d) d.initFocus.SisFWHMAtIso/2.3548, data)';

figure('Color','w','Position',[100 560 1150 420]);
subplot(1,3,1); hold on; box on; grid on;
plot(allE(~isInt), allR80(~isInt), 'ko', 'MarkerFaceColor','k', 'DisplayName','native');
plot(allE(isInt),  allR80(isInt),  'r.', 'MarkerSize',10, 'DisplayName','interpolated');
xlabel('Energy [MeV]'); ylabel('R80 [mm]'); legend('Location','northwest');
title('Range vs energy'); set(gca,'TickDir','out');

subplot(1,3,2); hold on; box on; grid on;
Rmid = 0.5*(allR80(1:end-1)+allR80(2:end));
plot(Rmid, diff(allR80)./widthModel(Rmid), 'r.-', 'DisplayName','this machine');
plot(0.5*(R80native(1:end-1)+R80native(2:end)), nativeRatio, 'ko-', 'DisplayName','native only');
yline(interpSpacingFraction,'k--','HandleVisibility','off');
xlabel('R80 [mm]'); ylabel('range step / 80-80 peak width');
legend('Location','northeast'); title('Layer spacing'); set(gca,'TickDir','out');

subplot(1,3,3); hold on; box on; grid on;
yyaxis left;  plot(allE, allCoreIso, '.-'); ylabel('\sigma_{core} at iso [mm]');
yyaxis right; plot(allE, 100*allW, '.-');   ylabel('halo weight [%]');
xlabel('Energy [MeV]'); title('Lateral parameters, all layers'); set(gca,'TickDir','out');

fprintf('\nmachine.data and machine.meta assembled for %d energies (doubleGauss, %d native).\n', ...
    numel(data), nE);
fprintf('Halo carries %.0f-%.0f%% of the fluence across the energy range.\n', ...
    100*min(qaTable.haloWeight2D), 100*max(qaTable.haloWeight2D));
fprintf('\nThen run:\n');
fprintf('  obj = matRad_MCemittanceBaseData(machine);\n');
fprintf('  obj.saveMatradMachine(''%s'');\n', machineName);
fprintf('\n-> writes protons_%s.mat, leaving protons_DWA_proton.mat\n', machineName);
fprintf('   and protons_DWA_proton_doubleGauss.mat in place for comparison.\n');
fprintf('   Plans must now set pln.machine = ''%s''.\n', machineName);

%% ========================================================================
%% Local functions
%% ========================================================================
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

function s = interpolateLayer(a, b, Ra, Rb, Rt)
% Builds a new energy layer with distal range Rt (R80, mm) between two
% native layers a and b with ranges Ra < Rt < Rb.
%
% Depth dose - range-scaled interpolation. Each neighbour is stretched in
% depth so that its range lands on Rt (z -> z*Ra/Rt and z*Rb/Rt), then the
% two are blended with weight t = (Rt-Ra)/(Rb-Ra). Stretching rather than
% shifting keeps the entrance region and the peak shape consistent, and
% because the two stretched peaks line up, the blend does not smear the
% peak into a double hump. Leave-one-out check on the 58 native layers
% (each rebuilt from neighbours two steps away): R80 within 0.2 mm, distal
% 80-20 falloff within 0.06 mm (median), 80-80 width +0.14 mm (median).
%
% Energy - local power law R = a*E^p through the two neighbours.
% Lateral parameters, energy spread - linear in t. The in-water MCS term is
% recomputed analytically for the new energy, exactly as for native layers.
%
% s starts as a copy of a so that the field ORDER is identical, which MATLAB
% requires when the struct is appended to the machine.data array.

    t = (Rt - Ra) / (Rb - Ra);
    s = a;

    p = log(Rb/Ra) / log(b.energy/a.energy);
    E = a.energy * (Rt/Ra)^(1/p);
    E = round(E*100)/100;                 % same precision as the native list
    s.energy = E;

    z  = a.depths(:);
    fa = stretchCurve(a.depths(:), a.Z(:), z * Ra/Rt);
    fb = stretchCurve(b.depths(:), b.Z(:), z * Rb/Rt);
    Z  = (1-t)*fa + t*fb;
    Z  = Z / max(Z) * ((1-t)*max(a.Z) + t*max(b.Z));

    s.depths = z;
    s.Z      = Z;
    s.offset = 0;
    s.peakPos = peakPosition(z, Z);

    % lateral: MCS recomputed, halo extra width / weight blended
    sigmaMCS = matRad_calcSigmaLatMCS_water(z, E);
    extraA = median(a.sigma2(:).^2 - a.sigma1(:).^2);
    extraB = median(b.sigma2(:).^2 - b.sigma1(:).^2);
    extraSq = max((1-t)*extraA + t*extraB, 0);
    wHalo   = (1-t)*a.weight(1) + t*b.weight(1);

    s.sigma1 = sigmaMCS(:);
    s.sigma2 = sqrt(extraSq + sigmaMCS(:).^2);
    s.weight = repmat(wHalo, numel(z), 1);

    if ~isequal(a.initFocus.dist, b.initFocus.dist)
        error('Neighbouring layers have different initFocus.dist - cannot blend.');
    end
    s.initFocus.sigma        = (1-t)*a.initFocus.sigma + t*b.initFocus.sigma;
    s.initFocus.SisFWHMAtIso = (1-t)*a.initFocus.SisFWHMAtIso + t*b.initFocus.SisFWHMAtIso;

    s.energySpectrum.type  = 'gaussian';
    s.energySpectrum.mean  = E;
    s.energySpectrum.sigma = (1-t)*a.energySpectrum.sigma + t*b.energySpectrum.sigma;  % percent

    s.interpolated = true;
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
    if isempty(iD)       % curve never falls below 80%: range undefined
        R80 = NaN; W80 = NaN; return;
    end
    R80 = zf(iMax + iD - 1);
    iP  = find(Df(1:iMax) < 0.8*M, 1, 'last');
    if isempty(iP)
        W80 = NaN;       % entrance above 80% of peak: width undefined
    else
        W80 = R80 - zf(iP);
    end
end

function pk = peakPosition(depths, Z)
% depth of the maximum, refined by a parabola through the three samples
% around it (vertex formula); falls back to the grid point at the edges
    depths = depths(:); Z = Z(:);
    [~, i] = max(Z);
    pk = depths(i);
    if i > 1 && i < numel(Z)
        den = Z(i-1) - 2*Z(i) + Z(i+1);
        if den < 0
            pk = depths(i) + 0.5*(Z(i-1) - Z(i+1))/den * (depths(i+1) - depths(i));
        end
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
