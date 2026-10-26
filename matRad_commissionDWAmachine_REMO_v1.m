%% matRad_commissionDWAmachine.m
%
% Parses the three raw commissioning exports for the DWA proton machine
% and assembles a matRad-format "pre-fit" machine struct:
%
%   DWA_PristineBraggPeaks.csv        -> depths / Z (IDD curves), 58 energies
%   DWA_SpotProfiles.csv              -> in-air lateral fluence profiles,
%                                         58 energies x 5 distances x X/Y
%   Alaina_sBeamSummary_-_Finalized.csv -> Twiss (alpha/beta/gamma/emittance)
%                                         and dE/E at the nozzle per energy
%
% This script does NOT run matRad_MCemittanceBaseData or save a final
% machine file - it builds and QA-checks machine.data/machine.meta so you
% can inspect it in MATLAB first, then run the emittance fit yourself:
%
%   obj = matRad_MCemittanceBaseData(machine);
%   obj.saveMatradMachine('DWA_proton');
%
% Copyright 2026 - DWA machine commissioning, built from RayStation/TOPAS
% raw exports per Bui et al 2026 Phys. Med. Biol. 71 175030.

%% ----------------------- USER CONFIGURATION ---------------------------

braggFile   = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD_local/DWA_Project/AlainaDWAFiles/TOPAS Output/DWA_PristineBraggPeaks.csv';
spotFile    = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD_local/DWA_Project/AlainaDWAFiles/TOPAS Output/DWA_SpotProfiles.csv';
twissFile   = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD_local/DWA_Project/AlainaDWAFiles/TOPAS Input/AlainasBeamSummaryFinalized.csv';

isoToSnout_mm = 400;   % physical isocenter-to-snout distance, from CSV headers
                       % (used below only as a cross-check reference value -
                       % NOT what goes into meta.BAMStoIsoDist, see note there)

SAD_mm = 10000;        % Fixed numerical convention (matches protons_generic_TOPAS.mat),
                       % NOT a physical distance - keeps the beam effectively
                       % parallel for the pencil-beam formalism.
                       % Confirmed with Remo Cristoforetti (DKFZ postdoc), Sept 2026.

sourceToIso_mm = 550;   % REAL source-to-isocenter distance (55cm) - per Remo
                        % Cristoforetti, this goes into meta.BAMStoIsoDist
                        % (not the physical snout distance above).

energyMatchTol_MeV = 0.05; % tolerance for matching energies across the 3 files

%% ========================================================================
%% STEP 1: Parse the Bragg peak (depth-dose) file
%% ========================================================================
fprintf('Parsing %s ...\n', braggFile);
braggBlocks = parseRayStationBlocks(braggFile);

nE = numel(braggBlocks);
fprintf('  found %d energy blocks\n', nE);

braggEnergy = nan(nE,1);
braggDepths = cell(nE,1);
braggZ      = cell(nE,1);

for k = 1:nE
    b = braggBlocks(k);
    braggEnergy(k) = str2double(b.meta('Nominal beam energy [MeV/A]:'));
    d = vertcat(b.data{:});         % [depth, dose] - vertcat forces row-stacking,
                                     % unlike cell2mat which horzcats a 1xN cell
    valid = ~any(isnan(d),2);       % drop the trailing 'End:' marker row
    braggDepths{k} = d(valid,1);
    braggZ{k}      = d(valid,2);
end

[braggEnergy, sortIdx] = sort(braggEnergy);
braggDepths = braggDepths(sortIdx);
braggZ      = braggZ(sortIdx);

%% ========================================================================
%% STEP 2: Parse the spot profile file and fit sigma at each (energy, axis, z)
%% ========================================================================
fprintf('Parsing %s ...\n', spotFile);
spotBlocks = parseRayStationBlocks(spotFile);
fprintf('  found %d profile blocks\n', numel(spotBlocks));

nSpot = numel(spotBlocks);
spotEnergy = nan(nSpot,1);
spotAxis   = strings(nSpot,1);
spotZpos   = nan(nSpot,1);
spotSigma  = nan(nSpot,1);

for k = 1:nSpot
    b = spotBlocks(k);
    spotEnergy(k) = str2double(b.meta('Nominal beam energy [MeV/A]:'));
    spotAxis(k)   = string(b.meta('Curve type [Depth/X/Y/AbsoluteDosimetry]:'));
    spotZpos(k)   = str2double(b.meta('Air profile Z pos [mm]:'));

    d = vertcat(b.data{:}); % [xpos, ypos, fluence] - vertcat, not cell2mat (see note above)
    valid = ~any(isnan(d),2);
    d = d(valid,:);
    if spotAxis(k) == "X"
        pos = d(:,1);
    else
        pos = d(:,2);
    end
    fluence = d(:,3);

    spotSigma(k) = weightedSigma(pos, fluence);
end

% distances from nozzle (confirmed convention: dist = isoToSnout + Zpos,
% cross-checked against Bui et al 2026, Fig 15 caption). NOTE: this is the
% physical distance for interpretation/cross-checking only - it is NOT
% what gets stored in initFocus.dist (see Step 4-6 below, which uses the
% SAD-referenced convention matRad_MCemittanceBaseData actually expects).
uniqueZ = unique(spotZpos);
distFromNozzle_mm = isoToSnout_mm + uniqueZ;   % e.g. 200,300,400,500,600
fprintf('  physical distances from nozzle (mm, for reference only): %s\n', mat2str(distFromNozzle_mm'));

%% ========================================================================
%% STEP 3: Parse the Twiss/beam-summary file
%% ========================================================================
fprintf('Parsing %s ...\n', twissFile);
twiss = readtable(twissFile, 'NumHeaderLines', 2, ...
    'ReadVariableNames', false);
twissHeader = readcell(twissFile, 'Range', '1:1');
twiss.Properties.VariableNames = matlab.lang.makeValidName(twissHeader);

fprintf('  found %d rows\n', height(twiss));

%% ========================================================================
%% STEP 4-6: Assemble machine.data, matching all three sources by energy
%% ========================================================================
fprintf('Assembling machine.data ...\n');

data = repmat(struct(), nE, 1);
qaTable = table('Size',[nE 6], ...
    'VariableTypes', {'double','double','double','double','double','double'}, ...
    'VariableNames', {'energy','sigmaX_iso_fit_mm','sigmaY_iso_fit_mm', ...
                       'sigma_iso_beamSummary_mm','dEoverE_pct','energySpread_MeV'});

for i = 1:nE
    E = braggEnergy(i);
    data(i).energy = E;
    data(i).depths = braggDepths{i};
    data(i).Z      = braggZ{i};
    data(i).sigma = matRad_calcSigmaLatMCS_water(data(i).depths, data(i).energy);
    [~, maxI] = max(data(i).Z);
    data(i).peakPos = data(i).depths(maxI);
    data(i).offset = 0;   % Isocenter to phantom surface distance = 0 mm

    % --- lateral envelope from measured spot profiles ---
    sigX = nan(size(distFromNozzle_mm));
    sigY = nan(size(distFromNozzle_mm));
    for j = 1:numel(uniqueZ)
        maskX = abs(spotEnergy-E)<energyMatchTol_MeV & spotAxis=="X" & spotZpos==uniqueZ(j);
        maskY = abs(spotEnergy-E)<energyMatchTol_MeV & spotAxis=="Y" & spotZpos==uniqueZ(j);
        if any(maskX), sigX(j) = spotSigma(find(maskX,1)); end
        if any(maskY), sigY(j) = spotSigma(find(maskY,1)); end
    end
    sigmaAvg = mean([sigX(:) sigY(:)], 2, 'omitnan');

    data(i).initFocus.dist  = (SAD_mm - uniqueZ(:))';   % 1 x 5, mm
                                     % SAD-referenced convention: z = SAD-dist
                                     % must recover the physical offset from
                                     % isocenter (uniqueZ), NOT distance from
                                     % the nozzle - see note below.
    data(i).initFocus.sigma = sigmaAvg(:)';             % 1 x 5, mm
    data(i).initFocus.unit  = 'mm';

    [data(i).initFocus.dist, sortIdx] = sort(data(i).initFocus.dist);
    data(i).initFocus.sigma = data(i).initFocus.sigma(sortIdx);

    isoIdx = find(uniqueZ==0, 1);
    sigmaAtIso = sigmaAvg(isoIdx);
    data(i).initFocus.SisFWHMAtIso = sigmaAtIso * 2.3548;

    % --- energy spectrum from the beam-summary dE/E (not back-fit from Z) ---
    [~, tRow] = min(abs(twiss.achieved_final_kinetic_energy - E));
    dEoverE_pct = twiss.dE_E_nozzle(tRow);
    data(i).energySpectrum.type  = 'gaussian';
    data(i).energySpectrum.mean  = E;
    data(i).energySpectrum.sigma = dEoverE_pct/100 * E;

    % --- QA row ---
    qaTable.energy(i)                   = E;
    qaTable.sigmaX_iso_fit_mm(i)        = sigX(isoIdx);
    qaTable.sigmaY_iso_fit_mm(i)        = sigY(isoIdx);
    qaTable.sigma_iso_beamSummary_mm(i) = twiss.x_iso_1rms_radius_mm(tRow);
    qaTable.dEoverE_pct(i)              = dEoverE_pct;
    qaTable.energySpread_MeV(i)         = data(i).energySpectrum.sigma;
end

%% ========================================================================
%% STEP 7: machine.meta
%% ========================================================================
meta.radiationMode = 'protons';
meta.machine       = 'DWA_proton';   % must match the filename used in saveMatradMachine
meta.SAD           = SAD_mm;             % fixed numerical convention (10000mm) - see note above
%meta.BAMStoIsoDist  = sourceToIso_mm;    % REAL source-to-isocenter distance (550mm) - per Remo Cristoforetti
meta.BAMStoIsoDist  = isoToSnout_mm;      % 400mm, snout-to-iso (matches real spot-profile measurement geometry)
meta.MCcode         = 'TOPAS';
meta.dataType        = 'singleGauss';
meta.fitAirOffset    = 0;
meta.LUTspotSize = [0 100; 0 0];   % flat LUT -> focusIx always 1 (we only fit one focus per energy)
meta.created_by      = 'matRad_buildDWAmachine.m (auto-generated)';
meta.created_on      = datestr(now, 'dd-mmm-yyyy');
meta.description = [ ...
    'DWA (dielectric wall accelerator) proton machine, built from raw ' ...
    'TOPAS/RayStation commissioning exports (Bui et al 2026 Phys. Med. Biol. ' ...
    '71 175030). Depth-dose from DWA_PristineBraggPeaks.csv; lateral spot ' ...
    'size fit from DWA_SpotProfiles.csv (weighted-variance sigma, X/Y ' ...
    'averaged); energy spread taken directly from dE/E_nozzle in ' ...
    'Alaina''sBeamSummary_-_Finalized.csv rather than back-fit from the ' ...
    'Bragg curve. Fitted sigma at isocenter (~4.5-5.5mm) matches the ' ...
    'paper''s 5mm virtual DWA machine, NOT the beam summary''s ' ...
    'x_iso_1rms_radius_mm column (~1.0 for all energies - likely an ' ...
    'achieved/requested ratio, not a literal mm value; unresolved, ' ...
    'flag for Alaina). SAD=', num2str(SAD_mm), 'mm is a FIXED NUMERICAL ' ...
    'CONVENTION (matches protons_generic_TOPAS.mat), not a physical ' ...
    'distance - keeps the beam parallel for the pencil-beam formalism; ' ...
    'BAMStoIsoDist=', num2str(sourceToIso_mm), 'mm is the REAL source-to-' ...
    'isocenter distance. Both confirmed with Remo Cristoforetti (DKFZ), ' ...
    'Sept 2026 - SAD/BAMStoIsoDist derivation from first principles is ' ...
    'no longer needed. In-water lateral broadening uses matRad''s built-in ' ...
    'analytical MCS (multiple Coulomb scattering) formula on top of the ' ...
    'fitted initFocus spot size, not a pre-measured per-depth table - ' ...
    'confirmed this is standard architecture, not a gap in this machine. ' ...
    'Intended for MC dose calculation only (no in-water depth-resolved ' ...
    'lateral data available for cross-checking the analytical MCS model)' ...
    'data(i).sigma (in-water lateral MCS broadening vs depth) is currently computed'...
    'from matRads own Highland/Fermi-Eyges formula (Gottschalk 1992, water' ...
    'radLength=36.3cm), NOT measured. None of the TOPAS/RayStation commissioning' ...
    'exports (Bragg peaks, spot profiles, beam summary) contained an in-water' ...
    'lateral-spread-vs-depth dataset. UNRESOLVED, flag for Remo - replace with' ...
    'measured data if/when available.'];

machine.meta = meta;
machine.data = data;

%% ========================================================================
%% QA output
%% ========================================================================
fprintf('\n--- QA: measured spot-profile fit vs. beam-summary achieved spot size at isocenter ---\n');
disp(qaTable);

figure('Color','w');
subplot(1,2,1); hold on; box on; grid on;
plot(qaTable.energy, qaTable.sigmaX_iso_fit_mm, 'o-', 'DisplayName','X (fit)');
plot(qaTable.energy, qaTable.sigmaY_iso_fit_mm, 's-', 'DisplayName','Y (fit)');
plot(qaTable.energy, qaTable.sigma_iso_beamSummary_mm, 'k--', 'DisplayName','beam summary');
xlabel('Energy [MeV]'); ylabel('\sigma at isocenter [mm]');
legend('Location','best'); title('Spot size cross-check');

subplot(1,2,2); hold on; box on; grid on;
plot(qaTable.energy, qaTable.energySpread_MeV, 'o-');
xlabel('Energy [MeV]'); ylabel('Energy spread \sigma [MeV]');
title('Energy spread from dE/E_{nozzle}');

fprintf('\nmachine.data and machine.meta assembled for %d energies.\n', nE);
fprintf('Inspect machine and the QA plot/table above, then run:\n');
fprintf('  obj = matRad_MCemittanceBaseData(machine);\n');
fprintf('  obj.saveMatradMachine(''DWA_proton'');\n');

%% ========================================================================
%% Local functions
%% ========================================================================
function blocks = parseRayStationBlocks(filepath)
    % Parses a RayStation-style commissioning CSV into an array of blocks,
    % each with a containers.Map of metadata and a numeric data matrix.
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Could not open file: %s', filepath);
    end
    rawLines = {};
    tline = fgetl(fid);
    while ischar(tline)
        rawLines{end+1} = tline; %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    blocks = struct('meta', {}, 'data', {});
    curMeta = containers.Map();
    curData = {};
    inData = false;

    for i = 1:numel(rawLines)
        line = rawLines{i};
        if startsWith(line, 'Origin of data')
            if ~isempty(curMeta) || ~isempty(curData)
                blocks(end+1) = struct('meta', curMeta, 'data', {curData}); %#ok<AGROW>
            end
            curMeta = containers.Map();
            curData = {};
            inData = false;
        end

        parts = strsplit(line, ';');
        key = strtrim(parts{1});

        if contains(key, 'curve') && (contains(key,'Measured') || contains(key,'curve ('))
            inData = true;
            continue;
        end

        if inData
            tokens = strtrim(parts);
            tokens = tokens(~cellfun(@isempty, tokens));
            if isempty(tokens)
                continue; % blank line within/after the data block - skip
            end
            nums = str2double(tokens);
            if any(isnan(nums))
                % first non-numeric line after the data rows (e.g. the
                % trailing 'End:' marker) - data block is finished
                inData = false;
                continue;
            end
            curData{end+1} = nums; %#ok<AGROW>
        else
            val = '';
            if numel(parts) > 1
                val = strtrim(parts{2});
            end
            curMeta(key) = val;
        end
    end
    if ~isempty(curMeta) || ~isempty(curData)
        blocks(end+1) = struct('meta', curMeta, 'data', {curData});
    end
end

function sigma = weightedSigma(pos, fluence)
    % Weighted-variance estimate of the lateral spread, treating the
    % measured fluence as a distribution over position. Robust and
    % optimizer-free; simple to swap for a Gaussian curve_fit if the
    % measured tails prove too noisy.
    w = max(fluence, 0);
    if sum(w) <= 0
        sigma = NaN;
        return;
    end
    x0 = sum(w.*pos) / sum(w);
    sigma = sqrt( sum(w.*(pos-x0).^2) / sum(w) );
end

function sigmaMCS = matRad_calcSigmaLatMCS_water(depthZ_mm, primaryEnergy_MeV)
alpha     = 2.2e-3;
p         = 1.77;
radLength = 36.3;

origSize = size(depthZ_mm);
depthZ_row = depthZ_mm(:)' ./ 10;            % work internally as a row
range = alpha * primaryEnergy_MeV^p;

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

sigmaMCS = reshape(sigmaMCS, origSize);      % match caller's shape (row or column)
end

