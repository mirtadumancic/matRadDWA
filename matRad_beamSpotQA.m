%% matRad_beamSpotQA.m
%
% Macro to load a matRad proton machine base data file and produce a set
% of publication-quality plots characterizing the beam spot size as a
% function of energy, focus setting and distance from the nozzle.
%
% Handles TWO machine.data(i).initFocus formats seen in matRad basedata:
%
%   FORMAT A - fitted emittance (e.g. protons_generic_TOPAS.mat):
%       .emittance(f).sigmaX/sigmaY/divX/divY/corrX/corrY   (per focus f)
%       .SisFWHMAtIso   [1 x nFoci]
%       .dist, .sigma   [nFoci x nSamples]
%     -> spot size/divergence at isocenter come directly from the fitted
%        Courant-Snyder optics (accurate).
%
%   FORMAT B - raw envelope only (e.g. protons_generic_MCsquare.mat):
%       .dist, .sigma   [1 x nSamples] (single focus, no emittance fit)
%       .SisFWHMAtIso   scalar
%     -> spot size/divergence at isocenter are ESTIMATED by interpolating
%        the raw sigma(distance) sampling at meta.SAD (approximate -
%        divergence here is a finite-difference dSigma/dDistance, not a
%        fitted emittance parameter, and is labeled as such in the plots).
%
% Plots produced:
%   1) Spot sigma at isocenter vs. nominal energy, one line per focus
%   2) Divergence at isocenter vs. nominal energy, one line per focus
%      (fitted, if available - otherwise apparent dSigma/dz)
%   3) Beam envelope (sigma vs. distance from nozzle) for a low- and a
%      high-energy layer, one line per focus setting
%   4) FWHM at isocenter vs. nominal energy, one line per focus
%
% All figures are saved as high-resolution PNGs (and .fig files) into
% outputDir, with filenames prefixed by the machine name so plots from
% multiple machine files can coexist in the same folder.
%
% Copyright 2026 - generated for spot-size commissioning QA plots.

%% ----------------------- USER CONFIGURATION ---------------------------
matFilePath = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD/matRad-master/matRad/basedata/protons_generic_MCsquare.mat';
outputDir   = '/Users/mirtadumancic/Work/MATLAB_Projects/matRAD/DWA_Project';

% Indices (into the sorted energy list) of the two energy layers used
% for the beam-envelope plot. 1 = lowest energy, end = highest energy.
% Leave as [] to auto-pick lowest/highest.
envelopeEnergyIdx = [];

%% ------------------------------------------------------------------
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

fprintf('Loading machine file: %s\n', matFilePath);
loaded = load(matFilePath, 'machine');
machine = loaded.machine;

data = machine.data;
nEnergies = numel(data);

machineName = 'machine';
if isfield(machine, 'meta') && isfield(machine.meta, 'machine')
    machineName = machine.meta.machine;
end
fileTag = regexprep(machineName, '[^a-zA-Z0-9_-]', '_');

SAD = NaN;
if isfield(machine, 'meta') && isfield(machine.meta, 'SAD')
    SAD = machine.meta.SAD;
end

fprintf('Loaded "%s": %d energy layers (%.2f - %.2f MeV)\n', ...
    machineName, nEnergies, min([data.energy]), max([data.energy]));

% Sort by energy just in case the file isn't already ordered
[energies, sortIdx] = sort([data.energy]);
data = data(sortIdx);

hasEmittance = isfield(data(1).initFocus, 'emittance');

if hasEmittance
    nFoci = numel(data(1).initFocus.emittance);
    fprintf('Detected FITTED emittance data, %d focus setting(s).\n', nFoci);
else
    d1 = data(1).initFocus.dist;
    nFoci = size(d1, 1);   % 1 if dist is a plain row vector
    fprintf(['Detected RAW envelope data only (no fitted emittance), ' ...
        '%d focus setting(s).\nSpot size/divergence at isocenter will ' ...
        'be estimated by interpolating sigma(distance) at SAD = %g mm.\n'], ...
        nFoci, SAD);
    if isnan(SAD)
        error(['machine.meta.SAD not found - cannot locate isocenter ' ...
            'along the sampled distance axis.']);
    end
end

focusColors = lines(nFoci);
focusLabels = arrayfun(@(f) sprintf('Focus %d', f), 1:nFoci, 'UniformOutput', false);

%% ----------------- Collect spot-size / divergence data ----------------
sigmaAtIso = nan(nEnergies, nFoci);
divAtIso   = nan(nEnergies, nFoci);   % rad (fitted) or apparent dSigma/dz
fwhmIso    = nan(nEnergies, nFoci);

for iE = 1:nEnergies
    focus = data(iE).initFocus;

    if hasEmittance
        for f = 1:nFoci
            em = focus.emittance(f);
            sigmaAtIso(iE,f) = em.sigmaX(1);
            divAtIso(iE,f)   = em.divX(1);
        end
        fwhmVals = focus.SisFWHMAtIso;
        fwhmIso(iE, 1:numel(fwhmVals)) = fwhmVals;
    else
        distMat  = focus.dist;
        sigmaMat = focus.sigma;
        if isvector(distMat)
            distMat  = distMat(:)';
            sigmaMat = sigmaMat(:)';
        end
        for f = 1:nFoci
            d = distMat(f,:);
            s = sigmaMat(f,:);
            sigmaAtIso(iE,f) = interp1(d, s, SAD, 'linear', 'extrap');
            divAtIso(iE,f)   = interp1(d, gradient(s, d), SAD, 'linear', 'extrap');
        end
        fwhmVals = focus.SisFWHMAtIso;   % scalar in format B
        fwhmIso(iE, 1:numel(fwhmVals)) = fwhmVals;
    end
end

%% ============ PLOT 1: Spot sigma at isocenter vs. energy ==============
fig1 = figure('Color', 'w', 'Position', [100 100 800 550]);
hold on; box on; grid on;
for f = 1:nFoci
    plot(energies, sigmaAtIso(:,f), '-', 'Color', focusColors(f,:), ...
        'LineWidth', 2, 'DisplayName', focusLabels{f});
end
xlabel('Nominal energy [MeV]', 'FontSize', 12);
ylabel('Spot size \sigma at isocenter [mm]', 'FontSize', 12);
title(sprintf('%s: Spot size at isocenter vs. energy', machineName), ...
    'Interpreter', 'none', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);
saveFigure(fig1, outputDir, fileTag, 'spotSize_vs_energy');

%% ============ PLOT 2: Divergence at isocenter vs. energy ==============
fig2 = figure('Color', 'w', 'Position', [100 100 800 550]);
hold on; box on; grid on;
for f = 1:nFoci
    plot(energies, divAtIso(:,f) * 1e3, '-', 'Color', focusColors(f,:), ...
        'LineWidth', 2, 'DisplayName', focusLabels{f});
end
xlabel('Nominal energy [MeV]', 'FontSize', 12);
if hasEmittance
    ylabel('Fitted beam divergence at isocenter [mrad]', 'FontSize', 12);
    titleStr = sprintf('%s: Fitted beam divergence at isocenter vs. energy', machineName);
else
    ylabel('Apparent divergence d\sigma/dz at isocenter [mrad]', 'FontSize', 12);
    titleStr = sprintf('%s: Apparent beam divergence (finite-difference) vs. energy', machineName);
end
title(titleStr, 'Interpreter', 'none', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);
saveFigure(fig2, outputDir, fileTag, 'divergence_vs_energy');

%% ============ PLOT 3: Beam envelope (sigma vs. distance) ==============
if isempty(envelopeEnergyIdx)
    envelopeEnergyIdx = [1, nEnergies];
end

fig3 = figure('Color', 'w', 'Position', [100 100 1100 500]);
for k = 1:numel(envelopeEnergyIdx)
    iE = envelopeEnergyIdx(k);
    focus = data(iE).initFocus;
    distMat  = focus.dist;
    sigmaMat = focus.sigma;
    if isvector(distMat)
        distMat  = distMat(:)';
        sigmaMat = sigmaMat(:)';
    end

    subplot(1, numel(envelopeEnergyIdx), k); hold on; box on; grid on;
    for f = 1:nFoci
        plot(distMat(f,:), sigmaMat(f,:), '-o', 'Color', focusColors(f,:), ...
            'LineWidth', 1.8, 'MarkerFaceColor', focusColors(f,:), ...
            'MarkerSize', 4, 'DisplayName', focusLabels{f});
    end
    if ~isnan(SAD)
        xline(SAD, '--k', 'SAD', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
    xlabel('Distance from nozzle [mm]', 'FontSize', 11);
    ylabel('Spot size \sigma [mm]', 'FontSize', 11);
    title(sprintf('E = %.1f MeV', energies(iE)), 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 9);
    set(gca, 'FontSize', 10);
end
sgtitle(sprintf('%s: Beam envelope for low/high energy layers', machineName), ...
    'Interpreter', 'none', 'FontSize', 13);
saveFigure(fig3, outputDir, fileTag, 'beam_envelope_lowhigh_energy');

%% ============ PLOT 4: FWHM at isocenter vs. energy ====================
fig4 = figure('Color', 'w', 'Position', [100 100 800 550]);
hold on; box on; grid on;
for f = 1:nFoci
    plot(energies, fwhmIso(:,f), '-', 'Color', focusColors(f,:), ...
        'LineWidth', 2, 'DisplayName', focusLabels{f});
end
xlabel('Nominal energy [MeV]', 'FontSize', 12);
ylabel('FWHM at isocenter [mm]', 'FontSize', 12);
title(sprintf('%s: Spot FWHM at isocenter vs. energy', machineName), ...
    'Interpreter', 'none', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);
saveFigure(fig4, outputDir, fileTag, 'fwhm_vs_energy');

fprintf('\nAll plots saved to: %s\n', outputDir);

%% ------------------------- helper function -----------------------------
function saveFigure(fig, outDir, fileTag, name)
    baseName = sprintf('%s_%s', fileTag, name);
    pngPath = fullfile(outDir, [baseName '.png']);
    figPath = fullfile(outDir, [baseName '.fig']);
    exportgraphics(fig, pngPath, 'Resolution', 300);
    savefig(fig, figPath);
    fprintf('  saved %s\n', pngPath);
end
