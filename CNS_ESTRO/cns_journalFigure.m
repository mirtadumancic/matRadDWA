function f = cns_journalFigure(widthCm, heightCm)
% cns_journalFigure  New figure sized for a Medical Physics / PMB / Green Journal column.
%
%   f = cns_journalFigure(8.5, 6.5)    % single column (8.5 cm)
%   f = cns_journalFigure(17.5, 10)    % double column (17.5 cm)
%
% Arial 8 pt, white background, axes lines 0.75 pt, no figure toolbar.
% Save with cns_saveJournalFigure (vector PDF + 300 dpi PNG + csv of the data).
if nargin < 1, widthCm = 8.5; end
if nargin < 2, heightCm = 0.75 * widthCm; end
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 widthCm heightCm], ...
    'PaperUnits', 'centimeters', 'PaperSize', [widthCm heightCm], ...
    'PaperPosition', [0 0 widthCm heightCm], 'Visible', get(0, 'DefaultFigureVisible'));
props = {'DefaultAxesFontName', 'Arial'; 'DefaultAxesFontSize', 8; 'DefaultTextFontName', 'Arial'; ...
    'DefaultTextFontSize', 8; 'DefaultAxesLineWidth', 0.75; 'DefaultLineLineWidth', 1.0; ...
    'DefaultAxesTickDir', 'out'; 'DefaultAxesBox', 'off'; 'DefaultLegendFontSize', 7; 'DefaultLegendBox', 'off'};
% MATLAB R2025a+ follows the desktop theme (dark theme = dark axes, light text):
% force the light theme so figures look the same everywhere and print correctly.
try, set(f, 'Theme', 'light'); catch, end          %#ok<NOCOM> % property exists from R2025a
props = [props; {'DefaultAxesColor', 'w'; 'DefaultAxesXColor', [0.1 0.1 0.1]; ...
    'DefaultAxesYColor', [0.1 0.1 0.1]; 'DefaultAxesGridColor', [0.15 0.15 0.15]; ...
    'DefaultTextColor', 'k'; 'DefaultLegendTextColor', 'k'; 'DefaultLegendColor', 'w'}];
for i = 1:size(props, 1)
    try, set(f, props{i, 1}, props{i, 2}); catch, end   %#ok<NOCOM> % not all exist in Octave
end
end
