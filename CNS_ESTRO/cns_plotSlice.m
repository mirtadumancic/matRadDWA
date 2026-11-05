function cns_plotSlice(ax, ct, cst, cube, z, window, cmap, names, labelText, showLegend, thresh)
% cns_plotSlice  Axial CT slice with a dose/LET colour wash and contours (no transparency needed).
%
%   cns_plotSlice(ax, ct, cst, cube, z, [0 77], jet(64), {'CTV','Brainstem'}, 'Gy(RBE)')
%
% The wash is blended into the grey CT in RGB, so it looks the same in
% MATLAB and Octave and prints reliably. Values below thresh (fraction of
% the window, default 0.05) are not coloured. Cube and CT share the CT grid ([y x z]). The view is
% cropped to the Body contour on that slice.
if nargin < 10, showLegend = true; end
if nargin < 11, thresh = 0.05; end
if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU)
    C = double(ct.cubeHU{1}(:, :, z)); cw = [-200 300];
else
    C = double(ct.cube{1}(:, :, z)); cw = [0.8 1.3];
end
g = min(max((C - cw(1)) / diff(cw), 0), 1);
rgb = repmat(g, [1 1 3]);

D = double(cube(:, :, z));
idx = round((min(max(D, window(1)), window(2)) - window(1)) / diff(window) * (size(cmap, 1) - 1)) + 1;
col = reshape(cmap(idx(:), :), [size(D) 3]);
show = repmat(D > window(1) + thresh * diff(window), [1 1 3]);
a = 0.6;
rgb(show) = (1 - a) * rgb(show) + a * col(show);

axes(ax); %#ok<LAXES>
image(rgb); axis image off; hold on;
colormap(ax, cmap); caxis(ax, window);
cb = colorbar('peer', ax); ylabel(cb, labelText);

cols = lines(numel(names));
h = []; lab = {};
for i = 1:numel(names)
    ix = find(strcmp(cst(:,2), names{i}), 1);
    if isempty(ix), continue; end
    m = false(ct.cubeDim); m(cst{ix,4}{1}) = true;
    s = double(m(:, :, z));
    if ~any(s(:)), continue; end
    [~, hh] = contour(ax, s, [0.5 0.5], 'LineColor', cols(i, :), 'LineWidth', 1.2);
    h(end+1) = hh; lab{end+1} = names{i}; %#ok<AGROW>
end
if showLegend && ~isempty(h), legend(h, lab, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'none'); end
ib = find(strcmp(cst(:,2), 'Body'), 1);
if ~isempty(ib)
    m = false(ct.cubeDim); m(cst{ib,4}{1}) = true; m = m(:, :, z);
    r = find(any(m, 2)); c = find(any(m, 1));
    if ~isempty(r), ylim(ax, [max(1, r(1) - 5), min(size(m, 1), r(end) + 5)]); end
    if ~isempty(c), xlim(ax, [max(1, c(1) - 5), min(size(m, 2), c(end) + 5)]); end
end
hold off;
end
