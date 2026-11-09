function cns_figSliceRS_v2(ct, cst, panels, opt, base)
% cns_figSliceRS_v2  (Nov 9 2026) same as cns_figSliceRS, with the layout set in centimetres so that
% nothing is cut off: a title strip above every panel, and the colour-bar label written along the
% bar (the label above the bar ran over the right edge of the figure in one-row figures).
%
% cns_figSliceRS  Publication-quality dose / LET slice panels in a treatment-planning-system style.
%
%   cns_figSliceRS(ct, cst, panels, opt, 'D:/.../fig_slices')
%
% panels  struct array (row-major), fields:
%   cube   3D array on the CT grid: total dose (Gy / Gy(RBE)) or LETd (keV/um) or a difference
%   type   'dose' | 'let' | 'diff'
%   title  panel title
%   row, col  position in the grid
%   doseForMask  (type 'let' / 'diff') total dose cube used to hide low-dose voxels
% opt (all optional)
%   plane        'sagittal' (default) | 'axial' | 'coronal'
%   slice        index along the plane normal; default = centroid of opt.centreOn
%   centreOn     structure for the default slice (default 'CTV')
%   rxGy         prescription (total), sets the % scale (default 54)
%   washMinPct   lowest dose shown in colour (% of Rx, default 50)
%   isoPct       isodose lines (default [50 95 107])
%   letWindow    [lo hi] keV/um (default [1 8]); letMinDoseGy (default 0.1*rx)
%   diffWindow   [-a a] for 'diff' panels (default [-5 5]), diffUnit label
%   contours     struct array: name, color (RGB), label ('' = none), labelOffsetMm [right down]
%   crop         structure name, or a cell list (box = union), default 'Brain' (else 'Body'); cropMarginMm (10)
%   ctWindowHU   (default [-40 120]: soft tissue / brain; bone saturates white)
%   pixelMm      resampled pixel size (default 0.25)
%   background   'dark' (TPS look, default) | 'light' (journal print)
%   widthCm      figure width (default 17.5)
%
% Rendering: CT and maps are resampled to pixelMm (cubic for CT, linear for
% dose), the colour wash is alpha-blended into the grey CT in RGB (identical
% in MATLAB and Octave, prints correctly), structures are smoothed outlines
% drawn as vector lines, and isodose lines are vector contours.
% Output: base.pdf, base.png (600 dpi) and base_slices.mat (the native-resolution
% 2D slices that are plotted, for re-plotting and checking).
% PATIENT ANATOMY: keep the output on hospital storage. A sagittal head slice
% can show the face profile: the default crop to the brain removes it; check
% consent / journal rules before publishing any patient image.

if nargin < 4 || isempty(opt), opt = struct(); end
g = @(f, d) cns_getfield(opt, f, d);
plane = g('plane', 'sagittal');
rx = g('rxGy', 54);
washMin = g('washMinPct', 50);
isoPct = g('isoPct', [50 95 107]);
letWin = g('letWindow', [1 8]);
letMinDose = g('letMinDoseGy', 0.1 * rx);
diffWin = g('diffWindow', [-5 5]);
ctWin = g('ctWindowHU', [-40 120]);
px = g('pixelMm', 0.25);
dark = strcmp(g('background', 'dark'), 'dark');
fg = [1 1 1] * dark + [0 0 0] * ~dark;
bg = [0 0 0] * dark + [1 1 1] * ~dark;
contoursOpt = g('contours', struct('name', {}, 'color', {}, 'label', {}, 'labelOffsetMm', {}));
contoursOpt = cns_toCell(contoursOpt);

% ---- slice index
if isempty(g('slice', []))
    iC = find(strcmp(cst(:,2), g('centreOn', 'CTV')), 1);
    [r, c, s] = ind2sub(ct.cubeDim, cst{iC,4}{1});
    switch plane
        case 'sagittal', k = round(median(c));
        case 'coronal',  k = round(median(r));
        otherwise,       k = round(median(s));
    end
else
    k = opt.slice;
end
HU = double(ct.cubeHU{1});
[ctS, dr, dc, ori] = getSlice(HU, k, plane, ct.resolution);

% ---- crop box (mm) from the crop structure on this slice
cropNames = cns_toCell(g('crop', 'Brain'));      % one name or a list: crop box = union of them
cropNames = cropNames(cellfun(@(n) any(strcmp(cst(:,2), n)), cropNames));
if isempty(cropNames), cropNames = {'Body'}; end
box = [];
m = false(ct.cubeDim);
for ic = 1:numel(cropNames)
    ib = find(strcmp(cst(:,2), cropNames{ic}), 1);
    if ~isempty(ib), m(cst{ib,4}{1}) = true; end
end
if any(m(:))
    mS = getSlice(double(m), k, plane, ct.resolution) > 0.5;
    [rr, cc] = find(mS);
    if ~isempty(rr)
        mg = g('cropMarginMm', 10);
        box = [(min(cc) - 1) * dc - mg, (max(cc) - 1) * dc + mg, (min(rr) - 1) * dr - mg, (max(rr) - 1) * dr + mg];
    end
end

% ---- resample CT
[ctU, xq, yq] = up(ctS, dr, dc, px, 'cubic');
grey = min(max((ctU - ctWin(1)) / diff(ctWin), 0), 1);

% ---- layout
nRow = max([panels.row]); nCol = max([panels.col]);
if isempty(box), box = [0 xq(end) 0 yq(end)]; end
aspect = (box(4) - box(3)) / (box(2) - box(1));
W = g('widthCm', 17.5);
% colour bars in a strip at the right: one per panel type in a row (e.g. Delta dose and LETd)
nTypeMax = 1;
for r = 1:max([panels.row])
    nTypeMax = max(nTypeMax, numel(unique({panels([panels.row] == r).type})));
end
cbW = (2.2 + 1.7 * (nTypeMax - 1)) / W;       % fraction of width for colour bars (2.2 cm: bar, ticks, label)
gapW = 0.15;                                  % gap between panels (cm)
panW = W * (1 - cbW) / nCol - gapW;           % panel width (cm)
panH = panW * aspect;                         % panel height (cm): the crop box fills the axes exactly
tH = 0.60; bH = 0.15;                         % title strip above and gap below every panel (cm)
rowH = panH + tH + bH;
H = nRow * rowH;
f = cns_journalFigure(W, H);
set(f, 'Color', bg, 'InvertHardcopy', 'off');
try, set(f, 'Theme', 'light'); catch, end %#ok<NOCOM>
saved = struct('plane', plane, 'slice', k, 'pixelMm', [dr dc], 'ct', ctS, 'panels', {{}});
rowExtra = struct('row', {}, 'col', {}, 'type', {});   % other types in the same row (Nov 8)
used = false(nRow, nCol);

for p = 1:numel(panels)
    P = panels(p);
    [S, ~, ~] = getSlice(double(P.cube), k, plane, ct.resolution);
    saved.panels{end+1} = struct('title', P.title, 'type', P.type, 'slice', S);
    U = up(S, dr, dc, px, 'linear');
    rgb = repmat(grey, [1 1 3]);
    switch P.type
        case 'dose'
            pct = 100 * U / rx;
            [col, a] = doseColour(pct, washMin);
            pType = 'dose';
        case 'let'
            Dm = up(getSlice(double(P.doseForMask), k, plane, ct.resolution), dr, dc, px, 'linear');
            t = (U - letWin(1)) / diff(letWin);
            col = mapColour(t, letMap());
            a = 0.70 * (Dm >= letMinDose & U > 0);
            pType = 'let';
        case 'diff'
            Dm = up(getSlice(double(P.doseForMask), k, plane, ct.resolution), dr, dc, px, 'linear');
            t = (U - diffWin(1)) / diff(diffWin);
            col = mapColour(t, divMap());
            a = 0.75 * min(1, abs(U) / (0.1 * max(abs(diffWin)))) .* (Dm >= letMinDose);
            pType = 'diff';
    end
    used(P.row, P.col) = true;
    rowExtra(end+1) = struct('row', P.row, 'col', P.col, 'type', pType); %#ok<AGROW>
    a3 = repmat(a, [1 1 3]);
    rgb = (1 - a3) .* rgb + a3 .* col;

    x0 = (0.10 + (P.col - 1) * (panW + gapW)) / W;
    y0 = (H - P.row * rowH + bH) / H;
    ax = axes('Parent', f, 'Position', [x0, y0, panW / W, panH / H]);
    image(ax, xq, yq, rgb); hold(ax, 'on');
    axis(ax, 'image'); set(ax, 'YDir', 'reverse', 'Visible', 'off');
    xlim(ax, box(1:2)); ylim(ax, box(3:4));

    if strcmp(P.type, 'dose')                 % isodose lines
        for L = isoPct(:)'
            drawLevel(ax, xq, yq, pct, L, doseColour(L, 0), 0.8);
        end
    end
    for j = 1:numel(contoursOpt)              % structure outlines + labels
        C = contoursOpt{j};
        i = find(strcmp(cst(:,2), C.name), 1);
        if isempty(i) || isempty(cst{i,4}{1}), continue; end
        m = false(ct.cubeDim); m(cst{i,4}{1}) = true;
        mS = getSlice(double(m), k, plane, ct.resolution);
        if ~any(mS(:)), continue; end
        mU = up(smooth2(mS, 0.7), dr, dc, px, 'linear');
        [xc, yc] = drawLevel(ax, xq, yq, mU, 0.5, C.color, 1.1);
        lab = cns_getfield(C, 'label', '');
        if ~isempty(lab) && P.col == 1 && ~isempty(xc)
            [~, iy] = min(yc); xa = xc(iy); ya = yc(iy);          % arrow to the top of the outline
            off = cns_getfield(C, 'labelOffsetMm', [15 -15]);
            xt = xa + off(1); yt = ya + off(2);
            plot(ax, [xt xa], [yt ya], '-', 'Color', C.color, 'LineWidth', 1.0);
            plot(ax, xa, ya, 'o', 'Color', C.color, 'MarkerFaceColor', C.color, 'MarkerSize', 2.5);
            text(ax, xt, yt, lab, 'Color', C.color, 'FontSize', 8, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
    % orientation letters and scale bar
    text(ax, box(1) + 0.03 * diff(box(1:2)), box(3) + 0.05 * diff(box(3:4)), ori{1}, ...
        'Color', fg, 'FontSize', 7, 'FontWeight', 'bold');
    text(ax, box(2) - 0.03 * diff(box(1:2)), box(3) + 0.05 * diff(box(3:4)), ori{2}, ...
        'Color', fg, 'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    if P.row == nRow && P.col == 1
        xb = box(1) + 0.06 * diff(box(1:2)); yb = box(4) - 0.05 * diff(box(3:4));
        plot(ax, [xb xb + 20], [yb yb], '-', 'Color', fg, 'LineWidth', 1.5);
        text(ax, xb + 10, yb - 0.015 * diff(box(3:4)), '2 cm', 'Color', fg, 'FontSize', 6.5, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    title(ax, P.title, 'Color', fg, 'FontSize', 8, 'FontWeight', 'normal', 'Visible', 'on', 'Interpreter', 'tex');
end

% ---- colour bars: one per panel type in each row, in the right strip, in the left-to-right
% order of the panels (rows with two bars label them along the bar instead of on top)
bars = struct('row', {}, 'pos', {}, 'type', {}, 'nInRow', {});
for r = 1:nRow
    E = rowExtra([rowExtra.row] == r);
    if isempty(E), continue; end
    [~, o] = sort([E.col]); E = E(o);
    types = {};
    for e = 1:numel(E), if ~any(strcmp(types, E(e).type)), types{end+1} = E(e).type; end, end %#ok<AGROW>
    y = (H - r * rowH + bH + 0.25) / H; h = (panH - 0.50) / H;
    for t = 1:numel(types)
        bars(end+1) = struct('row', r, 'pos', [1 - cbW + (0.30 + 1.7 * (t - 1)) / W, y, 0.30 / W, h], ...
            'type', types{t}, 'nInRow', numel(types)); %#ok<AGROW>
    end
end
for b = 1:numel(bars)
    ax = axes('Parent', f, 'Position', bars(b).pos);
    switch bars(b).type
        case 'dose'
            stops = [50 90 95 107 120 140 160];
            stops = stops(stops >= washMin);
            n = 64; img = zeros(n * (numel(stops) - 1), 1, 3);
            for s = 1:numel(stops) - 1
                v = linspace(stops(s), stops(s + 1), n);
                img((s - 1) * n + (1:n), 1, :) = reshape(doseColour(v(:), 0), [n 1 3]);
            end
            image(ax, [0 1], [0 numel(stops) - 1], img); set(ax, 'YDir', 'normal');
            set(ax, 'YTick', 0:numel(stops) - 1, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), stops, 'UniformOutput', false));
            lab = sprintf('%% of %g cGy (RBE)', 100 * rx);
        case 'let'
            img = reshape(mapColour(linspace(0, 1, 256)', letMap()), [256 1 3]);
            image(ax, [0 1], letWin, img); set(ax, 'YDir', 'normal');
            lab = 'LET_d (keV/\mum)';
        case 'diff'
            img = reshape(mapColour(linspace(0, 1, 256)', divMap()), [256 1 3]);
            image(ax, [0 1], diffWin, img); set(ax, 'YDir', 'normal');
            lab = g('diffUnit', '\Delta');
    end
    set(ax, 'XTick', [], 'YAxisLocation', 'right', 'FontSize', 7, 'XColor', fg, 'YColor', fg, ...
        'TickDir', 'out', 'Box', 'off', 'Color', bg);
    ylabel(ax, lab, 'Color', fg, 'FontSize', 7, 'Interpreter', 'tex');
end

d = fileparts(base);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
save([base '_slices.mat'], 'saved');
isOctave = exist('OCTAVE_VERSION', 'builtin') > 0;
try
    if ~isOctave && exist('exportgraphics', 'file')
        exportgraphics(f, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', bg);
        exportgraphics(f, [base '.png'], 'Resolution', 600, 'BackgroundColor', bg);
    else
        print(f, [base '.pdf'], '-dpdf', '-painters');
        print(f, [base '.png'], '-dpng', '-r600');
    end
catch err
    warning('cns:figure', 'Could not save %s: %s', base, err.message);
end
close(f);
end

% =========================================================================
function [S, dr, dc, ori] = getSlice(A, k, plane, res)
% 2D slice for display: rows top->bottom, columns left->right, spacing dr/dc (mm)
switch plane
    case 'axial'        % anterior at the top, patient left on the right (radiological)
        S = A(:, :, k); dr = res.y; dc = res.x; ori = {'R', 'L'};
    case 'sagittal'     % superior at the top, anterior on the left
        S = flipud(squeeze(A(:, k, :))'); dr = res.z; dc = res.y; ori = {'A', 'P'};
    case 'coronal'      % superior at the top, patient left on the right
        S = flipud(squeeze(A(k, :, :))'); dr = res.z; dc = res.x; ori = {'R', 'L'};
    otherwise
        error('plane must be axial, sagittal or coronal');
end
end

function [U, xq, yq] = up(S, dr, dc, px, method)
y = (0:size(S, 1) - 1) * dr; x = (0:size(S, 2) - 1) * dc;
yq = 0:px:y(end); xq = 0:px:x(end);
[Xq, Yq] = meshgrid(xq, yq);
U = interp2(x, y, double(S), Xq, Yq, method);
U(~isfinite(U)) = 0;
end

function M = smooth2(M, sigmaPx)
r = ceil(3 * sigmaPx);
k = exp(-0.5 * ((-r:r) / sigmaPx).^2); k = k / sum(k);
M = conv2(k, k, double(M), 'same');
end

function [xc, yc] = drawLevel(ax, xq, yq, Z, level, colour, lw)
xc = []; yc = [];
C = contourc(xq, yq, double(Z), [level level]);
i = 1;
while i < size(C, 2)
    n = C(2, i);
    xs = C(1, i + 1:i + n); ys = C(2, i + 1:i + n);
    if n > 8
        plot(ax, xs, ys, '-', 'Color', colour, 'LineWidth', lw);
        xc = [xc, xs]; yc = [yc, ys]; %#ok<AGROW>
    end
    i = i + n + 1;
end
end

function [col, a] = doseColour(pct, washMin)
% TPS-style colour wash: blue (50 %) - cyan - green - yellow (90 %) - red (95-100 %)
% - dark red (107 %) - purple (120 %) - magenta (140 %) - pale pink (160 %)
A = [ 30 0.05 0.10 0.60;  50 0.10 0.35 1.00;  70 0.00 0.80 0.95;  80 0.25 0.85 0.30;
      90 1.00 0.92 0.00;  95 1.00 0.45 0.00; 100 0.90 0.05 0.05; 107 0.62 0.00 0.25;
     120 0.50 0.10 0.70; 140 0.85 0.35 0.90; 160 1.00 0.80 1.00];
sz = size(pct);
v = min(max(pct(:), A(1, 1)), A(end, 1));
col = [interp1(A(:,1), A(:,2), v), interp1(A(:,1), A(:,3), v), interp1(A(:,1), A(:,4), v)];
if numel(sz) == 2 && sz(2) > 1
    col = reshape(col, [sz 3]);
end
if nargout > 1
    a = 0.30 + 0.35 * min(max((pct - washMin) / max(95 - washMin, 1), 0), 1);
    a(pct < washMin) = 0;
end
end

function M = letMap()
% perceptually ordered (inferno-like) for LETd
M = [0 0.00 0.00 0.30; 0.15 0.20 0.04 0.45; 0.35 0.55 0.10 0.50; 0.55 0.85 0.25 0.25;
     0.75 0.98 0.55 0.05; 0.90 0.99 0.80 0.20; 1 0.99 0.99 0.65];
end

function M = divMap()
% diverging blue - white - red
M = [0 0.13 0.30 0.70; 0.25 0.40 0.62 0.90; 0.5 0.97 0.97 0.97; 0.75 0.95 0.55 0.40; 1 0.70 0.08 0.10];
end

function col = mapColour(t, M)
sz = size(t);
t = min(max(t(:), 0), 1);
col = [interp1(M(:,1), M(:,2), t), interp1(M(:,1), M(:,3), t), interp1(M(:,1), M(:,4), t)];
if numel(sz) == 2 && sz(2) > 1
    col = reshape(col, [sz 3]);
end
end
