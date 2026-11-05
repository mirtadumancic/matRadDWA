function cns_figCurves(S, xField, yField, structures, series, xLabel, yLabel, xLimits, base, titleStr)
% cns_figCurves  Journal-style grid of DVH or LVH panels, one panel per structure.
%
%   cns_figCurves(dvh, 'dose', 'volPct', {'CTV','Chiasm'}, series, 'Dose (Gy(RBE))', ...
%                 'Volume (%)', [0 60], 'D:/runs/CNS_01/figures/fig_dvh', '')
%
% S        struct array with fields machine, plan, structure, [quantity], xField, yField
%          (the dvh / lvh output of cns_runCase / cns_evaluatePlan)
% series   struct array: machine, plan, quantity ('' = any), label, color, style
% xLimits  [lo hi], or a function handle @(structure) -> [lo hi]
% Writes base.pdf, base.png (300 dpi) and base.csv (every plotted point).
nS = numel(structures);
nCol = min(nS, 5);
nRow = ceil(nS / nCol);
nLegCols = 4;
nLegRows = max(1, ceil(numel(unique({series.label})) / nLegCols));
legH = 0.45 * nLegRows;                           % cm for the legend
hFig = 3.6 * nRow + 1.5 + legH;                   % cm
f = cns_journalFigure(17.5, hFig);
% explicit layout: room for y labels on the left and the legend at the bottom
mL = 0.07; mR = 0.015; mT = 0.6 / hFig; mB = (1.2 + legH) / hFig; gX = 0.035; gY = 1.25 / hFig;
wAx = (1 - mL - mR - (nCol - 1) * gX) / nCol;
hAx = (1 - mT - mB - (nRow - 1) * gY) / nRow;
csv = {};
hLeg = []; labLeg = {};
for k = 1:nS
    rr = ceil(k / nCol); cc = k - (rr - 1) * nCol;
    ax = axes('Parent', f, 'Position', [mL + (cc - 1) * (wAx + gX), 1 - mT - rr * hAx - (rr - 1) * gY, wAx, hAx]);
    hold(ax, 'on'); grid(ax, 'on');
    set(ax, 'FontName', 'Arial', 'FontSize', 7, 'TickDir', 'out', 'LineWidth', 0.6, ...
        'Color', 'w', 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    try, set(ax, 'GridColor', [0.15 0.15 0.15]); catch, end %#ok<NOCOM>
    try, set(ax, 'GridAlpha', 0.15); catch, end %#ok<NOCOM>
    for j = 1:numel(series)
        sr = series(j);
        sel = strcmp({S.machine}, sr.machine) & strcmp({S.plan}, sr.plan) & strcmp({S.structure}, structures{k});
        if ~isempty(sr.quantity) && isfield(S, 'quantity')
            sel = sel & strcmp({S.quantity}, sr.quantity);
        end
        i = find(sel, 1);
        if isempty(i), continue; end
        x = S(i).(xField); y = S(i).(yField);
        last = find(y > 0, 1, 'last'); if isempty(last), last = 1; end
        x = x(1:min(numel(x), last + 1)); y = y(1:min(numel(y), last + 1));
        h = plot(ax, x, y, sr.style, 'Color', sr.color, 'LineWidth', 1.0);
        if k == 1 || ~any(strcmp(labLeg, sr.label))
            if ~any(strcmp(labLeg, sr.label)), hLeg(end+1) = h; labLeg{end+1} = sr.label; end %#ok<AGROW>
        end
        for p = 1:numel(x)
            csv(end+1, :) = {structures{k}, sr.machine, sr.plan, sr.quantity, x(p), y(p)}; %#ok<AGROW>
        end
    end
    if isa(xLimits, 'function_handle'), xl = xLimits(structures{k}); else, xl = xLimits; end
    if ~isempty(xl), xlim(ax, xl); set(ax, 'XTick', niceTicks(xl)); end
    yl = ylim(ax); ylim(ax, [0 max(yl(2), eps)]);
    title(ax, strrep(structures{k}, '_', ' '), 'FontWeight', 'normal', 'FontSize', 7, 'Color', 'k');
    if cc == 1, ylabel(ax, yLabel, 'FontSize', 7); end
    if k > nS - nCol, xlabel(ax, xLabel, 'FontSize', 7); end
end
if ~isempty(hLeg)
    lg = legend(hLeg, labLeg, 'Interpreter', 'tex', 'FontSize', 7);
    try, set(lg, 'NumColumns', min(nLegCols, numel(hLeg))); catch, set(lg, 'Orientation', 'horizontal'); end %#ok<NOCOM>
    try, set(lg, 'Box', 'off', 'Units', 'normalized', 'Position', [0.05 0.01 0.9 legH / hFig]); catch, end %#ok<NOCOM>
    try, set(lg, 'TextColor', 'k', 'Color', 'w'); catch, end %#ok<NOCOM>
end
if nargin > 9 && ~isempty(titleStr)
    annotation(f, 'textbox', [0 0.96 1 0.04], 'String', titleStr, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Interpreter', 'none');
end
cns_saveJournalFigure(f, base, {'structure', 'machine', 'plan', 'quantity', xField, yField}, csv);
end

function t = niceTicks(xl)
% 3-5 round tick values inside xl
span = xl(2) - xl(1);
steps = [0.5 1 2 2.5 5 10 20 25 50 100];
st = steps(find(span ./ steps <= 5, 1));
if isempty(st), st = span / 4; end
t = (ceil(xl(1) / st):floor(xl(2) / st)) * st;
end
