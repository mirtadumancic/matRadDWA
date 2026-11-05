function cns_figDVHandBars(dvh, rows, machines, dvhStructs, barMetrics, rx, base, opt)
% cns_figDVHandBars  Journal figure: (a) DVHs, [(b) LVHs,] (c) dose-LET bars per plan.
%
% opt.lvh, opt.lvhStructs, opt.letdMinDoseGy: add the LET-volume histogram row
% (volume % of the structure with dose >= letdMinDoseGy AND LETd >= x).
%
%   cns_figDVHandBars(dvh, rows, {'DWA_proton_v2','Generic','IBA_cyclotron'}, ...
%       {'CTV','OpticNerve_L','Chiasm','Brainstem','Brain_healthy'}, barMetrics, 54, base)
%
% dvh        struct array (machine, plan, structure, quantity, dose, volPct); photon
%            rows have machine 'photon', plan 'clinical'
% rows       long metric rows {alias, machine, plan, structure, quantity, metric, value, unit}
% barMetrics struct array: structure, quantity, metric, label (tex)
% Plans per machine: 'doseOnly' (solid line / light bar) and 'LETxD_best'
% (dashed line / full-colour bar, with the % change vs dose-only printed above).
% Photon: black line in (a); dashed black reference line in (b) where it is
% defined (dose metrics only; LET metrics are not defined for photons).
% Saves base.pdf, base.png (300 dpi), base_DVH.csv and base_bars.csv.
if nargin < 8, opt = struct(); end
plans = {'doseOnly', 'LETxD_best'};
nA = numel(dvhStructs); nB = numel(barMetrics);
lvh = cns_getfield(opt, 'lvh', []);
lvhS = cns_toCell(cns_getfield(opt, 'lvhStructs', {}));
doLVH = ~isempty(lvh) && ~isempty(lvhS);
nC = max([nA, nB, numel(lvhS)]);
W = 17.5;
if doLVH
    H = 19; yA = 0.75; hA = 0.19; legY = 0.625; yL = 0.40; hL = 0.17; yB = 0.075; hB = 0.17;
    labY = [0.955 0.60 0.29];
else
    H = 13.5; yA = 0.64; hA = 0.29; legY = 0.46; yB = 0.10; hB = 0.25;
    labY = [0.94 NaN 0.40];
end
f = cns_journalFigure(W, H);
mL = 0.065; mR = 0.012; gX = 0.04;
wAx = (1 - mL - mR - (nC - 1) * gX) / nC;

% ---------------- (a) DVHs
csvA = {};
hLeg = []; labLeg = {};
for k = 1:nA
    ax = axes('Parent', f, 'Position', [mL + (k - 1) * (wAx + gX), yA, wAx, hA]);
    styleAxes(ax); hold(ax, 'on');
    sel = find(strcmp({dvh.machine}, 'photon') & strcmp({dvh.structure}, dvhStructs{k}) & strcmp({dvh.quantity}, 'RBE11'), 1);
    if ~isempty(sel)
        h = plot(ax, dvh(sel).dose, dvh(sel).volPct, '-', 'Color', [0 0 0], 'LineWidth', 1.3);
        if k == 1, hLeg(end+1) = h; labLeg{end+1} = 'Photon (delivered)'; end %#ok<AGROW>
        csvA = addCsv(csvA, dvhStructs{k}, 'photon', 'clinical', dvh(sel));
    end
    for m = 1:numel(machines)
        for p = 1:2
            sel = find(strcmp({dvh.machine}, machines{m}) & strcmp({dvh.plan}, plans{p}) & ...
                strcmp({dvh.structure}, dvhStructs{k}) & strcmp({dvh.quantity}, 'RBE11'), 1);
            if isempty(sel), continue; end
            ls = '-'; if p == 2, ls = '--'; end
            h = plot(ax, dvh(sel).dose, dvh(sel).volPct, ls, 'Color', cns_okabeIto(machines{m}), 'LineWidth', 1.0);
            if k == 1
                hLeg(end+1) = h; %#ok<AGROW>
                labLeg{end+1} = sprintf('%s %s', shortName(machines{m}), ifelse(p == 1, 'dose-only', 'LET\timesD')); %#ok<AGROW>
            end
            csvA = addCsv(csvA, dvhStructs{k}, machines{m}, plans{p}, dvh(sel));
        end
    end
    if any(strcmp(dvhStructs{k}, {'CTV', 'PTV', 'GTV'}))
        xlim(ax, [0.85 1.1] * rx);
    else
        xlim(ax, [0 1.1 * rx]);
    end
    ylim(ax, [0 102]);
    title(ax, strrep(dvhStructs{k}, '_', ' '), 'FontWeight', 'normal', 'FontSize', 7.5, 'Color', 'k');
    xlabel(ax, 'Dose (Gy(RBE))', 'FontSize', 7);
    if k == 1, ylabel(ax, 'Volume (%)', 'FontSize', 7); end
end
lg = legend(hLeg, labLeg, 'Interpreter', 'tex', 'FontSize', 6.5);
try, set(lg, 'NumColumns', 4); catch, set(lg, 'Orientation', 'horizontal'); end %#ok<NOCOM>
try, set(lg, 'Box', 'off', 'Units', 'normalized', 'Position', [0.08 legY 0.86 0.9 / H], 'TextColor', 'k', 'Color', 'w'); catch, end %#ok<NOCOM>
annotation(f, 'textbox', [0.005 labY(1) 0.04 0.04], 'String', '(a)', 'EdgeColor', 'none', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');

% ---------------- (b) LVHs (proton plans only)
csvL = {};
if doLVH
    for k = 1:numel(lvhS)
        ax = axes('Parent', f, 'Position', [mL + (k - 1) * (wAx + gX), yL, wAx, hL]);
        styleAxes(ax); hold(ax, 'on');
        for m = 1:numel(machines)
            for p = 1:2
                sel = find(strcmp({lvh.machine}, machines{m}) & strcmp({lvh.plan}, plans{p}) & ...
                    strcmp({lvh.structure}, lvhS{k}), 1);
                if isempty(sel), continue; end
                ls = '-'; if p == 2, ls = '--'; end
                x = lvh(sel).letd; y = lvh(sel).volPct;
                last = find(y > 0, 1, 'last'); if isempty(last), last = 1; end
                plot(ax, x(1:min(end, last + 1)), y(1:min(end, last + 1)), ls, 'Color', cns_okabeIto(machines{m}), 'LineWidth', 1.0);
                for i = 1:min(numel(x), last + 1)
                    csvL(end+1, :) = {lvhS{k}, machines{m}, plans{p}, x(i), y(i)}; %#ok<AGROW>
                end
            end
        end
        xlim(ax, [0 10]); yl = ylim(ax); ylim(ax, [0 max(yl(2), 1)]);
        title(ax, strrep(lvhS{k}, '_', ' '), 'FontWeight', 'normal', 'FontSize', 7.5, 'Color', 'k');
        xlabel(ax, 'LET_d (keV/\mum)', 'FontSize', 7);
        if k == 1, ylabel(ax, 'Volume (%)', 'FontSize', 7); end
    end
    annotation(f, 'textbox', [0.005 labY(2) 0.04 0.04], 'String', '(b)', 'EdgeColor', 'none', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
end

% ---------------- (c) bars
csvB = {};
for k = 1:nB
    B = barMetrics(k);
    ax = axes('Parent', f, 'Position', [mL + (k - 1) * (wAx + gX), yB, wAx, hB]);
    styleAxes(ax); hold(ax, 'on'); grid(ax, 'off'); set(ax, 'YGrid', 'on');
    vmax = 0;
    allV = [];
    for m = 1:numel(machines)
        allV = [allV, get1(rows, machines{m}, 'doseOnly', B), get1(rows, machines{m}, 'LETxD_best', B)]; %#ok<AGROW>
    end
    panelMax = max([allV(~isnan(allV)), eps]);
    for m = 1:numel(machines)
        v = [get1(rows, machines{m}, 'doseOnly', B), get1(rows, machines{m}, 'LETxD_best', B)];
        c = cns_okabeIto(machines{m});
        x = (m - 1) * 3 + [1 2];
        if ~isnan(v(1)), drawBar(ax, x(1), v(1), 1 - 0.45 * (1 - c), c); end
        if ~isnan(v(2)), drawBar(ax, x(2), v(2), c, c); end
        if all(~isnan(v))
            if v(1) >= 0.1 * panelMax          % % change only on a meaningful baseline
                lab = sprintf('%+.0f%%', 100 * (v(2) - v(1)) / v(1));
            else
                lab = sprintf('%+.2g', v(2) - v(1));
            end
            text(ax, x(2), v(2), lab, 'FontSize', 6, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
        end
        vmax = max([vmax, v(~isnan(v))]);
        csvB(end+1, :) = {B.structure, B.quantity, B.metric, machines{m}, v(1), v(2)}; %#ok<AGROW>
    end
    vp = get1(rows, 'photon', 'clinical', B);
    if ~isnan(vp) && ~any(strcmp(B.quantity, {'DxLETd', 'LETd', 'hotSpot', 'dirty'}))
        plot(ax, [0.3, 3 * numel(machines) - 0.3], [vp vp], 'k--', 'LineWidth', 0.9);
        vmax = max(vmax, vp);
        csvB(end+1, :) = {B.structure, B.quantity, B.metric, 'photon', vp, NaN}; %#ok<AGROW>
    end
    xlim(ax, [0.2, 3 * numel(machines) - 0.2]);
    if vmax > 0, ylim(ax, [0 1.18 * vmax]); end
    set(ax, 'XTick', (0:numel(machines) - 1) * 3 + 1.5, 'XTickLabel', cellfun(@shortName, machines, 'UniformOutput', false));
    title(ax, {B.label, cns_getfield(B, 'sublabel', '')}, 'FontWeight', 'normal', 'FontSize', 6.5, 'Color', 'k', 'Interpreter', 'tex');
end
annotation(f, 'textbox', [0.005 labY(3) 0.04 0.04], 'String', ifelse(doLVH, '(c)', '(b)'), 'EdgeColor', 'none', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
note = 'light bars: dose-only;  full bars: LET\timesD-optimised (% change vs dose-only; absolute change when dose-only < 10% of the panel maximum);  dashed line: photon (delivered)';
if doLVH
    note = sprintf('%s.  (b) LET_d in voxels with dose \\geq %g Gy(RBE)', note, cns_getfield(opt, 'letdMinDoseGy', 5.4));
end
annotation(f, 'textbox', [0.03 0.0 0.94 0.035], 'String', note, ...
    'EdgeColor', 'none', 'FontSize', 6.5, 'Color', 'k', 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

cns_saveJournalFigure(f, base, {}, {});
cns_writeCsvTable([base '_DVH.csv'], {'structure', 'machine', 'plan', 'doseGyRBE', 'volPct'}, csvA);
cns_writeCsvTable([base '_bars.csv'], {'structure', 'quantity', 'metric', 'machine', 'doseOnly', 'LETxD_best'}, csvB);
if doLVH
    cns_writeCsvTable([base '_LVH.csv'], {'structure', 'machine', 'plan', 'letdKeVum', 'volPct'}, csvL);
end
end

% -------------------------------------------------------------------------
function drawBar(ax, x, v, face, edge)
patch(ax, [x - 0.42, x + 0.42, x + 0.42, x - 0.42], [0 0 v v], face, 'EdgeColor', edge, 'LineWidth', 0.6);
end

function styleAxes(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 7, 'TickDir', 'out', 'LineWidth', 0.6, 'Color', 'w', ...
    'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1], 'Box', 'off');
grid(ax, 'on');
try, set(ax, 'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.15); catch, end %#ok<NOCOM>
end

function c = addCsv(c, s, m, p, d)
for i = 1:numel(d.dose)
    c(end+1, :) = {s, m, p, d.dose(i), d.volPct(i)}; %#ok<AGROW>
end
end

function v = get1(rows, machine, plan, B)
sel = strcmp(rows(:,2), machine) & strcmp(rows(:,3), plan) & strcmp(rows(:,4), B.structure) & ...
    strcmp(rows(:,5), B.quantity) & strcmp(rows(:,6), B.metric);
if any(sel)
    v = rows{find(sel, 1), 7};
    if ischar(v), v = str2double(v); end
else
    v = NaN;
end
end

function s = shortName(m)
switch m
    case 'DWA_proton_v2', s = 'DWA';
    case 'IBA_cyclotron', s = 'IBA';
    otherwise, s = strrep(m, '_', ' ');
end
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end
