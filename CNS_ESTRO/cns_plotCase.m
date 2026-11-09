function cns_plotCase(caseDir, alias, rows, dvh, lvh, slices, caseCfg, study, dataClass)
% cns_plotCase  DVH, LVH, hot-spot and (optionally) slice figures for one case.
%
% Figures go to <caseDir>/figures. DVH/LVH/bar figures contain no images.
% Dose/LET slice maps are drawn for public and synthetic cases; for
% patient cases only if study.figures.slicesForPatient is true (they stay
% on hospital storage either way).

figDir = fullfile(caseDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
set(0, 'DefaultFigureVisible', 'off');
rx = caseCfg.prescription.totalDoseGyRBE;
F = cns_getfield(study, 'figures', struct());
oars = cns_toCell(cns_getfield(F, 'oars', {'Brainstem', 'Chiasm', 'Brain_healthy', 'BrainRing_0_15'}));
machinesAll = unique({dvh.machine}, 'stable');
refLabel = cns_getfield(caseCfg, 'reference.label', 'photon');
hasRef = any(strcmp(machinesAll, refLabel));
machines = machinesAll(~strcmp(machinesAll, refLabel));   % proton machines; the photon plan is overlaid
plans = unique({dvh.plan}, 'stable');
plans = plans(~strcmp(plans, 'clinical'));
cols = lines(max(numel(plans), 1));
showMcN = logical(cns_getfield(F, 'showMcNamara', true));      % (Nov 8) off in the robust DWA study
if ~showMcN && isfield(dvh, 'quantity'), dvh = dvh(~strcmp({dvh.quantity}, 'McN')); end

% ---- 1) target DVH: CTV and PTV, RBE 1.1 (solid) and McNamara (dashed)
tg = {'CTV', 'PTV'};
f = figure('Position', [50 50 420 * numel(machines) 640]);
for r = 1:2
    for m = 1:numel(machines)
        subplot(2, numel(machines), (r - 1) * numel(machines) + m); hold on; box on; grid on;
        [h, lab] = drawCurves(dvh, machines{m}, plans, tg{r}, 'dose', cols);
        if hasRef, [h, lab] = addRef(h, lab, dvh, refLabel, tg{r}); end
        xline0(rx);
        xlim([0.8 * rx, 1.15 * rx]); ylim([0 102]);
        title(sprintf('%s  %s', tg{r}, machines{m}), 'Interpreter', 'none');
        xlabel('Dose [Gy(RBE)]'); ylabel('Volume [%]');
        if r == 1 && m == 1 && ~isempty(h), legend(h, lab, 'Interpreter', 'none', 'Location', 'southwest'); end
    end
end
annotation('textbox', [0 0.95 1 0.05], 'String', sprintf('%s: RBE 1.1%s%s', alias, ifelse(showMcN, ' solid, McNamara dashed', ''), ...
    ifelse(hasRef, sprintf(', black = %s clinical plan', refLabel), '')), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'none');
saveFig(f, fullfile(figDir, 'dvh_targets'));

% ---- 2) OAR DVH and 3) OAR LVH
% (Nov 8) no pencil-beam LET when the dij has no LET (doseCalc.calcLET = 0, e.g. the robust
% DWA study: LETd comes from the MCsquare recalculation): skip the empty LVH and hot-spot figures
hasLET = logical(cns_getfield(study, 'doseCalc.calcLET', 1));
kinds = {'dvh', 'lvh'}; if ~hasLET, kinds = {'dvh'}; end
for kind = kinds
    f = figure('Position', [50 50 380 * numel(machines) 300 * numel(oars)]);
    for r = 1:numel(oars)
        for m = 1:numel(machines)
            subplot(numel(oars), numel(machines), (r - 1) * numel(machines) + m); hold on; box on; grid on;
            if strcmp(kind{1}, 'dvh')
                [h, lab] = drawCurves(dvh, machines{m}, plans, oars{r}, 'dose', cols);
                if hasRef, [h, lab] = addRef(h, lab, dvh, refLabel, oars{r}); end
                xlabel('Dose [Gy(RBE)]'); ylabel('Volume [%]'); xlim([0 1.15 * rx]);
            else
                [h, lab] = drawCurves(lvh, machines{m}, plans, oars{r}, 'letd', cols);
                xlabel(sprintf('LETd [keV/\\mum] (dose \\geq %g Gy(RBE))', ...
                    cns_getfield(study, 'evaluation.letdMinDoseGy', 1)));
                ylabel('Volume [cc]'); xlim([0 10]);
            end
            title(sprintf('%s  %s', oars{r}, machines{m}), 'Interpreter', 'none');
            if r == 1 && m == 1 && ~isempty(h), legend(h, lab, 'Interpreter', 'none', 'Location', 'northeast'); end
        end
    end
    saveFig(f, fullfile(figDir, [kind{1} '_oars']));
end

% ---- 4) dose-LET hot spots
hs = cns_toCell(cns_getfield(F, 'hotSpotBars', {}));
if ~isempty(hs) && hasLET
    f = figure('Position', [50 50 360 * numel(hs) 360]);
    for k = 1:numel(hs)
        subplot(1, numel(hs), k);
        V = nan(numel(plans), numel(machines));
        for p = 1:numel(plans)
            for m = 1:numel(machines)
                V(p, m) = lookup(rows, machines{m}, plans{p}, hs{k}.structure, 'hotSpot', hs{k}.metric);
            end
        end
        bar(V); grid on;
        set(gca, 'XTick', 1:numel(plans), 'XTickLabel', plans, 'TickLabelInterpreter', 'none');
        try, xtickangle(45); catch, end %#ok<NOCOM>
        ylabel('Volume [cc]');
        title(sprintf('%s %s', hs{k}.structure, strrep(hs{k}.metric, '_', ' ')), 'Interpreter', 'none');
        if k == 1, legend(machines, 'Interpreter', 'none', 'Location', 'northeast'); end
    end
    saveFig(f, fullfile(figDir, 'hotspots'));
end

% ---- 5) slice maps (dose-only vs each LET plan, iso slice through the CTV)
doSlices = ~strcmp(dataClass, 'patient') || logical(cns_getfield(F, 'slicesForPatient', false));
if doSlices && ~isempty(slices) && isfield(slices(1).masks, 'CTV')
    for m = 1:numel(machinesAll)
        sl = slices(strcmp({slices.machine}, machinesAll{m}));
        if isempty(sl) || ~isfield(sl(1).masks, 'CTV'), continue; end
        [~, ~, zz] = ind2sub(size(sl(1).masks.CTV), find(sl(1).masks.CTV));   % own grid per machine
        z = round(median(zz));
        f = figure('Position', [50 50 300 * numel(sl) 560]);
        for p = 1:numel(sl)
            for r = 1:2
                subplot(2, numel(sl), (r - 1) * numel(sl) + p);
                if r == 1
                    img = sl(p).cubes.RBE11(:, :, z); clim0 = [0 1.1 * rx]; lab = 'Gy(RBE)';
                else
                    img = sl(p).cubes.LETd(:, :, z); clim0 = [0 8]; lab = 'LETd keV/\mum';
                end
                imagesc(img, clim0); axis image off; colormap(gca, jet); hold on;
                contourMask(sl(p).masks, 'CTV', z, 'w');
                contourMask(sl(p).masks, 'Brainstem', z, 'm');
                cb = colorbar; ylabel(cb, lab);
                title(sprintf('%s', sl(p).plan), 'Interpreter', 'none');
            end
        end
        saveFig(f, fullfile(figDir, ['slices_' machinesAll{m}]));
    end
end
set(0, 'DefaultFigureVisible', 'on');
end

% =========================================================================
function [h, lab] = drawCurves(S, machine, plans, structure, xField, cols)
% one colour per plan; RBE 1.1 solid, McNamara dashed (legend shows the solid line)
h = []; lab = {};
for p = 1:numel(plans)
    sel = find(strcmp({S.machine}, machine) & strcmp({S.plan}, plans{p}) & strcmp({S.structure}, structure));
    hp = [];
    for k = sel
        if isfield(S, 'quantity') && strcmp(S(k).quantity, 'McN'), ls = '--'; else, ls = '-'; end
        if isfield(S, 'volCc'), y = S(k).volCc; else, y = S(k).volPct; end
        hh = plot(S(k).(xField), y, ls, 'Color', cols(p, :), 'LineWidth', 1.3);
        if ~strcmp(ls, '--'), hp = hh; end
    end
    if ~isempty(hp), h(end+1) = hp; lab{end+1} = plans{p}; end %#ok<AGROW>
end
end

function [h, lab] = addRef(h, lab, S, refLabel, structure)
sel = find(strcmp({S.machine}, refLabel) & strcmp({S.structure}, structure) & strcmp({S.quantity}, 'RBE11'), 1);
if isempty(sel), return; end
h(end+1) = plot(S(sel).dose, S(sel).volPct, '-', 'Color', [0 0 0], 'LineWidth', 1.6);
lab{end+1} = [refLabel ' clinical'];
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end

function xline0(x)
yl = ylim;
plot([x x], [0 102], 'k:');
ylim(yl);
end

function v = lookup(rows, machine, plan, structure, quantity, metric)
v = NaN;
if isempty(rows), return; end
sel = strcmp(rows(:,2), machine) & strcmp(rows(:,3), plan) & strcmp(rows(:,4), structure) ...
    & strcmp(rows(:,5), quantity) & strcmp(rows(:,6), metric);
if any(sel), v = rows{find(sel, 1), 7}; end
end

function contourMask(masks, name, z, c)
if isfield(masks, name)
    s = double(masks.(name)(:, :, z));
    if any(s(:)), contour(s, [0.5 0.5], c, 'LineWidth', 1); end
end
end

function saveFig(f, base)
% macOS + invisible figures can report a zero-size canvas: give the figure a real
% size, then try print, exportgraphics, and finally a visible redraw.
try
    set(f, 'Units', 'pixels');
    p = get(f, 'Position');
    if p(3) < 10 || p(4) < 10, set(f, 'Position', [100 100 1200 800]); end
    set(f, 'PaperPositionMode', 'auto');
catch
end
ok = false; msg = '';
try
    print(f, [base '.png'], '-dpng', '-r150'); ok = true;
catch err
    msg = err.message;
end
if ~ok && exist('OCTAVE_VERSION', 'builtin') == 0 && exist('exportgraphics', 'file')
    try, exportgraphics(f, [base '.png'], 'Resolution', 150); ok = true; catch err, msg = err.message; end %#ok<NOCOM>
end
if ~ok
    try
        set(f, 'Visible', 'on'); drawnow;
        print(f, [base '.png'], '-dpng', '-r150'); ok = true;
        set(f, 'Visible', 'off');
    catch err
        msg = err.message;
    end
end
if ~ok, fprintf(2, '[cns] could not save %s: %s\n', base, msg); end
if exist('OCTAVE_VERSION', 'builtin') == 0
    try, savefig(f, [base '.fig']); catch, end %#ok<NOCOM>
end

close(f);
end
