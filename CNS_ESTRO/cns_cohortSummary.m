function T = cns_cohortSummary(caseDirs, studyFile, outDir)
% cns_cohortSummary  Combine metrics.csv of several cases into cohort tables.
%
%   T = cns_cohortSummary({'D:/runs/CNS_01', 'D:/runs/CNS_02'}, 'study_ESTRO2027.json', 'D:/runs/cohort')
%
% 1) Adds a pseudo-plan "LETselected" per case and machine: among the
%    LET-modified plans that keep target coverage (study.selection), the one
%    with the lowest selection metric (default: brainstem dose-LET hot spot).
% 2) For every structure / quantity / metric and machine x plan: n, median,
%    min, max over cases, and the paired difference to
%      - the reference plan of the same machine (study.summary.referencePlan)
%      - the reference machine with the same plan (study.referenceMachine)
%    as median absolute and median relative (%) difference, plus a Wilcoxon
%    signed-rank p value when n >= 5 and signrank is available.
% Writes summary_long.csv, metrics_all.csv and summary_abstract.md to outDir.
%
% Only metrics are read here (no images), but for UHN cases keep outDir on
% hospital storage like everything else. What leaves the hospital is your call
% under the REB and UHN policy.

study = cns_readJson(studyFile);
if ~exist(outDir, 'dir'), mkdir(outDir); end
caseDirs = cns_toCell(caseDirs);

R = {};
for i = 1:numel(caseDirs)
    f = fullfile(caseDirs{i}, 'metrics.csv');
    if ~exist(f, 'file'), warning('cns:summary', 'missing %s', f); continue; end
    R = [R; readMetrics(f)]; %#ok<AGROW>
end
if isempty(R), error('cns:summary', 'No metrics found.'); end
R = addSelected(R, study);
writeRows(fullfile(outDir, 'metrics_all.csv'), R);

refPlan = cns_getfield(study, 'summary.referencePlan', 'doseOnly');
refMach = cns_getfield(study, 'referenceMachine', '');
aliases = unique(R(:,1), 'stable');
keys = unique(strcat(R(:,4), '|', R(:,5), '|', R(:,6)), 'stable');
mp = unique(strcat(R(:,2), '|', R(:,3)), 'stable');

T = {};
for k = 1:numel(keys)
    kp = strsplit(keys{k}, '|');
    selK = strcmp(R(:,4), kp{1}) & strcmp(R(:,5), kp{2}) & strcmp(R(:,6), kp{3});
    unit = R{find(selK, 1), 8};
    for j = 1:numel(mp)
        mpj = strsplit(mp{j}, '|');
        x = vals(R, selK, mpj{1}, mpj{2}, aliases);
        if all(isnan(x)), continue; end
        xr = vals(R, selK, mpj{1}, refPlan, aliases);
        xm = vals(R, selK, refMach, mpj{2}, aliases);
        [dP, dPp, pP] = paired(x, xr);
        [dM, dMp, pM] = paired(x, xm);
        T(end+1, :) = {kp{1}, kp{2}, kp{3}, unit, mpj{1}, mpj{2}, nnz(~isnan(x)), ...
            med(x), min(x), max(x), dP, dPp, pP, dM, dMp, pM}; %#ok<AGROW>
    end
end

hdr = {'structure', 'quantity', 'metric', 'unit', 'machine', 'plan', 'n', 'median', 'min', 'max', ...
    'dVsRefPlan_median', 'dVsRefPlan_pct_median', 'p_vsRefPlan', ...
    'dVsRefMachine_median', 'dVsRefMachine_pct_median', 'p_vsRefMachine'};
fid = fopen(fullfile(outDir, 'summary_long.csv'), 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(T, 1)
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%d', T{i,1:7});
    fprintf(fid, ',%.5g', T{i,8:16});
    fprintf(fid, '\n');
end
fclose(fid);

writeAbstractTable(fullfile(outDir, 'summary_abstract.md'), T, study, refPlan, numel(aliases));
fprintf('[cns] cohort summary of %d cases in %s\n', numel(aliases), outDir);
end

% =========================================================================
function R = readMetrics(f)
fid = fopen(f, 'r');
fgetl(fid);
C = textscan(fid, '%s %s %s %s %s %s %f %s', 'Delimiter', ',');
fclose(fid);
R = [C{1:6}, num2cell(C{7}), C{8}];
end

function writeRows(f, R)
fid = fopen(f, 'w');
fprintf(fid, 'alias,machine,plan,structure,quantity,metric,value,unit\n');
for i = 1:size(R, 1)
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%.6g,%s\n', R{i,1:6}, R{i,7}, R{i,8});
end
fclose(fid);
end

function R = addSelected(R, study)
% Pick per case and machine the best LET plan that keeps coverage.
S = cns_getfield(study, 'selection', struct());
if isempty(fieldnames(S)), return; end
plans = cns_toCell(study.plans);
letIds = {};
for p = 1:numel(plans)
    if ~strcmp(cns_getfield(plans{p}, 'let.mode', 'none'), 'none'), letIds{end+1} = plans{p}.id; end %#ok<AGROW>
end
cov = S.coverage;  obj = S.objective;
aliases = unique(R(:,1), 'stable');
machines = unique(R(:,2), 'stable');
add = {};
for a = 1:numel(aliases)
    for m = 1:numel(machines)
        best = ''; bestV = Inf;
        refId = cns_getfield(study, 'summary.referencePlan', 'doseOnly');
        rx = get1(R, aliases{a}, machines{m}, refId, 'case', 'info', 'RxGyRBE');
        if isnan(rx), continue; end      % e.g. the clinical photon plan: nothing to select
        for p = 1:numel(letIds)
            ok = true;
            for c = 1:numel(cns_toCell(cov))
                C = cns_toCell(cov); C = C{c};
                v = get1(R, aliases{a}, machines{m}, letIds{p}, C.structure, C.quantity, C.metric);
                if isnan(v), ok = false; break; end
                if isfield(C, 'min') && v < C.min, ok = false; end
                if isfield(C, 'max') && v > C.max, ok = false; end
                if isfield(C, 'minFractionOfRx') && v < C.minFractionOfRx * rx, ok = false; end
                if isfield(C, 'maxFractionOfRx') && v > C.maxFractionOfRx * rx, ok = false; end
                ref = get1(R, aliases{a}, machines{m}, refId, C.structure, C.quantity, C.metric);
                if isfield(C, 'maxDecreasePctVsRef') && v < ref * (1 - C.maxDecreasePctVsRef / 100), ok = false; end
                if isfield(C, 'maxIncreasePctVsRef') && v > ref * (1 + C.maxIncreasePctVsRef / 100), ok = false; end
            end
            if ~ok, continue; end
            v = get1(R, aliases{a}, machines{m}, letIds{p}, obj.structure, obj.quantity, obj.metric);
            if ~isnan(v) && v < bestV, bestV = v; best = letIds{p}; end
        end
        if isempty(best)       % no LET plan keeps coverage: fall back to the reference
            best = refId;
            fprintf('[cns] %s / %s: no LET plan met the coverage rules\n', aliases{a}, machines{m});
        end
        sel = strcmp(R(:,1), aliases{a}) & strcmp(R(:,2), machines{m}) & strcmp(R(:,3), best);
        rows = R(sel, :);
        rows(:,3) = {'LETselected'};
        add = [add; rows; {aliases{a}, machines{m}, 'LETselected', 'case', 'info', ['selected_' best], 1, '-'}]; %#ok<AGROW>
        fprintf('[cns] %s / %s: LETselected = %s\n', aliases{a}, machines{m}, best);
    end
end
R = [R; add];
end

function v = get1(R, a, m, p, s, q, met)
sel = strcmp(R(:,1), a) & strcmp(R(:,2), m) & strcmp(R(:,3), p) & strcmp(R(:,4), s) ...
    & strcmp(R(:,5), q) & strcmp(R(:,6), met);
if any(sel), v = R{find(sel, 1), 7}; else, v = NaN; end
end

function x = vals(R, selK, machine, plan, aliases)
x = nan(numel(aliases), 1);
if isempty(machine) || isempty(plan), return; end
sel = selK & strcmp(R(:,2), machine) & strcmp(R(:,3), plan);
idx = find(sel);
for i = idx'
    a = find(strcmp(aliases, R{i,1}), 1);
    x(a) = R{i,7};
end
end

function [d, dp, p] = paired(x, r)
ok = ~isnan(x) & ~isnan(r);
d = NaN; dp = NaN; p = NaN;
if ~any(ok), return; end
d = med(x(ok) - r(ok));
rel = 100 * (x(ok) - r(ok)) ./ r(ok);
rel = rel(isfinite(rel));
if ~isempty(rel), dp = med(rel); end
if nnz(ok) >= 5 && exist('signrank', 'file') == 2 && any(x(ok) ~= r(ok))
    try, p = signrank(x(ok), r(ok)); catch, end %#ok<NOCOM>
end
end

function m = med(x)
x = x(~isnan(x));
if isempty(x), m = NaN; else, m = median(x); end
end

function writeAbstractTable(f, T, study, refPlan, nCases)
K = cns_toCell(cns_getfield(study, 'summary.keyMetrics', {}));
machines = unique(T(:,5), 'stable');
cmpPlan = cns_getfield(study, 'summary.comparePlan', 'LETselected');
fid = fopen(f, 'w');
fprintf(fid, '# Cohort summary (%d case%s)\n\n', nCases, ifelse(nCases > 1, 's', ''));
fprintf(fid, ['Median [min, max] over cases. Delta = paired median difference of %s vs %s, same machine. ' ...
    'Clinical (photon) plan: RBE = 1, so RBE11 and McN rows are its physical dose and its LET rows are zero.\n\n'], cmpPlan, refPlan);
refLabel = cns_getfield(study, 'summary.clinicalLabel', 'photon');
hasClin = any(strcmp(T(:,5), refLabel) & strcmp(T(:,6), 'clinical'));
machines = machines(~strcmp(machines, refLabel));
hdr = '| Structure | Quantity | Metric | Unit |';
sep = '|---|---|---|---|';
if hasClin
    hdr = [hdr sprintf(' %s clinical |', refLabel)]; sep = [sep '---|'];
end
for m = 1:numel(machines)
    hdr = [hdr sprintf(' %s %s | %s %s | Delta (%%) |', machines{m}, refPlan, machines{m}, cmpPlan)]; %#ok<AGROW>
    sep = [sep '---|---|---|']; %#ok<AGROW>
end
fprintf(fid, '%s\n%s\n', hdr, sep);
for k = 1:numel(K)
    s = K{k};
    sel = strcmp(T(:,1), s.structure) & strcmp(T(:,2), s.quantity) & strcmp(T(:,3), s.metric);
    if ~any(sel), continue; end
    line = sprintf('| %s | %s | %s | %s |', s.structure, s.quantity, strrep(s.metric, '_', ' '), T{find(sel, 1), 4});
    if hasClin
        line = [line, ' ', cellTxt(T, sel & strcmp(T(:,5), refLabel) & strcmp(T(:,6), 'clinical')), ' |']; %#ok<AGROW>
    end
    for m = 1:numel(machines)
        a = sel & strcmp(T(:,5), machines{m}) & strcmp(T(:,6), refPlan);
        b = sel & strcmp(T(:,5), machines{m}) & strcmp(T(:,6), cmpPlan);
        line = [line, ' ', cellTxt(T, a), ' | ', cellTxt(T, b), ' | ', ...
            deltaTxt(T, b), ' |']; %#ok<AGROW>
    end
    fprintf(fid, '%s\n', line);
end

% machine comparison, same plan
refMach = cns_getfield(study, 'referenceMachine', '');
others = machines(~strcmp(machines, refMach) & ~strcmp(machines, refLabel));
if ~isempty(refMach) && ~isempty(others)
    fprintf(fid, '\n## Machine comparison vs %s (same plan)\n\n', refMach);
    hdr = '| Structure | Quantity | Metric | Unit |'; sep = '|---|---|---|---|';
    for m = 1:numel(others)
        for pl = {refPlan, cmpPlan}
            hdr = [hdr sprintf(' %s %s: delta vs %s |', others{m}, pl{1}, refMach)]; %#ok<AGROW>
            sep = [sep '---|']; %#ok<AGROW>
        end
    end
    fprintf(fid, '%s\n%s\n', hdr, sep);
    for k = 1:numel(K)
        s = K{k};
        sel = strcmp(T(:,1), s.structure) & strcmp(T(:,2), s.quantity) & strcmp(T(:,3), s.metric);
        if ~any(sel), continue; end
        line = sprintf('| %s | %s | %s | %s |', s.structure, s.quantity, strrep(s.metric, '_', ' '), T{find(sel, 1), 4});
        for m = 1:numel(others)
            for pl = {refPlan, cmpPlan}
                b = sel & strcmp(T(:,5), others{m}) & strcmp(T(:,6), pl{1});
                line = [line, ' ', machineDeltaTxt(T, b), ' |']; %#ok<AGROW>
            end
        end
        fprintf(fid, '%s\n', line);
    end
end

% which LET plan was selected
sel = strcmp(T(:,1), 'case') & strcmp(T(:,2), 'info') & strncmp(T(:,3), 'selected_', 9) & strcmp(T(:,6), 'LETselected');
if any(sel)
    fprintf(fid, '\n## LETselected\n\n| Machine | Plan | Cases |\n|---|---|---|\n');
    for i = find(sel)'
        fprintf(fid, '| %s | %s | %d |\n', T{i,5}, strrep(T{i,3}, 'selected_', ''), T{i,7});
    end
end
fclose(fid);
end

function s = machineDeltaTxt(T, sel)
if ~any(sel), s = 'n/a'; return; end
i = find(sel, 1);
if isnan(T{i,15}), s = sprintf('%+.3g', T{i,14}); else, s = sprintf('%+.3g (%+.1f %%)', T{i,14}, T{i,15}); end
if ~isnan(T{i,16}), s = sprintf('%s, p=%.2g', s, T{i,16}); end
end

function s = cellTxt(T, sel)
if ~any(sel), s = 'n/a'; return; end
i = find(sel, 1);
if T{i,7} > 1
    s = sprintf('%.3g [%.3g, %.3g]', T{i,8}, T{i,9}, T{i,10});
else
    s = sprintf('%.3g', T{i,8});
end
end

function s = deltaTxt(T, sel)
if ~any(sel), s = 'n/a'; return; end
i = find(sel, 1);
if isnan(T{i,12}), s = sprintf('%+.3g', T{i,11}); else, s = sprintf('%+.3g (%+.1f %%)', T{i,11}, T{i,12}); end
if ~isnan(T{i,13}), s = sprintf('%s, p=%.2g', s, T{i,13}); end
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end
