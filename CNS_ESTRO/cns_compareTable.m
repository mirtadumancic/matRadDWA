function T = cns_compareTable(rows, keyMetrics, columns, base, constraints)
% cns_compareTable  One row per key metric, one column per plan; csv + markdown.
%
%   T = cns_compareTable(rows, study.summary.keyMetrics, columns, 'D:/.../table_doseOnly')
%   T = cns_compareTable(rows, keyMetrics, columns, base, study.constraintsCheck.list)
%
% rows     long metric rows {alias, machine, plan, structure, quantity, metric, value, unit}
% columns  struct array: machine, plan, label
% A single case gives values only (no statistics); for cohorts use cns_cohortSummary.
% With constraints, a second table lists pass/fail per plan.
K = cns_toCell(keyMetrics);
hdr = [{'structure', 'quantity', 'metric', 'unit'}, {columns.label}];
T = cell(0, numel(hdr));
for k = 1:numel(K)
    s = K{k};
    r = {s.structure, s.quantity, s.metric, unitOf(rows, s)};
    for c = 1:numel(columns)
        r{end+1} = get1(rows, columns(c).machine, columns(c).plan, s.structure, s.quantity, s.metric); %#ok<AGROW>
    end
    T(end+1, :) = r; %#ok<AGROW>
end
writeCsv([base '.csv'], hdr, T);

fid = fopen([base '.md'], 'w');
fprintf(fid, '| %s |\n|%s\n', strjoin(hdr, ' | '), repmat('---|', 1, numel(hdr)));
for i = 1:size(T, 1)
    v = cellfun(@fmt, T(i, 5:end), 'UniformOutput', false);
    fprintf(fid, '| %s | %s | %s | %s | %s |\n', T{i, 1}, T{i, 2}, strrep(T{i, 3}, '_', ' '), T{i, 4}, strjoin(v, ' | '));
end
if nargin > 4 && ~isempty(constraints)
    C = cns_toCell(constraints);
    fprintf(fid, '\n**Constraint check** (limits from study.constraintsCheck: replace with the UHN protocol)\n\n');
    fprintf(fid, '| Structure | Metric | Limit | %s |\n|%s\n', strjoin({columns.label}, ' | '), repmat('---|', 1, 3 + numel(columns)));
    Crows = {};
    for k = 1:numel(C)
        c = C{k};
        line = sprintf('| %s | %s %s | <= %g |', c.structure, c.quantity, strrep(c.metric, '_', ' '), c.max);
        r = {c.structure, c.quantity, c.metric, c.max};
        for j = 1:numel(columns)
            v = get1(rows, columns(j).machine, columns(j).plan, c.structure, c.quantity, c.metric);
            if isnan(v), t = 'n/a'; elseif v <= c.max, t = sprintf('%.1f ok', v); else, t = sprintf('**%.1f FAIL**', v); end
            line = [line ' ' t ' |']; %#ok<AGROW>
            r{end+1} = v; %#ok<AGROW>
        end
        fprintf(fid, '%s\n', line);
        Crows(end+1, :) = r; %#ok<AGROW>
    end
    writeCsv([base '_constraints.csv'], [{'structure', 'quantity', 'metric', 'max'}, {columns.label}], Crows);
end
fclose(fid);
end

function v = get1(rows, m, p, s, q, met)
sel = strcmp(rows(:,2), m) & strcmp(rows(:,3), p) & strcmp(rows(:,4), s) & strcmp(rows(:,5), q) & strcmp(rows(:,6), met);
if any(sel), v = rows{find(sel, 1), 7}; else, v = NaN; end
end

function u = unitOf(rows, s)
sel = strcmp(rows(:,4), s.structure) & strcmp(rows(:,5), s.quantity) & strcmp(rows(:,6), s.metric);
if any(sel), u = rows{find(sel, 1), 8}; else, u = ''; end
end

function t = fmt(v)
if ischar(v), t = v; elseif isnan(v), t = 'n/a'; elseif abs(v) >= 100, t = sprintf('%.0f', v);
elseif abs(v) >= 10, t = sprintf('%.1f', v); else, t = sprintf('%.2f', v); end
end

function writeCsv(f, hdr, T)
fid = fopen(f, 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(T, 1)
    parts = cellfun(@(v) ifc(ischar(v), v, sprintf('%.6g', v)), T(i, :), 'UniformOutput', false);
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end

function o = ifc(c, a, b)
if c, o = a; else, o = b; end
end
