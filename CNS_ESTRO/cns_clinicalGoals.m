function T = cns_clinicalGoals(dose, cstG, voxCc, goals, label)
% cns_clinicalGoals  Evaluate a list of clinical goals on one dose cube.
%
%   T = cns_clinicalGoals(dose, cstG, voxCc, goals, label)
%
% dose   total dose cube, Gy(RBE) (any grid)
% cstG   cst on the same grid as dose (canonical structure names)
% voxCc  voxel volume of that grid, cm^3
% goals  study.clinicalGoals.list: structure, type, param, comparator, limit, unit, label
%          type 'Vgy'   volume (%) receiving >= param Gy
%               'Vgycc' volume (cc) receiving >= param Gy
%               'Dpct'  dose (Gy) to the hottest param % of the volume
%               'Dcc'   dose (Gy) to the hottest param cc
%               'mean'  mean dose (Gy)
% label  text for the 'evaluation' column (e.g. 'PB nominal')
%
% T: table with goal, structure, value, limit, unit, pass (NaN if the structure is missing).
% No patient identifiers are written: only canonical structure names and numbers.
goals = cns_toCell(goals);
n = numel(goals);
goal = cell(n, 1); st = cell(n, 1); val = nan(n, 1); lim = nan(n, 1); unit = cell(n, 1);
pass = nan(n, 1); nVox = zeros(n, 1);
for k = 1:n
    g = goals{k};
    goal{k} = g.label; st{k} = g.structure; lim(k) = g.limit; unit{k} = g.unit;
    ix = find(strcmp(cstG(:,2), g.structure), 1);
    if isempty(ix) || isempty(cstG{ix,4}{1}), continue; end
    d = sort(dose(cstG{ix,4}{1}), 'descend');
    d = d(:); nVox(k) = numel(d);
    switch g.type
        case 'Vgy',   v = 100 * mean(d >= g.param);
        case 'Vgycc', v = nnz(d >= g.param) * voxCc;
        case 'Dpct',  v = d(max(1, min(numel(d), ceil(g.param / 100 * numel(d)))));
        case 'Dcc',   v = d(max(1, min(numel(d), ceil(g.param / voxCc))));
        case 'mean',  v = mean(d);
        otherwise, error('cns:goal', 'Unknown goal type %s', g.type);
    end
    val(k) = v;
    switch g.comparator
        case '>=', pass(k) = v >= g.limit - 1e-9;
        case '<=', pass(k) = v <= g.limit + 1e-9;
    end
end
evaluation = repmat({label}, n, 1);
T = table(evaluation, goal, st, val, lim, unit, pass, nVox, ...
    'VariableNames', {'evaluation', 'goal', 'structure', 'value', 'limit', 'unit', 'pass', 'nVoxels'});
end
