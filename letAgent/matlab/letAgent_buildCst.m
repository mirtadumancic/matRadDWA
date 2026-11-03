function [cst, notes] = letAgent_buildCst(cst, caseCfg, recipe, pln)
% letAgent_buildCst  Turn a search "recipe" into matRad objectives/constraints.
%
% The case file holds the clinical template (prescription, OAR limits,
% base penalties). The recipe holds what the search varies:
%   recipe.targetPenaltyScale   multiplies every target penalty (default 1)
%   recipe.oarDosePenaltyScale  multiplies every OAR dose penalty (default 1)
%   recipe.let.mode             'none' | 'LETxDose' | 'dirtyDose'
%   recipe.let.penalty          penalty of the LET objective
%   recipe.let.structures       OARs that get it (default: OARs with let=true)
%   recipe.let.refGy            reference of the squared overdosing (default 0)
%   recipe.robust               'none' | 'COWC' | 'STOCH' | 'VWWC' (objectives)
%   recipe.coverage.mode        'penalty' (default) | 'constraint'
%   recipe.coverage.minFraction dose level as fraction of Rx (default 0.95)
%   recipe.coverage.minVolumePct volume that must reach it (default 95)

notes = {};
for i = 1:size(cst,1)
    cst{i,6} = {};
end

rob = letAgent_getfield(recipe, 'robust', 'none');
if ~strcmp(rob, 'none') && pln.multScen.totNumScen < 2
    notes{end+1} = sprintf('robust=%s requested but the case has only the nominal scenario; optimised nominally', rob);
    rob = 'none';
end
tScale = letAgent_getfield(recipe, 'targetPenaltyScale', 1);
oScale = letAgent_getfield(recipe, 'oarDosePenaltyScale', 1);

% --- targets
targets = letAgent_toCell(caseCfg.structures.targets);
for i = 1:numel(targets)
    t = targets{i};
    ix = letAgent_findStructure(cst, t.name);
    obj = DoseObjectives.matRad_SquaredDeviation(letAgent_getfield(t, 'penalty', 800) * tScale, t.doseGy);
    obj.quantity = 'physicalDose';
    obj.robustness = rob;
    cst{ix,6}{end+1} = obj;

    if strcmp(letAgent_getfield(recipe, 'coverage.mode', 'penalty'), 'constraint')
        f = letAgent_getfield(recipe, 'coverage.minFraction', 0.95);
        v = letAgent_getfield(recipe, 'coverage.minVolumePct', 95);
        con = DoseConstraints.matRad_MinMaxDVH(f * t.doseGy, v, 100);
        con.quantity = 'physicalDose';
        cst{ix,6}{end+1} = con;
    end
end

% --- OAR dose objectives
oars = letAgent_toCell(letAgent_getfield(caseCfg, 'structures.oars', {}));
letDefault = {};
for i = 1:numel(oars)
    o = oars{i};
    ix = letAgent_findStructure(cst, o.name);
    if isfield(o, 'maxDoseGy') && ~isempty(o.maxDoseGy)
        obj = DoseObjectives.matRad_SquaredOverdosing(letAgent_getfield(o, 'penalty', 100) * oScale, o.maxDoseGy);
        obj.quantity = 'physicalDose';
        obj.robustness = rob;
        cst{ix,6}{end+1} = obj;
    end
    if isfield(o, 'meanDosePenalty') && ~isempty(o.meanDosePenalty)
        obj = DoseObjectives.matRad_MeanDose(o.meanDosePenalty, 0);
        obj.quantity = 'physicalDose';
        cst{ix,6}{end+1} = obj;
    end
    if letAgent_getfield(o, 'let', false)
        letDefault{end+1} = o.name; %#ok<AGROW>
    end
end

% --- LET-modifying objective
mode = letAgent_getfield(recipe, 'let.mode', 'none');
switch mode
    case 'none'
        quantity = '';
    case 'LETxDose'
        quantity = 'LETd';   % matRad_LETd optimises dij.mLETDose * w, i.e. LET x dose
    case 'dirtyDose'
        quantity = 'dirtyDose';
    otherwise
        error('letAgent:recipe', 'Unknown let.mode %s', mode);
end
if ~isempty(quantity)
    letStructs = letAgent_toCell(letAgent_getfield(recipe, 'let.structures', letDefault));
    if isempty(letStructs)
        notes{end+1} = 'LET objective requested but no OAR has let=true in the case file';
    end
    p = recipe.let.penalty;
    ref = letAgent_getfield(recipe, 'let.refGy', 0);
    for i = 1:numel(letStructs)
        ix = letAgent_findStructure(cst, letStructs{i});
        obj = DoseObjectives.matRad_SquaredOverdosing(p, ref);
        obj.quantity = quantity;
        obj.robustness = rob;
        cst{ix,6}{end+1} = obj;
    end
end
end
