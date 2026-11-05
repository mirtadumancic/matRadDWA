function [cst, notes] = cns_buildCst(cst, caseCfg, study, planDef, pln)
% cns_buildCst  Objectives for one plan variant.
%
% Dose objectives come from study.objectives (one list per canonical
% structure). Doses are total Gy(RBE) with RBE 1.1: either "doseGy"
% (absolute) or "doseFraction" (x prescription). A case file can override
% or add structures with caseCfg.objectives (same format).
%
% planDef.let.mode   'none' | 'LETxDose' | 'dirtyDose'
% planDef.let.penalty          weight of SquaredOverdosing(., 0) on the quantity
% planDef.let.structures       default: study.letStructures
% planDef.targetPenaltyScale   multiplies all PTV/CTV penalties (default 1)
%
% 'LETxDose' uses Remo's 'LETd' quantity, which is LET x dose
% (dij.mLETDose * w). 'dirtyDose' needs matRad_DirtyDose (letAgent).

notes = {};
rx = caseCfg.prescription.totalDoseGyRBE;
for i = 1:size(cst,1), cst{i,6} = {}; end

objs = cns_getfield(study, 'objectives', struct());
caseObjs = cns_getfield(caseCfg, 'objectives', struct());
fn = fieldnames(caseObjs);
for i = 1:numel(fn), objs.(fn{i}) = caseObjs.(fn{i}); end

rob = cns_getfield(study, 'robustness.objective', 'none');
if pln.multScen.totNumScen < 2, rob = 'none'; end
tScale = cns_getfield(planDef, 'targetPenaltyScale', 1);

names = fieldnames(objs);
for i = 1:numel(names)
    ix = find(strcmp(cst(:,2), names{i}), 1);
    if isempty(ix)
        notes{end+1} = sprintf('no structure %s: its objectives are skipped', names{i}); %#ok<AGROW>
        continue;
    end
    list = cns_toCell(objs.(names{i}));
    isTarget = any(strcmp(names{i}, {'PTV', 'CTV'}));
    for k = 1:numel(list)
        o = list{k};
        d = cns_getfield(o, 'doseGy', []);
        if isempty(d), d = cns_getfield(o, 'doseFraction', 1) * rx; end
        p = o.penalty;
        if isTarget, p = p * tScale; end
        obj = makeObjective(o.type, p, d, cns_getfield(o, 'volumePct', []));
        if ~strcmp(rob, 'none'), obj.robustness = rob; end
        cst{ix,6}{end+1} = obj;
    end
end

mode = cns_getfield(planDef, 'let.mode', 'none');
switch mode
    case 'none'
        return;
    case 'LETxDose'
        quantity = 'LETd';
    case 'dirtyDose'
        quantity = 'dirtyDose';
    otherwise
        error('cns:plan', 'Unknown let.mode %s', mode);
end
letStructs = cns_toCell(cns_getfield(planDef, 'let.structures', cns_getfield(study, 'letStructures', {})));
p = planDef.let.penalty;
ref = cns_getfield(planDef, 'let.ref', 0);
nUsed = 0;
for i = 1:numel(letStructs)
    ix = find(strcmp(cst(:,2), letStructs{i}), 1);
    if isempty(ix), continue; end
    obj = DoseObjectives.matRad_SquaredOverdosing(p, ref);
    obj.quantity = quantity;
    if ~strcmp(rob, 'none'), obj.robustness = rob; end
    cst{ix,6}{end+1} = obj;
    nUsed = nUsed + 1;
end
if nUsed == 0
    notes{end+1} = 'LET objective requested but none of letStructures exists';
end
end

function obj = makeObjective(type, p, d, vol)
switch type
    case 'SquaredDeviation',   obj = DoseObjectives.matRad_SquaredDeviation(p, d);
    case 'SquaredUnderdosing', obj = DoseObjectives.matRad_SquaredUnderdosing(p, d);
    case 'SquaredOverdosing',  obj = DoseObjectives.matRad_SquaredOverdosing(p, d);
    case 'MeanDose',           obj = DoseObjectives.matRad_MeanDose(p, d);
    case 'MinDVH',             obj = DoseObjectives.matRad_MinDVH(p, d, vol);
    case 'MaxDVH',             obj = DoseObjectives.matRad_MaxDVH(p, d, vol);
    otherwise
        error('cns:objective', 'Objective type %s not supported', type);
end
end
