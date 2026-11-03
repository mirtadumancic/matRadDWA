function multScen = letAgent_buildMultScen(ct, caseCfg)
% letAgent_buildMultScen  Robustness scenario model from the case file.
% Rebuilt on load because Octave stores classdef objects as plain structs.
rob = letAgent_getfield(caseCfg, 'robustness', struct('scenarios', 'nomScen'));
multScen = matRad_multScen(ct, rob.scenarios);
if ~strcmp(rob.scenarios, 'nomScen')
    if isfield(rob, 'shiftSDmm'),     multScen.shiftSD = rob.shiftSDmm(:)'; end
    if isfield(rob, 'rangeRelSDpct'), multScen.rangeRelSD = rob.rangeRelSDpct; end
    if isfield(rob, 'rangeAbsSDmm'),  multScen.rangeAbsSD = rob.rangeAbsSDmm; end
end
end
