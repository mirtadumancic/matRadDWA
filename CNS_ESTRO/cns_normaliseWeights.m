function [w, f, goal] = cns_normaliseWeights(w, dij, cstD, caseCfg, study, refDose, cstCt)
% cns_normaliseWeights  Scale spot weights so <structure> <Dxx> (RBE 1.1, total) hits a goal.
%
%   [w, f, goal] = cns_normaliseWeights(w, dij, cstD, caseCfg, study)
%   [w, f, goal] = cns_normaliseWeights(w, dij, cstD, caseCfg, study, refDose, cstCt)
%
% study.normalisation:
%   structure     'PTV' (or 'CTV'; 'none' switches normalisation off)
%   metric        'D95' (any Dxx)
%   mode          'fractionOfRx' : goal = fractionOfRx * prescription      (default)
%                 'matchReference': goal = the same Dxx of the clinical photon
%                                   plan (refDose), so protons and photons
%                                   are compared at equal target coverage
%   fractionOfRx  e.g. 0.95
% cstD is the cst on the dose grid; cstCt (CT grid) is needed for matchReference.
f = 1; goal = NaN;
N = cns_getfield(study, 'normalisation', struct());
s = cns_getfield(N, 'structure', 'none');
if strcmp(s, 'none'), return; end
ix = find(strcmp(cstD(:,2), s), 1);
if isempty(ix), return; end
pct = str2double(regexprep(cns_getfield(N, 'metric', 'D95'), '^D', ''));
rx = caseCfg.prescription.totalDoseGyRBE;
nFx = caseCfg.prescription.numOfFractions;

mode = cns_getfield(N, 'mode', 'fractionOfRx');
if strcmp(mode, 'matchReference') && nargin >= 7 && ~isempty(refDose)
    t = refDose.toTotal; if isnan(t), t = nFx; end
    D = double(refDose.cube) * t;
    j = find(strcmp(cstCt(:,2), s), 1);
    goal = dx(D(cstCt{j,4}{1}), pct);
else
    if strcmp(mode, 'matchReference')
        warning('cns:norm', 'normalisation.mode = matchReference but no reference dose: using fractionOfRx.');
    end
    goal = cns_getfield(N, 'fractionOfRx', 0.95) * rx;
end
d = full(dij.physicalDose{1} * w) * nFx * 1.1;
cur = dx(d(cstD{ix,4}{1}), pct);
if cur > 0
    f = goal / cur;
    w = w * f;
end
end

function v = dx(x, pct)
x = sort(x(:), 'descend');
v = x(max(1, min(numel(x), ceil(pct / 100 * numel(x)))));
end
