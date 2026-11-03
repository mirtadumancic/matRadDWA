function m = letAgent_metrics(dij, cst, w, caseCfg, pln)
% letAgent_metrics  Plan metrics for the agent: numbers only, no images.
%
% Doses are total doses in Gy (fraction dose x numOfFractions). LETd is the
% dose-averaged LET in keV/um, evaluated where the total dose exceeds
% metrics.letdMinDoseGy. Dirty dose uses the threshold the dij was built
% with. Every metric is given for the nominal scenario and as the worst
% case over all robustness scenarios in the dij.

nFx = pln.numOfFractions;
cfgM = letAgent_getfield(caseCfg, 'metrics', struct());
rx = caseCfg.prescription.totalDoseGy;
letdMinDose = letAgent_getfield(cfgM, 'letdMinDoseGy', 0.02 * rx);
vdlDose = letAgent_getfield(cfgM, 'letVolume.minDoseGy', 0.5 * rx);
vdlLet  = letAgent_getfield(cfgM, 'letVolume.minLETd', 2.5);

cstD = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
    dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

scens = find(~cellfun(@isempty, dij.physicalDose(:)))';
nS = numel(scens);
d = cell(1,nS); letd = cell(1,nS); dd = cell(1,nS);
for k = 1:nS
    s = scens(k);
    dFx = full(dij.physicalDose{s} * w);
    d{k} = dFx * nFx;
    letd{k} = zeros(size(dFx));
    if isfield(dij, 'mLETDose') && ~isempty(dij.mLETDose{s})
        L = full(dij.mLETDose{s} * w);
        ok = d{k} > letdMinDose;
        letd{k}(ok) = L(ok) ./ dFx(ok);
    end
    dd{k} = zeros(size(dFx));
    if isfield(dij, 'mDirtyDose') && ~isempty(dij.mDirtyDose{s})
        dd{k} = full(dij.mDirtyDose{s} * w) * nFx;
    end
end

m = struct();
m.numOfScenarios = nS;
m.prescriptionGy = rx;
m.dirtyDoseLETThreshold = letAgent_getfield(dij, 'mDirtyDoseThreshold', NaN);

targets = letAgent_toCell(caseCfg.structures.targets);
m.targets = {};
for i = 1:numel(targets)
    t = targets{i};
    v = cstD{letAgent_findStructure(cstD, t.name), 4}{1};
    D95 = zeros(1,nS); D98 = D95; D2 = D95; Dmean = D95; V95 = D95;
    for k = 1:nS
        x = d{k}(v);
        D95(k) = dx(x, 95); D98(k) = dx(x, 98); D2(k) = dx(x, 2);
        Dmean(k) = mean(x); V95(k) = 100 * mean(x >= 0.95 * t.doseGy);
    end
    r = struct('name', t.name, 'doseGy', t.doseGy);
    r.D95 = nomwc(D95, 'min'); r.D98 = nomwc(D98, 'min'); r.D2 = nomwc(D2, 'max');
    r.Dmean = nomwc(Dmean, 'min'); r.V95pct = nomwc(V95, 'min');
    r.LETdMean = doseWeightedLETd(d{1}(v), letd{1}(v));
    r.dirtyDoseMean = mean(dd{1}(v));
    m.targets{end+1} = r;
end

oars = letAgent_toCell(letAgent_getfield(caseCfg, 'structures.oars', {}));
m.oars = {};
for i = 1:numel(oars)
    o = oars{i};
    v = cstD{letAgent_findStructure(cstD, o.name), 4}{1};
    D2 = zeros(1,nS); Dmean = D2; Dmax = D2; LETdP98 = D2; VDL = D2; DD2 = D2; DDmean = D2;
    for k = 1:nS
        x = d{k}(v); l = letd{k}(v); y = dd{k}(v);
        D2(k) = dx(x, 2); Dmean(k) = mean(x); Dmax(k) = max(x);
        hot = x > letdMinDose;
        if any(hot), LETdP98(k) = dx(l(hot), 2); end
        VDL(k) = 100 * mean(x >= vdlDose & l >= vdlLet);
        DD2(k) = dx(y, 2); DDmean(k) = mean(y);
    end
    r = struct('name', o.name, 'maxDoseGy', letAgent_getfield(o, 'maxDoseGy', NaN), ...
        'let', logical(letAgent_getfield(o, 'let', false)));
    r.D2 = nomwc(D2, 'max'); r.Dmax = nomwc(Dmax, 'max'); r.Dmean = nomwc(Dmean, 'max');
    r.LETdP98 = nomwc(LETdP98, 'max');
    r.LETdMean = doseWeightedLETd(d{1}(v), letd{1}(v));
    r.VDLpct = nomwc(VDL, 'max');
    r.dirtyDoseD2 = nomwc(DD2, 'max'); r.dirtyDoseMean = nomwc(DDmean, 'max');
    m.oars{end+1} = r;
end
m.letVolumeDefinition = struct('minDoseGy', vdlDose, 'minLETd', vdlLet);
end

function v = dx(x, pct)
% Dose to the hottest pct % of the volume (D2 = near max, D95 = coverage).
if isempty(x), v = NaN; return; end
x = sort(x(:), 'descend');
v = x(max(1, min(numel(x), ceil(pct / 100 * numel(x)))));
end

function s = nomwc(vals, worst)
s.nom = vals(1);
if strcmp(worst, 'min'), s.wc = min(vals); else, s.wc = max(vals); end
end

function l = doseWeightedLETd(dose, letd)
if sum(dose) <= 0, l = 0; else, l = sum(dose .* letd) / sum(dose); end
end
