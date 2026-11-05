function E = cns_evaluatePlan(dij, cst, w, caseCfg, study, opts)
% cns_evaluatePlan  Target DVH and OAR dose-LET metrics for one plan.
%
%   E = cns_evaluatePlan(dij, cst, w, caseCfg, study)
%
% All doses are totals (fraction dose x numOfFractions) on the dose grid.
% Quantities:
%   phys    physical dose [Gy]
%   RBE11   RBE-weighted dose with RBE = 1.1 [Gy(RBE)]   (the clinical dose)
%   McN     McNamara variable-RBE dose [Gy(RBE)], alpha/beta 10 Gy in
%           CTV/PTV/GTV and 2 Gy elsewhere (study.evaluation.alphaBeta*)
%   LETd    dose-averaged LET [keV/um], only where RBE11 >= letdMinDoseGy
%   DxLETd  RBE11 x LETd smoothed with a Gaussian (Heinzelmann 2026), [Gy keV/um]
%   dirty   dirty dose [Gy], if the dij carries mDirtyDose
%
% E.rows   {structure, quantity, metric, value, unit}
% E.dvh    struct array: structure, quantity, dose (bin edges), volPct
% E.lvh    struct array: structure, letd (bin edges), volPct, volCc
% E.cubes  RBE11, McN, LETd (dose grid), kept for slice figures
%
% opts.photon = true evaluates a photon dose (cns_evaluateReferenceDose):
% RBE = 1, so phys = RBE11 = McN, and LETd and the LET metrics are zero.
if nargin < 6, opts = struct(); end
isPhoton = logical(cns_getfield(opts, 'photon', false));

ev = cns_getfield(study, 'evaluation', struct());
nFx = caseCfg.prescription.numOfFractions;
rx = caseCfg.prescription.totalDoseGyRBE;
abT = cns_getfield(ev, 'alphaBetaTarget', 10);
abN = cns_getfield(ev, 'alphaBetaNormal', 2);
letdMin = cns_getfield(ev, 'letdMinDoseGy', 1);
dvhBin = cns_getfield(ev, 'dvhBinGy', 0.5);
lvhBin = cns_getfield(ev, 'lvhBinKeVum', 0.1);
sigmaDL = cns_getfield(ev, 'dLetd.sigmaMm', 6);
minDL = cns_getfield(ev, 'dLetd.minDoseGy', 5);

dims = dij.doseGrid.dimensions;
res = dij.doseGrid.resolution;
vox = res.x * res.y * res.z / 1000;      % cc
if isPhoton
    cstD = cst;                          % reference dose is already on the CT grid
else
    cstD = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
        dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
end

[phys, rbe11, letd, mcn, dirty] = cubesForScenario(dij, 1, w, nFx, cstD, abT, abN, letdMin, isPhoton);
dl = cns_gaussSmooth3(reshape(rbe11 .* letd, dims), sigmaDL, res);
dl = dl(:);

targets = present(cstD, cns_toCell(cns_getfield(ev, 'targets', {'CTV', 'PTV', 'GTV'})));
oars = present(cstD, cns_toCell(cns_getfield(ev, 'oars', {})));
hot = cns_toCell(cns_getfield(ev, 'hotSpots', {}));

rows = {'case', 'info', 'RxGyRBE', rx, 'Gy(RBE)'; 'case', 'info', 'numOfFractions', nFx, '-'};
dvh = struct('structure', {}, 'quantity', {}, 'dose', {}, 'volPct', {});
lvh = struct('structure', {}, 'letd', {}, 'volPct', {}, 'volCc', {});
dBins = 0:dvhBin:ceil(1.3 * rx);
lBins = 0:lvhBin:15;

% ---- targets
for i = 1:numel(targets)
    v = idx(cstD, targets{i});
    rows(end+1, :) = {targets{i}, 'geom', 'volumeCc', numel(v) * vox, 'cc'}; %#ok<AGROW>
    for q = {'phys', 'RBE11', 'McN'}
        x = pick(q{1}, phys, rbe11, mcn, v);
        u = unitOf(q{1});
        rows = [rows; {targets{i}, q{1}, 'D98', dx(x, 98), u; targets{i}, q{1}, 'D95', dx(x, 95), u; ...
            targets{i}, q{1}, 'D50', dx(x, 50), u; targets{i}, q{1}, 'D2', dx(x, 2), u; ...
            targets{i}, q{1}, 'Dmean', mean(x), u; ...
            targets{i}, q{1}, 'V95pct', 100 * mean(x >= 0.95 * rx), '%'; ...
            targets{i}, q{1}, 'V107pct', 100 * mean(x >= 1.07 * rx), '%'}]; %#ok<AGROW>
        if ~strcmp(q{1}, 'phys')
            dvh(end+1) = struct('structure', targets{i}, 'quantity', q{1}, 'dose', dBins, ...
                'volPct', cumVol(x, dBins)); %#ok<AGROW>
        end
    end
    rows(end+1, :) = {targets{i}, 'LETd', 'LETdMean', dwMean(phys(v), letd(v)), 'keV/um'}; %#ok<AGROW>
end

% ---- robustness of the CTV (RBE 1.1), if the dij has error scenarios
scens = find(~cellfun(@isempty, dij.physicalDose(:)))';
if numel(scens) > 1 && any(strcmp(targets, 'CTV'))
    v = idx(cstD, 'CTV');
    D95 = zeros(1, numel(scens)); D98 = D95;
    for k = 1:numel(scens)
        dS = full(dij.physicalDose{scens(k)} * w) * nFx * 1.1;
        D95(k) = dx(dS(v), 95); D98(k) = dx(dS(v), 98);
    end
    rows = [rows; {'CTV', 'RBE11', 'D95_worst', min(D95), 'Gy(RBE)'; ...
        'CTV', 'RBE11', 'D98_worst', min(D98), 'Gy(RBE)'}];
end

% ---- OARs and helper structures
for i = 1:numel(oars)
    s = oars{i};
    v = idx(cstD, s);
    rows(end+1, :) = {s, 'geom', 'volumeCc', numel(v) * vox, 'cc'}; %#ok<AGROW>
    n03 = max(1, round(0.03 / vox));     % voxels in 0.03 cc
    if isempty(v)
        rows(end+1, :) = {s, 'geom', 'emptyOnDoseGrid', 1, '-'}; %#ok<AGROW>
        continue;                        % structure smaller than a dose voxel: no metrics
    end
    for q = {'RBE11', 'McN'}
        x = pick(q{1}, phys, rbe11, mcn, v);
        xs = sort(x, 'descend');
        u = unitOf(q{1});
        rows = [rows; {s, q{1}, 'Dmean', mean(x), u; s, q{1}, 'D2', dx(x, 2), u; ...
            s, q{1}, 'D0_03cc', xs(min(n03, numel(xs))), u; s, q{1}, 'Dmax', xs(1), u; ...
            s, q{1}, 'V20Gy_cc', nnz(x >= 20) * vox, 'cc'; s, q{1}, 'V40Gy_cc', nnz(x >= 40) * vox, 'cc'}]; %#ok<AGROW>
        dvh(end+1) = struct('structure', s, 'quantity', q{1}, 'dose', dBins, 'volPct', cumVol(x, dBins)); %#ok<AGROW>
    end
    % LETd where the dose is meaningful
    on = rbe11(v) >= letdMin;
    l = letd(v); l = l(on);
    rows = [rows; {s, 'LETd', 'LETdMean', dwMean(phys(v(on)), l), 'keV/um'; ...
        s, 'LETd', 'LETd_P50', dxOr0(l, 50), 'keV/um'; s, 'LETd', 'LETd_P98', dxOr0(l, 2), 'keV/um'; ...
        s, 'LETd', 'volumeAboveDoseCc', numel(l) * vox, 'cc'}]; %#ok<AGROW>
    lvh(end+1) = struct('structure', s, 'letd', lBins, 'volPct', 100 * cumCount(l, lBins) / max(1, numel(v)), ...
        'volCc', cumCount(l, lBins) * vox); %#ok<AGROW>
    % smoothed D x LETd (voxels >= minDL Gy(RBE))
    m = dl(v); m = m(rbe11(v) >= minDL);
    rows = [rows; {s, 'DxLETd', 'mean', meanOr0(m), 'Gy*keV/um'; s, 'DxLETd', 'D2', dxOr0(m, 2), 'Gy*keV/um'}]; %#ok<AGROW>
    % joint dose-LET hot spots
    for h = 1:numel(hot)
        H = hot{h};
        if ~(strcmp(H.structure, '*') || strcmp(H.structure, s)), continue; end
        D = pick(cns_getfield(H, 'quantity', 'RBE11'), phys, rbe11, mcn, v);
        in = D >= H.doseGy & letd(v) >= H.letd;
        tag = sprintf('V_D%g_L%g', H.doseGy, H.letd);
        rows = [rows; {s, 'hotSpot', [tag '_cc'], nnz(in) * vox, 'cc'; ...
            s, 'hotSpot', [tag '_pct'], 100 * mean(in), '%'}]; %#ok<AGROW>
    end
    if ~isempty(dirty)
        y = dirty(v);
        rows = [rows; {s, 'dirty', 'Dmean', mean(y), 'Gy'; s, 'dirty', 'D2', dx(y, 2), 'Gy'}]; %#ok<AGROW>
    end
end

E = struct('rows', {rows}, 'dvh', dvh, 'lvh', lvh);
E.cubes = struct('RBE11', reshape(rbe11, dims), 'McN', reshape(mcn, dims), 'LETd', reshape(letd, dims));
end

% =========================================================================
function [phys, rbe11, letd, mcn, dirty] = cubesForScenario(dij, s, w, nFx, cstD, abT, abN, letdMin, isPhoton)
dFx = full(dij.physicalDose{s} * w);
phys = dFx * nFx;
letd = zeros(size(dFx));
dirty = [];
if isPhoton
    rbe11 = phys; mcn = phys;            % photons: RBE = 1
    return;
end
rbe11 = 1.1 * phys;
if isfield(dij, 'mLETDose') && ~isempty(dij.mLETDose{s})
    L = full(dij.mLETDose{s} * w);
    ok = rbe11 >= letdMin & dFx > 0;
    letd(ok) = L(ok) ./ dFx(ok);
else
    warning('cns:LET', 'dij has no mLETDose: LET metrics are zero (set calcLET).');
end
ab = abN * ones(size(dFx));
for t = {'CTV', 'PTV', 'GTV'}
    ix = find(strcmp(cstD(:,2), t{1}), 1);
    if ~isempty(ix), ab(cstD{ix,4}{1}) = abT; end
end
mcn = nFx * cns_mcNamara(dFx, letd, ab);
if isfield(dij, 'mDirtyDose') && ~isempty(dij.mDirtyDose{s})
    dirty = full(dij.mDirtyDose{s} * w) * nFx;
end
end

function c = present(cstD, names)
c = names(cellfun(@(n) any(strcmp(cstD(:,2), n)), names));
end

function v = idx(cstD, name)
v = cstD{find(strcmp(cstD(:,2), name), 1), 4}{1};
end

function x = pick(q, phys, rbe11, mcn, v)
switch q
    case 'phys',  x = phys(v);
    case 'RBE11', x = rbe11(v);
    case 'McN',   x = mcn(v);
    otherwise, error('cns:quantity', 'Unknown quantity %s', q);
end
end

function u = unitOf(q)
if strcmp(q, 'phys'), u = 'Gy'; else, u = 'Gy(RBE)'; end
end

function v = dx(x, pct)
% value reached by the hottest pct % of the voxels (D95 = coverage, D2 = near max)
if isempty(x), v = NaN; return; end
x = sort(x(:), 'descend');
v = x(max(1, min(numel(x), ceil(pct / 100 * numel(x)))));
end

function v = cumVol(x, bins)
x = x(:);
if isempty(x), v = zeros(size(bins)); return; end
v = 100 * cumCount(x, bins) / numel(x);
end

function c = cumCount(x, bins)
% number of voxels with value >= each bin edge
x = sort(x(:));
n = numel(x);
c = zeros(size(bins));
for b = 1:numel(bins)
    k = find(x >= bins(b), 1);
    if ~isempty(k), c(b) = n - k + 1; end
end
end

function l = dwMean(dose, letd)
if isempty(dose) || sum(dose) <= 0, l = 0; else, l = sum(dose(:) .* letd(:)) / sum(dose(:)); end
end

function v = dxOr0(x, pct)
% as dx, but 0 when no voxel qualifies (no dose above the threshold = no burden)
if isempty(x), v = 0; else, v = dx(x, pct); end
end

function m = meanOr0(x)
if isempty(x), m = 0; else, m = mean(x); end
end
