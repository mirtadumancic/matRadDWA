function R = cns_letWeightSearch(dij, cst, pln, caseCfg, study, wDose, refDose)
% cns_letWeightSearch  Gradient descent on the LET x dose penalty weight for one machine.
%
%   R = cns_letWeightSearch(dij, cst, pln, caseCfg, study, wDoseOnly, refDose)
%
% Settings: study.letSearch (see study_UHN_CNS_v1.json for the definition).
%   Inner problem : matRad fluence optimisation with the study's dose objectives
%                   plus SquaredOverdosing(w, LxD_ref) on matRad's 'LETd'
%                   quantity (= LET x dose) in letStructures. Every run starts
%                   from the dose-only spot weights (wInit), so the outcome
%                   is a deterministic function of w.
%   Outer problem : minimise over u = log10(w)
%                     J(u) = mean_s B_s(u)/B_s(dose-only) + lambda * P_cov(u)
%                   B_s = burdenQuantity/burdenMetric of structure s (default
%                   smoothed D x LETd, D2), P_cov = squared coverage
%                   violations in percent (see letSearch.coverage).
%   Method        : steepest descent with a central finite-difference
%                   gradient (step h in decades), step -eta*g capped at
%                   maxStepDecades, Armijo backtracking, w kept in
%                   [wMin, wMax]. If a finite-difference neighbour is
%                   better than the line-search point, the search moves
%                   there (J is noisy at the optimiser tolerance). It stops
%                   when the current point beats its neighbours (local
%                   optimum at resolution h) or after maxEvals plans.
%                   Evaluations are cached; the reported optimum is the
%                   best J among ALL evaluated weights.
%   Caveat        : J(u) is not guaranteed convex or smooth (IPOPT stops at
%                   a tolerance, coverage terms are hinges), so this finds a
%                   local optimum at the resolution of h. The table of all
%                   evaluations is returned so the trade-off can be inspected.
%
% R.table  one row per evaluated weight (cell + header)
% R.best   w, J, weights (normalised), metric rows
% R.base   dose-only metric rows, normalised weights
% R.ref    LET x dose reference (total, keV/um*Gy)

C = cns_getfield(study, 'letSearch', struct());
X = struct();                                 % constants of the search
X.w0 = cns_getfield(C, 'w0', 0.1);
X.uMin = log10(cns_getfield(C, 'wMin', 1e-3)); X.uMax = log10(cns_getfield(C, 'wMax', 3));
h = cns_getfield(C, 'h', 0.15);
eta = cns_getfield(C, 'eta', 1);
maxStep = cns_getfield(C, 'maxStepDecades', 0.5);
tolStep = cns_getfield(C, 'tolStepDecades', 0.05);
cA = cns_getfield(C, 'armijoC', 1e-4);
maxBT = cns_getfield(C, 'maxBacktracks', 3);
X.maxEvals = cns_getfield(C, 'maxEvals', 14);
X.lambda = cns_getfield(C, 'lambda', 1);
X.burdenS = cns_toCell(cns_getfield(C, 'burdenStructures', {'Chiasm', 'Brainstem'}));
X.bQ = cns_getfield(C, 'burdenQuantity', 'DxLETd'); X.bM = cns_getfield(C, 'burdenMetric', 'D2');
X.cov = cns_getfield(C, 'coverage', struct());
X.covS = cns_getfield(X.cov, 'structure', 'CTV');
X.dispS = cns_toCell(cns_getfield(C, 'displacement', {}));
X.letS = cns_toCell(cns_getfield(C, 'letStructures', cns_getfield(study, 'letStructures', {})));
X.rx = caseCfg.prescription.totalDoseGyRBE;
nFx = caseCfg.prescription.numOfFractions;
X.dij = dij; X.cst = cst; X.pln = pln; X.caseCfg = caseCfg; X.study = study; X.wDose = wDose; X.refDose = refDose;
X.cstD = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
    dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

% ---- baseline: dose-only plan (normalised exactly like every LET plan)
[wB, fB] = cns_normaliseWeights(wDose, dij, X.cstD, caseCfg, study, refDose, cst);
EB = cns_evaluatePlan(dij, cst, wB, caseCfg, study);
X.B0 = burdenValues(EB.rows, X);
X.active = X.B0 > 0 & ~isnan(X.B0);
if ~any(X.active)
    warning('cns:search', 'No burden structure has a nonzero %s %s in the dose-only plan: nothing to optimise.', X.bQ, X.bM);
end
X.D98_0 = get1(EB.rows, X.covS, 'RBE11', 'D98');
X.D2_0 = get1(EB.rows, X.covS, 'RBE11', 'D2');

% ---- LET x dose reference from the dose-only plan (total course, physical dose)
X.ref = 0;
refMode = cns_getfield(C, 'refMode', 'targetP5');
if ischar(refMode) && strcmp(refMode, 'targetP5')
    iC = find(strcmp(X.cstD(:,2), X.covS), 1);
    LxD = full(dij.mLETDose{1} * wB) * nFx;          % LET x dose, total course
    x = sort(LxD(X.cstD{iC,4}{1}));
    X.ref = cns_getfield(C, 'refFactor', 1) * x(max(1, ceil(0.05 * numel(x))));
    clear LxD
elseif isnumeric(refMode)
    X.ref = refMode;
end
fprintf('[cns] LET x dose reference = %.1f keV/um*Gy (total course; refMode %s, factor %g)\n', ...
    X.ref, num2str(refMode), cns_getfield(C, 'refFactor', 1));

hdr = {'evaluation', 'w', 'log10w', 'J', 'burdenRatio', 'Pcov', 'CTV_D98', 'CTV_D2', 'CTV_D95', ...
    'optimizerStatus', 'seconds'};
for s = 1:numel(X.burdenS), hdr{end+1} = sprintf('%s_%s_%s', X.burdenS{s}, X.bQ, X.bM); end %#ok<AGROW>
for s = 1:numel(X.dispS)
    hdr{end+1} = [X.dispS{s} '_RBE11_Dmean']; hdr{end+1} = [X.dispS{s} '_DxLETd_mean']; %#ok<AGROW>
end
hdr{end+1} = 'BrainRing_0_15_V_D20_L4_cc';

% search state: cache of evaluated weights and the table
T = struct('cacheU', [], 'cacheJ', [], 'cacheW', {{}}, 'cacheRows', {{}}, 'tab', {{}}, 'nEval', 0);

% ---- baseline row (w = 0)
[J0, Br0, Pc0, m0] = outer(EB.rows, X);
T.tab(end+1, :) = [{0, 0, -Inf, J0, Br0, Pc0, m0(1), m0(2), m0(3), 'doseOnly', 0}, ...
    num2cell(X.B0), num2cell(dispValues(EB.rows, X)), {get1(EB.rows, 'BrainRing_0_15', 'hotSpot', 'V_D20_L4_cc')}];

% ---- steepest descent in u = log10(w)
u = min(max(log10(X.w0), X.uMin), X.uMax);
[Ju, T] = evalJ(u, X, T);
trace = {};
while T.nEval < X.maxEvals
    up = min(u + h, X.uMax); um = max(u - h, X.uMin);
    [Jp, T] = evalJ(up, X, T);
    [Jm, T] = evalJ(um, X, T);
    g = (Jp - Jm) / (up - um);
    if abs(g) < 1e-6
        trace{end+1} = sprintf('u = %.3f: gradient ~ 0, stop', u); %#ok<AGROW>
        break;
    end
    step = -eta * g;
    step = sign(step) * min(abs(step), maxStep);
    accepted = false; uN = u; JN = Inf;
    for bt = 0:maxBT
        if T.nEval >= X.maxEvals, break; end
        uTry = min(max(u + step, X.uMin), X.uMax);
        if abs(uTry - u) < 1e-9, break; end
        [JN, T] = evalJ(uTry, X, T);
        if JN <= Ju - cA * abs(step) * abs(g)
            accepted = true; uN = uTry; break;
        end
        step = step / 2;
    end
    % move to the best point evaluated so far if it beats the line-search result
    % (the finite-difference neighbours count: J is noisy at the optimiser tolerance)
    [Jbest, kb] = min(T.cacheJ);
    if (~accepted || Jbest < JN - 1e-9) && Jbest < Ju - 1e-9
        uN = T.cacheU(kb); accepted = true;
    end
    trace{end+1} = sprintf('u = %.3f (w = %.4g): J = %.4f, dJ/du = %.4f, move %+.3f decades, %s', ...
        u, 10^u, Ju, g, uN - u, ifelse(accepted, 'accepted', 'no better point: local optimum at step h')); %#ok<AGROW>
    fprintf('[cns]   %s\n', trace{end});
    if ~accepted, break; end
    small = abs(uN - u) < tolStep;
    u = uN; [Ju, T] = evalJ(u, X, T);
    if small && Ju <= min(T.cacheJ) + 1e-9
        trace{end+1} = sprintf('u = %.3f: move below %.3f decades, stop', u, tolStep); %#ok<AGROW>
        break;
    end
end
if T.nEval >= X.maxEvals, trace{end+1} = sprintf('stopped: maxEvals = %d reached', X.maxEvals); end

[~, kb] = min(T.cacheJ);
R = struct();
R.header = hdr;
R.table = T.tab;
R.trace = trace;
R.ref = X.ref;
R.base = struct('w', wB, 'normalisationFactor', fB, 'rows', {EB.rows}, 'J', J0);
R.best = struct('w', 10^T.cacheU(kb), 'log10w', T.cacheU(kb), 'J', T.cacheJ(kb), ...
    'weights', T.cacheW{kb}, 'rows', {T.cacheRows{kb}});
R.nEvaluations = T.nEval;
R.hitBound = abs(T.cacheU(kb) - X.uMin) < 1e-3 || abs(T.cacheU(kb) - X.uMax) < 1e-3;
if J0 <= T.cacheJ(kb)
    R.note = 'No LET weight improved J over the dose-only plan.';
elseif R.hitBound
    R.note = 'Best weight is at a search bound: widen wMin/wMax.';
else
    R.note = '';
end
fprintf('[cns] best w = %.4g (J = %.4f; dose-only J = %.4f) after %d evaluations. %s\n', ...
    R.best.w, R.best.J, J0, T.nEval, R.note);
end

% =========================================================================
function [J, T] = evalJ(u, X, T)
u = min(max(u, X.uMin), X.uMax);
k = find(abs(T.cacheU - u) < 1e-3, 1);
if ~isempty(k), J = T.cacheJ(k); return; end
T.nEval = T.nEval + 1;
w = 10^u;
t0 = tic;
P = struct('id', sprintf('LETxD_w%.4g', w), 'let', struct('mode', 'LETxDose', 'penalty', w, ...
    'ref', X.ref, 'structures', {X.letS}));
cstRun = cns_buildCst(X.cst, X.caseCfg, X.study, P, X.pln);
[res, ~, optUsed] = cns_optimize(X.dij, cstRun, X.pln, X.study, X.wDose);
if ~strcmpi(optUsed, 'IPOPT'), fprintf(2, '[cns] w = %.4g optimised with %s\n', w, optUsed); end
wO = cns_normaliseWeights(res.wUnsequenced, X.dij, X.cstD, X.caseCfg, X.study, X.refDose, X.cst);
E = cns_evaluatePlan(X.dij, X.cst, wO, X.caseCfg, X.study);
[J, Bratio, Pc, mets] = outer(E.rows, X);
st = '';
try, st = num2str(res.info.status); catch, end %#ok<NOCOM>
T.tab(end+1, :) = [{T.nEval, w, u, J, Bratio, Pc, mets(1), mets(2), mets(3), st, toc(t0)}, ...
    num2cell(burdenValues(E.rows, X)), num2cell(dispValues(E.rows, X)), ...
    {get1(E.rows, 'BrainRing_0_15', 'hotSpot', 'V_D20_L4_cc')}];
T.cacheU(end+1) = u; T.cacheJ(end+1) = J; T.cacheW{end+1} = wO; T.cacheRows{end+1} = E.rows;
fprintf('[cns]   eval %2d: w = %.4g  J = %.4f  (burden ratio %.3f, Pcov %.3f, %s D98 %.2f, D2 %.2f)\n', ...
    T.nEval, w, J, Bratio, Pc, X.covS, mets(1), mets(2));
end

function [J, Bratio, Pc, mets] = outer(rows, X)
b = burdenValues(rows, X);
if any(X.active), Bratio = mean(b(X.active) ./ X.B0(X.active)); else, Bratio = 1; end
D98 = get1(rows, X.covS, 'RBE11', 'D98'); D2 = get1(rows, X.covS, 'RBE11', 'D2');
D95 = get1(rows, X.covS, 'RBE11', 'D95');
v1 = max(0, 100 * (cns_getfield(X.cov, 'minFractionOfRx', 0.95) * X.rx - D98) / X.rx);
if isfield(X, 'D98_0')
    v2 = max(0, 100 * (X.D98_0 - D98) / X.D98_0 - cns_getfield(X.cov, 'maxD98DropPct', 1));
    v3 = max(0, 100 * (D2 - X.D2_0) / X.D2_0 - cns_getfield(X.cov, 'maxD2RisePct', 2));
else
    v2 = 0; v3 = 0;
end
Pc = v1^2 + v2^2 + v3^2;
J = Bratio + X.lambda * Pc;
mets = [D98 D2 D95];
end

function b = burdenValues(rows, X)
b = zeros(1, numel(X.burdenS));
for q = 1:numel(X.burdenS), b(q) = get1(rows, X.burdenS{q}, X.bQ, X.bM); end
end

function d = dispValues(rows, X)
d = zeros(1, 2 * numel(X.dispS));
for q = 1:numel(X.dispS)
    d(2*q - 1) = get1(rows, X.dispS{q}, 'RBE11', 'Dmean');
    d(2*q) = get1(rows, X.dispS{q}, 'DxLETd', 'mean');
end
end

function v = get1(rows, s, q, m)
sel = strcmp(rows(:,1), s) & strcmp(rows(:,2), q) & strcmp(rows(:,3), m);
if any(sel), v = rows{find(sel, 1), 4}; else, v = NaN; end
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end
