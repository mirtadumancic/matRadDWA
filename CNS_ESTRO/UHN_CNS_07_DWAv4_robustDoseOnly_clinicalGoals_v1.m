%% UHN_CNS_07_DWAv4_robustDoseOnly_clinicalGoals_v1.m
% Can matRad reproduce the approved RayStation proton plan?
% For one patient: robust dose-only DWA plan with the beam angles of the approved
% RayStation plan, evaluated against the UHN clinical goals.
%
%   1  robust dose-only plan, DWA_proton_v4, pencil beam (cns_runCase with
%      study_UHN_CNS_robustDWA_v1.json): CTV objectives composite worst case over 9
%      scenarios (setup +/-3 mm, range +/-3.5 %), OAR objectives nominal, no PTV,
%      nominal CTV D95 = 100 % Rx. The delivered photon plan is evaluated too.
%   2  clinical goals (study.clinicalGoals) in the nominal scenario and in the worst of
%      the 9 scenarios (pencil-beam dose of every scenario).
%   3  MCsquare recalculation of the nominal plan (same spots and weights): clinical
%      goals again, plus LETd (RayStation-style: every voxel with dose).
%   4  table: goal | limit | PB nominal | PB worst scenario | MCsquare nominal, pass/fail.
%      Fill in the RayStation column by hand from the approved plan's clinical goals.
%
% RUN ON THE UHN PC ONLY (patient data), from X:\matRadDWA_LET (MCsquare needs a path
% without spaces). Outputs go to the case work folder on hospital storage. The goal
% table contains only the study alias, canonical structure names and numbers.
% Time: robust dij ~20-60 min (3 beams x 9 scenarios), optimisation ~10-30 min,
% MCsquare ~10-30 min. Memory: 5 mm spots and lateral cut-off 0.99 (study file) for 16 GB.

%% 0  Settings ---------------------------------------------------------------
S = struct();
S.caseFile  = 'M:/My Documents/CNS_ESTRO_cases/CNS_01.json';  % per patient (made with cns_newCaseFile_v1)
S.studyFile = '';                         % '' = CNS_ESTRO/study_UHN_CNS_robustDWA_v1.json
S.gantryAngles = [0 0 0];                 % <-- the three beams of the approved RayStation plan
S.couchAngles  = [0 0 0];                 %     (matRad convention: see CNS_planning_pipeline / beamCandidates)
S.setLabel  = 'RS_approved_3beam';
S.runMC     = true;                       % MCsquare recalculation of the nominal plan
S.nHistMC   = 2e7;
S.forceDij  = false;                      % true: recompute the robust dij even if cached

if isequal(S.gantryAngles, [0 0 0]), error('Set S.gantryAngles / S.couchAngles to the approved RayStation beams.'); end
matRadRoot = getenv('MATRAD_ROOT');
if isempty(matRadRoot), matRadRoot = fileparts(fileparts(mfilename('fullpath'))); end
cd(matRadRoot); matRad_rc;
matRad_cfg = MatRad_Config.instance(); matRad_cfg.disableGUI = true;
addpath(fullfile(matRadRoot, 'CNS_ESTRO'));
if S.runMC && any(matRadRoot == ' ')
    error('MCsquare fails on paths with spaces: run from X:\\matRadDWA_LET (subst X: "M:\\My Documents").');
end
if isempty(S.studyFile), S.studyFile = fullfile(matRadRoot, 'CNS_ESTRO', 'study_UHN_CNS_robustDWA_v1.json'); end
caseCfg0 = cns_readJson(S.caseFile);
study    = cns_readJson(S.studyFile);
cns_checkPrivacy(caseCfg0, S.caseFile);
mc  = study.machines{1};
rx  = caseCfg0.prescription.totalDoseGyRBE;
nFx = caseCfg0.prescription.numOfFractions;
goals = study.clinicalGoals.list;

caseCfg = caseCfg0;
caseCfg.alias = [caseCfg0.alias '_DWAv4rob'];
caseCfg.beams = struct('setId', S.setLabel, 'gantryAngles', S.gantryAngles, 'couchAngles', S.couchAngles);
caseDir = fullfile(cns_expandPath(caseCfg.workDir), caseCfg.alias);
if ~exist(caseDir, 'dir'), mkdir(caseDir); end
caseFile  = fullfile(caseDir, [caseCfg.alias '.json']);
studyFile = fullfile(caseDir, 'study_robustDWA_used.json');
cns_writeJson(caseFile, caseCfg);
cns_writeJson(studyFile, study);
outDir = fullfile(caseDir, 'results_clinicalGoals_v1');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('[UHN_CNS_07] %s, beams gantry [%s] couch [%s]\n', mc, num2str(S.gantryAngles), num2str(S.couchAngles));

%% 1  Robust dose-only plan (+ photon plan evaluation) ------------------------------
out = cns_runCase(caseFile, studyFile, struct('forceDij', S.forceDij));
planId = study.plans{1}.id;
if ~any(strcmp(out.rows(:,2), mc) & strcmp(out.rows(:,3), planId))
    error('The robust plan failed: see the log above.');
end

%% 2  Clinical goals: nominal and worst scenario (pencil beam) ----------------------
P = load(fullfile(caseDir, 'prepared.mat'));
refDose = []; if isfield(P, 'refDose'), refDose = P.refDose; end
D = load(fullfile(caseDir, ['dij_' mc '.mat']), 'pln', 'stf', 'dij');
W = load(fullfile(caseDir, 'plans', [mc '__' planId '_w.mat']), 'w');
cstD = matRad_resizeCstToGrid(P.cst, D.dij.ctGrid.x, D.dij.ctGrid.y, D.dij.ctGrid.z, ...
    D.dij.doseGrid.x, D.dij.doseGrid.y, D.dij.doseGrid.z);
w = cns_normaliseWeights(W.w, D.dij, cstD, caseCfg, study, refDose, P.cst);
r = D.dij.doseGrid.resolution; voxCc = r.x * r.y * r.z / 1000;
ms = cns_buildPln(P.ct, P.cst, caseCfg, study, mc); ms = ms.multScen;
nS = ms.totNumScen;
Tn = []; vals = nan(numel(cns_toCell(goals)), nS);
for s = 1:nS
    lm = ms.linearMask(s, :);
    Ds = reshape(full(D.dij.physicalDose{lm(1), lm(2), lm(3)} * w), D.dij.doseGrid.dimensions) * 1.1 * nFx;
    T = cns_clinicalGoals(Ds, cstD, voxCc, goals, sprintf('PB scenario %d', s));
    vals(:, s) = T.value;
    if s == 1, Tn = T; Tn.evaluation(:) = {'PB nominal'}; Dnom = Ds; end
end
isMin = cellfun(@(g) strcmp(g.comparator, '>='), cns_toCell(goals));
worst = max(vals, [], 2); worst(isMin) = min(vals(isMin, :), [], 2);
Tw = Tn; Tw.evaluation(:) = {'PB worst scenario'}; Tw.value = worst;
Tw.pass = double((isMin & worst >= Tw.limit - 1e-9) | (~isMin & worst <= Tw.limit + 1e-9));
Tw.pass(isnan(worst)) = NaN;
fprintf('[UHN_CNS_07] pencil beam: %d / %d goals met nominally, %d / %d in the worst scenario\n', ...
    sum(Tn.pass, 'omitnan'), nnz(~isnan(Tn.pass)), sum(Tw.pass, 'omitnan'), nnz(~isnan(Tw.pass)));

% photon plan (delivered), same goals, on the CT grid
Tp = [];
if ~isempty(refDose)
    t = refDose.toTotal; if isnan(t), t = nFx; end
    rc = P.ct.resolution;
    Tp = cns_clinicalGoals(double(refDose.cube) * t, P.cst, rc.x * rc.y * rc.z / 1000, goals, 'Photon delivered');
end

%% 3  MCsquare recalculation of the nominal plan --------------------------------------
Tm = []; letRows = {};
if S.runMC
    plnMC = D.pln; plnMC.multScen = 'nomScen';
    plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', 1, 'numHistoriesDirect', S.nHistMC, ...
        'doseGrid', D.pln.propDoseCalc.doseGrid);
    stf = D.stf; for i = 1:numel(stf), stf(i).machine = mc; end
    matRad_cfg.defaults.doseWeightingThreshold = 0;          % LETd wherever dose > 0 (RayStation-like)
    t0 = tic;
    res = matRad_calcDoseForward(P.ct, P.cst, stf, plnMC, w);
    matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'doseWeightingThreshold');
    cd(matRadRoot);
    if isfield(res, 'RBExD'), Dm = res.RBExD; elseif isfield(res, 'RBExDose'), Dm = res.RBExDose; else, Dm = 1.1 * res.physicalDose; end
    Dm = Dm * nFx; Lm = res.LET; Lm(~isfinite(Lm) | Dm <= 0) = 0;
    g = D.dij;
    if isequal(size(Dm), P.ct.cubeDim) && ~isequal(size(Dm), g.doseGrid.dimensions)   % CT grid -> dose grid
        Dm = matRad_interp3(g.ctGrid.x, g.ctGrid.y', g.ctGrid.z, Dm, g.doseGrid.x, g.doseGrid.y', g.doseGrid.z, 'linear', 0);
        Lm = matRad_interp3(g.ctGrid.x, g.ctGrid.y', g.ctGrid.z, Lm, g.doseGrid.x, g.doseGrid.y', g.doseGrid.z, 'nearest', 0);
    end
    Tm = cns_clinicalGoals(Dm, cstD, voxCc, goals, 'MCsquare nominal');
    fprintf('[UHN_CNS_07] MCsquare (%.0e protons, %.0f s): %d / %d goals met\n', S.nHistMC, toc(t0), ...
        sum(Tm.pass, 'omitnan'), nnz(~isnan(Tm.pass)));
    % LETd (MCsquare, every voxel with dose) and PB-vs-MC dose agreement in the CTV
    for sName = {'CTV', 'GTV', 'Brainstem', 'Chiasm', 'OpticNerve_L', 'OpticNerve_R'}
        ix = find(strcmp(cstD(:,2), sName{1}), 1);
        if isempty(ix), continue; end
        v = cstD{ix,4}{1}; l = Lm(v); dd = Dm(v);
        letRows(end+1, :) = {sName{1}, mean(l(dd > 0)), prctile(l(dd > 0), 98), sum(dd .* l) / max(sum(dd), eps), ...
            mean(Dm(v)) / max(mean(Dnom(v)), eps)}; %#ok<SAGROW>
    end
end

%% 4  Tables ----------------------------------------------------------------------------
G = table(Tn.goal, Tn.limit, Tn.unit, 'VariableNames', {'goal', 'limit', 'unit'});
add = @(G, T, name) addvars(G, T.value, T.pass, 'NewVariableNames', {name, [name '_pass']});
if ~isempty(Tp), G = add(G, Tp, 'Photon'); end
G = add(G, Tn, 'PB_nominal');
G = add(G, Tw, 'PB_worst');
if ~isempty(Tm), G = add(G, Tm, 'MC_nominal'); end
G.RayStation = nan(height(G), 1);                        % fill in from the approved plan
alias = caseCfg0.alias;
writetable(G, fullfile(outDir, sprintf('clinicalGoals_%s_%s_v1.csv', alias, mc)));
if ~isempty(letRows)
    L = cell2table(letRows, 'VariableNames', {'structure', 'LETd_mean_MC', 'LETd_p98_MC', 'LETd_doseWeighted_MC', 'meanDose_MC_over_PB'});
    writetable(L, fullfile(outDir, sprintf('LETd_MCsquare_%s_%s_v1.csv', alias, mc)));
    disp(L);
end
info = struct('alias', alias, 'machine', mc, 'gantryAngles', S.gantryAngles, 'couchAngles', S.couchAngles, ...
    'robustness', study.robustness, 'normalisation', study.normalisation, 'nScenarios', nS, ...
    'bixelWidthMm', D.pln.propStf.bixelWidth, 'doseGridMm', r.x, 'nHistMC', S.nHistMC, 'finished', datestr(now, 31));
cns_writeJson(fullfile(outDir, sprintf('run_info_%s_v1.json', alias)), info);
disp(G);
fprintf('[UHN_CNS_07] done. Table: %s\n', fullfile(outDir, sprintf('clinicalGoals_%s_%s_v1.csv', alias, mc)));
