%% UHN_CNS_09_RayStation_vs_matRad_v1.m
% Compare the approved RayStation DWA proton plan of one patient (Alaina's RTPLAN + RTDOSE)
% with matRad, in two tests:
%   TEST 1 (beam model / dose engine): the RayStation plan itself (spots, energies, MU from
%          the RTPLAN) recalculated in matRad with the pencil beam and MCsquare on the same
%          CT, compared with the RayStation dose. Absolute MU calibration is not available, so
%          matRad is scaled to RayStation by ONE factor (median CTV dose), as agreed.
%   TEST 2 (plan level): the matRad robust dose-only plan (UHN_CNS_07) vs the RayStation plan:
%          clinical goals, DVHs and gamma.
%
% Run on the UHN PC AFTER UHN_CNS_07 has finished for the patient (uses its dij, weights).
%   clear all; matRad_rc; run('CNS_ESTRO/UHN_CNS_09_RayStation_vs_matRad_v1.m')
% Outputs (on hospital storage, <case>_DWAv4rob/results_RayStation_v1/):
%   rs_dicom_summary_v1.txt           what the RTDOSE/RTPLAN files contain (no identifiers)
%   clinicalGoals_<alias>_withRayStation_v1.csv   photon | matRad PB/worst/MC | RayStation | RS plan in matRad
%   comparison_metrics_<alias>_v1.csv  CTV D98/D95/D50/D2, scale factors, gamma pass rates
%   fig_DVH_RayStation_vs_matRad_<alias>_v1.(pdf|png|csv)
%   rs_doses_<alias>_v1.mat           dose cubes (stay on the PC)
% Only the csv/txt numbers may leave the PC.

%% 0  Settings -------------------------------------------------------------------------
S = struct();
S.caseFile = 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data/Patient1/CNS_01.json';
S.rsDir    = '';            % <-- folder with Alaina's RayStation RTPLAN (RP*) and RTDOSE (RD*) files
S.machine  = 'DWA_proton_v4';
S.runMC    = true;          % MCsquare for the RayStation plan (TEST 1) and the matRad plan (TEST 2)
S.nHistMC  = 2e7;           % matRad robust plan (TEST 2)
S.nHistMC_RSplan = 5e7;     % RayStation plan in matRad (TEST 1): ~1e5 spots, needs more protons
S.mcScanner = 'matRad_hlutConsistent';
S.gamma    = {[3 3], [2 2]};  % [% mm], global
S.lowDosePct = 10;          % gamma evaluated where either dose > 10 % of the RayStation max
if exist('S_override', 'var') && isstruct(S_override)
    fo = fieldnames(S_override); for k = 1:numel(fo), S.(fo{k}) = S_override.(fo{k}); end
end
if isempty(S.rsDir), error('Set S.rsDir to the folder with the RayStation RTPLAN and RTDOSE files.'); end

matRadRoot = 'X:\matRadDWA_LET';
if ispc && ~exist(matRadRoot, 'dir'), system('subst X: "M:\My Documents"'); end
if ~exist(matRadRoot, 'dir'), matRadRoot = fileparts(fileparts(mfilename('fullpath'))); end
cd(matRadRoot); matRad_rc; addpath(fullfile(matRadRoot, 'CNS_ESTRO'));
matRad_cfg = MatRad_Config.instance(); matRad_cfg.disableGUI = true;

caseCfg0 = cns_readJson(S.caseFile); cns_checkPrivacy(caseCfg0, S.caseFile);
study = cns_readJson(fullfile(matRadRoot, 'CNS_ESTRO', 'study_UHN_CNS_robustDWA_v1.json'));
alias = caseCfg0.alias; rx = caseCfg0.prescription.totalDoseGyRBE; nFx = caseCfg0.prescription.numOfFractions;
goals = study.clinicalGoals.list;
caseDir = fullfile(cns_expandPath(caseCfg0.workDir), [alias '_DWAv4rob']);
outDir = fullfile(caseDir, 'results_RayStation_v1'); if ~exist(outDir, 'dir'), mkdir(outDir); end
logTxt = fullfile(outDir, 'rs_dicom_summary_v1.txt');
fid = fopen(logTxt, 'w'); if fid < 0, error('Cannot write %s', logTxt); end
say = @(varargin) sayBoth(fid, varargin{:});        % print to the screen and to the summary txt

%% 1  What is in Alaina's folder (headers only, no identifiers) --------------------------
F = dir(fullfile(S.rsDir, '**', '*.dcm'));
rd = {}; rp = {};
for k = 1:numel(F)
    f = fullfile(F(k).folder, F(k).name);
    try, I = dicominfo(f, 'UseDictionaryVR', true); catch, continue; end
    switch I.Modality
        case 'RTDOSE', rd{end+1} = struct('file', f, 'info', I); %#ok<SAGROW>
        case 'RTPLAN', rp{end+1} = struct('file', f, 'info', I); %#ok<SAGROW>
    end
end
say('RayStation folder: %d RTDOSE, %d RTPLAN\n', numel(rd), numel(rp));
if isempty(rd), error('No RTDOSE in %s', S.rsDir); end
for k = 1:numel(rd)
    I = rd{k}.info; bn = NaN;
    try, bn = I.ReferencedRTPlanSequence.Item_1.ReferencedFractionGroupSequence.Item_1.ReferencedBeamSequence.Item_1.ReferencedBeamNumber; catch, end %#ok<NOCOM>
    rd{k}.beam = bn;
    say('RTDOSE %d: DoseType %s, Summation %s, beam %g, units %s, grid %dx%dx%d at %.2f x %.2f x %.2f mm, scaling %g, comment "%s"\n', ...
        k, I.DoseType, I.DoseSummationType, bn, I.DoseUnits, I.Columns, I.Rows, numel(I.GridFrameOffsetVector), ...
        I.PixelSpacing(2), I.PixelSpacing(1), abs(median(diff(I.GridFrameOffsetVector))), I.DoseGridScaling, getf(I, 'DoseComment', ''));
end
% frame of reference: must match the planning CT that matRad imported
forCT = '';
try
    ctF = dir(fullfile(cns_expandPath(caseCfg0.source.dicomDir), '**', '*.dcm'));
    for k = 1:numel(ctF)
        J = dicominfo(fullfile(ctF(k).folder, ctF(k).name));
        if strcmp(J.Modality, 'CT'), forCT = J.FrameOfReferenceUID; break; end
    end
catch
end
forRD = rd{1}.info.FrameOfReferenceUID;
say('Frame of reference RTDOSE vs planning CT: %s\n', ifelse(strcmp(forRD, forCT), 'SAME', ...
    ifelse(isempty(forCT), 'CT not checked (source.dicomDir not found)', 'DIFFERENT - stop and check')));

%% 2  RayStation plan dose on the CT grid ------------------------------------------------
P = load(fullfile(caseDir, 'prepared.mat'), 'ct', 'cst', 'refDose');
ct = P.ct; cst = P.cst; [cx, cy, cz] = ctAxes(ct);
rc = ct.resolution; voxCcCT = rc.x * rc.y * rc.z / 1000;
isPlan = cellfun(@(r) strcmpi(r.info.DoseSummationType, 'PLAN'), rd);
isEff  = cellfun(@(r) strncmpi(r.info.DoseType, 'EFFECTIVE', 6), rd);
if any(isPlan & isEff), use = find(isPlan & isEff, 1); Drs = readDose(rd{use}, cx, cy, cz); src = 'PLAN EFFECTIVE';
elseif any(isPlan),     use = find(isPlan, 1);         Drs = readDose(rd{use}, cx, cy, cz); src = ['PLAN ' rd{use}.info.DoseType];
else
    Drs = 0; for k = 1:numel(rd), Drs = Drs + readDose(rd{k}, cx, cy, cz); end; src = 'sum of BEAM doses';
end
ctvIx = idx(cst, 'CTV');
med = median(Drs(ctvIx));
if med < 0.2 * rx                      % per-fraction dose
    Drs = Drs * nFx; src = [src ' x ' num2str(nFx) ' fractions'];
end
if ~any(isEff), say('NOTE: RayStation dose is PHYSICAL: multiplied by 1.1 for Gy(RBE)\n'); Drs = 1.1 * Drs; end
say('RayStation dose used: %s; CTV median %.2f Gy(RBE), max %.2f Gy(RBE)\n', src, median(Drs(ctvIx)), max(Drs(:)));
% per-beam doses, if present, should add up to the plan dose
isBeam = cellfun(@(r) strcmpi(r.info.DoseSummationType, 'BEAM'), rd);
if any(isBeam) && any(isPlan)
    Db = 0; for k = find(isBeam), Db = Db + readDose(rd{k}, cx, cy, cz); end
    Dp = readDose(rd{use}, cx, cy, cz);
    say('Sum of %d BEAM doses / PLAN dose in the CTV: %.4f\n', nnz(isBeam), median(Db(ctvIx)) / median(Dp(ctvIx)));
    clear Db Dp
end

%% 3  matRad doses: robust plan (PB nominal + MCsquare) -----------------------------------
D = load(fullfile(caseDir, ['dij_' S.machine '.mat']), 'pln', 'stf', 'dij');
W = load(fullfile(caseDir, 'plans', [S.machine '__robustDoseOnly_w.mat']), 'w');
g = D.dij; gridD = struct('doseGrid', g.doseGrid, 'ctGrid', g.ctGrid);
cstD = matRad_resizeCstToGrid(cst, g.ctGrid.x, g.ctGrid.y, g.ctGrid.z, g.doseGrid.x, g.doseGrid.y, g.doseGrid.z);
w = cns_normaliseWeights(W.w, g, cstD, setfield(caseCfg0, 'alias', [alias '_DWAv4rob']), study, P.refDose, cst); %#ok<SFLD>
Dpb_d = reshape(double(full(g.physicalDose{1} * w)), g.doseGrid.dimensions) * 1.1 * nFx;   % dose grid
plnBase = D.pln; stfRob = D.stf; clear D g W
toCT = @(c) matRad_interp3(gridD.doseGrid.x, gridD.doseGrid.y', gridD.doseGrid.z, double(c), ...
    gridD.ctGrid.x, gridD.ctGrid.y', gridD.ctGrid.z, 'linear', 0);
toDG = @(c) matRad_interp3(gridD.ctGrid.x, gridD.ctGrid.y', gridD.ctGrid.z, double(c), ...
    gridD.doseGrid.x, gridD.doseGrid.y', gridD.doseGrid.z, 'linear', 0);
Dpb = toCT(Dpb_d);
Dmc = []; Dmc_d = [];
if S.runMC
    for i = 1:numel(stfRob), stfRob(i).machine = S.machine; end
    Dmc = mcDose(ct, cst, stfRob, plnBase, w, S, nFx, matRadRoot);
    Dmc_d = toDG(Dmc);
end

%% 4  TEST 1: the RayStation plan recalculated in matRad ---------------------------------
Drsm_pb = []; Drsm_mc = []; fPB = NaN; fMC = NaN;
if ~isempty(rp)
    M = load(machineFile(matRadRoot, S.machine), 'machine'); M = M.machine;
    [stfRS, wRS, info] = rsPlanToStf(rp{1}.info, M, S.machine);
    say('RTPLAN: %d beams, %d spots, %d energy layers; energies %.1f-%.1f MeV; max |E(RS) - E(matRad table)| = %.3f MeV\n', ...
        numel(stfRS), numel(wRS), info.nLayers, info.Emin, info.Emax, info.maxEdiff);
    for i = 1:numel(stfRS)
        say('  beam %d: gantry %.1f, couch %.1f, iso [%.1f %.1f %.1f] mm, snout %.1f mm, range shifter %s, spots %d, MU %.2f, spot spacing ~%.1f mm\n', ...
            i, stfRS(i).gantryAngle, stfRS(i).couchAngle, stfRS(i).isoCenter, info.snout(i), info.rashi{i}, ...
            stfRS(i).totalNumOfBixels, info.MU(i), info.spacing(i));
    end
    say('  fractions planned: %g\n', info.nFx);
    plnRS = plnBase; plnRS.multScen = 'nomScen';
    plnRS.propStf.gantryAngles = [stfRS.gantryAngle]; plnRS.propStf.couchAngles = [stfRS.couchAngle];
    plnRS.propStf.isoCenter = vertcat(stfRS.isoCenter); plnRS.propStf.numOfBeams = numel(stfRS);
    plnRS.propDoseCalc.calcLET = 0; plnRS.propDoseCalc.precision = 'double';
    r = matRad_calcDoseForward(ct, cst, stfRS, plnRS, wRS); cd(matRadRoot);
    Drsm_pb = 1.1 * double(r.physicalDose);                     % per fraction, arbitrary MU scale
    fPB = median(Drs(ctvIx)) / median(Drsm_pb(ctvIx)); Drsm_pb = fPB * Drsm_pb;
    say('TEST 1 pencil beam: scale factor matRad -> RayStation (median CTV) = %.4g\n', fPB);
    if S.runMC
        S1 = S; S1.nHistMC = S.nHistMC_RSplan;
        Drsm_mc = mcDose(ct, cst, stfRS, plnRS, wRS, S1, 1, matRadRoot);
        fMC = median(Drs(ctvIx)) / median(Drsm_mc(ctvIx)); Drsm_mc = fMC * Drsm_mc;
        say('TEST 1 MCsquare: scale factor = %.4g (MC/PB factor ratio %.4f)\n', fMC, fPB / fMC);
    end
else
    say('No RTPLAN: TEST 1 skipped\n');
end

%% 5  Clinical goals ----------------------------------------------------------------------
G0 = fullfile(caseDir, 'results_clinicalGoals_v1', sprintf('clinicalGoals_%s_%s_v1.csv', alias, S.machine));
T = readtable(G0);
add = @(T, X, name) addGoals(T, cns_clinicalGoals(X, cst, voxCcCT, goals, name), name);
T.RayStation = []; T = add(T, Drs, 'RayStation');
if ~isempty(Drsm_pb), T = add(T, Drsm_pb, 'RSplan_matRadPB_scaled'); end
if ~isempty(Drsm_mc), T = add(T, Drsm_mc, 'RSplan_matRadMC_scaled'); end
writetable(T, fullfile(outDir, sprintf('clinicalGoals_%s_withRayStation_v1.csv', alias)));
disp(T(:, [{'goal', 'limit', 'Photon', 'PB_nominal', 'MC_nominal', 'RayStation'}, ...
    T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'RSplan') & ~endsWith(T.Properties.VariableNames, '_pass'))]));

%% 6  Metrics and gamma (on the 3 mm dose grid) -------------------------------------------
Drs_d = toDG(Drs);
bodyD = false(gridD.doseGrid.dimensions); bodyD(idx(cstD, 'Body')) = true;
ctvD = false(gridD.doseGrid.dimensions); ctvD(idx(cstD, 'CTV')) = true;
res = [gridD.doseGrid.resolution.x gridD.doseGrid.resolution.y gridD.doseGrid.resolution.z];
pairs = {'RSplan_matRadMC_vs_RayStation', Drsm_mc; 'RSplan_matRadPB_vs_RayStation', Drsm_pb; ...
         'matRadRobustMC_vs_RayStation', Dmc; 'matRadRobustPB_vs_RayStation', Dpb};
rows = {};
dx = @(x, p) prctile(x, 100 - p);
doses = {'RayStation', Drs; 'matRadRobust_PB', Dpb; 'matRadRobust_MC', Dmc; 'RSplan_matRadPB', Drsm_pb; 'RSplan_matRadMC', Drsm_mc};
for k = 1:size(doses, 1)
    if isempty(doses{k, 2}), continue; end
    x = doses{k, 2}(ctvIx);
    rows(end+1, :) = {doses{k, 1}, 'CTV', dx(x, 98), dx(x, 95), dx(x, 50), dx(x, 2), mean(x), NaN, NaN, NaN, NaN}; %#ok<SAGROW>
end
for k = 1:size(pairs, 1)
    if isempty(pairs{k, 2}), continue; end
    A = Drs_d; B = toDG(pairs{k, 2});
    pr = nan(1, 4);
    for q = 1:numel(S.gamma)
        gc = matRad_gammaIndex(A, B, res, S.gamma{q});
        m = (A > S.lowDosePct / 100 * max(A(:)) | B > S.lowDosePct / 100 * max(A(:))) & bodyD;
        pr(2 * q - 1) = 100 * mean(gc(m) < 1);
        pr(2 * q)     = 100 * mean(gc(m & ctvD) < 1);
    end
    rows(end+1, :) = {pairs{k, 1}, 'gamma', NaN, NaN, NaN, NaN, NaN, pr(1), pr(2), pr(3), pr(4)}; %#ok<SAGROW>
    say('gamma %s: 3%%/3mm body %.1f %% (CTV %.1f %%), 2%%/2mm body %.1f %% (CTV %.1f %%)\n', pairs{k, 1}, pr);
end
rows(end+1, :) = {'scaleFactor_RSplan_PB', 'MU->RS', fPB, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
rows(end+1, :) = {'scaleFactor_RSplan_MC', 'MU->RS', fMC, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
Mt = cell2table(rows, 'VariableNames', {'dose', 'what', 'D98_or_value', 'D95', 'D50', 'D2', 'Dmean', ...
    'gamma33_body_pct', 'gamma33_CTV_pct', 'gamma22_body_pct', 'gamma22_CTV_pct'});
writetable(Mt, fullfile(outDir, sprintf('comparison_metrics_%s_v1.csv', alias))); disp(Mt);

%% 7  DVH figure (CT grid) ----------------------------------------------------------------
try
    tR = P.refDose.toTotal; if isnan(tR), tR = nFx; end
    Dph = double(P.refDose.cube) * tR;
    sD = {'CTV', 'PTV', 'GTV', 'Brainstem', 'Chiasm', 'OpticNerve_L', 'OpticNerve_R', 'Pituitary', 'Hippocampus_L', 'Hippocampus_R'};
    sD = sD(cellfun(@(n) ~isempty(idx(cst, n)), sD));
    curves = {'Photon (delivered)', Dph, [0 0 0], '-'; 'RayStation DWA (approved)', Drs, [0.835 0.369 0], '-'; ...
              'matRad DWA robust, MCsquare', Dmc, [0 0.447 0.698], '-'; ...
              'RayStation plan in matRad, MCsquare (scaled)', Drsm_mc, [0.337 0.706 0.914], '--'};
    curves = curves(~cellfun(@isempty, curves(:, 2)), :);
    gx = (0:0.1:65)'; nC = 4; nR = ceil(numel(sD) / nC);
    set(0, 'DefaultFigureVisible', 'off');
    fD = cns_journalFigure(17.5, 4.0 * nR + 0.8); rowsD = {};
    for i = 1:numel(sD)
        v = idx(cst, sD{i}); ax = subplot(nR, nC, i); hold(ax, 'on'); grid(ax, 'on'); h = gobjects(0);
        for c = 1:size(curves, 1)
            y = 100 * arrayfun(@(t) mean(curves{c, 2}(v) >= t), gx);
            h(end+1) = plot(ax, gx, y, curves{c, 4}, 'Color', curves{c, 3}, 'LineWidth', 1.1); %#ok<SAGROW>
            rowsD = [rowsD; [repmat({sD{i}, curves{c, 1}}, numel(gx), 1), num2cell(gx), num2cell(y)]]; %#ok<AGROW>
        end
        xline(ax, rx, ':', 'Color', [0.4 0.4 0.4]); xlim(ax, [0 62]); ylim(ax, [0 102]);
        title(ax, strrep(sD{i}, '_', ' '), 'FontWeight', 'normal');
        if mod(i - 1, nC) == 0, ylabel(ax, 'Volume (%)'); end
        if i > numel(sD) - nC, xlabel(ax, 'Dose (Gy(RBE))'); end
        if i == 1, legend(ax, h, curves(:, 1), 'Location', 'southwest', 'FontSize', 6); end
    end
    cns_saveJournalFigure(fD, fullfile(outDir, sprintf('fig_DVH_RayStation_vs_matRad_%s_v1', alias)), ...
        {'structure', 'dose', 'dose_GyRBE', 'volume_pct'}, rowsD);
    close(fD); set(0, 'DefaultFigureVisible', 'on');
catch err
    set(0, 'DefaultFigureVisible', 'on'); fprintf(2, 'DVH figure failed: %s\n', err.message);
end
save(fullfile(outDir, sprintf('rs_doses_%s_v1.mat', alias)), 'Drs', 'Dpb', 'Dmc', 'Drsm_pb', 'Drsm_mc', 'fPB', 'fMC', '-v7.3');
say('Done. Results in %s\n', outDir);
fclose(fid);

%% ===================================================================================
function sayBoth(fid, varargin)
fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end
function v = getf(s, f, d), if isfield(s, f), v = s.(f); else, v = d; end, end
function o = ifelse(c, a, b), if c, o = a; else, o = b; end, end
function v = idx(cst, n)
k = find(strcmp(cst(:,2), n), 1); if isempty(k), v = []; else, v = cst{k,4}{1}; end
end
function [x, y, z] = ctAxes(ct)
if all(isfield(ct, {'x', 'y', 'z'})) && ~isempty(ct.x)
    x = ct.x(:)'; y = ct.y(:)'; z = ct.z(:)';
else
    o = ct.dicomInfo.ImagePositionPatient; r = ct.resolution;
    x = o(1) + r.x * (0:ct.cubeDim(2) - 1); y = o(2) + r.y * (0:ct.cubeDim(1) - 1);
    z = ct.dicomInfo.SlicePositions(:)'; if numel(z) ~= ct.cubeDim(3), z = o(3) + r.z * (0:ct.cubeDim(3) - 1); end
end
end
function D = readDose(r, cx, cy, cz)
% RTDOSE -> CT grid (matRad cube order y, x, z), Gy as stored (no fraction scaling)
I = r.info; A = double(squeeze(dicomread(r.file))) * double(I.DoseGridScaling);
if ~isequal(round(I.ImageOrientationPatient(:)'), [1 0 0 0 1 0])
    error('RTDOSE orientation %s not supported (HFS expected)', mat2str(I.ImageOrientationPatient'));
end
x = I.ImagePositionPatient(1) + I.PixelSpacing(2) * (0:double(I.Columns) - 1);
y = I.ImagePositionPatient(2) + I.PixelSpacing(1) * (0:double(I.Rows) - 1);
z = I.ImagePositionPatient(3) + double(I.GridFrameOffsetVector(:)');
if numel(z) > 1 && z(2) < z(1), z = fliplr(z); A = flip(A, 3); end
Fi = griddedInterpolant({y, x, z}, A, 'linear', 'none');
D = Fi({cy, cx, cz}); D(isnan(D)) = 0;
end
function f = machineFile(root, m)
f = fullfile(root, 'matRad', 'basedata', ['protons_' m '.mat']);
if ~exist(f, 'file'), f = fullfile(root, 'userdata', 'machines', ['protons_' m '.mat']); end
end
function D = mcDose(ct, cst, stf, pln, w, S, nFx, root)
matRad_cfg = MatRad_Config.instance();
plnMC = pln; plnMC.multScen = 'nomScen';
plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', 0, 'numHistoriesDirect', S.nHistMC, ...
    'doseGrid', pln.propDoseCalc.doseGrid);
matRad_cfg.defaults.MCsquareScanner = S.mcScanner;
r = matRad_calcDoseForward(ct, cst, stf, plnMC, w);
matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'MCsquareScanner'); cd(root);
D = 1.1 * double(r.physicalDose) * nFx;
end
function T = addGoals(T, G, name)
T.(name) = G.value; T.([name '_pass']) = G.pass;
end
function [stf, w, info] = rsPlanToStf(I, machine, mName)
% RayStation ion RTPLAN -> matRad stf (one ray per spot position, energies mapped to the
% nearest machine energy) and w = spot MU (in stf order). Same geometry conventions as
% matRad_importDicomSteeringParticles (spot map (X,Y) -> rayPos_bev = [X 0 Y]).
Eall = [machine.data.energy]; SAD = machine.meta.SAD;
B = I.IonBeamSequence; bn = fieldnames(B);
stf = struct([]); w = []; info = struct('nLayers', 0, 'Emin', Inf, 'Emax', -Inf, 'maxEdiff', 0, ...
    'snout', [], 'rashi', {{}}, 'MU', [], 'spacing', [], 'nFx', NaN);
try, info.nFx = I.FractionGroupSequence.Item_1.NumberOfFractionsPlanned; catch, end %#ok<NOCOM>
muBeam = containers.Map('KeyType', 'double', 'ValueType', 'double');
try
    rb = I.FractionGroupSequence.Item_1.ReferencedBeamSequence; rn = fieldnames(rb);
    for k = 1:numel(rn), muBeam(double(rb.(rn{k}).ReferencedBeamNumber)) = double(rb.(rn{k}).BeamMeterset); end
catch
end
for b = 1:numel(bn)
    Bm = B.(bn{b});
    if isfield(Bm, 'TreatmentDeliveryType') && ~strcmpi(Bm.TreatmentDeliveryType, 'TREATMENT'), continue; end
    cp = Bm.IonControlPointSequence; cn = fieldnames(cp); c1 = cp.(cn{1});
    gantry = double(c1.GantryAngle); couch = double(c1.PatientSupportAngle); iso = double(c1.IsocenterPosition(:)');
    snout = NaN; if isfield(c1, 'SnoutPosition'), snout = double(c1.SnoutPosition); end
    rsTxt = 'none'; rsWET = 0; rsIsoDist = 0;
    if isfield(c1, 'RangeShifterSettingsSequence')
        rs = c1.RangeShifterSettingsSequence.Item_1;
        rsWET = double(getf(rs, 'RangeShifterWaterEquivalentThickness', 0));
        rsIsoDist = double(getf(rs, 'IsocenterToRangeShifterDistance', 0));
        rsTxt = sprintf('%s (WET %.1f mm)', getf(rs, 'RangeShifterSetting', 'IN'), rsWET);
    end
    % cumulative meterset weights -> MU per spot
    cumTot = double(getf(Bm, 'FinalCumulativeMetersetWeight', 1));
    beamMU = cumTot; if isKey(muBeam, double(Bm.BeamNumber)), beamMU = muBeam(double(Bm.BeamNumber)); end
    T = zeros(0, 4);
    for k = 1:numel(cn)
        c = cp.(cn{k});
        if ~isfield(c, 'ScanSpotMetersetWeights'), continue; end
        mw = double(c.ScanSpotMetersetWeights(:)); if sum(mw) == 0, continue; end
        if isfield(c, 'NominalBeamEnergy'), E = double(c.NominalBeamEnergy); end
        n = double(c.NumberOfScanSpotPositions); pm = double(c.ScanSpotPositionMap(:));
        T = [T; pm(1:2:2*n) pm(2:2:2*n) E * ones(n, 1) mw]; %#ok<AGROW>
    end
    T(:, 4) = T(:, 4) / cumTot * beamMU;                       % MU
    % energies -> nearest machine energy
    uE = unique(T(:, 3)); mapE = zeros(size(uE));
    for k = 1:numel(uE)
        [d, j] = min(abs(Eall - uE(k))); mapE(k) = Eall(j); info.maxEdiff = max(info.maxEdiff, d);
    end
    [~, ie] = ismember(T(:, 3), uE); T(:, 3) = mapE(ie);
    info.nLayers = info.nLayers + numel(uE); info.Emin = min(info.Emin, min(uE)); info.Emax = max(info.Emax, max(uE));
    [P, ~, ic] = unique(round(T(:, 1:2) * 1e3) / 1e3, 'rows');
    s = struct('gantryAngle', gantry, 'couchAngle', couch, 'bixelWidth', NaN, 'radiationMode', 'protons', ...
        'machine', mName, 'SAD', SAD, 'isoCenter', iso, 'sourcePoint_bev', [0 -SAD 0]);
    R = matRad_getRotationMatrix(gantry, couch)';
    s.sourcePoint = s.sourcePoint_bev * R;
    ray = struct('rayPos_bev', {}, 'targetPoint_bev', {}, 'rayPos', {}, 'targetPoint', {}, 'energy', {}, ...
        'focusIx', {}, 'rangeShifter', {}, 'weight', {});
    wb = [];
    for j = 1:size(P, 1)
        sel = find(ic == j); [Ej, o] = sort(T(sel, 3)); mu = T(sel(o), 4);
        rp = [P(j, 1) 0 P(j, 2)];
        rsh = repmat(struct('ID', double(rsWET > 0), 'eqThickness', rsWET, 'sourceRashiDistance', SAD - rsIsoDist), 1, numel(Ej));
        ray(j) = struct('rayPos_bev', rp, 'targetPoint_bev', [2 * rp(1) SAD 2 * rp(3)], 'rayPos', rp * R, ...
            'targetPoint', [2 * rp(1) SAD 2 * rp(3)] * R, 'energy', Ej(:)', 'focusIx', ones(1, numel(Ej)), ...
            'rangeShifter', rsh, 'weight', mu(:)'); %#ok<AGROW>
        wb = [wb; mu(:)]; %#ok<AGROW>
    end
    s.ray = ray; s.numOfRays = numel(ray); s.numOfBixelsPerRay = arrayfun(@(r) numel(r.energy), ray);
    s.totalNumOfBixels = sum(s.numOfBixelsPerRay);
    if isempty(stf), stf = s; else, stf(end+1) = s; end %#ok<AGROW>
    w = [w; wb]; %#ok<AGROW>
    info.snout(end+1) = snout; info.rashi{end+1} = rsTxt; info.MU(end+1) = beamMU;
    ux = unique(P(:, 1)); info.spacing(end+1) = median(diff(ux));
end
end
