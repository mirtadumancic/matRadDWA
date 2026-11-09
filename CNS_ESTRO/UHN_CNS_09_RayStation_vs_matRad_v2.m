%% UHN_CNS_09_RayStation_vs_matRad_v2.m
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
%   clear all; matRad_rc; run('CNS_ESTRO/UHN_CNS_09_RayStation_vs_matRad_v2.m')
% Outputs (on hospital storage, <case>_DWAv4rob/results_RayStation_v2/):
%   rs_dicom_summary_v2.txt           what the RTDOSE/RTPLAN files contain (no identifiers)
%   clinicalGoals_<alias>_withRayStation_v2.csv   photon | matRad PB/worst/MC | RayStation | RS plan in matRad
%   comparison_metrics_<alias>_v2.csv  CTV D98/D95/D50/D2, scale factors, gamma pass rates
%   MUfit_perEnergyBin_<alias>_v2.csv, MUfit_perBeam_<alias>_v2.csv   fitted MU calibration and its quality
%   fig_DVH_RayStation_vs_matRad_<alias>_v2.(pdf|png|csv)
%   rs_doses_<alias>_v2_PATIENT_DATA_keep_on_PC.mat   dose cubes: patient data, never leaves the PC
% Only the csv/txt numbers may leave the PC.

% v2 (Nov 9 2026), after the first CNS_02 run (v1: gamma 3%/3mm only 35 % for the RayStation plan
% recalculated in matRad, although OAR doses agreed, i.e. the geometry is right):
%   - v1 used the RayStation spot MU directly as matRad particle numbers. MU are not protons:
%     protons per MU depend on the energy (RayStation beam-model calibration). v2 FITS this
%     calibration: per beam, the plan is split into S.nBins energy bins, each recalculated with
%     the pencil beam, and one factor per bin is fitted (non-negative least squares) to that
%     beam's RayStation RTDOSE. Output: MUfit_perEnergyBin_*.csv (the calibration curve).
%     The factors absorb every energy-dependent difference (MU definition, beam-model output).
%   - v1 moved each RayStation energy to the nearest of the 200 machine energies (up to
%     0.65 MeV = ~1 mm range). v2 builds a machine at the exact RayStation energies.
%   - plan-level comparison also at the same CTV D95 as RayStation (CT grid).
%   - CT-grid cubes in single precision; legend outside the DVH panels.
%
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
S.exactEnergies = true;     % machine entries at the exact RayStation energies
S.fitMU    = true;          % fit protons-per-MU vs energy per beam (needs the per-beam RTDOSE files)
S.nBins    = 8;             % energy bins per beam for the fit
S.outsideBodyAir = false;   % true: everything outside the Body contour is air in ALL matRad calculations
                            % (RayStation ignores material outside the External ROI unless it is a
                            % support/fixation ROI; Alaina, Nov 9: mask and couch are not ROIs)
S.outTag   = '';            % appended to the output folder name, e.g. '_airOutside'
S.calcLVH  = false;         % true: MCsquare also scores LETd (no dose threshold) and LETd-volume histograms are written
S.reuseFit = false;         % true: take the MU calibration from results_RayStation_v2/MUfit_perEnergyBin_*.csv (skips the 24 pencil-beam runs)
S.skipUnfitMC = false;      % true: skip the MCsquare run with MU used as particle numbers
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
outDir = fullfile(caseDir, ['results_RayStation_v2' S.outTag]); if ~exist(outDir, 'dir'), mkdir(outDir); end
logTxt = fullfile(outDir, 'rs_dicom_summary_v2.txt');
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
if S.outsideBodyAir
    bodyCT = false(ct.cubeDim); bodyCT(idx(cst, 'Body')) = true;
    nOut = nnz(~bodyCT & ct.cubeHU{1} > -900);
    ct.cubeHU{1}(~bodyCT) = -1000;                     % matRad rebuilds the stopping-power cube from cubeHU
    say('Outside the Body contour set to air: %.0f cc of material (HU > -900) removed\n', ...
        nOut * ct.resolution.x * ct.resolution.y * ct.resolution.z / 1000);
    clear bodyCT
end
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
fxFac = 1; effFac = 1;                 % the same factors are applied to the per-beam doses
if med < 0.2 * rx                      % per-fraction dose
    fxFac = nFx; src = [src ' x ' num2str(nFx) ' fractions'];
end
if ~any(isEff), say('NOTE: RayStation dose is PHYSICAL: multiplied by 1.1 for Gy(RBE)\n'); effFac = 1.1; end
Drs = single(Drs * fxFac * effFac);
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
toCT = @(c) single(matRad_interp3(gridD.doseGrid.x, gridD.doseGrid.y', gridD.doseGrid.z, double(c), ...
    gridD.ctGrid.x, gridD.ctGrid.y', gridD.ctGrid.z, 'linear', 0));
toDG = @(c) matRad_interp3(gridD.ctGrid.x, gridD.ctGrid.y', gridD.ctGrid.z, double(c), ...
    gridD.doseGrid.x, gridD.doseGrid.y', gridD.doseGrid.z, 'linear', 0);
Dpb = toCT(Dpb_d);
Drs_d = toDG(Drs);
bodyD = false(gridD.doseGrid.dimensions); bodyD(idx(cstD, 'Body')) = true;
ctvD = false(gridD.doseGrid.dimensions); ctvD(idx(cstD, 'CTV')) = true;
res = [gridD.doseGrid.resolution.x gridD.doseGrid.resolution.y gridD.doseGrid.resolution.z];
dx = @(x, p) prctile(double(x), 100 - p);
Dmc = []; Dmc_d = []; Lmc = []; Lrs = [];
if S.runMC
    for i = 1:numel(stfRob), stfRob(i).machine = S.machine; end
    [Dmc, Lmc] = mcDose(ct, cst, stfRob, plnBase, w, S, nFx, matRadRoot, S.calcLVH);
    Dmc_d = toDG(Dmc);
end

%% 4  TEST 1: the RayStation plan recalculated in matRad ---------------------------------
Drsm_pb = []; Drsm_mc = []; Drsm_fit = []; Drsm_mcfit = []; fPB = NaN; fMC = NaN; fMCfit = NaN;
if ~isempty(rp)
    M = load(machineFile(matRadRoot, S.machine), 'machine'); M = M.machine;
    [stfRS, wRS, info] = rsPlanToStf(rp{1}.info, M, S.machine, S.exactEnergies);
    say('RTPLAN: %d beams, %d spots, %d energy layers; energies %.1f-%.1f MeV; max |E(RS) - E(matRad table)| = %.3f MeV\n', ...
        numel(stfRS), numel(wRS), info.nLayers, info.Emin, info.Emax, info.maxEdiff);
    mRS = S.machine;
    if S.exactEnergies
        mRS = [S.machine '_RSenergies'];
        machine = machineAtEnergies(M, info.energies); machine.meta.machine = mRS;
        mdir = fullfile(matRad_cfg.primaryUserFolder, 'machines'); if ~exist(mdir, 'dir'), mkdir(mdir); end
        save(fullfile(mdir, ['protons_' mRS '.mat']), 'machine'); clear machine
        for i = 1:numel(stfRS), stfRS(i).machine = mRS; end
        say('Exact energies: machine %s with %d entries at the RayStation energies (range-scaled from the nearest %s entry)\n', ...
            mRS, numel(info.energies), S.machine);
    end
    for i = 1:numel(stfRS)
        say('  beam %d: gantry %.1f, couch %.1f, iso [%.1f %.1f %.1f] mm, snout %.1f mm, range shifter %s, spots %d, MU %.2f, spot spacing ~%.2f mm\n', ...
            info.beamNumber(i), stfRS(i).gantryAngle, stfRS(i).couchAngle, stfRS(i).isoCenter, info.snout(i), info.rashi{i}, ...
            stfRS(i).totalNumOfBixels, info.MU(i), info.spacing(i));
    end
    say('  fractions planned: %g\n', info.nFx);
    plnRS = plnBase; plnRS.machine = mRS; plnRS.multScen = 'nomScen';
    plnRS.propStf.gantryAngles = [stfRS.gantryAngle]; plnRS.propStf.couchAngles = [stfRS.couchAngle];
    plnRS.propStf.isoCenter = vertcat(stfRS.isoCenter); plnRS.propStf.numOfBeams = numel(stfRS);
    plnRS.propDoseCalc.calcLET = 0; plnRS.propDoseCalc.precision = 'double';
    kOfSpot = ones(size(wRS)); reused = false;
    fitCsv = fullfile(caseDir, 'results_RayStation_v2', sprintf('MUfit_perEnergyBin_%s_v2.csv', alias));
    if S.fitMU && S.reuseFit && exist(fitCsv, 'file')
        Kt = readtable(fitCsv); c0 = cumsum([0 stfRS.totalNumOfBixels]);
        for b = 1:numel(stfRS)
            es = cell2mat(arrayfun(@(q) q.energy(:), stfRS(b).ray(:), 'UniformOutput', false));
            kb = Kt(Kt.beam == info.beamNumber(b), :);
            kOfSpot(c0(b) + 1:c0(b + 1)) = kInterp(kb.Emean_MeV, kb.k_fit, es);
        end
        reused = true; say('MU calibration reused from %s\n', fitCsv);
    end
    if S.fitMU && ~reused
        % ---- per beam: pencil-beam dose of each energy bin (dose grid), then fit to the RayStation beam dose
        dimsD = gridD.doseGrid.dimensions; nVox = prod(dimsD);
        gxD = gridD.doseGrid.x(:)'; gyD = gridD.doseGrid.y(:)'; gzD = gridD.doseGrid.z(:)';
        c0 = cumsum([0 stfRS.totalNumOfBixels]);
        Dk1_d = zeros(nVox, 1); Dfit_d = zeros(nVox, 1); kRows = {}; bRows = {};
        matRad_cfg.defaults.calcCubesInterpToCtGrid = false;         % results stay on the dose grid
        try
            for b = 1:numel(stfRS)
                wb = wRS(c0(b) + 1:c0(b + 1));
                es = cell2mat(arrayfun(@(q) q.energy(:), stfRS(b).ray(:), 'UniformOutput', false));
                Eb = unique(es); nBin = min(S.nBins, numel(Eb));
                binOf = ceil((1:numel(Eb))' / numel(Eb) * nBin);
                Bas = zeros(nVox, nBin, 'single'); Ebin = zeros(nBin, 4);
                for j = 1:nBin
                    Ej = Eb(binOf == j); inb = ismember(es, Ej);
                    [sj, wj] = subStf(stfRS(b), wb, Ej);
                    plnB = plnRS; plnB.propStf.gantryAngles = sj.gantryAngle; plnB.propStf.couchAngles = sj.couchAngle;
                    plnB.propStf.isoCenter = sj.isoCenter; plnB.propStf.numOfBeams = 1;
                    r = matRad_calcDoseForward(ct, cst, sj, plnB, wj); cd(matRadRoot);
                    d = double(r.physicalDose); if ~isequal(size(d), dimsD), d = toDG(d); end
                    Bas(:, j) = single(1.1 * d(:));
                    Ebin(j, :) = [min(Ej) max(Ej) sum(es(inb) .* wb(inb)) / sum(wb(inb)) sum(wb(inb))];
                    fprintf('[UHN_CNS_09] beam %d, energy bin %d / %d (%.1f-%.1f MeV, %d spots) done\n', ...
                        info.beamNumber(b), j, nBin, Ebin(j, 1), Ebin(j, 2), nnz(inb));
                end
                kb = find(isBeam & cellfun(@(q) isequal(double(q.beam), info.beamNumber(b)), rd), 1);
                if isempty(kb), error('No per-beam RTDOSE for beam %d: set S.fitMU = false.', info.beamNumber(b)); end
                y = readDose(rd{kb}, gxD, gyD, gzD); y = double(y(:)) * fxFac * effFac;
                m = y > S.lowDosePct / 100 * max(y) & bodyD(:);
                A = double(Bas(m, :));
                k = lsqnonneg(A, y(m));                              % one factor per energy bin
                s1 = sum(A, 2) \ y(m);                               % one factor for the whole beam
                D1 = s1 * double(sum(Bas, 2)); Df = double(Bas) * k;
                g1 = gpass(y, D1, dimsD, res, m); gf = gpass(y, Df, dimsD, res, m);
                Dk1_d = Dk1_d + double(sum(Bas, 2)); Dfit_d = Dfit_d + Df;
                kOfSpot(c0(b) + 1:c0(b + 1)) = kInterp(Ebin(:, 3), k, es);
                for j = 1:nBin
                    kRows(end+1, :) = {info.beamNumber(b), j, Ebin(j, 1), Ebin(j, 2), Ebin(j, 3), Ebin(j, 4), k(j), ...
                        k(j) / s1, k(j) / k(end), (Ebin(j, 3) / Ebin(end, 3))^0.77}; %#ok<SAGROW>
                end
                bRows(end+1, :) = {info.beamNumber(b), s1, 100 * norm(D1(m) - y(m)) / norm(y(m)), ...
                    100 * norm(Df(m) - y(m)) / norm(y(m)), g1(1), g1(2), gf(1), gf(2)}; %#ok<SAGROW>
                say('beam %d: one factor: residual %.1f %%, gamma 3%%/3mm %.1f %%, 2%%/2mm %.1f %% | per-energy factors: residual %.1f %%, gamma %.1f %% / %.1f %%\n', ...
                    info.beamNumber(b), bRows{end, 3}, g1(1), g1(2), bRows{end, 4}, gf(1), gf(2));
                clear Bas A
            end
        catch err
            matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'calcCubesInterpToCtGrid'); rethrow(err);
        end
        matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'calcCubesInterpToCtGrid');
        Dk1_d = reshape(Dk1_d, dimsD); Dfit_d = reshape(Dfit_d, dimsD);
        fPB = median(Drs_d(ctvD)) / median(Dk1_d(ctvD)); Drsm_pb = toCT(fPB * Dk1_d);
        Drsm_fit = toCT(Dfit_d);
        Kt = cell2table(kRows, 'VariableNames', {'beam', 'bin', 'Emin_MeV', 'Emax_MeV', 'Emean_MeV', 'MU', 'k_fit', ...
            'k_rel_to_single_factor', 'k_rel_to_top_bin', 'ref_Epow0p77_rel_to_top_bin'});
        Bt = cell2table(bRows, 'VariableNames', {'beam', 'singleFactor', 'residual_single_pct', 'residual_fit_pct', ...
            'gamma33_single', 'gamma22_single', 'gamma33_fit', 'gamma22_fit'});
        writetable(Kt, fullfile(outDir, sprintf('MUfit_perEnergyBin_%s_v2.csv', alias)));
        writetable(Bt, fullfile(outDir, sprintf('MUfit_perBeam_%s_v2.csv', alias)));
        disp(Kt); disp(Bt);
    elseif ~S.fitMU
        r = matRad_calcDoseForward(ct, cst, stfRS, plnRS, wRS); cd(matRadRoot);
        Drsm_pb = single(1.1 * double(r.physicalDose));
        fPB = median(Drs(ctvIx)) / median(Drsm_pb(ctvIx)); Drsm_pb = fPB * Drsm_pb;
    end
    say('TEST 1 pencil beam, MU used as particle numbers: one scale factor (median CTV) = %.4g\n', fPB);
    if S.runMC
        S1 = S; S1.nHistMC = S.nHistMC_RSplan;
        if ~S.skipUnfitMC
            Drsm_mc = mcDose(ct, cst, stfRS, plnRS, wRS, S1, 1, matRadRoot, false);
            fMC = median(Drs(ctvIx)) / median(Drsm_mc(ctvIx)); Drsm_mc = fMC * Drsm_mc;
            say('TEST 1 MCsquare, MU used as particle numbers: scale factor = %.4g (PB/MC factor ratio %.4f)\n', fMC, fPB / fMC);
        end
        if S.fitMU
            [Drsm_mcfit, Lrs] = mcDose(ct, cst, stfRS, plnRS, wRS .* kOfSpot, S1, 1, matRadRoot, S.calcLVH);
            fMCfit = median(Drs(ctvIx)) / median(Drsm_mcfit(ctvIx)); Drsm_mcfit = fMCfit * Drsm_mcfit;
            say('TEST 1 MCsquare with the fitted MU calibration: residual scale factor (median CTV) = %.4f\n', fMCfit);
        end
    end
else
    say('No RTPLAN: TEST 1 skipped\n');
end
% plan-level comparison at the same CTV D95 as RayStation (CT grid)
Dm0 = Dmc; m0 = 'MC'; if isempty(Dm0), Dm0 = Dpb; m0 = 'PB'; end
fMatch = dx(Drs(ctvIx), 95) / dx(Dm0(ctvIx), 95); Dmatch = fMatch * Dm0;
say('matRad robust plan (%s) scaled by %.4f to the RayStation CTV D95\n', m0, fMatch);

%% 5  Clinical goals ----------------------------------------------------------------------
G0 = fullfile(caseDir, 'results_clinicalGoals_v1', sprintf('clinicalGoals_%s_%s_v1.csv', alias, S.machine));
T = readtable(G0);
add = @(T, X, name) addGoals(T, cns_clinicalGoals(X, cst, voxCcCT, goals, name), name);
T.RayStation = []; T = add(T, Drs, 'RayStation');
if ~isempty(Drsm_pb), T = add(T, Drsm_pb, 'RSplan_matRadPB_scaled'); end
if ~isempty(Drsm_mc), T = add(T, Drsm_mc, 'RSplan_matRadMC_scaled'); end
if ~isempty(Drsm_fit), T = add(T, Drsm_fit, 'RSplan_matRadPB_MUfit'); end
if ~isempty(Drsm_mcfit), T = add(T, Drsm_mcfit, 'RSplan_matRadMC_MUfit'); end
T = add(T, Dmatch, ['matRadRobust' m0 '_sameD95asRS']);
writetable(T, fullfile(outDir, sprintf('clinicalGoals_%s_withRayStation_v2.csv', alias)));
disp(T(:, [{'goal', 'limit', 'Photon', 'PB_nominal', 'MC_nominal', 'RayStation'}, ...
    T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'RSplan') & ~endsWith(T.Properties.VariableNames, '_pass'))]));

%% 6  Metrics and gamma (on the 3 mm dose grid) -------------------------------------------
pairs = {'RSplan_matRadMC_MUfit_vs_RayStation', Drsm_mcfit; 'RSplan_matRadPB_MUfit_vs_RayStation', Drsm_fit; ...
         'RSplan_matRadMC_vs_RayStation', Drsm_mc; 'RSplan_matRadPB_vs_RayStation', Drsm_pb; ...
         ['matRadRobust' m0 '_sameD95asRS_vs_RayStation'], Dmatch; ...
         'matRadRobustMC_vs_RayStation', Dmc; 'matRadRobustPB_vs_RayStation', Dpb};
rows = {};
dx = @(x, p) prctile(x, 100 - p);
doses = {'RayStation', Drs; 'matRadRobust_PB', Dpb; 'matRadRobust_MC', Dmc; ['matRadRobust_' m0 '_sameD95asRS'], Dmatch; ...
         'RSplan_matRadPB', Drsm_pb; 'RSplan_matRadMC', Drsm_mc; 'RSplan_matRadPB_MUfit', Drsm_fit; 'RSplan_matRadMC_MUfit', Drsm_mcfit};
for k = 1:size(doses, 1)
    if isempty(doses{k, 2}), continue; end
    x = double(doses{k, 2}(ctvIx));
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
rows(end+1, :) = {'residualScale_RSplan_MC_MUfit', 'after fit', fMCfit, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
rows(end+1, :) = {'scale_matRadRobust_to_RS_D95', m0, fMatch, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
Mt = cell2table(rows, 'VariableNames', {'dose', 'what', 'D98_or_value', 'D95', 'D50', 'D2', 'Dmean', ...
    'gamma33_body_pct', 'gamma33_CTV_pct', 'gamma22_body_pct', 'gamma22_CTV_pct'});
writetable(Mt, fullfile(outDir, sprintf('comparison_metrics_%s_v2.csv', alias))); disp(Mt);

%% 6b  LETd-volume histograms (MCsquare, every voxel counts, no dose threshold; CT grid) ---------
if S.calcLVH
    gL = (0:0.1:25)'; rowsL = {}; mL = {};
    sL = {'CTV', 'PTV', 'GTV', 'Brainstem', 'Chiasm', 'OpticNerve_L', 'OpticNerve_R', 'Pituitary', ...
          'Hippocampus_L', 'Hippocampus_R', 'Cochlea_L', 'Cochlea_R', 'BrainMinusPTV', 'Brain'};
    sL = sL(cellfun(@(n) ~isempty(idx(cst, n)), sL));
    setsL = {'matRad robust plan MCsquare', Lmc, Dmc; 'RayStation plan in matRad MCsquare', Lrs, Drsm_mcfit};
    for q = 1:size(setsL, 1)
        if isempty(setsL{q, 2}), continue; end
        for i = 1:numel(sL)
            v = idx(cst, sL{i}); l = double(setsL{q, 2}(v)); d = double(setsL{q, 3}(v));
            vol = 100 * arrayfun(@(x) mean(l >= x), gL);
            rowsL = [rowsL; [repmat({setsL{q, 1}, sL{i}}, numel(gL), 1), num2cell(gL), num2cell(vol)]]; %#ok<AGROW>
            mL(end+1, :) = {setsL{q, 1}, sL{i}, mean(l), sum(l .* d) / max(sum(d), eps), prctile(l, 50), prctile(l, 98), mean(d)}; %#ok<SAGROW>
        end
    end
    if ~isempty(rowsL)
        writetable(cell2table(rowsL, 'VariableNames', {'plan', 'structure', 'LETd_keVum', 'volume_pct'}), ...
            fullfile(outDir, sprintf('LVH_MCsquare_%s_v2.csv', alias)));
        Lt = cell2table(mL, 'VariableNames', {'plan', 'structure', 'LETd_mean', 'LETd_doseWeighted', 'LETd_median', 'LETd_p98', 'meanDose_GyRBE'});
        writetable(Lt, fullfile(outDir, sprintf('LETd_metrics_MCsquare_%s_v2.csv', alias))); disp(Lt);
    end
end

%% 7  DVH figure (CT grid) ----------------------------------------------------------------
try
    tR = P.refDose.toTotal; if isnan(tR), tR = nFx; end
    Dph = single(double(P.refDose.cube) * tR);
    sD = {'CTV', 'PTV', 'GTV', 'Brainstem', 'Chiasm', 'OpticNerve_L', 'OpticNerve_R', 'Pituitary', 'Hippocampus_L', 'Hippocampus_R'};
    sD = sD(cellfun(@(n) ~isempty(idx(cst, n)), sD));
    Dre = Drsm_mcfit; lre = 'RayStation plan in matRad (MCsquare, fitted MU calibration)';
    if isempty(Dre), Dre = Drsm_fit; lre = 'RayStation plan in matRad (pencil beam, fitted MU calibration)'; end
    if isempty(Dre), Dre = Drsm_mc; lre = 'RayStation plan in matRad (MCsquare, MU as particles)'; end
    if isempty(Dre), Dre = Drsm_pb; lre = 'RayStation plan in matRad (pencil beam, MU as particles)'; end
    curves = {'Photon (delivered)', Dph, [0 0 0], '-'; 'RayStation DWA plan', Drs, [0.835 0.369 0], '-'; ...
              sprintf('matRad DWA robust (%s), same CTV D95', m0), Dmatch, [0 0.447 0.698], '-'; ...
              lre, Dre, [0.337 0.706 0.914], '--'};
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
        if i == numel(sD)                                  % legend in a free panel (or below the last one)
            if numel(sD) < nR * nC
                axL = subplot(nR, nC, nR * nC); hold(axL, 'on'); hl = gobjects(0);
                for c = 1:size(curves, 1), hl(end+1) = plot(axL, NaN, NaN, curves{c, 4}, 'Color', curves{c, 3}, 'LineWidth', 1.1); end %#ok<SAGROW>
                axis(axL, 'off'); legend(axL, hl, curves(:, 1), 'Location', 'west', 'FontSize', 6);
            else
                legend(ax, h, curves(:, 1), 'Location', 'southoutside', 'FontSize', 6);
            end
        end
    end
    cns_saveJournalFigure(fD, fullfile(outDir, sprintf('fig_DVH_RayStation_vs_matRad_%s_v2', alias)), ...
        {'structure', 'dose', 'dose_GyRBE', 'volume_pct'}, rowsD);
    close(fD); set(0, 'DefaultFigureVisible', 'on');
catch err
    set(0, 'DefaultFigureVisible', 'on'); fprintf(2, 'DVH figure failed: %s\n', err.message);
end
save(fullfile(outDir, sprintf('rs_doses_%s_v2_PATIENT_DATA_keep_on_PC.mat', alias)), 'Lmc', 'Lrs', 'Drs', 'Dpb', 'Dmc', 'Dmatch', 'Drsm_pb', 'Drsm_mc', 'Drsm_fit', 'Drsm_mcfit', 'fPB', 'fMC', 'fMCfit', 'fMatch', '-v7.3');
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
function [D, L] = mcDose(ct, cst, stf, pln, w, S, nFx, root, wantLET)
% MCsquare forward dose (Gy(RBE), x nFx) and, if wantLET, LETd in every voxel with dose
if nargin < 9, wantLET = false; end
matRad_cfg = MatRad_Config.instance();
plnMC = pln; plnMC.multScen = 'nomScen';
plnMC.propDoseCalc = struct('engine', 'MCsquare', 'calcLET', double(wantLET), 'numHistoriesDirect', S.nHistMC, ...
    'doseGrid', pln.propDoseCalc.doseGrid);
matRad_cfg.defaults.MCsquareScanner = S.mcScanner;
if wantLET, matRad_cfg.defaults.doseWeightingThreshold = 0; end
try
    r = matRad_calcDoseForward(ct, cst, stf, plnMC, w);
catch err
    matRad_cfg.defaults = rmfield(matRad_cfg.defaults, intersect(fieldnames(matRad_cfg.defaults), {'MCsquareScanner', 'doseWeightingThreshold'}));
    cd(root); rethrow(err);
end
matRad_cfg.defaults = rmfield(matRad_cfg.defaults, intersect(fieldnames(matRad_cfg.defaults), {'MCsquareScanner', 'doseWeightingThreshold'}));
cd(root);
D = single(1.1 * nFx * double(r.physicalDose)); L = [];
if wantLET
    if isfield(r, 'LET'), L = single(r.LET); L(~isfinite(L) | D <= 0) = 0;
    else, warning('MCsquare returned no LET cube.'); end
end
end
function T = addGoals(T, G, name)
T.(name) = G.value; T.([name '_pass']) = G.pass;
end
function p = gpass(a, b, dims, res, m)
% gamma pass rates [3%/3mm, 2%/2mm] (global) of b against a, in the voxels m
p = [NaN NaN]; crit = {[3 3], [2 2]}; m = reshape(m, dims);
for q = 1:2
    gc = matRad_gammaIndex(reshape(a, dims), reshape(b, dims), res, crit{q});
    p(q) = 100 * mean(gc(m) < 1);
end
end
function kq = kInterp(Ec, k, Eq)
% factor per spot: linear in energy between the bin centres, constant outside
if numel(Ec) < 2, kq = k(1) * ones(size(Eq)); return; end
kq = interp1(Ec, k, Eq, 'linear');
kq(Eq <= Ec(1)) = k(1); kq(Eq >= Ec(end)) = k(end); kq(isnan(kq)) = k(end);
end
function [s, w] = subStf(s0, w0, Eset)
% one beam restricted to the spots with energies in Eset (w0: weights of the beam, stf order)
c = 0; keepRay = false(1, s0.numOfRays); ray = s0.ray; w = [];
for j = 1:s0.numOfRays
    n = numel(s0.ray(j).energy); wj = w0(c + (1:n)); c = c + n;
    keep = ismember(s0.ray(j).energy, Eset);
    if any(keep)
        keepRay(j) = true;
        ray(j).energy = ray(j).energy(keep); ray(j).focusIx = ray(j).focusIx(keep);
        ray(j).rangeShifter = ray(j).rangeShifter(keep); ray(j).weight = ray(j).weight(keep);
        w = [w; wj(keep)]; %#ok<AGROW>
    end
end
s = s0; s.ray = ray(keepRay); s.numOfRays = nnz(keepRay);
s.numOfBixelsPerRay = arrayfun(@(q) numel(q.energy), s.ray);
s.totalNumOfBixels = sum(s.numOfBixelsPerRay);
end
function machine = machineAtEnergies(machine, Enew)
% machine entries at the energies Enew: copy of the nearest entry with the depth axis scaled
% by the range ratio (R ~ E^1.77) and Z scaled so that the integral dose stays ~ E
Eall = [machine.data.energy]; Enew = sort(Enew(:))'; data = machine.data(1:0);
for k = 1:numel(Enew)
    [~, j] = min(abs(Eall - Enew(k))); d = machine.data(j);
    fE = Enew(k) / d.energy; fR = fE^1.77;
    d.depths = d.depths * fR; d.peakPos = d.peakPos * fR;
    if isnumeric(d.Z), d.Z = d.Z * fE / fR; end
    if isfield(d, 'energySpectrum') && isstruct(d.energySpectrum) && isfield(d.energySpectrum, 'mean')
        d.energySpectrum.mean = Enew(k);
    end
    d.energy = Enew(k);
    data(k) = d;
end
machine.data = data;
end
function [stf, w, info] = rsPlanToStf(I, machine, mName, exact)
% RayStation ion RTPLAN -> matRad stf (one ray per spot position) and w = spot MU (stf order).
% exact = true: the RayStation energies are kept; false: nearest machine energy.
% Geometry as matRad_importDicomSteeringParticles (spot map (X,Y) -> rayPos_bev = [X 0 Y]).
Eall = [machine.data.energy]; SAD = machine.meta.SAD;
B = I.IonBeamSequence; bn = fieldnames(B);
stf = struct([]); w = [];
info = struct('nLayers', 0, 'Emin', Inf, 'Emax', -Inf, 'maxEdiff', 0, 'snout', [], 'MU', [], 'spacing', [], ...
    'nFx', NaN, 'beamNumber', [], 'energies', []);
info.rashi = {};
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
    cumTot = double(getf(Bm, 'FinalCumulativeMetersetWeight', 1));
    beamMU = cumTot; if isKey(muBeam, double(Bm.BeamNumber)), beamMU = muBeam(double(Bm.BeamNumber)); end
    T = zeros(0, 4); E = NaN;
    for k = 1:numel(cn)
        c = cp.(cn{k});
        if isfield(c, 'NominalBeamEnergy'), E = double(c.NominalBeamEnergy); end
        if ~isfield(c, 'ScanSpotMetersetWeights'), continue; end
        mw = double(c.ScanSpotMetersetWeights(:)); if sum(mw) == 0, continue; end
        n = double(c.NumberOfScanSpotPositions); pm = double(c.ScanSpotPositionMap(:));
        T = [T; pm(1:2:2*n) pm(2:2:2*n) E * ones(n, 1) mw]; %#ok<AGROW>
    end
    T(:, 4) = T(:, 4) / cumTot * beamMU;                       % MU per spot
    T(:, 3) = round(T(:, 3) * 1e4) / 1e4;
    uE = unique(T(:, 3)); mapE = uE;
    for k = 1:numel(uE)
        [d, j] = min(abs(Eall - uE(k))); info.maxEdiff = max(info.maxEdiff, d);
        if ~exact, mapE(k) = Eall(j); end
    end
    info.nLayers = info.nLayers + numel(uE); info.Emin = min(info.Emin, min(uE)); info.Emax = max(info.Emax, max(uE));
    [~, ie] = ismember(T(:, 3), uE); T(:, 3) = mapE(ie);
    info.energies = union(info.energies, unique(T(:, 3)));
    % spot spacing: median nearest-neighbour distance in the layer with the most spots
    [~, ~, il] = unique(T(:, 3)); [~, kmax] = max(accumarray(il, 1)); Pl = T(il == kmax, 1:2); sp = NaN;
    if size(Pl, 1) > 1
        Dm = sqrt((Pl(:, 1) - Pl(:, 1)').^2 + (Pl(:, 2) - Pl(:, 2)').^2); Dm(1:size(Dm, 1) + 1:end) = Inf;
        sp = median(min(Dm, [], 2));
    end
    [P, ~, ic] = unique(round(T(:, 1:2) * 1e3) / 1e3, 'rows');
    s = struct('gantryAngle', gantry, 'couchAngle', couch, 'bixelWidth', NaN, 'radiationMode', 'protons', ...
        'machine', mName, 'SAD', SAD, 'isoCenter', iso, 'sourcePoint_bev', [0 -SAD 0]);
    R = matRad_getRotationMatrix(gantry, couch)';
    s.sourcePoint = s.sourcePoint_bev * R;
    ray = struct('rayPos_bev', {}, 'targetPoint_bev', {}, 'rayPos', {}, 'targetPoint', {}, 'energy', {}, ...
        'focusIx', {}, 'rangeShifter', {}, 'weight', {});
    wb = [];
    for j = 1:size(P, 1)
        sel = find(ic == j); [Ej, ~, iu] = unique(T(sel, 3)); mu = accumarray(iu, T(sel, 4));
        rp = [P(j, 1) 0 P(j, 2)];
        rsh = repmat(struct('ID', double(rsWET > 0), 'eqThickness', rsWET, 'sourceRashiDistance', SAD - rsIsoDist), 1, numel(Ej));
        ray(j) = struct('rayPos_bev', rp, 'targetPoint_bev', [2 * rp(1) SAD 2 * rp(3)], 'rayPos', rp * R, ...
            'targetPoint', [2 * rp(1) SAD 2 * rp(3)] * R, 'energy', Ej(:)', 'focusIx', ones(1, numel(Ej)), ...
            'rangeShifter', rsh, 'weight', mu(:)'); %#ok<AGROW>
        wb = [wb; mu(:)]; %#ok<AGROW>
    end
    s.ray = ray; s.numOfRays = numel(ray); s.numOfBixelsPerRay = arrayfun(@(q) numel(q.energy), ray);
    s.totalNumOfBixels = sum(s.numOfBixelsPerRay);
    if isempty(stf), stf = s; else, stf(end+1) = s; end %#ok<AGROW>
    w = [w; wb]; %#ok<AGROW>
    info.snout(end+1) = snout; info.rashi{end+1} = rsTxt; info.MU(end+1) = beamMU;
    info.spacing(end+1) = sp; info.beamNumber(end+1) = double(Bm.BeamNumber);
end
end
