%% matRad_waterPhantom_distalDepthCheck.m
% Post-processing only (no dose calculation, no optimisation).
% Re-evaluates the saved water-phantom plans "behind the target" at several
% depths instead of only the 15 mm region used so far.
%
% WHY: the old distal region was everything downstream of the CTV within 15 mm
% of its surface (a half shell), including the lateral sides and the planning
% margin. For shallow targets the protons stop before 15 mm, and IBA's wider
% lateral penumbra adds many low-dose, low-LET voxels, so the mean LETd there
% is hard to interpret. This script reports, for 2, 5, 10 and 15 mm:
%   shell0toX : the old definition, cut at X mm (X = 15 reproduces the csv)
%   footSlabX : a 3 mm slab centred X mm beyond the distal CTV surface,
%               only inside the beam footprint (lateral distance <= target
%               radius), so lateral penumbra is excluded
% plus a central-axis profile of dose and LETd from 30 mm before to 30 mm
% beyond the distal CTV surface, with the planning-margin edge.
%
% Metrics per region: number of voxels, fraction with dose > 1 % Rx, mean
% dose, mean LETd (voxels with dose > 1 % Rx, as before), dose-weighted LETd,
% mean LET x dose, McNamara mean dose (alpha/beta 2 Gy), V(D>=20 & LETd>=4).
%
% Note: the dose grid is 3 mm, so 2 mm and 5 mm are only about one voxel apart.
% On axis the distal planning margin is 3 mm + 3.5 % of the depth: about
% 4.8 mm (anterior), 8.3 mm (center), 8.8 mm (oblique), 11.8 mm (posterior),
% so 2 and 5 mm behind the CTV are still inside the planning target.
%
% Run on each computer that has plan files (Mac: IBA and Generic in
% WaterPhantom_LETopt_v3; PC: DWA in WaterPhantom_LETopt_v2), from the
% matRadDWA folder:
%   matRad_rc; run('DWA_final_machine/matRad_waterPhantom_distalDepthCheck.m')
% Output: DWA_final_machine/WaterPhantom_distalDepthCheck/
%   distalDepth_<computer>.csv, distalProfile_<computer>.csv

clear; clc;
rootDir = fileparts(mfilename('fullpath'));
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
cd(rootDir);
matRad_rc;

%% ---------------- settings (must match the planning scripts) -------------
inDirs = {fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v2'), ...
          fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v3')};
filePattern   = 'plan_*_cov.mat';   % margin-protected runs (coverage weight 1000)
depths_mm     = [2 5 10 15];
slabHalf_mm   = 1.5;                % slab = X +/- 1.5 mm (one voxel)
res_mm        = 3;
targetDiam_mm = 50;
marginLat_mm  = 3;  rangeUnc = 0.035;
Rx_total      = 54; nFx = 30; abNormal = 2;
voxVol_cc     = (res_mm/10)^3;
if ispc, tagPC = 'PC'; else, tagPC = 'Mac'; end
outDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_distalDepthCheck');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- phantom grid (same as the planning scripts) -------------
semiLR = 250; semiAP = 150; halfSI = 150; airMargin = 30;   % mm
Ny = ceil((2*semiAP + 2*airMargin)/res_mm); Nx = ceil((2*semiLR + 2*airMargin)/res_mm);
Nz = ceil((2*halfSI + 2*airMargin)/res_mm);
ct = struct();
ct.resolution.x = res_mm; ct.resolution.y = res_mm; ct.resolution.z = res_mm;
ct.cubeDim = [Ny Nx Nz];
ct.numOfCtScen = 1;
ct.cube = {zeros(ct.cubeDim)};
ct = matRad_getWorldAxes(ct);
c0 = [mean(ct.x) mean(ct.y) mean(ct.z)];
[Yg, Xg, Zg] = ndgrid(ct.y, ct.x, ct.z);
body = ((Xg - c0(1))/semiLR).^2 + ((Yg - c0(2))/semiAP).^2 <= 1 & abs(Zg - c0(3)) <= halfSI;

%% ---------------- find plan files ----------------------------------------
files = {};
for k = 1:numel(inDirs)
    f = dir(fullfile(inDirs{k}, filePattern));
    files = [files, fullfile({f.folder}, {f.name})]; %#ok<AGROW>
end
if isempty(files), error('No %s files found in the WaterPhantom_LETopt_v2/_v3 folders.', filePattern); end
fprintf('%d plan files found.\n', numel(files));

T = table(); P = table();
R = targetDiam_mm/2;
sProf = (-30:0.5:30)';
for k = 1:numel(files)
    [~, nm] = fileparts(files{k});
    tok = regexp(nm, '^plan_(.+)_(anterior|center|posterior|oblique)_LETw([0-9.eE+-]+)_cov$', 'tokens', 'once');
    if isempty(tok), fprintf('  skipped (name not recognised): %s\n', nm); continue; end
    mName = tok{1}; loc = tok{2};
    S = load(files{k}, 'Dtot', 'LETc', 'tc', 'd', 'wLET');
    Dtot = double(S.Dtot); LET = double(S.LETc); LET(~isfinite(LET)) = 0;
    if ~isequal(size(Dtot), ct.cubeDim)
        fprintf('  skipped (grid %s does not match): %s\n', mat2str(size(Dtot)), nm); continue;
    end
    tc = S.tc(:)'; d = S.d(:)'; d = d / norm(d); wLET = S.wLET;

    % geometry relative to this target and beam
    rr    = sqrt((Xg - tc(1)).^2 + (Yg - tc(2)).^2 + (Zg - tc(3)).^2);
    tMask = rr <= R & body;
    dist  = rr - R;                                         % distance from the CTV surface
    ax    = (Xg - tc(1))*d(1) + (Yg - tc(2))*d(2) + (Zg - tc(3))*d(3);   % along the beam
    pp    = sqrt(max(rr.^2 - ax.^2, 0));                    % lateral distance from the beam axis
    down  = ax > 0;
    depthToCentre = norm(tc - entryPoint(ct, body, tc, d));
    mAx   = marginLat_mm + rangeUnc*depthToCentre;          % distal planning margin on axis
    foot  = body & ~tMask & down & pp <= R;                 % inside the beam footprint
    sBeyond = ax - sqrt(max(R^2 - pp.^2, 0));               % distance beyond the distal CTV surface

    for X = depths_mm
        regs = {'shell0toX', body & ~tMask & down & dist <= X; ...
                'footSlabX', foot & abs(sBeyond - X) <= slabHalf_mm};
        for q = 1:size(regs, 1)
            r = regionMetrics(regs{q,2}, Dtot, LET, Rx_total, nFx, abNormal, voxVol_cc);
            meta = table({mName}, {loc}, wLET, depthToCentre, mAx, {regs{q,1}}, X, ...
                'VariableNames', {'machine','location','LETweight','centreDepth_mm', ...
                                  'PTVdistalMargin_mm','region','depth_mm'});
            T = [T; [meta struct2table(r)]]; %#ok<AGROW>
        end
    end

    % central-axis profile, s = 0 at the distal CTV surface
    pts = tc + (R + sProf) * d;
    FD = griddedInterpolant({ct.y, ct.x, ct.z}, Dtot, 'linear', 'nearest');
    FL = griddedInterpolant({ct.y, ct.x, ct.z}, LET,  'linear', 'nearest');
    n = numel(sProf);
    P = [P; table(repmat({mName}, n, 1), repmat({loc}, n, 1), repmat(wLET, n, 1), repmat(mAx, n, 1), ...
                  sProf, FD(pts(:,2), pts(:,1), pts(:,3)), FL(pts(:,2), pts(:,1), pts(:,3)), ...
        'VariableNames', {'machine','location','LETweight','PTVdistalMargin_mm', ...
                          's_beyondCTV_mm','dose_Gy','LETd_keVum'})]; %#ok<AGROW>
    fprintf('  %-16s %-9s w=%-4g  LETd mean (shell 0-2/5/10/15 mm): %s   foot slab: %s\n', mName, loc, wLET, ...
        sprintf('%.2f ', T.LETd_mean(end-7:2:end)), sprintf('%.2f ', T.LETd_mean(end-6:2:end)));
end

writetable(T, fullfile(outDir, sprintf('distalDepth_%s.csv', tagPC)));
writetable(P, fullfile(outDir, sprintf('distalProfile_%s.csv', tagPC)));
fprintf('\nSaved to %s\n', outDir);

%% ---------------- helpers -------------------------------------------------
function r = regionMetrics(m, Dtot, LET, Rx, nFx, ab, voxVol)
    r = struct('nVox', nnz(m), 'fracDose1pct', NaN, 'Dmean', NaN, 'LETd_mean', NaN, ...
               'LETd_doseWeighted', NaN, 'LETxD_mean', NaN, 'DMcN_mean', NaN, 'V_D20_LET4_cc', NaN);
    if r.nVox == 0, return; end
    Dm = Dtot(m); Lm = LET(m); Pm = Dm / 1.1;              % Dtot is RBE 1.1 dose
    hasD = Dm > 0.01*Rx;
    r.fracDose1pct = mean(hasD);
    r.Dmean = mean(Dm);
    if any(hasD), r.LETd_mean = mean(Lm(hasD)); end
    if sum(Pm) > 0, r.LETd_doseWeighted = sum(Pm .* Lm) / sum(Pm); end
    r.LETxD_mean = mean(Pm .* Lm);
    r.DMcN_mean = mean(mcNamara(Pm, Lm, nFx, ab));
    r.V_D20_LET4_cc = nnz(Dm >= 20 & Lm >= 4) * voxVol;
end

function D = mcNamara(Dphys, LET, nFx, ab)
% McNamara et al. 2015 variable RBE, total RBE-weighted dose
    d = Dphys / nFx;
    RBEmax = 0.99064 + 0.35605 ./ ab .* LET;
    RBEmin = 1.1012 - 0.0038703 .* sqrt(ab) .* LET;
    RBE = (sqrt(ab^2 + 4*d*ab.*RBEmax + 4*d.^2.*RBEmin.^2) - ab) ./ (2*d);
    RBE(d <= 0) = 0;
    D = nFx * d .* RBE;
end

function e = entryPoint(ct, body, tc, d)
    F = griddedInterpolant({ct.y, ct.x, ct.z}, double(body), 'nearest', 'nearest');
    s = 0;
    while F(tc(2) - s*d(2), tc(1) - s*d(1), tc(3) - s*d(3)) > 0.5 && s < 1000
        s = s + 0.5;
    end
    e = tc - s * d;
end
