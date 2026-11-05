function [ct, cst, refDose] = cns_makeSyntheticUHNlike(opt)
% cns_makeSyntheticUHNlike  Synthetic head that mimics the UHN CNS case layout (no patient data).
%
%   [ct, cst, refDose] = cns_makeSyntheticUHNlike(struct('resolutionMm', 3))
%
% Target behind the LEFT eye (patient left = +x in matRad's LPS system,
% shown on the RIGHT of an axial image), abutting the left optic nerve.
% Structure names are those of the UHN RTSTRUCT (GTVp, CTVp_5400,
% PTVp_5400, OpticNrv_L, x_Ring1, CT SETUP, ...), and there is NO body
% contour, like the real case. refDose is an analytic photon-like dose
% (total, Gy) standing in for RD1: 54 Gy to the PTV with a smooth fall-off.
% For testing UHN_CNS_01/02/03 anywhere, including on an AI-connected machine.

if nargin < 1, opt = struct(); end
res = cns_getfield(opt, 'resolutionMm', 3);
rx  = cns_getfield(opt, 'prescriptionGy', 54);

dimsMm = [216 180 168];                 % y (AP), x (LR), z (SI)
n = round(dimsMm / res);
[Y, X, Z] = ndgrid(((1:n(1)) - (n(1)+1)/2) * res, ((1:n(2)) - (n(2)+1)/2) * res, ...
    ((1:n(3)) - (n(3)+1)/2) * res);     % mm from the centre; +x patient left, +y posterior, +z superior
ell = @(cy, cx, cz, ay, ax, az) ((Y-cy)/ay).^2 + ((X-cx)/ax).^2 + ((Z-cz)/az).^2 <= 1;
cyl = @(cy, cx, r, z0, z1) ((Y-cy).^2 + (X-cx).^2 <= r^2) & Z >= z0 & Z <= z1;

body  = ell(0, 0, 0, 95, 75, 80);
inner = ell(0, 0, 0, 88, 68, 73);
skull = body & ~inner;
eyeL = ell(-68, 32, -2, 12, 12, 12);   eyeR = ell(-68, -32, -2, 12, 12, 12);
lensL = ell(-78, 32, -2, 2.5, 4, 4);   lensR = ell(-78, -32, -2, 2.5, 4, 4);
onL = ell(-40, 14, 0, 20, 3, 3);        onR = ell(-40, -14, 0, 20, 3, 3);
chiasm = ell(-18, 0, 2, 5, 9, 4);
pit = ell(-12, 0, -10, 4, 5, 4);
brainstem = cyl(15, 0, 11, -60, -5);
cochL = ell(10, 30, -35, 4, 4, 4);      cochR = ell(10, -30, -35, 4, 4, 4);
hipL = ell(5, 24, -15, 8, 6, 15);       hipR = ell(5, -24, -15, 8, 6, 15);
parL = ell(5, 50, -50, 12, 8, 14);      parR = ell(5, -50, -50, 12, 8, 14);
cord = cyl(20, 0, 5, -80, -60);
c = [-44, 27, 0];                       % CTV centre behind the left eye (y, x, z)
ctv = ell(c(1), c(2), c(3), 13, 13, 13);
gtv = ell(c(1), c(2), c(3), 7, 7, 7);
ptv = ell(c(1), c(2), c(3), 16, 16, 16) & body;
brain = inner & ~(eyeL | eyeR | onL | onR) & Z > -60;
brainMinusPTV = brain & ~ptv;
ring1 = ell(c(1), c(2), c(3), 26, 26, 26) & ~ptv & body;
ring2 = ell(c(1), c(2), c(3), 40, 40, 40) & ~ell(c(1), c(2), c(3), 26, 26, 26) & body;

HU = -1000 * ones(n);
HU(body) = 20; HU(skull) = 900; HU(inner) = 35;
HU(eyeL | eyeR) = 10; HU(ctv) = 45;

ct = struct('cubeHU', {{HU}}, 'resolution', struct('x', res, 'y', res, 'z', res), ...
    'cubeDim', n, 'numOfCtScen', 1);

names = {'GTVp', 'CTVp_5400', 'PTVp_5400', 'Brain', 'BrainMinusPTV', 'Brainstem', 'Chiasm', ...
    'Cochlea_L', 'Cochlea_R', 'Eye_L', 'Eye_R', 'Hippocampus_L', 'Hippocampus_R', 'Lens_L', 'Lens_R', ...
    'OpticNrv_L', 'OpticNrv_R', 'Parotid_L', 'Parotid_R', 'Pituitary', 'SpinalCord', 'x_Ring1', 'x_Ring2', ...
    'CT SETUP', 'ISO_A', 'ICRU_A'};
empty = false(n);
masks = {gtv, ctv, ptv, brain, brainMinusPTV, brainstem, chiasm, cochL, cochR, eyeL, eyeR, hipL, hipR, ...
    lensL, lensR, onL, onR, parL, parR, pit, cord, ring1, ring2, empty, empty, empty};
cst = cell(numel(names), 6);
for i = 1:numel(names)
    cst{i,1} = i - 1; cst{i,2} = names{i};
    if i <= 3, cst{i,3} = 'TARGET'; else, cst{i,3} = 'OAR'; end
    cst{i,4} = {find(masks{i})};
    cst{i,5} = struct('Priority', 2, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1, 'visibleColor', [1 0 0]);
    cst{i,6} = {};
end

% analytic photon-like dose (total Gy): PTV at Rx, smooth fall-off over ~1.5 cm
d = cns_gaussSmooth3(double(ell(c(1), c(2), c(3), 24, 24, 24)), 4, ct.resolution);
x = sort(d(ptv), 'descend');
d = d / x(ceil(0.95 * numel(x))) * 0.97 * rx;            % PTV D95 = 97 % Rx, like a clinical plan
d = min(d, 1.05 * rx);
d(~body) = 0;
refDose = struct('cube', d, 'toTotal', 1, 'summationType', 'PLAN', 'doseType', 'physicalDose', ...
    'numOfFractions', 30, 'sourceFields', {{'physicalDose'}});
end
