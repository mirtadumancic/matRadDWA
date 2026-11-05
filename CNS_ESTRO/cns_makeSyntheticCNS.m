function [ct, cst] = cns_makeSyntheticCNS(opt)
% cns_makeSyntheticCNS  Synthetic head with a CNS target next to the brainstem.
%
%   [ct, cst] = cns_makeSyntheticCNS(struct('gapMm', 3, 'resolutionMm', 3))
%
% For testing the pipeline without any patient data. Water-like brain in a
% bone shell, a spherical CTV (with GTV) beside a cylindrical brainstem,
% chiasm, optic nerves, eyes, cochleae and hippocampi. The structure names
% are deliberately inconsistent ('Brain Stem', 'Optic Nrv R', ...) so the
% alias matching in cns_prepareStructures is exercised too.
%
% Options (all optional): gapMm (CTV surface to brainstem surface, 3),
% ctvRadiusMm (18), resolutionMm (3), includeBrain (false: the brain is
% then rebuilt from HU by cns_prepareStructures, as for GLIS-RT).

if nargin < 1, opt = struct(); end
gap  = cns_getfield(opt, 'gapMm', 3);
rCTV = cns_getfield(opt, 'ctvRadiusMm', 18);
res  = cns_getfield(opt, 'resolutionMm', 3);
incBrain = logical(cns_getfield(opt, 'includeBrain', false));

dimsMm = [216 180 168];                 % y (AP), x (LR), z (SI)
n = round(dimsMm / res);
[Y, X, Z] = ndgrid(((1:n(1)) - (n(1)+1)/2) * res, ((1:n(2)) - (n(2)+1)/2) * res, ...
    ((1:n(3)) - (n(3)+1)/2) * res);   % mm from the centre; +y posterior

ell = @(cy, cx, cz, ay, ax, az) ((Y-cy)/ay).^2 + ((X-cx)/ax).^2 + ((Z-cz)/az).^2 <= 1;
cyl = @(cy, cx, r, z0, z1) ((Y-cy).^2 + (X-cx).^2 <= r^2) & Z >= z0 & Z <= z1;

body  = ell(0, 0, 0, 95, 75, 80);
inner = ell(0, 0, 0, 88, 68, 73);
skull = body & ~inner;

rBS = 11;
brainstem = cyl(15, 0, rBS, -60, -5);
cx = -(rCTV + rBS + gap);             % CTV beside the brainstem (patient's right)
ctv = ell(15, cx, -20, rCTV, rCTV, rCTV);
gtv = ell(15, cx, -20, 0.55*rCTV, 0.55*rCTV, 0.55*rCTV);
chiasm = ell(-18, 0, 2, 5, 9, 4);
onL = cyl(-45, 16, 3, -1, 5) | ell(-40, 12, 2, 22, 4, 3.5);
onR = cyl(-45, -16, 3, -1, 5) | ell(-40, -12, 2, 22, 4, 3.5);
eyeL = ell(-68, 32, -2, 12, 12, 12);
eyeR = ell(-68, -32, -2, 12, 12, 12);
cochL = ell(10, 30, -35, 4, 4, 4);
cochR = ell(10, -30, -35, 4, 4, 4);
hipL = ell(5, 24, -15, 8, 6, 15);
hipR = ell(5, -24, -15, 8, 6, 15);

HU = -1000 * ones(n);
HU(body) = 20;            % scalp / soft tissue
HU(skull) = 900;
HU(inner) = 35;           % brain
HU(eyeL | eyeR) = 10;
HU(ctv) = 45;

ct = struct();
ct.cubeHU = {HU};
ct.resolution = struct('x', res, 'y', res, 'z', res);
ct.cubeDim = n;
ct.numOfCtScen = 1;

names = {'CTV_54', 'GTV', 'Brain Stem', 'Chiasm', 'OpticNerve_L', 'Optic Nrv R', 'Eye_L', 'Eye_R', ...
    'Cochlea L', 'Cochlea_R', 'Hippocampus_L', 'Hippo_R', 'External'};
masks = {ctv, gtv, brainstem, chiasm, onL, onR, eyeL, eyeR, cochL, cochR, hipL, hipR, body};
types = [{'TARGET', 'TARGET'}, repmat({'OAR'}, 1, 11)];
if incBrain
    names{end+1} = 'Brain'; masks{end+1} = inner & ~(eyeL | eyeR); types{end+1} = 'OAR';
end

cst = cell(numel(names), 6);
for i = 1:numel(names)
    cst{i,1} = i - 1;
    cst{i,2} = names{i};
    cst{i,3} = types{i};
    cst{i,4} = {find(masks{i})};
    cst{i,5} = struct('Priority', 1 + ~strcmp(types{i}, 'TARGET'), 'alphaX', 0.1, 'betaX', 0.05, ...
        'Visible', 1, 'visibleColor', [1 0 0]);
    cst{i,6} = {};
end
end
