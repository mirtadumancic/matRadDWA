function [ct, cst] = letAgent_makeCnsPhantom(opts)
% letAgent_makeCnsPhantom  Synthetic CNS stress case: a spherical target
% next to a brainstem-like box inside a water head. The gap between them
% (gapMm) sets how hard the dose / LET trade-off is, so the search can be
% tested on easy and hard geometries without any patient data.
%
% opts fields (all optional): gapMm (2), targetRadiusMm (15),
% headRadiusMm (75), brainstemSizeMm ([20 20 40]), resolutionMm (3)

if nargin < 1, opts = struct(); end
gap  = letAgent_getfield(opts, 'gapMm', 2);
rT   = letAgent_getfield(opts, 'targetRadiusMm', 15);
rH   = letAgent_getfield(opts, 'headRadiusMm', 75);
bs   = letAgent_getfield(opts, 'brainstemSizeMm', [20 20 40]);
res  = letAgent_getfield(opts, 'resolutionMm', 3);

n = 2 * ceil((rH + 15) / res);
nz = 2 * ceil(max(rH, bs(3)) / res) + 4;
builder = matRad_PhantomBuilder([n n nz], [res res res], 1);

% offsets are [i j k] in mm; the brainstem sits beside the target along j
offBS = rT + gap + bs(2) / 2;
builder.addSphericalTarget('PTV', rT, 'coordType', 'mm', 'offset', [0 -offBS/2 0]);
builder.addBoxOAR('Brainstem', bs(:)', 'coordType', 'mm', 'offset', [0 offBS/2 0]);
builder.addSphericalOAR('BODY', rH, 'coordType', 'mm');
[ct, cst] = builder.getctcst();
ct.cubeHU{1}(ct.cubeHU{1} < -1000) = -1000;
end
