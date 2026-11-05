function out = cns_gaussSmooth3(cube, sigmaMm, res)
% cns_gaussSmooth3  Separable 3D Gaussian smoothing (no toolbox needed).
%
%   out = cns_gaussSmooth3(cube, sigmaMm, struct('x',3,'y',3,'z',3))
%
% Cube dimensions are matRad's [y x z]. The kernel is truncated at 3 sigma.
% Used for the smoothed D x LETd predictor of Heinzelmann et al. 2026
% (6 mm Gaussian).
out = double(cube);
if sigmaMm <= 0, return; end
sig = [sigmaMm / res.y, sigmaMm / res.x, sigmaMm / res.z];   % voxels along dims 1,2,3
for d = 1:3
    r = ceil(3 * sig(d));
    if r < 1, continue; end
    k = exp(-0.5 * ((-r:r) / sig(d)).^2);
    k = k / sum(k);
    shape = ones(1, 3); shape(d) = numel(k);
    out = convn(out, reshape(k, shape), 'same');
end
end
