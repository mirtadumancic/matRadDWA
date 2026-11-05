function [rbeDoseFx, rbe] = cns_mcNamara(dFx, letd, ab)
% cns_mcNamara  McNamara et al. (PMB 2015) variable-RBE-weighted dose.
%
%   [rbeDoseFx, rbe] = cns_mcNamara(dFx, LETd, alphaBeta)
%
% dFx   physical dose per fraction [Gy]
% LETd  dose-averaged LET [keV/um]
% ab    (alpha/beta)_x [Gy], scalar or one value per voxel
%
%   RBEmax = 0.99064 + 0.35605/(a/b) * LETd
%   RBEmin = 1.1012  - 0.0038703 * sqrt(a/b) * LETd
%   RBE*d  = 0.5 * ( sqrt((a/b)^2 + 4 d (a/b) RBEmax + 4 d^2 RBEmin^2) - a/b )
%
% Returns the RBE-weighted dose per fraction (multiply by the number of
% fractions for the total) and the RBE (RBEmax where d = 0).
dFx = double(dFx); letd = double(letd);
if isscalar(ab), ab = ab * ones(size(dFx)); end
rbeMax = 0.99064 + 0.35605 ./ ab .* letd;
rbeMin = 1.1012 - 0.0038703 .* sqrt(ab) .* letd;
rbeDoseFx = 0.5 * (sqrt(ab.^2 + 4 * dFx .* ab .* rbeMax + 4 * dFx.^2 .* rbeMin.^2) - ab);
rbeDoseFx(dFx <= 0) = 0;
if nargout > 1
    rbe = rbeMax;
    pos = dFx > 0;
    rbe(pos) = rbeDoseFx(pos) ./ dFx(pos);
end
end
