function body = cns_bodyFromHU(ct, thresholdHU)
% cns_bodyFromHU  Patient outline from the CT, for cases without a BODY/External contour.
%
%   body = cns_bodyFromHU(ct)          % threshold -400 HU
%
% Per axial slice: voxels above the threshold, holes filled (sinuses, airways
% stay inside), then the largest 3D connected component (drops the couch and
% immobilisation mask if they are not attached). Check it in the GUI once.
if nargin < 2, thresholdHU = -400; end
if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU)
    HU = double(ct.cubeHU{1});
else
    HU = (double(ct.cube{1}) - 1) * 1000;
end
m = HU > thresholdHU;
for k = 1:size(m, 3)
    m(:, :, k) = fillHoles2D(m(:, :, k));
end
body = largestComponent(m);
end

function f = fillHoles2D(b)
if ~any(b(:)), f = b; return; end
if exist('imfill', 'file') == 2
    f = imfill(b, 'holes');
    return;
end
bg = false(size(b));
bg([1 end], :) = ~b([1 end], :);
bg(:, [1 end]) = ~b(:, [1 end]);
k = [0 1 0; 1 1 1; 0 1 0];
while true
    nb = conv2(double(bg), k, 'same') > 0 & ~b;
    if isequal(nb, bg), break; end
    bg = nb;
end
f = ~bg;
end

function m = largestComponent(m)
if ~any(m(:)), return; end
if exist('bwconncomp', 'file') == 2
    cc = bwconncomp(m, 6);
    [~, i] = max(cellfun(@numel, cc.PixelIdxList));
    keep = false(size(m)); keep(cc.PixelIdxList{i}) = true; m = keep;
elseif exist('bwlabeln', 'file') == 2
    L = bwlabeln(m, 6);
    n = accumarray(L(L > 0), 1);
    [~, i] = max(n);
    m = L == i;
end
end
