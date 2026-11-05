function brain = cns_brainProxy(ct, body)
% cns_brainProxy  Intracranial soft tissue from the CT, for cases without a Brain contour.
%
%   brain = cns_brainProxy(ct, bodyMask)
%
% Per axial slice: bone (HU > 250) is hole-filled; the enclosed region that
% is not bone and has brain/CSF-like HU (-30..100) is kept if it is larger
% than 10 cm^2. The largest 3D connected component is returned. Slices
% where the skull is open (skull base, foramen magnum) drop out, so the
% proxy is conservative inferiorly. GLIS-RT has no Brain contour; the
% Burdenko set does, and a contoured brain is always preferred.

if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU)
    HU = double(ct.cubeHU{1});
else
    HU = (double(ct.cube{1}) - 1) * 1000;     % rough rED -> HU, enough for thresholds
end
dims = size(HU);
bone = HU > 250 & body;
soft = HU > -30 & HU < 100 & body;
minVox = 10 * 100 / (ct.resolution.x * ct.resolution.y);   % 10 cm^2
brain = false(dims);
for k = 1:dims(3)
    b = bone(:,:,k);
    if nnz(b) < 20, continue; end
    filled = fillHoles2D(b);
    s = filled & ~b & soft(:,:,k);
    if nnz(s) >= minVox, brain(:,:,k) = s; end
end
brain = largestComponent(brain);
end

function f = fillHoles2D(b)
if exist('imfill', 'file') == 2
    f = imfill(b, 'holes');
    return;
end
% Flood the background in from the border; what it cannot reach is a hole.
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
