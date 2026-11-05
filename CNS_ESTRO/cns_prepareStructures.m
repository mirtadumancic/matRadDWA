function [cst, info] = cns_prepareStructures(ct, cst, caseCfg, study)
% cns_prepareStructures  Canonical structure names and helper structures.
%
%   [cst, info] = cns_prepareStructures(ct, cst, caseCfg, study)
%
% 1) Renames the structures the study needs to canonical names:
%      CTV PTV GTV Brain Brainstem Chiasm OpticNerve_L/R Eye_L/R Lens_L/R
%      Cochlea_L/R Hippocampus_L/R LacrimalGland_L/R TemporalLobe_L/R Cerebellum
%      Pituitary SpinalCord Body
%    caseCfg.structures.map (canonical -> name as in the RTSTRUCT) wins;
%    otherwise the regular expressions below (or study.structures.aliasRegex)
%    are matched against the lower-case name with spaces/underscores removed.
% 2) Builds what is missing or derived:
%      Body           = CT outline (HU > -400), if there is no BODY/External
%      PTV            = CTV + caseCfg.structures.ptvMarginMm (default 3 mm), if absent
%      Brain          = HU-based intracranial proxy, if absent (flagged in info)
%      Brain_healthy  = Brain minus CTV   (literature "healthy brain")
%      BrainMinusPTV  = Brain minus PTV
%      BrainRing_0_15 = Brain within 15 mm outside the PTV (LET objective shell)
%      Ring_0_5, Ring_5_10, Ring_10_20 = Body shells around the PTV
% 3) Sets TARGET/OAR types, overlap priorities and alpha/beta (10 Gy for
%    targets, 2 Gy otherwise) and clears all objectives.
%
% info.map lists canonical -> original name; info.notes the warnings. Both
% stay in the local case folder (on hospital storage for UHN cases).

info = struct('map', struct(), 'notes', {{}}, 'brainIsProxy', false);
map = cns_getfield(caseCfg, 'structures.map', struct());
regexTab = defaultRegex();
extra = cns_getfield(study, 'structures.aliasRegex', struct());
fn = fieldnames(extra);
for i = 1:numel(fn), regexTab.(fn{i}) = extra.(fn{i}); end

orig = cst(:,2);
norm = cellfun(@(s) regexprep(lower(s), '[^a-z0-9]', ''), orig, 'UniformOutput', false);
canon = fieldnames(regexTab);
used = false(size(orig));
for c = 1:numel(canon)
    name = canon{c};
    ix = [];
    if isfield(map, name)
        % compare letters and digits only: the matRad importer turns '_' into ' ',
        % and RTSTRUCT names often carry trailing spaces
        ix = find(strcmp(norm, regexprep(lower(map.(name)), '[^a-z0-9]', '')) & ~used, 1);
        if isempty(ix)
            error('cns:structure', 'structures.map.%s = "%s" not found. Available: %s', ...
                name, map.(name), strjoin(orig', ', '));
        end
    else
        hits = find(~cellfun(@isempty, regexp(norm, regexTab.(name), 'once')) & ~used);
        if numel(hits) > 1
            if any(strcmp(name, {'CTV', 'PTV', 'GTV'}))
                error('cns:structure', ['Several candidates for %s: %s. Say which one in ' ...
                    'caseCfg.structures.map.%s.'], name, strjoin(orig(hits)', ', '), name);
            end
            info.notes{end+1} = sprintf('%s: several matches (%s), using %s', name, ...
                strjoin(orig(hits)', ', '), orig{hits(1)});
        end
        if ~isempty(hits), ix = hits(1); end
    end
    if ~isempty(ix)
        used(ix) = true;
        cst{ix,2} = name;
        info.map.(name) = orig{ix};
    end
end

if ~hasStruct(cst, 'CTV'), error('cns:structure', 'No CTV found. Set caseCfg.structures.map.CTV.'); end
dims = ct.cubeDim;
if ~hasStruct(cst, 'Body')
    % no BODY/External in the RTSTRUCT: outline from the CT (HU > -400, holes filled)
    body = cns_bodyFromHU(ct);
    if nnz(body) == 0
        error('cns:structure', 'No body contour and the HU-based outline failed. Set structures.map.Body.');
    end
    cst = addStruct(cst, 'Body', 'OAR', body);
    info.notes{end+1} = sprintf('Body: no external contour, outline from CT HU > -400 (%.0f cc)', ...
        nnz(body) * voxCc(ct));
end
body = getMask(cst, 'Body', dims);
ctv  = getMask(cst, 'CTV', dims);

% ---- PTV
if ~hasStruct(cst, 'PTV')
    m = cns_getfield(caseCfg, 'structures.ptvMarginMm', cns_getfield(study, 'structures.ptvMarginMm', 3));
    ptv = grow(ctv, cst, ct, m) & body;
    cst = addStruct(cst, 'PTV', 'TARGET', ptv);
    info.notes{end+1} = sprintf('PTV built as CTV + %g mm', m);
end
ptv = getMask(cst, 'PTV', dims);

% ---- Brain
if ~hasStruct(cst, 'Brain')
    brain = cns_brainProxy(ct, body);
    if nnz(brain) == 0
        error('cns:structure', 'No Brain contour and the HU-based proxy failed. Add a Brain contour.');
    end
    cst = addStruct(cst, 'Brain', 'OAR', brain);
    info.brainIsProxy = true;
    info.notes{end+1} = sprintf('Brain: HU-based intracranial proxy (%.0f cc)', ...
        nnz(brain) * voxCc(ct));
end
brain = getMask(cst, 'Brain', dims);

% ---- derived structures
d5 = grow(ptv, cst, ct, 5);
d10 = grow(d5, cst, ct, 5);
d15 = grow(d10, cst, ct, 5);
d20 = grow(d15, cst, ct, 5);
cst = addStruct(cst, 'Brain_healthy', 'OAR', brain & ~ctv);
if hasStruct(cst, 'BrainMinusPTV')
    info.notes{end+1} = 'BrainMinusPTV: clinical contour kept';
else
    cst = addStruct(cst, 'BrainMinusPTV', 'OAR', brain & ~ptv);
end
cst = addStruct(cst, 'BrainRing_0_15', 'OAR', brain & d15 & ~ptv);
cst = addStruct(cst, 'Ring_0_5', 'OAR', body & d5 & ~ptv);
cst = addStruct(cst, 'Ring_5_10', 'OAR', body & d10 & ~d5);
cst = addStruct(cst, 'Ring_10_20', 'OAR', body & d20 & ~d10);

% ---- types, priorities, alpha/beta, objectives
targets = {'CTV', 'PTV', 'GTV'};
helpers = {'Brain_healthy', 'BrainMinusPTV', 'BrainRing_0_15', 'Ring_0_5', 'Ring_5_10', 'Ring_10_20', 'Body'};
for i = 1:size(cst,1)
    n = cst{i,2};
    p = cst{i,5};
    if ~isstruct(p), p = struct(); end
    if any(strcmp(n, targets))
        cst{i,3} = 'TARGET'; p.Priority = 1; p.alphaX = 0.1; p.betaX = 0.01;
    else
        cst{i,3} = 'OAR'; p.alphaX = 0.1; p.betaX = 0.05;
        if any(strcmp(n, helpers)), p.Priority = 3; else, p.Priority = 2; end
    end
    if ~isfield(p, 'Visible'), p.Visible = 1; end
    cst{i,5} = p;
    cst{i,6} = {};
end
% GTV lies inside the CTV: keep it for evaluation only (not a planning target)
ig = find(strcmp(cst(:,2), 'GTV'));
if ~isempty(ig), cst{ig,3} = 'OAR'; cst{ig,5}.Priority = 1; cst{ig,5}.betaX = 0.01; end

info.volumesCc = struct();
for i = 1:size(cst,1)
    key = matlab.lang.makeValidName(cst{i,2});
    info.volumesCc.(key) = numel(cst{i,4}{1}) * voxCc(ct);
end
for k = 1:numel(info.notes), fprintf('[cns] %s\n', info.notes{k}); end
end

% =========================================================================
function r = defaultRegex()
L = '(l|left|lt|li|sin)'; R = '(r|right|rt|re|dx)';
side = @(stem, s) sprintf('^(%s)%s$|^%s(%s)$', stem, s, s, stem);
r = struct();
r.CTV = '^ctv';
r.PTV = '^ptv';
r.GTV = '^gtv';
r.Brainstem = '^(brainstem|brainstm|bs|tronc|truncus)$';
r.Chiasm = '^(opticchiasm|chiasm|chiasma|chiasmopt)$';
r.OpticNerve_L = side('opticnerve|opticnrv|opticn|on|nopt', L);
r.OpticNerve_R = side('opticnerve|opticnrv|opticn|on|nopt', R);
r.Eye_L = side('eye|globe|retina|bulbus', L);
r.Eye_R = side('eye|globe|retina|bulbus', R);
r.Lens_L = side('lens', L);
r.Lens_R = side('lens', R);
r.Cochlea_L = side('cochlea|cochlear|innerear', L);
r.Cochlea_R = side('cochlea|cochlear|innerear', R);
r.Hippocampus_L = side('hippocampus|hippo|hc', L);
r.Hippocampus_R = side('hippocampus|hippo|hc', R);
r.LacrimalGland_L = side('lacrimalgland|lacrimal|lacgland', L);
r.LacrimalGland_R = side('lacrimalgland|lacrimal|lacgland', R);
r.TemporalLobe_L = side('temporallobe|templobe|temporal', L);
r.TemporalLobe_R = side('temporallobe|templobe|temporal', R);
r.Cerebellum = '^(cerebellum|cerebel)$';
r.Parotid_L = side('parotid|parot', L);
r.Parotid_R = side('parotid|parot', R);
r.Pituitary = '^(pituitary|hypophysis|pit)$';
r.SpinalCord = '^(spinalcord|cord|spinal|myelon)$';
r.Brain = '^(brain|wholebrain|cerebrum|encephalon)$';
r.Body = '^(body|external|ext|skin|outline|head|bodycontour|patient|outercontour)$';
end

function tf = hasStruct(cst, name)
tf = any(strcmp(cst(:,2), name));
end

function m = getMask(cst, name, dims)
m = false(dims);
m(cst{find(strcmp(cst(:,2), name), 1), 4}{1}) = true;
end

function cst = addStruct(cst, name, type, mask)
ix = find(strcmp(cst(:,2), name), 1);
if isempty(ix)
    ix = size(cst,1) + 1;
    cst{ix,1} = ix - 1;
    cst{ix,5} = struct('Priority', 3, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', 1, 'visibleColor', [0 0 1]);
end
cst{ix,2} = name;
cst{ix,3} = type;
cst{ix,4} = {find(mask)};
cst{ix,6} = {};
end

function out = grow(mask, cst, ct, mm)
out = logical(matRad_addMargin(mask, cst, ct.resolution, struct('x', mm, 'y', mm, 'z', mm), true));
end

function v = voxCc(ct)
v = ct.resolution.x * ct.resolution.y * ct.resolution.z / 1000;
end
