function v = cns_getfield(s, path, default)
% cns_getfield  Nested field access with a default: cns_getfield(s,'a.b',3).
parts = strsplit(path, '.');
v = s;
for i = 1:numel(parts)
    if isstruct(v) && isfield(v, parts{i}) && ~isempty(v.(parts{i}))
        v = v.(parts{i});
    else
        v = default;
        return;
    end
end
end
