function ix = letAgent_findStructure(cst, name)
% letAgent_findStructure  Row of cst whose name matches (case-insensitive).
ix = find(strcmpi(cst(:,2), name), 1);
if isempty(ix)
    error('letAgent:structure', 'Structure "%s" not found. Available: %s', name, strjoin(cst(:,2)', ', '));
end
end
