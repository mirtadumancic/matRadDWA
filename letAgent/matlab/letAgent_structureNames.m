function names = letAgent_structureNames(caseCfg)
% letAgent_structureNames  Structure names used by the case (targets, then OARs).
names = {};
t = letAgent_toCell(caseCfg.structures.targets);
o = letAgent_toCell(letAgent_getfield(caseCfg, 'structures.oars', {}));
for i = 1:numel(t), names{end+1} = t{i}.name; end %#ok<AGROW>
for i = 1:numel(o), names{end+1} = o{i}.name; end %#ok<AGROW>
end
