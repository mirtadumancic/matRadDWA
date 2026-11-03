function cst = letAgent_checkStructures(cst, caseCfg)
% letAgent_checkStructures  Check every structure named in the case exists
% and set its type (TARGET/OAR) as the case file says.
t = letAgent_toCell(caseCfg.structures.targets);
o = letAgent_toCell(letAgent_getfield(caseCfg, 'structures.oars', {}));
for i = 1:numel(t)
    ix = letAgent_findStructure(cst, t{i}.name);
    cst{ix,3} = 'TARGET';
end
for i = 1:numel(o)
    ix = letAgent_findStructure(cst, o{i}.name);
    cst{ix,3} = 'OAR';
end
for i = 1:size(cst,1)
    cst{i,6} = {};
end
end
