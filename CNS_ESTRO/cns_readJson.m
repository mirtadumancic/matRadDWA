function s = cns_readJson(f)
% cns_readJson  Read a JSON file (MATLAB R2016b+ / Octave 7+).
if ~exist(f, 'file')
    error('cns:file', 'File not found: %s', f);
end
s = jsondecode(fileread(f));
end
