function s = letAgent_readJson(file)
% letAgent_readJson  Read a JSON file into a struct (MATLAB R2016b+ / Octave 7+).
txt = fileread(file);
s = jsondecode(txt);
end
