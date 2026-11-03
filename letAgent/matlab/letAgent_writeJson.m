function letAgent_writeJson(file, s)
% letAgent_writeJson  Write a struct as JSON (NaN/Inf become null).
if exist('OCTAVE_VERSION', 'builtin')
    txt = jsonencode(s);
else
    txt = jsonencode(s, 'PrettyPrint', true);
end
fid = fopen(file, 'w');
if fid < 0, error('letAgent:io', 'Cannot write %s', file); end
fwrite(fid, txt, 'char');
fclose(fid);
end
