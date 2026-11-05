function cns_writeJson(f, s)
% cns_writeJson  Write a struct as JSON (pretty-printed where supported).
try
    txt = jsonencode(s, 'PrettyPrint', true);
catch
    txt = jsonencode(s);
end
fid = fopen(f, 'w');
if fid < 0, error('cns:file', 'Cannot write %s', f); end
fprintf(fid, '%s\n', txt);
fclose(fid);
end
