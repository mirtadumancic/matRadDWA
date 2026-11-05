function cns_writeCsvTable(f, hdr, T)
% cns_writeCsvTable  Write a header + cell table (numbers or text) as csv.
fid = fopen(f, 'w');
if fid < 0, error('cns:file', 'Cannot write %s', f); end
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(T, 1)
    parts = cell(1, size(T, 2));
    for j = 1:size(T, 2)
        v = T{i, j};
        if ischar(v), parts{j} = v; else, parts{j} = sprintf('%.6g', v); end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end
