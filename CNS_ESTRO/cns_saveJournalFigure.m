function cns_saveJournalFigure(f, base, dataHeader, dataRows)
% cns_saveJournalFigure  Save a figure as vector PDF + 300 dpi PNG, and its numbers as csv.
%
%   cns_saveJournalFigure(f, 'D:/runs/CNS_01/figures/fig_dvh', header, rows)
%
% header : 1xN cell of column names
% rows   : MxN cell (numbers or text) - the exact values that are plotted
% Every figure gets its csv (project rule: no chart without its numbers).
d = fileparts(base);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
isOctave = exist('OCTAVE_VERSION', 'builtin') > 0;
try
    if ~isOctave && exist('exportgraphics', 'file')
        exportgraphics(f, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
        exportgraphics(f, [base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
    else
        print(f, [base '.pdf'], '-dpdf', '-painters');
        print(f, [base '.png'], '-dpng', '-r300');
    end
catch err
    % fallback (e.g. macOS invisible figure reporting a zero-size canvas)
    try
        set(f, 'Visible', 'on'); drawnow;
        print(f, [base '.pdf'], '-dpdf', '-painters');
        print(f, [base '.png'], '-dpng', '-r300');
    catch err2
        warning('cns:figure', 'Could not save %s: %s / %s', base, err.message, err2.message);
    end
end
if nargin > 2 && ~isempty(dataHeader)
    fid = fopen([base '.csv'], 'w');
    fprintf(fid, '%s\n', strjoin(dataHeader, ','));
    for i = 1:size(dataRows, 1)
        parts = cell(1, size(dataRows, 2));
        for j = 1:size(dataRows, 2)
            v = dataRows{i, j};
            if ischar(v), parts{j} = v; else, parts{j} = sprintf('%.6g', v); end
        end
        fprintf(fid, '%s\n', strjoin(parts, ','));
    end
    fclose(fid);
end
close(f);
end
