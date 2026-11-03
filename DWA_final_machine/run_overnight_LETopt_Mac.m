function run_overnight_LETopt_Mac()
%% run_overnight_LETopt_Mac.m
% Mac overnight queue (Oct 28, 2026): LET x dose optimisation for IBA_cyclotron
% and Generic at all four targets, with the same settings as the DWA runs on
% the PC (coverage weight 1000, LET weights 0 / 0.1 / 0.3 / 1, tag _cov).
% Runs matRad_waterPhantom_LETopt_v3.m (range shifter for IBA anterior) once
% per target, so each target is logged separately and one failure does not
% stop the others. The v3 script itself is not changed: a temporary copy with
% the settings below is written next to it and run.
%
% Must sit in the same folder as matRad_waterPhantom_LETopt_v3.m
% (DWA_final_machine). Usage, in MATLAB:
%   cd('~/Work/MATLAB_Projects/matRadDWA'); matRad_rc;
%   addpath('DWA_final_machine'); run_overnight_LETopt_Mac
% Results and log: DWA_final_machine/WaterPhantom_LETopt_v3/ (overnight_log_Mac.txt)
% Expected time: about 3-5 h (32 optimisations).

here = fileparts(mfilename('fullpath'));
rootDir = here;
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
src = fullfile(here, 'matRad_waterPhantom_LETopt_v3.m');
if ~exist(src, 'file'), error('matRad_waterPhantom_LETopt_v3.m not found next to this file (%s)', here); end

W = '[0 0.1 0.3 1]';
LW = ['struct(''DWA_proton_v2'', 0, ''IBA_cyclotron'', ' W ', ''Generic'', ' W ')'];
locs = {'center', 'posterior', 'oblique', 'anterior'};     % most informative first

logDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v3');
if ~exist(logDir, 'dir'), mkdir(logDir); end
logFile = fullfile(logDir, 'overnight_log_Mac.txt');
for k = 1:numel(locs)
    writeLog(logFile, sprintf('%s  run %d/%d start: IBA_cyclotron + Generic, %s, LET weights %s, coverage 1000, tag _cov', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), k, numel(locs), locs{k}, W));
    try
        txt = fileread(src);
        txt = setLine(txt, 'machines', '{''IBA_cyclotron'', ''Generic''}');
        txt = setLine(txt, 'locations', ['{''' locs{k} '''}']);
        txt = regexprep(txt, '(?m)^(\s*letWeights\s*=\s*)struct\([\s\S]*?\);', ['$1' LW ';'], 'once');
        if isempty(strfind(txt, LW)), error('letWeights line not found in the v3 script'); end
        txt = setLine(txt, 'mcRecalcFor', '[]');
        txt = setLine(txt, 'mcOpt', 'false');
        txt = setLine(txt, 'targetUnderdoseWeight', '1000');
        txt = setLine(txt, 'runTag', '''_cov''');
        tmp = fullfile(here, sprintf('tmp_overnight_Mac_run%d.m', k));
        fid = fopen(tmp, 'w'); fwrite(fid, txt); fclose(fid);
        runVariant(tmp);
        writeLog(logFile, sprintf('%s  run %d done', datestr(now, 'yyyy-mm-dd HH:MM'), k));
    catch ME
        writeLog(logFile, sprintf('%s  run %d FAILED: %s', datestr(now, 'yyyy-mm-dd HH:MM'), k, ME.message));
    end
    close all force;
end
writeLog(logFile, sprintf('%s  queue finished', datestr(now, 'yyyy-mm-dd HH:MM')));
end

function txt = setLine(txt, name, value)
% replace "name = <anything>;" at the start of a line (first occurrence)
    pat = ['(?m)^(\s*' name '\s*=\s*)[^;\n]*;'];
    if isempty(regexp(txt, pat, 'once'))
        error('Setting "%s" not found in the script', name);
    end
    txt = regexprep(txt, pat, ['$1' strrep(value, '\', '\\') ';'], 'once');
end

function runVariant(tmp)
% separate function: the script's "clear" only clears this workspace
    run(tmp);
end

function writeLog(logFile, msg)
    fprintf('\n#### %s\n', msg);
    fid = fopen(logFile, 'a'); fprintf(fid, '%s\n', msg); fclose(fid);
end
