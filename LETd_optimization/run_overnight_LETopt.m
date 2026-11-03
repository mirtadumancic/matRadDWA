function run_overnight_LETopt()
%% run_overnight_LETopt.m
% Runs matRad_waterPhantom_LETopt_forMCSquare.m several times in a row, each time with
% different settings, so the queue below can run unattended overnight.
% Each run gets its own runTag, so nothing is overwritten. If one run fails,
% the error is logged and the next run starts.
%
% Must sit in the same folder as matRad_waterPhantom_LETopt_forMCSquare.m
% (LETd_optimization). Usage, in MATLAB from the matRadDWA_LET folder:
%   matRad_rc; addpath('LETd_optimization'); run_overnight_LETopt
% Results and log: DWA_final_machine/WaterPhantom_LETopt_v2/ (overnight_log.txt)

here = fileparts(mfilename('fullpath'));          % folder of this file (LETd_optimization)
rootDir = here;
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end
src  = fullfile(here, 'matRad_waterPhantom_LETopt_forMCSquare.m');

% ---- the queue: each row overrides settings lines of the v2 script -------
W = '[0 0.1 0.3 1]';
queue = {
  % 1. posterior, margin protected, 4 LET weights
  struct('locations', '{''posterior''}', 'w', W, 'cov', '1000', 'mc', '[]', 'tag', '''_cov''')
  % 2. oblique, margin protected, 4 LET weights
  struct('locations', '{''oblique''}',   'w', W, 'cov', '1000', 'mc', '[]', 'tag', '''_cov''')
  % 3. anterior, margin protected, 4 LET weights (completes DWA at all four targets)
  struct('locations', '{''anterior''}',  'w', W, 'cov', '1000', 'mc', '[]', 'tag', '''_cov''')
  % 4. center, MCsquare recalculation of the dose-only and LET-weight-0.3 plans
  struct('locations', '{''center''}',    'w', '[0 0.3]', 'cov', '1000', 'mc', '[0 0.3]', 'tag', '''_covMC''')
  % 5. posterior, MCsquare recalculation of the same two plans
  struct('locations', '{''posterior''}', 'w', '[0 0.3]', 'cov', '1000', 'mc', '[0 0.3]', 'tag', '''_covMC''')
};

logDir = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v2');   % same folder as the results
if ~exist(logDir, 'dir'), mkdir(logDir); end
logFile = fullfile(logDir, 'overnight_log.txt');
for k = 1:numel(queue)
    q = queue{k};
    msg = sprintf('%s  run %d/%d start: locations %s, LET weights %s, coverage %s, MC %s, tag %s', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), k, numel(queue), q.locations, q.w, q.cov, q.mc, q.tag);
    writeLog(logFile, msg);
    try
        tmp = makeVariant(src, here, k, q);
        runVariant(tmp);
        writeLog(logFile, sprintf('%s  run %d done', datestr(now, 'yyyy-mm-dd HH:MM'), k));
    catch ME
        writeLog(logFile, sprintf('%s  run %d FAILED: %s', datestr(now, 'yyyy-mm-dd HH:MM'), k, ME.message));
    end
    close all force;
end
writeLog(logFile, sprintf('%s  queue finished', datestr(now, 'yyyy-mm-dd HH:MM')));
end

function tmp = makeVariant(src, here, k, q)
% copy of the v2 script with the settings lines replaced
    txt = fileread(src);
    txt = setLine(txt, 'machines', '{''DWA_proton_v2''}');
    txt = setLine(txt, 'locations', q.locations);
    txt = regexprep(txt, 'letWeights = struct\(''DWA_proton_v2'', \[[^\]]*\]', ...
                    ['letWeights = struct(''DWA_proton_v2'', ' q.w], 'once');
    txt = setLine(txt, 'mcRecalcFor', q.mc);
    txt = setLine(txt, 'targetUnderdoseWeight', q.cov);
    txt = setLine(txt, 'runTag', q.tag);
    % the copy lives next to the original so its rootDir detection works
    tmp = fullfile(here, sprintf('tmp_overnight_run%d.m', k));
    fid = fopen(tmp, 'w'); fwrite(fid, txt); fclose(fid);
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
