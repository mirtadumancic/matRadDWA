function test_MCsquare_PC()
%% test_MCsquare_PC.m
% Quick (~10 min) check that MCsquare runs on this computer before the
% overnight queue. Runs matRad_waterPhantom_LETopt_v2.m once with test
% settings: DWA, center target, dose-only plan (LET weight 0), fast
% optimisation (coverage weight 10), then an MCsquare recalculation with
% 1e6 histories. matRad_waterPhantom_LETopt_v2.m itself is NOT changed: a
% temporary copy with the test settings is written next to it and run.
%
% Must sit in the same folder as matRad_waterPhantom_LETopt_v2.m
% (LETd_optimization). Usage, in MATLAB from the matRadDWA_LET folder:
%   matRad_rc; addpath('LETd_optimization'); test_MCsquare_PC
%
% Success: the Command Window ends with "MCsquare TEST PASSED" and the
% line "MC: D_in/D_out ... D95 ... distal LETd ...".
% Output (tag _MCtest): DWA_final_machine/WaterPhantom_LETopt_v2/

here = fileparts(mfilename('fullpath'));
src  = fullfile(here, 'matRad_waterPhantom_LETopt_v2.m');
if ~exist(src, 'file'), error('matRad_waterPhantom_LETopt_v2.m not found next to this file (%s)', here); end
rootDir = here;
if ~exist(fullfile(rootDir, 'matRad_rc.m'), 'file'), rootDir = fileparts(rootDir); end

% ---- 1. is the MCsquare program there?
binDir = fullfile(rootDir, 'thirdParty', 'MCsquare', 'bin');
if ispc, exe = 'MCSquare_windows.exe'; elseif ismac, exe = 'MCsquare_mac'; else, exe = 'MCsquare_linux'; end
fprintf('\nMCsquare program: %s ... ', fullfile(binDir, exe));
if ~exist(fullfile(binDir, exe), 'file')
    error('not found. MCsquare cannot run on this computer.');
end
fprintf('found.\n');

% ---- 2. test settings in a temporary copy of the v2 script
txt = fileread(src);
txt = setLine(txt, 'machines', '{''DWA_proton_v2''}');
txt = setLine(txt, 'locations', '{''center''}');
txt = regexprep(txt, 'letWeights = struct\(''DWA_proton_v2'', \[[^\]]*\]', 'letWeights = struct(''DWA_proton_v2'', [0]', 'once');
txt = setLine(txt, 'mcRecalcFor', '0');
txt = setLine(txt, 'nHistFinal', '1e6');
txt = setLine(txt, 'mcOpt', 'false');
txt = setLine(txt, 'targetUnderdoseWeight', '10');
txt = setLine(txt, 'runTag', '''_MCtest''');
tmp = fullfile(here, 'tmp_MCsquare_test.m');
fid = fopen(tmp, 'w'); fwrite(fid, txt); fclose(fid);

% ---- 3. run it
t0 = tic;
try
    runVariant(tmp);
catch ME
    fprintf('\n#### MCsquare TEST FAILED after %.1f min:\n%s\n', toc(t0)/60, ME.message);
    for k = 1:min(5, numel(ME.stack))
        fprintf('   in %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    fprintf('Send this message to Claude. The overnight queue can still run (MCsquare is last).\n');
    return
end

% ---- 4. check that an MCsquare row was written
csv = fullfile(rootDir, 'DWA_final_machine', 'WaterPhantom_LETopt_v2', 'metrics_DWA_proton_v2_center_MCtest.csv');
ok = false;
if exist(csv, 'file')
    T = readtable(csv);
    r = T(strcmp(T.engine, 'MCrecalc'), :);
    if ~isempty(r)
        ok = true;
        p = T(strcmp(T.engine, 'PB'), :);
        fprintf('\n#### MCsquare TEST PASSED (%.1f min)\n', toc(t0)/60);
        fprintf('   pencil beam: D95 %.1f Gy, D_in/D_out %.1f, LETd behind target %.2f keV/um\n', p.D95(1), p.DinDout(1), p.LETd_distal_mean(1));
        fprintf('   MCsquare   : D95 %.1f Gy, D_in/D_out %.1f, LETd behind target %.2f keV/um\n', r.D95(1), r.DinDout(1), r.LETd_distal_mean(1));
        fprintf('   (1e6 histories: MC numbers are noisy; the overnight runs use 1e7)\n');
    end
end
if ~ok
    fprintf('\n#### MCsquare TEST FAILED: the run finished but no MCsquare result was written (%s).\n', csv);
end
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
