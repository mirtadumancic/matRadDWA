%% run_UHN_CNS_07_patients_v1.m
% Runs UHN_CNS_07 (robust dose-only DWA_proton_v4 plan + clinical goals + MCsquare)
% for the two UHN CNS patients, one after the other (or one only: runOnly).
% Same as run_UHN_CNS_07_overnight_v1.m, plus runOnly and its own log/summary names.
% A failure in one patient does not stop the other. Everything is logged.
%
% RUN ON THE UHN PC (X: is created automatically if missing), e.g. from M:\My Documents\matRadDWA_LET:
%   cd X:\matRadDWA_LET; run('CNS_ESTRO/run_UHN_CNS_07_patients_v1.m')
% Afterwards: open the log (path printed at the start) and the summary file
%   UHN_CNS_07_patients_summary_<date>.txt in the same folder.

dataRoot = 'C:/Users/t146774uhn/Documents/UHN_Test_Patient_Data';   % hospital storage, outside the repo
patients = struct( ...
    'caseFile',     {fullfile(dataRoot, 'Patient1', 'CNS_01.json'), fullfile(dataRoot, 'Patient2', 'CNS_02.json')}, ...
    'gantryAngles', {[270 100 75], [270 290 60]}, ...
    'couchAngles',  {[0 0 300],    [0 55 270]}, ...
    'label',        {'CNS_01', 'CNS_02'});
runMC = true;
runOnly = {};            % e.g. {'CNS_01'} to run one patient only; {} = both
if ~isempty(runOnly), patients = patients(ismember({patients.label}, runOnly)); end

stamp = datestr(now, 'yyyy-mm-dd_HHMM');
logFile = fullfile(dataRoot, sprintf('UHN_CNS_07_patients_log_%s.txt', stamp));
sumFile = fullfile(dataRoot, sprintf('UHN_CNS_07_patients_summary_%s.txt', stamp));
diary(logFile); diary on
fprintf('Run started %s\nLog: %s\n', datestr(now), logFile);
here = fileparts(mfilename('fullpath'));
% (Nov 9) MCsquare fails on paths with spaces ("M:\My Documents"): work from the X: drive
% (subst X: "M:\My Documents"). subst is lost at every restart, so set it again if needed.
xRoot = 'X:\matRadDWA_LET';
if ispc && ~exist(xRoot, 'dir')
    [st, msg] = system('subst X: "M:\My Documents"');
    if st ~= 0 || ~exist(xRoot, 'dir')
        error('Could not create X: with subst (%s): %s', 'subst X: "M:\My Documents"', msg);
    end
    fprintf('%s\n', 'X: created with subst X: "M:\My Documents"');
end
if ispc, here = fullfile(xRoot, 'CNS_ESTRO'); end
script = fullfile(here, 'UHN_CNS_07_DWAv4_robustDoseOnly_clinicalGoals_v1.m');
summary = {};
for iPat = 1:numel(patients)
    pat = patients(iPat);
    fprintf('\n\n================ %s  (%s) ================\n', pat.label, datestr(now));
    t0 = tic;
    S_override = struct('caseFile', pat.caseFile, 'gantryAngles', pat.gantryAngles, ...
        'couchAngles', pat.couchAngles, 'setLabel', ['RS_approved_' pat.label], ...
        'geometryOnly', false, 'runMC', runMC);
    try
        run(script);
        msg = sprintf('%s: OK in %.0f min. PB nominal %d/%d, PB worst %d/%d goals met', pat.label, toc(t0)/60, ...
            sum(Tn.pass, 'omitnan'), nnz(~isnan(Tn.pass)), sum(Tw.pass, 'omitnan'), nnz(~isnan(Tw.pass)));
        if exist('Tm', 'var') && ~isempty(Tm)
            msg = sprintf('%s, MCsquare %d/%d', msg, sum(Tm.pass, 'omitnan'), nnz(~isnan(Tm.pass)));
        end
        msg = sprintf('%s. Table: %s', msg, outDir);
    catch err
        msg = sprintf('%s: FAILED after %.0f min: %s (line %d of %s)', pat.label, toc(t0)/60, err.message, ...
            err.stack(1).line, err.stack(1).name);
    end
    fprintf('\n>>> %s\n', msg);
    summary{end+1} = msg; %#ok<SAGROW>
    fid = fopen(sumFile, 'w'); fprintf(fid, '%s\n', summary{:}); fclose(fid);
    clearvars -except patients runMC runOnly stamp logFile sumFile here script summary iPat dataRoot
    close all force
end
fprintf('\nAll done %s\n', datestr(now));
fprintf('%s\n', summary{:});
diary off
