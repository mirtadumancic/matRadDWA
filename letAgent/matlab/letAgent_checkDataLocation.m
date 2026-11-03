function dataClass = letAgent_checkDataLocation(caseCfg, caseFile)
% letAgent_checkDataLocation  Refuse patient cases stored in AI-reachable folders.
%
% Mirrors privacy.check_data_location in the Python agent, so the guard also
% applies when letAgent_prepareCase / letAgent_runBatch are called directly
% from MATLAB. For a patient case (dataClass 'patient', which is the default
% for source.type 'matFile') the case file, source.path and workDir must all
% be outside the matRad repository and outside every folder listed in
% letAgent/config/ai_connected_folders.txt, and reb.protocol must be set.

dataClass = letAgent_getfield(caseCfg, 'dataClass', '');
srcType = letAgent_getfield(caseCfg, 'source.type', '');
if isempty(dataClass)
    if strcmp(srcType, 'matFile'), dataClass = 'patient'; else, dataClass = 'synthetic'; end
end
if ~any(strcmp(dataClass, {'patient', 'synthetic', 'public'}))
    error('letAgent:privacy', 'dataClass must be patient, synthetic or public.');
end
if ~strcmp(dataClass, 'patient')
    if strcmp(srcType, 'matFile') && isempty(letAgent_getfield(caseCfg, 'dataClassNote', ''))
        error('letAgent:privacy', ['A case loaded from a .mat file is treated as patient data unless ' ...
            'dataClassNote says where the data comes from.']);
    end
    return;
end

proto = strtrim(num2str(letAgent_getfield(caseCfg, 'reb.protocol', '')));
if isempty(proto) || ~isempty(strfind(proto, 'XX'))
    error('letAgent:privacy', 'Patient case: set reb.protocol to your REB protocol number.');
end
workDir = letAgent_getfield(caseCfg, 'workDir', '');
if isempty(workDir)
    error('letAgent:privacy', 'Patient case: set workDir explicitly to a folder on hospital storage.');
end

repoRoot = fileparts(which('matRad_rc'));
roots = [{repoRoot}, readAiFolders(repoRoot)];
paths = {caseFile, 'case file'; workDir, 'workDir'};
srcPath = letAgent_getfield(caseCfg, 'source.path', '');
if ~isempty(srcPath), paths(end+1, :) = {srcPath, 'source.path'}; end

for i = 1:size(paths, 1)
    for r = 1:numel(roots)
        if isInside(paths{i,1}, roots{r})
            error('letAgent:privacy', ['Patient case: the %s is inside the code repository or a folder ' ...
                'connected to an AI assistant. Move it to hospital storage outside these folders.'], paths{i,2});
        end
    end
end
end

function roots = readAiFolders(repoRoot)
roots = {};
f = fullfile(repoRoot, 'letAgent', 'config', 'ai_connected_folders.txt');
if exist(f, 'file')
    lines = strsplit(fileread(f), {'\r\n', '\n'});
    for i = 1:numel(lines)
        l = strtrim(regexprep(lines{i}, '#.*$', ''));
        if ~isempty(l), roots{end+1} = l; end %#ok<AGROW>
    end
end
env = getenv('LETAGENT_AI_FOLDERS');
if ~isempty(env)
    roots = [roots, strsplit(env, pathsep)];
end
end

function tf = isInside(p, root)
p = normPath(p); root = normPath(root);
tf = strcmp(p, root) || strncmp(p, [root '/'], numel(root) + 1);
end

function p = normPath(p)
p = strtrim(p);
if strncmp(p, '~', 1)
    home = getenv('HOME');
    if isempty(home), home = getenv('USERPROFILE'); end
    p = [home p(2:end)];
end
p = strrep(p, '\', '/');
isAbs = strncmp(p, '/', 1) || ~isempty(regexp(p, '^[A-Za-z]:/', 'once'));
if ~isAbs, p = [strrep(pwd, '\', '/') '/' p]; end
p = regexprep(p, '/+$', '');
p = regexprep(p, '/\./', '/');
if ispc, p = lower(p); end
end
