function dataClass = cns_checkPrivacy(caseCfg, caseFile)
% cns_checkPrivacy  Refuse to run a patient case from a folder an AI tool can reach.
%
%   dataClass = cns_checkPrivacy(caseCfg, caseFile)
%
% Same rules as letAgent (README, "Privacy protocol"):
%   * caseCfg.alias must be a study code (CNS_03), never a name, MRN or date.
%   * caseCfg.dataClass is 'patient', 'public' or 'synthetic'. A case read
%     from a .mat file or a DICOM folder counts as 'patient' unless
%     dataClass says otherwise AND dataClassNote says where the data is from.
%   * For a patient case, the case file, source.path / source.dir and
%     workDir must all be outside the code repository and outside every
%     folder listed in CNS_ESTRO/ai_connected_folders.txt (and in the
%     environment variable LETAGENT_AI_FOLDERS), workDir must be set, and
%     reb.protocol must hold the REB number.
%
% The check only looks at paths and settings. It never opens patient data.

alias = cns_getfield(caseCfg, 'alias', '');
checkAlias(alias);

srcType = cns_getfield(caseCfg, 'source.type', '');
dataClass = cns_getfield(caseCfg, 'dataClass', '');
if isempty(dataClass)
    if any(strcmp(srcType, {'matFile', 'dicom'})), dataClass = 'patient'; else, dataClass = 'synthetic'; end
end
if ~any(strcmp(dataClass, {'patient', 'public', 'synthetic'}))
    error('cns:privacy', 'dataClass must be patient, public or synthetic.');
end

if ~strcmp(dataClass, 'patient')
    if any(strcmp(srcType, {'matFile', 'dicom'})) && isempty(cns_getfield(caseCfg, 'dataClassNote', ''))
        error('cns:privacy', ['A case read from a .mat file or DICOM folder is treated as patient data ' ...
            'unless dataClassNote says where it comes from (e.g. "TCIA GLIS-RT, GLI_xxx").']);
    end
    return;
end

proto = strtrim(num2str(cns_getfield(caseCfg, 'reb.protocol', '')));
if isempty(proto) || ~isempty(strfind(upper(proto), 'XX'))
    error('cns:privacy', 'Patient case: set reb.protocol to your REB protocol number.');
end
workDir = cns_getfield(caseCfg, 'workDir', '');
if isempty(workDir)
    error('cns:privacy', 'Patient case: set workDir explicitly to a folder on hospital storage.');
end

roots = aiRoots();
paths = {caseFile, 'case file'; workDir, 'workDir'};
for f = {'source.path', 'source.dir', 'source.dicomDir'}
    p = cns_getfield(caseCfg, f{1}, '');
    if ~isempty(p), paths(end+1, :) = {p, f{1}}; end %#ok<AGROW>
end
for i = 1:size(paths, 1)
    for r = 1:numel(roots)
        if isInside(paths{i,1}, roots{r})
            error('cns:privacy', ['Patient case: the %s is inside the code repository or a folder ' ...
                'connected to an AI assistant (%s). Move it to hospital storage.'], paths{i,2}, roots{r});
        end
    end
end
end

% -------------------------------------------------------------------------
function checkAlias(alias)
if ~ischar(alias) || isempty(regexp(alias, '^[A-Za-z][A-Za-z0-9_-]{0,31}$', 'once'))
    error('cns:alias', ['Case alias must start with a letter and use only letters, digits, _ or - ' ...
        '(max 32 characters), e.g. CNS_03.']);
end
if ~isempty(regexp(alias, '\d{6,}', 'once')) || ...
        ~isempty(regexp(alias, '(19|20)\d{2}[-_]?\d{2}[-_]?\d{2}', 'once'))
    error('cns:alias', 'Case alias looks like an ID number or a date. Use a study code such as CNS_03.');
end
end

function roots = aiRoots()
roots = {};
repo = fileparts(which('matRad_rc'));
if ~isempty(repo), roots{end+1} = repo; end
here = fileparts(mfilename('fullpath'));
roots{end+1} = here;
f = fullfile(here, 'ai_connected_folders.txt');
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
