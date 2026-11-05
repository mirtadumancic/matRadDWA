function caseFile = cns_newCaseFile_v1(dicomDir, rebProtocol, alias)
% cns_newCaseFile_v1  Write the case file for a patient, with the paths filled in.
%
%   caseFile = cns_newCaseFile_v1()                                  % dialogs
%   caseFile = cns_newCaseFile_v1('/Users/.../CNS_01/DICOM', 'REB-26-1234', 'CNS_01')
%
% dicomDir : folder with CT + RS + RP + RD only (keep MR and REG elsewhere)
% Writes <parent of dicomDir>/<alias>.json from cases/UHN_CNS_ependymoma_TEMPLATE_v1.json:
%   source.dicomDir = dicomDir
%   source.path     = <parent>/<alias>_imported.mat   (written by UHN_CNS_01)
%   workDir         = <parent>/CNS_ESTRO_runs
% Keep the parent folder OUT of any folder connected to an AI tool: the
% scripts refuse to run a patient case from the repository or from the
% folders listed in CNS_ESTRO/ai_connected_folders.txt.
if nargin < 1 || isempty(dicomDir)
    dicomDir = uigetdir(pwd, 'Patient DICOM folder (CT + RS + RP + RD only)');
    if isequal(dicomDir, 0), error('No folder selected.'); end
end
if nargin < 2 || isempty(rebProtocol)
    a = inputdlg({'REB protocol number', 'Study code (alias), e.g. CNS_01'}, 'Case file', 1, {'', 'CNS_01'});
    if isempty(a), error('Cancelled.'); end
    rebProtocol = strtrim(a{1}); alias = strtrim(a{2});
end
if nargin < 3 || isempty(alias), alias = 'CNS_01'; end
dicomDir = regexprep(dicomDir, '[\\/]+$', '');
root = fileparts(dicomDir);
here = fileparts(mfilename('fullpath'));
c = cns_readJson(fullfile(here, 'cases', 'UHN_CNS_ependymoma_TEMPLATE_v1.json'));
c.alias = alias;
c.reb.protocol = rebProtocol;
c.workDir = fullfile(root, 'CNS_ESTRO_runs');
c.source.type = 'matFile';
c.source.path = fullfile(root, [alias '_imported.mat']);
c.source.dicomDir = dicomDir;
caseFile = fullfile(root, [alias '.json']);
if exist(caseFile, 'file')
    error('%s exists already; edit it or delete it first.', caseFile);
end
cns_writeJson(caseFile, c);
fprintf(['Case file written: %s\nCheck prescription (54 Gy, 30 fx?), structures.map and beams, ' ...
    'then set S.caseFile to it (or pick it when the script asks).\n'], caseFile);
end
