function caseFile = cns_makeTestCase_v1(outDir)
% cns_makeTestCase_v1  Synthetic stand-in for the UHN CNS case, for rehearsing UHN_CNS_01/02/03.
%
%   caseFile = cns_makeTestCase_v1('~/Work/CNS_ESTRO_test')
%
% Writes, in outDir:
%   SYN_UHNlike.mat       ct, cst (UHN structure names, target behind the LEFT eye,
%                         no BODY) and refDose (analytic photon-like dose, 54 Gy)
%   case_SYN_UHNlike.json case file (dataClass synthetic) pointing at it
% Then set S.caseFile = caseFile in the three scripts. No patient data.
if nargin < 1 || isempty(outDir), outDir = fullfile(pwd, 'CNS_ESTRO_test'); end
outDir = cns_expandPath(outDir);
if ~exist(outDir, 'dir'), mkdir(outDir); end
[ct, cst, refDose] = cns_makeSyntheticUHNlike(struct('resolutionMm', 3));
matFile = fullfile(outDir, 'SYN_UHNlike.mat');
save(matFile, 'ct', 'cst', 'refDose');
here = fileparts(mfilename('fullpath'));
c = cns_readJson(fullfile(here, 'cases', 'UHN_CNS_ependymoma_TEMPLATE_v1.json'));
c.alias = 'SYN_UHNlike';
c.dataClass = 'synthetic';
c.dataClassNote = 'cns_makeSyntheticUHNlike (no patient data)';
c = rmfield(c, 'reb');
c.workDir = fullfile(outDir, 'runs');
c.source = struct('type', 'matFile', 'path', matFile);
caseFile = fullfile(outDir, 'case_SYN_UHNlike.json');
cns_writeJson(caseFile, c);
fprintf('Test case written: %s\nSet S.caseFile = ''%s'' in UHN_CNS_01/02/03.\n', caseFile, caseFile);
end
