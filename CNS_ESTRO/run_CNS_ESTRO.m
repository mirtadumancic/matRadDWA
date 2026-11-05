%% run_CNS_ESTRO.m  CNS test-patient planning for the ESTRO 2027 abstract
%
% Plans each case with every machine x plan variant in study_ESTRO2027.json
% (DWA_proton_v2 and Generic; dose-only and LET x dose / dirty-dose plans),
% evaluates target DVH and OAR dose-LET metrics, and builds the cohort table.
%
% Run from the matRadDWA root on the let_opt branch:
%   matRad_rc; addpath('CNS_ESTRO'); run('CNS_ESTRO/run_CNS_ESTRO.m')
%
% WHERE TO RUN WHAT
%   Synthetic / public (TCIA) cases: any computer, including the Mac.
%   UHN patients: ONLY on the UHN machine, with case files, data and
%   workDir on hospital storage (cns_checkPrivacy enforces this). Never on a
%   computer that has an AI tool connected to the data folders.

matRad_rc;
here = fileparts(mfilename('fullpath'));
addpath(here);
studyFile = fullfile(here, 'study_ESTRO2027.json');

%% 1) Choose the cases  ---------------------------------------------------
% (a) smoke test, no patient data (about 10-20 min at 3 mm)
caseFiles = {fullfile(here, 'cases', 'CNS_synthetic.json')};

% (a2) matRad HEAD_AND_NECK rehearsal with a photon reference: run
%      matRad_example5_DWA_HN_training.m (sections 1-2 make the .mat), then
%   caseFiles = {fullfile(here, 'cases', 'HN_matRad_training.json')};

% (b) public TCIA case: import once, then plan
%   cns_importDicom('~/Work/TCIA/GLIS-RT/GLI_xxx_AAC', '~/Work/CNS_ESTRO_public/GLI_xxx_AAC.mat');
%   caseFiles = {'~/Work/CNS_ESTRO_public/GLI_xxx_AAC.json'};   % from cases/CNS_public_GLISRT_TEMPLATE.json

% (c) UHN patients, on the UHN machine only (case files on hospital storage)
%   cns_importDicom('D:/CNS_DICOM/CNS_01', 'D:/CNS_ESTRO_cases/CNS_01.mat');
%   caseFiles = {'D:/CNS_ESTRO_cases/CNS_01.json', 'D:/CNS_ESTRO_cases/CNS_02.json'};

%% 2) Plan and evaluate every case  ---------------------------------------
caseDirs = {};
for i = 1:numel(caseFiles)
    out = cns_runCase(cns_expandPath(caseFiles{i}), studyFile);
    caseDirs{end+1} = out.caseDir; %#ok<SAGROW>
end

%% 3) Cohort tables (median [range], paired differences)  -----------------
cohortDir = fullfile(fileparts(caseDirs{1}), 'cohort');
cns_cohortSummary(caseDirs, studyFile, cohortDir);
type(fullfile(cohortDir, 'summary_abstract.md'));
