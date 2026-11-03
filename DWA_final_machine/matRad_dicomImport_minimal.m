%% matRad_dicomImport_minimal.m
% Minimal DICOM import test for matRad (script version of the import GUI).
% matRad ships no sample DICOM, so step 1 exports a bundled phantom to DICOM
% (CT slices + RTSTRUCT), step 2 imports that folder back with
% matRad_DicomImporter, and step 3 checks that the CT and structures match.
% For a real patient: skip step 1 and set dicomDir to the folder holding all
% .dcm files (CT + RTSTRUCT [+ RTPLAN, RTDOSE]) in one flat folder.
%
% Run from the matRadDWA folder:  matRad_rc; run('DWA_final_machine/matRad_dicomImport_minimal.m')

matRad_cfg = matRad_rc;

phantom  = 'HEAD_AND_NECK';                         % or 'TG119' (smaller, faster)
%phantom  = 'TG119';  
dicomDir = fullfile(matRad_cfg.primaryUserFolder, 'patients', ['dicomTest_' phantom]);  % userdata: ignored by git

%% 1) Make sample DICOM from a bundled phantom (skip for real patient data)
orig = load([phantom '.mat']);                       % gives ct, cst
if exist(dicomDir, 'dir'), delete(fullfile(dicomDir, '*.dcm')); else, mkdir(dicomDir); end   % start clean
% Pass ct and cst explicitly: with no inputs the exporter reads ct/cst/pln/stf/resultGUI
% from the base workspace and would also export any leftover dose (RTDose).
dcmExp = matRad_DicomExporter(orig.ct, orig.cst);
dcmExp.dicomDir = dicomDir;
tic; dcmExp.matRad_exportDicom(); fprintf('Export: %.1f s, %d .dcm files in %s\n', toc, numel(dir(fullfile(dicomDir,'*.dcm'))), dicomDir);

%% 2) Import the DICOM folder (this is the part you use for patients)
tic;
imp = matRad_DicomImporter(dicomDir);                % scans the folder: CT, RTSTRUCT, RTPLAN, RTDOSE
fprintf('Found: %d CT slices, %d RTSTRUCT, %d RTPLAN, %d RTDOSE\n', size(imp.importFiles.ct,1), ...
    size(imp.importFiles.rtss,1), size(imp.importFiles.rtplan,1), size(imp.importFiles.rtdose,1));
imp.matRad_importDicom();                            % also puts ct, cst (pln, stf, resultGUI) in the workspace
ct = imp.ct; cst = imp.cst;
fprintf('Import: %.1f s\n', toc);

%% 3) Check the round trip
fprintf('\nCT  original %s @ %.2f/%.2f/%.2f mm | imported %s @ %.2f/%.2f/%.2f mm\n', ...
    mat2str(orig.ct.cubeDim), orig.ct.resolution.x, orig.ct.resolution.y, orig.ct.resolution.z, ...
    mat2str(ct.cubeDim), ct.resolution.x, ct.resolution.y, ct.resolution.z);
fprintf('%-25s %-8s %12s %12s\n', 'Structure', 'Type', 'orig voxels', 'imported');
for i = 1:size(orig.cst,1)
    j = find(strcmpi(cst(:,2), orig.cst{i,2}), 1);
    if isempty(j), nImp = NaN; else, nImp = numel(cst{j,4}{1}); end
    fprintf('%-25s %-8s %12d %12d\n', orig.cst{i,2}, orig.cst{i,3}, numel(orig.cst{i,4}{1}), nImp);
end

%% 4) Look at it and save as a matRad .mat
save(fullfile(dicomDir, [phantom '_imported.mat']), 'ct', 'cst');
matRadGUI;
