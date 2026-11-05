function [ct, cst, refDose] = cns_importDicom(dicomDir, outMat)
% cns_importDicom  Import CT + RTSTRUCT (+ RTPLAN, RTDOSE) from a DICOM folder.
%
%   [ct, cst, refDose] = cns_importDicom(dicomDir)            % in memory only
%   [ct, cst, refDose] = cns_importDicom(dicomDir, outMat)    % also save ct, cst, refDose
%
% Uses matRad_DicomImporter (same calls as matRad_dicomImport_minimal.m).
% Put all .dcm files of one patient in one folder; series in subfolders
% (TCIA layout) are collected too. MR images and REG files are ignored:
% planning uses the CT, and the RT objects are on the CT frame.
%
% refDose is the clinical dose from the RTDOSE (e.g. the delivered photon
% plan), interpolated to the CT grid by matRad, or [] if there is none:
%   refDose.cube             dose cube on the CT grid, as matRad imported it
%   refDose.toTotal          factor giving the total dose (NaN: multiply by
%                            the number of fractions of the case file)
%   refDose.summationType    DICOM DoseSummationType (PLAN, FRACTION, ...)
%   refDose.doseType         'physicalDose' or 'RBExDose'
%   refDose.numOfFractions   from the RTPLAN, if it was imported
% matRad divides a PLAN dose by the RTPLAN fraction number when it imports
% the RTPLAN and leaves it as the total otherwise; toTotal undoes either.
%
% PRIVACY. Use this on public data on any computer. For UHN patients run it
% ONLY on the hospital machine, with dicomDir and outMat on hospital storage.

if ~exist(dicomDir, 'dir')
    error('cns:dicom', 'DICOM folder not found: %s', dicomDir);
end
flatDir = dicomDir;
if isempty(dir(fullfile(dicomDir, '*.dcm')))
    list = findDcm(dicomDir);
    if isempty(list), error('cns:dicom', 'No .dcm files under %s', dicomDir); end
    flatDir = tempname();
    mkdir(flatDir);
    for i = 1:numel(list)
        [~, n, e] = fileparts(list{i});
        copyfile(list{i}, fullfile(flatDir, sprintf('%05d_%s%s', i, n, e)));
    end
    cleanup = onCleanup(@() rmdir(flatDir, 's')); %#ok<NASGU>
end

imp = matRad_DicomImporter(flatDir);
nCt = size(imp.importFiles.ct, 1); nRs = size(imp.importFiles.rtss, 1);
nRp = size(imp.importFiles.rtplan, 1); nRd = size(imp.importFiles.rtdose, 1);
fprintf('[cns] DICOM: %d CT slices, %d RTSTRUCT, %d RTPLAN, %d RTDOSE\n', nCt, nRs, nRp, nRd);
if nRs == 0, error('cns:dicom', 'No RTSTRUCT found in %s', dicomDir); end
if nRs > 1
    warning('cns:dicom', 'Several RTSTRUCTs found; matRad uses the first one.');
end

try
    imp.matRad_importDicom();
catch err
    if nRp == 0, rethrow(err); end
    % Clinical photon plans (VMAT, many control points) can fail to import as
    % matRad steering. The plan is not needed here, only the dose.
    warning('cns:dicom', 'RTPLAN import failed (%s); importing without it.', err.message);
    imp = matRad_DicomImporter(flatDir);
    imp.importFiles.rtplan = cell(0, size(imp.importFiles.rtplan, 2));
    imp.matRad_importDicom();
end
ct = imp.ct;
cst = imp.cst;
fprintf('[cns] imported CT %s at %.2f x %.2f x %.2f mm, %d structures:\n', mat2str(ct.cubeDim), ...
    ct.resolution.x, ct.resolution.y, ct.resolution.z, size(cst, 1));
fprintf('      %s\n', strjoin(cst(:,2)', ' | '));

refDose = pickDose(imp);

if nargin > 1 && ~isempty(outMat)
    save(outMat, 'ct', 'cst', 'refDose', '-v7');
    fprintf('[cns] saved %s\n', outMat);
end
end

% -------------------------------------------------------------------------
function refDose = pickDose(imp)
refDose = [];
R = imp.resultGUI;
if isempty(R) || ~isstruct(R), return; end
fn = fieldnames(R);
isCube = cellfun(@(f) isnumeric(R.(f)) && ndims(R.(f)) == 3, fn);
fn = fn(isCube);
cand = fn(~cellfun(@isempty, regexp(fn, '^(physicalDose|RBExDose)', 'once')));
if isempty(cand), return; end
whole = cand(cellfun(@isempty, regexp(cand, '_beam\d+', 'once')));
if ~isempty(whole)
    use = whole(1);
    cube = double(R.(use{1}));
else
    use = cand;                          % only per-beam doses: add them up
    cube = zeros(size(R.(cand{1})));
    for i = 1:numel(cand), cube = cube + double(R.(cand{i})); end
end
meta = struct();
try, meta = R.doseMetaInfo.(use{1}); catch, end %#ok<NOCOM>
sumType = upper(cns_getfield(meta, 'DoseSummationType', 'PLAN'));
nFxPlan = NaN;
if ~isempty(imp.pln) && isstruct(imp.pln) && isfield(imp.pln, 'numOfFractions')
    nFxPlan = imp.pln.numOfFractions;
end
if any(strcmp(sumType, {'PLAN', 'BEAM'}))
    if isnan(nFxPlan), toTotal = 1; else, toTotal = nFxPlan; end   % matRad divided by nFxPlan
else
    toTotal = nFxPlan;                    % per-fraction dose; NaN -> case file fractions
end
refDose = struct('cube', cube, 'toTotal', toTotal, 'summationType', sumType, ...
    'doseType', regexprep(use{1}, '_.*$', ''), 'numOfFractions', nFxPlan, 'sourceFields', {use});
if numel(whole) > 1
    warning('cns:dicom', 'Several dose cubes (%s); using %s.', strjoin(whole', ', '), use{1});
end
fprintf('[cns] reference dose: %s (%s), max %.2f Gy as imported, factor to total %g\n', ...
    use{1}, sumType, max(cube(:)), toTotal);
end

function list = findDcm(d)
list = {};
e = dir(d);
for i = 1:numel(e)
    if any(strcmp(e(i).name, {'.', '..'})), continue; end
    p = fullfile(d, e(i).name);
    if e(i).isdir
        list = [list, findDcm(p)]; %#ok<AGROW>
    elseif ~isempty(regexpi(e(i).name, '\.dcm$', 'once'))
        list{end+1} = p; %#ok<AGROW>
    end
end
end
