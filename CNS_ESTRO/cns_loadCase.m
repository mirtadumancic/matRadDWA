function [ct, cst, refDose] = cns_loadCase(caseCfg)
% cns_loadCase  Load ct, cst and (if available) the clinical reference dose.
%
%   source.type = 'matFile'        source.path = '<file>.mat' with ct, cst
%                                  [, refDose] (from cns_importDicom)
%   source.type = 'dicom'          source.dir  = folder with CT + RTSTRUCT
%                                  [+ RTPLAN + RTDOSE], imported in memory
%   source.type = 'matRadPhantom'  source.name = 'HEAD_AND_NECK.mat'
%   source.type = 'synthetic'      other fields go to cns_makeSyntheticCNS
%
% refDose is [] when there is no RTDOSE, or when caseCfg.reference.use is false.
src = caseCfg.source;
refDose = [];
switch src.type
    case 'matFile'
        S = load(cns_expandPath(src.path));
        if isfield(S, 'refDose'), refDose = S.refDose; end
    case 'dicom'
        [S.ct, S.cst, refDose] = cns_importDicom(cns_expandPath(src.dir));
    case 'matRadPhantom'
        S = load(src.name, 'ct', 'cst');
    case 'synthetic'
        [S.ct, S.cst] = cns_makeSyntheticCNS(src);
    otherwise
        error('cns:source', 'Unknown source.type %s', src.type);
end
ct = S.ct;
cst = S.cst;
if ~logical(cns_getfield(caseCfg, 'reference.use', true))
    refDose = [];
end
end
