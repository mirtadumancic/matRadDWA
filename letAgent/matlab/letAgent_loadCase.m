function [ct, cst] = letAgent_loadCase(caseCfg)
% letAgent_loadCase  Load ct and cst for a case.
%   source.type = 'matRadPhantom'  source.name = 'TG119.mat'  (ships with matRad)
%   source.type = 'matFile'        source.path = 'D:\local\CNS_03.mat'
%                 (a matRad .mat with ct and cst, e.g. from matRad_importDicom;
%                  keep it on the hospital machine)
%   source.type = 'synthetic'      source.name = 'cnsStress', other fields
%                 are passed to letAgent_makeCnsPhantom
src = caseCfg.source;
switch src.type
    case 'matRadPhantom'
        S = load(src.name, 'ct', 'cst');
    case 'matFile'
        S = load(src.path, 'ct', 'cst');
    case 'synthetic'
        switch src.name
            case 'cnsStress'
                [S.ct, S.cst] = letAgent_makeCnsPhantom(src);
            otherwise
                error('letAgent:source', 'Unknown synthetic case %s', src.name);
        end
    otherwise
        error('letAgent:source', 'Unknown source type %s', src.type);
end
ct = S.ct;
cst = S.cst;
end
