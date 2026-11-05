function c = cns_okabeIto(name)
% cns_okabeIto  Colourblind-safe Okabe-Ito colours (RGB 0-1).
%
%   c = cns_okabeIto()          % 8x3 palette
%   c = cns_okabeIto('DWA_proton_v2')   % fixed colour per machine / plan type
%
% Fixed assignments keep every figure of the study consistent:
%   photon (clinical)  black      DWA_proton_v2  blue      Generic  vermillion
%   IBA_cyclotron      bluish green   extra: orange, sky blue, reddish purple, yellow
P = [0 0 0; 0 114 178; 213 94 0; 0 158 115; 230 159 0; 86 180 233; 204 121 167; 240 228 66] / 255;
if nargin == 0, c = P; return; end
switch name
    case {'photon', 'clinical'},        c = P(1, :);
    case {'DWA_proton_v2', 'DWA'},      c = P(2, :);
    case {'Generic'},                   c = P(3, :);
    case {'IBA_cyclotron', 'IBA'},      c = P(4, :);
    otherwise,                          c = P(5, :);
end
end
