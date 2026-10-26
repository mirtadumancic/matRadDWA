classdef matRad_EffectProjection < matRad_BackProjection
% matRad_EffectProjection Deprecated. Use matRad_BackProjection with defaultQuantity 'effect'.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019-2026 the matRad development team.
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function obj = matRad_EffectProjection()
            obj.defaultQuantity = 'effect';
            obj.warnDeprecatedProjection();
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc, alphaX, betaX)
            optiFunc = matRad_Effect.setBiologicalDosePrescriptions(optiFunc, alphaX, betaX);
        end
    end
end
