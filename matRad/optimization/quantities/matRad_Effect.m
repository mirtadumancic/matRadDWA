classdef matRad_Effect < matRad_DistributionQuantity
% matRad_Effect LQ effect: alphaDose + sqrtBetaDose.^2
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        quantityName = 'effect'
        requiredSubquantities = {'alphaDose', 'sqrtBetaDose'}
    end

    methods
        function e = computeQuantity(this, dij, w, scen)
            alphaDose = this.subQuantities.alphaDose.getResult(dij, w, scen);
            sqrtBetaDose = this.subQuantities.sqrtBetaDose.getResult(dij, w, scen);
            e = alphaDose + sqrtBetaDose.^2;
        end

        function wGrad = projectGradient(this, dij, fGrad, w, scen)
            sqrtBetaDose = this.subQuantities.sqrtBetaDose.getResult(dij, w, scen);
            wGrad = this.subQuantities.alphaDose.projectGradient(dij, fGrad, w, scen) + ...
                    this.subQuantities.sqrtBetaDose.projectGradient(dij, 2 * fGrad .* sqrtBetaDose, w, scen);
        end

        function dij = prepareDij(~, dij, cst)
            dij = matRad_prepareBiologicalOptimizationDij(dij, cst);
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc, alphaX, betaX)
            doses = optiFunc.getDoseParameters();
            optiFunc = optiFunc.setDoseParameters(alphaX * doses + betaX * doses.^2);
        end
    end
end
