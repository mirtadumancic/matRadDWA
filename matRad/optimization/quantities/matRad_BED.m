classdef matRad_BED < matRad_DistributionQuantity
% matRad_BED biologically effective dose: effect / ax
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
        quantityName = 'BED'
        requiredSubquantities = {'effect'}
    end

    methods
        function BED = computeQuantity(this, dij, w, scen)
            effect = this.subQuantities.effect.getResult(dij, w, scen);
            BED = effect .* this.getInverseAlpha(dij, scen);
        end

        function wGrad = projectGradient(this, dij, fGrad, w, scen)
            wGrad = this.subQuantities.effect.projectGradient(dij, fGrad .* this.getInverseAlpha(dij, scen), w, scen);
        end

        function dij = prepareDij(~, dij, cst)
            dij = matRad_prepareBiologicalOptimizationDij(dij, cst);
        end
    end

    methods (Access = protected)
        function invAlpha = getInverseAlpha(this, dij, scen)
            ctScen = this.getCtScenario(dij, scen);
            ix = dij.ixDose{ctScen};
            invAlpha = zeros(dij.doseGrid.numOfVoxels, 1);
            invAlpha(ix) = 1 ./ dij.ax{ctScen}(ix);
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc, alphaX, betaX)
            doses = optiFunc.getDoseParameters();
            optiFunc = optiFunc.setDoseParameters(doses .* (1 + doses ./ (alphaX / betaX)));
        end
    end
end
