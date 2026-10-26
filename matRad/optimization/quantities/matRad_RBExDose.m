classdef matRad_RBExDose < matRad_DistributionQuantity
% matRad_RBExDose variable RBE-weighted dose from LQ effect:
% sqrt(effect/bx + gamma^2) - gamma, with gamma = ax / (2 bx)
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
        quantityName = 'RBExDose'
        requiredSubquantities = {'effect'}
    end

    methods
        function RBExDose = computeQuantity(this, dij, w, scen)
            effect = this.subQuantities.effect.getResult(dij, w, scen);
            ctScen = this.getCtScenario(dij, scen);
            ix = dij.ixDose{ctScen};
            gamma = dij.gamma{ctScen}(ix);

            RBExDose = zeros(dij.doseGrid.numOfVoxels, 1);
            RBExDose(ix) = sqrt(effect(ix) ./ dij.bx{ctScen}(ix) + gamma.^2) - gamma;
        end

        function wGrad = projectGradient(this, dij, fGrad, w, scen)
            RBExDose = this.getResult(dij, w, scen);
            ctScen = this.getCtScenario(dij, scen);
            ix = dij.ixDose{ctScen};

            dRBExDose_dEffect = zeros(dij.doseGrid.numOfVoxels, 1);
            dRBExDose_dEffect(ix) = 1 ./ (2 * dij.bx{ctScen}(ix) .* (RBExDose(ix) + dij.gamma{ctScen}(ix)));

            wGrad = this.subQuantities.effect.projectGradient(dij, fGrad .* dRBExDose_dEffect, w, scen);
        end

        function dij = prepareDij(~, dij, cst)
            dij = matRad_prepareBiologicalOptimizationDij(dij, cst, true);
        end
    end
end
