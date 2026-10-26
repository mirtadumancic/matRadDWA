classdef matRad_ConstantRBExDose < matRad_DistributionQuantity
% matRad_ConstantRBExDose RBE-weighted dose with constant RBE: dij.RBE .* physicalDose
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
        quantityName = 'constantRBExDose'
        requiredSubquantities = {'physicalDose'}
    end

    methods
        function RBExDose = computeQuantity(this, dij, w, scen)
            RBExDose = dij.RBE .* this.subQuantities.physicalDose.getResult(dij, w, scen);
        end

        function wGrad = projectGradient(this, dij, fGrad, w, scen)
            wGrad = this.subQuantities.physicalDose.projectGradient(dij, dij.RBE .* fGrad, w, scen);
        end

        function dij = prepareDij(~, dij, ~)
            if ~isfield(dij, 'RBE')
                dij.RBE = 1.1;
            end
        end
    end
end
