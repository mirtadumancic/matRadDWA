classdef matRad_AlphaDose < matRad_DijQuantity
% matRad_AlphaDose alpha x dose. Uses dij.mAlphaDose if available,
% otherwise ax .* RBE .* (physicalDose * w).
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
        quantityName = 'alphaDose'
        requiredSubquantities = {}
        dijField = 'mAlphaDose'
    end

    methods
        function d = computeQuantity(this, dij, w, scen)
            if isfield(dij, this.dijField)
                d = computeQuantity@matRad_DijQuantity(this, dij, w, scen);
            else
                d = this.getVoxelScaling(dij, scen) .* (dij.physicalDose{scen} * w);
            end
        end

        function wGrad = projectGradient(this, dij, fGrad, w, scen)
            if isfield(dij, this.dijField)
                wGrad = projectGradient@matRad_DijQuantity(this, dij, fGrad, w, scen);
            else
                wGrad = ((fGrad .* this.getVoxelScaling(dij, scen))' * dij.physicalDose{scen})';
            end
        end

        function dij = prepareDij(this, dij, cst)
            if ~isfield(dij, this.dijField)
                dij = matRad_prepareBiologicalOptimizationDij(dij, cst);
            end
        end
    end

    methods (Access = protected)
        function s = getVoxelScaling(this, dij, scen)
            s = dij.ax{this.getCtScenario(dij, scen)} .* this.getConstantRBE(dij, scen);
        end
    end
end
