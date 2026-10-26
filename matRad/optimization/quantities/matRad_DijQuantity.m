classdef (Abstract) matRad_DijQuantity < matRad_DistributionQuantity
% matRad_DijQuantity Voxel-wise quantity linear in w: d = dij.(dijField){scen} * w
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

    properties (Abstract, Constant)
        dijField
    end

    methods
        function d = computeQuantity(this, dij, w, scen)
            D = this.getInfluenceMatrix(dij, scen);
            if isempty(D)
                d = [];
            else
                d = D * w;
            end
        end

        function wGrad = projectGradient(this, dij, fGrad, ~, scen)
            D = this.getInfluenceMatrix(dij, scen);
            if isempty(D)
                wGrad = [];
            else
                wGrad = (fGrad' * D)';
            end
        end
    end

    methods (Access = protected)
        function D = getInfluenceMatrix(this, dij, scen)
            if ~isfield(dij, this.dijField)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('dij.%s required for quantity ''%s'' not found!', this.dijField, this.quantityName);
            end
            D = dij.(this.dijField){scen};
            if isempty(D)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario %d in dij.%s detected! This should not happen...\n', scen, this.dijField);
            end
        end
    end
end
