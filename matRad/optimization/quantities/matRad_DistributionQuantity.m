classdef (Abstract) matRad_DistributionQuantity < matRad_OptimizationQuantity
% matRad_DistributionQuantity Superclass for voxel-wise quantities.
% computeQuantity returns one value per dose grid voxel.
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

    methods (Static, Access = protected)
        function ctScen = getCtScenario(dij, scen)
            [ctScen, ~, ~] = ind2sub(size(dij.physicalDose), scen);
        end

        function RBE = getConstantRBE(dij, scen)
            % dij.RBE if valid (scalar or voxel-wise), 1 otherwise
            if isfield(dij, 'RBE') && (isscalar(dij.RBE) || numel(dij.RBE) == size(dij.physicalDose{scen}, 1)) && all(isfinite(dij.RBE))
                RBE = dij.RBE;
            else
                RBE = 1;
            end
        end
    end
end
