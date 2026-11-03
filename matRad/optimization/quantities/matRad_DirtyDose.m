classdef matRad_DirtyDose < matRad_DijQuantity
% matRad_DirtyDose Dirty dose: the part of the physical dose delivered by
% bixel contributions whose LET exceeds a threshold k (Seckler et al.,
% Med Phys 2026, Eq. 4; concept from Heuchel et al.):
%
%   D^dirty_ij = D_ij   if LET_ij > k
%              = 0      otherwise
%
% With L_ij = LET_ij * D_ij stored in dij.mLETDose, LET_ij > k is evaluated
% as L_ij > k * D_ij, so no division is needed. The quantity stays linear in
% w (d = dij.mDirtyDose{scen} * w), so it is as cheap as physical dose.
%
% Threshold: dij.dirtyDoseLETThreshold [keV/um], default 2 keV/um.
% Use in cst:  obj = DoseObjectives.matRad_SquaredOverdosing(100,0);
%              obj.quantity = 'dirtyDose';
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
        quantityName = 'dirtyDose'
        requiredSubquantities = {}
        dijField = 'mDirtyDose'
    end

    methods
        function dij = prepareDij(~, dij, ~)
            dij = matRad_DirtyDose.addDirtyDoseDij(dij);
        end
    end

    methods (Static)
        function dij = addDirtyDoseDij(dij, threshold)
            % Adds dij.mDirtyDose (one matrix per scenario) if missing or
            % if it was built with a different threshold.
            if nargin < 2 || isempty(threshold)
                if isfield(dij, 'dirtyDoseLETThreshold') && ~isempty(dij.dirtyDoseLETThreshold)
                    threshold = dij.dirtyDoseLETThreshold;
                else
                    threshold = 2;
                end
            end

            if isfield(dij, 'mDirtyDose') && isfield(dij, 'mDirtyDoseThreshold') ...
                    && isequal(dij.mDirtyDoseThreshold, threshold)
                return;
            end

            if ~isfield(dij, 'mLETDose')
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Dirty dose needs dij.mLETDose. Calculate the dij with LET (engine.calcLET = true).');
            end

            dij.mDirtyDose = cell(size(dij.physicalDose));
            for s = 1:numel(dij.physicalDose)
                D = dij.physicalDose{s};
                if isempty(D) || s > numel(dij.mLETDose) || isempty(dij.mLETDose{s})
                    continue;
                end
                L = dij.mLETDose{s};
                % LET_ij > k  <=>  L_ij > k * D_ij  (for D_ij > 0)
                highLET = (L - threshold .* D) > 0;
                dij.mDirtyDose{s} = D .* highLET;
            end
            dij.mDirtyDoseThreshold = threshold;
            dij.dirtyDoseLETThreshold = threshold;
        end
    end
end
