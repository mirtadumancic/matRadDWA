classdef (Abstract) matRad_OptimizationQuantity < handle
% matRad_OptimizationQuantity Superclass for quantities used in optimization.
% Results are cached per scenario and recomputed only when w changes.
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
        quantityName            % identifier used in cst objectives / constraints
        requiredSubquantities   % cell array of quantityNames this quantity depends on
    end

    properties
        subQuantities = struct()   % quantityName -> instance
    end

    properties (SetAccess = protected)
        dCache = {}
        wCache = []
    end

    methods
        function d = getResult(this, dij, w, scen)
            if ~isequal(this.wCache, w)
                this.dCache = {};
                this.wCache = w;
            end
            if scen > numel(this.dCache) || isempty(this.dCache{scen})
                this.dCache{scen} = this.computeQuantity(dij, w, scen);
            end
            d = this.dCache{scen};
        end

        function reset(this)
            this.dCache = {};
            this.wCache = [];
        end

        function dij = prepareDij(~, dij, ~)
            % Adds required fields to dij before optimization
        end

        % Should be abstract. Octave can't parse abstract methods.
        function d = computeQuantity(this, dij, w, scen) %#ok<STOUT,INUSD>
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('computeQuantity not implemented for %s!', class(this));
        end

        % fGrad: voxel gradient (nVox x k). Returns gradient wrt w (nBixel x k).
        function wGrad = projectGradient(this, dij, fGrad, w, scen) %#ok<STOUT,INUSD>
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('projectGradient not implemented for %s!', class(this));
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc, ~, ~)
        end

        function classList = getAvailableQuantities()
            matRad_cfg = MatRad_Config.instance();

            folders = [{fileparts(mfilename('fullpath'))} matRad_cfg.userfolders];

            persistent metaQuantities lastFolders

            if ~matRad_cfg.isOctave && ~isempty(metaQuantities) && ~all(cellfun(@isvalid, metaQuantities))
                metaQuantities = [];
            end

            if isempty(metaQuantities) || ~isequal(lastFolders, folders)
                lastFolders = folders;
                metaQuantities = matRad_findSubclasses(meta.class.fromName(mfilename('class')), 'folders', folders, 'includeSubfolders', true);
            end

            classList = matRad_identifyClassesByConstantProperties(metaQuantities, 'quantityName');
        end

        function quantity = create(quantityName)
            classList = matRad_OptimizationQuantity.getAvailableQuantities();
            ix = find(strcmp({classList.quantityName}, quantityName), 1);
            if isempty(ix)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Optimization quantity ''%s'' not available!', quantityName);
            end
            quantity = classList(ix).handle();
        end
    end
end
