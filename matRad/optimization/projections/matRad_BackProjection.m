classdef matRad_BackProjection < handle
% matRad_BackProjection manages the optimization quantities found in the
% cst and the scenarios used during optimization
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

    properties (SetAccess = protected)
        quantities = struct()           % quantityName -> instance, including subquantities
        optimizationQuantities = {}     % quantityNames used by objectives
        constraintQuantities = {}       % quantityNames used by constraints
        d                               % d.(quantityName){scen}
        wGrad                           % objective gradient wrt w, summed over scenarios
    end

    properties
        defaultQuantity = 'physicalDose' % used for objectives / constraints without quantity
        scenarios    = 1        %Scenario indices to evaluate (used for 4D & robust/stochastic optimization)
        scenarioProb = 1        %Probability associated with scenario (for stochastic optimization)
        nominalCtScenarios = 1; %nominal ct scenario (no shift, no range error) indices to evaluate (used for 4D & robust/stochastic optimization, when at least one cst structure does not have robustness)
    end

    methods
        function obj = matRad_BackProjection()
            obj.d = struct();
            obj.wGrad = [];
        end

        function instantiateQuantities(this, cst)
            this.quantities = struct();
            optNames = {};
            constrNames = {};
            for i = 1:size(cst, 1)
                for j = 1:numel(cst{i, 6})
                    optiFunc = cst{i, 6}{j};
                    if isa(optiFunc, 'DoseObjectives.matRad_DoseObjective')
                        optNames{end + 1} = this.getQuantityName(optiFunc); %#ok<AGROW>
                    elseif isa(optiFunc, 'DoseConstraints.matRad_DoseConstraint')
                        constrNames{end + 1} = this.getQuantityName(optiFunc); %#ok<AGROW>
                    end
                end
            end
            this.optimizationQuantities = unique(optNames);
            this.constraintQuantities = unique(constrNames);

            cellfun(@(name) this.addQuantity(name), [this.optimizationQuantities this.constraintQuantities], 'UniformOutput', false);
        end

        function quantity = addQuantity(this, quantityName)
            % Instantiates a quantity and its required subquantities
            if ~isfield(this.quantities, quantityName)
                quantity = matRad_OptimizationQuantity.create(quantityName);
                if ~isa(quantity, 'matRad_DistributionQuantity')
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Only distribution quantities are supported, ''%s'' is not!', quantityName);
                end
                for sub = quantity.requiredSubquantities
                    quantity.subQuantities.(sub{1}) = this.addQuantity(sub{1});
                end
                this.quantities.(quantityName) = quantity;
            end
            quantity = this.quantities.(quantityName);
        end

        function quantityName = getQuantityName(this, optiFunc)
            quantityName = optiFunc.quantity;
            if isempty(quantityName)
                quantityName = this.defaultQuantity;
            end
        end

        function quantity = getQuantity(this, quantityNameOrOptiFunc)
            if ischar(quantityNameOrOptiFunc)
                quantityName = quantityNameOrOptiFunc;
            else
                quantityName = this.getQuantityName(quantityNameOrOptiFunc);
            end
            quantity = this.addQuantity(quantityName);
        end

        function dij = prepareDij(this, dij, cst)
            for name = fieldnames(this.quantities)'
                dij = this.quantities.(name{1}).prepareDij(dij, cst);
            end
        end

        function obj = compute(obj, dij, w)
            obj.d = struct();
            for name = unique([obj.optimizationQuantities obj.constraintQuantities])
                quantity = obj.quantities.(name{1});
                dQuantity = cell(size(dij.physicalDose));
                for scen = obj.scenarios
                    dQuantity{scen} = quantity.getResult(dij, w, scen);
                end
                obj.d.(name{1}) = dQuantity;
            end
        end

        function obj = computeGradient(obj, dij, fGrad, w)
            % fGrad.(quantityName){scen}: objective gradient wrt the quantity
            obj.wGrad = zeros(dij.totalNumOfBixels, 1);
            for name = fieldnames(fGrad)'
                quantity = obj.quantities.(name{1});
                for scen = obj.scenarios
                    if scen <= numel(fGrad.(name{1})) && ~isempty(fGrad.(name{1}){scen})
                        obj.wGrad = obj.wGrad + quantity.projectGradient(dij, fGrad.(name{1}){scen}, w, scen);
                    end
                end
            end
        end

        function d = GetResult(obj)
            d = obj.d;
        end

        function wGrad = GetGradient(obj)
            wGrad = obj.wGrad;
        end

        % Deprecated single-scenario interface, evaluates defaultQuantity
        function d = computeSingleScenario(obj, dij, scen, w)
            d = obj.getQuantity(obj.defaultQuantity).computeQuantity(dij, w, scen);
        end

        function wGrad = projectSingleScenarioGradient(obj, dij, doseGrad, scen, w)
            wGrad = obj.getQuantity(obj.defaultQuantity).projectGradient(dij, doseGrad{scen}, w, scen);
        end
    end

    methods (Access = protected)
        function warnDeprecatedProjection(obj)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('%s is deprecated. Use matRad_BackProjection with defaultQuantity ''%s'' or set the quantity of the objectives.', class(obj), obj.defaultQuantity);
        end
    end

    methods (Static)
        % Deprecated, use the setBiologicalDosePrescriptions of the quantity
        function optiFunc = setBiologicalDosePrescriptions(optiFunc, ~, ~)
        end
    end
end
