function weightGradient = matRad_objectiveGradient(optiProb,w,dij,cst)
% matRad IPOPT callback: gradient function for inverse planning
% supporting mean dose objectives, EUD objectives, squared overdosage,
% squared underdosage, squared deviation and DVH objectives
%
% call:
%   g = matRad_gradFuncWrapper(optiProb,w,dij,cst)
%
% input:
%   optiProb: option struct defining the type of optimization
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%
% output:
%   g: gradient of objective function
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

% get current quantities, d.(quantityName){scen}
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;
useNominalCtScen = optiProb.BP.nominalCtScenarios;

% retrieve matching 4D scenarios
fullScen      = cell(ndims(dij.physicalDose),1);
[fullScen{:}] = ind2sub(size(dij.physicalDose),useScen);
contourScen   = fullScen{1};

% objective gradient wrt each quantity, fGrad.(quantityName){scen}
fGrad = struct();
emptyGrad          = cell(size(dij.physicalDose));
emptyGrad(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};

% min/max cache for voxel-wise worst case, per quantity
dTmp = struct();

%For COWC
f_COWC = zeros(size(dij.physicalDose));

% compute objective function for every VOI.
for  i = 1:size(cst,1)
   
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && any(strcmp(cst{i,3},{'OAR','TARGET','EXTERNAL'}))
        
        % loop over the number of constraints and objectives for the current VOI
        for j = 1:numel(cst{i,6})
            
            %Get current optimization function
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective')
                
                % retrieve the robustness type
                robustness = objective.robustness;
                
                quantityName = optiProb.BP.getQuantityName(objective);
                dQ = d.(quantityName);
                if ~isfield(fGrad,quantityName)
                    fGrad.(quantityName) = emptyGrad;
                end

                % rescale dose parameters to the optimization quantity if required
                objective = optiProb.BP.getQuantity(quantityName).setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);

                switch robustness
                    case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                        for s = useNominalCtScen
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            d_i = dQ{ixScen}(cst{i,4}{ixContour});
                            %add to dose gradient
                            fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) = fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        end
                    case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = dQ{ixScen}(cst{i,4}{ixContour});
                            
                            fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) = fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) + ...
                                (objective.penalty*objective.computeDoseObjectiveGradient(d_i) * scenProb(s));
                            
                        end
                        
                    case 'PROB'
                        matRad_cfg.dispError('PROB robustness is currently not supported!');

                    case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                        end
                        
                        % prepare min/max dose vector for voxel-wise worst case
                        if ~isfield(dTmp,quantityName)
                            dTmp.(quantityName) = [dQ{useScen}];
                        end
                        
                        d_Scen = dTmp.(quantityName)(cst{i,4}{contourIx},:);
                        [d_max,max_ix] = max(d_Scen,[],2);
                        [d_min,min_ix] = min(d_Scen,[],2);
                        
                        if isequal(cst{i,3},'OAR')
                            d_i = d_max;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_min;
                        end
                        
                        if any(isnan(d_i))
                            matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                        end
                        
                        deltaTmp = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            if isequal(cst{i,3},'OAR')
                                currWcIx = double(max_ix == s);
                            elseif isequal(cst{i,3},'TARGET')
                                currWcIx = double(min_ix == s);
                            end
                            
                            fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) = fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                        end
                        
                    case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                        end
                        
                        % prepare min/max dose vector for voxel-wise worst case
                        if ~isfield(dTmp,quantityName)
                            dTmp.(quantityName) = [dQ{useScen}];
                        end
                        
                        d_Scen = dTmp.(quantityName)(cst{i,4}{1},:);
                        [d_max,max_ix] = max(d_Scen,[],2);
                        [d_min,min_ix] = min(d_Scen,[],2);
                        
                        if isequal(cst{i,3},'OAR')
                            d_i = d_min;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_max;
                        end
                        
                        if any(isnan(d_i))
                            matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                        end
                        
                        deltaTmp = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            if isequal(cst{i,3},'OAR')
                                currWcIx = double(min_ix == s);
                            elseif isequal(cst{i,3},'TARGET')
                                currWcIx = double(max_ix == s);
                            end
                            
                            fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) = fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                        end
                        
                    case 'COWC' % composite worst case consideres ovarall the worst objective function value
                        %First check the speficic cache for COWC
                        if ~exist('delta_COWC','var')
                            delta_COWC = struct();
                        end
                        if ~isfield(delta_COWC,quantityName)
                            delta_COWC.(quantityName) = emptyGrad;
                        end

                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);

                            d_i = dQ{ixScen}(cst{i,4}{ixContour});

                            f_COWC(ixScen) = f_COWC(ixScen) + objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                            delta_COWC.(quantityName){ixScen}(cst{i,4}{ixContour}) = delta_COWC.(quantityName){ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        end

                    case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                        f_OWC = zeros(size(dij.physicalDose));
                        delta_OWC = emptyGrad;

                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = dQ{ixScen}(cst{i,4}{ixContour});
                            
                            f_OWC(ixScen) = objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                            
                            delta_OWC{ixScen}(cst{i,4}{ixContour}) = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            
                        end
                          
                        switch optiProb.useMaxApprox
                            case 'logsumexp'
                                [~,fMaxGrad] = optiProb.logSumExp(f_OWC);
                            case 'pnorm'
                                [~,fMaxGrad] = optiProb.pNorm(f_OWC,numel(useScen));
                            case 'none'
                                [~,ix] = max(f_OWC(:));
                                fMaxGrad = zeros(size(f_OWC));
                                fMaxGrad(ix) = 1;
                            case 'otherwise'
                                matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                [~,ix] = max(f_OWC(:));
                                fMaxGrad = zeros(size(f_OWC));
                                fMaxGrad(ix) = 1;
                        end
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            if fMaxGrad(ixScen ) ~= 0
                                fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) = fGrad.(quantityName){ixScen}(cst{i,4}{ixContour}) + fMaxGrad(ixScen)*delta_OWC{ixScen}(cst{i,4}{ixContour});
                            end
                        end
                        
                    otherwise
                        matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);
                        
                end  %robustness type                              
            end  % objective check         
        end %objective loop       
    end %empty check    
end %cst structure loop

if exist('delta_COWC','var')   
    switch optiProb.useMaxApprox
        case 'logsumexp'
            [~,fMaxGrad] = optiProb.logSumExp(f_COWC);
        case 'pnorm'
            [~,fMaxGrad] = optiProb.pNorm(f_COWC,numel(useScen));
        case 'none'
            [~,ixCurrWC] = max(f_COWC(:));
            fMaxGrad = zeros(size(f_COWC));
            fMaxGrad(ixCurrWC) = 1;
        case 'otherwise'
            matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
            [~,ixCurrWC] = max(f_COWC(:));
            fMaxGrad = zeros(size(f_COWC));
            fMaxGrad(ixCurrWC) = 1;
    end
    
    for name = fieldnames(delta_COWC)'
        for s = 1:numel(useScen)
            ixScen = useScen(s);
            if fMaxGrad(ixScen) ~= 0
                fGrad.(name{1}){ixScen} = fGrad.(name{1}){ixScen} + fMaxGrad(ixScen)*delta_COWC.(name{1}){ixScen};
            end
        end
    end
end

optiProb.BP.computeGradient(dij,fGrad,w);
weightGradient = optiProb.BP.GetGradient();
