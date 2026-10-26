function c = matRad_constraintFunctions(optiProb,w,dij,cst)
% matRad IPOPT callback: constraint function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, 
% max DVH constraint, min DVH constraint 
% 
% call:
%   c = matRad_constraintFunctions(optiProb,w,dij,cst)
%
% input:
%   optiProb:   option struct defining the type of optimization
%   w:          bixel weight vector
%   dij:        dose influence matrix
%   cst:        matRad cst struct
%
% output:
%   c:          value of constraints
%
% References
%   -
%
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

% retrieve matching 4D scenarios
fullScen = cell(ndims(dij.physicalDose),1);
[fullScen{:}] = ind2sub(size(dij.physicalDose),useScen);
contourScen = fullScen{1};

% Initializes constraints
c = [];

% min/max cache for voxel-wise worst case, per quantity
dTmp = struct();

% compute objective function for every VOI.
for  i = 1:size(cst,1)
   
   % Only take OAR or target VOI.
   if ~isempty(cst{i,4}{1}) && any(strcmp(cst{i,3},{'OAR','TARGET','EXTERNAL'}))
      
      % loop over the number of constraints for the current VOI
      for j = 1:numel(cst{i,6})
         
         constraint = cst{i,6}{j};
         
         % only perform computations for constraints
         if isa(constraint,'DoseConstraints.matRad_DoseConstraint')
            
            quantityName = optiProb.BP.getQuantityName(constraint);
            dQ = d.(quantityName);

            % rescale dose parameters to the optimization quantity if required
            constraint = optiProb.BP.getQuantity(quantityName).setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            % retrieve the robustness type
            robustness = constraint.robustness;
            
            switch robustness
               case 'none' % if conventional opt: just sum objectives of nominal dose
                   d_i = dQ{1}(cst{i,4}{1});
                   c = [c; constraint.computeDoseConstraintFunction(d_i)];
                  
               case 'PROB'
                  matRad_cfg.dispError('PROB robustness is currently not supported!');
                  
               case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~isfield(dTmp,quantityName)
                     dTmp.(quantityName) = [dQ{useScen}];
                  end
                  
                  d_Scen = dTmp.(quantityName)(cst{i,4}{contourIx},:);
                  
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_max;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_min;
                  end
                  
                  c = [c; constraint.computeDoseConstraintFunction(d_i)];
                  
               case 'VWWC_INV'  %inverse voxel-wise conformitiy - takes maximum dose in TARGET and minimum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D inverted VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~isfield(dTmp,quantityName)
                     dTmp.(quantityName) = [dQ{useScen}];
                  end
                  
                  d_Scen = dTmp.(quantityName)(cst{i,4}{contourIx},:);
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_min;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_max;
                  end
                  
                  c = [c; constraint.computeDoseConstraintFunction(d_i)];   
               otherwise
                  matRad_cfg.dispError('Robustness setting %s not yet supported!',constraint.robustness);
            end
            
            
         end
         
      end % if we are a constraint
      
   end % over all defined constraints & objectives
   
end % if structure not empty and oar or target

end % over all structures
