function jacob = matRad_constraintJacobian(optiProb,w,dij,cst)
% matRad IPOPT callback: jacobian function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 
% 
% call:
%   jacob = matRad_jacobFunc(optiProb,w,dij,cst)
%
% input:
%   optiProb: option struct defining the type of optimization
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%
% output:
%   jacob: jacobian of constraint function
%
% References
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

% voxel-wise constraint jacobians and their rows, per quantity (only first scenario supported)
fJacob = struct();
rowIx  = struct();
numOfConstraints = 0;

% get the used scenarios
useScen  = optiProb.BP.scenarios;

% retrieve matching 4D scenarios
fullScen      = cell(ndims(dij.physicalDose),1);
[fullScen{:}] = ind2sub(size(dij.physicalDose),useScen);
contourScen   = fullScen{1};

% min/max cache for voxel-wise worst case, per quantity
dTmp = struct();

% compute objective function for every VOI.
for i = 1:size(cst,1)
   
   % Only take OAR or target VOI.
   if ~isempty(cst{i,4}{1}) && any(strcmp(cst{i,3},{'OAR','TARGET','EXTERNAL'}))
      
      % loop over the number of constraints for the current VOI
      for j = 1:numel(cst{i,6})
         
         constraint = cst{i,6}{j}; %Get the Optimization Object
         
         % only perform computations for constraints
         if isa(constraint,'DoseConstraints.matRad_DoseConstraint')
            
            % retrieve the robustness type
            robustness = constraint.robustness;

            quantityName = optiProb.BP.getQuantityName(constraint);
            dQ = d.(quantityName);
            
            % rescale dose parameters to the optimization quantity if required
            constraint = optiProb.BP.getQuantity(quantityName).setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            switch robustness
               
               case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                     d_i = dQ{1}(cst{i,4}{1});
                     jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               case 'PROB'
                  matRad_cfg.dispError('PROB robustness is currently not supported!');
                  
               case {'VWWC','VWWC_INV'}  % voxel-wise worst case / inverted voxel-wise worst case
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
                  
                  if xor(isequal(cst{i,3},'OAR'),strcmp(robustness,'VWWC_INV'))
                     d_i = d_max;
                  else
                     d_i = d_min;
                  end
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               otherwise
                  matRad_cfg.dispError('Robustness setting %s not yet supported!',constraint.robustness);
            end
            
            nConst = size(jacobSub,2);
            
            if nConst > 0
               if ~isfield(fJacob,quantityName)
                  fJacob.(quantityName) = sparse(dij.doseGrid.numOfVoxels,0);
                  rowIx.(quantityName)  = [];
               end

               fJacobSub = sparse(dij.doseGrid.numOfVoxels,nConst);
               fJacobSub(cst{i,4}{1},:) = jacobSub;

               fJacob.(quantityName) = [fJacob.(quantityName), fJacobSub];
               rowIx.(quantityName)  = [rowIx.(quantityName), numOfConstraints + (1:nConst)];
               numOfConstraints = numOfConstraints + nConst;
            end
         end
      end
   end
end

if numOfConstraints == 0
   jacob = sparse([]);
   return;
end

% project the voxel-wise jacobians onto the weights
jacob = sparse(numOfConstraints,dij.totalNumOfBixels);
for name = fieldnames(fJacob)'
   jacob(rowIx.(name{1}),:) = optiProb.BP.getQuantity(name{1}).projectGradient(dij,fJacob.(name{1}),w,1)';
end
