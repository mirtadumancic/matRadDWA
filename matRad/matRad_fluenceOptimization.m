function [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
%
% call:
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
%
% input:
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   wInit:      (optional) custom weights to initialize problems
%
% output:
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
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

if exist('wInit','var')
    [dij,cst,pln,wInit,optiProb] = matRad_initOptimization(dij,cst,pln,wInit);
else
    [dij,cst,pln,wInit,optiProb] = matRad_initOptimization(dij,cst,pln);
end

if ~isfield(pln.propOpt,'optimizer')
    %While the default optimizer is IPOPT, we can try to fallback to
    %fmincon in case it does not work for some reason
    if ~matRad_OptimizerIPOPT.isAvailable()
        pln.propOpt.optimizer = 'fmincon';
    else
        pln.propOpt.optimizer = 'IPOPT';
    end    
end

if isstring(pln.propOpt.optimizer) || ischar(pln.propOpt.optimizer)
    pln.propOpt.optimizer = struct('name',pln.propOpt.optimizer);
end    

if isstruct(pln.propOpt.optimizer) && isfield(pln.propOpt.optimizer,'name')
    optimizerName = char(pln.propOpt.optimizer.name);
    optimizerOptions = rmfield(pln.propOpt.optimizer,'name');
else
    matRad_cfg.dispError('Could not identify optimizer! Please provide a valid optimizer name or optimizer struct with field ''name''!');
end

switch optimizerName
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    case 'simulannealbnd'
        optimizer = matRad_OptimizerSimulannealbnd;
    otherwise
        warning(['Optimizer ''' optimizerName ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end

matRad_assignPropertiesFromStruct(optimizer,optimizerOptions);

if ~optimizer.isAvailable()
    matRad_cfg.dispError(['Optimizer ''' optimizerName ''' not available!']);
end

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;

resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;

%Robust quantities
if pln.multScen.totNumScen > 1
    for i = 1:pln.multScen.totNumScen
        scenSubIx = pln.multScen.linearMask(i,:);
        resultGUItmp = matRad_calcCubes(wOpt,dij,pln.multScen.sub2scenIx(scenSubIx(1),scenSubIx(2),scenSubIx(3)));
        resultGUI = matRad_appendResultGUI(resultGUI,resultGUItmp,false,sprintf('scen%d',i));
    end
end

% unblock mex files
clear mex
