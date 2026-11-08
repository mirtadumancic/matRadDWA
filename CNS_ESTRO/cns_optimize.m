function [resultGUI, pln, usedOptimizer] = cns_optimize(dij, cst, pln, study, w0)
% cns_optimize  matRad_fluenceOptimization with a robust check of the IPOPT mex.
%
%   [resultGUI, pln, used] = cns_optimize(dij, cst, pln, study)
%   [resultGUI, pln, used] = cns_optimize(dij, cst, pln, study, w0)
%
% matRad builds a new IPOPT optimiser for every plan and refuses to run when its
% tiny test problem fails ("IPOPT mex interface not available"). On macOS this can
% happen intermittently (the mex is still loaded in a bad state after an earlier run
% or an interrupted one). This wrapper:
%   1. tests IPOPT itself and prints the real reason when it fails;
%   2. unloads the mex (clear ipopt), rehashes the path and tests again;
%   3. if IPOPT still fails, uses study.optimisation.fallbackOptimizer when it is set
%      (e.g. "fmincon"), otherwise stops with a clear message.
% No silent switch: the optimiser used is returned and printed, so mixed-optimiser
% comparisons can be spotted.
if nargin < 5, w0 = []; end
used = cns_getfield(pln, 'propOpt.optimizer', 'IPOPT');

if strcmpi(used, 'IPOPT')
    [ok, why] = ipoptWorks();
    if ~ok
        fprintf(2, '[cns] IPOPT test failed (%s). Unloading the mex and retrying...\n', why);
        clear ipopt
        rehash path
        [ok, why] = ipoptWorks();
    end
    if ~ok
        fb = cns_getfield(study, 'optimisation.fallbackOptimizer', '');
        diag = sprintf('arch %s, mexext %s, ipopt found at: %s, reason: %s', ...
            computer('arch'), mexext, strtrim(evalc('which ipopt')), why);
        if isempty(fb)
            error(['IPOPT is not usable in this MATLAB (%s).\n' ...
                'Fix: restart MATLAB and run again; if it persists, see the README (IPOPT on Mac), ' ...
                'or set "fallbackOptimizer": "fmincon" under "optimisation" in the study JSON.'], diag);
        end
        fprintf(2, '[cns] IPOPT still not usable (%s): using %s for this plan.\n', diag, fb);
        used = fb;
        pln.propOpt.optimizer = fb;
    end
end

% (Nov 8) keep the result cubes on the dose grid: matRad_fluenceOptimization computes
% cubes for every robust scenario and would interpolate all of them to the sub-mm CT
% (out of memory in matRad_calcCubes). Only the weights are used downstream.
matRad_cfg = MatRad_Config.instance();
hadFlag = isstruct(matRad_cfg.defaults) && isfield(matRad_cfg.defaults, 'calcCubesInterpToCtGrid');
oldFlag = true; if hadFlag, oldFlag = matRad_cfg.defaults.calcCubesInterpToCtGrid; end
matRad_cfg.defaults.calcCubesInterpToCtGrid = false;
try
    if isempty(w0)
        resultGUI = matRad_fluenceOptimization(dij, cst, pln);
    else
        resultGUI = matRad_fluenceOptimization(dij, cst, pln, w0);
    end
catch err
    restoreFlag(matRad_cfg, hadFlag, oldFlag);
    rethrow(err);
end
restoreFlag(matRad_cfg, hadFlag, oldFlag);
usedOptimizer = used;
end

function restoreFlag(matRad_cfg, had, oldFlag)
if had
    matRad_cfg.defaults.calcCubesInterpToCtGrid = oldFlag;
elseif isfield(matRad_cfg.defaults, 'calcCubesInterpToCtGrid')
    matRad_cfg.defaults = rmfield(matRad_cfg.defaults, 'calcCubesInterpToCtGrid');
end
end

function [ok, why] = ipoptWorks()
% Same test problem as matRad_OptimizerIPOPT.isAvailable, but keeps the error message.
ok = false; why = '';
if exist(['ipopt.' mexext], 'file') ~= 3
    why = sprintf('no ipopt.%s on the path', mexext);
    return;
end
funcs.objective        = @(x) x.^2;
funcs.gradient         = @(x) 2 * x;
funcs.hessian          = @(x, sigma, lambda) sigma * sparse(2);
funcs.hessianstructure = @() sparse(true);
s.ipopt.tol = 1e-5;
s.ipopt.print_level = 0;
s.ipopt.print_user_options = 'no';
s.ipopt.print_options_documentation = 'no';
try
    [x, ~] = ipopt(1, funcs, s);   % the mex requires two outputs
    ok = abs(x) < 1e-5;
    if ~ok, why = sprintf('test problem returned x = %g', x); end
catch ME
    why = ME.message;
end
end
