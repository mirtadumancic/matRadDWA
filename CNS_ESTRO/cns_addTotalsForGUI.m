function R = cns_addTotalsForGUI(R, ct, cst, nFx, study)
% cns_addTotalsForGUI  Add total-dose cubes to a matRad resultGUI for display.
%
%   resultGUI = cns_addTotalsForGUI(resultGUI, ct, cst, nFx, study)
%
% Adds RBExDose_total (RBE 1.1 x fractions) and McN_total (McNamara,
% alpha/beta 10 Gy in CTV/PTV, 2 Gy elsewhere), both on the CT grid, so the
% matRad GUI can show them next to LET. matRad itself shows per-fraction dose.
if isfield(R, 'RBExDose'), D = R.RBExDose; else, D = 1.1 * R.physicalDose; end
R.RBExDose_total = D * nFx;
L = zeros(ct.cubeDim);
if isfield(R, 'LET'), L = R.LET; L(~isfinite(L)) = 0; end
ab = cns_getfield(study, 'evaluation.alphaBetaNormal', 2) * ones(ct.cubeDim);
for t = {'CTV', 'PTV'}
    ix = find(strcmp(cst(:,2), t{1}), 1);
    if ~isempty(ix), ab(cst{ix,4}{1}) = cns_getfield(study, 'evaluation.alphaBetaTarget', 10); end
end
R.McN_total = nFx * cns_mcNamara(R.physicalDose, L, ab);
end
