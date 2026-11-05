function E = cns_evaluateReferenceDose(ct, cst, refDose, caseCfg, study)
% cns_evaluateReferenceDose  Same metrics as the proton plans, for the clinical (photon) dose.
%
%   E = cns_evaluateReferenceDose(ct, cst, refDose, caseCfg, study)
%
% refDose comes from cns_importDicom (RTDOSE on the CT grid). The total dose
% is refDose.cube * refDose.toTotal, or * the case file's numOfFractions
% when toTotal is NaN (per-fraction RTDOSE without an imported RTPLAN).
% Photons: RBE = 1, so the RBE11 and McN rows equal the physical dose and
% the LET rows are zero. That lets the cohort table put the photon plan in
% the same columns as the proton plans.

nFx = caseCfg.prescription.numOfFractions;
f = refDose.toTotal;
if isnan(f), f = nFx; end
total = double(refDose.cube) * f;

if isfield(refDose, 'numOfFractions') && ~isnan(refDose.numOfFractions) && refDose.numOfFractions ~= nFx
    warning('cns:reference', 'RTPLAN has %d fractions, the case file %d. Check the prescription.', ...
        refDose.numOfFractions, nFx);
end

% pseudo-dij: one "spot" whose dose per fraction is the whole cube
dij = struct();
dij.physicalDose = {total(:) / nFx};
dij.doseGrid = struct('dimensions', ct.cubeDim, 'resolution', ct.resolution);
dij.ctGrid = dij.doseGrid;
E = cns_evaluatePlan(dij, cst, 1, caseCfg, study, struct('photon', true));

% sanity check against the prescription
ix = find(strcmp(cst(:,2), 'PTV'), 1);
if isempty(ix), ix = find(strcmp(cst(:,2), 'CTV'), 1); end
if ~isempty(ix)
    med = median(total(cst{ix,4}{1}));
    rx = caseCfg.prescription.totalDoseGyRBE;
    fprintf('[cns] reference dose: median %s dose %.1f Gy (prescription %.1f)\n', cst{ix,2}, med, rx);
    if med < 0.8 * rx || med > 1.2 * rx
        warning('cns:reference', ['Reference %s median dose is %.0f %% of the prescription: check ' ...
            'the fraction handling (refDose.toTotal) or the prescription.'], cst{ix,2}, 100 * med / rx);
    end
end
end
