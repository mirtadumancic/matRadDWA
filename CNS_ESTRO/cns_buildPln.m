function pln = cns_buildPln(ct, cst, caseCfg, study, machine)
% cns_buildPln  matRad pln for one case and machine (same settings as the
% LETd water-phantom scripts: HongPB, double Gaussian, LET on, constant RBE 1.1).
nFx = caseCfg.prescription.numOfFractions;

pln = struct();
pln.radiationMode = 'protons';
pln.machine = machine;
pln.bioModel = 'constRBE';
pln.numOfFractions = nFx;

g = cns_getfield(caseCfg, 'beams.gantryAngles', []);
if isempty(g), error('cns:beams', 'Set beams.gantryAngles in the case file.'); end
pln.propStf.gantryAngles = g(:)';
c = cns_getfield(caseCfg, 'beams.couchAngles', zeros(size(g)));
pln.propStf.couchAngles = c(:)';
% safety: if the case names a beam set of this study, its angles must match (setId is only a label)
sid = cns_getfield(caseCfg, 'beams.setId', '');
sets = cns_toCell(cns_getfield(study, 'beamCandidates.sets', {}));
for i = 1:numel(sets)
    if strcmp(sets{i}.id, sid) && (~isequal(sets{i}.gantryAngles(:)', g(:)') || ~isequal(sets{i}.couchAngles(:)', c(:)'))
        error('cns:beams', ['Case file says beams.setId = %s but gantry [%s] / couch [%s]; set %s is gantry [%s] / couch [%s]. ' ...
            'Copy the angles of the set into the case file.'], sid, num2str(g(:)'), num2str(c(:)'), sid, ...
            num2str(sets{i}.gantryAngles(:)'), num2str(sets{i}.couchAngles(:)'));
    end
end

bw = cns_getfield(study, ['stf.bixelWidthMm.' machine], cns_getfield(study, 'stf.bixelWidthMm.default', 5));
pln.propStf.bixelWidth = bw;
pln.propStf.longitudinalSpotSpacing = cns_getfield(study, 'stf.longitudinalSpotSpacingMm', 3);
iso = cns_getfield(caseCfg, 'beams.isoCenter', []);
if isempty(iso)
    iso = matRad_getIsoCenter(cst, ct, 0);
end
pln.propStf.isoCenter = repmat(iso(:)', numel(g), 1);
rs = cns_getfield(study, ['stf.rangeShifterEqDmm.' machine], []);
if ~isempty(rs) && rs > 0
    pln.propStf.useRangeShifter = true;
    pln.propStf.rangeShifterEqD = rs;
end

pln.propDoseCalc.engine = cns_getfield(study, 'doseCalc.engine', 'HongPB');
pln.propDoseCalc.lateralModel = cns_getfield(study, 'doseCalc.lateralModel', 'double');
cut = cns_getfield(study, 'doseCalc.geometricLateralCutOffMm', []);
if ~isempty(cut), pln.propDoseCalc.geometricLateralCutOff = cut; end
pln.propDoseCalc.calcLET = 1;
r = cns_getfield(study, 'doseCalc.doseGridMm', 3);
pln.propDoseCalc.doseGrid.resolution = struct('x', r, 'y', r, 'z', r);

if exist('matRad_OptimizationQuantity', 'class') == 8      % let_opt branch
    pln.propOpt.quantityOpt = 'RBExDose';
else
    pln.propOpt.quantityOpt = 'RBExD';
end
pln.propOpt.runDAO = 0;
pln.propSeq.runSequencing = 0;
pln.propOpt.optimizer = cns_getfield(study, 'optimisation.optimizer', 'IPOPT');

rob = cns_getfield(study, 'robustness', struct('scenarios', 'nomScen'));
pln.multScen = matRad_multScen(ct, rob.scenarios);
if ~strcmp(rob.scenarios, 'nomScen')
    if isfield(rob, 'shiftSDmm'),     pln.multScen.shiftSD = rob.shiftSDmm(:)'; end
    if isfield(rob, 'rangeRelSDpct'), pln.multScen.rangeRelSD = rob.rangeRelSDpct; end
    if isfield(rob, 'rangeAbsSDmm'),  pln.multScen.rangeAbsSD = rob.rangeAbsSDmm; end
end
end
