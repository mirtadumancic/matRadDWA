clear all;
matRad_cfg = MatRad_Config.instance();

load('TG119.mat');
%%

pln.radiationMode = 'protons';

pln.propDoseCalc.calcLET = 0;

pln.numOfFractions       = 30;
pln.propStf.gantryAngles = [45, 315];
pln.propStf.couchAngles   = zeros(numel(pln.propStf.gantryAngles),1);
pln.propStf.bixelWidth    = 3;

pln.propStf.isoCenter     = matRad_getIsoCenter(cst,ct,0);
pln.propStf.longitudinalSpotSpacing = 2;
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

pln.multScen = matRad_multScen(ct, 'nomScen');

pln.propStf.bixelWidth = 5;

%% DOSE CALCULATION WITH DWA MACHINE
pln.machine       = 'DWA_proton_v2';
stf_DWA = matRad_generateStf(ct,cst,pln);

engine = DoseEngines.matRad_ParticleHongPencilBeamEngine(pln);
engine.calcLET = true;
dij_DWA = engine.calcDoseInfluence(ct,cst,stf_DWA);

% resultGUI_DWA = matRad_fluenceOptimization(dij_DWA,cst,pln);
% %% DOSE CALCULATION WITH GENRIC MACHINE
% pln.machine       = 'Generic';
% stf_generic = matRad_generateStf(ct,cst,pln);
% 
% engine = DoseEngines.matRad_ParticleHongPencilBeamEngine(pln);
% dij_generic = engine.calcDoseInfluence(ct,cst,stf_generic);
% % resultGUI_generic = matRad_fluenceOptimization(dij_generic,cst,pln);
%% PHYSICAL DOSE PLAN OPTIMIZATION WITH DWA MACHINE

for i=1:size(cst,1)
    cst{i,6} = {};
end

% TARGET
cst{2,6}{1} = DoseObjectives.matRad_SquaredDeviation(100,42);
cst{2,6}{1}.quantity = 'physicalDose';

% OAR
cst{1,6}{1} = DoseObjectives.matRad_SquaredOverdosing(100,42);
cst{1,6}{1}.quantity = 'physicalDose';

% OAR
cst{3,6}{1} = DoseObjectives.matRad_SquaredOverdosing(100,42);
cst{3,6}{1}.quantity = 'physicalDose';

resultGUI_physicalDose = matRad_fluenceOptimization(dij_DWA,cst,pln);

%% PHYSICAL DOSE + LETd PLAN OPTIMIZATION WITH DWA MACHINE

for i=1:size(cst,1)
    cst{i,6} = {};
end

% TARGET
cst{2,6}{1} = DoseObjectives.matRad_SquaredDeviation(100,42);
cst{2,6}{1}.quantity = 'physicalDose';

% OAR
cst{1,6}{1} = DoseObjectives.matRad_SquaredOverdosing(100,42);
cst{1,6}{1}.quantity = 'physicalDose';

cst{1,6}{2} = DoseObjectives.matRad_SquaredOverdosing(1,0);
cst{1,6}{2}.quantity = 'LETd';

% OAR
cst{3,6}{1} = DoseObjectives.matRad_SquaredOverdosing(100,42);
cst{3,6}{1}.quantity = 'physicalDose';

resultGUI_physicalDoseLETd = matRad_fluenceOptimization(dij_DWA,cst,pln);



%%
figure,
nexttile;
matRad_plotSlice(ct,'cst', cst, 'dose', resultGUI_physicalDose.physicalDose, 'axesHandle', gca(), 'plane', 3);
title('Physical Dose Optimization - Physical Dose');

nexttile;
matRad_plotSlice(ct,'cst', cst, 'dose', resultGUI_physicalDoseLETd.physicalDose, 'axesHandle', gca(), 'plane', 3);
title('LETd Optimization - Physical Dose');

nexttile;

matRad_plotSlice(ct,'cst', cst, 'dose', resultGUI_physicalDose.LET, 'axesHandle', gca(), 'plane', 3);
title('Physical Dose Optimization - LETd');

nexttile;
matRad_plotSlice(ct,'cst', cst, 'dose', resultGUI_physicalDoseLETd.LET, 'axesHandle', gca(), 'plane', 3);
title('LETd Optimization - LETd');
%% DVHS


pdOpt_pd   = matRad_calcDVH(cst,resultGUI_physicalDose.physicalDose);
pdOpt_let  = matRad_calcDVH(cst,resultGUI_physicalDose.LET);
letOpt_pd  = matRad_calcDVH(cst,resultGUI_physicalDoseLETd.physicalDose);
letOpt_let = matRad_calcDVH(cst,resultGUI_physicalDoseLETd.LET);


figure;

nexttile;
structIdx = 2;

plot(pdOpt_pd(structIdx).doseGrid, pdOpt_pd(structIdx).volumePoints, '-k', 'DisplayName','Physical Dose Optimization');
hold on;
plot(letOpt_pd(structIdx).doseGrid, letOpt_pd(structIdx).volumePoints, '-r', 'DisplayName','LETd Optimization');
xlabel('Dose [Gy]');
ylabel('Volume [%]');
grid on;
title('Physical Dose - PTV');

nexttile;
structIdx = 1;

plot(pdOpt_let(structIdx).doseGrid, pdOpt_let(structIdx).volumePoints, '-k', 'DisplayName','Physical Dose Optimization');
hold on;
plot(letOpt_let(structIdx).doseGrid, letOpt_let(structIdx).volumePoints, '-r', 'DisplayName','LETd Optimization');
xlabel('LETd [keV/\mu m]');
ylabel('Volume [%]');
grid on;
title('LETd - Core');
legend;
