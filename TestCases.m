%%%%%%%%%%%%%%%%%%%%
%TC0001 - addPenaltiesToModel
%%%%%%%%%%%%%%%%%%%%
cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ecModel = load('data/MinModel.mat').minModel;

%Check that there is no penalty in the no penalty case
m = addPenaltiesToModel(ecModel, 1, 1, 1, 1, 1, true);
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_other')) == [0;0])
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_cyto')) == [0;0])
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_mito')) == [0;0])
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_mt')) == [0;0])

%Now with penalty
%model, B, TM, TG, UCyto, UMito, dontPenalize
m = addPenaltiesToModel(ecModel, 1, 0.1, 0.2, 0.3, 0.4, false);
%check that we consume ATP
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_other')) == [-1;0])
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_other')) == [-1;0])
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_other')) == [-1;0])
all(m.S(strcmp(m.metNames, 'ATP'), strcmp(m.rxns, 'prot_maint_other')) == [-1;0])
%check that the protein penalty is right (according to the formulas in the paper)
%there is a roundoff problem here since the function converts the numbers to text in addRxns, which is ok.
abs(m.S(strcmp(m.metNames, 'other_prot_pool'), strcmp(m.rxns, 'prot_maint_other')) - 1/((1+0.2)/0.3)) < 10^-5 %same as cyto
abs(m.S(strcmp(m.metNames, 'cyto_prot_pool'), strcmp(m.rxns, 'prot_maint_cyto')) - 1/((1+0.2)/0.3)) < 10^-5
abs(m.S(strcmp(m.metNames, 'mito_prot_pool'), strcmp(m.rxns, 'prot_maint_mito')) - 1/((1+0.1)/0.4)) < 10^-5
abs(m.S(strcmp(m.metNames, 'mt_prot_pool'), strcmp(m.rxns, 'prot_maint_mt')) - 1/((1)/0.4)) < 10^-5

%All ok!

%%%%%%%%%%%%%%%%%%%%
%TC0002 - buildFullUtilBrainModel
%%%%%%%%%%%%%%%%%%%%

cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ecModel = load('data/MinModel.mat').minModel;
%buildFullUtilBrainModel(baseModel, fracN, fracA, TNM, TNG, TAM, TAG, mobUtilN, mobUtilA)
fracN = 0.17;
fracA = 0.03;
TNM = 0.1;
TNG = 0.2;
TAM = 0.3;
TAG = 0.4;
mobUtilN = 0.5;
mobUtilA = 0.6;
fm = buildFullUtilBrainModel(ecModel, fracN, fracA, TNM, TNG, TAM, TAG, mobUtilN, mobUtilA);

%Step 1 - Create the combined model. Check that the maintenance costs are right.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%look at slice 37
%neurons
us = 0.37;
uMito = us + (1-us)*mobUtilN;
uCyto = us;
nComp = find(strcmp(fm.comps, 'N_37_c'));
abs(fm.S(strcmp(fm.metNames, 'other_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'N_37_prot_maint_other')) - 1/((1+TNG)/uCyto)) < 10^-5 %same as cyto
abs(fm.S(strcmp(fm.metNames, 'cyto_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'N_37_prot_maint_cyto')) - 1/((1+TNG)/uCyto)) < 10^-5 
abs(fm.S(strcmp(fm.metNames, 'mito_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'N_37_prot_maint_mito')) - 1/((1+TNM)/uMito)) < 10^-5 
abs(fm.S(strcmp(fm.metNames, 'mt_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'N_37_prot_maint_mt')) - 1/(1/uMito)) < 10^-5 

%astrocytes
us = 0.37;
uMito = us + (1-us)*mobUtilA;
uCyto = us;
nComp = find(strcmp(fm.comps, 'A_37_c'));
abs(fm.S(strcmp(fm.metNames, 'other_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'A_37_prot_maint_other')) - 1/((1+TAG)/uCyto)) < 10^-5 %same as cyto
abs(fm.S(strcmp(fm.metNames, 'cyto_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'A_37_prot_maint_cyto')) - 1/((1+TAG)/uCyto)) < 10^-5 
abs(fm.S(strcmp(fm.metNames, 'mito_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'A_37_prot_maint_mito')) - 1/((1+TAM)/uMito)) < 10^-5 
abs(fm.S(strcmp(fm.metNames, 'mt_prot_pool') & (fm.metComps == nComp), strcmp(fm.rxns, 'A_37_prot_maint_mt')) - 1/(1/uMito)) < 10^-5 

%all ok

%Step 2 - Add glucose uptake to slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%It is obvious from the modeling results that this works reasonably well, no need to test it further

%Step 3 - Add an uptake of lactate from rob in all slices - but only say 1/5 of the glucose uptake. Also transport in the opposite direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%It is obvious from the modeling results lactate import works reasonably well (see for example region I in Fig. 2F, neurons take up lactate)
%It is also obvious from the modeling results that lactate export works reasonably well (see for example region II in Fig. 2F, ANLS is present)

%Step 4 - Add communication of lactate from slices with higher utilization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check that communication of lactate is available from higher to lower utilization, but not the other way around.
constructEquations(fm, 'lac_slice_transp_R_79_80_') %{'L-lactate[N_80_s] => L-lactate[N_79_s]'}, ok
sum(strcmp(fm.rxns, 'lac_slice_transp_R_80_79_')) == 0  %doesn't exist, ok

%all ok

%Step 5 - Create a common ATP hydrolysis function that consumes ATP from all slices in the right amounts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check that the total ATP amount sums to 1 from all submodels
sel = strcmp(fm.rxns,'tot_atp_hydr');
subS = fm.S(:, sel);
abs(sum(subS(subS > 0)) - 3) < 10^-5  %3, as it should be (ADP + Pi + H+), ok
abs(sum(subS(subS < 0)) + 2) < 10^-5  %-2, as it should be (H2O + ATP), ok
%check that the slice with static utilization 0.8 has twice as high stoiciometric coeff as the one with 0.4
comp1 = find(strcmp(fm.comps, 'N_80_c'));
comp2 = find(strcmp(fm.comps, 'N_40_c'));
c80 = fm.S(strcmp(fm.metNames, 'ATP') & (fm.metComps == comp1), sel);
c40 = fm.S(strcmp(fm.metNames, 'ATP') & (fm.metComps == comp2), sel);
abs(c80 - 2*c40) < 10^-5 %ok (c80 = -0.0026931, so no numerical issues)

%all ok

%Step 6 - Constrain the lactate uptake in the Rob to ~25% of the expected glucose uptake (then assuming no lactate uptake)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs(fm.ub(strcmp(fm.rxns, 'MAR08515'))  - 0.8 * 2 * 0.25) < 10^-5 %ok

%all ok

%all steps ok


%%%%%%%%%%%%%%%%%%%%
%TC0003 - buildSimpleAstroNeuronModel
%%%%%%%%%%%%%%%%%%%%

cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ecModel = load('data/MinModel.mat').minModel;
%buildSimpleAstroNeuronModel(baseModel, fracA, TNM, TNG, TAM, TAG, statUtil, mobUtilN, mobUtilA)
fracA = 0.15;
TNM = 0.1;
TNG = 0.2;
TAM = 0.3;
TAG = 0.4;
statUtil = 0.4;
mobUtilN = 0.5;
mobUtilA = 0.6;
sm = buildSimpleAstroNeuronModel(ecModel, fracA, TNM, TNG, TAM, TAG, statUtil, mobUtilN, mobUtilA);

%Step 1 - Create the combined model. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check that the maintenance costs are right

uMito = statUtil + (1-statUtil)*mobUtilN;
uCyto = statUtil;
nComp = find(strcmp(sm.comps, 'c'));
abs(sm.S(strcmp(sm.metNames, 'other_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'prot_maint_other')) - 1/((1+TNG)/uCyto)) < 10^-5 %same as cyto
abs(sm.S(strcmp(sm.metNames, 'cyto_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'prot_maint_cyto')) - 1/((1+TNG)/uCyto)) < 10^-5 
abs(sm.S(strcmp(sm.metNames, 'mito_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'prot_maint_mito')) - 1/((1+TNM)/uMito)) < 10^-5 
abs(sm.S(strcmp(sm.metNames, 'mt_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'prot_maint_mt')) - 1/(1/uMito)) < 10^-5 

%astrocytes
uMito = statUtil + (1-statUtil)*mobUtilA;
uCyto = statUtil;
nComp = find(strcmp(sm.comps, 'A_c'));
abs(sm.S(strcmp(sm.metNames, 'other_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'A_prot_maint_other')) - 1/((1+TAG)/uCyto)) < 10^-5 %same as cyto
abs(sm.S(strcmp(sm.metNames, 'cyto_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'A_prot_maint_cyto')) - 1/((1+TAG)/uCyto)) < 10^-5 
abs(sm.S(strcmp(sm.metNames, 'mito_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'A_prot_maint_mito')) - 1/((1+TAM)/uMito)) < 10^-5 
abs(sm.S(strcmp(sm.metNames, 'mt_prot_pool') & (sm.metComps == nComp), strcmp(sm.rxns, 'A_prot_maint_mt')) - 1/(1/uMito)) < 10^-5 

%check that the astrocytes cannot import lactate/glucose directly, but need to get them from the neurons
constructEquations(sm, {'A_MAR09034_REV';'A_MAR09135_REV'}) %these are the correct reactions
sm.ub(strcmp(sm.rxns, 'A_MAR09034_REV')) == 0 %%ok
sm.ub(strcmp(sm.rxns, 'A_MAR09135_REV')) == 0 %%ok

%all ok

%Step 2 - Add communication of lactate and glucose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constructEquations(sm, {'lactate_transport'; 'glucose_transport'}) %looks ok, they exist

%Step 3 - Create a common ATP hydrolysis function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check that they sum to 1
constructEquations(sm, 'tot_atp_hydr')
sel = strcmp(sm.rxns,'tot_atp_hydr');
subS = sm.S(:, sel);
abs(sum(subS(subS > 0)) - 3) < 10^-5  %3, as it should be (ADP + Pi + H+), ok
abs(sum(subS(subS < 0)) + 2) < 10^-5  %-2, as it should be (H2O + ATP), ok
%check that the slice with static utilization 0.8 has twice as high stoiciometric coeff as the one with 0.4
compN = find(strcmp(sm.comps, 'c'));
compA = find(strcmp(sm.comps, 'A_c'));
cN = sm.S(strcmp(sm.metNames, 'ATP') & (sm.metComps == compN), sel);
cA = sm.S(strcmp(sm.metNames, 'ATP') & (sm.metComps == compA), sel);
cN == -0.85
cA == -0.15

%all ok


%all steps ok

%%%%%%%%%%%%%%%%%%%%
%TC0004 - CreateCompartmentECGEM and GetMWAndKcatsNeurons
%%%%%%%%%%%%%%%%%%%%

%The code is mainly copied from GECKO Light and tested there. However,
%we should test the split up of protein usage

cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ecModel = load('data/MinModel.mat').minModel;

%we look at the min model, and see if this behaves as expected
%We know that complex I should have both MT and mito costs, but no other:
constructEquations(ecModel, 'MAR06921')
%{'5 H+[m] + NADH[m] + ubiquinone[m] + 0.081101 mito_prot_pool[c] + 0.024257 mt_prot_pool[c] => NAD+[m] + ubiquinol[m] + 4 H+[i]'}
%looks reasonable

%We know that the first step in glycolysis should be cytosolic proteins only:
constructEquations(ecModel, 'MAR04394')
% {'ATP[c] + glucose[c] + 0.00021195 cyto_prot_pool[c] => ADP[c] + glucose-6-phosphate[c] + H+[c]'}
%ok

%table(ecModel.rxns,constructEquations(ecModel))

%Also test that the new function does the same as GECKO Light - this requires GECKO Light to run
%We test that the sum of the 4 enzyme category costs are the same as the single one in GECKO Light.
ihuman = load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat').ihuman;

ecLightModel = CreateECLightModel(ihuman, true, 1);

minS = ecLightModel.S(strcmp(ecLightModel.mets, 'prot_pool'), ismember(ecLightModel.rxns, ecModel.rxns));

minS2 = sum(ecModel.S((length(ecModel.mets)-3):length(ecModel.mets), :),1)

diff = abs(minS(:) - minS2(:));
all (diff < 10^-14) %ok

%All ok


