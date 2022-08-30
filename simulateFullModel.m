function simulateFullModel(baseModel, outFilename, fracN, fracA, TNM, TNG, TAM, TAG, mobUtilN, mobUtilA)
% simulateFullModel
%
% Generates a combined model from the params supplied and runs a simulation.
% The fluxes of glycolysis, mitochondrial respiration, and both are then extracted
% for both neurons and astrocytes in all slices and written to file.
%
% Input:
%
%   baseModel       Should be an ec model generated with CreateCompartmentECGEM, 
%                   and addPenaltiesToModel should not be used on it.
%
%   outFilename     Filename to write the results to. A .mat file name is expected.
%
%   fracN           Fraction of neurons in the model
%
%   fracA           Fraction of astrocytes in the model
%
%   TNM             Transportation penalty neurons, mitochondria - typically a number between 0 and 2, such as 0.10
%
%   TNG             Transportation penalty neurons, glycolysis - typically a number between 0 and 2, such as 0.10
%
%   TAM             Transportation penalty astrocytes, mitochondria - typically a number between 0 and 2, such as 0.10
%
%   TAG             Transportation penalty astrocytes, glycolysis - typically a number between 0 and 2, such as 0.10
%
%   mobUtilN        Fraction of unused time in neurons that mitochondria can be used somewhere else
%
%   mobUtilA        Fraction of unused time in astrocytes that mitochondria can be used somewhere else
%
% Output:
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests with the new model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fullUtilModel = buildFullUtilBrainModel(baseModel, fracN, fracA, TNM, TNG, TAM, TAG, mobUtilN, mobUtilA);

mu = fullUtilModel;
%glucose + unl oxygen
mu.ub(strcmp(mu.rxns, 'MAR09034_REV')) = 1; %Supply glucose - oxygen should already be open
mu.ub(strcmp(mu.rxns, 'MAR09135_REV')) = 0; %Do not supply lactate from the outside
mu.c = double(strcmp(mu.rxns, 'tot_atp_hydr'));

%mu.ub(strcmp(mu.rxns, 'MAR09048_REV')) = Inf;
%run simulation
res = solveLP(mu,1);
res %-31.2660 - looks good

%look at lactate output and export in each slice
%

%table(baseModel.rxns,constructEquations(baseModel))

% Lactate export: MAR06048, import: MAR08515

netImpNeur = nan(100,1);
netImpAstr = nan(100,1);
ATPConsNeur = nan(100,1);
ATPConsAstr = nan(100,1);
glycNeur = nan(100,1);
glycAstr = nan(100,1);
mitoNeur = nan(100,1);
mitoAstr = nan(100,1);


ATPConsTot = -res.f;
fracA = 0.03;
fracN = 0.17;

for i = 1:100
    neurExpSel = strcmp(mu.rxns, ['N_' num2str(i) '_MAR06048']);
    neurImpSel = strcmp(mu.rxns, ['N_' num2str(i) '_MAR08515']);
    astrExpSel = strcmp(mu.rxns, ['A_' num2str(i) '_MAR06048']);
    astrImpSel = strcmp(mu.rxns, ['A_' num2str(i) '_MAR08515']);

    neurGlycSel = strcmp(mu.rxns, ['N_' num2str(i) '_MAR04394']);
    astrGlycSel = strcmp(mu.rxns, ['A_' num2str(i) '_MAR04394']);
    neurMitSel = strcmp(mu.rxns, ['N_' num2str(i) '_MAR04152']);
    astrMitSel = strcmp(mu.rxns, ['A_' num2str(i) '_MAR04152']);

    netImpNeur(i) = res.x(neurImpSel) - res.x(neurExpSel);
    netImpAstr(i) = res.x(astrImpSel) - res.x(astrExpSel);
    
    scale = i/(100*50.5);%The average slice takes up 50.5 times more than the 0.01 util slice (mean of 1:100), and we have 100 slices
    ATPConsNeur(i) = ATPConsTot * scale*fracN;
    ATPConsAstr(i) = ATPConsTot * scale*fracA;

    glycNeur(i) = res.x(neurGlycSel);
    glycAstr(i) = res.x(astrGlycSel);
    mitoNeur(i) = res.x(neurMitSel);
    mitoAstr(i) = res.x(astrMitSel);

end

%ATPProdNeur = glycNeur.*2 + mitoNeur.*14.75;
%ATPProdAstr = glycAstr.*2 + mitoAstr.*14.75;


d = struct();
d.netImpNeur = netImpNeur;
d.netImpAstr = netImpAstr;
d.ATPConsNeur = ATPConsNeur;
d.ATPConsAstr = ATPConsAstr;
d.ATPProdGlycNeur = glycNeur.*2;
d.ATPProdGlycAstr = glycAstr.*2;
d.ATPProdMitoNeur = mitoNeur.*14.75;
d.ATPProdMitoAstr = mitoAstr.*14.75;

save(outFilename, 'd')

end

