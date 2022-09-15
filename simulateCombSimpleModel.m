function simulateCombSimpleModel(baseModel, outFilename, fracA, TNM, TNG, TAM, TAG, statUtil, mobUtilN, mobUtilA)
% simulateCombSimpleModel
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
%   statUtil        Static utilization
%
%   mobUtilN        Fraction of unused time in neurons that mitochondria can be used somewhere else
%
%   mobUtilA        Fraction of unused time in astrocytes that mitochondria can be used somewhere else
%
% Output:
%

mu = buildSimpleAstroNeuronModel(baseModel, fracA, TNM, TNG, TAM, TAG, statUtil, mobUtilN, mobUtilA);

%glucose + unl oxygen
mu.ub(strcmp(mu.rxns, 'MAR09034_REV')) = 1; %Supply glucose - oxygen should already be open
mu.ub(strcmp(mu.rxns, 'MAR09135_REV')) = 0; %Do not supply lactate from the outside
mu.c = double(strcmp(mu.rxns, 'tot_atp_hydr'));
%run simulation
res = solveLP(mu,1);
res %-29.5036 - looks good

%look at lactate output and export in each slice
%

%table(baseModel.rxns,constructEquations(baseModel))

% Lactate export: MAR06048, import: MAR08515

neurExpSel = strcmp(mu.rxns, 'MAR06048');
neurImpSel = strcmp(mu.rxns, 'MAR08515');
astrExpSel = strcmp(mu.rxns, 'A_MAR06048');
astrImpSel = strcmp(mu.rxns, 'A_MAR08515');

neurGlycSel = strcmp(mu.rxns, 'MAR04394');
astrGlycSel = strcmp(mu.rxns, 'A_MAR04394');
neurMitSel = strcmp(mu.rxns, 'MAR04152');
astrMitSel = strcmp(mu.rxns, 'A_MAR04152');

glycNeur = res.x(neurGlycSel);
glycAstr = res.x(astrGlycSel);
mitoNeur = res.x(neurMitSel);
mitoAstr = res.x(astrMitSel);

glycATPConsNeur = glycNeur*2;
glycATPConsAstr = glycAstr*2;
mitoATPConsNeur = mitoNeur*14.75;
mitoATPConsAstr = mitoAstr*14.75;

%test
%(glycATPConsNeur + mitoATPConsNeur)/(1-fracA)%31.5974
%(glycATPConsAstr + mitoATPConsAstr)/fracA%30.9482 - lower - this is expected since the astrocytes use glycolysis, and therefore less enzyme mass


d = struct();
d.glycATPGenNeur = glycATPConsNeur;
d.glycATPGenAstr = glycATPConsAstr;
d.mitoATPGenNeur = mitoATPConsNeur;
d.mitoATPGenAstr = mitoATPConsAstr;

save(outFilename, 'd')

end

