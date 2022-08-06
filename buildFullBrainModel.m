function [outModel] = buildFullBrainModel(robModel, neuronModel, astroModel, fracN, fracA)
% buildFullBrainModel
%
% Generates a combined model with neurons, astrocytes, and 'the rest of the body' (rob).
% The model supports a total tumor biomass objective function which builds
% ECM and tumor cells. The other cell types are assumed to be imported into
% the tumor.
%
% Input:
%
%   robModel     Should be an ec model derived from Human-GEM
%
%   fracC           fracCion of cancer cells in the model
%
%   fracF           fracCion of fibroblasts in the model. fracCion other cells 
%                   is calculated as 1 - fracC - fracF
%
%   fracECM         fracCion extracellular matrix in the tumor. The rest is
%                   cells
%
%   fracGAGsInECM   fracCion of GAGs in the ECM. The rest is protein.
%
% Output:
%
%   outModel        The full model
%

%for debugging
%load('C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/ecHumanGEM_batch.mat')
%ecModelOrig = ecModel_batch;

%length(getExchangeRxns(robModel)) %3347
%constructEquations(robModel, getExchangeRxns(robModel))
%sum(getMetsInComp(robModel, 'S'))
%robModel.metNames(getMetsInComp(robModel, 'S'))

%robModel.metComps

%figure out which metabolite is the protein pool, i.e. in the
%'prot_pool_exchange' reaction
%constructEquations(robModel, robModel.rxns(length(robModel.rxns))) %{' => prot_pool[c]'}
%So, this is not an exchange flux from the s compartment, which means we do not need to have any
%special treatment of this when copying the model

%The strategy for creating the three cell types combined model is to
%1. create 3 copies of all metabolites
%2. create 3 copies of all reactions
%3. In the matrix, the second copy of the reactions should have the second
%copy of the metabolites. Same for the third copy.





%Identify metabolites to copy
metsToCopy = true(length(robModel.mets),1);
%identify rxns to copy
rxnsToCopy = true(1,length(robModel.rxns));

model2 = robModel;

%generate new compartments
nmsN = strcat('n_', neuronModel.comps);
nmsA = strcat('a_', astroModel.comps);
model2.comps = [robModel.comps;nmsN;nmsA];
model2.compNames = [robModel.compNames;strcat('n_',robModel.compNames);strcat('a_',robModel.compNames)];

%extend the metabolites
nms = robModel.mets(metsToCopy);
nmsN = strcat('n_', neuronModel.mets(metsToCopy));
nmsA = strcat('a_', astroModel.mets(metsToCopy));
model2.mets = [robModel.mets; nmsN; nmsA];
model2.metNames = [robModel.metNames; neuronModel.metNames(metsToCopy); astroModel.metNames(metsToCopy)];
%met compartments
metComps = robModel.metComps(metsToCopy);

toAdd = length(robModel.comps);
compsN = neuronModel.metComps(metsToCopy) + toAdd;
compsA = astroModel.metComps(metsToCopy) + toAdd.*2;
%unique(compsN)
%unique(compsA) %looks good
model2.metComps = [robModel.metComps;compsN;compsA];
%b
model2.b = [robModel.b;neuronModel.b(metsToCopy);astroModel.b(metsToCopy)];

%now, the reactions
nmsN = strcat('n_', neuronModel.rxns(rxnsToCopy));
nmsA = strcat('a_', astroModel.rxns(rxnsToCopy));
model2.rxns = [robModel.rxns; nmsN; nmsA];
nmsN = strcat('n_', neuronModel.rxnNames(rxnsToCopy));
nmsA = strcat('a_', astroModel.rxnNames(rxnsToCopy));
model2.rxnNames = [robModel.rxnNames; nmsN; nmsA];

%lb, ub, rev and c
model2.lb = [robModel.lb;neuronModel.lb(rxnsToCopy);astroModel.lb(rxnsToCopy)];
model2.ub = [robModel.ub;neuronModel.ub(rxnsToCopy);astroModel.ub(rxnsToCopy)];
model2.rev = [robModel.rev;neuronModel.rev(rxnsToCopy);astroModel.rev(rxnsToCopy)];
model2.c = [robModel.c;neuronModel.c(rxnsToCopy);astroModel.c(rxnsToCopy)];

%grRules subSystems and eccodes
model2.grRules = [robModel.grRules;neuronModel.grRules(rxnsToCopy);astroModel.grRules(rxnsToCopy)];
model2.subSystems = [robModel.subSystems;neuronModel.subSystems(rxnsToCopy);astroModel.subSystems(rxnsToCopy)];
model2.eccodes = [robModel.eccodes;neuronModel.eccodes(rxnsToCopy);astroModel.eccodes(rxnsToCopy)];

%rxnGeneMat
model2.rxnGeneMat = [robModel.rxnGeneMat;neuronModel.rxnGeneMat(rxnsToCopy,:);astroModel.rxnGeneMat(rxnsToCopy,:)];

%And now finally the S matrix:
%first just add zeros in both directions
model2.S(length(model2.mets), length(model2.rxns)) = 0;
[nm,nr] = size(robModel.S);
model2.S(1:nm,1:nr) = robModel.S;
model2.S((1:nm)+nm,(1:nr)+nr) = neuronModel.S;
model2.S((1:nm)+2*nm,(1:nr)+2*nr) = astroModel.S;

%Now, redirect the exchange reactions. The uptake reactions of both a and n
%should take up everything from the S compartment, not from nothing. They
%should also be unconstrained
rxnSelRob = 1:nr;
rxnSelN = (1:nr) + nr;
rxnSelA = (1:nr) + nr*2;
metSelRob = (1:nm).';
metSelN = (1:nm).' + nm;
metSelA = (1:nm).' + nm*2;

%now, we change the exchange rxns of glucose and lactate in neurons and astrocytes to instead
%go to the rob model - it is then in the rob model that we manipulate the outputs.
sComp = find(strcmp(model2.comps, 's'));
%glucose
%import
glucoseSel = find(strcmp(model2.metNames,'glucose') & model2.metComps == sComp);
model2.S(glucoseSel,strcmp(model2.rxns, 'n_MAR09034_REV')) = -1;
model2.ub(strcmp(model2.rxns, 'n_MAR09034_REV')) = Inf;
%constructEquations(model2,'n_MAR09034_REV')%looks good
model2.S(glucoseSel,strcmp(model2.rxns, 'a_MAR09034_REV')) = -1;
model2.ub(strcmp(model2.rxns, 'a_MAR09034_REV')) = Inf;
%export
model2.S(glucoseSel,strcmp(model2.rxns, 'n_MAR09034')) = 1;
model2.ub(strcmp(model2.rxns, 'n_MAR09034')) = Inf;
%constructEquations(model2,'n_MAR09034')%looks good
model2.S(glucoseSel,strcmp(model2.rxns, 'a_MAR09034')) = 1;
model2.ub(strcmp(model2.rxns, 'a_MAR09034')) = Inf;

%lactate
%import
lactateSel = find(strcmp(model2.metNames,'L-lactate') & model2.metComps == sComp);
model2.S(lactateSel,strcmp(model2.rxns, 'n_MAR09135_REV')) = -1;
model2.ub(strcmp(model2.rxns, 'n_MAR09135_REV')) = Inf;
%constructEquations(model2,'n_MAR09135_REV')%looks good
model2.S(lactateSel,strcmp(model2.rxns, 'a_MAR09135_REV')) = -1;
model2.ub(strcmp(model2.rxns, 'a_MAR09135_REV')) = Inf;
%export
model2.S(lactateSel,strcmp(model2.rxns, 'n_MAR09135')) = 1;
model2.ub(strcmp(model2.rxns, 'n_MAR09135')) = Inf;
%constructEquations(model2,'n_MAR09135')%looks good
model2.S(lactateSel,strcmp(model2.rxns, 'a_MAR09135')) = 1;
model2.ub(strcmp(model2.rxns, 'a_MAR09135')) = Inf;



%% Now add the specialized objective function

fracRob = 1 - fracN - fracA;

rxn = [num2str(fracRob), ' ATP[c] + ', num2str(fracRob), ' H2O[c] + ', ...
       num2str(fracN), ' ATP[n_c] + ', num2str(fracN), ' H2O[n_c] + ', ...
       num2str(fracA), ' ATP[a_c] + ', num2str(fracA), ' H2O[a_c] => ', ...
       num2str(fracRob), ' ADP[c] + ', num2str(fracRob), ' Pi[c] + ', num2str(fracRob), ' H+[c] + ', ...
       num2str(fracN), ' ADP[n_c] + ', num2str(fracN), ' Pi[n_c] + ', num2str(fracN), ' H+[n_c] + ', ...
       num2str(fracA), ' ADP[a_c] + ', num2str(fracA), ' Pi[a_c] + ', num2str(fracA), ' H+[a_c]'];
   

rxnsToAdd = struct();   
rxnsToAdd.rxns = {'tot_ATP_hydr'};
rxnsToAdd.equations = {rxn};
model2 = addRxns(model2,rxnsToAdd, 3);

outModel = model2;

