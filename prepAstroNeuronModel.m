cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ihuman = load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat').ihuman;
modSimp = simplifyModel(ihuman,false,false,true,true); %remove dead-end reactions


%figure out which reactions that can carry flux

rxnsToKeepOpen = {'MAR09135'; ... Lactate
              'MAR09034'; ... %Glucose
              'MAR09047'; ... %H2O
              'MAR09072'; ... %Pi
              'MAR09079'; ... %H+
              'MAR09048'; ... %O2
};

[inExchRxns,inExchRxnsInd] = getExchangeRxns(modSimp);
sel = ~ismember(inExchRxns, rxnsToKeepOpen);
%sum(~sel) %6, ok
constructEquations(modSimp, modSimp.rxns(inExchRxnsInd(~sel)))%looks good
modSimpSim = modSimp;
modSimpSim.lb(inExchRxnsInd(sel)) = 0; %remove import of everything except glucose etc.

%Check which reactions cannot carry flux - these are later to be removed. This takes a while to run (a few hours)
hfRes = haveFlux(modSimpSim,10^-14); %12053 reactions before
save('data/hfres.mat', 'hfRes')
sum(hfRes)%10844
minModel = removeReactions(modSimp,~hfRes, true, true, false); %don't remove compartments, could perhaps lead to problems later.

ecModel = CreateCompartmentECGEM(minModel);%takes ~2 min


ecModel.lb(ecModel.lb == -1000) = -Inf; %These operations help the solver, it runs faster with inf
ecModel.ub(ecModel.ub == 1000) = Inf;

save('data/ecModel.mat', 'ecModel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now make a smaller model only including the central carbon metabolism reactions that will be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
load('data/ecModel.mat')


%We run 3 simulations where we optimize for growth:
% 1. Glucose only - lactate export, pure glycolysis
% 2. Glucose + oxygen - full oxidation of glucose, both pathways
% 3. Lactate + oxygen - lactate import, mitochondrial respiration
%
% We then only keep the reactions that carry flux in either of these 3 simulations.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start with glycolysis, i.e. feed on glucose, no oxygen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%open up the flux for protein usage, but don't penalize it
modUnpen = addPenaltiesToModel(ecModel, 1, 1, 1, 1, 1, true);


%now feed the model with some basics + glucose, no oxygen:

rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
              'MAR09072_REV'; ... %Pi
              'MAR09079_REV'; ...%H+
              %'MAR09048_REV'; %O2
              'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

[inExchRxns,inExchRxnsInd] = getExchangeRxns(modUnpen, 'in');
sel = ~ismember(inExchRxns, rxnsToKeepOpen);
%sum(~sel) %7, ok
%constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
modUnpen.ub(inExchRxnsInd(sel)) = 0;

%give only glucose as input
modUnpen.ub(strcmp(modUnpen.rxns, 'MAR09034_REV')) = 1;

%optimize for ATP hydrolysis, and give a small penalty to using protein
%to avoid meaningless fluxes
modUnpen.c = double(strcmp(modUnpen.rxns, 'MAR03964'));
modUnpen.c(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;


%run simulation
res = solveLP(modUnpen,1);
res

selGlyc = res.x ~= 0;

table(modUnpen.rxns(selGlyc), res.x(selGlyc), constructEquations(modUnpen, modUnpen.rxns(selGlyc)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Full oxidation of glucose, i.e. feed on glucose + oxygen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open up the flux for protein usage, but don't penalize it
modUnpen = addPenaltiesToModel(ecModel, 1, 1, 1, 1, 1, true);

rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
              'MAR09072_REV'; ... %Pi
              'MAR09079_REV'; ...%H+
              'MAR09048_REV'; %O2
              'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

[inExchRxns,inExchRxnsInd] = getExchangeRxns(modUnpen, 'in');
sel = ~ismember(inExchRxns, rxnsToKeepOpen);
%sum(~sel) %8, ok
%constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
modUnpen.ub(inExchRxnsInd(sel)) = 0;

%give glucose as input (oxygen + some other is already open)
modUnpen.ub(strcmp(modUnpen.rxns, 'MAR09034_REV')) = 1;

%optimize for ATP hydrolysis, and give a small penalty to using protein
%to avoid meaningless fluxes
modUnpen.c = double(strcmp(modUnpen.rxns, 'MAR03964'));
modUnpen.c(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;


%run simulation
res = solveLP(modUnpen,1);
res

selBoth = res.x ~= 0;

table(modUnpen.rxns(selBoth), res.x(selBoth))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. Now for mitochondria, i.e. feed on lactate + oxygen, no glucose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open up the flux for protein usage, but don't penalize it
modUnpen = addPenaltiesToModel(ecModel, 1, 1, 1, 1, 1, true);

rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
              'MAR09072_REV'; ... %Pi
              'MAR09079_REV'; ...%H+
              'MAR09048_REV'; %O2
              'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

[inExchRxns,inExchRxnsInd] = getExchangeRxns(modUnpen, 'in');
sel = ~ismember(inExchRxns, rxnsToKeepOpen);
%sum(~sel) %8, ok
%constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
modUnpen.ub(inExchRxnsInd(sel)) = 0;

%give only lactate as input (oxygen + some other is already open)
modUnpen.ub(strcmp(modUnpen.rxns, 'MAR09135_REV')) = 2; %two comparable to one glucose, doesn't really matter though

%optimize for ATP hydrolysis, and give a small penalty to using protein
%to avoid meaningless fluxes
modUnpen.c = double(strcmp(modUnpen.rxns, 'MAR03964'));
modUnpen.c(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;


%run simulation
res = solveLP(modUnpen,1);
res

selMito = res.x ~= 0;

table(modUnpen.rxns(selMito), res.x(selMito))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now keep all reactions that are used in any of the cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selAll = selGlyc | selBoth | selMito;
sum(selAll) %59

%so, we don't want to keep the unpenalized model, because it already has reactions
%for supplying enzymatic capacity. We instead minimize the ecModel directly

rxnIds = modUnpen.rxns(selAll);

ecModelSel = ismember(ecModel.rxns, rxnIds);
sum(ecModelSel) %55, as expected - the 4 enzymatic capacity reactions from penalizeModel are not there, so 4 less here

minModelTmp = removeReactions(ecModel, ~ecModelSel, true, true, true);

%For some reason the OriginalmetNames is still large and there - we don't need that, remove it.
minModel = rmfield(minModelTmp, 'OriginalmetNames')

save('data/MinModel.mat', 'minModel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test minModel with the three conditions again, just to be sure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


modUnpen = addPenaltiesToModel(minModel, 1, 1, 1, 1, 1, true);

%now feed the model with some basics + glucose, no oxygen:

rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
              'MAR09072_REV'; ... %Pi
              'MAR09079_REV'; ...%H+
              %'MAR09048_REV'; %O2
              'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

[inExchRxns,inExchRxnsInd] = getExchangeRxns(modUnpen, 'in');
sel = ~ismember(inExchRxns, rxnsToKeepOpen);
sum(~sel) %4, ok, only the prot pool rxns
constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(sel)))%looks good, glucose, O2, lactate
modUnpen.ub(inExchRxnsInd(sel)) = 0;

%optimize for ATP hydrolysis, and give a small penalty to using protein
%to avoid meaningless fluxes
modUnpen.c = double(strcmp(modUnpen.rxns, 'MAR03964'));
modUnpen.c(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;

%glucose only:
%%%%%%%%%%%%%%

mu1 = modUnpen;
%give only glucose as input
mu1.ub(strcmp(mu1.rxns, 'MAR09034_REV')) = 1;
%run simulation
res = solveLP(mu1,1);
res %2, ok

%both:
%%%%%%%%%%%%%%

mu2 = modUnpen;
%glucose + unl oxygen
mu2.ub(strcmp(mu2.rxns, 'MAR09034_REV')) = 1;
mu2.ub(strcmp(mu2.rxns, 'MAR09048_REV')) = Inf;
%run simulation
res = solveLP(mu2,1);
res % ~31.5, ok

%lactate:
%%%%%%%%%%%%%%

mu3 = modUnpen;
%2 lactate + unl oxygen
mu3.ub(strcmp(mu3.rxns, 'MAR09135_REV')) = 2;
mu3.ub(strcmp(mu3.rxns, 'MAR09048_REV')) = Inf;
%run simulation
res = solveLP(mu3,1);
res % ~29.5, ok

%all ok!