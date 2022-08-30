cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ecModel = load('data/ecModel.mat').ecModel;
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

sel = res.x ~= 0;

table(modUnpen.rxns(sel), res.x(sel))

ATPProdGlyc = res.x(strcmp(modUnpen.rxns, 'MAR03964'));
totProtUsageGlyc = sum(res.x(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito'})));%we don't count the mt genes
fluxesGlyc = res.x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now for mitochondria, i.e. feed on lactate + oxygen, no glucose
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

ATPProdMito = res.x(strcmp(modUnpen.rxns, 'MAR03964'));
totProtUsageMito = sum(res.x(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito'})));%we don't count the mt genes
totProtUsageMitoInclMT = sum(res.x(ismember(modUnpen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})));
fluxesMito = res.x;

%%%%%%%%%%%%%%%%%
% Analyze
%%%%%%%%%%%%%%%%%

glycProtPerATP = totProtUsageGlyc/ATPProdGlyc;%0.0031
mitoProtPerATP = totProtUsageMito/ATPProdMito;%0.0329
mitoProtPerATPWithMT = totProtUsageMitoInclMT/ATPProdMito;%0.0434
ATPProdGlyc % 2
ATPProdMito % 29.5

ratio = mitoProtPerATP / glycProtPerATP %10.5520
ratioWithMT = mitoProtPerATPWithMT / glycProtPerATP %13.9088

%%%%%%%%%%%%%%%%%
% Investigate half-lives of the enzymes used for the fluxes
%%%%%%%%%%%%%%%%%
selGlyc = fluxesGlyc ~= 0 & ~strcmp(modUnpen.grRules,'');

[convGrRules, convGenes, convGeneRxnMatrix] = translateGrRules(modUnpen.grRules, 'Name');
grRulesGlyc = convGrRules(selGlyc);
grRulesGlyc
selGenesGlyc = sum(convGeneRxnMatrix(selGlyc,:) ~= 0, 1) > 0;
glycGenes = convGenes(selGenesGlyc);

selMito = fluxesMito ~= 0 & ~strcmp(modUnpen.grRules,'');
grRulesMito = convGrRules(selMito);
grRulesMito
selGenesMito = sum(convGeneRxnMatrix(selMito,:) ~= 0, 1) > 0;
mitoGenesTmp = convGenes(selGenesMito);
%filter out MT-genes
mitoGenes = mitoGenesTmp(~startsWith(mitoGenesTmp, 'MT-'));
hlData = ReadHalflifeData('C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Single-cell Modeling/NeuronModeling/data/pr101183k_si_002_HELA.txt');
hlData.GeneNames
split('A;B;APA',';')
splitGeneNames = cell(length(hlData.GeneNames),1);
for i = 1:length(hlData.GeneNames)
    splitGeneNames{i} = split(hlData.GeneNames(i), ';');
end

%just extract all rows where the gene name exists
glycRows = false(length(hlData.GeneNames),1);
for i = 1:length(hlData.GeneNames)
    glycRows(i) = any(ismember(glycGenes, splitGeneNames{i}));
end
sum(glycRows)
hlData.GeneNames(glycRows)
glycVals = hlData.half_lifeT1_2InH(glycRows);

mitoRows = false(length(hlData.GeneNames),1);
for i = 1:length(hlData.GeneNames)
    mitoRows(i) = any(ismember(mitoGenes, splitGeneNames{i}));
end
sum(mitoRows)
hlData.GeneNames(mitoRows)
mitoVals = hlData.half_lifeT1_2InH(mitoRows);

s = struct();
s.glycVals = glycVals;
s.mitoVals = mitoVals;
cd 'C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Single-cell Modeling/NeuronModeling/'
save('data/HalfLives.mat', 's');



