cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
ecModel = load('data/MinModel.mat').minModel;
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

glycProtPerATP = totProtUsageGlyc/ATPProdGlyc;%0.00311949499922544
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
hlData = ReadHalflifeData('C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling/data/pr101183k_si_002_HELA.txt');
hlData.GeneNames
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
cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling/'
save('data/HalfLives.mat', 's');



%%%%%%%%%%%%%%%%%%%%%%%
% Generate data for fig 2 B-D with the model
%%%%%%%%%%%%%%%%%%%%%%%

% Fig 2C first - need these results for Fig 2B
%%%%%%%%%%%%%%%%

util = (1:100)*0.01;
prodGlyc = NaN(length(util),1);
prodMito = NaN(length(util),1);
prodMitoMob = NaN(length(util),1);
EAMCAGlyc = NaN(length(util),1);
EAMCAMito = NaN(length(util),1);
EAMCAMitoMob = NaN(length(util),1);

%run glycolysis
for i = 1:length(util)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, 0, 0, util(i), util(i), false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  %'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
        };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only glucose as input
    modPen.ub(strcmp(modPen.rxns, 'MAR09034_REV')) = 1;

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;


    %run simulation
    res = solveLP(modPen,1);
    prodGlyc(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(util)
    EAMCAGlyc(i) = (prodGlyc(100) - prodGlyc(i))/prodGlyc(100);
end
log2(EAMCAGlyc)

%run mitochondrial respiration, static utilization
for i = 1:length(util)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, 0, 0, util(i), util(i), false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only lactate as input (oxygen + some other is already open)
    modPen.ub(strcmp(modUnpen.rxns, 'MAR09135_REV')) = 2; %two comparable to one glucose, doesn't really matter though

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;

    %run simulation
    res = solveLP(modPen,1);
    prodMito(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(util)
    EAMCAMito(i) = (prodMito(100) - prodMito(i))/prodMito(100);
end
log2(EAMCAMito)

util2 = util

util = (1:100)/100;

%run mitochondrial respiration, mitochondrial mobility
for i = 1:length(util)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, 0, 0, util(i), util(i) + 0.4*(1-util(i)), false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only lactate as input (oxygen + some other is already open)
    modPen.ub(strcmp(modUnpen.rxns, 'MAR09135_REV')) = 2; %two comparable to one glucose, doesn't really matter though

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;

    %run simulation
    res = solveLP(modPen,1);
    prodMitoMob(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(util)
    EAMCAMitoMob(i) = (prodMitoMob(100) - prodMitoMob(i))/prodMitoMob(100);
end
log2(EAMCAMitoMob)

% Fig 2B
%%%%%%%%%%%%%%%%

transp = (1:200)*0.01;
prodTrGlyc = NaN(length(transp),1);
prodTrMito = NaN(length(transp),1);
EAMCATrGlyc = NaN(length(transp),1);
EAMCATrMito = NaN(length(transp),1);

%run glycolysis
for i = 1:length(transp)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, transp(i), transp(i), 1, 1, false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  %'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
        };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only glucose as input
    modPen.ub(strcmp(modPen.rxns, 'MAR09034_REV')) = 1;

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;


    %run simulation
    res = solveLP(modPen,1);
    prodTrGlyc(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(transp)
    EAMCATrGlyc(i) = (prodGlyc(100) - prodTrGlyc(i))/prodGlyc(100); %We use the data prodGlyc(100) from 2C here, it has no transport or util costs
end
log2(EAMCATrGlyc)

%run mitochondrial respiration, static utilization
for i = 1:length(transp)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, transp(i), transp(i), 1, 1, false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only lactate as input (oxygen + some other is already open)
    modPen.ub(strcmp(modUnpen.rxns, 'MAR09135_REV')) = 2; %two comparable to one glucose, doesn't really matter though

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;

    %run simulation
    res = solveLP(modPen,1);
    prodTrMito(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(transp)
    EAMCATrMito(i) = (prodMito(100) - prodTrMito(i))/prodMito(100);%We use the data prodGlyc(100) from 2C here, it has no transport or util costs
end
log2(EAMCATrMito)

% Fig 2C first - need these results for Fig 2B
%%%%%%%%%%%%%%%%

util2 = (1:1000)*0.001;
prodGlycD = NaN(length(util2),1);
prodMitoMobD = NaN(length(util2),1);
EAMCAGlycD = NaN(length(util2),1);
EAMCAMitoMobD = NaN(length(util2),1);

%run glycolysis
for i = 1:length(util2)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, 0.1, 0.1, util2(i), util2(i), false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  %'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
        };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only glucose as input
    modPen.ub(strcmp(modPen.rxns, 'MAR09034_REV')) = 1;

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;


    %run simulation
    res = solveLP(modPen,1);
    prodGlycD(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(util2)
    EAMCAGlycD(i) = (prodGlycD(1000) - prodGlycD(i))/prodGlycD(1000);
end
log2(EAMCAGlycD)

%run mitochondrial respiration, mitochondrial mobility
for i = 1:length(util2)
    %model, B, TM, TG, UCyto, UMito, dontPenalize
    modPen = addPenaltiesToModel(ecModel, 0.1, 0.1, 0.1, util2(i), util2(i) + 0.4*(1-util2(i)), false);

    %now feed the model with some basics + glucose, no oxygen:

    rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
                  'MAR09072_REV'; ... %Pi
                  'MAR09079_REV'; ...%H+
                  'MAR09048_REV'; %O2
                  'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt' ...
    };%we close the rest

    [inExchRxns,inExchRxnsInd] = getExchangeRxns(modPen, 'in');
    sel = ~ismember(inExchRxns, rxnsToKeepOpen);
    %sum(~sel) %7, ok
    %constructEquations(modUnpen, modUnpen.rxns(inExchRxnsInd(~sel)))%looks good
    modPen.ub(inExchRxnsInd(sel)) = 0;

    %give only lactate as input (oxygen + some other is already open)
    modPen.ub(strcmp(modUnpen.rxns, 'MAR09135_REV')) = 2; %two comparable to one glucose, doesn't really matter though

    %optimize for ATP hydrolysis, and give a small penalty to using protein
    %to avoid meaningless fluxes
    modPen.c = double(strcmp(modPen.rxns, 'MAR03964'));
    modPen.c(ismember(modPen.rxns, {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'})) = -0.001;

    %run simulation
    res = solveLP(modPen,1);
    prodMitoMobD(i) = res.x(strcmp(modPen.rxns, 'MAR03964'));
end

%the one to compare with is the one with full util, so subtract from the util(100) to get the diff
for i = 1:length(util2)
    EAMCAMitoMobD(i) = (prodMitoMobD(1000) - prodMitoMobD(i))/prodMitoMobD(1000);
end
log2(EAMCAMitoMobD)

%For practical reasons, we had to reduce the B to 0.1 in these simulations, since the EAMCA
%sometimes become larger than 1 (for mitochondrial respiration for small utilizations), and
%we therefore cannot produce any ATP. This is a modeling issue only, most of the maintenance
%energy is not spent in the same cell in that sense. It doesn't happen in the other modeling
%scenarios either, since we don't use static utilization with mito.
s = struct();
s.util = util;
s.util2 = util2;
s.transp = transp;
s.b_EAMCAGlyc = EAMCATrGlyc*10;
s.b_EAMCAMito = EAMCATrMito*10;
s.c_EAMCAGlyc = EAMCAGlyc*10;
s.c_EAMCAMito = EAMCAMito*10;
s.c_EAMCAMitoMob = EAMCAMitoMob*10;
s.d_EAMCAGlyc = EAMCAGlycD*10;
s.d_EAMCAMitoMob = EAMCAMitoMobD*10;
save('data/2B-DData.mat','s')







