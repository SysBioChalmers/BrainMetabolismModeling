cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
load('data/MinModel.mat')

%some temp test code:
%
%ihuman = load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat').ihuman;
%ecModel2 = CreateCompartmentECGEM(ihuman);%takes ~2 min



%neuronModelTempl = ecModel;
%neuronModelTempl = ecModel2;
neuronModelTempl = minModel;

%neuronModelTempl.S((length(neuronModelTempl.mets)-3):length(neuronModelTempl.mets),:)
%neuronModel = neuronModelTempl;

%We have a basic maintenance cost for all proteins, i.e. protein turnover. Call this B.
%At synapses, there is a transport cost associated with maintenance for all except the MT proteins. The total is then T*B, where T is typically 1.xx, it should be 
%roughly proportional to the amount of protein (we ignore that different proteins have different half-life times).
%The cytosolic and other proteins are also penalized by a utilization factor. To simplify the calculations, we assume that these cannot move (which is not entirely true)
%and that mitochondria have a perfect mobility and are always where they are needed, and thereby can be used up to 100% (not entirely true either). There must therefore be a spare 
%capacity for cytosolic and other, while not for mitochondria. The cost is therefore T*B/U, where U is the utilization factor (which is 1 for mitochondria)

%We ignore astrocytes for the first simulation, and instead assume that we have 'the rest of the body' that consumes 80% of the ATP and the neurons that consume 20% of the ATP. 
%The rest of the body only have the turnover maintenance B. 

%Add maintenance cost for all proteins (transport not included)

%block all uptake except what we need:
rxnsToKeepOpen = {'MAR09047_REV'; ... %H2O
              'MAR09072_REV'; ... %Pi
              'MAR09079_REV'; ... %H+
              'MAR09048_REV' %O2
    };%we close the rest

[inExchRxns,inExchRxnsInd] = getExchangeRxns(neuronModelTempl, 'in');
sel = ~ismember(inExchRxns, rxnsToKeepOpen);
%sum(~sel) %4, ok
%constructEquations(neuronModelTempl, neuronModelTempl.rxns(inExchRxnsInd(~sel)))%looks good
neuronModelTempl.ub(inExchRxnsInd(sel)) = 0;
%mBase.c = double(strcmp(mBase.rxns,'MAR03964'));%set objective function to ATP consumption

UCytos = [0.1;0.05;0.01];

B = 1;
TNeuron = 1.10;
UNeuronMito = 0.6;

TAstro = 1.10;
UAstroMito = 0.3;

%we make these small and not proportional to the energy consumption,
%since they should be able to produce only lactate at a certain utilization
%level
fracN = 0.02;
fracA = 0.005;

%rest of the body, rob:
robModel = addPenaltiesToModel(neuronModelTempl, B, 1, 1, 1);


%Loop through the different static utilization points and check if the 
%neurons and astrocytes are producing or consuming lactate
neuronLactateImport = nan(length(UCytos),1);
astroLactateImport = nan(length(UCytos),1);
%this take a few minutes to run
results = cell(length(UCytos),0);
for i = 1:length(UCytos)
    disp(i)
    neuronModel = addPenaltiesToModel(neuronModelTempl, B, TNeuron, UCytos(i), UNeuronMito);
    astroModel = addPenaltiesToModel(neuronModelTempl, B, TAstro, UCytos(i), UAstroMito);

    combModel = buildFullBrainModel(robModel,neuronModel, astroModel, fracN, fracA);
    %constructEquations(combModel, 'tot_ATP_hydr')%looks good
    %set objective
    combModel.c = double(strcmp(combModel.rxns, 'tot_ATP_hydr'));
    %give only glucose as input (oxygen + some other is already open)
    combModel.ub(strcmp(combModel.rxns, 'MAR09034_REV')) = 1;
    
    %run simulation
    res = solveLP(combModel,1);
    listExchBetweenCT(combModel, res)
    %extract lactate import as import - export
    neuronLactateImport(i) = res.x(strcmp(combModel.rxns, 'n_MAR09135_REV')) - res.x(strcmp(combModel.rxns, 'n_MAR09135'));
    astroLactateImport(i) = res.x(strcmp(combModel.rxns, 'a_MAR09135_REV')) - res.x(strcmp(combModel.rxns, 'a_MAR09135'));
    results{i} = res;
end
neuronLactateImport
%   -0.6211
%    0.0444
%    0.0445
astroLactateImport
%   -0.1553
%   -0.1597
%    0.0121

%investigate



%Now do the same again for point 1, but limit lactate uptake in Rob
%this take a few minutes to run
for i = 1:1
    disp(i)
    neuronModel = addPenaltiesToModel(neuronModelTempl, B, TNeuron, UCytos(i), UNeuronMito);
    astroModel = addPenaltiesToModel(neuronModelTempl, B, TAstro, UCytos(i), UAstroMito);

    combModel = buildFullBrainModel(robModel,neuronModel, astroModel, fracN, fracA);
    %constructEquations(combModel, 'tot_ATP_hydr')%looks good
    %set objective
    combModel.c = double(strcmp(combModel.rxns, 'tot_ATP_hydr'));
    %limit the lactate uptake in rob - this is a bit complicated, since there is 
    %not one reaction. The easiest is to add an extra metabolite to the export function
    %and add an exchange rxn on that which we can limit.
    metsToAdd.mets = {'LactateExportControl'};
    metsToAdd.metNames = {'LactateExportControl'};
    metsToAdd.compartments = {'s'};
    combModel = addMets(combModel, metsToAdd, false);
    combModel.S(numel(combModel.mets), ismember(combModel.rxns, {'n_MAR09135';'a_MAR09135'})) = 1;
    %constructEquations(combModel, 'n_MAR09135') %looks good
    rxnsToAdd.rxns = {'ExpLactateExportControl'};
    rxnsToAdd.equations = {'LactateExportControl[s] =>'};
    combModel = addRxns(combModel,rxnsToAdd, 3);
    %constructEquations(combModel, 'ExpLactateExportControl') %looks good
    %finally constrain it
    combModel.ub(strcmp(combModel.rxns, 'ExpLactateExportControl')) = 0.2;
    
    %give only glucose as input (oxygen + some other is already open)
    combModel.ub(strcmp(combModel.rxns, 'MAR09034_REV')) = 1;
    
    %run simulation
    res = solveLP(combModel,1);
    
    %extract lactate import as import - export
    neuronLactateImport2 = res.x(strcmp(combModel.rxns, 'n_MAR09135_REV')) - res.x(strcmp(combModel.rxns, 'n_MAR09135'));
    astroLactateImport2 = res.x(strcmp(combModel.rxns, 'a_MAR09135_REV')) - res.x(strcmp(combModel.rxns, 'a_MAR09135'));
    neuronGlucoseImport2 = res.x(strcmp(combModel.rxns, 'n_MAR09034_REV')) - res.x(strcmp(combModel.rxns, 'n_MAR09034'));
    astroGlucoseImport2 = res.x(strcmp(combModel.rxns, 'a_MAR09034_REV')) - res.x(strcmp(combModel.rxns, 'a_MAR09034'));
    %constructEquations(combModel, 'n_MAR09034_REV')
end

fracNeuronGlucExpAsLact = -neuronLactateImport2/(2*neuronGlucoseImport2)%0.5380
fracAstroGlucExpAsLact = -astroLactateImport2/(2*astroGlucoseImport2)%1
%{
neuronModel = addPenaltiesToModel(neuronModelTempl, B, TNeuron, UNeuronCyto, UNeuronMito);
astroModel = addPenaltiesToModel(neuronModelTempl, B, TAstro, UAstroCyto, UAstroMito);


combModel = buildFullBrainModel(robModel,neuronModel, astroModel, fracN, fracA);
constructEquations(combModel, 'tot_ATP_hydr')%looks good

%set objective
combModel.c = double(strcmp(combModel.rxns, 'tot_ATP_hydr'));
%give only glucose as input (oxygen + some other is already open)
combModel.ub(strcmp(combModel.rxns, 'MAR09034_REV')) = 1;

%combModel.c = double(strcmp(combModel.rxns, 'n_MAR03964'));


%run simulation
res = solveLP(combModel,1);
res

res.x(strcmp(combModel.rxns, 'n_MAR09135_REV'))%0.3650 - so, the neurons run on lactate here
res.x(strcmp(combModel.rxns, 'a_MAR09135'))%the astrocytes do export lactate in this setting
res.x(strcmp(combModel.rxns, 'a_MAR09135_REV'))%0 - the astrocytes do not import lactate

constructEquations(combModel, 'n_MAR03964')
constructEquations(combModel,'n_MAR09034_REV')










%In astrocytes, we assume the situation is similar. However, we assume that the mobility of the .


%}
