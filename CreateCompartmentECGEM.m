function outModel = CreateCompartmentECGEM(model, fillInMissingGPRs, minAcceptableKCat, speciesAdapter)

%This code was copied from CreateECLightModel in GeckoLight and modified. What we do is to divide the costs into 
%4 categories: cytosolic, mitochondria, MT-genes, and other.
%Cytosolic reactions are identified as reactions only handling cytosolic metabolites.
%Mitochondrial as having any metabolite in the m compartment
%The MT part is a bit trickier. For selection of which isozyme that have the lowest cost, we ignore MT.
%Then, we just take the fraction of the modelcular weight that belongs to MT genes.
%So, we produce 4 metabolites instead of 1 (compared to Gecko), but we don't create any 
%reaction to deal with them, that is to be customized outside this function

if nargin < 4
	speciesAdapter = HumanGEMAdapter();
end
if nargin < 3
	minAcceptableKCat = 1; %This is what we ususally use
end
if nargin < 2
	fillInMissingGPRs = true;
end

%pr101183k_si_002.xls
%for test only:
%load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat')
%model = ihuman;

%Get model-specific parameters
parameters = speciesAdapter.getParameters();

oldMetNames = model.metNames;

%Remove blocked rxns + correct model.rev:
[model,name,modelVer] = preprocessModelLt(model,'ecLight','1');

fprintf('\n==================')
fprintf('\nGenerating light ecModel:')
fprintf('\n==================\n')

model_data = getEnzymeCodesOpt(model, speciesAdapter); %Elapsed time is 38.600622 seconds, I think this could be acceptable.
fprintf('\n')


%Load BRENDA data:
[KCATcell, SAcell] = loadBRENDAdataLt(speciesAdapter);
kcatsRes      = matchKcatsOpt(model_data,parameters.org_name, KCATcell, SAcell, speciesAdapter, minAcceptableKCat);
fprintf('\n')
%code here copied from readKcatData

%Get kcat value for both directions:
Fkcat = kcatsRes.forw.kcats;
Bkcat = kcatsRes.back.kcats;
rev   = logical(model_data.model.rev);
kcats = [Fkcat;Bkcat(rev,:)];
wcLevels = [kcatsRes.forw.wcLevel;kcatsRes.back.wcLevel(rev,:)];

%Update uniprots with both directions:
uniprots = [model_data.uniprots; model_data.uniprots(rev,:)];

%Update matched genes with both directions:
matchedGenes = [model_data.matchedGenes; model_data.matchedGenes(rev,:)];

%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
model = convertToIrrev(model_data.model);


%What we want to find for each reaction
%is the min value of MW/kcat when comparing all parallel reactions (i.e. OR in the GR rules)
%The plan is then to add one reaction with prot_pool_exchange, one
%metabolite called prot_pool, and a stochiometric coefficient of MW/kcat of
%that for each reaction. For complexes, the MW are summed up, since the kcat is the same.

%Load databases:
data      = load(speciesAdapter.getFilePath('ProtDatabase.mat')); %this is loaded twice, (also in a function above), a bit unnecessary
swissprot = data.swissprot;

standardMW = median([swissprot{:,5}])/1000;%the median of all the proteins

%profile on -history
%profile off
%p = profile('info')
%[srt,I] = sort([p.FunctionTable.TotalTime],2, 'ascend');
%fnksTmp = {p.FunctionTable.FunctionName};
%fnks = fnksTmp(I);
%fnks{54}
%srt(54)

fprintf('Creating cytosol-penalized model')

mtGenes = {'ENSG00000198888'; ... %MT-ND1
           'ENSG00000198763'; ... %MT-ND2
           'ENSG00000198840'; ... %MT-ND3
           'ENSG00000212907'; ... %MT-ND4L
           'ENSG00000198886'; ... %MT-ND4
           'ENSG00000198786'; ... %MT-ND5
           'ENSG00000198695'; ... %MT-ND6
           'ENSG00000198727'; ... %MT-CYB
           'ENSG00000198804'; ... %MT-CO1
           'ENSG00000198712'; ... %MT-CO2
           'ENSG00000198938'; ... %MT-CO3
           'ENSG00000198899'; ... %MT-ATP6
           'ENSG00000228253'; ... %MT-ATP8
           'ENSG00000210082' ... %MT-RNR2
    };


%Figure out which reactions are purely cytosolic, i.e. only contains metabolites
%in the c compartment
cytComp = find(strcmp(model.comps,'c'));
cytosolicMets = model.metComps == cytComp;
onlyCytS = model.S;
onlyCytS(~cytosolicMets, :) = 0;
diffS = model.S - onlyCytS;
pureCytRxns = full(sum(diffS ~= 0, 1) == 0).';
%sum(pureCytRxns)
%constructEquations(model, model.rxns(pureCytRxns))

mitoComp = find(strcmp(model.comps,'m'));
mitoMets = model.metComps == mitoComp;
mitoRxns  = full(sum(model.S(mitoMets,:) ~= 0, 1) ~= 0).';
%constructEquations(model, model.rxns(mitoRxns))%looks good

[MWDivKcatsPerCategory,rxnWcLevels] = GetMWAndKcatsNeurons(uniprots,kcats,wcLevels, swissprot, standardMW, mtGenes, pureCytRxns, mitoRxns);
%test the mt part
sel = ~isnan(MWDivKcatsPerCategory(:,4)) & MWDivKcatsPerCategory(:,4) > 0;
sum(sel)%5, looks good, probably the complexes
constructEquations(model, model.rxns(sel))%the complexes, one is duplicated (with ROS), so complex II is not here, ok.
MWDivKcatsPerCategory(sel,:) %ok, all of these are in MT or mito.



fprintf('\n');

%now restore the metabolite names - they were temporarily changed in preprocessModelLt to match the entries in the enzyme databases
%this messes up some things, for example the names of the metabolites in the biomass function
model.metNames = oldMetNames;


%so, what we do now is to add the reaction prot_pool_exchange and the metabolite prot_pool
%rxnsToAdd = struct();
%rxnsToAdd.rxns = {'prot_pool_exchange'};
%rxnsToAdd.equations = {[num2str(ATPPenalty) ' ATP[c] => ' num2str(ATPPenalty) ' ADP[c] + prot_pool[c]']};
%rxnsToAdd.rxnNames = {'prot_pool_exchange'};
%rxnsToAdd.subSystems = {'GeckoLight Rxns'};
metsToAdd = struct();
metsToAdd.mets = {'other_prot_pool';'cyto_prot_pool';'mito_prot_pool';'mt_prot_pool'};
metsToAdd.metNames = {'other_prot_pool';'cyto_prot_pool';'mito_prot_pool';'mt_prot_pool'};
metsToAdd.compartments = {'c';'c';'c';'c'};

model = addMets(model, metsToAdd);
otherMetIndex = length(model.mets) - 3;
cytoMetIndex = length(model.mets) - 2;
mitoMetIndex = length(model.mets) - 1;
 





%Now add the protein metabolites to the reactions
metRows = MWDivKcatsPerCategory.';
%metRows(~pureCytRxns) = metRow(~pureCytRxns)*(fixedPenalty + otherTranspPen);
%metRows(pureCytRxns) = metRow(pureCytRxns)*(fixedPenalty + cytTranspPen);
metRows(isnan(metRows)) = 0;
model.S((length(model.mets)-3):length(model.mets),:) = -metRows;


%set the protein pool constraint
%model.ub(strcmp(model.rxns, 'prot_pool_exchange')) = poolConstraint;%;2.238315e-02; %the ecModels value was 0.0505717472, but this was fit using a model with the bug in the biomass reaction.

%standardRxnProtCost = median(MWDivKcats(~isnan(MWDivKcats)));%1.2894e-04

%standardRxnProtCostOther = median(metRows(~pureCytRxns));
%standardRxnProtCostCyt = median(metRow(pureCytRxns));


%The enzymes based on wildcard matches are uncertain - so make sure the
%protein cost is not unrealistically high - a single such value could dominate a whole
%simulation.
%The strategy is to assume that none of those reactions should have a
%protein cost higher than that of complex I. Replace all of those with a
%standard protein cost

%look at the complexes in oxphos:
%MWDivKcats(strcmp(model.rxns,'HMR_6921')) %complex I: 0.1054
%MWDivKcats(strcmp(model.rxns,'HMR_6918')) %complex III: 0.0011
%MWDivKcats(strcmp(model.rxns,'HMR_6914')) %complex IV: 0.0264
%MWDivKcats(strcmp(model.rxns,'HMR_6916')) %complex V: 7.0032e-04

%we assume that none of the reactions based on a wildcard should have a
%higher protein cost than human complex I - use a value of 0.1 as max
%realistic protein cost - if above that it means that it is totally
%unrealiable - we therefore take a standard value (and do not assign 0.1)
%sel = MWDivKcats > 0.1 & rxnWcLevels > 0;
%model.rxns(sel)
%constructEquations(model, model.rxns(sel))
%if (sum(sel > 0))
%    model.S(length(model.mets),[sel;false]) = -standardRxnProtCost;
%end

%So, this step was replaced by a change in matchKcatsOpt, where all kcats
%smaller than 1 s^-1 (3600 h^-1) are set to 1 s^-1. It works better with
%the stoichiometry.

if fillInMissingGPRs
    %Now comes the task of filling in a standard protein cost for reactions
    %with missing GPRs

    tmpMWDiv = MWDivKcatsPerCategory.';
    standardRxnProtCost = median(sum(metRows(:,~isnan(tmpMWDiv(1,:))), 1));

    %get info about spontanoeus rxns
    [spont,spontRxnNames] = speciesAdapter.getSpontaneousReactions(model);


    [~,exchangeRxnsIndexes] = getExchangeRxns(model);

    numToFix = 0;
%    protPoolIndex = find(strcmp(model.mets, 'prot_pool'));

    for i = 1:(length(MWDivKcatsPerCategory(:,1))-1)%avoid the protein_pool_exchange
       if isnan(MWDivKcatsPerCategory(i,1))
           %Step 1: Skip exchange reactions
           if ismember (i,exchangeRxnsIndexes)
               %disp(strcat(num2str(i),': skipping exchange: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
               continue;
           end
           %Step 2: Skip transport reactions (with missing GPRs)
           %We really only skip the reactions with two metabolites
           %where they are the same but different compartments
           numMets = sum(model.S(:,i) ~= 0);
           if numMets == 2
              mets = model.metComps(model.S(:,i) ~= 0);
              metNames = model.metNames(model.S(:,i) ~= 0);
              if (~strcmp(mets(1), mets(2))) %different compartments
                  if strcmp(metNames{1}, metNames{2}) %the same metabolite, it is a transport reaction
                      %disp(strcat(num2str(i),': skipping transport: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
                      continue;
                  end
              end
           end
           %Step 3: check if this is a spontaneous reaction (last because it is a bit heavier than the other checks)
           rxnName = model.rxns{i};
           if endsWith(rxnName, '_REV')
              rxnName = extractBefore(rxnName, strlength(rxnName)-3);
           end
           spontaneous = spont(strcmp(spontRxnNames, rxnName)); %could potentially be optimized
           if ~isempty(spontaneous)%it shouldn't be
               if ~isnan(spontaneous)%we treat NaN as that it is unknown if the reaction is spontaneous, which is seen as a no
                   if spontaneous
                      %disp(strcat(num2str(i),': skipping spontaneous: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
                      continue; %leave spontaneous reactions without protein cost
                   end
               end
           end
           numToFix = numToFix + 1;
           %disp(strcat(num2str(i),': Assigning std val: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
           %fill in a standard value
           if pureCytRxns(i)
               model.S(cytoMetIndex, i) = -standardRxnProtCost;
           elseif mitoRxns(i)
               model.S(mitoMetIndex, i) = -standardRxnProtCost;
           else
               model.S(otherMetIndex, i) = -standardRxnProtCost;
           end
       end
    end
    
    disp(['Filled in ' num2str(numToFix) ' protein usage costs with a standard value due to missing GPRs or missing kcats.'])
end


outModel = model;

end