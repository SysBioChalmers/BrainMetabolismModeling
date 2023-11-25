function totModel = buildFullUtilBrainModel(baseModel, fracN, fracA, TNM, TNG, TAM, TAG, mobUtilN, mobUtilA, TMT, addExtraATP)
% buildFullUtilBrainModel
%
% Generates a combined model with neurons, astrocytes, and 'the rest of the body' (rob).
% 100 model pairs of astrocytes and neurons with varying utilization (0.01 - 1) is created.
% The rest of the body (Rob) is just one model.
% We assume that each static utilization occurs equally often across time.
% The combined model thus contains 201 models. It is therefore important to have
% a minimized model.
% The model is assumed to be run with a total glucose uptake bound of 1.
%
% Input:
%
%   baseModel       Should be an ec model generated with CreateCompartmentECGEM, 
%                   and addPenaltiesToModel should not be used on it.
%
%   fracN           fraction of neurons in the model
%
%   fracA           fraction of astrocytes in the model. 
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
%   TMT             Transportation penalty for MT enzymes, - optional, default 0
%
%   addExtraATP     Adds extra ATP to ensure that we don't run out of ATP in slices - optional, default false
%
% Output:
%
%   outModel        The full model
%

if nargin < 10
    TMT = 0;
end
if nargin < 11
    addExtraATP = false;
end

%for debugging
%baseModel = minModel;
%fracN = 0.17;
%fracA = 0.03;
%TA = 0.1;
%TN = 0.1;
%mobUtilN = 0.4;
%mobUtilA = 0.2;

%Model tagging
%The Rob model is not tagged and has no penalites from utilization
%The N and A models are tagged n_[x]_, where [x] is a number between 1-100 stating the static utilization in percent

%Steps for building the model:
% 1. Create penalized models for Rob and all slices, and merge them. Lactate can be transported
%    freely within each slice, i.e. between astrocytes and neurons. Glucose uptake is blocked in the slices.
% 2. Add free communication of glucose from [s] to all slices.
% 3. Add an uptake of lactate from rob in all slices - but only 25% of the glucose uptake (times 2 since it is lactate). 
%    Also transport in the opposite direction (which is even more important).
% 4. Add communication of lactate from slices with higher utilization - these will always be there 
%    when a slice is active. For simplicity, we don't constrain the flux - it will not matter practically.
% 5. Create a common ATP hydrolysis function that consumes ATP from all slices in the right amounts
% 6. Constrain the lactate uptake in the Rob to ~25% of the expected glucose uptake (then assuming no lactate uptake, also times 2 since it is lactate)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Create the combined model
%%%%%%%%%%%%%%%%%%%%%%%%%%
rob = addPenaltiesToModel(baseModel, 1, 1, 1, 1, 1, true); %untagged, no penalties

totModel = rob; %add Rob

%remove in exch rxns for glucose and lactate for the neurons and astrocytes - this is to be distributed from rob
ANBase = baseModel;
sel = ismember(ANBase.rxns, {'MAR09034_REV';'MAR09135_REV'});
%constructEquations(ANBase, ANBase.rxns(sel)) %test, looks ok
ANBase.ub(sel) = 0;

%Now create the slices, both N and A
for i = 1:100
    NTag = ['N_' num2str(i) '_'];
    ATag = ['A_' num2str(i) '_'];
    N = addPenaltiesToModel(ANBase, 1, TNM, TNG, i/100, i/100 + (1 - i/100)*mobUtilN, false, TMT);
    A = addPenaltiesToModel(ANBase, 1, TAM, TAG, i/100, i/100 + (1 - i/100)*mobUtilA, false, TMT);
    N = tagModel(N, NTag);
    A = tagModel(A, ATag);
    AN = joinModels(N,A);
    %Add free communication of lactate between A and S, since they are close spatially
    rxnsToAdd = struct();   
    rxnsToAdd.rxns = {['internal_AN_lactate_' num2str(i)]};
    rxnsToAdd.equations = {['L-lactate[' NTag 's] <=> L-lactate[' ATag 's]']};
    AN = addRxns(AN,rxnsToAdd, 3);

    totModel = joinModels(totModel, AN);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Add glucose uptake to slices
%%%%%%%%%%%%%%%%%%%%%%%%%%

rxnsToAdd = struct();   
rxnsToAdd.rxns = cell(200,1);
rxnsToAdd.equations = cell(200,1);

for i = 1:100
    NTag = ['N_' num2str(i) '_'];
    ATag = ['A_' num2str(i) '_'];
    
    rxnsToAdd.rxns{i*2-1} = ['gluc_transp_' NTag];
    rxnsToAdd.rxns{i*2} = ['gluc_transp_' ATag];
    rxnsToAdd.equations{i*2-1} = ['glucose[s] <=> glucose[' NTag 's]'];
    rxnsToAdd.equations{i*2} = ['glucose[s] <=> glucose[' ATag 's]'];
end

totModel = addRxns(totModel,rxnsToAdd, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. Add an uptake of lactate from rob in all slices - but only say 1/5 of the glucose uptake. Also transport in the opposite direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%

%We only need to supply it to the neurons - there is already communication to the astrocytes from the neurons
rxnsToAdd = struct();   
rxnsToAdd.rxns = cell(100,1);
rxnsToAdd.equations = cell(100,1);
rxnsToAdd.ub = nan(100,1);

for i = 1:100
    NTag = ['N_' num2str(i) '_'];
    
    rxnsToAdd.rxns{i} = ['lac_transp_' NTag];
    rxnsToAdd.equations{i} = ['L-lactate[s] <=> L-lactate[' NTag 's]'];
    
    %We calculate the expected uptake of glucose if it is fully oxidized everywhere:
    fracB = fracA + fracN;
    totGlucBrain = fracB * 1;
    totGlucHighestSlice = totGlucBrain /(100*50.5); %The average slice takes up 50.5 times more than the 0.01 util slice (mean of 1:100), and we have 100 slices
    
    rxnsToAdd.ub(i) = totGlucHighestSlice * i / 2; %expected glucose uptake in a slice is i * totGlucHighestSlice, the 2 is to say that there is less lactate transported in the blood, max ~25 perc of the energy.
end

totModel = addRxns(totModel,rxnsToAdd, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. Add communication of lactate from slices with higher utilization
%%%%%%%%%%%%%%%%%%%%%%%%%%

%We only need to supply it to the neurons - there is already communication to the astrocytes from the neurons
numRxns = 100*(100-1)/2;
rxnsToAdd = struct();   
rxnsToAdd.rxns = cell(numRxns,1);
rxnsToAdd.equations = cell(numRxns,1);
nNextInd = 1;
for i = 1:100
    for j = (i+1):100
        NTag1 = ['N_' num2str(i) '_'];
        NTag2 = ['N_' num2str(j) '_'];
        RTag = ['R_' num2str(i) '_' num2str(j) '_'];

        rxnsToAdd.rxns{nNextInd} = ['lac_slice_transp_' RTag];
        rxnsToAdd.equations{nNextInd} = ['L-lactate[' NTag2 's] => L-lactate[' NTag1 's]'];
        nNextInd = nNextInd + 1;
    end
end

%test:
%numRxns == nNextInd - 1 %ok

totModel = addRxns(totModel,rxnsToAdd, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%5. Create a common ATP hydrolysis function that consumes ATP from all slices in the right amounts
%%%%%%%%%%%%%%%%%%%%%%%%%%

fracRob = 1 - fracN - fracA;

%First add rob consumption
rxnLeft = [num2str(fracRob) ' ATP[c] + ' num2str(fracRob), ' H2O[c]'];
rxnRight = [num2str(fracRob) ' ADP[c] + ' num2str(fracRob) ' Pi[c] + ' num2str(fracRob) ' H+[c] + '];

%Now all slices
for i = 1:100
    NTag = ['N_' num2str(i) '_'];
    ATag = ['A_' num2str(i) '_'];

    scale = i/(100*50.5);%The average slice takes up 50.5 times more than the 0.01 util slice (mean of 1:100), and we have 100 slices
    leftA = [num2str(scale*fracA) ' ATP[' ATag 'c] + ' num2str(scale*fracA), ' H2O[' ATag 'c]'];
    leftN = [num2str(scale*fracN) ' ATP[' NTag 'c] + ' num2str(scale*fracN), ' H2O[' NTag 'c]'];
    
    rightA = [num2str(scale*fracA) ' ADP[' ATag 'c] + ' num2str(scale*fracA) ' Pi[' ATag 'c] + ' num2str(scale*fracA) ' H+[' ATag 'c]'];
    rightN = [num2str(scale*fracN) ' ADP[' NTag 'c] + ' num2str(scale*fracN) ' Pi[' NTag 'c] + ' num2str(scale*fracN) ' H+[' NTag 'c]'];

    rxnLeft = [rxnLeft ' + ' leftA ' + ' leftN];
    rxnRight = [rxnRight ' + ' rightA ' + ' rightN];
end

rxnsToAdd = struct();   
rxnsToAdd.rxns = {'tot_atp_hydr'};
rxnsToAdd.equations = {[rxnLeft ' => ' rxnRight]};

totModel = addRxns(totModel,rxnsToAdd, 3);

%test
%sel = strcmp(totModel.rxns,'tot_atp_hydr');
%subS = totModel.S(:, sel);
%sum(subS(subS > 0)) %3, as it should be (ADP + Pi + H+), ok
%sum(subS(subS < 0)) % -2, as it should be (H2O + ATP), ok

%%%%%%%%%%%%%%%%%%%%%%%%%%
%6. Constrain the lactate uptake in the Rob to ~25% of the expected glucose uptake (then assuming no lactate uptake)
%%%%%%%%%%%%%%%%%%%%%%%%%%
totModel.ub(strcmp(totModel.rxns, 'MAR08515')) = 2*fracRob*0.25;%total glucos supply is 1, 2 is for that 2 lactate corresponds to 1 glucose


%%%%%%%%%%%%%%%%%%%%%%%%%%
%7. For some simulations, the ATP maintenance costs become higher than the ATP production,
%   which messes up the simulation. To solve this, we move the maintenance costs to the ROB
%   model, and add extra ATP there.
%%%%%%%%%%%%%%%%%%%%%%%%%%
if addExtraATP
    for i = 1:100
        NTag = ['N_' num2str(i) '_'];
        ATag = ['A_' num2str(i) '_'];
        rxnNamesN = {[NTag 'prot_maint_other'];[NTag 'prot_maint_cyto'];[NTag 'prot_maint_mito'];[NTag 'prot_maint_mt']};
        rxnNamesA = {[ATag 'prot_maint_other'];[ATag 'prot_maint_cyto'];[ATag 'prot_maint_mito'];[ATag 'prot_maint_mt']};
        rxnIndN = find(ismember(totModel.rxns,rxnNamesN));
        rxnIndA = find(ismember(totModel.rxns,rxnNamesA));
        
        %Get indices of metabolites in the S matrix
        %H2O[c] + ATP[c] => ADP[c] + Pi[c] + H+[c] + X other_prot_pool[c]
        sliceNCompInd = find(strcmp(totModel.comps, [NTag 'c']));
        sliceACompInd = find(strcmp(totModel.comps, [ATag 'c']));
        robCompInd = find(strcmp(totModel.comps, 'c'));
        metH2OSliceNInd = find(totModel.metComps == sliceNCompInd & strcmp(totModel.metNames, 'H2O'));
        metH2OSliceAInd = find(totModel.metComps == sliceACompInd & strcmp(totModel.metNames, 'H2O'));
        metH2ORobInd = find(totModel.metComps == robCompInd & strcmp(totModel.metNames, 'H2O'));

        metATPSliceNInd = find(totModel.metComps == sliceNCompInd & strcmp(totModel.metNames, 'ATP'));
        metATPSliceAInd = find(totModel.metComps == sliceACompInd & strcmp(totModel.metNames, 'ATP'));
        metATPRobInd = find(totModel.metComps == robCompInd & strcmp(totModel.metNames, 'ATP'));

        metADPSliceNInd = find(totModel.metComps == sliceNCompInd & strcmp(totModel.metNames, 'ADP'));
        metADPSliceAInd = find(totModel.metComps == sliceACompInd & strcmp(totModel.metNames, 'ADP'));
        metADPRobInd = find(totModel.metComps == robCompInd & strcmp(totModel.metNames, 'ADP'));

        metPiSliceNInd = find(totModel.metComps == sliceNCompInd & strcmp(totModel.metNames, 'Pi'));
        metPiSliceAInd = find(totModel.metComps == sliceACompInd & strcmp(totModel.metNames, 'Pi'));
        metPiRobInd = find(totModel.metComps == robCompInd & strcmp(totModel.metNames, 'Pi'));

        metHSliceNInd = find(totModel.metComps == sliceNCompInd & strcmp(totModel.metNames, 'H+'));
        metHSliceAInd = find(totModel.metComps == sliceACompInd & strcmp(totModel.metNames, 'H+'));
        metHRobInd = find(totModel.metComps == robCompInd & strcmp(totModel.metNames, 'H+'));
        
        %now change the prot pool reactions so we consume ATP in Rob instead of in slices
        robMetInd = [metH2ORobInd;metATPRobInd;metADPRobInd;metPiRobInd;metHRobInd];
        sliceMetInd = [metH2OSliceNInd;metATPSliceNInd;metADPSliceNInd;metPiSliceNInd;metHSliceNInd];
        totModel.S(robMetInd,rxnIndN) = totModel.S(sliceMetInd,rxnIndN);
        totModel.S(sliceMetInd,rxnIndN) = 0;
        sliceMetInd = [metH2OSliceAInd;metATPSliceAInd;metADPSliceAInd;metPiSliceAInd;metHSliceAInd];
        totModel.S(robMetInd,rxnIndA) = totModel.S(sliceMetInd,rxnIndA);
        totModel.S(sliceMetInd,rxnIndA) = 0;
    end
    %And, Add some extra ATP to the ROB model
    %Perhaps not needed?
end


end
