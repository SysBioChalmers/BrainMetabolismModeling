function [totModel] = buildSimpleAstroNeuronModel(baseModel, fracA, TNM, TNG, TAM, TAG, statUtil, mobUtilN, mobUtilA)
% buildSimpleAstroNeuronModel
%
% Generates a combined model having one astrocyte model and one neuron model.
% The models can exchange lactate. Glucose is administered via the neuron model (untagged)
%
% Input:
%
%   baseModel       Should be an ec model generated with CreateCompartmentECGEM, 
%                   and addPenaltiesToModel should not be used on it.
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
%   statUtil        Static utilization
%
%   mobUtilN        Fraction of unused time in neurons that mitochondria can be used somewhere else
%
%   mobUtilA        Fraction of unused time in astrocytes that mitochondria can be used somewhere else
%
% Output:
%
%   outModel        The full model
%

%for debugging
%baseModel = minModel;
%fracA = 0.15;
%TNM = 0.1;
%TNG = 0.1;
%TAM = 0.1;
%TAG = 0.1;
%statUtil = 0.5;
%mobUtilN = 0.4;
%mobUtilA = 0.2;

%Model tagging
%The A model is tagged with 'A_'
ATag = 'A_';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Create the combined model
%%%%%%%%%%%%%%%%%%%%%%%%%%
N = addPenaltiesToModel(baseModel, 1, TNM, TNG, statUtil, statUtil + (1 - statUtil)*mobUtilN);
A = addPenaltiesToModel(baseModel, 1, TAM, TAG, statUtil, statUtil + (1 - statUtil)*mobUtilA);
A = tagModel(A, ATag);

%remove in exch rxns for glucose and lactate for the astrocytes - this is to be distributed from the neurons
sel = ismember(baseModel.rxns, {'MAR09034_REV';'MAR09135_REV'});
%constructEquations(A, A.rxns(sel)) %test, looks ok
A.ub(sel) = 0;

totModel = joinModels(N,A);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Add communication of lactate and glucose
%%%%%%%%%%%%%%%%%%%%%%%%%%

rxnsToAdd = struct();   
rxnsToAdd.rxns = {'lactate_transport'; 'glucose_transport'};
rxnsToAdd.equations = {['L-lactate[s] <=> L-lactate[' ATag 's]']; ['glucose[s] <=> glucose[' ATag 's]']};
totModel = addRxns(totModel,rxnsToAdd, 3);
%constructEquations(totModel, rxnsToAdd.rxns) %test, looks ok

%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. Create a common ATP hydrolysis function that consumes ATP from both models
%%%%%%%%%%%%%%%%%%%%%%%%%%

fracN = 1 - fracA;

rxnLeft = [num2str(fracN) ' ATP[c] + ' num2str(fracN) ' H2O[c] + ' num2str(fracA) ' ATP[' ATag 'c] + '  num2str(fracA) ' H2O[' ATag 'c]'];
rxnRight = [num2str(fracN) ' ADP[c] + ' num2str(fracN) ' Pi[c] + ' num2str(fracN) ' H+[c] + ' num2str(fracA) ' ADP[' ATag 'c] + ' num2str(fracA) ' Pi[' ATag 'c] + ' num2str(fracA) ' H+[' ATag 'c]'];

rxnsToAdd = struct();   
rxnsToAdd.rxns = {'tot_atp_hydr'};
rxnsToAdd.equations = {[rxnLeft ' => ' rxnRight]};
totModel = addRxns(totModel,rxnsToAdd, 3);
%constructEquations(totModel, rxnsToAdd.rxns) %test, looks ok

%test
%sel = strcmp(totModel.rxns,'tot_atp_hydr');
%subS = totModel.S(:, sel);
%sum(subS(subS > 0)) %3, as it should be (ADP + Pi + H+), ok
%sum(subS(subS < 0)) % -2, as it should be (H2O + ATP), ok

end
