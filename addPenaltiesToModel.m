function model = addPenaltiesToModel(model, B, TM, TG, UCyto, UMito, dontPenalize)
% addPenaltiesToModel
%
% Adds EAMCA maintenance costs to the model
%
% Input:
%
%   model           The model to modify
%
%   B               Base cost - maintenance cost per protein unit without transportation or utilization penalties
%
%   TM              Transportation penalty, mitochondria - typically a number between 0 and 2, such as 0.10
%
%   TG              Transportation penalty, glycolysis - typically a number between 0 and 2, such as 0.10
%
%   UCyto           Utilization of cytosolic proteins, such as glycolysis
%
%   UMito           Utilization of mitochondrial proteins, such as mitochondrial respiration
%
%   dontPenalize    If true, we add reactions to consume protein pool, but at no cost. Default: false
%
% Output:
%
%   outModel        The penalized model
%


if nargin < 7
    dontPenalize = false;
end
costOther = B*(1+TG)/UCyto;
costCyto = B*(1+TG)/UCyto;
costMito = B*(1+TM)/UMito;
costMT = B/UMito;

rxnsToAdd = {};
rxnsToAdd.rxns = {'prot_maint_other';'prot_maint_cyto';'prot_maint_mito';'prot_maint_mt'};
rxnsToAdd.rxnNames = rxnsToAdd.rxns;
rxnsToAdd.lb = [0;0;0;0];
rxnsToAdd.ub = [1000;1000;1000;1000];
%rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnData{5});
if dontPenalize
    rxnsToAdd.equations = {'=> other_prot_pool[c]';
                           '=> cyto_prot_pool[c]';
                           '=> mito_prot_pool[c]';
                           '=> mt_prot_pool[c]'};
else
    rxnsToAdd.equations = {['H2O[c] + ATP[c] => ADP[c] + Pi[c] + H+[c] + ' num2str(1/costOther) ' other_prot_pool[c]'];
                           ['H2O[c] + ATP[c] => ADP[c] + Pi[c] + H+[c] + ' num2str(1/costCyto) ' cyto_prot_pool[c]'];
                           ['H2O[c] + ATP[c] => ADP[c] + Pi[c] + H+[c] + ' num2str(1/costMito) ' mito_prot_pool[c]'];
                           ['H2O[c] + ATP[c] => ADP[c] + Pi[c] + H+[c] + ' num2str(1/costMT) ' mt_prot_pool[c]']};
end
%rxnsToAdd.equations

model = addRxns(model,rxnsToAdd,3);
