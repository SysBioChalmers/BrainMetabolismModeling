function model = addPenaltiesToModel(model, B, T, UCyto, UMito, dontPenalize)
%dontPenalize means that we add reactions to consume protein pool, but no cost
%B is the maintenance cost per protein unit with no transport or utilization penalties

if nargin < 6
    dontPenalize = false;
end
costOther = B*T/UCyto;
costCyto = B*T/UCyto;
costMito = B*T/UMito;
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
