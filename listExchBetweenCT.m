function res = listExchBetweenCT(model, sol, lowerLimit)
if nargin < 3
    lowerLimit = 0;
end

%just list all communication with the [s] compartment
sComp = find(strcmp(model.comps,'s'));%1
sMets = model.metComps == sComp;

rxnsWithsMets = sum(model.S(sMets,:) ~= 0, 1) > 0;
sum(rxnsWithsMets)%11529


fluxSel = sol.x > lowerLimit;
totSel = rxnsWithsMets.' & fluxSel;

fracFluxOfUb = sol.x./model.ub;

tbl = table(model.rxns(totSel), sol.x(totSel), fracFluxOfUb(totSel), constructEquations(model, model.rxns(totSel)));
tbl.Properties.VariableNames = {'Rxn', 'Flux', 'FracFluxOfUb', 'Formula'};
tbl
res = tbl;

end
