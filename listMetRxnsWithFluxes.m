function res = listMetRxnsWithFluxes(model, sol, met, consumption, lowerLimit)
if nargin < 5
    lowerLimit = 0;
end
mets = strcmp(model.metNames,met);
SMultFlux = model.S(mets,:) .* sol.x.';
if (consumption)
    rxns = sum(SMultFlux < -lowerLimit, 1) > 0;
    SMultFluxNeg = SMultFlux(:,rxns);
    SMultFluxNeg(SMultFluxNeg > 0) = 0;
    fluxes = sum(-SMultFluxNeg,1);
else
    rxns = sum(SMultFlux > lowerLimit, 1) > 0;
    SMultFluxPos = SMultFlux(:,rxns);
    SMultFluxPos(SMultFluxPos < 0) = 0;
    fluxes = sum(SMultFluxPos,1);
end
%select those that carry flux at a[15]
%rxnsWithFlux = rxns.' & (sol.x ~= 0);
tbl = table(model.rxns(rxns), fluxes.', constructEquations(model, model.rxns(rxns)));
tbl.Properties.VariableNames = {'Rxn', 'Metabolite Flux', 'Formula'};
tbl
res = tbl;

end
