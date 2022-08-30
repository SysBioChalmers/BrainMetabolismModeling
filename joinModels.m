function [outModel] = joinModels(inModel1, inModel2)

outModel = inModel1;
outModel.rxns = [inModel1.rxns;inModel2.rxns];
outModel.mets = [inModel1.mets;inModel2.mets];
outModel.lb = [inModel1.lb;inModel2.lb];
outModel.ub = [inModel1.ub;inModel2.ub];
outModel.c = [inModel1.c;inModel2.c];
outModel.b = [inModel1.b;inModel2.b];
outModel.comps = [inModel1.comps;inModel2.comps];
outModel.compNames = [inModel1.compNames;inModel2.compNames];
outModel.rxnNames = [inModel1.rxnNames;inModel2.rxnNames];
outModel.grRules = [inModel1.grRules;inModel2.grRules];
outModel.subSystems = [inModel1.subSystems;inModel2.subSystems];
outModel.eccodes = [inModel1.eccodes;inModel2.eccodes];
outModel.rxnNotes = [inModel1.rxnNotes;inModel2.rxnNotes];
%outModel.genes - we assume the genes are the same for both models
outModel.metNames = [inModel1.metNames;inModel2.metNames];
outModel.metComps = [inModel1.metComps;(inModel2.metComps + length(inModel1.comps))];%so, if we have 9 compartments, the second model starts at index 10
outModel.rxnGeneMat = [inModel1.rxnGeneMat;inModel2.rxnGeneMat];%so, we assume that genes are the same for both datasets
zeroMat1 = zeros(length(inModel1.mets), length(inModel2.rxns));
zeroMat2 = zeros(length(inModel2.mets), length(inModel1.rxns));
outModel.S = [inModel1.S zeroMat1;zeroMat2 inModel2.S];
outModel.rev = [inModel1.rev;inModel2.rev];

%we skip some new fields that I don't recognize...
    