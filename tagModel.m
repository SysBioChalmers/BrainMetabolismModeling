function [outModel] = tagModel(inModel, tag)

outModel = inModel;
outModel.rxns = strcat(tag,inModel.rxns);
outModel.mets = strcat(tag,inModel.mets);
outModel.comps = strcat(tag,inModel.comps);
outModel.compNames = strcat(tag,inModel.compNames);
outModel.subSystems = strcat(tag,inModel.mets);
