cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
load('data/MinModel.mat')

%Data for Fig. S2
simulateCombSimpleModel(minModel, 'data/simpleModelOutputRedMitMob.mat', 0.15, 0.1, 0.1, 0.1, 0.1, 0.5, 0.4, 0.2);
simulateCombSimpleModel(minModel, 'data/simpleModelOutputTransport.mat', 0.15, 0.1, 0.2, 0.2, 0.2, 0.5, 0.4, 0.4);

%Data for Fig. 2F:
simulateFullModel(minModel, 'data/FullModelOutputRedMitMob.mat', 0.17, 0.03, 0.1, 0.1, 0.1, 0.1, 0.4, 0.2);

%Data for Fig. S3:
simulateFullModel(minModel, 'data/FullModelOutputTransport.mat', 0.17, 0.03, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4);

