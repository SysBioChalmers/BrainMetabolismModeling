cd 'C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling'
load('data/MinModel.mat')

%Data for Fig. 4E, S2
simulateCombSimpleModel(minModel, 'data/simpleModelOutputRedMitMob.mat', 0.15, 0.1, 0.1, 0.1, 0.1, 0.5, 0.4, 0.2);
simulateCombSimpleModel(minModel, 'data/simpleModelOutputTransport.mat', 0.15, 0.1, 0.2, 0.2, 0.2, 0.5, 0.4, 0.4);

%Data for Fig. 4G:
simulateFullModel(minModel, 'data/FullModelOutputRedMitMob.mat', 0.17, 0.03, 0.1, 0.1, 0.1, 0.1, 0.4, 0.2);

%Data for Fig. S3:
simulateFullModel(minModel, 'data/FullModelOutputTransport.mat', 0.17, 0.03, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4);

%Data for Fig. S4:
simulateFullModel(minModel, 'data/FullModelOutputRedMitMobTMT.mat', 0.17, 0.03, 0.1, 0.1, 0.1, 0.1, 0.4, 0.2, 0.1);

%Sensitivity analysis for transport cost
%Simulate between 0.01 and 100
log10Vals = -2:0.1:2;
TVals = 10 .^ log10Vals;
TVals
save('data/Sens/TVals.mat', 'TVals')
for i = 1:length(TVals)
    disp(num2str(i))
    simulateFullModel(minModel, ['data/Sens/' num2str(i) '.mat'], 0.17, 0.03, TVals(i), TVals(i), TVals(i), TVals(i), 0.4, 0.2, 0, true);
end

x = load(['data/Sens/' num2str(10) '.mat']);
x.d
x.d.ATPConsNeur(100)
x.d.ATPProdGlycNeur
x.d.ATPProdMitoNeur
[x.d.ATPProdGlycNeur x.d.ATPProdMitoNeur]
%Assemble the data for astrocytes and neurons
%What we do here is basically to just determine if we have glyc, mito, or both for each slice
grid = zeros(100,length(TVals));
for i = 1:length(TVals)
    disp(num2str(i))
    fn = ['data/Sens/' num2str(i) '.mat'];
    x = load(fn);
    
    simulateFullModel(minModel, ['data/Sens/' num2str(i) '.mat'], 0.17, 0.03, TVals(i), TVals(i), TVals(i), TVals(i), 0.4, 0.2, 0, true);
end



