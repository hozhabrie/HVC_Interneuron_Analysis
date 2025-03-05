clear all; close all; clc; warning off;

% Define directories
dirJux = 'C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/';
dirSi = 'C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/';

% Load data
load([dirJux 'allCJuxta.mat']);
load([dirSi 'allCSi.mat']);

% Read goodCells data
opts = detectImportOptions([dirJux 'goodCellsnew.xlsx']);
opts = setvartype(opts, 'stdmax', 'string'); 
goodCells = readtable([dirJux 'goodCellsnew.xlsx'], opts);
toPrint = goodCells{:,5};

set(0, 'DefaultFigureVisible', 'on');

%% Fix Juxta Labels and Fields to Match Si
identityJuxta = allCJuxta.identity;
whichJuxta = (identityJuxta ~= 'x')';
identityJuxta(contains(identityJuxta, 'nothing')) = "NOLAB";
identityJuxta = upper(identityJuxta);
identityJuxta = identityJuxta' + '-' + allCJuxta.clusterID';
allCJuxta.identity = identityJuxta';
allCJuxta.wfRatio = allCJuxta.wfFWHM ./ allCJuxta.wfRise;
allCJuxta.songFR = nan(1,59);
allCJuxta.songChange = nan(1,59);
allCJuxta.songBurst = nan(1,59);
allCJuxta.songSparse = nan(1,59);
allCJuxta.songCorr = nan(1,59);

%% Fix Si Labels to Match Juxta
identitySi = allCSi.identity';
whichSi = contains(identitySi, 'i');
allCSi.spontSleepFR = allCSi.spontSleepEarlyFR;
allCSi.spontSleepCV = allCSi.spontSleepEarlyCV;
allCSi.spontSleepBurst = allCSi.spontSleepEarlyBurst;
allCSi.wfRatio = allCSi.wfFWHM ./ allCSi.wfRise;
allCSi.spontStateCorr = allCSi.spontStateCorrE;
allCSi.spontRatioVigilance = allCSi.spontRatioVigilanceE;
allCSi.pbComboChange = allCSi.pbAwakeChange;

%% Match Parameters
samesies = {'wfRatio', 'wfRise', 'wfFWHM', 'spontAwakeCV', 'spontAwakeFR', 'spontStateCorr', 'spontSleepCV', 'spontSleepFR', 'spontRatioVigilance', 'pbAwakeCorr', 'pbAwakeChange', 'pbAwakeFR', 'pbSleepCorr', 'pbSleepChange', 'pbSleepFR', 'songFR', 'songBurst', 'songChange', 'songCorr', 'songSparse', 'pbComboChange'};

goodVars = samesies;
useOG = [];
allC = struct;

for i = 1:length(goodVars)
    if isfield(allCSi, goodVars{i})
       field = string(goodVars{i});
       a = getfield(allCSi, field);
       a = a(whichSi);
       b = getfield(allCJuxta, field);
       allC.(field) = [a'; b'];
       useOG(:,i) = [a'; b'];
    end
end

% Concatenate and label which is si or juxta
allC = orderfields(allC);
which = [whichSi'; whichJuxta];
identity = [identitySi'; identityJuxta];
ints = identity(which);
method = logical([ones(sum(whichSi),1); zeros(sum(whichJuxta),1)]);

%% Plot Parameters
set(0, 'DefaultFigureVisible', 'on');
paramsFolder = 'C:/Users/User/Dropbox (NYU Langone Health)/JointParams';
cd(paramsFolder);

for i = 1:length(goodVars)
    field = string(goodVars{i});
    close;
    plotVars(allC, field, method);
    set(gca, 'TickDir', 'out', 'color', 'none', 'FontSize', 15, 'box', 'off');
    set(gcf, 'PaperOrientation', 'landscape');
    print(field + '_paperfig', '-dpdf', '-fillpage');
end

%% Combined k-means clustering and visualization
state = matchRatio(allC.spontStateCorr);
modul = abs(allC.pbComboChange);
wf = allC.wfRatio;
zuse = [zscorenan(wf) zscorenan(modul) state];

meas = 10;
clus = zeros(size(zuse,1), meas);
for iter = 1:meas
    clus(:, iter) = kmeans(zuse, iter, 'Replicates', 1000);
end
clusData = evalclusters(zuse, clus, 'CalinskiHarabasz');
k = clusData.OptimalK;
whichClus = clus(:, k);

% Plot cluster data with optimal k
figure; plot(clusData);
title(k);

% Plot clusters and differentiate two data types
figure; plot3(wf, modul, state, 'ko');
hold on;
scatter3(wf(which), modul(which), state, 'go', 'MarkerFaceColor', 'g');
scatter3(wf(~which), modul(~which), state, 'mo', 'MarkerFaceColor', 'm');
legend('Si', 'Jux');
xlabel('wfRatio');
ylabel('abs pb change: awake (si); combo (jux)');
zlabel('corrected isi corr');

%% Hierarchical Clustering of combined data
useGoodZ = zscorenan(useOG);
tree = linkage(useGoodZ(:,[1,5,13,18]),'ward',@naneucdist);
figure;
dendrogram(tree,0,'Orientation','top','ColorThreshold',8,'Labels',ints);
set(gca, 'TickDir', 'out','color','none','box','off','FontSize',9);
shg;

%% Plot combined data using dim red and top 3 PCs
[coeff,score,latent,~,explained] = pca(useGoodZ(:,1:14));
figure;
scatter3(score(:,1), score(:,2), score(:,3), 'k');
hold on;
scatter3(score(pv,1), score(pv,2), score(pv,3), 'filled', 'b');
scatter3(score(som,1), score(som,2), score(som,3), 'filled', 'r');
scatter3(score(reln,1), score(reln,2), score(reln,3), 'filled', 'g');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); shg;
