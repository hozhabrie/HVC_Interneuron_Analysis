clear all; close all; clc; warning off;

%% Define Paths
path = '/Users/ellie/NYU Langone Health Dropbox/Elnaz Hozhabri/Int Juxta/';
datafile = fullfile(path, 'goodCellsnew.xlsx');

%% Load Data
opts = detectImportOptions(datafile);
opts = setvartype(opts, 'stdmax', 'string'); 
goodCells = readtable(datafile, opts);

% Extract relevant variables
birds = string(goodCells{:,1});
cells = string(goodCells{:,2});
toPrint = goodCells{:,5};
identity = string(goodCells{:,9});
badx = identity == 'x';
fontsize = 15;
allC = [];

% Filter identities of interest
refound = identity(ismember(identity, {'pv', 'som', 'reln', 'pvsom', 'nothing'}));

% Fix the positive deflecting juxta spikes
deflect = contains(string(goodCells{toPrint == 1,4}), '+');

%% Load and Process Cluster Data
for i = 1:height(goodCells)
    if toPrint(i) == 1 && ~badx(i)
        changedir = fullfile(path, birds(i), cells(i));
        clustersFile = fullfile(changedir, 'clusters.mat');
        
        if exist(clustersFile, 'file')
            load(clustersFile);
            if isfield(clusters, 'spontStateBrian')
                clusters = rmfield(clusters, 'spontStateBrian');
                disp(i);
            end
            allC = [allC clusters];
        end
    end
end

%% Extract Raksin Variables
raksinVars = {'pbAwakeFRAll', 'pbSleepFRAll', 'prePbAwakeFRAll', 'prePbSleepFRAll'};
for n = 1:length(allC)
    for m = 1:length(raksinVars)
        raksin(n).(raksinVars{m}) = allC(n).(raksinVars{m});
    end
end

%% Prepare Juxta Structure
allC = orderfields(allC);
strOG = fieldnames(allC);
str = strOG;

allCJuxta.clusterID = string([allC.bird] + [allC.cell]);
allCJuxta.identity = identity(toPrint == 1 & ~badx)';

% Manually fix values that are nan because empty but existent trials
allCJuxta.pbAwakeCorr([2 24 27 31]) = 0;
allCJuxta.pbSleepCorr([2 20 32 50 51]) = 0;
allCJuxta.pbAwakeChange(2) = 0;
allCJuxta.pbSleepChange(32) = 0; %2


%% Calculate pbComboChange
comboChange = arrayfun(@(x) nanmean(([x.pbAwakeFRAll; x.pbSleepFRAll] - [x.prePbAwakeFRAll; x.prePbSleepFRAll]) ./ ([x.pbAwakeFRAll; x.pbSleepFRAll] + [x.prePbAwakeFRAll; x.prePbSleepFRAll])), raksin);
allCJuxta.pbComboChange = comboChange;

%% Identify Found Subtypes
load(fullfile(path, 'allCJuxta.mat'));
pv = allCJuxta.identity == 'pv';
som = allCJuxta.identity == 'som';
reln = allCJuxta.identity == 'reln';
pvsom = allCJuxta.identity == 'pvsom';
none = allCJuxta.identity == 'nothing';
unfound = ~ismember(allCJuxta.identity, {'pv', 'som', 'reln', 'pvsom', 'nothing'});

%% 3D Plot of Variables
wf = allCJuxta.wfRise;
modul = abs(allCJuxta.pbSleepChange);
state = allCJuxta.spontStateCorr;
id = allCJuxta.clusterID;

fig1 = figure;
hold on;
colors = {[.6 .6 .6], 'b', 'r', 'g', 'k', 'c'};
groups = {unfound, pv, som, reln, none, pvsom};
labels = {'Unfound', 'PV', 'SOM', 'RELN', 'None', 'PVSOM'};

for i = 1:numel(groups)
    plot3(wf(groups{i}), modul(groups{i}), state(groups{i}), 'o', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', [.6 .6 .6]);
end
text(wf(~unfound), modul(~unfound), state(~unfound), id(~unfound));

xlabel('Rise'); ylabel('PB Sleep Change'); zlabel('ISI Corr');
set(gca, 'TickDir', 'out', 'Color', 'none', 'Box', 'off');
fig1.Renderer = 'Painters';
print(fullfile(path, '3vars_everything'), '-dpdf');

%% Save
cd(path);
save(fullfile(path, 'allCJuxtaFull.mat'), 'allCJuxta');
disp('Saved allCJuxta');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extra Plots and Quantifications (not always run)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster 3 Variables
clc;
addpath('/Users/ellie/NYU Langone Health Dropbox/Elnaz Hozhabri/EllieScripts/');

use = [wf;modul;state]';
zuse = [zscorenan(wf);zscorenan(modul);matchRatio(state);zscorenan(allCJuxta.spontAwakeCV);zscorenan(sqrt(allCJuxta.spontAwakeBurst))]';

symbols = 'svd*o.|*.|sv';
colors = 'gbrmkcgbmkrc';
meas = 12;

clear clus;
for iter = 1:meas
    clus(:,iter) = kmeans(zuse, iter, 'Replicates', 1000);
end

clusData = evalclusters(zuse, clus, 'CalinskiHarabasz');
k = clusData.OptimalK;
whichClus = clus(:,k);

fig2 = figure; plot(clusData);
set(gca, 'TickDir','out','color','none','FontSize',15,'TickDir','out','box','off');

fig3 = figure; 
hold on;
for i = 1:k
    pts = use(whichClus==i,:);
    plot3(pts(:,1),pts(:,2),pts(:,3),[colors(i) symbols(i)]);
end
set(gca, 'TickDir','out','color','none','FontSize',15,'TickDir','out','box','off');
fig3.Renderer='Painters';
print(fullfile(path, '3d_clusters'), '-dpdf');

%% Compute Subtype Distances
close all; clear randDist; clear randMat; clear pvDist; clear pvMat;

pvMat = use(pvsom,:);
for i = 1:length(pvMat)
    for j = 1:length(pvMat)
        pvDist(i,j) = norm(pvMat(i,:) - pvMat(j,:));
    end
end
whichVals = triu(logical(ones(length(pvMat),length(pvMat))),1);
pvVals = pvDist(whichVals);

numRand = 1e4;
randVals = zeros([numRand,1]);

for a = 1:10000
    randMat = use(randperm(length(allC),length(pvMat)),:);
    for i = 1:length(randMat)
        for j = 1:length(randMat)
            randDist(i,j) = norm(randMat(i,:) - randMat(j,:));
        end
    end
    randVals(a) = mean(randDist(whichVals));
end

pvPrctile = sum(randVals < mean(pvVals)) / numRand;

fig1 = figure;
fig1.Renderer='Painters';
histogram(randVals)
hold on;
line([mean(pvVals),mean(pvVals)],[0 max(ylim)],'Color','r')
title(['pvsstPrctile: ' num2str((1 - pvPrctile) * 100)])
set(gca, 'TickDir','out','color','none','FontSize',15,'TickDir','out','box','off')
xlabel('mean distance')
ylabel('Count')
cd('C:\Users\User\Dropbox (NYU Langone Health)\Paper\fig4\')
print ('pvsom dist','-dpdf')

%% Raksin comparison for thesis

%wake
change = allCJuxta.pbAwakeChange;
for i = 1:length(raksin)
    pre = raksin(i).pbAwakeFRAll;
    post = raksin(i).prePbAwakeFRAll;
    if ~isempty(pre) & ~isempty(post)
        [pval(i) sig(i)] = signrank(pre,post);
    end
end

unmod = ~sig;
pos = (change > 0) & sig;
neg = (change < 0) & sig;
mod = sig;

clear whichpb
whichpb(pos) = {'up'};
whichpb(neg) = {'down'};
whichpb(unmod) = {'same'};

% sleep
% change = allCJuxta.pbSleepChange;
% for i = 1:length(raksin)
%     pre = raksin(i).prePbSleepFRAll;
%     post = raksin(i).pbSleepFRAll;
%     if ~isempty(pre) && ~isempty(post)
%         [pval(i) sig(i)] = signrank(pre,post);
%     else
%        pval(i) = nan;
%     end
% end
% unmod = ~sig;
% pos = (change > 0) & sig;
% neg = (change < 0) & sig;
% mod = sig;
% 
% clear whichpb
% whichpb(pos) = {'up'};
% whichpb(neg) = {'down'};
% whichpb(unmod) = {'same'};

%% plot params

cd('C:/Users/User/Dropbox (NYU Langone Health)/Thesis/Dissertation/figures/')

samesies = {'wfRise','wfFWHM','pbAwakeFR','pbAwakeChange','pbAwakeCorr','spontAwakeFR','spontAwakeCV',...
    'spontStateCorr','spontRatioVigilance','pbSleepFR','pbSleepChange','pbSleepCorr','spontSleepFR','spontSleepCV'};

thesisWords = {'WF Rise','WF FWHM', 'Playback FR_O_N','Playback FR Change_O_N','Playback Precision_O_N','Spont FR_O_N','Spont CV_O_N',...
    'ISI Corr','Vigilance Change','Playback FR_O_F_F','Playback FR Change_O_F_F','Playback Precision_O_F_F', 'Spont FR_O_F_F','Spont CV_O_F_F'};

use = nan(length(allC),length(samesies));

for i = 1:length(samesies)
    field = string(samesies(i));
    figure('Visible','off');
    plotVarJuxta(allCJuxta,field)
    print(['juxta_'+ field], '-dpdf','-fillpage');
    set(gcf,'PaperOrientation','landscape')
    pars = getfield(allCJuxta,field);  
    if i == 4 | i == 11
       pars = abs(pars);
    elseif i == 3 | i == 6 | i == 10 | i == 13
       pars = log(pars);
    elseif i == 5 | i == 8 | i == 12
        pars = matchRatio(pars);
    end
    pars(isinf(pars)) = nan;
    use(:,i) = pars;
end

%normalize variables
useGood = zscorenan(use);

histOrder = [1,6,7,13,14,8,9,2,3,4,5,10,11,12];

figure;     
for i = 1:length(samesies)
    v = useGood(:,histOrder(i));
    subplot(7,2,i)
    histogram(v,30,'DisplayStyle','stairs','EdgeColor','k')
    xlim([-5 5])
    xticks(-4:2:4)
    set(gca,'box','off','TickDir','out','color','none')
    [h p] = swtest(v);
    pass = p < .05; 
    if pass == 1
        verdict = 'not normal';
    else 
        verdict = 'normal';
    end
    title(thesisWords(histOrder(i)))
end

%% Make corr matrix of params

useGoodFull = useGood([1,3:26,28:30,33:35,37:40,42:49,52:end],:);
corrOrder = [1 2 6 7 13 14 8 9 3 4 5 10 11 12];
useCorr = useGoodFull(:,corrOrder);
[R,P] = corrcoef(useCorr);

figure; imagesc(R); hold on;
sig = string(nan(size(P)));
sig(P<.05) = '*';
nvars = length(samesies);
for a = 1:nvars
    i = a-1;
    text([1:nvars],a*ones(nvars,1),sig((1:nvars)+nvars*i),'Color','k')
    hold on;
    shg;
end

xticks(1:length(samesies))
xticklabels(thesisWords(corrOrder))
yticks(1:length(samesies))
yticklabels(thesisWords(corrOrder))
set(gca,'box','off','TickDir','out')
colorbar 

useGood = useGood(:,corrOrder);

%% PCA/dim red plots
 
[coeff,score,latent,~,explained] = pca(useGood);
% % 
x = score(:,1);
y = score(:,2);
z = score(:,3);
d = .04;

figure;
scatter3(x,y,z,'MarkerEdgeColor',[.5 .5 .5]);
hold on;
scatter3(x(pv),y(pv),z(pv),'filled','b');
hold on; 
scatter3(x(som),y(som),z(som),'filled','r');
hold on
scatter3(x(reln),y(reln),z(reln),'filled','g');
hold on
scatter3(x(none),y(none),z(none),'k');
scatter3(x(pvsom),y(pvsom),z(pvsom),'filled','c');
xlabel('pc1')
ylabel('pc2')
zlabel('pc3')
grid off;
xticks(-4:2:4)
yticks(-4:2:4)
zticks(-4:2:4)
set(gca,'box','off','TickDir','out','color','none')

cumulsum = cumsum(explained);
figure; plot(cumulsum,'k-o')
set(gca,'box','off','TickDir','out','color','none')

elbow = find(cumulsum > 90, 1);

%sum(explained(1:elbow))
figure; plot(mean(abs(coeff(:,1:elbow)),2),'k-o')
set(gca,'box','off','TickDir','out','color','none')

xticks(1:length(samesies(1:end)))
xticklabels(thesisWords(corrOrder(1:end)))
shg

%% Hierarchical clustering plots

%account for pos deflecting spikes!
transform = struct2cell(allCJuxta);
use = cell2mat(transform(2:end))';
clear transform;
use(isinf(use)) = nan;
tree = linkage(use(~deflect,:),'ward',@naneucdist);

figure;
subplot(1,2,1)
[~,~,outperm] = dendrogram(tree,0,'Orientation','left','Labels',identity(~deflect),'ColorThreshold','default');
set(gca, 'TickDir', 'out','color','none','box','off','FontSize',9)
title('all params')

params2 = erase(params,{'spontRatioVigilance';'spontAwakeBurst';'spontAwakeCV';...
'pbAwakeBurst';'pbAwakePrecInd10';'pbAwakeMaxIFR';'pbAwakeChange';'pbAwakePostChange';...
'pbSleepChange';'pbSleepPostChange'});
params2(ismissing(params2)) = [];

params3 = {'wfRise'}%,'wfRise70','pbSleepFR'%'pbSleepFR','pbSleepFRPost','pbSleepPrecInd50','spontAwakeFR','spontSleepFR','pbAwakeFR'}

clear trunk;
for k=1:length(params3)
    trunk.(params3{k}) = [allCJuxta.(params3{k})];
end
use = cell2mat(struct2cell(trunk))';
use(isinf(use)) = nan;
tree = linkage(use(~deflect,:),'ward',@naneucdist);

subplot(1,2,2)
[~,~,outperm] = dendrogram(tree,0,'Orientation','left','Labels',identity(~deflect),'ColorThreshold','default');
set(gca, 'TickDir', 'out','color','none','box','off','FontSize',9)
title('rise time only')

