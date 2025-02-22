clear all; close all; clc; warning off;

path = 'C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/';
opts = detectImportOptions([path 'goodCellsnew.xlsx']);
opts = setvartype(opts, 'stdmax', 'string'); 
goodCells = readtable([path 'goodCellsnew.xlsx'],opts);

birds = string([goodCells{:,1}]);
cells = string([goodCells{:,2}]);
toPrint = goodCells{:,5};
identity = string([goodCells{:,9}]);
refound = identity(identity=='pv' | identity=='som' | identity=='reln' | identity=='pvsom' | identity=='nothing');
deflect = string(goodCells{:,4});
deflect = deflect(toPrint == 1);
deflect = contains(deflect, '+');
badx = identity == 'x';

fontsize = 15;
allC = [];

for i = 1:height(goodCells)
    if toPrint(i) == 1 && ~badx(i)
        changedir = [path + birds(i) + '/' + cells(i) + '/'];
        load([changedir + 'clusters.mat']);
        if sum(ismember(fieldnames(clusters),'spontStateBrian')) > 0
        clusters = rmfield(clusters, 'spontStateBrian');
        disp(i)
        end
        allC = [allC clusters];
    end
end

raksinVars = {'pbAwakeFRAll','pbSleepFRAll','prePbAwakeFRAll','prePbSleepFRAll'};
for n = 1:length(allC)
    for m = 1:length(raksinVars)
        raksin(n).(raksinVars{m}) = allC(n).(raksinVars{m});
    end
end

%mult =(([allC.wfRise]-[allC.wfRise70]).*([allC.wfRise]./[allC.wfRise70]))';

allC = orderfields(allC);

% for i = 1:length(allC)
%     on = allC(i).onIsis;
%     on(isnan(on)) = [];
%     off = allC(i).offIsis;
%     off(isnan(off)) = [];
%     if ~isempty(on) & ~isempty(off)
%         distance(i) = ws_distance(log10(on),log10(off),2);
%         correlation(i) = allC(i).spontStateCorr;
%     else
%         distance(i) = nan;
%         correlation(i) = nan;
%     end
%     allC(i).spontStateBrian = distance(i);
% end

shitfields = {'spikeTimes','isis','onIsis','offIsis','pbIsis','spontIsis','onPbIsis','offPbIsis','meanWF','wfNorm','recLength'};
allC = rmfield(allC,shitfields);
allC = rmfield(allC,raksinVars);
allC = orderfields(allC);
strOG = fieldnames(allC);
str = strOG;

for k = 3:length(str)
    allCJuxta.(str{k}) = [allC.(strOG{k})];
    disp(k)
end
%allCJuxta.wfMult = mult;
allCJuxta.clusterID = string([allC.bird] + [allC.cell]);
allCJuxta.identity = identity(toPrint == 1 & ~badx)';

%allCJuxta.pbAwakeMatchScore = 
%allCJuxta.pbSleepMatchScore = 
%headers = fieldnames(allCFixed);
%doubles = cell2mat(struct2cell(allCFixed(2:end)));
% excel = array2table(doubles','VariableNames',headers);
% excel = rows2vars(excel);
% matrix = doubles(2:end,:)';
% matrix(isinf(matrix)) = nan;

% fix values that are nan because empty but existent trials
 allCJuxta.pbAwakeCorr([2 24 27 31]) = 0; 
 allCJuxta.pbSleepCorr([2 20 32 50 51]) = 0;
 allCJuxta.pbAwakeChange(2) = 0;
 allCJuxta.pbSleepChange(32) = 0; %2
 allCJuxta.wfRatio = allCJuxta.wfFWHM./allCJuxta.wfRise;

% calculate pb as combined state
for i = 1:length(allC)
    wake = raksin(i).pbAwakeFRAll;
    sleep = raksin(i).pbSleepFRAll;
    pb = [wake;sleep];
    prewake = raksin(i).prePbAwakeFRAll;
    presleep = raksin(i).prePbSleepFRAll;
    pre = [prewake;presleep];
    comboChange(i) = nanmean((pb-pre)./(pb+pre)); 
    other(i) = ( length(wake) * nanmean((wake-prewake)./(wake+prewake)) ...
        + length(sleep) * nanmean((sleep-presleep)./(sleep+presleep)) ) ...
        ./ (length(wake)+length(sleep));
end
allCJuxta.pbComboChange = comboChange;

cd(path)
allCJuxta = orderfields(allCJuxta);
save allCJuxta.mat allCJuxta
disp('saved allCJuxta')


%%

load([path 'allCJuxta.mat'])
pv = allCJuxta.identity == 'pv';
som = allCJuxta.identity == 'som';
reln = allCJuxta.identity == 'reln';
pvsom = allCJuxta.identity == 'pvsom';
none = allCJuxta.identity == 'nothing';
unfound = allCJuxta.identity ~= 'pv' & allCJuxta.identity ~= 'som' & allCJuxta.identity ~= 'reln' & allCJuxta.identity ~= 'pvsom' & allCJuxta.identity ~= 'nothing';


%% plot 3 vars & vid

wf = allCJuxta.wfRise;
modul = abs(allCJuxta.pbSleepChange); 
state = allCJuxta.spontStateCorr;
id = allCJuxta.clusterID;

fig1 = figure; 
%subplot(5,5,[1 2 6 7])
p = plot3(wf(unfound),modul(unfound),state(unfound),'o','MarkerFaceColor','none','MarkerEdgeColor',[.6 .6 .6]);
%alpha(p,.7)
hold on; plot3(wf(pv),modul(pv),state(pv),'o','MarkerFaceColor','b','MarkerEdgeColor',[.6 .6 .6])
hold on; plot3(wf(som),modul(som),state(som),'o','MarkerFaceColor','r','MarkerEdgeColor',[.6 .6 .6])
hold on; plot3(wf(reln),modul(reln),state(reln),'o','MarkerFaceColor','g','MarkerEdgeColor',[.6 .6 .6])
hold on; plot3(wf(none),modul(none),state(none),'o','MarkerFaceColor', 'k','MarkerEdgeColor',[.6 .6 .6])
hold on; plot3(wf(pvsom),modul(pvsom),state(pvsom),'o','MarkerFaceColor','c','MarkerEdgeColor',[.6 .6 .6])
hold on; text(wf(~unfound),modul(~unfound),state(~unfound),allCJuxta.clusterID(~unfound))
xlabel('rise')
ylabel('pb sleep change')
zlabel('ISIcorr')
% cd('C:\Users\User\Dropbox (NYU Langone Health)\Committee Meetings\Committee Meeting 7\')
% OptionZ.FrameRate=50;OptionZ.Duration=3;OptionZ.Periodic=true;
% CaptureFigVid([0,-10;-50,-20;-100,-30], 'juxta',OptionZ)
set(gca, 'TickDir','out','color','none','TickDir','out','box','off')
fig1.Renderer='Painters';
%print(['3vars everything'], '-dpdf');


%% cluster 3 vars

clc; 

% figure; 
% subplot(3,2,1)
% histogram(wf,20)
% subplot(3,2,3)
% histogram(modul,20)
% subplot(3,2,5)
% histogram(state,20)
% subplot(3,2,2)
% histogram(sqrt(wf),20)
% subplot(3,2,4)
% histogram(sqrt(modul),20)
% subplot(3,2,6)
% histogram(matchRatio(state),20)

use = [wf;modul;state]';
zuse = [zscorenan(wf);zscorenan(modul);matchRatio(state);zscorenan(allCJuxta.spontAwakeCV);zscorenan(sqrt(allCJuxta.spontAwakeBurst))]'; %[zscorenan(wf);zscorenan(modul);matchRatio(state)]';

symbols = 'svd*o.|*.|sv';
colors = 'gbrmkcgbmkrc';
meas = 12;

clear clus
for iter = 1:meas
    clus(:,iter) = kmeans(zuse,iter,'Replicates',1000);
end

clusData = evalclusters(zuse,clus,'CalinskiHarabasz');
k = clusData.OptimalK
whichClus = clus(:,k);

fig2 = figure; plot(clusData);
set(gca, 'TickDir','out','color','none','FontSize',15,'TickDir','out','box','off')

fig3 = figure; 
for i = 1:k
    pts = use(whichClus==i,:);
    plot3(pts(:,1),pts(:,2),pts(:,3),[colors(i) symbols(i)])
    hold on
end
set(gca, 'TickDir','out','color','none','FontSize',15,'TickDir','out','box','off')
fig3.Renderer='Painters';
% %print(['3d clusters'], '-dpdf');


%% compute subtype distances

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
    % rng('shuffle')
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


%% raksin

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
%whichpb(mod) = {'diff'};
whichpb(unmod) = {'same'};


%% plot params

% % cd('C:/Users/User/Dropbox (NYU Langone Health)/Thesis/Dissertation/figures/')

% % samesies = {'wfRise','wfFWHM','pbAwakeFR','pbAwakeChange','pbAwakeCorr','spontAwakeFR','spontAwakeCV',...
% %     'spontStateCorr','spontRatioVigilance','pbSleepFR','pbSleepChange','pbSleepCorr','spontSleepFR','spontSleepCV'};
% % 
% % thesisWords = {'WF Rise','WF FWHM', 'Playback FR_O_N','Playback FR Change_O_N','Playback Precision_O_N','Spont FR_O_N','Spont CV_O_N',...
% %     'ISI Corr','Vigilance Change','Playback FR_O_F_F','Playback FR Change_O_F_F','Playback Precision_O_F_F', 'Spont FR_O_F_F','Spont CV_O_F_F'};

% % use = nan(length(allC),length(samesies));

% % for i = 1:length(samesies)
% %     field = string(samesies(i));
% % %     figure('Visible','off');
% % %     plotVarJuxta(allCJuxta,field)
% % %     print(['juxta_'+ field], '-dpdf','-fillpage');
% % %     set(gcf,'PaperOrientation','landscape')
% %     pars = getfield(allCJuxta,field);  
% %     %if i == 1 | i == 2
% %     %   pars = log(pars);%nthroot(pars,3);
% %     if i == 4 | i == 11
% %        pars = abs(pars);
% % %        pars = log(pars);%nthroot(pars,3);
% %     elseif i == 3 | i == 6 | i == 10 | i == 13
% %        pars = log(pars);%nthroot(pars,3);
% %     elseif i == 5 | i == 8 | i == 12
% %         pars = matchRatio(pars);
% %     end
% %     pars(isinf(pars)) = nan;
% %     use(:,i) = pars;
% % end
% % 
% % useGood = zscorenan(use);

% % %use(2,7) = -1;
% % 
% % %useGood(2,7) = nan;
% % histOrder = [1,6,7,13,14,8,9,2,3,4,5,10,11,12];
% % 
% % figure;     
% % for i = 1:length(samesies)
% %     v = useGood(:,histOrder(i));
% %     subplot(7,2,i)
% %     histogram(v,30,'DisplayStyle','stairs','EdgeColor','k')
% %     xlim([-5 5])
% %     xticks(-4:2:4)
% %     set(gca,'box','off','TickDir','out','color','none')
% %     [h p] = swtest(v);
% %     pass = p < .05; 
% %     if pass == 1
% %         verdict = 'not normal';
% %     else 
% %         verdict = 'normal';
% %     end
% %     title(thesisWords(histOrder(i)))
% % end


%% make corr matrix of params

% % useGoodFull = useGood([1,3:26,28:30,33:35,37:40,42:49,52:end],:);
% % corrOrder = [1 2 6 7 13 14 8 9 3 4 5 10 11 12];
% % useCorr = useGoodFull(:,corrOrder);
% % [R,P] = corrcoef(useCorr);
% % 
% % % figure; imagesc(R); hold on;
% % % sig = string(nan(size(P)));
% % % sig(P<.05) = '*';
% % % nvars = length(samesies);
% % % for a = 1:nvars
% % %     i = a-1;
% % %     text([1:nvars],a*ones(nvars,1),sig((1:nvars)+nvars*i),'Color','k')
% % %     hold on;
% % %     shg;
% % % end
% % % 
% % % xticks(1:length(samesies))
% % % xticklabels(thesisWords(corrOrder))
% % % yticks(1:length(samesies))
% % % yticklabels(thesisWords(corrOrder))
% % % set(gca,'box','off','TickDir','out')
% % % colorbar 
% % 
% % useGood = useGood(:,corrOrder);


%% raksin comparison

% % change = allCJuxta.pbSleepChange;
% % for i = 1:length(raksin)
% %     pre = raksin(i).prePbSleepFRAll;
% %     post = raksin(i).pbSleepFRAll;
% %     if ~isempty(pre) && ~isempty(post)
% %         [pval(i) sig(i)] = signrank(pre,post);
% %     else
% %        pval(i) = nan;
% %      %  sig(i) = logical(nan);
% %     end
% % end
% % unmod = ~sig;
% % pos = (change > 0) & sig;
% % neg = (change < 0) & sig;
% % mod = sig;
% % 
% % clear whichpb
% % whichpb(pos) = {'up'};
% % whichpb(neg) = {'down'};
% % %whichpb(mod) = {'diff'};
% % whichpb(unmod) = {'same'};


%% pca
 
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
% % 
% % cumulsum = cumsum(explained);
% % figure; plot(cumulsum,'k-o')
% % set(gca,'box','off','TickDir','out','color','none')
% % 
% % elbow = find(cumulsum > 90, 1);
% % 
% % %sum(explained(1:elbow))
% % figure; plot(mean(abs(coeff(:,1:elbow)),2),'k-o')
% % set(gca,'box','off','TickDir','out','color','none')
% % 
% % xticks(1:length(samesies(1:end)))
% % xticklabels(thesisWords(corrOrder(1:end)))
% % shg


%% each pca
%  
% [coeff1,score1,latent1,~,explained1] = pca(useGood(:,1:2));
% [coeff2,score2,latent2,~,explained2] = pca(useGood(:,3:8));
% [coeff3,score3,latent3,~,explained3] = pca(useGood(:,9:14));
% 
% x = score1(:,1);
% y = score2(:,1);
% z = score3(:,1);
% d = .04;
% 
% pv = allCJuxta.identity == 'pv';
% som = allCJuxta.identity == 'som';
% reln = allCJuxta.identity == 'reln';
% pvsom = allCJuxta.identity == 'pvsom';
% none = allCJuxta.identity == 'nothing';
% unfound = allCJuxta.identity ~= 'pv' & allCJuxta.identity ~= 'som' & allCJuxta.identity ~= 'reln' & allCJuxta.identity ~= 'pvsom' & allCJuxta.identity ~= 'nothing';
% 
% figure;
% scatter3(x(unfound),y(unfound),z(unfound),'MarkerEdgeColor',[.5 .5 .5]);
% hold on;
% scatter3(x(pv),y(pv),z(pv),'filled','c');
% hold on; 
% scatter3(x(som),y(som),z(som),'filled','r');
% hold on
% scatter3(x(reln),y(reln),z(reln),'filled','g');
% hold on
% scatter3(x(none),y(none),z(none),'k');
% %scatter3(x(pvsom),y(pvsom),z(pvsom),'filled','m');
% xlabel('pc1')
% ylabel('pc2')
% zlabel('pc3')
% grid off;
% xticks(-4:2:4)
% yticks(-4:2:4)
% zticks(-4:2:4)
% set(gca,'box','off','TickDir','out','color','none')


%% introduce x projectors

% Xp = [];
% for i = 1:height(goodCells)
%     if toPrint(i) == 1 && badx(i)
%         changedir = [path + birds(i) + '/' + cells(i) + '/'];
%         load([changedir + 'clusters.mat']);
%         Xp = [Xp clusters];
%     end
% end
% 
% xuse = [Xp.wfFWHM; abs([Xp.pbSleepChange]); Xp.spontStateCorr]';
% xuse(1,2) = 0;
% 
% hold on; plot3(xuse(:,1),xuse(:,2),xuse(:,3),'ko','MarkerFaceColor','k')
% 

%% cluster graphs

% options = statset('MaxIter',1000);
% gmm = fitgmdist(use,4,'RegularizationValue', .05, 'Options', options);
% clusterIDs = cluster(gmm, zuse);
% dens = dbscan(zuse,2,5);
% use = zuse;
% figure;
% %for i = 1:
% scatter3(use(unfound,1),use(unfound,2),use(unfound,3),'k');
% hold on
% scatter3(use(pv,1),use(pv,2),use(pv,3),'b','filled')
% hold on
% scatter3(use(pvsom,1),use(pvsom,2),use(pvsom,3),'m','filled')
% hold on
% scatter3(use(som,1),use(som,2),use(som,3),'r','filled')
% hold on
% scatter3(use(reln,1),use(reln,2),use(reln,3),'g','filled')
% hold on
% d = .01;
% text(use(:,1)+d,use(:,2)+d,use(:,3)+d,string(whichClus))
% 

%% pca 3 vars
% 
% [coeff,score,latent,~,explained] = pca(use);
% figure; scatter(score(:,1),score(:,2))
% hold on; scatter(score(pv,1),score(pv,2),'b','filled')
% hold on; scatter(score(som,1),score(som,2),'r','filled')
% hold on; scatter(score(reln,1),score(reln,2),'g','filled')
% 
% groups = kmeans(use,6);
% figure;
% gscatter(score(:,1),score(:,2),groups);
% hold on
% %text(x,y,groups,ints)
% set(gca, 'TickDir','out','color','none','FontSize',9,'TickDir','out')
% xlabel('pc1')
% ylabel('pc2')


%% hierarchical clustering

% %account for pos deflecting spikes!
% 
% transform = struct2cell(allCJuxta);
% use = cell2mat(transform(2:end))';
% clear transform;
% use(isinf(use)) = nan;
% tree = linkage(use(~deflect,:),'ward',@naneucdist);
% 
% figure;
% subplot(1,2,1)
% [~,~,outperm] = dendrogram(tree,0,'Orientation','left','Labels',identity(~deflect),'ColorThreshold','default');
% set(gca, 'TickDir', 'out','color','none','box','off','FontSize',9)
% title('all params')
% 
% params2 = erase(params,{'spontRatioVigilance';'spontAwakeBurst';'spontAwakeCV';...
% 'pbAwakeBurst';'pbAwakePrecInd10';'pbAwakeMaxIFR';'pbAwakeChange';'pbAwakePostChange';...
% 'pbSleepChange';'pbSleepPostChange'});
% params2(ismissing(params2)) = [];
% 
% params3 = {'wfRise'}%,'wfRise70','pbSleepFR'%'pbSleepFR','pbSleepFRPost','pbSleepPrecInd50','spontAwakeFR','spontSleepFR','pbAwakeFR'}
% 
% clear trunk;
% for k=1:length(params3)
%     trunk.(params3{k}) = [allCJuxta.(params3{k})];
% end
% use = cell2mat(struct2cell(trunk))';
% use(isinf(use)) = nan;
% tree = linkage(use(~deflect,:),'ward',@naneucdist);
% 
% subplot(1,2,2)
% [~,~,outperm] = dendrogram(tree,0,'Orientation','left','Labels',identity(~deflect),'ColorThreshold','default');
% set(gca, 'TickDir', 'out','color','none','box','off','FontSize',9)
% title('rise time only')


%% fix pos spks
% 
% load([path 'wfs/posvals.mat'])
% diffRise(isnan(diffRise))=[];
% %diffRise(deflect) = diffRise(deflect)*2;
% diffRise70(isnan(diffRise70))=[];
% diffWidth(isnan(diffWidth))=[];
% 
%
%% mean normalized waveforms
% 
% for j = 1:length(allC)
%     if deflect(j)
%         allC(j).meanwf = allC(j).meanwf*-1;
%        % allC(j).meanwf = nan;
%     end
%     normalized = allC(j).meanwf/-min(allC(j).meanwf);
%     troughPt = find(normalized == min(normalized));
%     allC(j).normwf = normalized(troughPt-40:troughPt+40);
% end
% 
% pv = allC(identity == 'pv');
% pvsom = allC(identity == 'pvsom');
% som = allC(identity == 'som');
% 
% nothing = allC(identity == 'nothing');
% x = allC(identity == 'x');
% 
% 
% %% isi histogram comparisons
% 
% pvn = 0;
% pvsomn = 0;
% nothingn = 0;
% xn = 0;
% somn = 0;
% for n = 1:length(allC)
%     ons = allC(n).onIsis(allC(n).onIsis>=1e-3 & allC(n).onIsis<=10);
%     offs = allC(n).offIsis(allC(n).offIsis>=1e-3 & allC(n).offIsis<=10);
%     alloni{n} = ons;
%     alloffi{n} = offs;
%     pbOni{n} = allC(n).onPbIsis(allC(n).onPbIsis>=1e-3 & allC(n).onPbIsis<=1.5);
%     pbOffi{n} = allC(n).offPbIsis(allC(n).offPbIsis>=1e-3 & allC(n).offPbIsis<=1.5);
%     if identity(n) == 'pv'
%         pvn = pvn+1;
%         pvoni{pvn,:} = ons;
%         pvoffi{pvn,:} = offs;
%     end
%     if identity(n) == 'som'
%         somn = somn+1;
%         somoni{somn,:} = ons;
%         somoffi{somn,:} = offs;
%     end
%     if identity(n) == 'pvsom'
%         pvsomn = pvsomn+1;
%         pvsomoni{pvsomn,:} = ons;
%         pvsomoffi{pvsomn,:} = offs;
%     end
%     if identity(n) == 'nothing'
%         nothingn = nothingn + 1;
%         nothingoni{nothingn,:} = ons;
%         nothingoffi{nothingn,:} = offs;
%     end
%     if identity(n) == 'x'
%         xn = xn + 1;
%         xoni{xn,:} = ons;
%         xoffi{xn,:} = offs;
%     end
% end
% 
% bw = .01;
% fa = .5;
% edges = -3:.1:1;
% 
% 
% %% corr scores
% 
% for i = 1:length(allC)
%     allonih{i} = histcounts(log10(alloni{i}),edges,'Normalization','probability');
%     alloffih{i} = histcounts(log10(alloffi{i}),edges,'Normalization','probability');
%     c = corrcoef(allonih{i}, alloffih{i});
%     stateCorr(i) = c(2,1);
% 
%     pbonih{i} = histcounts(log10(pbOni{i}),edges,'Normalization','probability');
%     d = corrcoef(allonih{i},pbonih{i});
%     stateCorrOnPb(i) = d(2,1);
%     
%     pboffih{i} = histcounts(log10(pbOffi{i}),edges,'Normalization','probability');
%     e = corrcoef(alloffih{i},pboffih{i});
%     stateCorrOffPb(i) = e(2,1);
% end
% pvCorr = stateCorr(identity=='pv')';
% somCorr = stateCorr(identity=='som')';
% xCorr = stateCorr(identity=='x')';
% pvsomCorr = stateCorr(identity=='pvsom')';
% nothingCorr = stateCorr(identity=='nothing')';
% 
% 
% %% rise times
% 
% fs = 50000;
% for i = 1:length(allC)
%   %  if ~deflect(i)
%         wave = allC(i).normwf;
%         trPos = find(wave == min(wave));
%         pkPos = find(wave == max(wave(trPos:end)));
%         trPk(i) = (pkPos-trPos)/fs*1000; 
%         allC(i).rise = trPk(i);
%   %  end
% end
% 
% pvRise = trPk(identity == 'pv')';
% somRise = trPk(identity == 'som')';
% xRise = trPk(identity == 'x')';
% pvsomRise = trPk(identity == 'pvsom')';
% nothingRise = trPk(identity == 'nothing')';
% %trPk(trPk == 0) = nan;
% allRise = trPk;
% 
% 
% %% pb change
% 
% pvPbFR = [pv.pbOnFR];
% pvSpontFR = [pv.spontOnFR];
% pvChange = abs(pvPbFR-pvSpontFR)./pvSpontFR;
% 
% somPbFR = [som.pbOnFR];
% somSpontFR = [som.spontOnFR];
% somChange = abs(somPbFR-somSpontFR)./somSpontFR;
% 
% xPbFR = [x.pbOnFR];
% xSpontFR = [x.spontOnFR];
% xChange = abs(xPbFR-xSpontFR)./xSpontFR;
% 
% pvsomPbFR = [pvsom.pbOnFR];
% pvsomSpontFR = [pvsom.spontOnFR];
% pvsomChange = abs(pvsomPbFR-pvsomSpontFR)./pvsomSpontFR;
% 
% nothingPbFR = [nothing.pbOnFR];
% nothingSpontFR = [nothing.spontOnFR];
% nothingChange = abs(nothingPbFR-nothingSpontFR)./nothingSpontFR;
% 
% allPbFR = [allC.pbOffFR];
% allSpontFR = [allC.spontOffFR];
% pbChange = abs(allPbFR - allSpontFR)./allSpontFR;
% pvChange = pbChange(identity=='pv');
% 
% 
% %% pb psth
% 
% pts = 6;
% pbOnFR = [];
% for i = 1:length(allC)
%     onPsth = allC(i).fullPsthOn;
%     len = floor(length(onPsth)/2);
%     onPre{i,:} = onPsth(1:len-pts);
%     onPost{i,:} = onPsth(len+pts:length(onPsth));
%     if isempty(onPre{i,:})
%         onBaseFR = 0;
%     else
%         onBaseFR = mean(onPre{i,:});
%     end
%     if isempty(onPost{i,:})
%         onPBFR = 0;
%     else
%         onPbFR = mean(onPost{i,:});
%     end
%     pbOnFR = [pbOnFR; onPbFR];
%     allChangeOn(i) = abs(onPbFR-onBaseFR)/onBaseFR;
%     normPostOn{i,:} = onPsth-onBaseFR;
% end
% 
% pbOffFR = [];
% for i = 1:length(allC)
%     offPsth = allC(i).fullPsthOff;
%     len = floor(length(offPsth)/2);
%     offPre{i,:} = offPsth(1:len-pts);
%     offPost{i,:} = offPsth(len+pts:length(offPsth));
%     if isempty(offPre{i,:})
%         offBaseFR = 0;
%     else
%         offBaseFR = mean(offPre{i,:});
%     end
%     if isempty(offPost{i,:})
%         offPBFR = 0;
%     else
%         offPbFR = mean(offPost{i,:});
%     end
%     pbOffFR = [pbOffFR; offPbFR];
%     allChangeOff(i) = abs(offPbFR-offBaseFR)/offBaseFR;
%     normPostOff{i,:} = offPsth-offBaseFR;
% end
% 
% pvChangeOn = allChangeOn(identity == 'pv');
% pvsomChangeOn = allChangeOn(identity == 'pvsom');
% somChangeOn = allChangeOn(identity == 'som');
% nothingChangeOn = allChangeOn(identity == 'nothing');
% xChangeOn = allChangeOn(identity == 'x');
% 
% pvChangeOff = allChangeOff(identity == 'pv');
% pvsomChangeOff = allChangeOff(identity == 'pvsom');
% somChangeOff = allChangeOff(identity == 'som');
% nothingChangeOff = allChangeOff(identity == 'nothing');
% xChangeOff = allChangeOff(identity == 'x');
% 
% pbChange = nanmean([allChangeOn; allChangeOff]);
% pvChange = nanmean([pvChangeOn; pvChangeOff]);
% pvsomChange = nanmean([pvsomChangeOn; pvsomChangeOff]);
% somChange = nanmean([somChangeOn; somChangeOff]);
% nothingChange = nanmean([nothingChangeOn; nothingChangeOff]);
% xChange =  nanmean([xChangeOn; xChangeOff]);
% 
% 
% %% normalize data
% 
% %zstateCorr(isnan(zstateCorr)) = nanmean(stateCorr);
% %spontOnFR = [allC.spontOnFR];
% %spontOnFR(isnan(spontOnFR))=0;
% %spontOffFR = [allC.spontOffFR];
% %spontOffFR(isnan(spontOffFR))=0;
% spontOnFR = nan(length(allC),1);
% spontOnFR = nan(length(allC),1);
% for i = 1:length(allC)
%     spontOnFR(i) = 1./nanmean([allC(i).onIsis]);
%     spontOffFR(i) = 1./nanmean([allC(i).offIsis]);
% end
% zspontOnFR = (spontOnFR-nanmean(spontOnFR))/nanstd(spontOnFR);
% zspontOffFR = (spontOffFR-nanmean(spontOffFR))/nanstd(spontOffFR);
% 
% spontOnBurst = [allC.spontOnBurst];
% %spontOnBurst(isnan(spontOnBurst))=0;
% zspontOnBurst = (spontOnBurst-nanmean(spontOnBurst))/nanstd(spontOnBurst);
% 
% spontOffBurst = [allC.spontOffBurst];
% %spontOffBurst(isnan(spontOffBurst))=0;
% zspontOffBurst = (spontOffBurst-nanmean(spontOffBurst))/nanstd(spontOffBurst);
% 
% pbOnBurst = [allC.pbOnBurst];
% %pbOnBurst(isnan(pbOnBurst))=0;
% zpbOnBurst = (pbOnBurst-nanmean(pbOnBurst))/nanstd(pbOnBurst);
% 
% pbOffBurst = [allC.pbOffBurst];
% %pbOffBurst(isnan(pbOffBurst))=0;
% zpbOffBurst = (pbOffBurst-nanmean(pbOffBurst))/nanstd(pbOffBurst);
% 
% % spontOnCV = [allC.spontOnCV];
% % zspontOnCV = (spontOnCV-nanmean(spontOnCV))/nanstd(spontOnCV);
% % spontOffCV = [allC.spontOffCV];
% % zspontOffCV = (spontOffCV-nanmean(spontOffCV))/nanstd(spontOffCV);
% % 
% % pbOnCV = [allC.pbOnCV];
% % zpbOnCV = (pbOnCV-nanmean(pbOnCV))/nanstd(pbOnCV);
% % pbOffCV = [allC.pbOffCV];
% % zpbOffCV = (pbOffCV-nanmean(pbOffCV))/nanstd(pbOffCV);
% 
% pbOnPrec = [allC.pbOnPrec];
% %pbOnPrec(isnan(pbOnPrec))=0;
% pbOnPrec(isinf(pbOnPrec))=0;
% zpbOnPrec = (pbOnPrec-nanmean(pbOnPrec))/nanstd(pbOnPrec);
% zpbOnPrec(isnan(zpbOnPrec)) = nanmean(zpbOnPrec);
% pbOffPrec = [allC.pbOffPrec];
% %pbOffPrec(isnan(pbOffPrec))=0;
% pbOffPrec(isinf(pbOffPrec))=0;
% zpbOffPrec = (pbOffPrec-nanmean(pbOffPrec))/nanstd(pbOffPrec);
% zpbOffPrec(isnan(zpbOffPrec)) = nanmean(zpbOffPrec);
% 
% ratioVigi = [allC.ratioVigilance];
% zratioVigi = zscore(ratioVigi);
% 
% spontFR = [allC.spontWtFr];
% zspontFR = (spontFR-nanmean(spontFR))/nanstd(spontFR);
% zspontFR(isnan(zspontFR)) = nanmean(spontFR);
% 
% pbFR = [allC.pbWtFr];
% zpbFR = (spontFR-nanmean(pbFR))/nanstd(pbFR);
% zpbFR(isnan(zpbFR)) = nanmean(pbFR);
% 
% pbOnPost = [allC.pbOnPost];
% zpbOnPost = zscore(pbOnPost);
% 
% pbOffPost = [allC.pbOffPost];
% zpbOffPost = zscore(pbOffPost);
% 
% pbPost = nanmean([pbOffPost; pbOnPost]);
% zpbPost = zscore(pbPost);
% 
% ratioOn = [allC.ratioPbSpontOn];
% ratioOn(isinf(ratioOn)) = nan;
% 
% ratioOff = [allC.ratioPbSpontOff];
% ratioOff(isinf(ratioOff)) = nan;
% 
% ratioMean = nanmean([ratioOn; ratioOff]);
% zratioMean = (ratioMean-nanmean(ratioMean))/nanstd(ratioMean);
% zratioMean(isnan(zratioMean)) = nanmean(zratioMean);
% 
% zallRise = zscore(allRise);
% zallChange = zscore(pbChange);
% zstateCorr = (stateCorr-nanmean(stateCorr))/nanstd(stateCorr);
% zstateCorr(isnan(zstateCorr)) = nanmean(zstateCorr);
% 
% 
% %% 3var plot NEW!
% 
% ratioMean(isnan(ratioMean)) = nan;
% 
% close all;
% 
% figure;
% l = 40;
% for i = 1:length(allC)
%     scatter3(allRise(i),ratioMean(i),stateCorr(i),'k');
%     hold on
%     %text(allRise(i),pbChange(i),stateCorr(i),string(i),'FontSize',fontsize,'VerticalAlignment', 'top','HorizontalAlignment','left');
%     hold on
% end
% 
% %xlabel('Trough-to-Peak Time (ms)')
% %ylabel('\Delta FR during playback')
% %zlabel('ISI Corr: On vs. Off')
% %allParams = [rise; change; corr]';
% set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% grid off
% %title('3 variable space')
% 
% hold on
% scatter3(allRise(identity == 'pv'),ratioMean(identity == 'pv'),stateCorr(identity == 'pv'),l,'b','filled')
% hold on
% scatter3(allRise(identity == 'pvsom'),ratioMean(identity == 'pvsom'),stateCorr(identity == 'pvsom'),l,'m','filled')
% hold on
% scatter3(allRise(identity == 'som'),ratioMean(identity == 'som'), stateCorr(identity == 'som'),l,'r','filled')
% hold on
% scatter3(allRise(identity == 'x'),ratioMean(identity == 'x'),stateCorr(identity == 'x'),l,[.4 .4 .4],'filled')
% hold on
% 
% 
% %% lda
% 
% use = [zallRise; zratioMean; zstateCorr; zspontFR; zratioVigi; zpbOnPrec; zpbOffPrec; zpbPost]'; 
% %use = [zallRise; zratioMean; zstateCorr; zspontFR]'; 
% params = {'riseTime';'pbMod';'stateCorr';'spontFR';'stateRatio';'precisionOn';'precisionOff'; 'pbPost'};
% which = identity;
% which(which == '') = '?';
% lin = fitcdiscr(use,which,'ClassNames',{'pv','som','pvsom','x','nothing'},'OptimizeHyperparameters','auto','DiscrimType','auto','DiscrimType','linear')
% [label,score,cost] = predict(lin,use);
% result = [identity label]
% 
% 
%% pca
% 
% % [coeff,score,latent,~,explained] = pca(use);
% % 
% % % figure;
% % % scatter3(score(:,1),score(:,2),score(:,3),'ko')
% % % hold on
% % % for i = 1:length(allC)
% % %     hold on
% % %     if identity(i) == 'pv'
% % %         scatter3(score(i,1),score(i,2),score(i,3),'b','filled')
% % %         hold on
% % %     elseif identity(i) == 'som'
% % %         scatter3(score(i,1),score(i,2),score(i,3),'r','filled')
% % %         hold on
% % % %     elseif identity(i) == 'pvsom'
% % % %         scatter3(score(i,1),score(i,2),score(i,3),'m','filled')
% % % %         hold on
% % %     elseif identity(i) == 'x'
% % %         scatter3(score(i,1),score(i,2),score(i,3),'k','filled')
% % %         hold on
% % %     end
% % %     hold on
% % % end
% % % grid off
% % % xlabel('pc1')
% % % ylabel('pc2')
% % % zlabel('pc3')
% % 
% % % figure;
% % % plot(explained,'k-o')
% % 
% % options = statset('MaxIter',1000);
% % gmm = fitgmdist(use,5,'RegularizationValue', .05, 'Options', options);
% % clusterIDs = cluster(gmm, use);
% % 
% % close all;
% % 
% % subplot(1,2,1)
% % gscatter(score(:,1),score(:,2),clusterIDs);
% % legend('Cluster1','Cluster2','Cluster3','Cluster4','Location','best');
% % 
% % subplot(1,2,2)
% % scatter(score(:,1),score(:,2),'ko')
% % for i = 1:length(allC)
% %     hold on
% %     if identity(i) == 'pv'
% %         scatter(score(i,1),score(i,2),'b','filled')
% %         hold on
% %     elseif identity(i) == 'som'
% %         scatter(score(i,1),score(i,2),'r','filled')
% %         hold on
% %     elseif identity(i) == 'pvsom'
% %         scatter(score(i,1),score(i,2),'m','filled')
% %         hold on
% %     elseif identity(i) == 'x'
% %         scatter(score(i,1),score(i,2),'k','filled')
% %         hold on
% %     end
% %     hold on 
% % end
% 
% 
% %% tsne plot
% 
% % close all;
% % 
% % Y2 = tsne(use,'NumDimensions',3)%,'Perplexity',2, 'Exaggeration',2);
% % figure; scatter3(Y2(:,1),Y2(:,2),Y2(:,3),15,'k')
% % hold on
% % scatter3(Y2(identity == 'pv',1),Y2(identity == 'pv',2),Y2(identity == 'pv',3),15,'b','filled')
% % hold on
% % %scatter3(Y2(identity == 'pvsom',1),Y2(identity == 'pvsom',2),Y2(identity == 'pvsom',3),15,'m','filled')
% % hold on
% % scatter3(Y2(identity == 'som',1),Y2(identity == 'som',2),Y2(identity == 'som',3),15,'r','filled')
% % hold on
% % scatter3(Y2(identity == 'x',1),Y2(identity == 'x',2),Y2(identity == 'x',3),15,'k','filled')
% % hold on
% 
% 
% %% plot individual variables
% 
% % plotVarJuxta(allRise,.02,identity)
% % title('allRise')
% % xlabel('ms')
% % 
% % plotVarJuxta(pbChange,.1,identity)
% % title('pb modulation')
% % xlabel('deltaFR/spontFR')
% % 
% % plotVarJuxta(stateCorr,.02,identity)
% % title('corr of spont on/off ISIs')
% % xlabel('r value')
% 
% % plotVarJuxta(spontOnFR,.5,identity)
% % title('on spont FR')
% % xlabel('Hz')
% % 
% % plotVarJuxta(spontOffFR,.5,identity)
% % title('off spont FR')
% % xlabel('Hz')
% % 
% % plotVarJuxta(spontOnBurst,1,identity)
% % title('spont on burstiness')
% % xlabel('percent bursty ISIs')
% % 
% % plotVarJuxta(spontOffBurst,1,identity)
% % title('spont off burstiness')
% % xlabel('percent bursty ISIs')
% % 
% % plotVarJuxta(pbOnBurst,1,identity)
% % title('pb on burstiness')
% % xlabel('percent bursty ISIs')
% % 
% % plotVarJuxta(pbOffBurst,1,identity)
% % title('pb off burstiness')
% % xlabel('percent bursty ISIs')
% % 
% % plotVarJuxta(pbOnFR,1,identity)
% % title('on playback FR')
% % xlabel('Hz')
% % 
% % plotVarJuxta(pbOffFR,1,identity)
% % title('off playback FR')
% % xlabel('Hz')
% % 
% % plotVarJuxta(pbOnPrec,.1,identity)
% % title('precision on')
% % xlabel('prec index on')
% % 
% % plotVarJuxta(pbOffPrec,.1,identity)
% % title('precision off')
% % xlabel('prec index off')
% % 
% % plotVarJuxta(stateCorrOnPb,.05,identity)
% % title('corr pb vs spont')
% % xlabel('r value')
% % 
% % plotVarJuxta(stateCorrOffPb,.05,identity)
% % title('corr pb vs spont')
% % xlabel('r value')
% % 
% % plotVarJuxta(ratioOn,.05,identity)
% % title('ratio playback to spont: on')
% % xlim([0 15])
% % 
% % plotVarJuxta(ratioOff,.05,identity)
% % title('ratio playback to spont: off')
% % 
% % plotVarJuxta(ratioMean,.05,identity)
% % title('mean ratio playback to spont')
% % 
% % plotVarJuxta(ratioVigi,.05,identity)
% % title('ratio vigilance')
% % xlim([0 6])
% % 
% % plotVarJuxta(spontFR,.05,identity)
% % title('spont wt fr')
% % 
% % plotVarJuxta(pbFR,.05,identity)
% % title('pb wt fr')
% % 
% % plotVarJuxta(pbOffPost,.05,identity)
% % title('pb off post')
% % 
% % plotVarJuxta(spontOnCV,.05,identity)
% % title('pb off cv')
% 
% 
% %% 3var plot
% % 
% % figure;
% % l = 40;
% % for i = 1:length(allC)
% %     scatter3(allRise(i),pbChange(i),stateCorr(i),'k');
% %     hold on
% %       %text(allRise(i),pbChange(i),stateCorr(i),string(i),'FontSize',fontsize,'VerticalAlignment', 'top','HorizontalAlignment','left');
% %     hold on
% % end
% % hold on
% % scatter3(allRise(identity == 'pv'),pbChange(identity == 'pv'),stateCorr(identity == 'pv'),l,'b','filled')
% % hold on
% % scatter3(allRise(identity == 'pvsom'),pbChange(identity == 'pvsom'),stateCorr(identity == 'pvsom'),l,'m','filled')
% % hold on
% % scatter3(allRise(identity == 'som'),pbChange(identity == 'som'), stateCorr(identity == 'som'),l,'r','filled')
% % hold on
% % scatter3(allRise(identity == 'x'),pbChange(identity == 'x'),stateCorr(identity == 'x'),l,'k','filled')
% % hold on
% % 
% % xlabel('Trough-to-Peak Time (ms)')
% % ylabel('\Delta FR during playback')
% % zlabel('ISI Corr: On vs. Off')
% % %allParams = [rise; change; corr]';
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% % grid off
% % title('3 variable space')
% % 
% % 
% %% plot mean normalized waveforms
% %
% % figure;
% % subplot(5,1,1)
% % for i = 1:length(pv)
% %     plot(pv(i).normwf,'k')
% %     hold on
% % end
% % hold on
% % plot(mean([pv.normwf],2),'b','MarkerSize',3)
% % title('pv')
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% %
% % subplot(5,1,2)
% % for i = 1:length(pvsom)
% %     plot(pvsom(i).normwf,'k')
% %     hold on
% % end
% % hold on
% % plot(mean([pvsom.normwf],2),'m','MarkerSize',3)
% % title('pvsom')
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% %
% % subplot(5,1,3)
% % for i = 1:length(som)
% %     plot(som(i).normwf,'k')
% %     hold on
% % end
% % hold on
% % plot(mean([som.normwf],2),'r','MarkerSize',3)
% % title('som')
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% %
% % subplot(5,1,4)
% % for i = 1:length(nothing)
% %     plot(nothing(i).normwf,'k')
% %     hold on
% % end
% % hold on
% % plot(mean([nothing.normwf],2),'g','MarkerSize',3)
% % title('nothing')
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% %
% % subplot(5,1,5)
% % for i = 1:length(x)
% %     plot(x(i).normwf,'k')
% %     hold on
% % end
% % hold on
% % plot(mean([x.normwf],2),'Color','[.5. .5 .5]','MarkerSize',3)
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% %
% %
% %% spont on vs off
% %
% % figure;
% % title('spontOnFR vs. spontOffFR')
% % onoff = [ones(length(pv),1),ones(length(pv),1)*2]';
% % plot(onoff,[pv.spontOnCV; pv.spontOffCV],'b-o')
% % xlim([0 11])
% % ylim([-1 20])
% % hold on
% %
% % onoff = [ones(length(pvsom),1)*3,ones(length(pvsom),1)*4]';
% % plot(onoff,[pvsom.spontOnBurst; pvsom.spontOffBurst],'m-o')
% % hold on
% %
% % onoff = [ones(length(som),1)*5,ones(length(som),1)*6]';
% % plot(onoff,[som.spontOnCV; som.spontOffCV],'r-o')
% % hold on
% %
% % onoff = [ones(length(nothing),1)*7,ones(length(nothing),1)*8]';
% % plot(onoff,[nothing.spontOnFR; nothing.spontOffFR],'g-o')
% % hold on
% %
% % onoff = [ones(length(x),1)*9,ones(length(x),1)*10]';
% % plot(onoff,[x.spontOnCV; x.spontOffCV],'k-o')
% % shg
% %
% %
% %% spont on vs pb on
% %
% % figure;
% % title('spontOnFR vs. pbOnFR')
% % onoff = [ones(length(pv),1),ones(length(pv),1)*2]';
% % plot(onoff,[pv.spontOnFR; pv.pbOnFR],'b-o')
% % xlim([0 11])
% % ylim([-1 20])
% % hold on
% %
% % onoff = [ones(length(pvsom),1)*3,ones(length(pvsom),1)*4]';
% % plot(onoff,[pvsom.spontOnFR; pvsom.pbOnFR],'m-o')
% % hold on
% %
% % onoff = [ones(length(som),1)*5,ones(length(som),1)*6]';
% % plot(onoff,[som.spontOnFR; som.pbOnFR],'r-o')
% % hold on
% %
% % onoff = [ones(length(nothing),1)*7,ones(length(nothing),1)*8]';
% % plot(onoff,[nothing.spontOnFR; nothing.pbOnFR],'g-o')
% % hold on
% %
% % onoff = [ones(length(x),1)*9,ones(length(x),1)*10]';
% % plot(onoff,[x.spontOnFR; x.pbOnFR],'k-o')
% % shg
% %
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% %
% %
% %% plot isi histograms of subtypes
% %
% % figure;
% % title('ISIs')
% % h(1) = subplot(5,1,1);
% % pvih = histogram(log10(pvi),10.^edges);
% % pvih.Normalization = 'probability';
% % %pvih.edgesWidth = bw;
% % pvih.FaceAlpha = fa;
% % pvih.FaceColor = 'b';
% % title('ISIs')
% %
% % h(2) = subplot(5,1,2);
% % pvsomih = histogram(log10(pvsomi),10.^edges);
% % pvsomih.Normalization = 'probability';
% % %pvsomih.edgesWidth = bw;
% % pvsomih.FaceAlpha = fa;
% % pvsomih.FaceColor = 'm';
% %
% % h(3) = subplot(5,1,3);
% % somih = histogram(log10(somi),10.^edges);
% % somih.Normalization = 'probability';
% % %somih.edgesWidth = bw;
% % somih.FaceAlpha = fa;
% % somih.FaceColor = 'r';
% %
% % h(4) = subplot(5,1,4);
% % nothingih = histogram(log10(nothingi),10.^edges);
% % nothingih.Normalization = 'probability';
% % %nothingih.edgesWidth = bw;
% % nothingih.FaceAlpha = fa;
% % nothingih.FaceColor = 'g';
% %
% % h(5) = subplot(5,1,5);
% % xih = histogram(log10(xi),10.^edges);
% % xih.Normalization = 'probability';
% % %xih.edgesWidth = bw;
% % xih.FaceAlpha = fa;
% % xih.FaceColor = 'k';
% % xlabel('s')
% % linkaxes(h,'xy')
% % xlim([-bw .5])
% %
% %
% %% pca on pb psths
% %
% % close all;
% %
% % for i = 1:length(allC)
% %     onPsth = allC(i).fullPsthOn;
% %     offPsth = allC(i).fullPsthOff;
% %     smoothOn = smoothdata(onPsth,'movmean',5);
% %     smoothOff = smoothdata(offPsth,'movmean',5);
% %     figure(11)
% %     clf;
% %     subplot(2,1,1)
% %     plot(onPsth,'k')
% %     hold on
% %     plot(smoothOn,'r')
% %     subplot(2,1,2)
% %     plot(offPsth,'k')
% %     hold on
% %     plot(smoothOff,'r')
% %     pause(.5)
% % end
% %
% %
% % %% pca plot
% % %
% % use = [zallRise; zRatioMean; zstateCorr; zpbOnPrec];
% % whos use
% % [coeff,score,latent,~,explained] = pca(use);
% % 

%scatter3(score(41,1),score(41,2),score(41,3),l,'c','filled');
% 
%  
% %% zscored 3var plot
% % 
% % figure;
% % l = 40;
% % for i = 1:length(allC)
% %     scatter3(zallRise(i),zallChange(i),zstateCorr(i),'k');
% %     hold on
% %     %  text(zallRise(i),zallChange(i),zstateCorr(i),string(i),'FontSize',fontsize,'VerticalAlignment', 'top','HorizontalAlignment','left');
% %     hold on
% % end
% % hold on
% % scatter3(zallRise(identity == 'pv'),zallChange(identity == 'pv'),zstateCorr(identity == 'pv'),l,'b','filled')
% % hold on
% % %scatter3(zallRise(identity == 'pvsom'),zallChange(identity == 'pvsom'),zstateCorr(identity == 'pvsom'),l,'m','filled')
% % hold on
% % scatter3(zallRise(identity == 'som'),zallChange(identity == 'som'),zstateCorr(identity == 'som'),l,'r','filled')
% % hold on
% % scatter3(zallRise(identity == 'x'),zallChange(identity == 'x'),zstateCorr(identity == 'x'),l,'k','filled')
% % hold on
% % %scatter3(zallRise(41),zallChange(41),zstateCorr(41),l,'c','filled');
% % 
% % xlabel('Trough-to-Peak Time (ms)')
% % ylabel('\Delta FR during playback')
% % zlabel('ISI Corr: On vs. Off')
% % %allParams = [rise; change; corr]';
% % set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize)
% % grid off
% % title('zscored 3 variable space')
% %
% %
% 
% %%
% 
% % % % ind = [ones(length(pvRise),1);2*ones(length(somRise),1);3*ones(length(xRise),1);4*ones(length(pvsomRise),1);5*ones(length(nothingRise),1)];
% % % % 
% % % % foundRise = [pvRise;somRise;xRise;pvsomRise;nothingRise];
% % % % [p,tbl,stats] = anova1(foundRise,ind)
% % % % multcompare(stats)
% % % % 
% % % % 
% % % % foundCorr = [pvCorr;somCorr;xCorr;pvsomCorr;nothingCorr];
% % % % [p,tbl,stats] = anova1(foundCorr,ind)
% % % % multcompare(stats)
% % % % 
% % % % 
% % % % foundPbMod = [pvChange';somChange';xChange';pvsomChange';nothingChange'];
% % % % [p,tbl,stats] = anova1(foundPbMod,ind)
% % % % multcompare(stats)
% % % % 
% 
% %%


%% stats

clc;
var = (allCJuxta.pbSleepCorr); 
disp(['pv: ' + string(mean(var(pv)))]);
disp(['som: ' + string(mean(var(som)))]);
disp(['reln: ' + string(var(reln))]);
disp(['pvsom: ' + string(mean(var(pvsom)))]);


