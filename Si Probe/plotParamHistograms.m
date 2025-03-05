clear all; close all; clc; warning off; 
set(0,'DefaultFigureVisible','off');

fs = 30000;
fontsize = 15;


%% save 'em
allbirds = {'C22'; 'C24'; 'C50'; 'C52'};
allGoodCells = [];
allC = [];

for i = 1:length(allbirds)
    whichBird = string(allbirds(i));
    path = ['C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/' + whichBird + '/'];
    goodCells = readtable([path + whichBird + '_goodSiUnits.xlsx']);
    take = logical(goodCells{:,3});
    allGoodCells = [allGoodCells; goodCells(take,1:5)];
    load([path + '/data/clusters.mat'])
    allC = [allC clusters(take)];
end

allC = orderfields(allC);

birds = string([allGoodCells{:,1}]);

shitfields = {'spikeTimes','coordinates','shank','maxChannelTemplate','maxChannelWF',...
    'isis','ogTimes','spontIsis','spontSpks','spontOnSpks','spontOnIsis',...
    'spontOffSpks','spontOffIsis','smIsis','pmIsis','pmsIsis','meanWF'...
    'wfNorm','spontOffDeepIsis','spontOffEarlyIsis','amps'};

allC = rmfield(allC,shitfields);
allC = orderfields(allC);
strOG = fieldnames(allC);
str = strOG;

for k=1:length(str)
  allCSi.(str{k}) = [allC.(strOG{k})];
end
allCSi.clusterID = string(allCSi.clusterID);
allCSi.bird = birds';
allCSi = orderfields(allCSi);

cd('C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/')
save allCSi.mat allCSi
disp('saved allCSi')


%% plot 'em

numSi = length(allCSi.clusterID);

params = {'spontRatioVigilance';'spontStateCorr';...
'spontAwakeFR';'spontAwakeBurst';'spontAwakeCV';...
'spontSleepEarlyFR';'spontSleepEarlyBurst';'spontSleepEarlyCV';...
'pbAwakeFR';'pbAwakeBurst';'pbAwakePrecInd10';'pbAwakePrecInd50';'pbAwakeFRPost';'pbAwakeMaxIFR';'pbAwakeChange';'pbAwakePostChange';...
'pbSleepFR';'pbSleepBurst';'pbSleepMaxIFR';'pbSleepPrecInd10';'pbSleepPrecInd50';'pbSleepFRPost';'pbSleepChange';'pbSleepPostChange';...
'wfRise';'wfRise70';'wfFWHM';...
'songSparse';'songPrecInd10';'songPrecInd50';'songFR';'songFRPost';'songBurst';'songMaxIFR';'songChange';'songPostChange';...
'replayFR';'replayBurst';'replayMaxIFR';'nonReplayFR';'nonReplayBurst';'nonReplayMaxIFR';'replayChange'};

% % cd('C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/params/')
% % for i = 41%[1:40, 42:length(params)]
% %    field = string(params(i));
% %    close; 
% %    plotVarSi(allCSi,field)
% %    set(gcf,'PaperOrientation','landscape')
% %    print([field + ' SiOnlyHists'], '-dpdf','-fillpage');
% %    disp(i)
% % end


%% explain 'em

use = rmfield(allCSi,{'allSampWF'; 'meanSampWF';'bird';'clusterID'});
use = cell2mat(struct2cell(use))';
use(isinf(use)) = nan;

[coeff]