clc; warning off;
%set(0,'DefaultFigureVisible','off');

fs = 30000;
fontsize = 15;

% Save data
allbirds = {'C22'; 'C24'; 'C50'; 'C52'; 'C57'};
allGoodCells = [];
allC = [];
dir = '/Users/ellie/NYU Langone Health Dropbox/Elnaz Hozhabri/SiProbe/';

for w = 1:length(allbirds)
    clear clusters;
    bird = char(allbirds(w));
    path = [dir bird '/'];

    opts = detectImportOptions([path bird '_goodSiUnits.xlsx']);
    opts.DataRange = 'A2';
    opts.VariableNamesRange = 'A1';
    goodCells = readtable([path bird '_goodSiUnits.xlsx'], opts);

    take = logical(goodCells{:,3});
    allGoodCells = [allGoodCells; goodCells(take,[1:5,6])];
end

birds = string([allGoodCells{:,1}]);
birds = erase(birds,'eh');

load([dir 'Connectivity/kmID.mat']);
identity = [kmID{:,3}]';
ids = string([kmID{:,2}])';
raksinVars = {'pbAwakeFRAll','pbSleepFRAll','prePbAwakeFRAll','prePbSleepFRAll','songFRAll'};

for n = 1:length(allC)
    for m = 1:length(raksinVars)
        raksin(n).(raksinVars{m}) = allC(n).(raksinVars{m});
    end
    allC(n).clusterID = ids(n);
    allC(n).identity = identity(n);
    allC(n).bird = birds(n);
end

allC = orderfields(allC);
disp('saved allC');

% Remove unnecessary fields
shitfields = {'spikeTimes','coordinates','shank','maxChannelTemplate','maxChannelWF',
    'ogTimes','spontAwakeIsis','smIsis','pmIsis','pmsIsis',
    'spontSleepDeepIsis','spontSleepEarlyIsis','pbAwakeFRAll','pbSleepFRAll',
    'songFRAll','prePbAwakeFRAll','prePbSleepFRAll','meanWF','wfNorm'};
allC = rmfield(allC, shitfields);

strOG = fieldnames(allC);
str = strOG;
for p = 1:length(str)
    allCSi.(str{p}) = [allC.(strOG{p})];
end

fullLabel = [string(identity) + '-' + birds + '-' + ids];
allCSi.identity = fullLabel;
allCSi = orderfields(allCSi);

which = identity == 'i';
ints = fullLabel(which);
disp('saved allCSi');

%good = which' & allCSi.bird ~= 'eh50';
good = which;
raksin = raksin(good);

%% Plot parameters
params = {'spontRatioVigilanceE','spontRatioVigilanceD','spontStateCorrE','spontStateCorrD',
    'spontAwakeFR','spontAwakeBurst','spontAwakeCV','spontSongRate','spontSongBurst','spontSongCV',
    'spontSleepEarlyFR','spontSleepEarlyBurst','spontSleepEarlyCV','spontSleepDeepFR','spontSleepDeepBurst',
    'spontSleepDeepCV','pbAwakeFR','pbAwakeBurst','pbAwakePrecInd10','pbAwakePrecInd50',
    'pbAwakeFRPost','pbAwakeMaxIFR','pbAwakeChange','pbAwakePostChange','pbSleepFR','pbSleepBurst',
    'pbSleepMaxIFR','pbSleepPrecInd10','pbSleepPrecInd50','pbSleepFRPost','pbSleepChange',
    'pbSleepPostChange','wfRise','wfRise70','wfFWHM','mult70','songSparse','songPrecInd10',
    'songFR','songFRPost','songBurst','songChange','songPostChange','songPrecInd50','songMaxIFR',
    'songCorr','songRtoZ','songMatch','replayFR','replayBurst','replayMaxIFR','nonReplayFR',
    'nonReplayBurst','nonReplayMaxIFR','replayChange','pbAwakeMatch','pbSleepMatch'};

%% Normalize data
addpath('/Users/ellie/NYU Langone Health Dropbox/Elnaz Hozhabri/EllieScripts/');

useOG = rmfield(allCSi,{'bird','clusterID','identity'});

goodVars = {'wfRise','wfFWHM','spontAwakeFR','spontAwakeCV','spontSleepEarlyFR','spontSleepEarlyCV',
    'spontSleepDeepFR','spontSleepDeepCV','spontStateCorrE','spontRatioVigilanceE','spontStateCorrD',
    'spontRatioVigilanceD','pbAwakeFR','pbAwakeChange','pbAwakeCorr','pbSleepFR','pbSleepChange',
    'pbSleepCorr','songFR','songSparse','songChange','songBurst','songCorr'}';

logs = [3,5,7,13,16,19];
atanhs = [9,11,15,18,23];
abss = [14,17];

useGood = [];
for i = 1:length(goodVars)
    pars = getfield(useOG, goodVars{i});
    pars = pars(good);
    if ismember(i, logs)
        pars = log(pars);
    elseif ismember(i, atanhs)
        pars = matchRatio(pars);
    elseif ismember(i, abss)
        pars = abs(pars);
    end
    pars(isinf(pars)) = nan;
    useGood(:,i) = pars;
end
useGood = zscorenan(useGood);

useGood = useGood(:,[1:4,7:8,11:end]);
goodVars = goodVars([1:4,7:8,11:end]);

%% Correlation Matrix
[R, P] = corrcoef(useGood, 'rows', 'complete');
figure; imagesc(R); hold on;
for a = 1:length(goodVars)
    text(1:length(goodVars), a * ones(length(goodVars), 1), string(P(a,:) < .05) * '*', 'Color', 'k');
end
xticks(1:length(goodVars));
yticks(1:length(goodVars));
set(gca, 'TickDir', 'out', 'color', 'none', 'box', 'off', 'FontSize', 9);
colorbar;
set(gcf, 'PaperSize', [11 10]);

%% Cluster 3D Plots
symbols = 'svd*o.|*';
colors = 'gbrmkcgb';
state = matchRatio(allCSi.spontStateCorrE(good));
modul = allCSi.pbAwakeChange(good);

figure;
hold on;
for i = 2
    zuse = [zscorenan(song{i}); zscorenan(modul); state]';
    clus = kmeans(zuse, 8, 'Replicates', 500);
    for j = 1:8
        plot3(zuse(clus == j, 1), zuse(clus == j, 2), zuse(clus == j, 3), [colors(j) symbols(j)]);
    end
    title(['Song Param Clusters, K = ' num2str(8)]);
end
xlabel('Song'); ylabel('Playback'); zlabel('State');
set(gca, 'box', 'off', 'TickDir', 'out');
