clear all; close all; clc; warning off;

fs = 30000;
birds = {'C22'; 'C24'; 'C50'; 'C52'};
addpath 'C:/Users/User/Dropbox (NYU Langone Health)/EllieScripts/SiClustering/Fcns/';
w = 4;
bird = char(birds(w));
disp(bird)
path = ['C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/' bird '/'];
load([path 'data/clusters.mat']);
audio = audioread([path 'audio/audioFull.wav']);
goodCells = readtable([path bird '_goodSiUnits.xlsx']);
toPrint = (goodCells{:,3});
%turnTakeTime = [2680 2710];
turnTakeTime.start = 5249;
turnTakeTime.stop = 5255;
turnTakeTime.warp = 1;
audio = audio(turnTakeTime.start*fs:turnTakeTime.stop*fs);


%%
[ttSpks, ~] = alignSpikeTimesSi(clusters, [0 6], turnTakeTime, turnTakeTime.start);

figure;

ax(1) = subplot(5,1,1);
plot(linspace(0, 6, length(audio)),audio,'k')

hold on
set(gca,'box','off','color','none','xcolor','none','ycolor','none')

ax(2) = subplot(5,1,2:5);
count = 0;
for i = 1:length(clusters)
        if toPrint(i) == 1
            count = count + 1;
            id(count) = string(clusters(i).clusterID); 
            for k = 1:length(turnTakeTime.start)
                trial = ttSpks(i).eventTrials == k;
                trialSpks = ttSpks(i).eventTimesAligned(trial);
                for i = 1:length(trialSpks)
                    plot(trialSpks,zeros(length(trialSpks),1)+count, 'k|',"MarkerSize",1.2)
                    hold on
                end
                hold on
            end 
        end
 end
%xlim([0 6])
ylim([0 count+1])
yticks(1:count)
yticklabels(id)
set(gca,'box','off','color','none','TickDir','out', 'xTick',[])
set(gcf,'PaperOrientation','landscape')
linkaxes(ax,'x')

cd([path 'unitsModulation']);
print(['turntaking'], '-dpdf','-bestfit');
