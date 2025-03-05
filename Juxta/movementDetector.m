%% load shit
clear all; close all; clc; warning off;
baseFolder = 'V:\Ellie\Int Juxta\091218_dlx15\';
subFolder = 'boutFinder\';
f = 'c9_0';
load([baseFolder subFolder f '\goodDataTimes.mat']);
load([baseFolder subFolder f '\threshY.mat']);
load([baseFolder subFolder f '\artifactThresh.mat']);
load([baseFolder subFolder f '\windBack.mat']);
load([baseFolder subFolder f '\windFor.mat']);
spikedata = abfload([baseFolder subFolder f '\' f '_filteredspikes.abf']);
fs = 50000;
smoothWindow = .002 * fs;
% threshY = 2e-2;
% artifactThresh = -5;
% windBack = .5 * fs;
% windFor = .5 * fs;
% goodDataTimes = 1:(611*fs);

spikes = spikedata(:,1);
if goodDataTimes(1) ~= 1
    spikes(1:goodDataTimes(1)) = nan;
elseif goodDataTimes(end) < length(spikes) 
    spikes(goodDataTimes(end):end) = nan;
end

clear spikedata;    
disp([goodDataTimes(1)/fs, goodDataTimes(end)/fs])

%% extract mvmt envelopes
mvmtTime = 1:length(spikes);
mvmtSpks = spikes(mvmtTime);
xAxis = (linspace(0,length(mvmtSpks)/fs, length(mvmtSpks)))';
smoothed = smoothdata(mvmtSpks, 'gaussian', smoothWindow);
smoothed = smoothed.^2;
smoothed = normalize(smoothed, 'range');

%calculate mvmt env timepoints
env = movmean(smoothed, [windBack windFor]);
[locs, ~] = find(env >= threshY);

%if ~isempty(locs)
startMove = [locs(1)];
stopMove = [];
for i = 1:(length(locs)-1)
    if locs(i+1)-locs(i) ~= 1
        stopMove = [stopMove; locs(i)];
        startMove = [startMove; locs(i+1)];
    end
end
stopMove = [stopMove; locs(end)];

%concatenate close together mvmt bouts
for i = 2:length(startMove)
    if startMove(i)-stopMove(i-1) < .1 * fs
        startMove(i) = nan;
        stopMove(i-1) = nan;
    end
end
startMove = startMove(~isnan(startMove));
stopMove = stopMove(~isnan(stopMove));

%remove short mvmt timestamps
for i = 1:length(startMove)
    if stopMove(i)-startMove(i) < .1 * fs
        startMove(i) = nan;
        stopMove(i) = nan;
    end
end
startMove = startMove(~isnan(startMove));
stopMove = stopMove(~isnan(stopMove));

%% remove mvmt bouts and blips from spikedata
fixedSpikes = spikes;
for i = 1:length(startMove)
    for j = startMove(i):stopMove(i)
        fixedSpikes(j) = 0;
    end
end
%remove other artifacts
[artlocs, ~] = find(fixedSpikes < artifactThresh);
for i = 1:length(artlocs)
    artifact = [artlocs(i)-smoothWindow:artlocs(i)+smoothWindow];
    fixedSpikes(artifact) = 0;
end

[artlocs, ~] = find(fixedSpikes > -artifactThresh);
for i = 1:length(artlocs)
    artifact = [artlocs(i)-smoothWindow:artlocs(i)+smoothWindow];
    fixedSpikes(artifact) = 0;
end

%% plot it
figure; 
ax(1) = subplot(4,1,1);
plot(xAxis, mvmtSpks, 'Color', 'k'); 
hold on;
plot(startMove/fs, threshY, 'go')
hold on;
plot(stopMove/fs, threshY, 'ro')
xlim([xAxis(1) xAxis(end)]);
ylim([-1 1.2]);
ylabel('mV');
set(gca,'TickDir','out')    
set(gca, 'xtick', [])
set(gca, 'box', 'off')

ax(2) = subplot(4,1,2);
plot(xAxis, smoothed, 'Color', 'k');
hold on;
plot(startMove/fs, threshY, 'go')
hold on;
plot(stopMove/fs, threshY, 'ro')
xlim([xAxis(1) xAxis(end)]);
set(gca,'TickDir','out')   
set(gca, 'ytick', [])
set(gca, 'xtick', [])
set(gca, 'box', 'off')
ylim([0 threshY*20])

ax(3) = subplot(4,1,3);
plot(xAxis, env, 'Color', 'k')
xlim([xAxis(1) xAxis(end)]);
hold on;
plot(startMove/fs, threshY, 'go')
hold on;
plot(stopMove/fs, threshY, 'ro')
hold on;
line([xAxis(1) xAxis(end)],[threshY threshY])
set(gca, 'ytick', [])
set(gca, 'xtick', [])
set(gca,'TickDir','out')
set(gca, 'box', 'off')


ax(4) = subplot(4,1,4);
plot(xAxis, fixedSpikes, 'Color', 'k');
set(gca,'TickDir','out')
set(gca, 'box', 'off')
xlabel('Time (s)')
ylim([-1 1.2]);
ylabel('mV');
linkaxes(ax, 'x');

% figure; 
% %ax2(1) = subplot(2,1,1);
% plot(xAxis, spikes, 'Color', 'k');
% %ax2(2) = subplot(2,1,2);
% %plot(xAxis, fixedSpikes, 'Color', 'k'); 
% ylim([-1.5 1.5])
% xlim([291 296])
% xlabel('Time (s)')
% ylabel('mV')
% set(gca,'TickDir','out')  
% set(gca,'box','off')
% %linkaxes(ax2, 'x');

% %% plot movement modulation
% thresh = .2;
% window = .002 * fs; 
% wing = 10 * fs;
% bin = 10;
% move = startMove;
% temp = fixedSpikes; 
% temp(temp==0) = nan;
% [allVals, allLocs] = findpeaks(-1*temp, 'MinPeakHeight', thresh, 'MinPeakDistance', window); 
% for i = 1:length(allLocs)
%      tempwind = allLocs(i):allLocs(i)+window;
%      [topVals, topLocs] = find(fixedSpikes(tempwind) > .15);
%      if isempty(topLocs)
%         allLocs(i) = nan;
%         allVals(i) = nan;
%      end
% end
% allLocs = allLocs(~isnan(allLocs));
% rasterAll = zeros(length(temp), 1);
% rasterAll(allLocs) = 1;
% rasterAll(isnan(temp)) = nan;
% binRast = zeros(round(length(rasterAll)/fs/bin), 1)';
% count = 1;
% for i = 1:bin*fs:length(rasterAll)-bin*fs
%     binTimes = i:i+bin*fs;
%     binSpikes = rasterAll(binTimes);
%     sumNans = sum(isnan(binSpikes));
%     binRast(count) = nansum(binSpikes)/(bin-sumNans/fs);
%     count = count + 1;
% end
% 
% figure;
% plot(linspace(0, length(rasterAll)/fs, length(binRast)), binRast, 'Color', 'k');
% hold on
% for i = 1:length(startMove)
%     line([startMove(i)/fs, stopMove(i)/fs], [.7 .7], 'LineWidth', 15, 'Color', 'r')
%     hold on 
% end
% set(gca,'TickDir','out')  
% set(gca,'box','off')
% ylabel('Avg Firing Rate (Hz)')
% xlabel('Time (s)')
% title('Modulation of firing rate in relation to movement')
% % figure;
% % for i = 1:length(move)
% %    if move(i)-wing > 1 & move(i)+wing < length(fixedSpikes)
% %        temp = fixedSpikes(move(i)-wing:move(i)+wing);
% %        raster = zeros(length(temp), 1);
% %        [allVals, allLocs] = findpeaks(-1*temp, 'MinPeakHeight', thresh, 'MinPeakDistance', window); 
% %        for i = 1:length(allLocs)
% %             tempwind = allLocs(i):allLocs(i)+window;
% %             [topVals, topLocs] = find(fixedSpikes(tempwind) > .1);
% %             if isempty(topLocs)
% %                 allLocs(i) = nan;
% %                 allVals(i) = nan;
% %             end
% %        end
% %        allLocs = allLocs(~isnan(allLocs));
% %        allVals = allVals(~isnan(allVals));
% %        raster(allLocs) = 1;
% %        for k = 1:length(raster)
% %           if raster(k) == 1 
% %               plot(([k, k])/fs, [0, 1]-1*i,'Color', 'k')
% %               hold on
% %           end
% %        end
% %    end
% %    hold on
% % % end

%else
%fixedSpikes = spikes;
%end
%% save shit
cd([baseFolder subFolder f])
% pause(1)
% save fixedSpikes.mat fixedSpikes
% save threshY.mat threshY
% save artifactThresh.mat artifactThresh
% save goodDataTimes.mat goodDataTimes
data = fixedSpikes;
save clustSpikes2.mat data