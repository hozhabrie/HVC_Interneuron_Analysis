clc; warning off; clear all; close all;

% Define bird and data path
bird = 'C50';
path = ['C:\Users\User\Dropbox (NYU Langone Health)\SiProbe/' bird '/'];

% Load cluster and recording data
load([path 'data/clusters.mat']);
fs = 30000;
numUnits = length(clusters);

% Check if it's the first pass
pass = input('First pass? 1/0: ');
if pass
    toPrint = ones(length(clusters),1);
    startRec = zeros(length(clusters),1);
    stopRec = string(repmat('end',length(clusters),1));
    mkdir([path 'amps/']);
else
    goodCells = readtable([path bird '_goodSiUnits.xlsx']);
    toPrint = (goodCells{:,3});
    startRec = (goodCells{:,4});
    stopRec = (goodCells{:,5});
end

% Load spike cluster, amplitude, and time data
spkClus = readNPY([path 'spike_clusters.npy']);
amps = readNPY([path 'amplitudes.npy']);
spkTimes = readNPY([path 'spike_times.npy']);
load([path 'data/lenRec.mat']);
lenRec = lenRec / fs;
fontsize = 15;

% Load light signal and audio labels
load([path 'data/lightSig.mat']);
audioPath = [path 'audio/'];
load([audioPath 'audioLabels.mat']);
template = audioread([audioPath 'singing/template.wav']);
motif = audioread([audioPath 'playback/motif.wav']);

labels = {'sm'; 'pm'; 'pms'};
for i = 1:length(labels)
    if strcmp(labels{i}, 'sm')
        audio = template;
    else
        audio = motif;
    end
    Event(i) = getEventTimes(audioPath, fs, audio, audioLabels, labels{i});
end

%% Plot Amplitudes
for n = 45 % Specific unit index for visualization
    unit = clusters(n).clusterID;
    close all;
    if toPrint(n) == 1
        figure('Visible','off','WindowState','maximized');

        Tstart = startRec(n);
        if strcmp(stopRec{n}, 'end')
            Tstop = lenRec;
        else
            Tstop = str2num(stopRec{n});
        end

        % Extract and filter spike data
        goodSpks = spkTimes(spkClus == unit) / fs;
        goods = goodSpks > Tstart & goodSpks < Tstop;
        goodSpks = goodSpks(goods);
        goodAmps = amps(spkClus == unit);
        goodAmps = goodAmps(goods);

        % Plot spike amplitudes
        hold on;
        plot(goodSpks(1:10:end), goodAmps(1:10:end), 'k.');
        title(['Unit ' num2str(unit)]);
        xlim([0 lenRec+100]);
        ylabel('MicroVolts');
        xlabel('Seconds');
        set(gca,'color','none','box','off','TickDir','out','FontSize', fontsize);
        ylim([0 quantile(goodAmps, .999) + 10]);

        % Overlay event times
        for j = 1:length(Event)
            for i = 1:length(Event(j).start)
                patch([Event(j).start(i) Event(j).stop(i) Event(j).stop(i) Event(j).start(i)], 
                      [max(ylim)*.9 max(ylim)*.9 max(ylim)*.8 max(ylim)*.8], 
                      'k', 'EdgeColor', 'k', 'EdgeAlpha', .2, 'FaceAlpha', .1);
                hold on;
            end
        end

        % Set figure properties
        set(gcf, 'PaperOrientation', 'landscape', 'PaperUnits', 'inches', 'PaperSize', [8 11]);
        disp(['Plotted amplitude data for unit ' num2str(unit)]);
    end
end