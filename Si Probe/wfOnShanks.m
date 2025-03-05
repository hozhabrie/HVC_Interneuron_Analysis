%% Setup
F = "C:\Users\User\Dropbox (NYU Langone Health)\SiProbe\C50\";
G = "E:\SiProbe\C50\";

% Load necessary data
load(fullfile(F, "data", "clusters.mat"));
addpath("C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Extra\chanMaps\");
load(fullfile(F, "KS_Output.mat"));
load(ops.chanMap, 'AnatGroup', 'connected', 'xcoords', 'ycoords');

% Constants
SiProbeSamplingInterval = 1 / 3e4;
downSampleFactor = 100;
timeAxis = [];
meanWF = [];

% Process selected clusters
selectedClusters = [55]; % Modify as needed (e.g., 34-2, 20-3, 45-4)
for n = selectedClusters
    shank = clusters(n).shank;
    sampleTimes = round(clusters(n).spikeTimes(1:downSampleFactor:end) / SiProbeSamplingInterval);
    sampleWindow = [-2 / SiProbeSamplingInterval / 1000, 2 / SiProbeSamplingInterval / 1000];
    deltaWindow = sampleWindow(2) - sampleWindow(1);
    medianEndBin = floor(deltaWindow / 3);

    % Find corresponding channel
    lookUpChannels = [AnatGroup{:}];
    lookUpChannels = lookUpChannels(connected);
    channel = lookUpChannels(clusters(n).maxChannelTemplate + 1);

    % Read waveforms
    [tmpWF, ~] = readWaveformsFromDat(fullfile(G, 'amplifier.dat'), 128, 1:128, sampleTimes, sampleWindow, []);
    tmpWF = tmpWF - median(tmpWF(1:medianEndBin));

    % Store waveform and time axis
    meanWF = [meanWF; tmpWF'];
    timeAxis = [timeAxis; (1:length(tmpWF)) * 1000 * SiProbeSamplingInterval];
end

%% Plot Waveforms
close all;
figure; hold on;

for cH = 33:64
    waveForm = meanWF(:, cH);
    x = xcoords(cH);
    y = ycoords(cH);
    plot(((1:length(waveForm)) / 30) + x / 2, (waveForm + y * 5), 'LineWidth', 2, 'Color', 'k');
end

title(['Cluster No. ', num2str(clusters(n).clusterID)]);
set(gca, 'TickDir', 'out', 'Box', 'off');

%% Autocorrelation Analysis
close all;
selectedClusters = [55, 34, 20, 45];
count = 0;

for n = selectedClusters
    count = count + 1;
    spks = clusters(n).spikeTimes;
    isi = [diff(-spks); diff(spks)];

    figure(count);
    histogram(isi, -0.02:0.0001:0.02, 'Normalization', 'probability', 'EdgeColor', 'none');
    title(clusters(n).clusterID);
    set(gca, 'TickDir', 'out', 'Color', 'none', 'Box', 'off', 'FontSize', 15);

    % Save figure
    saveDir = "C:\Users\User\Dropbox (NYU Langone Health)\Thesis\figures\";
    print(fullfile(saveDir, string(clusters(n).clusterID) + "_autocorr_zoom"), '-dpdf');
end
