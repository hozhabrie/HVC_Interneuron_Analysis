clear all; close all; warning('off', 'all'); clc;

fs = 30000;
wing = 1;
wind = [-wing wing]; 
labels = {'sm', 'pm', 'pms'};
samplingInterval = 1/fs;
wavWind = [-fs/1000 fs/1000];
timeAxis = linspace(0, (2 * wavWind(2) + 1) / fs * 1000, (2 * wavWind(2) + 1));
ticklen = 3;
spontWind = [0 5];

% Hide figures by default
set(0, 'DefaultFigureVisible', 'off');

birds = {'C22', 'C24', 'C50', 'C52', 'C57'};
dir = 'C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/';

for w = 1:length(birds)
    bird = birds{w};
    path = fullfile(dir, bird);
    audioPath = fullfile(path, 'audio', 'warp', 'warpAudioLabels.mat');
    
    % Check if file exists before loading
    if ~isfile(audioPath)
        warning('File not found: %s', audioPath);
        continue; % Skip to next bird if file is missing
    end
    
    load(audioPath, 'warpAudioLabels');
    
    if ~exist('warpAudioLabels', 'var') || isempty(warpAudioLabels)
        warning('warpAudioLabels variable is missing or empty in: %s', audioPath);
        continue;
    end
    
    times = diff([warpAudioLabels.stop warpAudioLabels.start]');

    % Ensure figure visibility inside loop
    figure('Visible', 'on');
    plot(times);
    title(['Time Differences for Bird: ', bird]);
end
