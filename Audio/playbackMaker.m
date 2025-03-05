%% Playback Maker: Short or Long?
clear all; clc; close all;
warning('off', 'all');  % Suppress all warnings

% Prompt user for mode selection
modeChoice = input('Enter 1 for Short/Juxta Playback, 2 for Long/Si Probe Playback: ');

% Prompt user for base folder
baseFolder = input('Enter the base folder path: ', 's');

% Validate the base folder path
if ~isfolder(baseFolder)
    error('Invalid directory! Please enter a valid base folder path.');
end

% Define Sampling Rate Based on Mode
if modeChoice == 1
    fs = 50000;
    playbackType = 'Short/Juxta';
elseif modeChoice == 2
    fs = 30000;  % Ensure correct sampling rate
    playbackType = 'Long/Si Probe';
else
    error('Invalid choice! Please enter 1 (Short/Juxta) or 2 (Long/Si Probe).');
end

% Load motif file
motifPath = fullfile(baseFolder, 'motif.wav');
if ~isfile(motifPath)
    error('Motif file not found: %s', motifPath);
end
[motifFile, ~] = audioread(motifPath);

% Spectrogram Visualization
figure;
vigiSpec(motifFile, fs, [], 0.65);
title([playbackType, ' Playback Motif']);
shg;

% Playback Parameters
numMotifs = 50;
itiMotifStandard = zeros(15 * fs, 1);
itiJitter = randi([fs, 5 * fs], 1, numMotifs);  % Random jitter between 1-5 sec

% Generate Playback File Based on Mode
if modeChoice == 1
    % Short Playback (Juxta)
    playbackFile = zeros(15 * fs, 1);
    for i = 1:numMotifs
        playbackFile = [playbackFile; motifFile; itiMotifStandard; zeros(itiJitter(i), 1)];
    end
    wavName = fullfile(baseFolder, 'playbackFile.wav');
else
    % Long Playback (Si Probe)
    playbackFile = zeros(10800 * fs, 1);  % 3-hour silence at start
    for i = 1:2  % Two playback sessions
        for j = 1:numMotifs
            playbackFile = [playbackFile; motifFile; itiMotifStandard; zeros(itiJitter(j), 1)];
        end
        playbackFile = [playbackFile; zeros(7200 * fs, 1)];  % 2-hour silence
    end
    wavName = fullfile(baseFolder, 'playbackFileLong.wav');
end

% Save Playback File
audiowrite(wavName, playbackFile, fs);
disp(['Playback file saved: ', wavName]);

% Plot Playback Waveform
figure;
timeAxis = linspace(0, length(playbackFile) / fs / 60, length(playbackFile));
plot(timeAxis, playbackFile);
xlabel('Time (minutes)');
ylabel('Amplitude');
title([playbackType, ' Playback File']);
shg;
