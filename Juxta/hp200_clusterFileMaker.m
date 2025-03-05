clear all; close all; clc; warning('off', 'all');

% Set parameters
fs = 50000;   % Sampling frequency
hp = 200;     % High-pass filter cutoff frequency
baseFolder = 'C:\Users\User\Dropbox (NYU Langone Health)\Int Juxta\032819_dlx32\';
cellName = 'c1_0001';

% Construct file path for ABF file
abfFile = fullfile(baseFolder, [cellName, '.abf']);

% Check if the file exists before attempting to load
if ~isfile(abfFile)
    error('ABF file not found: %s', abfFile);
end

% Load ABF file
fileMat = abfload(abfFile);

% Process cell name to remove extra zeros
cellName(regexp(cellName, '0')) = [];
if cellName(end) == '_'
    cellName = [cellName, '0'];
end

% Define output folder paths
newFolder = fullfile(baseFolder, 'hp200');
cellFolder = fullfile(newFolder, cellName);

% Create necessary directories if they do not exist
if ~isfolder(newFolder)
    mkdir(newFolder);
end
if ~isfolder(cellFolder)
    mkdir(cellFolder);
end

% Extract spike data (channel 1) and apply high-pass filter
rawSpikes = fileMat(:,1);
filteredData = highpass(rawSpikes, hp, fs);

% Construct output file path
hpFile = fullfile(cellFolder, [cellName, '_hp200.mat']);

% Save filtered spikes
save(hpFile, 'filteredData');

disp(['Filtered spikes saved: ', hpFile]);

% Uncomment to visualize raw vs. filtered spikes
% figure;
% ax(1) = subplot(2,1,1);
% plot(rawSpikes);
% title('Raw Spikes');
% ax(2) = subplot(2,1,2);
% plot(filteredData);
% title('High-Pass Filtered Spikes');
% linkaxes(ax, 'x');

% Run wave_clus for clustering (ensure it's installed and in MATLAB path)
wave_clus
