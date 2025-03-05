% This script handles loading abf files from ephys data and extracting the audio file from it for findMotifs


clear all; close all; clc;

% Set parameters
fs = 50000;  % Sampling frequency
baseFolder = 'V:\Ellie\Int Juxta\103118_dlx18\';
cellName = 'c3_0000';

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

% Extract raw audio (channel 3)
rawAudio = fileMat(:,3);

% Create output folder
boutFinderFolder = fullfile(baseFolder, 'boutFinder');
cellFolder = fullfile(boutFinderFolder, cellName);

if ~isfolder(cellFolder)
    mkdir(cellFolder);
end

% Construct output file path
audioFile = fullfile(cellFolder, [cellName, '_playback.wav']);

% Save audio file
audiowrite(audioFile, rawAudio, fs);

disp(['Audio saved: ', audioFile]);
