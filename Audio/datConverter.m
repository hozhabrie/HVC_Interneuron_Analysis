clear all; close all; clc; warning('off', 'all');

baseFolder = 'V:\Ellie\BOS\eh76\eh76\071322\';
fs = 50000;
fileNames = dir(fullfile(baseFolder, '*.dat'));

if isempty(fileNames)
    error('No .dat files found in the specified directory: %s', baseFolder);
end

for i = 1:length(fileNames)
    fname = fullfile(baseFolder, fileNames(i).name);
    
    % Load audio file
    try
        [audioFileRaw, samplingRate, dateandtime, label, props] = LoadEGUI_daq(fname, 1);
    catch ME
        warning('Error loading file: %s\n%s', fname, ME.message);
        continue; % Skip this file and move to the next
    end
    
    % Normalize audio
    if ~isempty(audioFileRaw)
        audioFileRaw = audioFileRaw ./ max(abs(audioFileRaw));
    else
        warning('Skipping empty audio file: %s', fname);
        continue;
    end

    % Save as WAV file
    outputFilename = fullfile(baseFolder, [fileNames(i).name(1:end-4), '.wav']);
    audiowrite(outputFilename, audioFileRaw, fs);
    
    % Plot spectrogram
    figure;
    vigiSpec(audioFileRaw, samplingRate, [], 0.6);
    title(sprintf('File #%d: %s', i, fileNames(i).name));
    shg;
end

disp('Processing complete.');
