close all; clear all; clc; warning('off', 'all');

% Define base directory for Juxta data
baseDir = 'C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/';

% Load goodCells information
goodCellsFile = fullfile(baseDir, 'goodCellsnew.xlsx');
if ~isfile(goodCellsFile)
    error('File not found: %s', goodCellsFile);
end

goodCells = readtable(goodCellsFile);
birds = string(goodCells{:,1});
cells = string(goodCells{:,2});
toPrint = logical(goodCells{:,5}); % Logical array for cells to process

% Set parameters
ds = 50;  % Downsampling factor
fs = 50000; % Sampling frequency

for number = 1:height(goodCells)
    if toPrint(number)  % Process only marked cells
        cellName = char(cells(number));
        birdFolder = fullfile(baseDir, birds(number));
        cellFolder = fullfile(birdFolder, cellName);

        % Construct file paths
        audioFile = fullfile(cellFolder, [cellName '_playback.wav']);
        motifFile = fullfile(cellFolder, 'MotifTimes.mat');

        % Display processing info
        disp(['Processing: ', birds(number), ' - ', cellName]);

        % Check if required files exist
        if ~isfile(audioFile)
            warning('Audio file not found: %s', audioFile);
            continue;
        end
        if ~isfile(motifFile)
            warning('MotifTimes.mat not found: %s', motifFile);
            continue;
        end

        % Load audio and motif data
        audio = audioread(audioFile);
        audioDs = downsample(audio, ds);
        load(motifFile, 'Motif');

        % Plot full downsampled waveform with motif regions
        figure;
        plot(audioDs * fs);
        hold on;
        for i = 1:length(Motif.start)
            patch([Motif.start(i) * fs/ds, Motif.stop(i) * fs/ds, Motif.stop(i) * fs/ds, Motif.start(i) * fs/ds], ...
                  [-1000 -1000 1000 1000], 'red', 'FaceAlpha', .3);
            hold on;
        end
        title(['Motif Regions: ', birds(number), ' - ', cellName]);
        xlabel('Time (samples)');
        ylabel('Amplitude');

        % Plot individual motifs
        figure;
        for i = 1:length(Motif.start)
            plot(audio(Motif.start(i)*fs : Motif.stop(i)*fs) - i);
            hold on;
        end
        title(['Individual Motifs: ', birds(number), ' - ', cellName]);
        xlabel('Time (samples)');
        ylabel('Amplitude');
    end
end

