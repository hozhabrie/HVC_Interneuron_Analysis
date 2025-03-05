close all; clear all; clc; warning('off', 'all');

% Define bird IDs and sampling parameters
birds = {'C22', 'C24', 'C50', 'C52', 'C57'};
ds = 300; 
fs = 30000;

for w = 1:length(birds)
    whichBird = string(birds{w});
    path = fullfile('E:\SiProbe\', whichBird);
    smTimesFile = fullfile(path, 'data', 'smTimesAll.mat');
    audioFile = fullfile(path, 'audio', 'audioFull.wav');

    % Check if required files exist before proceeding
    if ~isfile(smTimesFile)
        warning('File not found: %s', smTimesFile);
        continue;
    end
    if ~isfile(audioFile)
        warning('File not found: %s', audioFile);
        continue;
    end

    % Load motif timing data
    load(smTimesFile, 'smTimesAll');
    structElem = fieldnames(smTimesAll);

    % Load audio file
    audio = audioread(audioFile);
    
    % Plot motifs for the bird
    figure;
    for i = 1:length(smTimesAll(5).start)
        motifRange = round(smTimesAll(5).start(i) * fs : smTimesAll(5).stop(i) * fs);
        plot(audio(motifRange) - i);
        set(gca, 'Color', 'none', 'Box', 'off', 'TickDir', 'out');
        title(whichBird);
        ylim([-length(smTimesAll(5).start) - 1, 0]);
        hold on;
    end

    % Save figure as a PDF
    outputDir = 'C:\Users\User\Dropbox (NYU Langone Health)\SiProbe\song motifs\';
    if ~isfolder(outputDir)
        warning('Output directory does not exist: %s', outputDir);
    else
        print(fullfile(outputDir, whichBird + " song motifs fixed"), '-dpdf', '-fillpage');
    end
    
    shg; % Show figure
end
