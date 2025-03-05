clear all; close all; clc; warning('off', 'all');

pause(2); % Allow user to review output
whichTech = input('Enter 1 for Juxta, 2 for Silicon Probe: ');

%% Juxta Processing
if whichTech == 1
    baseFolder = input('Enter the base folder path for Juxta: ', 's');
    
    % Validate the folder
    if ~isfolder(baseFolder)
        error('Invalid directory! Please enter a valid base folder path.');
    end

    cellFile = input('Enter the cell file name (e.g., x5_0): ', 's');
    fullRecordPath = fullfile(baseFolder, cellFile, [cellFile '_filteredspikes.abf']);
    
    % Check if the file exists before proceeding
    if ~isfile(fullRecordPath)
        error('Filtered spikes file not found: %s', fullRecordPath);
    end

    fullRecord = abfload(fullRecordPath);
    playback = fullRecord(:,4);
    fs = 50000;
    thresh = input('Enter threshold (0-1): ');

    BorM = 5; % Default to motif processing
    templateF = fullfile(baseFolder, 'motif.wav');
    outputFolder = fullfile(baseFolder, cellFile);
    audioF = fullfile(outputFolder, [cellFile '_playback.wav']);
    
    % Save playback audio
    audiowrite(audioF, playback, fs);
    
    % Find motifs
    [starts, stops, centers, warps] = findMotifsJuxta(audioF, templateF, thresh);
    
    % Save results
    Motif.file = audioF;
    Motif.start = starts;
    Motif.stop = stops;
    Motif.center = centers;
    Motif.warp = warps;
    Motif.thresh = thresh;
    
    cd(outputFolder);
    save('MotifTimes.mat', 'Motif');

%% Silicon Probe Processing
elseif whichTech == 2   
    baseFolder = input('Enter the base folder path for Silicon Probe: ', 's');

    % Validate the folder
    if ~isfolder(baseFolder)
        error('Invalid directory! Please enter a valid base folder path.');
    end

    addpath('C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\audio\MotifFinder-master\fcns\');
    fs = 30000;
    thresh = 0.7;
    audioF = fullfile(baseFolder, 'audioFull.wav');

    BorM = input('Enter 3 for Singing, 4 for Bout Playback, 5 for Motif Playback: ');

    if BorM == 3
        templateF = fullfile(baseFolder, 'singing', 'template.wav');
        outputFolder = fullfile(baseFolder, 'singing');
    elseif BorM == 5
        templateF = fullfile(baseFolder, 'playback', 'motif.wav');
        outputFolder = fullfile(baseFolder, 'playback', 'motifFolder');
    else
        error('Invalid selection. Only options 3 and 5 are supported.');
    end

    % Find motifs
    [starts, stops, centers, warps] = findMotifsSi(audioF, templateF, thresh);

    % Save results
    Motif.file = audioF;
    Motif.start = starts;
    Motif.stop = stops;
    Motif.center = centers;
    Motif.warp = warps;
    Motif.thresh = thresh;

    cd(outputFolder);
    save('MotifTimes.mat', 'Motif');

    % Display start/stop times
    disp([Motif.start; Motif.stop]');

else
    error('Invalid selection. Enter 1 for Juxta or 2 for Silicon Probe.');
end
