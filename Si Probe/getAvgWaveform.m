%% Set up environment

clear all; close all; warning off; clc;

% Define list of birds and recording window
birds = {'C22'; 'C24'; 'C50'; 'C52'; 'C57'};
start = 36;
stop = 56;

% Loop through selected birds
for w = 1:length(birds)
    bird = char(birds(w));
    disp(['Processing bird: ', bird]);
    
    % Define file paths
    path = ['E:/SiProbe/'  bird '/'];
    fs = 30000;
    
    % Load necessary data files
    disp('Loading data...');
    load([path 'KS_Output.mat']);
    cd('C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Extra\chanMaps\');
    load(ops.chanMap, 'AnatGroup', 'connected', 'xcoords', 'ycoords');
    clear rez;
    load([path 'data/clusters.mat']);
    disp('Data loaded.');
    
    % Load good units from Excel
    opts = detectImportOptions([path bird '_goodSiUnits.xlsx']);
    opts.DataRange = 'A2';
    opts.VariableNamesRange = 'A1';
    goodCells = readtable([path bird '_goodSiUnits.xlsx'], opts);
    toPrint = (goodCells{:,3});
    startRec = (goodCells{:,4});
    stopRec = (goodCells{:,5});
    
    %% Extract mean waveforms
    
    samplingInterval = 1 / fs;
    wind = [-1.5 * fs / 1000, 1.5 * fs / 1000];
    numUnits = length(clusters);
    waveforms = nan(121, numUnits);
    timeAxis = linspace(0, (2 * wind(2) + 1) / fs * 1000, (2 * wind(2) + 1));
    
    for n = 1:numUnits
        if toPrint(n) == 1
            % Select the maximum channel
            channel = clusters(n).maxChannelWF + 1;
            x = xcoords(channel);
            y = ycoords(channel);
            
            % Sample a subset of spike times for waveform extraction
            samples = randperm(length(clusters(n).spikeTimes), min(length(clusters(n).spikeTimes), 2000));
            sampleTimes = round(clusters(n).spikeTimes(samples) * fs);
            
            % Read waveforms from .dat file
            [meanWF, allWF] = readWaveformsFromDat(['E:/SiProbe/' bird '/amplifier.dat'], ops.NchanTOT, channel, sampleTimes, wind, []);
            
            % Filter waveforms
            [b, a] = butter(3, 500 / (fs / 2), 'high');
            allWF2 = filtfilt(b, a, allWF);
            clear allWF;
            
            % Align waveforms by local minima
            numWavs = size(allWF2, 2);
            fixedWF = zeros(61, numWavs);
            center = find(islocalmin(meanWF, 'MaxNumExtrema', 1));
            
            for i = 1:numWavs
                eachCenter = start + find(islocalmin(allWF2(start:stop, i), 'MaxNumExtrema', 1));
                jitter = eachCenter - center;
                if isempty(jitter)
                    jitter = 0;
                end
                fixedWF(:, i) = allWF2(16 + jitter:76 + jitter, i);
            end
            
            % Compute mean waveform
            meanWF2 = mean(fixedWF, 2);
            clusters(n).meanWF = meanWF2;
            
            disp(['Processed unit: ', num2str(n)]);
        end
    end
    
    % Save processed clusters
    cd([path 'data/']);
    save clusters.mat clusters;
    disp(['Saved updated clusters for ', bird]);
end
