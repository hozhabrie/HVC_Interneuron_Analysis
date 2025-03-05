clc; warning off;

% Add necessary script functions directory to path
addpath 'C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Fcns\'

% Define sampling frequency
fs = 30000;

% Define analysis window range
wing = 1;
wind = [-wing wing];

% Define labels for different event types
labels = {'sm'; 'pm'; 'pms'};

% Set sampling interval
samplingInterval = 1/fs;

% Define waveform window for alignment
wavWind = [-1*fs/1000 1*fs/1000];
timeAxis = linspace(0, (2*wavWind(2)+1)/fs*1000, (2*wavWind(2)+1));

% Define tick length for plots
ticklen = 3;

% Define spontaneous activity analysis window
spontWind = [0 5];

% Disable figure visibility
set(0,'DefaultFigureVisible','off');

% List of bird identifiers
birds = {'C22'; 'C24'; 'C50'; 'C52'; 'C57'};

dir = 'C:\Users\User\Dropbox (NYU Langone Health)\SiProbe/';

% Loop through birds for analysis
for w = length(birds)
    bird = char(birds(w));
    path = [dir bird '/'];
    
    % Load audio labels and templates
    audioPath = [path 'audio/'];
    load([audioPath 'audioLabels.mat']);
    template = audioread([audioPath 'singing/template.wav']);
    motif = audioread([audioPath 'playback/motif.wav']);
    
    % Load cluster data and recording length
    dataPath = [path 'data/'];
    load([dataPath 'clusters.mat']);
    load([dataPath 'lenRec.mat']);
    
    % Read good units from Excel file
    opts = detectImportOptions([path bird '_goodSiUnits.xlsx']);
    opts.DataRange = 'A2';
    opts.VariableNamesRange = 'A1';
    goodCells = readtable([path bird '_goodSiUnits.xlsx'], opts);
    
    % Extract necessary data from table
    toPrint = (goodCells{:,3});
    startRec = (goodCells{:,4});
    stopRec = (goodCells{:,5});
    
    numUnits = length(clusters);
    disp(bird);
    
    %% Process data for plotting
    for n = 1:numUnits
        if toPrint(n) == 1
            % Adjust start and stop times based on recording length
            startRecFixed = startRec(n) * fs;
            startRecFixed(startRecFixed == 0) = 1;
            if strcmp(stopRec{n}, 'end')
                stopRecFixed = lenRec;
            else
                stopRecFixed = str2num(stopRec{n}) * fs;
            end
            
            % Define motif range
            motifRange = [floor(startRecFixed/fs) ceil(stopRecFixed/fs)];
            
            % Process each event label
            for i = 1:length(labels)
                if strcmp(labels{i}, 'sm')
                    audio = template;
                else
                    audio = motif;
                end
                
                % Retrieve event times
                Event = getEventTimes(audioPath, fs, audio, audioLabels, labels{i});
                takeMotif = Event.start >= motifRange(1) & Event.start <= motifRange(2);
                structElem = fieldnames(Event);
                for j = 2:5
                    Event.(structElem{j}) = Event.(structElem{j})(takeMotif);
                end
                assignin('base', [labels{i} 'Times'], Event);
                
                % Align spike times and compute histograms
                if ~isempty(Event.start)
                    [outStruct, countSp] = alignSpikeTimesSi(clusters(n), wind, Event, Event.center);
                else
                    outStruct = [];
                    countSp = [];
                end
                
                % Save aligned spike times and histograms
                assignin('base', [labels{i} 'Hist'], countSp);
                assignin('base', [labels{i} 'Spks'], outStruct);
            end
        end
    end
    disp(['Done processing data for ' bird]);
end
