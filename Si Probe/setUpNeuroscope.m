clc; warning off;

% Add necessary script functions directory to path
addpath 'C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Fcns\'

% Define sampling frequency
fs = 30000;

% List of bird identifiers
allbirds = {'C22'; 'C24'; 'C50'; 'C52'; 'C57'};

% Select which bird from the list
w = 1;
disp(allbirds(w))  % Display selected bird

% Define file path for the selected bird
bird = char(allbirds(w));
path = ['C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/' bird '/'];

% Load cluster data
load([path 'data/clusters.mat']);

% Read the list of good units from Excel
opts = detectImportOptions([path bird '_goodSiUnits.xlsx']);
opts.DataRange = 'A2';
opts.VariableNamesRange = 'A1';
goodCells = readtable([path bird '_goodSiUnits.xlsx'], opts);

% Extract the logical mask for good clusters
toPrint = logical(goodCells{:,3});

% Initialize counter and clear previous structure
count = 0;
clear goodClusters;

% Filter clusters based on waveform and song sparsity properties
for i = 1:length(clusters)
    duo = clusters(i).wfRise * clusters(i).songSparse;  % Compute filtering metric
    split = duo < 0.25;  % Define threshold condition
    if toPrint(i) && ~split  % Retain only good clusters that meet the condition
        count = count + 1;
        goodClusters(count).clusterID = clusters(i).clusterID;
        goodClusters(count).class = count;
        goodClusters(count).spikeTimes = clusters(i).spikeTimes;
    end
end

% Create directory for storing neuroscope project output
mkdir([path 'data/neuroscopeProj/']);
cd([path 'data/neuroscopeProj/']);

% Write the filtered clusters to neuroscope-compatible format
writeNeuroscopeClass(goodClusters, fs, bird);

disp('done');  % Indicate script completion
