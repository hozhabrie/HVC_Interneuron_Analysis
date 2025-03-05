function [starts, stops, centers, warps] = findMotifsSi(audioF, templateF, thresh)
% findMotifsSi - Identifies motifs in an audio file based on a template.
%
% Inputs:
%   audioF    - String: Filename of the .wav file containing the trial (~25 s)
%   templateF - String: Filename of the .wav file containing the template (~1 s)
%   thresh    - Scalar: Threshold for correlation (between 0 and 1)
%
% Outputs:
%   starts  - Start times of detected motifs
%   stops   - Stop times of detected motifs
%   centers - Center points of motifs
%   warps   - Warping factors (length of detected motif / length of template)

% Validate inputs
if ~isfile(audioF)
    error('Audio file not found: %s', audioF);
end
if ~isfile(templateF)
    error('Template file not found: %s', templateF);
end
if thresh < 0 || thresh > 1
    error('Threshold must be between 0 and 1.');
end

% Set parameters
ds = 100;        % Downsampling factor
cutoutAddon = 0.2; % Additional time added to cutouts (seconds)
samplingRate = 30000; % Hz
cutoutAddon = cutoutAddon * samplingRate / ds;

% Load template and preprocess
motifTemplateFull = audioread(templateF);
[audioMin] = audioBounds(motifTemplateFull, ds, 1);
motifTemplate = audioMin; % Final motif template

% Generate warped templates (default range: 98% to 102%)
warpingRange = 980:2:1020;
numWarps = numel(warpingRange);
Warping.templates = cell(1, numWarps);
Warping.warp_factor = cell(1, numWarps);

for k = 1:numWarps
    Warping.templates{k} = resample(motifTemplate, warpingRange(k), 1000);
    Warping.warp_factor{k} = warpingRange(k) / 1000;
end

% Load trial audio
currTrialFull = audioread(audioF);
[audioMin] = audioBounds(currTrialFull, ds, 1);
currTrial = audioMin;

% Perform sliding cross-correlation
fprintf('\nFinding motifs...\n');
offset = length(motifTemplate);
slidingCrossCorr = zeros(1, length(currTrial) - offset);

for i = 1:length(slidingCrossCorr)
    slidingCrossCorr(i) = corr(motifTemplate', currTrial(i:i+offset-1)');
end

% Detect peaks above threshold
peakStart = find(slidingCrossCorr(1:end-1) <= thresh & slidingCrossCorr(2:end) >= thresh);
peakEnd = find(slidingCrossCorr(1:end-1) >= thresh & slidingCrossCorr(2:end) <= thresh);
peaks = [peakStart(:), peakEnd(:)];

% Identify motifs within peak regions
motifPos = [];
for y = 1:size(peaks, 1)
    [peakValue, peakIdx] = max(slidingCrossCorr(peaks(y,1) - offset/2 : peaks(y,2) - offset/2));
    
    if y < 3
        motifPos = [motifPos; peakValue, peakIdx + peaks(y,1) - 1];
    else
        if peakIdx + peaks(y,1) - 1 > motifPos(end,2) + cutoutAddon
            motifPos = [motifPos; peakValue, peakIdx + peaks(y,1) - 1];
        elseif peakValue > motifPos(end,1)
            motifPos(end, :) = [];
            motifPos = [motifPos; peakValue, peakIdx + peaks(y,1) - 1];
        end
    end
end

% Output detected motifs
numMotifs = size(motifPos, 1);
fprintf('\nMotifs detected: %d\n', numMotifs);

if numMotifs == 0
    fprintf('Try lowering the threshold if motifs are expected.\n');
    starts = [];
    stops = [];
    centers = [];
    warps = [];
    return;
end

motifOnset = motifPos(:, 2)' - offset/2;
motifOffset = motifPos(:, 2)' + offset/2;

% Determine best warping factor for each motif
warpFactors = zeros(1, numMotifs);
motifCenters = zeros(1, numMotifs);
motifStarts = zeros(1, numMotifs);
motifEnds = zeros(1, numMotifs);

for e = 1:numMotifs
    motifCutout = currTrial(motifOnset(e) - offset/2 - cutoutAddon : motifOffset(e) + offset/2 + cutoutAddon);

    maxCross = zeros(numWarps, 3);
    for f = 1:numWarps
        warpOffset = length(Warping.templates{f});
        slidingCrossCorrWarp = zeros(1, length(motifCutout) - warpOffset);

        for i = 1:length(slidingCrossCorrWarp)
            slidingCrossCorrWarp(i) = corr(Warping.templates{f}', motifCutout(i:i+warpOffset-1)');
        end

        [maxCross(f,2), maxCross(f,3)] = max(slidingCrossCorrWarp);
        maxCross(f,1) = warpingRange(f) / 1000;
    end

    % Select the best warping factor
    [~, maxWarpIdx] = max(maxCross(:,2));
    warpFactors(e) = maxCross(maxWarpIdx, 1);
    motifCenters(e) = (motifPos(e,2) + maxCross(maxWarpIdx,3) + length(Warping.templates{maxWarpIdx})/2) * ds;
    motifStarts(e) = motifCenters(e) - (length(Warping.templates{maxWarpIdx}) * ds) / 2;
    motifEnds(e) = motifCenters(e) + (length(Warping.templates{maxWarpIdx}) * ds) / 2;
end

% Convert results to time (seconds)
starts = motifStarts / samplingRate;
stops = motifEnds / samplingRate;
centers = motifCenters / samplingRate;
warps = warpFactors;
