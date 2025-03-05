function siQuantifications()
    %% Setup and Initialization
    clear; clc; close all; warning off;
    
    % Parameters
    fs = 30000;
    fontsize = 15;
    birds = {'C22','C24','C50','C52','C57'};
    addpath('C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Fcns\');
    labels = {'sm','pm','pms'};
    sigwin10 = 0.01; 
    sigwin50 = 0.05;
    calcBuffer = 0.03; 
    buffer = 0.5;
    lim = [1e-3 1e1];
    edges = -3:0.1:1;
    spontTrials = 50;
    windSpont = [0 5];
    binSparse = 0.005;
    
    badFields = {'ifrBottom','isis','spontIsis','spontOnIsis','spontOffIsis','spontSpks',...
        'spontOnSpks','spontOffSpks','spontAwakePrePbFR','spontAwakePrePbCV',...
        'spontOffEarlyIsis','spontOffDeepIsis','spontAwakePrePbBurst',...
        'nonReplayBurstProp','replayBurstProp','spontAwakeFR','spontAwakeCV',...
        'spontAwakeBurst','spontSleepDeepFR','songPrecInd','pbAwakePrecInd',...
        'pbSleepPrecInd','spontRatioVigilanceE','nonReplayOverReplay','maxChannel'};
    
    baseDir = 'C:\Users\User\Dropbox (NYU Langone Health)\SiProbe\';
    
    %% Process each bird
    for w = 1:length(birds)
        birdID = birds{w};
        fprintf('Processing bird: %s\n', birdID);
        
        % Define directories
        birdDir   = fullfile(baseDir, birdID);
        dataPath  = fullfile(birdDir, 'data');
        audioPath = fullfile(birdDir, 'audio');
        spontSongDir = fullfile(baseDir, 'spontSong');
        if ~exist(spontSongDir, 'dir')
            mkdir(spontSongDir);
        end

        % Load data files (assumes variables like lightSeg exist in the mat files)
        load(fullfile(dataPath, 'lightSig.mat'), 'lightSeg'); 
        load(fullfile(dataPath, 'clusters.mat'));
        load(fullfile(dataPath, 'lenRec.mat'));
        load(fullfile(audioPath, 'audioLabels.mat'), 'audioLabels');
        template = audioread(fullfile(audioPath, 'singing', 'template.wav'));
        motif    = audioread(fullfile(audioPath, 'playback', 'motif.wav'));

        % Import good unit info from Excel
        opts = detectImportOptions(fullfile(birdDir, [birdID '_goodSiUnits.xlsx']));
        opts.DataRange = 'A2';
        opts.VariableNamesRange = 'A1';
        goodCells = readtable(fullfile(birdDir, [birdID '_goodSiUnits.xlsx']), opts);
        toPrint = goodCells{:,3};
        startRec = goodCells{:,4};
        stopRec  = goodCells{:,5};

        % Remove unwanted fields from clusters
        clusters = removeBadFields(clusters, badFields);

        %% Process each unit
        numUnits = length(clusters);
        for n = 1:numUnits
            if toPrint(n) ~= 1
                continue;
            end
            fprintf('  Processing unit %d\n', n);
            clusters(n) = processUnit(clusters(n), startRec(n), stopRec{n}, lenRec, fs, ...
                template, motif, audioPath, audioLabels, labels, calcBuffer, buffer, lim, ...
                binSparse, sigwin10, sigwin50, edges, spontTrials, windSpont, birdID, n);
            
            % --- Additional Sections ---
            % 1. Before motif stats for pm and pms events:
            clusters(n) = computeBeforeMotifStats(clusters(n), motif, calcBuffer, fs);
            
            % 2. Extract spontaneous segments (auditory segments)
            [SpontAwake, SpontSleepEarly, SpontSleepDeep, audTrials] = ...
                computeSpontSnips(startRec, stopRec, lightSeg, buffer, spontTrials, fs, n);
            
            % 3. Juxta/spontaneous snippet analysis
            clusters(n) = computeJuxtaStats(clusters(n), windSpont, SpontAwake, SpontSleepEarly, SpontSleepDeep, fs, n);
            
            % 4. Compute post-event firing rates (post singing and playback)
            clusters(n) = computePostEventFRs(clusters(n), fs, template, motif, calcBuffer);
            
            % 5. Combinatorial ISI and ratio statistics
            clusters(n) = computeCombinatorialStats(clusters(n), lim, edges, SpontAwake, SpontSleepDeep, lightSeg);
            
            % 6. Spontaneous song metrics (using pre-defined snip times per bird)
            clusters(n) = computeSpontSongMetrics(clusters(n), w, fs, lightSeg);
        end
        
        % Order fields and save outputs
        clusters = orderfields(clusters);
        cd(dataPath);
        save('clusters.mat', 'clusters');
    end
end

%% --- Helper Functions ---

% Remove unwanted fields from the clusters structure
function clusters = removeBadFields(clusters, badFields)
    for i = 1:length(badFields)
        if isfield(clusters, badFields{i})
            clusters = rmfield(clusters, badFields{i});
        end
    end
end

% Process an individual unit (basic filtering, event alignment, waveform and binary precision)
function unitStruct = processUnit(unitStruct, startRecVal, stopRecVal, lenRec, fs, ...
    template, motif, audioPath, audioLabels, labels, calcBuffer, buffer, lim, ...
    binbin, sigwin10, sigwin50, edges, spontTrials, windSpont, birdID, unitIndex)

    %% Fix recording times and filter spikes
    startRecFixed = startRecVal * fs;
    if startRecFixed == 0
        startRecFixed = 1;
    end
    if strcmp(stopRecVal, 'end')
        stopRecFixed = lenRec;
    else
        stopRecFixed = str2double(stopRecVal) * fs;
    end
    goodDataTimes = [startRecFixed/fs, stopRecFixed/fs];
    
    spkTimes = unitStruct.spikeTimes;
    spkTimes = spkTimes(spkTimes >= startRecFixed/fs & spkTimes <= stopRecFixed/fs);
    isis = diff(spkTimes);
    validIdx = isis > 0.0015;
    spkTimes = spkTimes([true; validIdx]); % Keep first spike plus valid ones
    unitStruct.ogTimes = unitStruct.spikeTimes;
    unitStruct.spikeTimes = spkTimes;
    
    % Define motif range for event selection
    motifRange = [ceil(startRecFixed/fs), floor(stopRecFixed/fs)];
    
    %% Process Event Alignments for each label (sm, pm, pms)
    for i = 1:length(labels)
        currentLabel = labels{i};
        if strcmp(currentLabel, 'sm')
            currentAudio = template;
        else
            currentAudio = motif;
        end
        EventOG = getEventTimes(audioPath, fs, currentAudio, audioLabels, currentLabel);
        % Restrict events to motifRange
        takeMotif = EventOG.start >= motifRange(1) & EventOG.stop <= motifRange(2);
        structElem = fieldnames(EventOG);
        for j = 2:5
            Event.(structElem{j}) = EventOG.(structElem{j})(takeMotif);
        end
        if isempty(Event.start)
            for j = 2:5
                Event.(structElem{j}) = nan;
            end
        end
        
        % Align spike times to events
        if strcmp(currentLabel, 'sm')
            [structs(i).outStruct, ~] = alignSpikeTimesSi(unitStruct, [0, length(currentAudio)/fs], Event, Event.center);
            whichTrials = 1:length(Event.center);
            for c = 2:5
                EventRand.(structElem{c}) = Event.(structElem{c})(whichTrials);
            end
            [~, psth] = alignSpikeTimesSi(unitStruct, [0, length(currentAudio)/fs], EventRand, EventRand.center);
            assignin('base', [currentLabel 'Psth'], psth);  % Keeping for compatibility
        else
            [structs(i).outStruct, ~] = alignSpikeTimesSi(unitStruct, [calcBuffer, length(currentAudio)/fs + calcBuffer], Event, Event.start);
        end
        assignin('base', [currentLabel 'Times'], Event);
        assignin('base', [currentLabel 'Spks'], structs(i).outStruct);
        
        % Compute trial metrics: rate, max IFR, burstiness
        [meanRate, burstyMean, meanMaxIFR, rates] = computeTrialMetrics(structs(i).outStruct, length(currentAudio)/fs);
        assignin('base', [currentLabel 'FRAll'], rates);
        assignin('base', [currentLabel 'FR'], meanRate);
        assignin('base', [currentLabel 'Burst'], burstyMean);
        assignin('base', [currentLabel 'MaxIFR'], meanMaxIFR);
    end
    
    %% Binary Song Precision (vectorized correlation)
    [songCorr, songRtoZ, songMatch] = computeBinaryPrecision(unitStruct, fs, template, buffer, binbin, structs(1).outStruct, Event);
    unitStruct.songCorr = songCorr;
    unitStruct.songRtoZ = songRtoZ;
    unitStruct.songMatch = round(songMatch, 4);
    
    %% Waveform Metrics
    [normWav, wfRise, wfRise70, wfFWHM] = computeWaveformMetrics(unitStruct.meanWF, fs);
    unitStruct.wfNorm = normWav;
    unitStruct.wfRise = wfRise;
    unitStruct.wfRise70 = wfRise70;
    unitStruct.wfFWHM = wfFWHM;
    
    % Additional basic processing can be added here if needed.
end

% Compute trial metrics for aligned spike times
function [meanRate, burstyMean, meanMaxIFR, rates] = computeTrialMetrics(alignedStruct, eventDuration)
    numTrials = length(alignedStruct.eventTrials);
    rates = nan(numTrials,1);
    maxIFRs = nan(numTrials,1);
    burstiness = nan(numTrials,1);
    for k = 1:numTrials
        trialIdx = alignedStruct.eventTrials == k;
        trialSpks = alignedStruct.eventTimesAligned(trialIdx);
        rates(k) = length(trialSpks) / eventDuration;
        trialIsis = diff(trialSpks);
        if ~isempty(trialIsis)
            ifrs = 1 ./ trialIsis;
            maxIFRs(k) = max(ifrs);
            burstiness(k) = sum(trialIsis < 0.01) / length(trialIsis);
        else
            maxIFRs(k) = nan;
            burstiness(k) = 0;
            if isempty(trialSpks)
                rates(k) = 0;
            end
        end
    end
    meanRate = round(mean(rates), 2);
    burstyMean = round(nanmean(burstiness), 3);
    if isnan(burstyMean)
        burstyMean = 0;
    end
    meanMaxIFR = round(nanmean(maxIFRs), 2);
end

% Compute binary precision metrics using vectorized correlation
function [rVal, zVal, matchVal] = computeBinaryPrecision(unitStruct, fs, audioTemplate, buff, binbin, alignedStruct, Event)
    songHalf = length(audioTemplate) / fs / 2;
    windSong = [-(buff + songHalf), (buff + songHalf)];
    fullSpec = diff(windSong);
    allBins = ceil(fullSpec / binbin);
    numTrials = length(Event.center);
    
    series = zeros(numTrials, allBins);
    poss = ceil(-allBins/2) : ceil(allBins/2)-1;
    
    for k = 1:numTrials
        trialIdx = alignedStruct.eventTrials == k;
        trialSpks = alignedStruct.eventTimesAligned(trialIdx);
        whichBin = ceil(trialSpks ./ binbin);
        binVec = zeros(1, length(poss));
        binVec(ismember(poss, whichBin)) = 1;
        series(k,:) = binVec;
    end
    songHalfBins = ceil(songHalf / binbin);
    midSeries = find(poss == 0);
    songBins = midSeries - songHalfBins : midSeries + songHalfBins;
    series = series(:, songBins);
    
    r = 1 - squareform(pdist(series, 'correlation'));
    r(isnan(r)) = 0;
    rLower = tril(r, -1);
    rVal = nanmean(rLower(rLower ~= 0));
    
    z = atanh(r);
    z(isinf(z)) = atanh(0.9999);
    zLower = tril(z, -1);
    zVal = nanmean(zLower(zLower ~= 0));
    
    % "All-or-nothing" match metric (loop retained for clarity)
    match = nan(size(series,1));
    for i = 1:size(series,1)
        for j = 1:size(series,1)
            comp1 = series(i,:);
            comp2 = series(j,:);
            same = nan(size(comp1));
            for k = 1:length(comp1)
                if comp1(k) == comp2(k)
                    if comp1(k) == 1
                        same(k) = 1;
                    else
                        same(k) = nan;
                    end
                else
                    same(k) = 0;
                end
            end
            match(i,j) = nanmean(same);
        end
    end
    matchLower = tril(match, -1);
    matchVal = nanmean(matchLower(:));
end

% Compute waveform metrics from the mean waveform
function [normWav, trPkTime, trPk70, width] = computeWaveformMetrics(wav, fs)
    meanSubWav = wav - mean(wav);
    maxAmp = max(abs(meanSubWav));
    normWav = meanSubWav / maxAmp;
    
    factor = 1000;
    normWavInterp = interp(normWav, factor);
    normWavInterp = round(normWavInterp, 5);
    
    trVal = min(normWavInterp);
    trPos = find(islocalmin(normWavInterp, 'MaxNumExtrema', 1), 1, 'first');
    if isempty(trPos)
        trPkTime = nan;
        trPk70   = nan;
        width    = nan;
        return;
    end
    pkIdxRange = trPos:length(normWavInterp);
    pkVal = max(normWavInterp(pkIdxRange));
    pkPos = trPos - 1 + find(islocalmax(normWavInterp(pkIdxRange), 'MaxNumExtrema', 1), 1, 'first');
    trPkTime = (pkPos - trPos) / fs * 1000 / factor;
    
    pkVal70 = trVal + 0.7 * (pkVal - trVal);
    closeEnuf = find(abs(normWavInterp - pkVal70) < 0.0015);
    inBet = closeEnuf(closeEnuf < pkPos & closeEnuf > trPos);
    if isempty(inBet)
        trPk70 = nan;
    else
        pkPos70 = inBet(find(islocalmin(abs(normWavInterp(inBet) - pkVal70), 'MaxNumExtrema', 1), 1, 'first'));
        trPk70 = (pkPos70 - trPos) / fs * 1000 / factor;
    end
    
    pkVal50 = trVal + 0.5 * (pkVal - trVal);
    closeEnuf50 = find(abs(normWavInterp - pkVal50) < 0.0015);
    beforePk = closeEnuf50(closeEnuf50 < pkPos);
    if length(beforePk) < 2
        width = nan;
    else
        pkPos50 = beforePk(find(islocalmin(abs(normWavInterp(beforePk) - pkVal50), 'MaxNumExtrema', 2), 1, 'first'));
        width = diff(pkPos50) / fs * 1000 / factor;
    end
end

%% Additional Sections

% 1. Compute "Before Motif" statistics for pm and pms events
function unitStruct = computeBeforeMotifStats(unitStruct, motif, calcBuffer, fs)
    % These variables are retrieved from the base workspace to preserve the current structure.
    % For further optimization, consider passing pmTimes and pmsTimes as arguments.
    pmTimes = evalin('base','pmTimes');
    pmsTimes = evalin('base','pmsTimes');
    structElem = fieldnames(pmTimes);
    
    beforePmTimes = struct();
    beforePmsTimes = struct();
    for i = 2:5
        if i ~= 5
            beforePmTimes.(structElem{i}) = pmTimes.(structElem{i}) - length(motif)/fs;
            beforePmsTimes.(structElem{i}) = pmsTimes.(structElem{i}) - length(motif)/fs;
        else
            beforePmTimes.(structElem{i}) = pmTimes.(structElem{i});
            beforePmsTimes.(structElem{i}) = pmsTimes.(structElem{i});
        end
    end
    [beforePmSpks, ~] = alignSpikeTimesSi(unitStruct, [calcBuffer, length(motif)/fs+calcBuffer], beforePmTimes, beforePmTimes.start);
    [beforePmsSpks, ~] = alignSpikeTimesSi(unitStruct, [calcBuffer, length(motif)/fs+calcBuffer], beforePmsTimes, beforePmTimes.start);
    
    % Compute firing rates for each trial with preallocated arrays
    for idx = 1:2
        if idx == 1
            Event = beforePmTimes;
            outStruct = beforePmSpks;
        else
            Event = beforePmsTimes;
            outStruct = beforePmsSpks;
        end
        numTrials = length(Event.center);
        rates = nan(numTrials,1);
        for k = 1:numTrials
            trial = outStruct.eventTrials == k;
            trialSpks = outStruct.eventTimesAligned(trial);
            rates(k) = length(trialSpks) / (length(motif)/fs);
            if isempty(trialSpks)
                rates(k) = 0;
            end
        end
        meanRate = round(mean(rates), 2);
        if idx == 1
            unitStruct.prePbAwakeFR = meanRate;
            unitStruct.prePbAwakeFRAll = rates;
        else
            unitStruct.prePbSleepFR = meanRate;
            unitStruct.prePbSleepFRAll = rates;
        end
    end
end

% 2. Extract spontaneous segments (spont snips) from non-auditory intervals
function [SpontAwake, SpontSleepEarly, SpontSleepDeep, audTrials] = computeSpontSnips(startRec, stopRec, lightSeg, buffer, spontTrials, fs, n)
    % Placeholder: Replace audTrials with actual auditory trial times from your pipeline.
    audTrials = [200, 210; 400, 410; 800, 810];
    
    % Build non-auditory segments between auditory trials
    nonAudTrials = [];
    for i = 1:size(audTrials,1)-1
        nonAudTrials = [nonAudTrials; audTrials(i,2)+buffer, audTrials(i+1,1)-buffer];
    end
    segs = sortrows([startRec/fs, audTrials(1,1); nonAudTrials; audTrials(end,2), stopRec/fs]);
    
    % Split segments at lightSeg (in seconds)
    splitInd = find(segs(:,2) > lightSeg & segs(:,1) < lightSeg, 1);
    segs = [segs(1:splitInd-1, :); [segs(splitInd,1), lightSeg]; [lightSeg, segs(splitInd,2)]; segs(splitInd+1:end,:)];
    
    necLen = 10;
    segLengths = segs(:,2) - segs(:,1);
    longSegs = segLengths >= necLen;
    longEnuf = segLengths >= 5 & segLengths < necLen;
    spontSegs = segs(longEnuf, :);
    for i = 1:length(segLengths)
        if longSegs(i)
            left = floor(segLengths(i)/necLen);
            shortSegs = [(segs(i,1):necLen:segs(i,2))', (segs(i,1):necLen:segs(i,2))' + necLen];
            if size(shortSegs,1) > left
                shortSegs(end,:) = [shortSegs(end-1,2), segs(i,2)];
            end
            spontSegs = [spontSegs; shortSegs];
        end
    end
    spontSegs = sortrows(spontSegs(spontSegs(:,2)-spontSegs(:,1) >= 5, :));
    
    % Split segments into on (before lightSeg) and off (after lightSeg)
    split = find(spontSegs(:,1) == lightSeg, 1);
    if isempty(split)
        segsOn = spontSegs;
        segsOff = [];
    else
        segsOn  = spontSegs(1:split-1, :);
        segsOff = spontSegs(split:end, :);
    end
    deepCut = lightSeg + 600; % seconds
    whichEarly = segsOff(:,2) <= deepCut;
    segsOffEarly = segsOff(whichEarly, :);
    segsOffDeep  = segsOff(~whichEarly, :);
    
    % Generate random 5-second snips from each group using preallocation
    SpontAwake = generateSnips(segsOn, spontTrials, 5, n);
    SpontSleepDeep = generateSnips(segsOffDeep, spontTrials, 5, n);
    SpontSleepEarly = generateSnips(segsOffEarly, spontTrials, 5, n);
end

% Helper: generate random snips from given segments with preallocation
function SnipStruct = generateSnips(segsGroup, spontTrials, snipLength, seedVal)
    rng(seedVal);
    numSegs = size(segsGroup,1);
    numSelect = min(numSegs, spontTrials);
    SnipStruct.start = zeros(numSelect,1);
    SnipStruct.stop = zeros(numSelect,1);
    SnipStruct.center = zeros(numSelect,1);
    count = 1;
    idxSelect = randperm(numSegs, numSelect);
    for i = idxSelect
        a = ceil(segsGroup(i,1));
        b = ceil(segsGroup(i,2)) - snipLength;
        snip = a + (b - a) * rand;
        SnipStruct.start(count,1) = snip;
        SnipStruct.stop(count,1) = snip + snipLength;
        SnipStruct.center(count,1) = snip + snipLength/2;
        count = count + 1;
    end
    SnipStruct.warp = ones(numSelect,1);
end

% 3. Compute juxta (spontaneous) metrics from spontaneous segments
function unitStruct = computeJuxtaStats(unitStruct, windSpont, SpontAwake, SpontSleepEarly, SpontSleepDeep, fs, n)
    spontLabels = {'spontAwake','spontSleepEarly','spontSleepDeep'};
    spontTimes = {SpontAwake, SpontSleepEarly, SpontSleepDeep};
    % Use goodDataTimes from unitStruct (from original spike times)
    goodDataTimes = [min(unitStruct.ogTimes), max(unitStruct.ogTimes)];
    
    for s = 1:length(spontLabels)
        Sevent = alignSpikeTimesSi(unitStruct, windSpont, spontTimes{s}, spontTimes{s}.start);
        numTrials = length(spontTimes{s}.start);
        bursty = nan(numTrials,1);
        rates  = nan(numTrials,1);
        cv     = nan(numTrials,1);
        for k = 1:numTrials
            if spontTimes{s}.start(k) >= goodDataTimes(1) && spontTimes{s}.stop(k) <= goodDataTimes(2)
                within = 1;
            else
                within = 0;
            end
            trial = Sevent.eventTrials == k;
            trialSpks = Sevent.eventTimesAligned(trial);
            trialIsi = diff(trialSpks);
            if isempty(trialIsi) && within
                bursty(k) = 0;
                cv(k) = nan;
                rates(k) = 0;
                if length(trialSpks) == 1
                    rates(k) = length(trialSpks) / range(windSpont);
                end
            elseif isempty(trialIsi) && ~within
                bursty(k) = nan;
                cv(k) = nan;
                rates(k) = nan;
            else
                bursty(k) = nansum(trialIsi < 0.01) / length(trialIsi);
                cv(k) = nanstd(trialIsi) / nanmean(trialIsi);
                rates(k) = length(trialSpks) / range(windSpont);
            end
        end
        unitStruct.(sprintf('%sFR',spontLabels{s})) = round(nanmean(rates),2);
        unitStruct.(sprintf('%sCV',spontLabels{s})) = nanmean(cv);
        unitStruct.(sprintf('%sBurst',spontLabels{s})) = round(nanmean(bursty),3);
        % Retain Sevent in base workspace if needed (could be passed as an output instead)
        assignin('base',[spontLabels{s} 'Spks'], Sevent);
    end
end

% 4. Compute post-event firing rates (for song and playback)
function unitStruct = computePostEventFRs(unitStruct, fs, template, motif, calcBuffer)
    % Post singing rates (using smTimes from base workspace)
    smTimes = evalin('base','smTimes');
    bouts = find(diff(smTimes.stop) > 2);
    smTimesLast.start = smTimes.start(bouts);
    smTimesLast.stop  = smTimes.stop(bouts);
    smTimesLast.warp  = ones(length(bouts),1);
    [smPost, ~] = alignSpikeTimesSi(unitStruct, [length(template)/fs, buffer+length(template)/fs], smTimesLast, smTimesLast.start);
    smPostRates = nan(length(smTimesLast.start),1);
    for k = 1:length(smTimesLast.start)
        trial = smPost.eventTrials == k;
        trialSpks = smPost.eventTimesAligned(trial);
        if isempty(trialSpks)
            smPostRates(k) = 0;
        else
            smPostRates(k) = length(trialSpks) / (length(template)/fs);
        end
    end
    unitStruct.songFRPost = round(mean(smPostRates),2);
    
    % Post playback rates for pm and pms events (using pmTimes, pmsTimes from base)
    pmTimes = evalin('base','pmTimes');
    [pmPost, ~] = alignSpikeTimesSi(unitStruct, [length(motif)/fs+calcBuffer, buffer+length(motif)/fs+calcBuffer], pmTimes, pmTimes.start);
    pmPostRates = nan(length(pmTimes.center),1);
    for k = 1:length(pmTimes.center)
        trial = pmPost.eventTrials == k;
        trialSpks = pmPost.eventTimesAligned(trial);
        if isempty(trialSpks)
            pmPostRates(k) = 0;
        else
            pmPostRates(k) = length(trialSpks) / (length(motif)/fs);
        end
    end
    unitStruct.pbAwakeFRPost = round(mean(pmPostRates),2);
    
    pmsTimes = evalin('base','pmsTimes');
    [pmsPost, ~] = alignSpikeTimesSi(unitStruct, [length(motif)/fs+calcBuffer, buffer+length(motif)/fs+calcBuffer], pmsTimes, pmsTimes.start);
    pmsPostRates = nan(length(pmsTimes.center),1);
    for k = 1:length(pmsTimes.center)
        trial = pmsPost.eventTrials == k;
        trialSpks = pmsPost.eventTimesAligned(trial);
        if isempty(trialSpks)
            pmsPostRates(k) = 0;
        else
            pmsPostRates(k) = length(trialSpks) / (length(motif)/fs);
        end
    end
    unitStruct.pbSleepFRPost = round(mean(pmsPostRates),2);
end

% 5. Compute combinatorial ISI and state ratio statistics
function unitStruct = computeCombinatorialStats(unitStruct, lim, edges, SpontAwake, SpontSleepDeep, lightSeg)
    % Use the current spike times from unitStruct and split into on/off segments.
    spontSpks = unitStruct.spikeTimes;
    spontOnSpks = spontSpks(spontSpks <= lightSeg);
    spontOffSpks = spontSpks(spontSpks > lightSeg);
    
    deepCut = lightSeg + 600;
    offEisi = diff(spontOffSpks(spontOffSpks <= deepCut));
    offDisi = diff(spontOffSpks(spontOffSpks > deepCut));
    spontAwakeIsis = diff(spontOnSpks);
    
    ons = spontAwakeIsis(spontAwakeIsis >= lim(1) & spontAwakeIsis <= lim(2));
    offEs = offEisi(offEisi >= lim(1) & offEisi <= lim(2));
    offDs = offDisi(offDisi >= lim(1) & offDisi <= lim(2));
    
    allonih = histcounts(log10(ons), edges, 'Normalization', 'probability');
    alloffihe = histcounts(log10(offEs), edges, 'Normalization', 'probability');
    alloffihd = histcounts(log10(offDs), edges, 'Normalization', 'probability');
    
    ec = corrcoef(allonih, alloffihe);
    dc = corrcoef(allonih, alloffihd);
    unitStruct.spontStateCorrE = round(ec(2,1),2);
    unitStruct.spontStateCorrD = round(dc(2,1),2);
    
    % Compute vigilance ratios (assumes corresponding FR fields exist)
    unitStruct.spontRatioVigilanceE = round((unitStruct.spontAwakeFR - unitStruct.spontSleepEarlyFR) / ...
        (unitStruct.spontAwakeFR + unitStruct.spontSleepEarlyFR), 2);
    unitStruct.spontRatioVigilanceD = round((unitStruct.spontAwakeFR - unitStruct.spontSleepDeepFR) / ...
        (unitStruct.spontAwakeFR + unitStruct.spontSleepDeepFR), 2);
    
    unitStruct.spontAwakeIsis = spontAwakeIsis;
    unitStruct.spontSleepEarlyIsis = offEisi;
    unitStruct.spontSleepDeepIsis = offDisi;
    
    % Compute change ratios between song and spontaneous, and playback states
    unitStruct.songChange = round((unitStruct.songFR - unitStruct.spontAwakeFR) / (unitStruct.songFR + unitStruct.spontAwakeFR), 2);
    unitStruct.songPostChange = round((unitStruct.songFRPost - unitStruct.spontAwakeFR) / (unitStruct.songFRPost + unitStruct.spontAwakeFR), 2);
    unitStruct.pbAwakeChange = round((unitStruct.pbAwakeFR - unitStruct.prePbAwakeFR) / (unitStruct.pbAwakeFR + unitStruct.prePbAwakeFR), 2);
    unitStruct.pbAwakeChangeOld = round((unitStruct.pbAwakeFR - unitStruct.spontAwakeFR) / (unitStruct.pbAwakeFR + unitStruct.spontAwakeFR), 2);
    unitStruct.pbAwakePostChange = round((unitStruct.pbAwakeFRPost - unitStruct.spontAwakeFR) / (unitStruct.pbAwakeFRPost + unitStruct.spontAwakeFR), 2);
    unitStruct.pbSleepChange = round((unitStruct.pbSleepFR - unitStruct.prePbSleepFR) / (unitStruct.pbSleepFR + unitStruct.prePbSleepFR), 2);
    unitStruct.pbSleepChangeOld = round((unitStruct.pbSleepFR - unitStruct.spontSleepDeepFR) / (unitStruct.pbSleepFR + unitStruct.spontSleepDeepFR), 2);
    unitStruct.pbSleepPostChange = round((unitStruct.pbSleepFRPost - unitStruct.spontSleepDeepFR) / (unitStruct.pbSleepFRPost + unitStruct.spontSleepDeepFR), 2);
end

% 6. Compute spontaneous song metrics (snips around song)
function unitStruct = computeSpontSongMetrics(unitStruct, birdIndex, fs, lightSeg)
    % Define spontaneous song snip times per bird (example values)
    switch birdIndex
        case 1
            spontSongSnips = [250; 280; 5875; 5970; 6040];
        case 2
            spontSongSnips = [1060; 1235; 1530; 1600; 6870];
        case 3
            spontSongSnips = [1815; 2285; 2660; 3880; 4560];
        case 4
            spontSongSnips = [1485; 1700; 4330; 4933; 5216];
        case 5
            spontSongSnips = [118; 1614; 2945; 3540; 6480];
        otherwise
            spontSongSnips = [];
    end
    spontSnipLength = 10;
    numSpontSnips = length(spontSongSnips);
    spontSong.start = spontSongSnips;
    spontSong.stop  = spontSongSnips + spontSnipLength;
    spontSong.center = spontSong.start + spontSnipLength/2;
    spontSong.warp = ones(numSpontSnips,1);
    
    % Align spikes for spontaneous song segments
    [spontSongSpks, ~] = alignSpikeTimesSi(unitStruct, [0, spontSnipLength], spontSong, spontSong.start);
    bursty = nan(numSpontSnips,1);
    rates  = nan(numSpontSnips,1);
    cv     = nan(numSpontSnips,1);
    for k = 1:numSpontSnips
        if spontSong.start(k) >= min(unitStruct.ogTimes) && spontSong.stop(k) <= max(unitStruct.ogTimes)
            within = 1;
        else
            within = 0;
        end
        trial = spontSongSpks.eventTrials == k;
        trialSpks = spontSongSpks.eventTimesAligned(trial);
        trialIsi = diff(trialSpks);
        if isempty(trialIsi) && within
            bursty(k) = 0;
            cv(k) = nan;
            rates(k) = 0;
            if length(trialSpks) == 1
                rates(k) = length(trialSpks) / spontSnipLength;
            end
        elseif isempty(trialIsi) && ~within
            bursty(k) = nan;
            cv(k) = nan;
            rates(k) = nan;
        else
            bursty(k) = nansum(trialIsi < 0.01) / length(trialIsi);
            cv(k) = nanstd(trialIsi) / nanmean(trialIsi);
            rates(k) = length(trialSpks) / spontSnipLength;
        end
    end
    unitStruct.spontSongRate = nanmean(rates);
    unitStruct.spontSongBurst = nanmean(bursty);
    unitStruct.spontSongCV   = nanmean(cv);
end
