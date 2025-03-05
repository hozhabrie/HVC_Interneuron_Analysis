%% Load everything
clear all; close all; clc; warning('off', 'all');

% Define base directory and load goodCells table
dirBase = 'C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/';
opts = detectImportOptions(fullfile(dirBase, 'goodCellsnew.xlsx'));
opts = setvartype(opts, 'stdmax', 'string');
goodCells = readtable(fullfile(dirBase, 'goodCellsnew.xlsx'), opts);

birds    = string(goodCells{:,1});
cells    = string(goodCells{:,2});
toPrint  = goodCells{:,5};
startRec = goodCells{:,6};
stopRec  = goodCells{:,7};
std      = string(goodCells{:,4});
identity = string(goodCells{:,9});

deflect = contains(std, '+');

fs = 50000;
wind = 0.002 * fs;
fontsize = 15;

set(0, 'DefaultFigureVisible', 'on');
counting = 0;

%% Loop through each row in goodCells
for number = 1:height(goodCells)
    if toPrint(number) == 1 && identity(number) ~= 'x'
        % Clear variables except the ones we need across iterations
        clearvars -except diffRise diffRise70 diffWidth fontsize fs wind number identity goodCells birds cells toPrint dirBase startRec stopRec deflect counting SpontsOnTimes SpontsOffTimes AwakePB AwakeSpont SleepPB SleepSpont MotifsOnTimes MotifsOffTimes AwakeHist SleepHist

        close all;
        % Build path for current bird
        path = fullfile(dirBase, birds(number));
        f = char(cells(number));
        % Get bird name from path
        out = regexp(path, '/', 'split');
        name = string(out(end-1));
        disp(name)
        disp(f)

        % Load high-pass times file
        load(fullfile(path, 'hp200', f, ['times_' f '_hp200.mat']));

        % Load fullRecord if it exists; otherwise, load from ABF and save it
        if isfile(fullfile(path, f, 'fullRecord.mat'))
            load(fullfile(path, f, 'fullRecord.mat'));
        else
            fullRecord = abfload(fullfile(path, f, [f '_filteredspikes.abf']));
            fullRecord(:,2:3) = [];
            cd(fullfile(path, f));
            save('fullRecord.mat', 'fullRecord');
        end

        % Determine good data times based on startRec and stopRec
        startRecFixed = startRec(number) * fs;
        if startRecFixed == 0, startRecFixed = 1; end
        if strcmp(stopRec{number}, 'end')
            stopRecFixed = size(fullRecord, 1);
        else
            stopRecFixed = str2num(stopRec{number}) * fs;
        end
        goodDataTimes = startRecFixed:stopRecFixed;

        % Load high-pass filtered data
        load(fullfile(path, 'hp200', f, [f '_hp200.mat']));
        spikes = data(:,1);
        clear data;
        eyeTrig = fullRecord(:, 2);
        clear fullRecord;

        lengthRecording = length(goodDataTimes) / fs;
        motif = audioread(fullfile(path, 'motif.wav')); % Change if needed
        load(fullfile(path, f, 'MotifTimes.mat'));

        % Filter motif times within recording range
        structElem = fieldnames(Motif);
        motifRange = [floor(startRecFixed/fs), stopRecFixed/fs];
        elem = Motif.start;
        takeMotif = elem >= motifRange(1) & elem <= motifRange(2);
        for i = 2:5
            Motif.(structElem{i}) = Motif.(structElem{i})(takeMotif);
        end

        % Store basic cluster info
        clusters.bird = name;
        clusters.cell = string(f);
        clusters.recLength = lengthRecording;

        %% Get spike times (from ms)
        spkTimesAll = round(cluster_class(:,2) / 1000 * fs);
        spkTimes = spkTimesAll(cluster_class(:,1) == 1);
        spkTimes = spkTimes(spkTimes >= startRecFixed & spkTimes <= stopRecFixed);
        isis = diff(spkTimes) / fs;
        fixed = isis > 0.0015;
        spkTimes = spkTimes(fixed == 1);
        clusters.spikeTimes = spkTimes / fs;

        %% Plot spike shape
        timeWindow = wind + 20;
        xAxisWaveForm = linspace(-timeWindow/(fs*2)*1000, timeWindow/(fs*2)*1000, timeWindow+1);
        allWF = zeros(timeWindow+1, length(spkTimes));
        startInd = 55;
        stopInd  = 65;
        for i = 1:length(spkTimes)
            span = spkTimes(i) - timeWindow/2 : spkTimes(i) + timeWindow/2;
            allWF(:, i) = spikes(span);
        end

        rng(number)
        randChoose = randperm(length(spkTimes), min(length(spkTimes), 200));
        randWindow = allWF(:, randChoose);
        meanWF = mean(allWF, 2);

        [h, g] = butter(3, 500/(fs/2), 'high');
        if ~deflect(number)
            center = find(islocalmin(meanWF, 'MaxNumExtrema', 1));
        else
            center = find(islocalmax(meanWF, 'MaxNumExtrema', 1));
        end
        allWF2 = filtfilt(h, g, allWF);
        numWavs = size(allWF2, 2);

        fixedWF = zeros(101, numWavs);
        for i = 1:numWavs
            if ~deflect(number)
                eachCenter = startInd + find(islocalmin(allWF2(startInd:stopInd, i), 'MaxNumExtrema', 1));
            else
                eachCenter = startInd + find(islocalmax(allWF2(startInd:stopInd, i), 'MaxNumExtrema', 1));
            end
            if ~isempty(eachCenter)
                jitter = eachCenter - center;
            else
                jitter = 0;
            end
            fixedWF(:, i) = allWF2(10+jitter:110+jitter, i);
        end

        meanWF2 = mean(fixedWF, 2);
        clusters.meanWF = meanWF2;

        %plot norm waveform
        plotWaveform(xAxisWaveForm, randWindow, meanWF, fontsize, name, f, path);

        %% Get spike parameters
        [wfNorm, trPkTime, trPk70, width, wfMult70] = computeWaveformParams(meanWF, fs, deflect, number);
        clusters.wfNorm = wfNorm;
        clusters.wfRise = trPkTime;
        clusters.wfRise70 = trPk70;
        clusters.wfFWHM = width;
        clusters.wfMult70 = wfMult70;

        %% Get ISIs
        isi_sort = sort(diff(clusters.spikeTimes));
        clusters.isis = isi_sort;
        isi = diff(clusters.spikeTimes);
        ifr = 1 ./ isi;
        isiTimes = (spkTimes(1:end-1) + 0.5 * isi) / fs;

        %% Split pb and spont times
        buffer = 0.5;
        pbTrial = [Motif.start, Motif.stop + buffer];
        pbTimes = cell(size(pbTrial, 1), 1);

        for i = 1:size(pbTrial, 1)
            pbTimes{i} = clusters.spikeTimes(clusters.spikeTimes > pbTrial(i, 1) & clusters.spikeTimes < pbTrial(i, 2));
        end

        % Flatten pbTimes and add NaNs
        pbTimes = vertcat(pbTimes{:});
        pbTimes = [pbTimes; nan];

        % Remove playback spikes from spontaneous times
        spontTimes = clusters.spikeTimes;
        spontTimes(ismember(spontTimes, pbTimes(~isnan(pbTimes)))) = nan;

        % Compute ISIs
        clusters.pbIsis = sort(diff(pbTimes));
        clusters.spontIsis = sort(diff(spontTimes));

        %% Open or make eye stack
        eyeRate = 4;
        if ~contains(name, 'xbeads')
            if isfile(fullfile(path, f, 'stack.mat'))
                load(fullfile(path, f, 'stack.mat'));
            else
                fnameTif = fullfile(path, f, ['eyes_' f '.tif']);
                info = imfinfo(fnameTif);
                numImg = numel(info);
                stack = [];
                for k = 1:numImg
                    currentImage = imread(fnameTif, k, 'Info', info);
                    stack(:, :, k) = currentImage;
                end
                cd(fullfile(path, f));
                save('stack.mat', 'stack');
            end
            lum = mean(mean(stack));
            lum = lum(:);

            %% Check eye vid against eye trigs
            binEyes = diff(eyeTrig);
            frameLoc = find(binEyes > 0.5);
            clear binEyes eyeTrig;
            extra = find(diff(frameLoc) < 10);
            frameLoc(extra) = [];
            if abs(length(frameLoc) - length(lum)) < 2
                disp('checks out')
            else
                disp('nope nope nope')
                disp(string(length(frameLoc)) + ' trigs')
                disp(string(length(lum)) + ' frames')
            end
        else
            switchOff = goodCells{number, 10};
            lum = [ones(switchOff * eyeRate, 1); zeros(round((lengthRecording - switchOff) * eyeRate), 1)];
        end

        %% Get on/off segments from lum
        onThresh = mean(lum);
        eyePts = zeros(length(lum), 1);
        trigs_on = (lum >= onThresh);
        eyePts(trigs_on) = 1;
        lightSeg = find(diff(eyePts ~= 0));
        lightSeg = lightSeg(lightSeg >= floor(startRecFixed / fs * eyeRate) & lightSeg <= floor(stopRecFixed / fs * eyeRate));
        lightSeg = [floor(startRecFixed / fs * eyeRate); lightSeg; floor(stopRecFixed / fs * eyeRate)];

        segs = [];
        for i = 1:length(lightSeg)-1
            segs(i, :) = [lightSeg(i), lightSeg(i+1)];
        end

        onTimes = [];
        offTimes = [];
        for i = 1:size(segs, 1)
            [y, ~] = find(spkTimes / fs >= segs(i, 1) / eyeRate & spkTimes / fs <= segs(i, 2) / eyeRate);
            if eyePts(segs(i, 2))
                onTimes = [onTimes; spkTimes(y) / fs; nan];
            else
                offTimes = [offTimes; spkTimes(y) / fs; nan];
            end
        end

        %extract pb wake ISIs
        pbAwake = ismember(onTimes, pbTimes);
        pbAwakeTimes = onTimes(pbAwake);
        pbAwakeIsis = sort(diff(pbAwakeTimes));
        clusters.onPbIsis = pbAwakeIsis;

        %extract non-pb wake ISIs
        onTimes(pbAwake) = nan;
        onIsis = diff(onTimes);
        on_isi_sort = sort(onIsis);
        clusters.onIsis = on_isi_sort;

        %extract pb sleep ISIs
        pbSleep = ismember(offTimes, pbTimes);
        pbSleepTimes = offTimes(pbSleep);
        pbSleepIsis = sort(diff(pbSleepTimes));
        clusters.offPbIsis = pbSleepIsis;

        %extract non-pb sleep ISIs
        offTimes(pbSleep) = nan;
        offIsis = diff(offTimes);
        off_isi_sort = sort(offIsis);
        clusters.offIsis = off_isi_sort;

        %calculate state modulation index by correlating On and Off ISIs
        edges = -3:.1:1;
        alloni = clusters.onIsis(clusters.onIsis>=1e-3 & clusters.onIsis<=10);
        alloffi = clusters.offIsis(clusters.offIsis>=1e-3 & clusters.offIsis<=10);
        allonih = histcounts(log10(alloni),edges,'Normalization','probability');
        alloffih = histcounts(log10(alloffi),edges,'Normalization','probability');
        c = corrcoef(allonih, alloffih);
        clusters.spontStateCorr = round(c(2,1),2);

        %% Align spikes to On/Off motifs and baseline comparison
        wing = 1 * fs;
        windJux = [-wing/fs (wing+length(motif))/fs];
        %[outStruct, countSp] = alignSpikeTimesJuxta(clusters, windJux, Motif);

        paddedAudio = [zeros(abs(windJux(1))*fs,1); motif; zeros(abs(windJux(1))*fs,1)];
        pbEye = round(Motif.center*eyeRate,0);
        for i = 2:5
            MotifOn.(structElem{i}) = Motif.(structElem{i})(eyePts(pbEye) == 1);
            if i ~= 5
                BeforeOn.(structElem{i}) = MotifOn.(structElem{i})-length(motif)/fs;
            else
                BeforeOn.(structElem{i}) = MotifOn.(structElem{i});
            end
        end
        for i = 2:5
            MotifOff.(structElem{i}) = Motif.(structElem{i})(eyePts(pbEye) == 0);
            if i ~= 5
                BeforeOff.(structElem{i}) = MotifOff.(structElem{i})-length(motif)/fs;
            else
                BeforeOff.(structElem{i}) = MotifOff.(structElem{i});
            end
        end
        [outStructOn, countSpOn] = alignSpikeTimesJuxta(clusters, windJux, MotifOn);
        [outStructOff, countSpOff] = alignSpikeTimesJuxta(clusters, windJux, MotifOff);


        %% Get spontaneous activity snippets
        spontNumTrials = 64;
        segsS = segs / eyeRate;
        fixedSegsOn = [];
        fixedSegsOff = [];

        for i = 1:size(segs,1)
            seg = segsS(i,:);
            if eyePts(segs(i,2)) == 1
                for j = 1:length(pbTrial)
                    while seg(2) > pbTrial(j,2) && seg(1) < pbTrial(j,1)
                        fixedSegsOn = [fixedSegsOn; seg(1), pbTrial(j,1)];
                        seg = [pbTrial(j,2), seg(2)];
                    end
                end
                fixedSegsOn = [fixedSegsOn; seg];
            elseif eyePts(segs(i,2)) == 0
                for j = 1:length(pbTrial)
                    while seg(2) > pbTrial(j,2) && seg(1) < pbTrial(j,1)
                        fixedSegsOff = [fixedSegsOff; seg(1), pbTrial(j,1)];
                        seg = [pbTrial(j,2), seg(2)];
                    end
                end
                fixedSegsOff = [fixedSegsOff; seg];
            end
        end

        % Process fixed segments for spontaneous activity (lights on)
        if ~isempty(fixedSegsOn)
            spontSegsOn = processSpontaneousSegments(fixedSegsOn);
            SpontOn = extractSpontaneousTrials(spontSegsOn, spontNumTrials, number);
        else
            SpontOn = emptySpontaneousStruct();
        end

        % Process fixed segments for spontaneous activity (lights off)
        if ~isempty(fixedSegsOff)
            spontSegsOff = processSpontaneousSegments(fixedSegsOff);
            SpontOff = extractSpontaneousTrials(spontSegsOff, spontNumTrials, number);
        else
            SpontOff = emptySpontaneousStruct();
        end

        % clear spikes
        % clear eyeTrig
        % clear stack

        %% Plot Spontaneous Snippets

        ticklen = 3;
        windSpont = [0 5];
        sigwin10 = .01;
        sigwin50 = .05;
        Trials = [length(SpontOn.start); length(SpontOff.start); length(MotifOn.start); length(MotifOff.start)];


        onTri = length(SpontOn.start);
        offTri = length(SpontOff.start);
        totTri = onTri + offTri + 2;

        if ~isempty(fixedSegsOn)
            clear trialIsi trialSpks
            [spontStructOn, psthOn] = alignSpikeTimesJuxta(clusters, windSpont, SpontOn);

            % Compute burst, CV, and rate for lights on
            for k = 1:length(SpontOn.start)
                trial = spontStructOn.eventTrials == k;
                trialSpks = spontStructOn.eventTimesAligned(trial);
                trialIsi = diff(trialSpks);

                if isempty(trialIsi)
                    burstOn(k) = 0;
                    cvOn(k) = nan;
                    if isempty(trialSpks)
                        onRates(k) = 0;
                    else
                        onRates(k) = length(trialSpks) / windSpont(2);
                    end
                else
                    burstOn(k) = nansum(trialIsi < .01) / length(trialIsi);
                    cvOn(k) = nanstd(trialIsi) / nanmean(trialIsi);
                    onRates(k) = length(trialSpks) / windSpont(2);
                end
            end

            rangeOn = 1:Trials(3);
            if number == 79, rangeOn = 1:11; end
            for i = 2:5
                subsetOn.(structElem{i}) = SpontOn.(structElem{i})(rangeOn);
            end
        else
            spontStructOn.eventTimes = [];
            spontStructOn.eventTimesAligned = [];
            spontStructOn.eventTrials = [];
        end

        if ~isempty(fixedSegsOff)
            clear trialIsi trialSpks
            [spontStructOff, psthOff] = alignSpikeTimesJuxta(clusters, windSpont, SpontOff);

            % Compute burst, CV, and rate for lights off
            for k = 1:length(SpontOff.start)
                trial = spontStructOff.eventTrials == k;
                trialSpks = spontStructOff.eventTimesAligned(trial);
                trialIsi = diff(trialSpks);

                if isempty(trialIsi)
                    burstOff(k) = 0;
                    cvOff(k) = nan;
                    if isempty(trialSpks)
                        offRates(k) = 0;
                    else
                        offRates(k) = length(trialSpks) / windSpont(2);
                    end
                else
                    burstOff(k) = nansum(trialIsi < .01) / length(trialIsi);
                    cvOff(k) = nanstd(trialIsi) / nanmean(trialIsi);
                    offRates(k) = length(trialSpks) / windSpont(2);
                end
            end

            rangeOff = 1:Trials(4);
            for i = 2:5
                subsetOff.(structElem{i}) = SpontOff.(structElem{i})(rangeOff);
            end
        else
            spontStructOff.eventTimes = [];
            spontStructOff.eventTimesAligned = [];
            spontStructOff.eventTrials = [];
        end

        % Make spont plot
        plotSpontaneousActivity(spontStructOn, spontStructOff, SpontOn, SpontOff, windSpont, onTri, offTri, totTri, path, f, name, fontsize);

        %% Compute spontaneous stats

        if ~isempty(fixedSegsOn)
            cvSpontOn = nanmean(cvOn);
            burstSpontOn = nanmean(burstOn);
            rateSpontOn = mean(onRates); % nanmean(psthOn)
            clusters.spontAwakeFR = round(rateSpontOn, 2);
            clusters.spontAwakeBurst = round(burstSpontOn, 2);
            clusters.spontAwakeCV = round(cvSpontOn, 2);
        else
            clusters.spontAwakeFR = nan;
            clusters.spontAwakeBurst = nan;
            clusters.spontAwakeCV = nan;
        end

        if ~isempty(fixedSegsOff)
            cvSpontOff = nanmean(cvOff);
            burstSpontOff = nanmean(burstOff);
            rateSpontOff = mean(offRates); % nanmean(psthOff)
            clusters.spontSleepFR = round(rateSpontOff, 2);
            clusters.spontSleepBurst = round(burstSpontOff, 2);
            clusters.spontSleepCV = round(cvSpontOff, 2);
        else
            clusters.spontSleepFR = nan;
            clusters.spontSleepBurst = nan;
            clusters.spontSleepCV = nan;
        end

        % plot spont stats
        plotMetricsStats(false, fixedSegsOn, fixedSegsOff, clusters, SpontOn, SpontOff, path, f, name, fontsize);

        %% Plot on vs off ISIs
        edges = -3:0.1:1;
        ticks = -3:1;
        labels = 10.^ticks;
        nameLabels = labels * 1000;
        plotISIHistogram(on_isi_sort, off_isi_sort, edges, labels, nameLabels, fontsize, path, f, name);

        %% Plot all pb trials raster
        plotPlaybackActivity(windJux, paddedAudio, fontsize, MotifOn, MotifOff, ...
            outStructOn, outStructOff, ticklen, countSpOn, countSpOff, path, f, name);

        %% Precision index spike alignment
        calcBuffer = .03;

        % precision on
        [onMO,~] = alignSpikeTimesJuxta(clusters, [0+calcBuffer length(motif)/fs+calcBuffer], MotifOn);

        % precision off
        [offMO,~] = alignSpikeTimesJuxta(clusters, [0+calcBuffer length(motif)/fs+calcBuffer], MotifOff);

        % Pre-pb prep & stats
        [beforeOnSpks,~] = alignSpikeTimesJuxta(clusters, [0+calcBuffer length(motif)/fs+calcBuffer], BeforeOn);
        [beforeOffSpks, ~] = alignSpikeTimesJuxta(clusters, [0+calcBuffer length(motif)/fs+calcBuffer], BeforeOff);

        for i = 1:2
            if i == 1
                Event = BeforeOn;
                outStruct = beforeOnSpks;
            elseif i == 2
                Event = BeforeOff;
                outStruct = beforeOffSpks;
            end

            rates = nan(length(Event.center),1);
            trialSpks = [];
            trialIsis = [];
            for k = 1:length(Event.center)
                trial = outStruct.eventTrials == k;
                trialSpks = outStruct.eventTimesAligned(trial);
                rates(k) = length(trialSpks)/(length(motif)/fs);
                trialIsis = diff(trialSpks);
                if isempty(trialSpks)
                    rates(k) = 0;
                end
            end
            meanRate =  round(mean(rates),2);
            if i == 1
                clusters.prePbAwakeFR = meanRate;
                clusters.prePbAwakeFRAll = rates;
            elseif i == 2
                clusters.prePbSleepFR = meanRate;
                clusters.prePbSleepFRAll = rates;
            end
        end

        %% All or nothing match corr

        binbin = 0.005;
        pbHalf = length(motif) / fs / 2;
        buff = 0.5;
        windPb = [-buff, pbHalf * 2 + calcBuffer + buff];
        fullSpec = range(windPb);
        allBins = ceil(fullSpec / binbin);

        conditions = {MotifOn, MotifOff, BeforeOn, BeforeOff};
        clusterFields = {
            'pbAwakeMatch', 'pbAwakeCorr', 'pbAwakeZCorr';
            'pbSleepMatch', 'pbSleepCorr', 'pbSleepZCorr';
            'pbAwakeBeforeMatch', 'pbAwakeBeforeCorr', '';
            'pbSleepBeforeMatch', 'pbSleepBeforeCorr', ''
            };

        for a = 1:4
            Times = conditions{a};

            if ~isempty(Times.center)
                [matchVal, rVal, zVal] = computeBinnedCorrelation(clusters, windPb, Times, binbin, allBins, buff);

                % Assign computed values to clusters structure
                clusters.(clusterFields{a, 1}) = matchVal;
                clusters.(clusterFields{a, 2}) = rVal;
                if ~isempty(clusterFields{a, 3})
                    clusters.(clusterFields{a, 3}) = zVal;
                end
            else
                clusters.(clusterFields{a, 1}) = nan;
                clusters.(clusterFields{a, 2}) = nan;
                if ~isempty(clusterFields{a, 3})
                    clusters.(clusterFields{a, 3}) = nan;
                end
            end
        end

        %% Playback Stats
        onMoRates = computeSpikeMetrics('rate', MotifOn, onMO, motif, fs);
        burstpbAwake = computeSpikeMetrics('burst', MotifOn, onMO, motif, fs);
        cvpbAwake = computeSpikeMetrics('cv', MotifOn, onMO, motif, fs);
        maxIFRAwake = computeSpikeMetrics('maxIFR', MotifOn, onMO, motif, fs);
        meanPostOn = computeSpikeMetrics('postRate', MotifOn, onMO, motif, fs, calcBuffer, clusters, buffer);

        offMoRates = computeSpikeMetrics('rate', MotifOff, offMO, motif, fs);
        burstpbSleep = computeSpikeMetrics('burst', MotifOff, offMO, motif, fs);
        cvpbSleep = computeSpikeMetrics('cv', MotifOff, offMO, motif, fs);
        maxIFRSleep = computeSpikeMetrics('maxIFR', MotifOff, offMO, motif, fs);
        meanPostOff = computeSpikeMetrics('postRate', MotifOff, offMO, motif, fs, calcBuffer, clusters, buffer);

        onMoRatesMean = round(nanmean(onMoRates), 2);
        offMoRatesMean = round(nanmean(offMoRates), 2);
        burstpbAwakeMean = round(nanmean(burstpbAwake), 3);
        burstpbSleepMean = round(nanmean(burstpbSleep), 3);

        % Update clusters structure
        clusters.pbAwakeFR = onMoRatesMean;
        clusters.pbAwakeFRAll = onMoRates;
        clusters.pbAwakeBurst = burstpbAwakeMean;
        clusters.pbAwakeFRPost = meanPostOn;
        clusters.pbAwakeMaxIFR = meanMaxIFRAwake;

        clusters.pbSleepFR = offMoRatesMean;
        clusters.pbSleepFRAll = offMoRates;
        clusters.pbSleepBurst = burstpbSleepMean;
        clusters.pbSleepFRPost = meanPostOff;
        clusters.pbSleepMaxIFR = meanMaxIFRSleep;

        % Plot Playback Stats
        plotMetricsStats(true, fixedSegsOn, fixedSegsOff, clusters, SpontOn, SpontOff, path, f, name, fontsize);

        %% Save other params

        clusters.pbAwakeChangeOld = round((clusters.pbAwakeFR-clusters.spontAwakeFR)./(clusters.pbAwakeFR+clusters.spontAwakeFR),2);
        clusters.pbAwakeChange = round((clusters.pbAwakeFR-clusters.prePbAwakeFR)./(clusters.pbAwakeFR+clusters.prePbAwakeFR),2);
        clusters.pbAwakePostChange = round((clusters.pbAwakeFRPost-clusters.spontAwakeFR)./(clusters.pbAwakeFRPost+clusters.spontAwakeFR),2);
        clusters.pbSleepChangeOld = round((clusters.pbSleepFR-clusters.spontSleepFR)./(clusters.pbSleepFR+clusters.spontSleepFR),2);
        clusters.pbSleepChange = round((clusters.pbSleepFR-clusters.prePbSleepFR)./(clusters.pbSleepFR+clusters.prePbSleepFR),2);
        clusters.pbSleepPostChange = round((clusters.pbSleepFRPost-clusters.spontSleepFR)./(clusters.pbSleepFRPost+clusters.spontSleepFR),2);
        clusters.spontRatioVigilance = round((clusters.spontAwakeFR-clusters.spontSleepFR)./(clusters.spontAwakeFR+clusters.spontSleepFR),2);

        %% Save everything
        cd(path + f)
        save clusters.mat clusters
        counting = counting + 1;
    end

    %Future structures for plotting figures
    AwakePB(counting) = outStructOn;
    SleepPB(counting) = outStructOff;
    AwakeSpont(counting) = spontStructOn;
    AwakeHist{counting,:} = countSpOn;
    SleepSpont(counting) = spontStructOff;
    SleepHist{counting,:} = countSpOff;
    MotifsOnTimes(counting) = MotifOn;
    MotifsOffTimes(counting) = MotifOff;
    SpontsOnTimes(counting) = SpontOn;
    SpontsOffTimes(counting) = SpontOff;
end

%% Save all structures for future plotting

cd(dir)
save AwakePB.mat AwakePB;
save SleepPB.mat SleepPB;
save AwakeSpont.mat AwakeSpont;
save SleepSpont.mat SleepSpont;
save MotifsOnTimes.mat MotifsOnTimes;
save MotifsOffTimes.mat MotifsOffTimes;
save SpontsOnTimes.mat SpontsOnTimes;
save SpontsOffTimes.mat SpontsOffTimes;
save AwakeHist.mat AwakeHist;
save SleepHist.mat SleepHist;
