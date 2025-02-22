%% load everything

clear all; clc; warning off; close all;

fs = 30000;
% wing = 1;
% wind = [-wing wing];
fontsize = 15;
birds = {'C22'; 'C24'; 'C50'; 'C52'; 'C57'};
addpath 'C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Fcns\';
labels = {'sm'; 'pm'; 'pms'};
sigwin10 = .01;
sigwin50 = .05;
calcBuffer = .03;
buffer = .5;
lim = [1e-3 1e1];
edges = -3:.1:1;
spontTrials = 50;
windSpont = [0 5];
binSparse = .005; %same as alignspiketimessi

shitFields = {'ifrBottom','isis','spontIsis','spontOnIsis','spontOffIsis','spontSpks',...
    'spontOnSpks','spontOffSpks','spontAwakePrePbFR','spontAwakePrePbCV',...
    'spontOffEarlyIsis','spontOffDeepIsis','spontAwakePrePbBurst',...
    'nonReplayBurstProp','replayBurstProp','spontAwakeFR','spontAwakeCV',...
    'spontAwakeBurst','spontSleepDeepFR','songPrecInd','pbAwakePrecInd',...
    'pbSleepPrecInd','spontRatioVigilanceE','nonReplayOverReplay','maxChannel'};

dir = 'C:\Users\User\Dropbox (NYU Langone Health)\SiProbe/';

for w = 1:length(birds)
    bird = char(birds(w));
    disp(bird)
    path = [dir bird '/'];
    mkdir([dir 'spontSong/'])
    dataPath = [path 'data/'];
    load([dataPath 'lightSig.mat']);
    load([dataPath 'clusters.mat']);
    load([dataPath 'lenRec.mat']);

    audioPath = [path 'audio/'];
    load([audioPath 'audioLabels.mat']);
    template = audioread([audioPath 'singing/template.wav']);
    motif = audioread([audioPath 'playback/motif.wav']);

    opts = detectImportOptions([path bird '_goodSiUnits.xlsx']);
    opts.DataRange = 'A2';
    opts.VariableNamesRange = 'A1';
    goodCells = readtable([path bird '_goodSiUnits.xlsx'],opts);

    toPrint = (goodCells{:,3});
    startRec = (goodCells{:,4});
    stopRec = (goodCells{:,5});


    %% get good data times

    numUnits = length(clusters);
    clear struct;

    for i = 1:length(shitFields)
        if isfield(clusters,shitFields(i))
            clusters = rmfield(clusters,shitFields(i));
        end
    end
    counter = 1;

    for n = 1:numUnits %[1,3,4,5,6,7,14,20,23,24,26,35,39,40,45,47,48,58,60]
        close all;

        if toPrint(n) == 1

            disp(n)

            startRecFixed = startRec(n)*fs;
            startRecFixed(startRecFixed == 0) = 1;
            if strcmp(stopRec{n},'end')
                stopRecFixed = lenRec;
            else
                stopRecFixed = str2double(stopRec{n})*fs;
            end
            goodDataTimes = [startRecFixed/fs stopRecFixed/fs];
            spkTimes = clusters(n).spikeTimes(clusters(n).spikeTimes >= startRecFixed/fs & clusters(n).spikeTimes <= stopRecFixed/fs);
            isis = diff(spkTimes);
            fixed = isis > .0015; %1.5 ms cutoff
            spkTimes = spkTimes(fixed);
            isis = isis(fixed);
            %clusters(n).isis = isis;
            %  if ~isfield(clusters,'ogTimes')
            % if isempty(clusters(n).ogTimes)
            clusters(n).ogTimes = clusters(n).spikeTimes;
            %    end
            %    end
            clusters(n).spikeTimes = spkTimes;


            %% get spike times for events

            motifRange = [ceil(startRecFixed/fs) floor(stopRecFixed/fs)];
            clear outStruct Event EventOG
            for i = 1:length(labels)
                if strcmp(labels{i},'sm')
                    audio = template;
                else
                    audio = motif;
                end
                EventOG = getEventTimes(audioPath, fs, audio, audioLabels, labels{i});
                %   times = [dataPath labels{i} 'Times.mat'];
                takeMotif = EventOG.start>=motifRange(1) & EventOG.stop<=motifRange(2);
                structElem = fieldnames(EventOG);
                for j = 2:5
                    Event.(structElem{j}) = EventOG.(structElem{j})(takeMotif);
                end
                assignin('base',[labels{i} 'Times'],Event);

                if isempty(Event.start)
                    for j = 2:5
                        Event.(structElem{j}) = nan;
                    end
                end
                assignin('base',[labels{i} 'Times'],Event);

                if strcmp(labels{i},'sm')
                    [structs(i).outStruct, ~] = alignSpikeTimesSi(clusters(n), [0 length(audio)/fs], Event, Event.center);
                    rng(n);
                    whichSongTrials = 1:length(Event.center);% randperm(length(Event.center),7);
                    for c = 2:5
                        EventRand.(structElem{c}) = Event.(structElem{c})(whichSongTrials);
                    end
                    [~, psth] = alignSpikeTimesSi(clusters(n), [0 length(audio)/fs], EventRand, EventRand.center);
                    assignin('base',[labels{i} 'Psth'],psth);
                else
                    [structs(i).outStruct, ~] = alignSpikeTimesSi(clusters(n), [calcBuffer length(audio)/fs+calcBuffer], Event, Event.start);
                end
                assignin('base',[labels{i} 'Spks'],structs(i).outStruct);

                rates = nan(length(Event.center),1);
                maxIFR = nan(length(Event.center),1);
                bursty = nan(length(Event.center),1);
                trialSpks = [];
                trialIsis = [];
                whenIFR = [];
                maxTime = nan(length(Event.center),1);

                for k = 1:length(Event.center)
                    trial = structs(i).outStruct.eventTrials == k;
                    trialSpks = structs(i).outStruct.eventTimesAligned(trial);
                    rates(k) = length(trialSpks)/(length(audio)/fs);
                    trialIsis = diff(trialSpks);
                    if ~isempty(trialIsis)
                        ifr = 1./trialIsis;
                        maxIFR(k) = max(ifr);
                        whenIFR = find(ifr == maxIFR(k));
%                         if length(whenIFR>1)
%                             disp('fucked')
%                         end
                        maxTime(k) = mean(trialSpks(whenIFR:whenIFR+1));
                        bursty(k) = sum(trialIsis < .01)/length(trialIsis);
                        %     elseif length(trialSpks)==1
                        %     maxIFR(k) = 1;
                    elseif isempty(trialIsis)
                        maxIFR(k) = nan;
                        bursty(k) = 0;
                        maxTime(k) = 0;
                        if isempty(trialSpks)
                            rates(k) = 0;
                        else
                            rates(k) = length(trialSpks)/(length(audio)/fs);
                        end
                    end
                end
                
%                 ifrMat = nan(length(maxTime),length(maxTime));
%                 for m = 1:length(maxTime)
%                     t1 = maxTime(m);
%                     for v = 1:length(maxTime)
%                         t2 = maxTime(v);
%                         if abs(t2-t1) <= .05
%                             ifrMat(m,v) = 2;
%                         else
%                             ifrMat(m,v) = 1;
%                         end
%                     end
%                 end
%                 
%                 ifrBottom = tril(ifrMat,1);
%             %    imagesc(ifrBottom);
%                 ifrBottom(ifrBottom == triu(ifrBottom)) = nan;
%                 ifrBottom = ifrBottom - 1;
%                 ifrScore = nanmean(ifrBottom(:));
%                 if isempty(ifrScore)
%                     disp('empty')
%                 end
%                 if i == 2
%                     clusters(n).pbAwakeIFRScore = ifrScore;
%                     clusters(n).ifrBottom = ifrBottom;
%                 elseif i == 3
%                     clusters(n).pbSleepIFRScore = ifrScore;
%                     clusters(n).ifrBottom = ifrBottom;
%                 end

                meanRate =  round(mean(rates),2);
                burstyMean = round(nanmean(bursty),3);
                burstyMean(isnan(burstyMean)) = 0;
                meanMaxIFR = round(nanmean(maxIFR),2);
                
                assignin('base',[labels{i} 'FRAll'], rates);
                assignin('base',[labels{i} 'FR'], meanRate);
                assignin('base',[labels{i} 'Burst'], burstyMean);
                assignin('base',[labels{i} 'MaxIFR'], meanMaxIFR);
            end
            
            
            %% prepare before Motif stats

        for i = 2:5
            if i ~= 5
                beforePmTimes.(structElem{i}) = pmTimes.(structElem{i})-length(motif)/fs;
                beforePmsTimes.(structElem{i}) = pmsTimes.(structElem{i})-length(motif)/fs;
            else
                beforePmTimes.(structElem{i}) = pmTimes.(structElem{i});
                beforePmsTimes.(structElem{i}) = pmsTimes.(structElem{i});
            end
        end

        [beforePmSpks,~] = alignSpikeTimesSi(clusters(n), [0+calcBuffer length(motif)/fs+calcBuffer], beforePmTimes,beforePmTimes.start);
        [beforePmsSpks, ~] = alignSpikeTimesSi(clusters(n), [0+calcBuffer length(motif)/fs+calcBuffer], beforePmsTimes,beforePmTimes.start);
        
        for i = 1:2
            if i == 1
                Event = beforePmTimes;
                outStruct = beforePmSpks;
            elseif i == 2
                Event = beforePmsTimes;
                outStruct = beforePmsSpks;
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
                clusters(n).prePbAwakeFR = meanRate;
                clusters(n).prePbAwakeFRAll = rates;
            elseif i == 2
                clusters(n).prePbSleepFR = meanRate; 
                clusters(n).prePbSleepFRAll = rates;
            end
        end
         
            Motif = [smTimes, pmTimes, pmsTimes, beforePmTimes, beforePmsTimes];


            %% binary song precision

            songHalf = length(template)/fs/2;
            buff = .5;
            binbin = .005;
            windSong = [-(buff + songHalf) (buff + songHalf)];
            fullSpec = windSong(2)*2;
            allBins = ceil(fullSpec/binbin);
            series = nan(length(smTimes.center),allBins);
            outStruct = alignSpikeTimesSi(clusters(n), windSong, smTimes, smTimes.center);
            poss = ceil(-allBins/2):ceil(allBins/2)-1;
            for k = 1:length(smTimes.center)
                trial = outStruct.eventTrials == k;
                trialSpks = outStruct.eventTimesAligned(trial);
                whichBin = ceil(trialSpks./binbin);
                time = zeros(length(poss),1);
                time(ismember(poss,whichBin)) = 1;
                series(k,:) = time;
            end

            songHalfBins = ceil(songHalf/binbin);
            midSeries = find(poss==0);
            songBins = midSeries-songHalfBins:midSeries+songHalfBins;
            series = series(:,songBins);

%                         %plot binarization
%                         close all;
%                         figure('Visible','off','WindowState','maximized');
%                         for i = 1:size(series,1)
%                             plot(series(i,:)-i);
%                             hold on;
%                         end
%             
%                         hold on

            r = nan(size(series,1),size(series,1));
            t = nan(size(series,1),size(series,1));
            for i = 1:size(series,1)
                for j = 1:size(series,1)
                    score = corrcoef(series(i,:),series(j,:));
                    r(i,j) = score(1,2);
                    if isnan(r(i,j))
                        r(i,j) = 0;
                    end
                end
            end

            rBottom = tril(r,1);
            rBottom(rBottom == triu(rBottom)) = nan;
            rVal = nanmean(rBottom(:));

            z = atanh(r);
            z(isinf(z)) = atanh(.9999);
            zBottom = tril(z,1);
            zBottom(zBottom == triu(zBottom)) = nan;
            zVal = nanmean(zBottom(:));
            clusters(n).songCorr = rVal;
            clusters(n).songRtoZ = zVal;

            % all or nothing corr
            match = nan(size(series,1),size(series,1));
            for i = 1:size(series,1)
                for j = 1:size(series,1)
                    comp1 = series(i,:);
                    comp2 = series(j,:);
                     for k = 1:length(comp1)
                            if comp1(k) == comp2(k)
                              if comp1(k) == 1
                                    same(k) = 1;
                               elseif comp1(k) == 0
                                   same(k) = nan;
                               end
                            elseif comp1(k) ~= comp2(k)
                                same(k) = 0;
                            end
                            match(i,j) = nanmean(same);
                        end
                end
            end

            matchBottom = tril(match,1);
            matchBottom(matchBottom == triu(matchBottom)) = nan;
            matchVal = nanmean(matchBottom(:));
            clusters(n).songMatch = round(matchVal,4);


            %% binary pb precision
            
            binbin = .005;
            for a = 2:3
                pbHalf = length(motif)/fs/2;
                windPb = [-buff pbHalf*2+calcBuffer+buff];
                fullSpec = range(windPb);
                allBins = ceil(fullSpec/binbin);
                series = nan(length(Motif(a).center),allBins);
                outStruct = alignSpikeTimesSi(clusters(n), windPb, Motif(a), Motif(a).start);
                poss = ceil(windPb(1)/binbin):ceil(windPb(2)/binbin)-1;
                for k = 1:length(Motif(a).center)
                    trial = outStruct.eventTrials == k;
                    trialSpks = outStruct.eventTimesAligned(trial);
                    whichBin = ceil(trialSpks./binbin);
                    time = zeros(length(poss),1);
                    time(ismember(poss,whichBin)) = 1;
                    series(k,:) = time;
                end
           
                pbFlankBins = ceil(buff/binbin);
                pbBins = pbFlankBins+1:allBins-pbFlankBins;
                series = series(:,pbBins);

            r = nan(size(series,1),size(series,1));
            t = nan(size(series,1),size(series,1));
            for i = 1:size(series,1)
                for j = 1:size(series,1)
                    score = corrcoef(series(i,:),series(j,:));
                    r(i,j) = score(1,2);
                    if isnan(r(i,j))
                        r(i,j) = 0;
                    end
                end
            end

            rBottom = tril(r,1);
            rBottom(rBottom == triu(rBottom)) = nan;
            rVal = nanmean(rBottom(:));

            z = atanh(r);
            z(isinf(z)) = atanh(.9999);
            zBottom = tril(z,1);
            zBottom(zBottom == triu(zBottom)) = nan;
            zVal = nanmean(zBottom(:));

                match = nan(size(series,1),size(series,1));
                for i = 1:size(series,1)
                    for j = 1:size(series,1)
                        comp1 = series(i,:);
                        comp2 = series(j,:);
                        for k = 1:length(comp1)
                            if comp1(k) == comp2(k)
                              if comp1(k) == 1
                                    same(k) = 1;
                               elseif comp1(k) == 0
                                   same(k) = nan;
                               end
                            elseif comp1(k) ~= comp2(k)
                                same(k) = 0;
                            end
                            match(i,j) = nanmean(same);
                        end
                    end
                end

                matchBottom = tril(match,1);
                matchBottom(matchBottom == triu(matchBottom)) = nan;
                matchVal = round(nanmean(matchBottom(:)),4);
                if a == 2
                    clusters(n).pbAwakeMatch = matchVal;
                    clusters(n).pbAwakeCorr = rVal;
                    clusters(n).pbAwakeZCorr = zVal;
                elseif a == 3
                    clusters(n).pbSleepMatch = matchVal;
                    clusters(n).pbSleepCorr = rVal;
                    clusters(n).pbSleepZCorr = zVal;
                elseif a == 4
                    clusters(n).pbAwakeBeforeMatch = matchVal;
                    clusters(n).pbAwakeBeforeCorr = rVal;

                elseif a == 5
                    clusters(n).pbSleepBeforeMatch = matchVal;
                    clusters(n).pbSleepBeforeCorr = rVal;
                end
            end

            %             title([string(birds(w)) + "_" + string(clusters(n).clusterID) + "   r = " + round(rVal,3) + '; z = ' + round(zVal,3) + '; match = ' + round(matchVal,3)])
            %
            %             cd('C:\Users\User\Dropbox (NYU Langone Health)\SiProbe\precision\')
            %             print([string(birds(w)) + "_" + string(clusters(n).clusterID) + ' adjusted'], '-dpdf','-fillpage');


            %% now do everything

            %wf params
            wav = clusters(n).meanWF;

            meanSubWav = wav-mean(wav);
            maxAmp = round(max(abs(meanSubWav)));
            normWav = meanSubWav./maxAmp;
            %normWavShort = normWav(40:94);
            clusters(n).wfNorm = normWav;

            factor = 1000;
            normWavInterp = round(interp(normWav,factor),5);

            trVal = min(normWavInterp);
            trPos = find(islocalmin(normWavInterp,'MaxNumExtrema',1));
            pkVal = max(normWavInterp(trPos:end));
            pkPos = (trPos-1)+find(islocalmax(normWavInterp(trPos:end),'MaxNumExtrema',1));
            trPkTime = (pkPos-trPos)/fs*1000/factor;
            clusters(n).wfRise = trPkTime;

            pkVal70 = trVal+.7*range([pkVal trVal]);
            closeEnuf = find(abs(normWavInterp - pkVal70) < .0015);
            inBet = closeEnuf(closeEnuf < pkPos & closeEnuf >  trPos);
            pkPos70 = inBet(islocalmin(abs(normWavInterp(inBet) - pkVal70),'MaxNumExtrema',1));
            trPk70 = (pkPos70-trPos)/fs*1000/factor;
            clusters(n).wfRise70 = trPk70;

            pkVal50 = trVal+.5*range([pkVal trVal]);
            closeEnuf = find(abs(normWavInterp - pkVal50) < .0015);
            beforePk =  closeEnuf(closeEnuf < pkPos);
            pkPos50 = beforePk(islocalmin(abs(normWavInterp(beforePk) - pkVal50),'MaxNumExtrema',2));
            width = diff(pkPos50)/fs*1000/factor;
            if isempty(width)
                width = nan;
            end
            clusters(n).wfFWHM = width;

            %             figure(w)
            %             counter = counter + 1;
            %             subplot(4,16,counter)
            %             plot([normWav(2:end)],'k')
            %             hold on
            %             plot(pkPos70/factor,pkVal70,'b*')
            %             hold on
            %             line([trPos pkPos70]/factor,[trVal pkVal70],'Color','b')
            %             hold on
            %             plot(pkPos/factor,pkVal,'r*')
            %             hold on
            %             plot(trPos/factor,trVal,'r*')
            %             hold on
            %             line([pkPos50 pkPos50]/factor,[pkVal50 pkVal50],'Color','b')
            %             title(clusters(n).clusterID)
            %             set(gca,'color','none')
            %             text(30,-.8,string(round(trPk70,3)))
            %             text(30,-.9,string(round(width,3)))

            % split pb and spont and song
            spontSpks = clusters(n).spikeTimes;
            audSegs = [];
            audTrials = [];
            audLabels = {'sm'; 'pm';'pms'};

            for i = 1:3
                if i == 1
                    %so that buffer doesnt overlap with next motif
                    audTrial = [Motif(i).start Motif(i).stop];
                else
                    audTrial = [Motif(i).start Motif(i).stop+buffer];
                end
                audTrials = [audTrials; audTrial];
                audSpks = [];
                for j = 1:length(audTrial)
                    audSpks = [audSpks; spontSpks(spontSpks >= audTrial(j,1) & spontSpks <= audTrial(j,2))];
                end
                audSpks = sort(audSpks);
                Motif(i).audSpks = audSpks;
                isi = sort(diff(audSpks));
                Motif(i).audIsis = isi(isi <= lim(2) & isi >= lim(1));
                [is loc] = ismember(audSpks, spontSpks);
                spontSpks(loc) = nan;

                % for precision
                struct = structs(i).outStruct;
                audLat = [];
                if ~isempty(struct.eventTimes)

                    for x = 1%:length(struct.eventTimesAligned)
                        t = struct.eventTrials(x);
                        compare = struct.eventTimesAligned(struct.eventTrials~=t);
                        lats = compare-struct.eventTimesAligned(x);
                        audLat = [audLat; lats'];
                        clear lats compare
                    end
                else
                    audLat = 0;
                end

                audProp10(i) = length(audLat(audLat >= -sigwin10 & audLat <= sigwin10))/(length(audLat));
                audProp50(i) = length(audLat(audLat >= -sigwin50 & audLat <= sigwin50))/(length(audLat));
                assignin('base',[audLabels{i} 'Lat'],audLat);
                % for spont segs
                audSegs = [audSegs; Motif(i).start(1) Motif(i).stop(end)+buffer];
            end


            %% deal with c22 and c24 pmb and pmbs

            if strcmp(bird, 'C22') | strcmp(bird, 'C24')
                badLabels = {'pmb'; 'pmbs'};
                bout = audioread([path '\audio\playback\bout.wav']);
                for a = 1:length(badLabels)
                    clear junk;
                    junk = getEventTimes(audioPath, fs, bout, audioLabels, badLabels(a));
                    takeJunk = junk.start>=motifRange(1) & junk.stop<=motifRange(2);
                    for j = 2:5
                        junk.(structElem{j}) = junk.(structElem{j})(takeJunk);
                    end
                    %[junkStruct, ~] = alignSpikeTimesSi(clusters(n), [0 length(bout)/fs+calcBuffer], junk, junk.start);
                    junkTrial = [junk.start junk.stop+buffer];
                    junkSpks = [];
                    for j = 1:length(junkTrial)
                        junkSpks = [junkSpks; spontSpks(spontSpks >= junkTrial(j,1) & spontSpks <= junkTrial(j,2))];
                    end
                    [is2 loc2] = ismember(junkSpks,spontSpks);
                    spontSpks(loc2) = nan;
                end
            end


            %%
            % spont snips
            audSegs = sort(audSegs);
            audTrials = sort(audTrials);
            nonAudTrials = [];
            for i = 1:length(audTrials)-1
                nonAudTrials = [nonAudTrials; audTrials(i,2)+buffer audTrials(i+1,1)-buffer];
            end
            segs = nonAudTrials;
            segs = sort(segs);
            segs = [startRecFixed/fs audTrials(1,1); segs];
            segs = [segs; segs(end,2) stopRecFixed/fs];

            splitInd = find(segs(:,2)>lightSeg & segs(:,1)<lightSeg);
            segs = [segs(1:splitInd-1, :); segs(splitInd,1) lightSeg; lightSeg segs(splitInd,2); segs(splitInd+1:end,:)];

            necLen = 10;
            segLengths = segs(:,2)-segs(:,1);
            longSegs = segLengths >= necLen;
            longEnuf = segLengths >= 5 & segLengths < necLen;
            spontSegs = segs(longEnuf,:);
            for i = 1:length(segLengths)
                if longSegs(i)
                    shortSegs = [];
                    left = floor(segLengths(i)/necLen);
                    shortSegs(:,1) = [segs(i,1):necLen:segs(i,2)];
                    shortSegs(:,2) = shortSegs(:,1)+necLen;
                    shortSegs(left+1,:) = [shortSegs(left,2) segs(i,2)];
                    spontSegs = [spontSegs; shortSegs];
                end
            end
            spontSegs = sort(spontSegs(spontSegs(:,2)-spontSegs(:,1)>=5,:));

            split = find(spontSegs(:,1) == lightSeg);
            segsOn  = spontSegs(1:split-1,:);
            segsOff = spontSegs(split:end,:);
            deepCut = lightSeg + 600; %THIS WAS 300 UNTIL 10/19/23
            whichEarly = segsOff(:,2)<= deepCut;
            segsOffEarly = segsOff(whichEarly,:);
            segsOffDeep = segsOff(~whichEarly,:);

            rng(n); 
            whichSpont = randperm(length(segsOn),min(length(segsOn),spontTrials)); 
            count = 1; 
            for i = whichSpont
                a = ceil(segsOn(i,1)); 
                b = ceil(segsOn (i,2))-5; 
                snip = a + (b-a).*rand(1,1);
                SpontAwake(n).start(count) = snip;
                SpontAwake(n).stop(count) = snip+5;
                SpontAwake(n).center(count) = snip+2.5; 
                count = count + 1;
            end 
            SpontAwake(n).warp = ones(length(SpontAwake(n).start),1);

            rng(n);
            whichSpont = randperm(length(segsOffDeep),min(length(segsOffDeep),spontTrials));
            count = 1;
            for i = whichSpont
                a = ceil(segsOffDeep(i,1));
                b = ceil(segsOffDeep(i,2))-5;
                snip = a + (b-a).*rand(1,1);
                SpontSleepDeep(n).start(count) = snip;
                SpontSleepDeep(n).stop(count) = snip+5;
                SpontSleepDeep(n).center(count) = snip+2.5;
                count = count + 1;
            end
            SpontSleepDeep(n).warp = ones(length(SpontSleepDeep(n).start),1);

            rng(n);
            whichSpont = randperm(length(segsOffEarly),min(length(segsOffEarly),spontTrials));
            count = 1;
            for i = whichSpont
                a = ceil(segsOffEarly(i,1));
                b = ceil(segsOffEarly (i,2))-5;
                snip = a + (b-a).*rand(1,1);
                SpontSleepEarly(n).start(count) = snip;
                SpontSleepEarly(n).stop(count) = snip+5;
                SpontSleepEarly(n).center(count) = snip+2.5;
                count = count + 1;
            end
            SpontSleepEarly(n).warp = ones(length(SpontSleepEarly(n).start),1);

            %  split on and off spikes/isis
            spontIsis = diff(spontSpks);
            spont_isi_sort = sort(spontIsis);
            %  clusters(n).spontIsis = spont_isi_sort;
            %   clusters(n).spontSpks = spontSpks;

            spontOnSpks = spontSpks(spontSpks<=lightSeg);
            spontOnIsis = sort(diff(spontOnSpks));
            spontOnIsis = spontOnIsis (spontOnIsis  <= lim(2) & spontOnIsis >= lim(1));
            %   clusters(n).spontOnSpks = spontOnSpks;
            %    clusters(n).spontOnIsis = spontOnIsis;

            spontOffSpks = spontSpks(spontSpks>lightSeg);
            spontOffIsis = sort(diff(spontOffSpks));
            spontOffIsis = spontOffIsis (spontOffIsis  <= lim(2) & spontOffIsis >= lim(1));
            %       clusters(n).spontOffSpks = spontOffSpks;
            %      clusters(n).spontOffIsis = spontOffIsis;

            smIsis = Motif(1).audIsis;
            clusters(n).smIsis = smIsis;
            pmIsis = Motif(2).audIsis;
            clusters(n).pmIsis = pmIsis;
            pmsIsis = Motif(3).audIsis;
            clusters(n).pmsIsis = pmsIsis;

            % sparseness during singing
            %    points = length(smPsth);
            templateLength = ceil(length(template)/fs/binSparse);
            %    smPsthWing = floor((points-templateLength)/2);
            %    smPsthInterval = smPsthWing+1:points-smPsthWing-1;
            %    smZeroInd = find(smPsth(smPsthInterval) == 0);
            smZeroInd = find(smPsth == 0);
            sparse = round(length(smZeroInd)/templateLength,3);
           % sparse(isempty(sparse)) = 0;
            clusters(n).songSparse = sparse;


            %% juxta

            %SpontOff = SpontOffEarly;
            rng(n)
            whichSpont = randperm(spontTrials);

%             for i = 2:5
%                 if i == 5
%                     SpontAwake.(structElem{i}) = pmTimes.(structElem{i})(whichSpont);
%                 else
%                     SpontAwake.(structElem{i})= (pmTimes.(structElem{i})(whichSpont))-6;
%                 end
%                                 %for a call during this epoch in c52; c50 all good
%                                 if w == 4
%                                     SpontAwake.(structElem{i})(13) = [];
%                                 end
%             end

            Trials = [min(length(segsOn),spontTrials); min(length(segsOffEarly),spontTrials); min(length(segsOffDeep),spontTrials)];
            spontLabels = {'spontAwake'; 'spontSleepEarly'; 'spontSleepDeep'};
            spontTimes = [SpontAwake(n); SpontSleepEarly(n); SpontSleepDeep(n)];

            for s = 1:length(spontLabels)
                [Sevent] = alignSpikeTimesSi(clusters(n), windSpont, spontTimes(s), spontTimes(s).start);
                assignin('base',[spontLabels{s} 'Times'],Sevent);

                trialSpks = [];
                trialIsi = [];
                numTimes = length(spontTimes(s).start);
                bursty = nan(numTimes,1);
                rates =  nan(numTimes,1);
                cv =  nan(numTimes,1);
                for k = 1:numTimes

                    if spontTimes(s).start(k)>= goodDataTimes(1) && spontTimes(s).stop(k)<= goodDataTimes(2)
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
                            rates(k) = length(trialSpks)/range(windSpont);
                        end
                    elseif isempty(trialIsi) && ~within
                        bursty(k) = nan;
                        cv(k) = nan;
                        rates(k) = nan;
                    else
                        bursty(k) = nansum(trialIsi < .01)./length(trialIsi);
                        cv(k) = nanstd(trialIsi)./nanmean(trialIsi);
                        rates(k) = length(trialSpks)/range(windSpont);
                    end
                end

                cvMean = nanmean(cv);
                burstyMean = nanmean(bursty);
                burstyMean(isnan(burstyMean)) = 0;
                ratesMean = nanmean(rates);
                assignin('base',[spontLabels{s} 'FR'],ratesMean);
                assignin('base',[spontLabels{s} 'CV'],cvMean);
                assignin('base',[spontLabels{s} 'Burst'],burstyMean);

                %spontaneous precision
                numMotifs = min(Trials(s),length(Motif(s).start));
                rangeT = 1:Trials(s);
                for i = 2:5
                    subset.(structElem{i}) = spontTimes(s).(structElem{i})(rangeT);
                end

%                 [EventPrec] = alignSpikeTimesSi(clusters(n), [0 length(motif)/fs], subset, subset.start);
%                 spontLat = [];
%                 if ~isempty(EventPrec.eventTimes)
% 
%                     for i = 1%:length(EventPrec.eventTimesAligned)
%                         t = EventPrec.eventTrials(i);
%                         compare = EventPrec.eventTimesAligned(EventPrec.eventTrials~=t);
%                         precSpont = compare-EventPrec.eventTimesAligned(i);
%                         spontLat = [spontLat; precSpont'];
%                         clear precSpont compare;
%                     end
%                 else
%                     spontLat = 0;
%                 end
% 
%                 propSpont1 = sum(spontLat >= -sigwin10 & spontLat <= sigwin10)/length(spontLat);
%                 propSpont2 = sum(spontLat >= -sigwin50 & spontLat <= sigwin50)/length(spontLat);
% 
%                 assignin('base',[spontLabels{s} 'Lat'], spontLat);
%                 assignin('base',[spontLabels{s} 'Prec10'], propSpont1);
%                 assignin('base',[spontLabels{s} 'Prec50'], propSpont2);
                assignin('base',[spontLabels{s} 'Spks'], Sevent);

                clear Sevent
                clear EventPrec
            end
            %             figure;
            %             plotEventRasters(spontAwakeSpks,SpontAwake,[0 5])

            %  clusters(n).spontAwakePrec = round(spontAwakePrec,4);
            %  clusters(n).spontSleepPrec = round(spontSleepDeepPrec,4);
            %  clusters(n).songPrec = round(audProp(1),4);
            %  clusters(n).pbAwakePrec = round(audProp(2),4);
            %  clusters(n).pbSleepPrec = round(audProp(3),4);

            clusters(n).spontAwakeFR = round(spontAwakeFR,2);
            clusters(n).spontAwakeCV = round(spontAwakeCV,2);
            clusters(n).spontAwakeBurst = round(spontAwakeBurst,3);
%             clusters(n).songPrecInd10 = round((audProp10(1)-spontAwakePrec10)/(audProp10(1)+spontAwakePrec10),4);
%             clusters(n).songPrecInd50 = round((audProp50(1)-spontAwakePrec50)/(audProp50(1)+spontAwakePrec50),4);
%             clusters(n).pbAwakePrecInd10 = round((audProp10(2)-spontAwakePrec10)/(audProp10(2)+spontAwakePrec10),4);
%             clusters(n).pbAwakePrecInd50 = round((audProp50(2)-spontAwakePrec50)/(audProp50(2)+spontAwakePrec50),4);
%             clusters(n).pbSleepPrecInd10 = round((audProp10(3)-spontSleepDeepPrec10)/(audProp10(3)+spontSleepDeepPrec10),4);
%             clusters(n).pbSleepPrecInd50 = round((audProp50(3)-spontSleepDeepPrec50)/(audProp50(3)+spontSleepDeepPrec50),4);

            % clusters(n).spontAwakeFR = round(spontAwakeFR,2);
            % clusters(n).spontAwakeBurst = round(spontAwakeBurstProp,3);
            % clusters(n).spontAwakeCV = round(spontAwakeCV,2);

            clusters(n).spontSleepDeepFR = round(spontSleepDeepFR,2);
            clusters(n).spontSleepDeepBurst = round(spontSleepDeepBurst,2);
            clusters(n).spontSleepDeepCV = round(spontSleepDeepCV,2);
            clusters(n).spontSleepEarlyFR = round(spontSleepEarlyFR,2);
            clusters(n).spontSleepEarlyBurst = round(spontSleepEarlyBurst,3);
            clusters(n).spontSleepEarlyCV = round(spontSleepEarlyCV,2);

            % %             %get awake spont pre-pm
            % %             [prePmSpont] = alignSpikeTimesSi(clusters(n), [-6 -1], pmTimes, pmTimes.start);
            % %             trialSpks = [];
            % %             trialIsi = [];
            % %             numTimes = length(pmTimes.center);
            % %             rates = nan(numTimes,1);
            % %             bursty = nan(numTimes,1);
            % %             cv =  nan(numTimes,1);
            % %             for k = 1:numTimes
            % %                 trial = prePmSpont.eventTrials == k;
            % %                 trialSpks = prePmSpont.eventTimesAligned(trial);
            % %                 trialIsi = diff(trialSpks);
            % %
            % %                 if isempty(trialIsi)
            % %                     rates(k) = 0;
            % %                     bursty(k) = 0;
            % %                     cv(k) = nan;
            % %                 %for a call during this epoch in c52
            % %                 elseif k == 13 && w == 4
            % %                     rates(k) = nan;
            % %                 else
            % %                     rates(k) = length(trialSpks)/windSpont(2);
            % %                     bursty(k) = nansum(trialIsi < .01)./length(trialIsi);
            % %                     cv(k) = nanstd(trialIsi)./nanmean(trialIsi);
            % %                 end
            % %             end
            % %             allPrePm(n) = prePmSpont;
            % %             cvMean = nanmean(cv);
            % %             burstyMean = nanmean(bursty);
            % %             burstyMean(isnan(burstyMean)) = 0;
            % %             prePmRatesMean = mean(rates);
            % %


            %get post sm frs
            bouts = find(diff(smTimes.stop)>2);
            smTimesLast.start = smTimes.start(bouts);
            smTimesLast.stop = smTimes.stop(bouts);
            smTimesLast.warp = ones(length(bouts),1);
            [smPost,psthSmPost] = alignSpikeTimesSi(clusters(n), [length(template)/fs buffer+length(template)/fs], smTimesLast, smTimesLast.start);
            smPostRates = nan(length(smTimesLast.start),1);
            trialSpks = [];
            for k = 1:length(smTimesLast.start)
                trial = smPost.eventTrials == k;
                trialSpks = smPost.eventTimesAligned(trial);
                if isempty(trialSpks)
                    smPostRates(k) = 0;
                else
                    smPostRates(k) = length(trialSpks)/(length(template)/fs);
                end
            end
            smPostRatesMean = round(mean(smPostRates),2);


            %get post pb frs
            [pmPost,psthPmPost] = alignSpikeTimesSi(clusters(n), [length(motif)/fs+calcBuffer buffer+length(motif)/fs+calcBuffer], pmTimes, pmTimes.start);
            pmPostRates = nan(length(pmTimes.center),1);
            trialSpks = [];
            for k = 1:length(pmTimes.center)
                trial = pmPost.eventTrials == k;
                trialSpks = pmPost.eventTimesAligned(trial);
                if isempty(trialSpks)
                    pmPostRates(k) = 0;
                else
                    pmPostRates(k) = length(trialSpks)/(length(motif)/fs);
                end
            end
            pmPostRatesMean = round(mean(pmPostRates),2);

            [pmsPost,psthPmsPost] = alignSpikeTimesSi(clusters(n), [length(motif)/fs+calcBuffer buffer+length(motif)/fs+calcBuffer], pmsTimes, pmsTimes.start);
            trialSpks = [];
            pmsPostRates = nan(length(pmsTimes.center),1);
            for k = 1:length(pmsTimes.center)
                trial = pmsPost.eventTrials == k;
                trialSpks = pmsPost.eventTimesAligned(trial);
                if isempty(trialSpks)
                    pmsPostRates(k) = 0;
                else
                    pmsPostRates(k) = length(trialSpks)/(length(motif)/fs);
                end
            end
            pmsPostRatesMean = round(mean(pmsPostRates),2);

            clusters(n).songFR = smFR;
            clusters(n).songFRAll = smFRAll;
            clusters(n).songBurst = smBurst;
            clusters(n).songMaxIFR = smMaxIFR;
            clusters(n).songFRPost = smPostRatesMean;
            clusters(n).pbAwakeFR = pmFR;
            clusters(n).pbAwakeFRAll = pmFRAll;
            clusters(n).pbAwakeFRPost = pmPostRatesMean;
            clusters(n).pbAwakeBurst = pmBurst;
            clusters(n).pbAwakeMaxIFR = pmMaxIFR;
            clusters(n).pbSleepFR = pmsFR;
            clusters(n).pbSleepFRAll = pmsFRAll;
            clusters(n).pbSleepFRPost = pmsPostRatesMean;
            clusters(n).pbSleepBurst = pmsBurst;
            clusters(n).pbSleepMaxIFR = pmsMaxIFR;


            %% next is combinatorial stats for pb ratios and statecorr ratios given isis

            spontAwakeIsis = diff(spontOnSpks);
            offEisi = diff(spontOffSpks(spontOffSpks<=deepCut));
            offDisi = diff(spontOffSpks(spontOffSpks>deepCut));

            ons = spontAwakeIsis(spontAwakeIsis>=lim(1) & spontAwakeIsis<=lim(2));
            offEs = offEisi(offEisi>=lim(1) & offEisi<=lim(2));
            offDs = offDisi(offDisi>=lim(1) & offDisi<=lim(2));

            alloni = ons;
            alloffie = offEs;
            alloffid = offDs;
            allonih = histcounts(log10(alloni),edges,'Normalization','probability');
            alloffihe = histcounts(log10(alloffie),edges,'Normalization','probability');
            alloffihd = histcounts(log10(alloffid),edges,'Normalization','probability');

            ec = corrcoef(allonih, alloffihe);
            dc = corrcoef(allonih, alloffihd);
            clusters(n).spontStateCorrE = round(ec(2,1),2);
            clusters(n).spontStateCorrD = round(dc(2,1),2);

            clusters(n).spontRatioVigilanceE = round((clusters(n).spontAwakeFR-clusters(n).spontSleepEarlyFR)/(clusters(n).spontAwakeFR+clusters(n).spontSleepEarlyFR),2);
            clusters(n).spontRatioVigilanceD = round((clusters(n).spontAwakeFR-clusters(n).spontSleepDeepFR)/(clusters(n).spontAwakeFR+clusters(n).spontSleepDeepFR),2);
            clusters(n).spontAwakeIsis = spontAwakeIsis;
            clusters(n).spontSleepEarlyIsis = offEisi;
            clusters(n).spontSleepDeepIsis = offDisi;

            clusters(n).songChange = round((clusters(n).songFR-clusters(n).spontAwakeFR)/(clusters(n).songFR+clusters(n).spontAwakeFR),2);
            clusters(n).songPostChange = round((clusters(n).songFRPost-clusters(n).spontAwakeFR)/(clusters(n).songFRPost+clusters(n).spontAwakeFR),2);
            clusters(n).pbAwakeChange = round((clusters(n).pbAwakeFR-clusters(n).prePbAwakeFR)/(clusters(n).pbAwakeFR+clusters(n).prePbAwakeFR),2);
            clusters(n).pbAwakeChangeOld = round((clusters(n).pbAwakeFR-clusters(n).spontAwakeFR)/(clusters(n).pbAwakeFR+clusters(n).spontAwakeFR),2);
            clusters(n).pbAwakePostChange = round((clusters(n).pbAwakeFRPost-clusters(n).spontAwakeFR)/(clusters(n).pbAwakeFRPost+clusters(n).spontAwakeFR),2);
            clusters(n).pbSleepChange = round((clusters(n).pbSleepFR-clusters(n).prePbSleepFR)/(clusters(n).pbSleepFR+clusters(n).prePbSleepFR),2);
            clusters(n).pbSleepChangeOld = round((clusters(n).pbSleepFR-spontSleepDeepFR)/(clusters(n).pbSleepFR+spontSleepDeepFR),2);
            clusters(n).pbSleepPostChange = round((clusters(n).pbSleepFRPost-spontSleepDeepFR)/(clusters(n).pbSleepFRPost+spontSleepDeepFR),2);


            %% replay code
% 
%             if w == 1
%                 nonReplaySnips = sort([15619; 15740; 15957; 16495; 16685]);
%                 replaySnips =  sort([15614.22; 15745.7; 15953.2; 16421.9; 16632.1; 16664.6]);
%             elseif w == 2
%                 nonReplaySnips = sort([28217; 28291; 28645; 28799; 29645]);
%                 replaySnips =  sort([28225.2; 28341; 28639.4; 28953; 29582]);
%             elseif w == 3
%                 nonReplaySnips = sort([5264; 5462; 5540; 6003; 6206; 6421; 7002.5; 7385; 8077; 8882; 8998; 9175; 9264; 9309; 9405; 9517; 9599; 9735; 10034.5; 10694]);
%                 replaySnips =  sort([5970.8; 6368; 6559; 6700.1; 6964.1; 7259.16; 7388.4; 7650.3; 7924.2; 8232.5; 8512.1; 9067.1; 9179.55; 9284.5; 9401.8; 9505.63;9801.7; 10195.3; 10402.2; 10767.6]);
%             elseif w == 4
%                 nonReplaySnips = sort([7592.5; 7739.7; 8001; 8123; 8214.5; 8333; 8440; 8616; 9137.5; 9213; 9309; 9455; 9623; 11282; 11161; 11065; 10988; 10338; 10450.6; 10567.8]);
%                 replaySnips = sort([7730.7; 7971.3; 8036.7; 8361.2; 8528.1; 8788.1; 9220.5; 9302.6; 9509.6; 11288.5; 11200; 11099; 10952.13; 10846.9; 10688.7; 10107.3; 10284; 10382.65; 10471.5; 10541.9; ]);
%             elseif w == 5
%                 nonReplaySnips = sort([11000;12000;13000;14000;15000;16000]);
%                 replaySnips = nonReplaySnips + 50;
%             end
% 
%             snipLength = 5;
%             numSnips = length(nonReplaySnips);
%             totalTime = snipLength*numSnips;
% 
%             nonReplay.start = nonReplaySnips;
%             nonReplay.stop = nonReplaySnips+snipLength;
%             nonReplay.center = nonReplay.start + snipLength/2;
%             nonReplay.warp = ones(length(nonReplaySnips),1);
%             [nonReplaySpks, ~] = alignSpikeTimesSi(clusters, [0 snipLength], nonReplay, nonReplay.start);
% 
%             replay.start = replaySnips;
%             replay.stop = replaySnips+snipLength;
%             replay.warp = ones(length(replaySnips),1);
%             replay.center = replay.start + snipLength/2;
%             [replaySpks, ~] = alignSpikeTimesSi(clusters, [0 snipLength], replay, replay.start);
% 
% 
%             trialSpks = [];
%             maxIFR = [];
%             trialIsis = [];
% 
%             for k = 1:length(replay.center)
%                 trial = replaySpks(n).eventTrials == k;
%                 trialSpks = replaySpks(n).eventTimesAligned(trial);
%                 replayRates(k) = length(trialSpks)/snipLength;
%                 trialIsis = diff(trialSpks);
%                 if ~isempty(trialIsis)
%                     maxIFR(k) = max(1./trialIsis);
%                     %                 elseif length(trialSpks)==1
%                     %                     maxIFR(k) = 1;
%                 elseif isempty(trialSpks)
%                     maxIFR(k) = nan;
%                     replayRates(k) = 0;
%                 end
%                 replayBurstNum(k) = length(trialIsis(trialIsis < .01));
%                 replayBursty(k) = replayBurstNum(k)/length(trialIsis);
%             end
%             replayBursty(isnan(replayBursty)) = 0;
%             clusters(n).replayFR = round(mean(replayRates),2);
%             clusters(n).replayBurst = round(nanmean(replayBursty),3);
%             clusters(n).replayMaxIFR = round(nanmean(maxIFR),2);
%             if isnan(clusters(n).replayMaxIFR)
%                 clusters(n).replayMaxIFR = 0;
%             end
% 
%             trialSpks = [];
%             maxIFR = [];
%             trialIsis = [];
% 
%             for k = 1:length(nonReplay.center)
%                 trial = nonReplaySpks(n).eventTrials == k;
%                 trialSpks = nonReplaySpks(n).eventTimesAligned(trial);
%                 nonReplayRates(k) = length(trialSpks)/snipLength;
%                 trialIsis = diff(trialSpks);
%                 if ~isempty(trialIsis)
%                     maxIFR(k) = max(1./trialIsis);
%                     %                 elseif length(trialSpks)==1
%                     %                     maxIFR(k) = 1;
%                 elseif isempty(trialSpks)
%                     maxIFR(k) = nan;
%                     nonReplayRates(k) = 0;
%                 end
%                 nonReplayBurstNum(k) = length(trialIsis(trialIsis < .01));
%                 nonReplayBursty(k) = nonReplayBurstNum(k)/length(trialIsis);
%             end
%             nonReplayBursty(isnan(nonReplayBursty)) = 0;
%             clusters(n).nonReplayFR = round(mean(nonReplayRates),2);
%             clusters(n).nonReplayBurst = round(nanmean(nonReplayBursty),3);
%             clusters(n).nonReplayMaxIFR = round(nanmean(maxIFR),2);
%             if isnan(clusters(n).nonReplayMaxIFR)
%                 clusters(n).nonReplayMaxIFR = 0;
%             end
%             clusters(n).replayChange = round((clusters(n).replayFR-clusters(n).nonReplayFR)/(clusters(n).nonReplayFR+clusters(n).replayFR),2);

            smSpksAll(n) = smSpks;
            smTimesAll(n) = smTimes;
            pmSpksAll(n) = pmSpks;
            pmTimesAll(n) = pmTimes;
            pmsSpksAll(n) = pmsSpks;
            pmsTimesAll(n) = pmsTimes;
            spontTimesAll(n).spont = spontTimes;

            %disp('fix spont saving line 1063')


            %% spont awake around song

            if w == 1
                spontSongSnips = [250; 280; 5875; 5970; 6040]; %5940;
            elseif w == 2
                spontSongSnips = [1060; 1235; 1530; 1600; 6870]; %5960;
            elseif w == 3
                spontSongSnips = [1815; 2285; 2660; 3880; 4560]; %1035;
            elseif w == 4
                spontSongSnips = [1485; 1700; 4330; 4933; 5216]; %3230;
            elseif w == 5
                spontSongSnips = [118; 1614; 2945; 3540; 6480]; %6415;
            end

            spontSnipLength = 10;
            numSpontSnips = length(spontSongSnips);

            spontSong.start = spontSongSnips;
            spontSong.stop = spontSongSnips+spontSnipLength;
            spontSong.center = spontSong.start + spontSnipLength/2;
            spontSong.warp = ones(length(spontSongSnips),1);
            [spontSongSpks, ~] = alignSpikeTimesSi(clusters, [0 spontSnipLength], spontSong, spontSong.start);

            trialSpks = [];
            trialIsi = [];
            bursty = nan(numSpontSnips,1);
            rates =  nan(numSpontSnips,1);
            cv =  nan(numSpontSnips,1);
            counting = numSpontSnips;

            for k = 1:numSpontSnips
                if spontSong.start(k)>= goodDataTimes(1) && spontSong.stop(k)<= goodDataTimes(2)
                    within = 1;
                else
                    within = 0;
                    counting = counting-1;
                end
                trial = spontSongSpks(n).eventTrials == k;
                trialSpks = spontSongSpks(n).eventTimesAligned(trial);
                trialIsi = diff(trialSpks);
                if isempty(trialIsi) && within
                    bursty(k) = 0;
                    cv(k) = nan;
                    rates(k) = 0;
                    if length(trialSpks) == 1
                        rates(k) = length(trialSpks)/spontSnipLength;
                    end
                elseif isempty(trialIsi) && ~within
                    bursty(k) = nan;
                    cv(k) = nan;
                    rates(k) = nan;
                else
                    bursty(k) = nansum(trialIsi < .01)./length(trialIsi);
                    cv(k) = nanstd(trialIsi)./nanmean(trialIsi);
                    rates(k) = length(trialSpks)/spontSnipLength;
                end
            end
            cvMean = nanmean(cv);
            burstyMean = nanmean(bursty);
            burstyMean(isnan(burstyMean)) = 0;
            ratesMean = nanmean(rates);
            clusters(n).spontSongRate = ratesMean;
            clusters(n).spontSongBurst = burstyMean;
            clusters(n).spontSongCV = cvMean;

%             figure('Visible','off');
%             subplot(3,1,[1:2])
%             plotEventRasters(spontSongSpks(n),spontSong,[0 spontSnipLength])
%             ylim([-8 1])
%             set(gca,'ycolor','none')
%             title([string(bird) + '  _  ' + string(clusters(n).clusterID)])
%             subplot(3,1,3)
%             text(0,.5,['GoodTrials: ' + string(counting) + '    Rate: ' + string(round(ratesMean,2)) + '    Burst: ' + string(round(burstyMean,4)) + '    CV: ' + string(round(cvMean,2))])
%             set(gca,'xtick',[],'ytick',[],'box','off','xcolor','none','color','none','ycolor','none')
% 
%             cd([dir 'spontSong/'])
%             print([string(bird) + '_' + string(clusters(n).clusterID) + '_spontSongRasters'], '-djpeg');


            %%
            % figure;
            % subplot(3,1,1)
            % histogram(spontAwakeLat,-1:.005:1,'Normalization','Probability','EdgeAlpha',.5,'EdgeColor','b');
            % hold on;
            % histogram(smLat,-1:.005:1,'Normalization','Probability','EdgeAlpha',.5,'EdgeColor','r')
            % legend('spont','song')
            % title('song vs spont awake')
            %
            % subplot(3,1,2)
            % histogram(spontAwakeLat,-1:.005:1,'Normalization','Probability','EdgeAlpha',.5,'EdgeColor','b');
            % hold on;
            % histogram(pmLat,-1:.005:1,'Normalization','Probability','EdgeAlpha',.5,'EdgeColor','r')
            % title('pb awake vs spont awake')
            % legend('spont','pb')
            %
            % subplot(3,1,3)
            % histogram(spontSleepDeepLat,-1:.005:1,'Normalization','Probability','EdgeAlpha',.5,'EdgeColor','b');
            % hold on;
            % histogram(pmsLat,-1:.005:1,'Normalization','Probability','EdgeAlpha',.5,'EdgeColor','r')
            % title('pb sleep vs spont deep sleep')
            % legend('spont','pb')
            %
            % cd([path 'precision/'])
            % print([string(clusters(n).clusterID) + 'fixed'], '-dpdf','-fillpage');


         end

    counter = counter+1;

    end
    clusters = orderfields(clusters);

    
    %% save stuff

    cd(dataPath)
    save clusters.mat clusters;
    save spontTimesAll.mat spontTimesAll;
    save smSpksAll.mat smSpksAll;
    save smTimesAll.mat smTimesAll;
    save pmSpksAll.mat pmSpksAll;
    save pmTimesAll.mat pmTimesAll;
    save pmsSpksAll.mat pmsSpksAll;
    save pmsTimesAll.mat pmsTimesAll;

%disp('not saving shit')


end
