%% load everything
close all; clear all; clc; warning off;

dir = 'C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/';
opts = detectImportOptions([dir 'goodCellsnew.xlsx']);
opts = setvartype(opts, 'stdmax', 'string');
goodCells = readtable([dir 'goodCellsnew.xlsx'],opts);
birds = string(goodCells{:,1});
cells = string(goodCells{:,2});
toPrint = [(goodCells{:,5})];
startRec = [(goodCells{:,6})];
stopRec = [(goodCells{:,7})];
std = string(goodCells{:,4});
identity = string([goodCells{:,9}]);

deflect = contains(std, '+');
fs = 50000;
wind = .002 * fs;
fontsize = 15;
diffRise = nan(height(goodCells),1);
diffRise70 = nan(height(goodCells),1);
diffWidth = nan(height(goodCells),1);
set(0,'DefaultFigureVisible','on');
counting = 0;


%%
for number = 6%1:height(goodCells)
    if toPrint(number) == 1 && identity(number) ~= 'x'
        clearvars -except diffRise diffRise70 diffWidth fontsize fs wind number identity goodCells birds cells toPrint dir startRec stopRec deflect counting SpontsOnTimes SpontsOffTimes AwakePB AwakeSpont SleepPB SleepSpont MotifsOnTimes MotifsOffTimes AwakeHist SleepHist
        close all;
        path = dir + birds(number) + '/';
        f = char(cells(number));
        out=regexp(path,'/','split');
        name = string(out(end-1));
        disp(name)
        disp(f)
        %new
        % load([path + f + '/clusters.mat'])
        % spkTimes = round(clusters.spikeTimes*fs);

            load(path + 'hp200/' + f + '/times_' + f + '_hp200.mat');        if isfile(path + f + '/fullRecord.mat')
            load(path + f + '/fullRecord.mat');
        else
            fullRecord = abfload(path + f + '/' + f + '_filteredspikes.abf');
            fullRecord(:,2:3) = [];
            cd(path + f);
            save fullRecord.mat fullRecord;
        end

        %goodDataTimes = 1:length(fullRecord(:,1)); %som is 1:795*fs
        startRecFixed = startRec(number)*fs;
        startRecFixed(startRecFixed == 0) = 1;
        if strcmp(stopRec{number},'end')
            stopRecFixed = length(fullRecord(:,1));
        else
            stopRecFixed = str2num(stopRec{number})*fs;
        end
        goodDataTimes = startRecFixed:stopRecFixed;
        load(path + 'hp200/' + f + '/' + f + '_hp200.mat');
        spikes = data(:,1);
        clear data;
        eyeTrig = fullRecord(:, 2);

        clear fullRecord;
        lengthRecording = length(goodDataTimes)/fs;
        motif = audioread(path + 'motif.wav'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%change this
        load(path + f + '/MotifTimes.mat');

        structElem = fieldnames(Motif);
        motifRange = [floor(startRecFixed/fs) stopRecFixed/fs];
        elem = Motif.start;
        takeMotif = elem>=motifRange(1) & elem<=motifRange(2);
        for i = 2:5
            Motif.(structElem{i}) = Motif.(structElem{i})(takeMotif);
        end

        clusters.bird = name;
        clusters.cell = string(f);
        clusters.recLength = lengthRecording;


        %% get spike times (from ms)

        spkTimesAll = round(cluster_class(:,2) / 1000 * fs);
        spkTimes = spkTimesAll(cluster_class(:,1)==1);
        spkTimes = spkTimes(spkTimes >= startRecFixed & spkTimes <= stopRecFixed);
        isis = diff(spkTimes)/fs;
        fixed = isis > .0015;
        spkTimes = spkTimes(fixed == 1);
        clusters.spikeTimes = spkTimes/fs;


        %% plot spike shape

        timeWindow = wind+20;
        xAxisWaveForm = linspace(-timeWindow/fs/2*1000, timeWindow/fs/2*1000, timeWindow+1);
        allWF = zeros(timeWindow+1, length(spkTimes));
        start = 55;
        stop = 65;

        for i = 1:length(spkTimes)
            span = spkTimes(i)-timeWindow/2 : spkTimes(i)+timeWindow/2;
            allWF(:,i) = spikes(span);
        end

        rng(number)
        randChoose = randperm(length(spkTimes),min(length(spkTimes),200));
        randWindow = allWF(:,randChoose);
        meanWF = mean(allWF,2);

        [h,g] = butter(3, 500/(fs/2),'high');

        if ~deflect(number)
            center = find(islocalmin(meanWF,'MaxNumExtrema',1));
        else
            center = find(islocalmax(meanWF,'MaxNumExtrema',1));
        end
        allWF2 = filtfilt(h, g, allWF);
        numWavs = size(allWF2,2);

        fixedWF = zeros(101,numWavs);
        for i = 1:numWavs
            if ~deflect(number)
                eachCenter = start+find(islocalmin(allWF2(start:stop,i),'MaxNumExtrema',1));
            else
                eachCenter = start+find(islocalmax(allWF2(start:stop,i),'MaxNumExtrema',1));
            end

            if ~isempty(eachCenter)
                jitter = eachCenter-center;
            else
                jitter = 0;
            end
            fixedWF(:,i) = allWF2(10+jitter:110+jitter,i);
        end

        meanWF2 = mean(fixedWF, 2);
        clusters.meanWF = meanWF2;

        %clusters.riseTime = abs(min(meanWF)-max(meanWF))
        f1 = figure('WindowState','maximized');
        clf;
        hold on;
        %subplot(1,2,1)
        for i = 1:length(randChoose)
            plot(xAxisWaveForm, randWindow(:,i),'Color',[.5 .5 .5]);
            hold on;
        end
        hold on;
        plot(xAxisWaveForm, meanWF, 'LineWidth', 4, 'Color', 'k');
        set(gca,'box','off','TickDir','out','FontSize',fontsize,'Color','none');
        xlabel('ms');
        ylabel('mV');
        title([name + ' ' + f]);
        cd([path + f])
        print([name + ' ' + f + ' wf'], '-dpdf');


        %% get spike params

        meanSubWav = meanWF-mean(meanWF);
        maxAmp = max(abs(meanSubWav));
        normWav = meanSubWav./maxAmp;
        %     normWavShort = normWav(10:end);
        clusters.wfNorm = normWav;
        factor = 1000;
        normWavInterp = round(interp(normWav,factor),5);
        %    clusters.interp = normWavInterp;

        %if ~deflect(number)
        if deflect(number)
            normWavInterp = normWavInterp * (-1);
        end
        trVal = min(normWavInterp);
        trPos = find(islocalmin(normWavInterp,'MaxNumExtrema',1));
        pkVal = max(normWavInterp(trPos:end));
        pkPos = (trPos-1)+find(islocalmax(normWavInterp(trPos:end),'MaxNumExtrema',1));
        trPkTime = (pkPos-trPos)/fs*1000/factor;

        pkVal70 = trVal+.7*range([pkVal trVal]);
        closeEnuf = find(abs(normWavInterp - pkVal70) < .0015);
        inBet = closeEnuf(closeEnuf < pkPos & closeEnuf >  trPos);
        pkPos70 = inBet(islocalmin(abs(normWavInterp(inBet) - pkVal70),'MaxNumExtrema',1));
        trPk70 = (pkPos70-trPos)/fs*1000/factor;

        pkVal50 = trVal+.5*range([pkVal trVal]);
        closeEnuf = find(abs(normWavInterp - pkVal50) < .0015);
        beforePk =  closeEnuf(closeEnuf < pkPos);
        pkPos50 = beforePk(islocalmin(abs(normWavInterp(beforePk) - pkVal50),'MaxNumExtrema',2));
        width = diff(pkPos50)/fs*1000/factor;

        clusters.wfRise = trPkTime;
        clusters.wfRise70 = trPk70;
        clusters.wfFWHM = width;
        clusters.wfMult70 = (clusters.wfRise-clusters.wfRise70).*(clusters.wfRise./clusters.wfRise70);

        %         subplot(1,2,2);
        %         plot([normWavShort(2:end)],'k')
        %         hold on
        %         set(gca,'box','off','TickDir','out','FontSize',fontsize,'Color','none');
        %
        %         hold on
        %         plot(pkPos70/factor,pkVal70,'b*')
        %         hold on
        %         line([trPos pkPos70]/factor,[trVal pkVal70],'Color','b')
        %         hold on
        %         plot(pkPos/factor,pkVal,'r*')
        %         hold on
        %         plot(trPos/factor,trVal,'r*')
        %         hold on
        %         line([pkPos50 pkPos50]/factor,[pkVal50 pkVal50],'Color','b')
        %
        %         text(60,-.4,['rise: ' + string(round(trPkTime,3))],'FontSize',fontsize)
        %         text(60,-.5,['rise70: ' + string(round(trPk70,3))],'FontSize',fontsize)
        %         text(60,-.6,['width: ' + string(round(width,3))],'FontSize',fontsize)

        %         else
        %             clusters.wfFWHM = nan;
        %             clusters.wfRise = nan;
        %             clusters.wfRise70 = nan;
        %             clusters.wfMult70 = nan;
        %         end


        %% derivative
        %
        %            %     if deflect(number)
        %                     clear trVal trPos pkVal pkPostrPkTime pkVal70 closeEnuf inBet pkPos70 trPk70 pkVal50 beforePk pkPos50 width
        %                     normWavShort2 = diff(normWavShort);
        %                     normWavInterp = round(interp(normWavShort2,factor),5);
        %                     trVal = min(downsample(normWavInterp,factor));
        %                     trPos = find(islocalmin(normWavInterp,'MaxNumExtrema',1));
        %                     pkVal = max(downsample(normWavInterp(trPos:end),factor));
        %                     pkPos = (trPos-1)+find(islocalmax(normWavInterp(trPos:end),'MaxNumExtrema',1));
        %                     trPkTime = (pkPos-trPos)/fs*1000/factor;
        %                     clusters.derivRise = trPkTime;
        %
        %                     pkVal70 = trVal+.7*range([pkVal trVal]);
        %                     closeEnuf = find(abs(normWavInterp - pkVal70) < .0015);
        %                     inBet = closeEnuf(closeEnuf < pkPos & closeEnuf >  trPos);
        %                     pkPos70 = inBet(islocalmin(abs(normWavInterp(inBet) - pkVal70),'MaxNumExtrema',1));
        %                     trPk70 = (pkPos70-trPos)/fs*1000/factor;
        %                     clusters.derivRise70 = trPk70;
        %
        %                     pkVal50 = trVal+.5*range([pkVal trVal]);
        %                     closeEnuf = find(abs(normWavInterp - pkVal50) < .0015);
        %                     beforePk =  closeEnuf(closeEnuf < pkPos);
        %                     pkPos50 = beforePk(islocalmin(abs(normWavInterp(beforePk) - pkVal50),'MaxNumExtrema',2));
        %                     width = diff(pkPos50)/fs*1000/factor;
        %                     clusters.derivWidth = width;
        %
        %                     subplot(1,3,3);
        %                     plot(normWavShort2(2:end),'k')
        %                     hold on
        %
        %                     set(gca,'box','off','TickDir','out','FontSize',fontsize,'Color','none');
        %                     hold on
        %                     plot(pkPos70/factor,pkVal70,'b*')
        %                     hold on
        %                     line([trPos pkPos70]/factor,[trVal pkVal70],'Color','b')
        %                     hold on
        %                     plot(pkPos/factor,pkVal,'r*')
        %                     hold on
        %                     plot(trPos/factor,trVal,'r*')
        %                     hold on
        %                     line([pkPos50 pkPos50]/factor,[pkVal50 pkVal50],'Color','b')
        %                     %   title(clusters.clusterID)
        %                     hold on
        %                     text(60,-.08,['rise: ' + string(round(trPkTime,3))],'FontSize',fontsize)
        %                     text(60,-.1,['rise70: ' + string(round(trPk70,3))],'FontSize',fontsize)
        %                     text(60,-.12,['width: ' + string(round(width,3))],'FontSize',fontsize)
        %          %       end
        %                 set(gcf,'PaperOrientation','landscape');
        %                 print(['C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/wfs/' + name + f + ' wf_deriv'], '-bestfit', '-dpdf');


        %% get isis

        isi_sort = sort(diff(clusters.spikeTimes));
        clusters.isis = isi_sort;
        isi = diff(clusters.spikeTimes);
        ifr = 1./isi;
        isiTimes = (spkTimes(1:end-1) + .5*isi)/fs;


        %% split pb and spont times

        buffer = .5;
        pbTrial = [Motif.start; Motif.stop+buffer]';
        pbTimes = [];
        for i = 1:length(pbTrial)
            clear pbSpks;
            pbSpks = clusters.spikeTimes(clusters.spikeTimes>pbTrial(i,1) & clusters.spikeTimes < pbTrial(i,2));
            pbSpks = [pbSpks; nan];
            pbTimes = [pbTimes; pbSpks];
        end
        spontTimes = clusters.spikeTimes;
        [~,b] = ismember(pbTimes(~isnan(pbTimes)),spontTimes);
        spontTimes(b)=nan;

        pbIsis = diff(pbTimes);
        pb_isi_sort = sort(pbIsis);
        clusters.pbIsis = pb_isi_sort;

        spontIsis = diff(spontTimes);
        spont_isi_sort = sort(spontIsis);
        clusters.spontIsis = spont_isi_sort;


        %% open or make eye stack

        eyeRate = 4;

        if ~contains(name,'xbeads')
            %lum = eyes.lum;
            %onThresh = eyes.onThresh;
            if isfile(path + f + '/stack.mat')
                load(path + f + '/stack.mat');
            else
                fname = char(path + f + '/eyes_' + f + '.tif');
                info = imfinfo(fname);
                numImg = length(info);
                stack = [];
                for k = 1:numImg
                    currentImage = imread(fname, k, 'Info', info);
                    stack(:,:,k) = currentImage;
                end
                cd(path + f)
                save stack.mat stack
            end
            lum = [mean(mean(stack))];
            lum = lum(:);


            %% check eye vid against eye trigs

            binEyes = diff(eyeTrig);
            [frameLoc] = find(binEyes > .5);
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
            switchOff = [(goodCells{number,10})];
            lum = [ones(switchOff*eyeRate,1); zeros(round((lengthRecording-switchOff)*eyeRate),1)];
        end


        %% get segs from lum

        onThresh = mean(lum); %input('what threshold? \n \n');
        eyePts = zeros(length(lum),1);
        trigs_on = (lum>=onThresh);
        eyePts(trigs_on) = 1;

        lightSeg = find(diff(eyePts ~= 0));
        lightSeg = lightSeg(lightSeg>=floor(startRecFixed/fs*eyeRate) & lightSeg<=floor(stopRecFixed/fs*eyeRate));
        lightSeg = [floor(startRecFixed/fs*eyeRate); lightSeg; floor(stopRecFixed/fs*eyeRate)];

        segs = [];
        for i = 1:length(lightSeg)-1
            segs(i,:) = [lightSeg(i) lightSeg(i+1)];
        end

        onTimes = [];
        offTimes = [];
        for i = 1:size(segs,1)
            [y, z] = find(spkTimes/fs>=segs(i,1)/eyeRate & spkTimes/fs<=segs(i,2)/eyeRate);
            if eyePts(segs(i,2))
                onTimes = [onTimes; spkTimes(y)/fs; nan];
            else
                offTimes = [offTimes; spkTimes(y)/fs; nan];
            end
        end

        pbAwake = ismember(onTimes,pbTimes);
        pbAwakeTimes = onTimes(pbAwake);
        pbAwakeIsis = sort(diff(pbAwakeTimes));
        clusters.onPbIsis = pbAwakeIsis;

        onTimes(pbAwake) = nan;
        onIsis = diff(onTimes);
        on_isi_sort = sort(onIsis);
        clusters.onIsis = on_isi_sort;

        pbSleep = ismember(offTimes,pbTimes);
        pbSleepTimes = offTimes(pbSleep);
        pbSleepIsis = sort(diff(pbSleepTimes));
        clusters.offPbIsis = pbSleepIsis;

        offTimes(pbSleep) = nan;
        offIsis = diff(offTimes);
        off_isi_sort = sort(offIsis);
        clusters.offIsis = off_isi_sort;

        edges = -3:.1:1;
        alloni = clusters.onIsis(clusters.onIsis>=1e-3 & clusters.onIsis<=10);
        alloffi = clusters.offIsis(clusters.offIsis>=1e-3 & clusters.offIsis<=10);
        allonih = histcounts(log10(alloni),edges,'Normalization','probability');
        alloffih = histcounts(log10(alloffi),edges,'Normalization','probability');
        c = corrcoef(allonih, alloffih);
        clusters.spontStateCorr = round(c(2,1),2);

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

        %         gm = [22 27 9];
        %         if ismember(number, gm)
        %             [~,countSpOn] = alignSpikeTimesJuxta(clusters, windJux, Motif);
        %             countSpOff = countSpOn;
        %         end


        %% get spontaneous activity snippets

        spontNumTrials = 64;
        segsS = segs/eyeRate;
        fixedSegsOn = [];
        fixedSegsOff = [];
        for i = 1:size(segs,1)
            if eyePts(segs(i,2)) == 1
                seg = segsS(i,:);
                for j = 1:length(pbTrial)
                    while seg(2) > pbTrial(j,2) && seg(1) < pbTrial(j,1)
                        fixedSegsOn = [fixedSegsOn; seg(1) pbTrial(j,1)];
                        seg = [pbTrial(j,2) seg(2)];
                    end
                end
                fixedSegsOn = [fixedSegsOn; seg];
            elseif eyePts(segs(i,2)) == 0
                seg = segsS(i,:);
                for j = 1:length(pbTrial)
                    while seg(2) > pbTrial(j,2) && seg(1) < pbTrial(j,1)
                        fixedSegsOff = [fixedSegsOff; seg(1) pbTrial(j,1)];
                        seg = [pbTrial(j,2) seg(2)];
                    end
                end
                fixedSegsOff = [fixedSegsOff; seg];
            end
        end

        if ~isempty(fixedSegsOn)
            segLengthsOn = fixedSegsOn(:,2)-fixedSegsOn(:,1);
            longSegsOn = segLengthsOn >= 30;
            longEnufOn = segLengthsOn >= 5 & segLengthsOn < 30;
            segsOn = fixedSegsOn(longEnufOn,:);
            for i = 1:length(segLengthsOn)
                if longSegsOn(i)
                    shortSegsOn = [];
                    left = floor(segLengthsOn(i)/20);
                    shortSegsOn(:,1) = fixedSegsOn(i,1):20:fixedSegsOn(i,2);
                    shortSegsOn(:,2) = shortSegsOn(:,1)+20;
                    shortSegsOn(left+1,:) = [shortSegsOn(left,2) fixedSegsOn(i,2)];
                    segsOn = [segsOn; shortSegsOn];
                end
            end
            spontSegsOn = sort(segsOn(segsOn(:,2)-segsOn(:,1)>=5,:));

            %spontSegsOn(spontSegsOn(:,1)>=motifRange(1) & spontSegsOn(:,2) <= motifRange(2));
            rng(number)
            whichSpontOn = randperm(min(size(spontSegsOn,1)));
            count = 1;
            for i = whichSpontOn(1:(min(spontNumTrials,length(whichSpontOn))))
                a = ceil(spontSegsOn(i,1));
                b = ceil(spontSegsOn(i,2))-5;
                rng(number)
                snip = a + (b-a).*rand(1,1);
                SpontOn.start(count) = snip;
                SpontOn.stop(count) = snip+5;
                SpontOn.center(count) = snip+2.5;
                count = count + 1;
            end
            SpontOn.warp = ones(length(SpontOn.start),1);
        else
            SpontOn.start = [];
            SpontOn.stop = [];
            SpontOn.center = [];
            SpontOn.warp = [];
        end

        if ~isempty(fixedSegsOff)
            segLengthsOff = fixedSegsOff(:,2)-fixedSegsOff(:,1);
            longSegsOff = segLengthsOff >= 30;
            longEnufOff = segLengthsOff >= 5 & segLengthsOff < 30;
            segsOff = fixedSegsOff(longEnufOff,:);
            for i = 1:length(segLengthsOff)
                if longSegsOff(i)
                    shortSegsOff = [];
                    left = floor(segLengthsOff(i)/20);
                    shortSegsOff(:,1) = fixedSegsOff(i,1):20:fixedSegsOff(i,2);
                    shortSegsOff(:,2) = shortSegsOff(:,1)+20;
                    shortSegsOff(left+1,:) = [shortSegsOff(left,2) fixedSegsOff(i,2)];
                    segsOff = [segsOff; shortSegsOff];
                end
            end
            spontSegsOff = sort(segsOff(segsOff(:,2)-segsOff(:,1)>=5,:));

            %spontSegsOff(spontSegsOff(:,1)>=motifRange(1) & spontSegsOff(:,2) <= motifRange(2));
            rng(number)
            whichSpontOff = randperm(min(size(spontSegsOff,1)));

            count = 1;
            for i = whichSpontOff(1:(min(spontNumTrials,length(whichSpontOff))))
                a = ceil(spontSegsOff(i,1));
                b = ceil(spontSegsOff(i,2))-5;
                rng(number)
                snip = a + (b-a).*rand(1,1);
                SpontOff.start(count) = snip;
                SpontOff.stop(count) = snip+5;
                SpontOff.center(count) = snip+2.5;
                count = count + 1;
            end
            SpontOff.warp = ones(length(SpontOff.start),1);
        else
            SpontOff.start = [];
            SpontOff.stop = [];
            SpontOff.center = [];
            SpontOff.warp = [];
        end


        %% amplitude distributions and ifr

        % amp = spikes(spkTimes);

        numPts = round(length(spikes)/fs*1000,0);
        plotPts = zeros(numPts,1);
        roundTimes = round(spkTimes(2:end)/fs*1000,0);
        plotPts(roundTimes) = ifr;
        eyePtsNAN = eyePts;
        eyePtsNAN(eyePtsNAN == 0) = nan;
        indEyePts = 1:length(eyePtsNAN);
        eyePtsNAN(indEyePts < motifRange(1)*eyeRate) = nan;
        eyePtsNAN(indEyePts > motifRange(2)*eyeRate) = nan;

        % a(1) = subplot(2,1,1);
        % plot(spkTimes/fs,amp,'k.')
        % hold on;
        % whichTick = xticks;
        % xax = whichTick/fs;
        %
        % hold on
        % ylabel('mV')
        % if par.detection == 'neg'
        %     ylim([min(amp)-.1 .1])
        % elseif par.detection == 'pos'
        %     ylim([-.1 max(amp)+.1])
        % end
        % title('amp')
        % set(gca, 'box', 'off', 'color', 'none', 'FontSize', fontsize, 'TickDir','out');
        % xlim([0 lengthRecording+5])
        %
        % a(2) = subplot(2,1,2);

        %         f7 = figure('WindowState','maximized');
        %         clf;
        %         plot(linspace(0,length(spikes)/fs,numPts),plotPts,'k')
        %         title('ifr')
        %         ylabel('Hz')
        %         xlabel('sec')
        %         set(gca, 'box', 'off', 'color', 'none', 'FontSize', fontsize, 'TickDir','out');
        %         xlim([0 length(spikes)/fs+5])
        %         hold on
        %         if ~isempty(fixedSegsOn)
        %             for i = 1:length(SpontOn.start)
        %                 patch([SpontOn.start(i) SpontOn.stop(i) SpontOn.stop(i) SpontOn.start(i)], [max(ylim)*.9 max(ylim)*.9 max(ylim)*.8 max(ylim)*.8], [1 0.55 0.125], 'EdgeColor', 'k', 'EdgeAlpha',.2,'FaceAlpha', .1)
        %                 hold on
        %             end
        %         end
        %         hold on
        %         if ~isempty(fixedSegsOff)
        %             for i = 1:length(SpontOff.start)
        %                 patch([SpontOff.start(i) SpontOff.stop(i) SpontOff.stop(i) SpontOff.start(i)], [max(ylim)*.9 max(ylim)*.9 max(ylim)*.8 max(ylim)*.8], 'k', 'EdgeColor', 'k', 'EdgeAlpha',.2,'FaceAlpha', .1)
        %                 hold on
        %             end
        %         end

        %         hold on
        %         if ~isempty(MotifOn.start)
        %             for i = 1:length(MotifOn.start)
        %                 patch([MotifOn.start(i) MotifOn.stop(i) MotifOn.stop(i) MotifOn.start(i)], [max(ylim)*.9 max(ylim)*.9 max(ylim) max(ylim)], [1 0.55 0.125],'EdgeColor', 'k', 'EdgeAlpha',.2,'FaceAlpha', .1)
        %                 hold on
        %             end
        %         end
        %         hold on
        %         if ~isempty(MotifOff.start)
        %             for i = 1:length(MotifOff.start)
        %                 patch([MotifOff.start(i) MotifOff.stop(i) MotifOff.stop(i) MotifOff.start(i)], [max(ylim)*.9 max(ylim)*.9 max(ylim) max(ylim)], 'k','EdgeColor', 'k', 'EdgeAlpha',.2,'FaceAlpha', .1)
        %                 hold on
        %             end
        %         end
        %         hold on
        %         plot(linspace(0, length(spikes)/fs, length(eyePtsNAN)), eyePtsNAN*max(ylim), '.', 'Color', 'k', 'LineWidth',20)
        %         hold onmatch
        %         line([motifRange(1) motifRange(1)],[0 max(ylim)], 'Color','r')
        %         hold on
        %         line([motifRange(2) motifRange(2)],[0 max(ylim)], 'Color','r')
        %         cd([path + f])
        %         set(gcf,'PaperOrientation','landscape', 'PaperUnits', 'inches', 'PaperSize', [20 30]);
        %
        %         print([name + ' ' + f + ' annotated ifr'], '-dpdf');
        %         cd(['C:/Users\User\Dropbox (NYU Langone Health)\Int Juxta\amps final\'])
        %         print([name + ' ' + f + ' amps and ifr'], '-dpdf');
        clear spikes
        clear eyeTrig
        clear stack


        %% plot spont snippets

        ticklen = 3;
        windSpont = [0 5];
        sigwin10 = .01;
        sigwin50 = .05;
        Trials = [length(SpontOn.start); length(SpontOff.start); length(MotifOn.start); length(MotifOff.start)];

        f6 = figure;
        clf;

        %   subplot(20,1,[1:16])
        %    count = 1;
        onTri = length(SpontOn.start);
        offTri = length(SpontOff.start);
        totTri = onTri+offTri+2;

        if ~isempty(fixedSegsOn)
            clear trialIsi trialSpks
            [spontStructOn,psthOn] = alignSpikeTimesJuxta(clusters, windSpont, SpontOn);
            subplot(totTri,1,1:onTri+1)

            plotEventRasters(spontStructOn,SpontOn,windSpont)
            title('spontaneous activity')
            ylabel([num2str(onTri) +'Trials'])
            set(gca,'box','off','color','none','TickDir','out','xcolor',[1 0.55 0.125],'ycolor',[1 0.55 0.125],'FontSize',fontsize)

            for k = 1:length(SpontOn.start)

                trial = spontStructOn.eventTrials == k;
                trialSpks = spontStructOn.eventTimesAligned(trial);
                trialIsi = diff(trialSpks);
                %                 for i = 1:length(trialSpks)
                %                     line([trialSpks(i) trialSpks(i)], [0 ticklen]-ticklen*count, 'Color',[1 0.55 0.125], 'LineWidth', 1)
                %                     hold on
                %                 end
                %                 hold on
                if isempty(trialIsi)
                    burstOn(k) = 0;
                    cvOn(k) = nan;
                    if isempty(trialSpks)
                        onRates(k) = 0;
                    else
                        onRates(k) = length(trialSpks)/windSpont(2);
                    end
                else
                    burstOn(k) = nansum(trialIsi < .01)./length(trialIsi);
                    cvOn(k) = nanstd(trialIsi)/nanmean(trialIsi);
                    onRates(k) = length(trialSpks)/windSpont(2);
                end
                count = count + 1;
            end
            rangeOn = 1:Trials(3);
            if number == 79
                rangeOn = 1:11;
            end
            for i = 2:5
                subsetOn.(structElem{i}) = SpontOn.(structElem{i})(rangeOn);
            end
        else
            spontStructOn.eventTimes = [];
            spontStructOn.eventTimesAligned = [];
            spontStructOn.eventTrials = [];
        end

        %             [SpontOnTimes] = alignSpikeTimesJuxta(clusters, [0 length(motif)/fs], subsetOn);
        %             spontLatOn = [];
        %             if ~isempty(SpontOnTimes.eventTimes)
        %                 for i = 1:length(SpontOnTimes.eventTimesAligned)
        %                     t = SpontOnTimes.eventTrials(i);
        %                     compare = SpontOnTimes.eventTimesAligned(SpontOnTimes.eventTrials~=t);
        %                     latencies = compare-SpontOnTimes.eventTimesAligned(i);
        %                     spontLatOn = [spontLatOn; latencies'];
        %                     clear latencies
        %                 end
        %                 propSpontOn10 = sum(spontLatOn >= -sigwin10 & spontLatOn <= sigwin10)/length(spontLatOn);
        %                 propSpontOn50 = sum(spontLatOn >= -sigwin50 & spontLatOn <= sigwin50)/length(spontLatOn);
        %             else
        %                 spontLatOn = [];
        %                 propSpontOn10 = nan;
        %                 propSpontOn50 = nan;
        %             end
        %         else
        %             propSpontOn10 = nan;
        %             propSpontOn50 = nan;
        %         end

        if ~isempty(fixedSegsOff)
            clear trialIsi trialSpks

            [spontStructOff,psthOff] = alignSpikeTimesJuxta(clusters, windSpont, SpontOff);
            subplot(totTri,1,onTri+2:totTri)
            plotEventRasters(spontStructOff,SpontOff,windSpont)
            set(gca,'box','off','color','none','TickDir','out','xcolor',[0.3010, 0.7450, 0.9330],'ycolor',[0.3010, 0.7450, 0.9330],'FontSize',fontsize)

            for k = 1:length(SpontOff.start)
                trial = spontStructOff.eventTrials == k;
                trialSpks = spontStructOff.eventTimesAligned(trial);
                trialIsi = diff(trialSpks);
                %                 for i = 1:length(trialSpks)
                %                     line([trialSpks(i) trialSpks(i)], [0 ticklen]-ticklen*count, 'Color','k', 'LineWidth', 1)
                %                     hold on
                %                 end
                %                 hold on
                if isempty(trialIsi)
                    burstOff(k) = 0;
                    cvOff(k) = nan;
                    offRates(k) = 0;
                    if isempty(trialSpks)
                        offRates(k) = 0;
                    else
                        offRates(k) = length(trialSpks)/windSpont(2);
                    end
                else
                    burstOff(k) = nansum(trialIsi < .01)./length(trialIsi);
                    cvOff(k) = nanstd(trialIsi)/nanmean(trialIsi);
                    offRates(k) = length(trialSpks)/windSpont(2);
                end
                %                 count = count + 1;
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
        %             [SpontOffTimes] = alignSpikeTimesJuxta(clusters, [0 length(motif)/fs], subsetOff);
        %             spontLatOff = [];
        %             if ~isempty(SpontOffTimes.eventTimes)
        %                 for i = 1:length(SpontOffTimes.eventTimesAligned)
        %                     t = SpontOffTimes.eventTrials(i);
        %                     compare = SpontOffTimes.eventTimesAligned(SpontOffTimes.eventTrials~=t);
        %                     latencies = compare-SpontOffTimes.eventTimesAligned(i);
        %                     spontLatOff = [spontLatOff; latencies'];
        %                     clear latencies;
        %                 end
        %                 propSpontOff10 = sum(spontLatOff >= -sigwin10 & spontLatOff <= sigwin10)/length(spontLatOff);
        %                 propSpontOff50 = sum(spontLatOff >= -sigwin50 & spontLatOff <= sigwin50)/length(spontLatOff);
        %             else
        %                 spontLatOff = [];
        %                 propSpontOff10 = nan;
        %                 propSpontOff50 = nan;
        %             end
        %         else
        %             propSpontOff10 = nan;
        %             propSpontOff50 = nan;
        %end

        xlim(windSpont)
        xticks(windSpont(1):windSpont(2));
        spontxLabels = string(0:5);
        xticklabels(spontxLabels);
        %  snips = count-1;
        %  ylim([-(snips+1)*ticklen 0])
        xlabel('sec')
        ylabel([num2str(offTri) +'Trials'])


        cd(path + f)
        print(name + ' ' + f + ' spont', '-dpdf');


        %% spont stats

        pbBuffer = length(motif)/fs+buffer;
        % lengthNonPb = lengthRecording - pbBuffer*length(Motif.center);

        %cvSpontOn = nanstd(onIsis)/nanmean(onIsis);
        %cvSpontOff = nanstd(offIsis)/nanmean(offIsis);
        %burstSpontOn = nansum(onIsis < .01)./lengthNonPb;
        %burstSpontOff = nansum(offIsis < .01)./lengthNonPb;
        %numSpksOnSpont = length(onTimes) - sum(isnan(onTimes));
        %rateSpontOn = numSpksOnSpont/lengthNonPb;
        %numSpksOffSpont = length(offTimes) - sum(isnan(offTimes));
        %rateSpontOff = numSpksOffSpont/lengthNonPb;

        f7 = figure('WindowState','normal');
        clf;
        txty = .9;
        txtx1 = .6;
        txtx2 = 1;
        xlim([0 1.6])
        ylim([-.4 1])
        set(gca,'box','off','color','none','xcolor','none','ycolor','none','FontSize',fontsize)
        text(txtx1-.15,txty+.2,'SPONTANEOUS STATS','FontSize',fontsize)

        text(txtx1,txty,'On','FontSize',fontsize)
        text(txtx2,txty,'Off','FontSize',fontsize)

        % num trials
        text(0,txty-.2,'Trials','FontSize',fontsize)

        % spont mean fr
        text(0,txty-.4,'Rate (Hz)','FontSize',fontsize)

        % spont cv
        text(0,txty-.8,'CV','FontSize',fontsize)

        % spont bursty
        text(0,txty-.6,'Burst %','FontSize',fontsize)
        if ~isempty(fixedSegsOn)
            cvSpontOn = nanmean(cvOn);
            burstSpontOn = nanmean(burstOn);
            rateSpontOn = mean(onRates); %nanmean(psthOn);
            clusters.spontAwakeFR = round(rateSpontOn,2);
            clusters.spontAwakeBurst = round(burstSpontOn,2);
            clusters.spontAwakeCV = round(cvSpontOn,2);
            text(txtx1,txty-.2,num2str(length(SpontOn.start)),'FontSize',fontsize)
            text(txtx1,txty-.4,num2str(clusters.spontAwakeFR),'FontSize',fontsize)
            text(txtx1,txty-.8,num2str(clusters.spontAwakeCV),'FontSize',fontsize)
            text(txtx1,txty-.6,num2str(clusters.spontAwakeBurst),'FontSize',fontsize)
        else
            clusters.spontAwakeFR = nan;
            clusters.spontAwakeBurst = nan;
            clusters.spontAwakeCV = nan;
        end
        if ~isempty(fixedSegsOff)
            cvSpontOff = nanmean(cvOff);
            burstSpontOff = nanmean(burstOff);
            rateSpontOff = mean(offRates);%nanmean(psthOff);
            clusters.spontSleepFR = round(rateSpontOff,2);
            clusters.spontSleepBurst = round(burstSpontOff,2);
            clusters.spontSleepCV = round(cvSpontOff,2);
            text(txtx2,txty-.6,num2str(clusters.spontSleepBurst), 'FontSize',fontsize)
            text(txtx2,txty-.4,num2str(clusters.spontSleepFR),'FontSize',fontsize)
            text(txtx2,txty-.2,num2str(length(SpontOff.start)),'FontSize',fontsize)
            text(txtx2,txty-.8,num2str(clusters.spontSleepCV),'FontSize',fontsize)
        else
            clusters.spontSleepFR = nan;
            clusters.spontSleepBurst = nan;
            clusters.spontSleepCV = nan;
        end

        % spont precise
        %text(0,txty-1,'Precise','FontSize',fontsize)
        %text(txtx1,txty-1,string(round(propSpontOn,3))', 'FontSize',fontsize)
        %text(txtx2,txty-1,string(round(propSpontOff,3)), 'FontSize',fontsize)

        cd(path + f)
        print(name + ' ' + f + ' spont stats', '-dpdf');


        %% plot on vs off ISIs

        edges = -3:.1:1;
        ticks = -3:1;
        labels = 10.^ticks;
        nameLabels = labels*1000;

        f5 = figure('WindowState','normal');
        clf;
        % hold on;
        % subplot(2,1,2)
        %[~,edges,bins] = histcounts(log10(on_isi_sort),20);
        histogram(on_isi_sort,10.^edges,'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 0.55 0.125]);
        set(gca, 'xscale', 'log','TickDir','out');

        hold on
        %[~,edges] = histcounts(log10(off_isi_sort),20);
        histogram(off_isi_sort,10.^edges,'Normalization','probability','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[0.3010, 0.7450, 0.9330]);
        set(gca, 'xscale', 'log', 'box', 'off', 'color', 'none', 'FontSize', fontsize,'TickDir','out');
        title('spontaneous ISIs');
        %xlabel('ISI (sec)');
        ylabel('ISI counts (prob.)');
        legend('lights on','lights off','Location','northeast');
        %xlim([1e-3 10]);
        %xticks([1e-3 1e-2 1e-1 1 1e1 1e2])

        xticks(labels);
        set(gca,'XTickLabels',nameLabels)
        xlabel('ms')

        cd(path + f)
        print(name + ' ' + f + ' isis', '-dpdf');


        %% plot all pb trials raster

        f2 = figure('WindowState','normal');

        clf;
        axis(1) = subplot(17,1,1:3);
        plot(linspace(windJux(1),windJux(2), length(paddedAudio)), paddedAudio, 'k')
        set(gca,'xtick',[],'ytick',[],'box','off','xcolor','none','ycolor','none','color','none','FontSize',fontsize)
        xlim(windJux)
        title('playback activity')
        hold on
        axis(2) = subplot(17,1,4:11);
        count = 1;
        for k = 1:length(MotifOn.center)
            trial = outStructOn.eventTrials == k;
            trialSpks = outStructOn.eventTimesAligned(trial);
            hold on
            for i = 1:length(trialSpks)
                line([trialSpks(i) trialSpks(i)], [0 ticklen]-ticklen*count, 'Color',[1 0.55 0.125], 'LineWidth', 1)
                hold on
            end
            count = count + 1;
            hold on
        end
        count = length(MotifOn.center)+1;
        for m = 1:length(MotifOff.center)
            trial = outStructOff.eventTrials == m;
            trialSpks = outStructOff.eventTimesAligned(trial);
            hold on
            for n = 1:length(trialSpks)
                line([trialSpks(n) trialSpks(n)], [0 ticklen]-ticklen*count, 'Color',[0.3010, 0.7450, 0.9330], 'LineWidth', 1)
                hold on
            end
            count = count + 1;
            hold on
        end
        ylim([-count*ticklen 0])
        set(gca,'xtick',[],'ytick',[],'box','off','xcolor','none','color','none','FontSize',fontsize)
        xlim(windJux)
        ylabel(string(length(Motif.center)) + ' trials')

        axis(3) = subplot(17,1,12:16);
        plot(linspace(windJux(1),windJux(2), length(countSpOff(:))), countSpOff(:), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth', 2)
        hold on
        plot(linspace(windJux(1),windJux(2), length(countSpOn(:))), countSpOn(:), 'Color', [1 0.55 0.125], 'LineWidth', 2)

        ylabel('Hz')
        xlabel('sec')
        set(gca,'box','off','color','none','TickDir','out','FontSize',fontsize)
        xlim(windJux)
        linkaxes(axis, 'x');

        cd(path + f)
        print(name + ' ' + f + ' pb trials', '-dpdf');


        %% precision on

        calcBuffer = .03;

        [onMO,psthOn] = alignSpikeTimesJuxta(clusters, [0+calcBuffer length(motif)/fs+calcBuffer], MotifOn);

        %
        %         if ~isempty(onMO.eventTimes)
        %             allLatOn = [];
        %             for i = 1:length(onMO.eventTimesAligned)
        %                 t = onMO.eventTrials(i);
        %                 compare = onMO.eventTimesAligned(onMO.eventTrials~=t);
        %                 latencies = compare-onMO.eventTimesAligned(i);
        %                 allLatOn = [allLatOn; latencies'];
        %                 clear latencies;
        %             end
        %             propOn10 = length(allLatOn(allLatOn >= -sigwin10 & allLatOn <= sigwin10))/length(allLatOn);
        %             propOn50 = length(allLatOn(allLatOn >= -sigwin50 & allLatOn <= sigwin50))/length(allLatOn);
        %         else
        %             propOn10 = nan;
        %             propOn50 = nan;
        %             allLatOn = [];
        %         end
        %
        %         precOnidx10 = (propOn10-propSpontOn10)/propSpontOn10;
        %         precOnidx50 = (propOn50-propSpontOn50)/propSpontOn50;


        %% precision off

        [offMO,psthOff] = alignSpikeTimesJuxta(clusters, [0+calcBuffer length(motif)/fs+calcBuffer], MotifOff);

        %         if ~isempty(offMO.eventTimes)
        %             allLatOff = [];
        %             for i = 1:length(offMO.eventTimesAligned)
        %                 t = offMO.eventTrials(i);
        %                 compare = offMO.eventTimesAligned(offMO.eventTrials~=t);
        %                 latencies = compare-offMO.eventTimesAligned(i);
        %                 allLatOff = [allLatOff; latencies'];
        %                 clear latencies;
        %             end
        %             propOff10 = length(allLatOff(allLatOff >= -sigwin10 & allLatOff <= sigwin10))/length(allLatOff);
        %             propOff50 = length(allLatOff(allLatOff >= -sigwin50 & allLatOff <= sigwin50))/length(allLatOff);
        %         else
        %             propOff10 = nan;
        %             propOff50 = nan;
        %             allLatOff = [];
        %         end
        %
        %         precOffidx10 = (propOff10-propSpontOff10)/propSpontOff10;
        %         precOffidx50 = (propOff50-propSpontOff50)/propSpontOff50;


        %% pre-pb prep & stats

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


        %% all or nothing corr

        binbin = .005;
        pbHalf = length(motif)/fs/2;
        buff = .5;
        windPb = [-buff pbHalf*2+calcBuffer+buff];
        fullSpec = range(windPb);
        allBins = ceil(fullSpec/binbin);

        for a = 1:4
            if a == 1
                Times = MotifOn;
            elseif a == 2
                Times = MotifOff;
            elseif a == 3
                Times = BeforeOn;
            elseif a == 4
                Times = BeforeOff;
            end
            if ~isempty(Times.center)
                series = nan(length(Times.center),allBins);
                outStruct = alignSpikeTimesJuxta(clusters, windPb, Times);
                poss = ceil(windPb(1)/binbin):ceil(windPb(2)/binbin)-1;
                for k = 1:length(Times.center)
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
                                same(k) = 1;
                            elseif comp1(k) ~= comp2(k)
                                same(k) = 0;
                            end
                            match(i,j) = mean(same);
                        end
                    end
                end
                matchBottom = tril(match,1);
                matchBottom(matchBottom == triu(matchBottom)) = nan;
                matchVal = nanmean(matchBottom(:));
            else
                matchVal = nan;
                rVal = nan;
                zVal = nan;
            end

            if a == 1
                clusters.pbAwakeMatch = matchVal;
                clusters.pbAwakeCorr = rVal;
                clusters.pbAwakeZCorr = zVal;

            elseif a == 2
                clusters.pbSleepMatch = matchVal;
                clusters.pbSleepCorr = rVal;
                clusters.pbSleepZCorr = zVal;
            elseif a == 3
                clusters.pbAwakeBeforeMatch = matchVal;
                clusters.pbAwakeBeforeCorr = rVal;
            elseif a == 4
                clusters.pbSleepBeforeMatch = matchVal;
                clusters.pbSleepBeforeCorr = rVal;
            end
        end

        %% pb stats

        f4 = figure('WindowState','normal');
        clf;
        txty = .9;
        txtx1 = .6;
        txtx2 = 1.4;
        xlim([0 2])
        ylim([-.4 1])

        text(txtx1+.1,txty+.2,'PLAYBACK STATS','FontSize',fontsize)
        text(txtx1,txty,'On','FontSize',fontsize)
        text(txtx2,txty,'Off','FontSize',fontsize)

        text(txtx1,txty-.2,num2str(length(MotifOn.start)),'FontSize',fontsize)

        onMoRates = nan(length(MotifOn.center),1);
        burstpbAwake = nan(length(MotifOn.center),1);
        cvpbAwake = nan(length(MotifOn.center),1);
        maxIFRAwake = nan(length(MotifOn.center),1);

        trialIsi = [];
        trialSpks = [];
        roundMOspks = [];
        for k = 1:length(MotifOn.center)
            trial = onMO.eventTrials == k;
            trialSpks = onMO.eventTimesAligned(trial);
            trialIsi = diff(trialSpks);
            roundMOspks = round(trialSpks(2:end)*1000,0);
            if isempty(trialIsi)
                burstpbAwake(k) = 0;
                cvpbAwake(k) = nan;
                if isempty(trialSpks)
                    onMoRates(k) = 0;
                else
                    onMoRates(k) = length(trialSpks)/(length(motif)/fs);
                end
                %  precMatOn(roundMOspks,k) = 0;
                maxIFRAwake(k) = nan;
            else
                maxIFRAwake(k) = max(1./trialIsi);
                burstpbAwake(k) = nansum(trialIsi <= .01)./length(trialIsi);
                cvpbAwake(k) = nanstd(trialIsi)/nanmean(trialIsi);
                %  precMatOn(roundMOspks,k) = 1./trialIsi;
                onMoRates(k) = length(trialSpks)/(length(motif)/fs);
            end
        end

        [onPost,psthOnPost] = alignSpikeTimesJuxta(clusters, [length(motif)/fs+calcBuffer buffer+length(motif)/fs+calcBuffer], MotifOn);
        onPostRates = nan(length(MotifOn.center),1);
        for k = 1:length(MotifOn.center)
            trial = onPost.eventTrials == k;
            trialSpks = onPost.eventTimesAligned(trial);
            trialIsi = diff(trialSpks);
            onPostRates(k) = length(trialSpks)/(length(motif)/fs);
        end
        meanMaxIFRAwake = round(nanmean(maxIFRAwake),2);
        meanPostOn = round(nanmean(onPostRates),2);
        onMoRatesMean =  round(nanmean(onMoRates),2); %round(mean(psthOnPb),2);
        onMoRatesStd =  round(nanstd(onMoRates),2); %round(std(psthOnPb),2);
        burstpbAwakeMean =  round(nanmean(burstpbAwake),3);
        cvpbAwakeMean = round(nanmean(cvpbAwake),2);
        cvpbAwakeStd = round(nanstd(cvpbAwake),2);

        text(txtx1,txty-.4,num2str(onMoRatesMean) ,'FontSize',fontsize) %+ ' +/- ' + string(onMoRatesStd)
        text(txtx1,txty-.6,num2str(meanPostOn),'FontSize',fontsize)
        text(txtx1,txty-.8,num2str(burstpbAwakeMean) ,'FontSize',fontsize) %+ ' +/- ' + string(burstpbAwakeStd)
        %   text(txtx1,txty-1,num2str(cvpbAwakeMean),'FontSize',fontsize) % + ' +/- ' + string(cvpbAwakeStd)
        %     text(txtx1,txty-1.2,num2str(round(precOnidx,3)),'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)
        text(txtx1,txty-1,num2str(round(clusters.pbAwakeMatch,3)),'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)
        %text(txtx1,txty-1.2,num2str(precOnidx50,4),'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)
        text(txtx1,txty-1.2,num2str(round(clusters.pbAwakeBeforeMatch),3),'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)

        %pb on trials
        % precMatOn = zeros(numPrecPts,length(MotifOn.center));

        % pb off trials
        %  precMatOff = zeros(numPrecPts,length(MotifOff.center));

        offMoRates = nan(length(MotifOff.center),1);
        burstpbSleep = nan(length(MotifOff.center),1);
        cvpbSleep = nan(length(MotifOff.center),1);
        maxIFRSleep = nan(length(MotifOff.center),1);

        trialIsi = [];
        trialSpks = [];
        roundMOspks = [];
        for k = 1:length(MotifOff.center)
            trial = offMO.eventTrials == k;
            trialSpks = offMO.eventTimesAligned(trial);
            trialIsi = diff(trialSpks);
            roundMOspks = round(trialSpks(2:end)*1000,0);
            if isempty(trialIsi)
                burstpbSleep(k) = 0;
                cvpbSleep(k) = nan;
                %        offMoRates(k) = 0;
                maxIFRSleep(k) = nan;
                if isempty(trialSpks)
                    offMoRates(k) = 0;
                else
                    offMoRates(k) = length(trialSpks)/(length(motif)/fs);
                end
                %      precMatOff(roundMOspks,k) = 0;
            else
                burstpbSleep(k) = nansum(trialIsi <= .01)./length(trialIsi);
                cvpbSleep(k) = nanstd(trialIsi)/nanmean(trialIsi);
                offMoRates(k) = length(trialSpks)/(length(motif)/fs);
                maxIFRSleep(k) = max(1./trialIsi);

                %    precMatOff(roundMOspks,k) = 1./trialIsi;
            end
        end

        [offPost,psthOffPost] = alignSpikeTimesJuxta(clusters, [length(motif)/fs+calcBuffer buffer+length(motif)/fs+calcBuffer], MotifOff);
        for k = 1:length(MotifOff.center)
            trial = offPost.eventTrials == k;
            trialSpks = offPost.eventTimesAligned(trial);
            trialIsi = diff(trialSpks);
            offPostRates(k) = length(trialSpks)/(length(motif)/fs);
        end
        meanMaxIFRSleep = round(nanmean(maxIFRSleep),2);
        meanPostOff = round(nanmean(offPostRates),2);
        burstpbSleepMean =  round(nanmean(burstpbSleep),3);
        cvpbSleepMean = round(nanmean(cvpbSleep),2);
        cvpbSleepStd = round(nanstd(cvpbSleep),2);
        offMoRatesMean = round(nanmean(offMoRates),2); %round(mean(psthOffPb),2);
        offMoRatesStd = round(nanstd(offMoRates),2); %round(std(psthOffPb),2);
        %         precOnidx10 = round(precOnidx10,4);
        %         precOffidx10 = round(precOffidx10,4);
        %         precOnidx50 = round(precOnidx50,4);
        %         precOffidx50 = round(precOffidx50,4);

        text(txtx2,txty-.2,num2str(length(MotifOff.start)),'FontSize',fontsize)
        text(txtx2,txty-.4,num2str(offMoRatesMean) ,'FontSize',fontsize)%+ ' +/- ' + string(offMoRatesStd)
        text(txtx2,txty-.6,num2str(meanPostOff),'FontSize',fontsize)
        text(txtx2,txty-.8,num2str(burstpbSleepMean) ,'FontSize',fontsize) %+ ' +/- ' + string(burstpbSleepStd)
        %      text(txtx2,txty-1,num2str(precOnidx10) ,'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)
        %      text(txtx2,txty-1.2,num2str(precOnidx50),'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)
        text(txtx2,txty-1,num2str(round(clusters.pbSleepMatch,3)),'FontSize',fontsize) %+ ' +/- ' + string(cvpbSleepStd)
        text(txtx2,txty-1.2,num2str(round(clusters.pbSleepBeforeMatch,3)),'FontSize',fontsize)

        % num trials
        text(0,txty-.2,'Trials','FontSize',fontsize)
        % pb fr
        text(0,txty-.4,'Rate Hz','FontSize',fontsize)
        % post rates
        text(0,txty-.6,'Post','FontSize',fontsize)
        % burstiness
        text(0,txty-.8,'Burst %','FontSize',fontsize)
        % index10
        text(0,txty-1,'Match Prec','FontSize',fontsize)
        %         % index50
        %         text(0,txty-1.2,'PrecInd50','FontSize',fontsize)
        % maxIFR
        text(0,txty-1.2,'Match Base','FontSize',fontsize)
        set(gca,'box','off','color','none','xcolor','none','ycolor','none','FontSize',fontsize)

        cd(path + f)
        print(name + ' ' + f + ' pb stats', '-dpdf');

        clusters.pbAwakeFR = onMoRatesMean;
        clusters.pbAwakeFRAll = onMoRates;

        clusters.pbAwakeBurst = burstpbAwakeMean;
        %  clusters.pbAwakeCV = cvpbAwakeMean;
        %    clusters.pbAwakePrecInd10 = precOnidx10;
        %    clusters.pbAwakePrecInd50 = precOnidx50;

        clusters.pbAwakeFRPost = meanPostOn;
        clusters.pbAwakeMaxIFR = meanMaxIFRAwake;

        clusters.pbSleepFR = offMoRatesMean;
        clusters.pbSleepFRAll = offMoRates;

        clusters.pbSleepBurst = burstpbSleepMean;
        %    clusters.pbSleepCV = cvpbSleepMean;
        %        clusters.pbSleepPrecInd10 = precOffidx10;
        %        clusters.pbSleepPrecInd50 = precOffidx50;
        clusters.pbSleepFRPost = meanPostOff;
        clusters.pbSleepMaxIFR = meanMaxIFRSleep;


        %         %% draw latency graphs
        %
        %         precAx = -length(motif)/fs:.005:length(motif)/fs;
        %         figure('WindowState','maximized');
        %         subplot(2,1,1)
        %         histogram(allLatOn,precAx,'Normalization','probability','FaceAlpha',.4);
        %         hold on;
        %         histogram(spontLatOn,precAx,'Normalization','probability','FaceAlpha',.4);
        %         xticks(-1:.05:1)
        %         xticklabels(string(-1000:50:1000))
        %         set(gca,'TickDir','out','color','none','box','off');
        %         legend('motif','spont');
        %         xlabel('milliseconds');
        %         ylabel('bin size = 5 ms');
        %         title([name + f + ' lights on = ' + string(round(precOnidx,2))])
        %
        %         subplot(2,1,2)
        %         histogram(allLatOff,precAx,'Normalization','probability','FaceAlpha',.4);
        %         hold on;
        %         histogram(spontLatOff,precAx,'Normalization','probability','FaceAlpha',.4);
        %         xticks(-1:.05:1)
        %         xticklabels(string(-1000:50:1000))
        %         set(gca,'TickDir','out','color','none','box','off');
        %         legend('motif','spont');
        %         xlabel('milliseconds');
        %         ylabel('bin size = 5 ms');
        %         title([name + f + ' lights off = ' + string(round(precOffidx,2))])
        %
        %         cd(['C:\Users\User\Dropbox (NYU Langone Health)\Int Juxta\precision\'])
        %         print([name + ' ' + f + ' precision histograms'], '-dpdf');


        %% psth significance

        %   sizePsth = length(motif)/fs+calcBuffer;
        %   [~,fullpsthOn] = alignSpikeTimesJuxta(clusters, [-sizePsth sizePsth], MotifOn);
        %  [~,fullpsthOff] = alignSpikeTimesJuxta(clusters, [-sizePsth sizePsth], MotifOff);
        %  clusters.fullPsthOn = fullpsthOn;
        %  clusters.fullPsthOff = fullpsthOff;


        %% save other params

        % clusters.spontWtFr = (clusters.spontAwakeFR*Trials(1) + clusters.spontSleepFR*Trials(2))/(Trials(1)+Trials(2));
        % clusters.pbWtFr = (clusters.pbAwakeFR*Trials(3) + clusters.pbSleepFR*Trials(4))/(Trials(3)+Trials(4));

        % clusters.ratioPbSpontOn = max(clusters.spontAwakeFR, clusters.pbAwakeFR)/min(clusters.spontAwakeFR, clusters.pbAwakeFR);
        % clusters.ratioPbSpontOff = max(clusters.spontSleepFR, clusters.pbSleepFR)/min(clusters.spontSleepFR, clusters.pbSleepFR);

        % clusters.ratioVigilance = max(clusters.spontAwakeFR, clusters.spontSleepFR)/min(clusters.spontAwakeFR, clusters.spontSleepFR);
        % clusters.numTrials = Trials;

        clusters.pbAwakeChangeOld = round((clusters.pbAwakeFR-clusters.spontAwakeFR)./(clusters.pbAwakeFR+clusters.spontAwakeFR),2);
        clusters.pbAwakeChange = round((clusters.pbAwakeFR-clusters.prePbAwakeFR)./(clusters.pbAwakeFR+clusters.prePbAwakeFR),2);
        clusters.pbAwakePostChange = round((clusters.pbAwakeFRPost-clusters.spontAwakeFR)./(clusters.pbAwakeFRPost+clusters.spontAwakeFR),2);
        clusters.pbSleepChangeOld = round((clusters.pbSleepFR-clusters.spontSleepFR)./(clusters.pbSleepFR+clusters.spontSleepFR),2);
        clusters.pbSleepChange = round((clusters.pbSleepFR-clusters.prePbSleepFR)./(clusters.pbSleepFR+clusters.prePbSleepFR),2);
        clusters.pbSleepPostChange = round((clusters.pbSleepFRPost-clusters.spontSleepFR)./(clusters.pbSleepFRPost+clusters.spontSleepFR),2);
        clusters.spontRatioVigilance = round((clusters.spontAwakeFR-clusters.spontSleepFR)./(clusters.spontAwakeFR+clusters.spontSleepFR),2);


        %% save errthang
        cd(path + f)

        save clusters.mat clusters
        %         diffRise(number) = trPkTime;
        %         diffRise70(number) = trPk70;
        %         diffWidth(number) = width;
        counting = counting + 1;
    end
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


%%

% cd(dir)
% save AwakePB.mat AwakePB;
% save SleepPB.mat SleepPB;
% save AwakeSpont.mat AwakeSpont;
% save SleepSpont.mat SleepSpont;
% save MotifsOnTimes.mat MotifsOnTimes;
% save MotifsOffTimes.mat MotifsOffTimes;
% save SpontsOnTimes.mat SpontsOnTimes;
% save SpontsOffTimes.mat SpontsOffTimes;
% save AwakeHist.mat AwakeHist;
% save SleepHist.mat SleepHist;


%% thesis eye shit for n = 1

% times = [6750 7650]/eyeRate; 
% which = spkTimes >= times(1)*fs & spkTimes <= times(2)*fs;
% spikes = spkTimes(which)/fs;
% spikes = spikes-spikes(1);
% onSpikes = spikes < 114;
% onRate = sum(onSpikes)/114;
% offSpikes = spikes > 114;
% offRate = sum(offSpikes)/(diff(times)-114);
% (onRate-offRate)/(onRate + offRate)
% 

