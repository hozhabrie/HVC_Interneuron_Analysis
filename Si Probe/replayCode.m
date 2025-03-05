close all; clear all; clc; warning off;

fs = 30000;
birds = {'C22'; 'C24'; 'C50'; 'C52'};
addpath 'C:/Users/User/Dropbox (NYU Langone Health)/EllieScripts/SiClustering/Fcns/';

for w = 4%1:length(birds)
    bird = char(birds(w));
    disp(bird)
    path = ['C:/Users/User/Dropbox (NYU Langone Health)/SiProbe/' bird '/'];
    load([path 'data/clusters.mat']);
    load([path 'data/lightSig.mat']);
    disp(lightSeg)
    if w == 1
        spontSnips = sort([15619 15620; 15740 15741; 15957 15958; 16495 16496; 16685 16686]);
        replaySnips =  sort([15614.22 15615.22; 15745.7 15746.7; 15953.2 15934.2; 16421.9 16422.9; 16632.1 16633.1]); %16664.6 16665.6 another good one
    elseif w == 2
        spontSnips = sort([28217 28218; 28291 28292; 28645 28646; 28799  28800; 29645 29646]);
        replaySnips =  sort([28225.2 28226.2; 28341 28342; 28639.4 28640.4; 28953 28954; 29582 29583]);
    elseif w == 30
        spontSnips = sort([9309 9310; 9405 9406; 9264 9265; 9175 9176; 8998 8999]); % 8882 8883 another good one
        replaySnips =  sort([9401.8 9402.8; 9505.63 9506.63; 9284.5 9285.5; 9179.55 9180.55; 9067.1  9068.1]);
    elseif w == 4
        spontSnips = sort([7592.5 7593.5; 7739.7 7740.7; 8001 8002; 8123 8124; 8214.5 8215.5; 8333 8334; 8440 8441; 8616 8617; 9137.5 9138.5; 9213 9214; 9309 9310; 9455 9456; 9623 9624; 11282 11283; 11161 11162; 11065 11066; 10988 10989; 10338 10339; 10450.6 10451.6; 10567.8 10568.8]);
        replaySnips = sort([7730.7 7731.7; 7971.3 7972.3; 8036.7 8037.7; 8361.2 8362.2; 8528.1 8529.1; 8788.1 8789.1; 9220.5 9221.5; 9302.6 9303.6; 9509.6 9510.6; 11288.5 11289.5; 11200 11201; 11099 11100; 10952.13 10952.14; 10846.9 10847.9; 10688.7 10688.1; 10107.3 10108.3; 10284 10285; 10382.65 10383.65; 10471.5 10472.5; 10541.9 10542.9]);
    end

    goodCells = readtable([path bird '_goodSiUnits.xlsx']);
    toPrint = (goodCells{:,3});

    %%%%%%%%%%%%don't forget to account for start and stop recs

    snipLength = spontSnips(1,2)-spontSnips(1,1);
    numSnips = length(spontSnips);
    totalTime = snipLength*numSnips;

    spont.start = spontSnips(:,1);
    spont.stop = spontSnips(:,2);
    spont.center = spont.start + .5;
    spont.warp = ones(length(spontSnips),1);
    [spontSpks, ~] = alignSpikeTimesSi(clusters, [0 snipLength], spont, spont.start);

    replay.start = replaySnips(:,1);
    replay.stop = replaySnips(:,2);
    replay.warp = ones(length(replaySnips),1);
    replay.center = replay.start + .5;
    [replaySpks, ~] = alignSpikeTimesSi(clusters, [0 1], replay, replay.start);
    i = 0;
    for a = 1:length(clusters)
        if toPrint(a) == 1                  
            i = i + 1;            
            clus_sleep(i).ID = clusters(a).clusterID;
            
            trialSpks = [];
            maxIFR = [];
            trialIsis = [];
            
            for k = 1:length(replay.center)
                trial = replaySpks(a).eventTrials == k;
                trialSpks = replaySpks(a).eventTimesAligned(trial);
                replayRates(k) = length(trialSpks)/snipLength;
                trialIsis = diff(trialSpks);
                if ~isempty(trialIsis)
                    maxIFR(k) = max(1./trialIsis);
%                 elseif length(trialSpks)==1
%                     maxIFR(k) = 1;
                elseif isempty(trialSpks)
                    maxIFR(k) = nan;
                    replayRates(k) = 0;
                end
                replayBurstNum(k) = length(trialIsis(trialIsis < .01));
                replayBursty(k) = replayBurstNum(k)/length(trialIsis);
            end
            replayBursty(isnan(replayBursty)) = 0;
        %   clus_sleep(i).replayNumSpks = length(replaySpks);
            clus_sleep(i).replayFR = round(mean(replayRates),2);%length(replaySpks(i).eventTimesAligned)/totalTime;
            clus_sleep(i).replayBurstProp = round(nanmean(replayBursty),2);
            clus_sleep(i).replayMaxIFR = round(nanmean(maxIFR),2);
            if isnan(clus_sleep(i).replayMaxIFR)
                clus_sleep(i).replayMaxIFR = 0;
            end
            
            %    clus_sleep(i).replayBurstMeanNum = round(nanmean(replayBurstNum,2));

            
            trialSpks = [];
            maxIFR = [];
            trialIsis = [];
            
            for k = 1:length(spont.center)
                trial = spontSpks(a).eventTrials == k;
                trialSpks = spontSpks(a).eventTimesAligned(trial);                
                spontRates(k) = length(trialSpks)/snipLength;
                trialIsis = diff(trialSpks);
                if ~isempty(trialIsis)
                    maxIFR(k) = max(1./trialIsis);
%                 elseif length(trialSpks)==1
%                     maxIFR(k) = 1;
                elseif isempty(trialSpks)
                    maxIFR(k) = nan;
                    spontRates(k) = 0;
                end
                spontBurstNum(k) = length(trialIsis(trialIsis < .01));
                spontBursty(k) = spontBurstNum(k)/length(trialIsis);
            end
            spontBursty(isnan(spontBursty)) = 0;
      %     clus_sleep(i).spontNumSpks = length(spontSpks);
            clus_sleep(i).nonReplayFR = round(mean(spontRates),2);%length(spontSpks(i).eventTimesAligned)/totalTime;
            clus_sleep(i).nonReplayBurstProp = round(nanmean(spontBursty),2);
       %     clus_sleep(i).spontBurstMeanNum = round(nanmean(spontBurstNum,2));
            clus_sleep(i).nonReplayMaxIFR = round(nanmean(maxIFR),2);            
            if isnan(clus_sleep(i).nonReplayMaxIFR)
                clus_sleep(i).nonReplayMaxIFR = 0;
            end
            clus_sleep(i).nonReplayOverReplay = round(clus_sleep(i).nonReplayFR/clus_sleep(i).replayFR,2);
        end
    end

    cd([path 'data/'])
    save clus_sleep.mat clus_sleep
end
