function [outStruct,countSp] = alignSpikeTimesSi(clusters,wind,Event,events)

    % align spikes to middle of motifs
    %   INPUT:
    %       clusters: spike times in seconds
    %       events: middle of motif times
    %       wind: window for time lock loop (before and after 0 time of events
    %       Event: output of motif finder with start and stop times
    %   OUTPUT:
    %       outStruct: output of time lock loop with fields - 
    %               eventTimes - indexing of trials
    %               eventTimes - times relative to middle of event 
    %               eventTimesAligned - times warped with recalculated warping
    %                   factor

   %events = Event.start;
    [outStruct] = timeLockLoop(clusters,events,wind);

    %% Warp re-calc

    warp = Event.warp;

    for k = 1:length(clusters)
        for n = 1:length(events)
            trial = (outStruct(k).eventTrials); 
            times = (outStruct(k).eventTimes); 
            indN = trial == n;
            times = times(indN); 
            times = times / warp(n);
            outStruct(k).eventTimesAligned(indN) = times';
        end
    end
    
    
    %% PSTH

    bw = 0.005; % sec
    start = wind(1);
    stop = wind(2);
    iter = start:bw:stop;
    countSp = nan(length(clusters),length(iter));

    for k = 1:length(clusters)
        sT = outStruct(k).eventTimesAligned;
        start = wind(1);
        for i = 1:length(iter)
            binEnd = start + bw;
            indST = (sT > start) & (sT < binEnd);
            countSp(k,i) = sum(indST)/bw/length(events);
            start = start + bw;
        end
    end

end

