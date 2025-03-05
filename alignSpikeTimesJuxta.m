function [outStruct, countSp] = alignSpikeTimesJuxta(clusters, wind, Event)
% alignSpikeTimesJuxta - Aligns spike times to the middle of detected motifs.
%
% Inputs:
%   clusters - Array: Spike times in seconds
%   wind     - Array: Window for time-lock loop (before and after 0 time of events)
%   Event    - Struct: Output of motif finder containing start, stop, and warp times
%
% Outputs:
%   outStruct - Struct array with the following fields:
%               * eventTrials: Index of trials
%               * eventTimes: Times relative to middle of event
%               * eventTimesAligned: Times warped with recalculated warp factor
%   countSp   - Matrix: Peri-stimulus time histogram (PSTH) spike counts

% Extract event start times
events = Event.start;

% Perform time-locking
outStruct = timeLockLoop(clusters, events, wind);

% Ensure eventTimesAligned exists even if no events are detected
if isempty(events)
    for k = 1:length(clusters)
        outStruct(k).eventTimesAligned = [];
    end
end

%% Warp Recalculation
warp = Event.warp;

for k = 1:length(clusters)
    for n = 1:length(events)
        trial = outStruct(k).eventTrials; 
        times = outStruct(k).eventTimes; 
        
        % Get indices of trials matching event n
        indN = trial == n;
        
        if any(indN)  % Ensure there are matching trials
            times = times(indN) / warp(n);
            outStruct(k).eventTimesAligned(indN) = times';
        end
    end
end

%% Peri-Stimulus Time Histogram (PSTH)
bw = 0.005; % Bin width in seconds
start = wind(1);
stop = wind(2);
inter = start:bw:stop;
countSp = nan(length(clusters), length(inter));

if ~isempty(outStruct)
    sT = outStruct(1).eventTimes; 
    start = wind(1);

    for i = 1:length(inter)
        binEnd = start + bw;
        indST = (sT > start) & (sT < binEnd);
        countSp(i) = sum(indST) / bw / length(Event.center);
        start = start + bw;
    end
end

% countSp = movmean(countSp, 4); % Uncomment if smoothing is needed
end
