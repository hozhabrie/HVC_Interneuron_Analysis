function [structOut] = timeLockLoop(structIn,events,win)
% Runs function timeLock.m in a loop using input structure structInput and produces new structure with aligned spike times 
%   see timeLock.m
%       structInp: define struct as structure containing spikeTimes for each cluster of interest
%       structOut: structOutput name
%       win: range over which time lock extracts spikes (with event time = 0)
%       events: event times defined in same time scale as spike times

for k = 1:length(structIn)
    spikes = (structIn(k).spikeTimes); % spike times in seconds (or at least consistent with event times)
    [times,trial]=timeLock(spikes,events,win);
    structOut(k).eventTimes = times;
    structOut(k).eventTrials = trial;
end

end

