function [structOut] = timeLockLoop(structIn, events, win)
% timeLockLoop - Aligns spike times relative to event times for multiple clusters.
%
% Inputs:
%   structIn - Struct array: Contains spikeTimes for each cluster
%   events   - Array: Event times (same time scale as spike times)
%   win      - Array: Time window [start, stop] for time-locking (event time = 0)
%
% Output:
%   structOut - Struct array with the following fields:
%               * eventTimes  - Spike times relative to event occurrences
%               * eventTrials - Trial indices corresponding to each spike time

% Preallocate structOut for efficiency
structOut = struct('eventTimes', [], 'eventTrials', []);

% Loop through each cluster in structIn
for k = 1:length(structIn)
    spikes = structIn(k).spikeTimes; % Extract spike times (consistent with event times)
    [times, trial] = timeLock(spikes, events, win); % Call timeLock function

    % Store aligned spike times in structOut
    structOut(k).eventTimes = times;
    structOut(k).eventTrials = trial;
end

end
