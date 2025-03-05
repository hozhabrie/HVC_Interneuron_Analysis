function [times, trial] = timeLock(spikes, events, range)
% timeLock - Aligns spike times relative to event times within a specified range.
%
% Inputs:
%   spikes - Array: Spike timestamps (in seconds)
%   events - Array: Event timestamps (in seconds)
%   range  - Array: Time window [start, stop] relative to events
%
% Outputs:
%   times  - Array: Spike times relative to event times
%   trial  - Array: Trial indices corresponding to each spike time

% Preallocate arrays for efficiency
times = [];
trial = [];

% Iterate over events and align spikes
for e = 1:length(events)
    tN = spikes - events(e); % Compute spike times relative to the event
    validIdx = (tN >= range(1)) & (tN <= range(2)); % Keep spikes within range
    
    if any(validIdx)
        times = [times; tN(validIdx)];
        trial = [trial; repmat(e, sum(validIdx), 1)];
    end
end

% Sort spike times and reorder trial indices accordingly
[times, inds] = sort(times);
trial = trial(inds);
end
