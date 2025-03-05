% Ensure the audio file exists before attempting to read
if ~isfile(audioF)
    error('Audio file not found: %s', audioF);
end

% Read the full audio file
audioFull = audioread(audioF);

% Copy motif event structure
Event = Motif;

% Plot individual event waveforms
figure;
for i = 1:length(Event.center)
    % Ensure indices are within valid range
    durStart = round(Event(i).start * fs);
    durStop = round(Event(i).stop * fs);
    
    if durStart < 1 || durStop > length(audioFull)
        warning('Skipping event %d: Index out of bounds.', i);
        continue;
    end

    % Extract duration range
    dur = durStart:durStop;
    
    % Plot waveform with vertical offset
    plot(audioFull(dur) - 2 * i);
    hold on;
end
title('Motif Events');
xlabel('Time (samples)');
ylabel('Amplitude');
grid on;

% Calculate event duration
duration = Event.stop - Event.start;

% Recalculate warp factor
warpRecalc = duration ./ (length(audioFull) / fs);

% Display warp recalculation results
disp([Event.warp, warpRecalc]);
