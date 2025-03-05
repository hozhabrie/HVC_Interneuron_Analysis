%% FUNCTION SCRIPT: juxFunctions.m
% This script contains multiple utility functions for processing and analysis of juxtacellular recordings of single cells.


%% Waveform functions
%Create normalized spike wf plot
function plotWaveform(xAxisWaveForm, randWindow, meanWF, fontsize, name, f, path)
    % Plot waveform figure
    f1 = figure('WindowState', 'maximized');
    clf;
    hold on;
    
    % Plot all waveforms in gray
    for i = 1:size(randWindow, 2)
        plot(xAxisWaveForm, randWindow(:, i), 'Color', [0.5 0.5 0.5]);
    end
    
    % Plot mean waveform in black
    plot(xAxisWaveForm, meanWF, 'LineWidth', 4, 'Color', 'k');
    
    % Formatting
    set(gca, 'box', 'off', 'TickDir', 'out', 'FontSize', fontsize, 'Color', 'none');
    xlabel('ms'); 
    ylabel('mV');
    title(name + " " + f);

    % Save plot
    cd(fullfile(path, f));
    print(name + " " + f + " wf", '-dpdf');
end

% Calculate waveform metrics
function [wfNorm, trPkTime, trPk70, width, wfMult70] = computeWaveformParams(meanWF, fs, deflect, number)
    meanSubWav = meanWF - mean(meanWF);
    maxAmp = max(abs(meanSubWav));
    wfNorm = meanSubWav ./ maxAmp;
    
    factor = 1000;
    x = 1:length(wfNorm);
    xi = linspace(1, length(wfNorm), length(wfNorm) * factor);
    normWavInterp = round(interp1(x, wfNorm, xi, 'spline'), 5);

    % Adjust polarity if necessary
    if deflect(number)
        normWavInterp = normWavInterp * (-1);
    end
    
    % Find minimum and peak values
    trVal = min(normWavInterp);
    trPos = find(islocalmin(normWavInterp, 'MaxNumExtrema', 1));
    pkVal = max(normWavInterp(trPos:end));
    pkPos = (trPos - 1) + find(islocalmax(normWavInterp(trPos:end), 'MaxNumExtrema', 1));
    trPkTime = (pkPos - trPos) / fs * 1000 / factor;

    % Compute 70% peak rise time
    pkVal70 = trVal + 0.7 * range([pkVal trVal]);
    closeEnuf = find(abs(normWavInterp - pkVal70) < 0.0015);
    inBet = closeEnuf(closeEnuf < pkPos & closeEnuf > trPos);
    pkPos70 = inBet(islocalmin(abs(normWavInterp(inBet) - pkVal70), 'MaxNumExtrema', 1));
    trPk70 = (pkPos70 - trPos) / fs * 1000 / factor;

    % Compute waveform width at 50% peak
    pkVal50 = trVal + 0.5 * range([pkVal trVal]);
    closeEnuf = find(abs(normWavInterp - pkVal50) < 0.0015);
    beforePk = closeEnuf(closeEnuf < pkPos);
    pkPos50 = beforePk(islocalmin(abs(normWavInterp(beforePk) - pkVal50), 'MaxNumExtrema', 2));
    width = diff(pkPos50) / fs * 1000 / factor;

    % Compute wfMult70
    wfMult70 = (trPkTime - trPk70) .* (trPkTime ./ trPk70);
end



%% Spontaneous activity functions

% Adjust spont segments to avoid playback trials
function fixedSegs = adjustSegment(seg, pbTrial, fixedSegs)
    for j = 1:size(pbTrial,1)
        while seg(2) > pbTrial(j,2) && seg(1) < pbTrial(j,1)
            fixedSegs = [fixedSegs; seg(1) pbTrial(j,1)];
            seg = [pbTrial(j,2) seg(2)];
        end
    end
    % Only append if the segment remains valid
    if seg(2) > seg(1)
        fixedSegs = [fixedSegs; seg];
    end
end

% Process spontaneous segments: separate long segments, create short ones
function spontSegs = processSpontaneousSegments(fixedSegs)
    % Calculate segment lengths
    segLengths = fixedSegs(:, 2) - fixedSegs(:, 1);
    
    % Identify long and short segments
    longSegs = segLengths >= 30;  % Segments >= 30 seconds
    shortSegs = fixedSegs(segLengths >= 5 & segLengths < 30, :); % Segments between 5 and 30 sec

    % Process long segments by dividing them into smaller chunks
    for i = find(longSegs)'  % Loop over long segments
        left = floor(segLengths(i) / 20); % Number of 20-sec chunks
        shortSeg = zeros(left+1, 2);  % Preallocate for speed
        shortSeg(:, 1) = (fixedSegs(i, 1):20:fixedSegs(i, 2))'; % Ensure column vector

        % Assign segment endpoints
        if left >= 1
            shortSeg(1:left, 2) = shortSeg(2:left+1, 1);
        end
        shortSeg(end, 2) = fixedSegs(i, 2);  % Last segment reaches end of original
        
        % Append processed short segments
        shortSegs = [shortSegs; shortSeg];
    end

    % Sort segments and ensure all are at least 5 sec long
    spontSegs = sortrows(shortSegs(shortSegs(:, 2) - shortSegs(:, 1) >= 5, :));
end

% Randomly select spontaneous trials from valid segments
function Spont = extractSpontaneousTrials(spontSegs, spontNumTrials, number)
    rng(number);
    whichSpont = randperm(min(size(spontSegs,1)));

    Spont.start = zeros(min(spontNumTrials, length(whichSpont)), 1);
    Spont.stop = zeros(size(Spont.start));
    Spont.center = zeros(size(Spont.start));

    for i = 1:length(Spont.start)
        a = ceil(spontSegs(whichSpont(i), 1));
        b = ceil(spontSegs(whichSpont(i), 2)) - 5;
        rng(number);
        snip = a + (b - a) * rand(1,1);
        Spont.start(i) = snip;
        Spont.stop(i) = snip + 5;
        Spont.center(i) = snip + 2.5;
    end
    Spont.warp = ones(size(Spont.start));
end

%Handle empty spont struct
function Spont = emptySpontaneousStruct()
    Spont.start = [];
    Spont.stop = [];
    Spont.center = [];
    Spont.warp = [];
end
       
%Plot spontaneous rasters for on and off conditions
function plotSpontaneousActivity(spontStructOn, spontStructOff, SpontOn, SpontOff, windSpont, onTri, offTri, totTri, path, f, name, fontsize)
    f6 = figure;
    clf;

    % Plot spontaneous activity (lights on)
    if ~isempty(spontStructOn.eventTimes)
        subplot(totTri, 1, 1:onTri+1)
        plotEventRasters(spontStructOn, SpontOn, windSpont)
        title('spontaneous activity')
        ylabel([num2str(onTri) ' Trials'])
        set(gca, 'box', 'off', 'color', 'none', 'TickDir', 'out', ...
                 'xcolor', [1 0.55 0.125], 'ycolor', [1 0.55 0.125], 'FontSize', fontsize)
    end

    % Plot spontaneous activity (lights off)
    if ~isempty(spontStructOff.eventTimes)
        subplot(totTri, 1, onTri+2:totTri)
        plotEventRasters(spontStructOff, SpontOff, windSpont)
        ylabel([num2str(offTri) ' Trials'])
        set(gca, 'box', 'off', 'color', 'none', 'TickDir', 'out', ...
                 'xcolor', [0.3010, 0.7450, 0.9330], 'ycolor', [0.3010, 0.7450, 0.9330], 'FontSize', fontsize)
    end

    % Set X-axis limits
    xlim(windSpont)
    xticks(windSpont(1):windSpont(2));
    xticklabels(string(0:5));
    xlabel('sec')

    % Save figure
    cd(fullfile(path, f))
    print(name + " " + f + " spont", '-dpdf');
end

%Plot spontaneous ISIs
function plotISIHistogram(on_isi_sort, off_isi_sort, edges, labels, nameLabels, fontsize, path, f, name)
    f5 = figure('WindowState', 'normal');
    
    histogram(on_isi_sort, 10.^edges, 'Normalization', 'probability', ...
              'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', [1 0.55 0.125]);
    set(gca, 'xscale', 'log', 'TickDir', 'out');
    
    hold on;
    histogram(off_isi_sort, 10.^edges, 'Normalization', 'probability', ...
              'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', [0.3010, 0.7450, 0.9330]);
    
    set(gca, 'xscale', 'log', 'box', 'off', 'color', 'none', 'FontSize', fontsize, 'TickDir', 'out');
    title('spontaneous ISIs');
    ylabel('ISI counts (prob.)');
    legend('lights on', 'lights off', 'Location', 'northeast');
    
    xticks(labels);
    set(gca, 'XTickLabels', nameLabels);
    xlabel('ms');

    % Save figure
    cd(fullfile(path, f));
    print(name + " " + f + " isis", '-dpdf');
end



%% Playback activity functions

% Plot playback rasters
function plotPlaybackActivity(windJux, paddedAudio, fontsize, MotifOn, MotifOff, ...
                              outStructOn, outStructOff, ticklen, countSpOn, countSpOff, path, f, name)
    f2 = figure('WindowState', 'normal');

    % Plot playback activity waveform
    axis(1) = subplot(17, 1, 1:3);
    plot(linspace(windJux(1), windJux(2), length(paddedAudio)), paddedAudio, 'k');
    set(gca, 'xtick', [], 'ytick', [], 'box', 'off', 'xcolor', 'none', ...
             'ycolor', 'none', 'color', 'none', 'FontSize', fontsize);
    xlim(windJux);
    title('playback activity');
    hold on;

    % Raster plot for Motif On
    axis(2) = subplot(17, 1, 4:11);
    count = 1;
    for k = 1:length(MotifOn.center)
        trial = outStructOn.eventTrials == k;
        trialSpks = outStructOn.eventTimesAligned(trial);
        for i = 1:length(trialSpks)
            line([trialSpks(i) trialSpks(i)], [0 ticklen] - ticklen * count, 'Color', [1 0.55 0.125], 'LineWidth', 1);
        end
        count = count + 1;
    end

    % Raster plot for Motif Off
    count = length(MotifOn.center) + 1;
    for m = 1:length(MotifOff.center)
        trial = outStructOff.eventTrials == m;
        trialSpks = outStructOff.eventTimesAligned(trial);
        for n = 1:length(trialSpks)
            line([trialSpks(n) trialSpks(n)], [0 ticklen] - ticklen * count, 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth', 1);
        end
        count = count + 1;
    end
    ylim([-count * ticklen 0]);
    set(gca, 'xtick', [], 'ytick', [], 'box', 'off', 'xcolor', 'none', 'color', 'none', 'FontSize', fontsize);
    xlim(windJux);
    ylabel(string(length(MotifOn.center) + length(MotifOff.center)) + " trials");

    % PSTH plot
    axis(3) = subplot(17, 1, 12:16);
    plot(linspace(windJux(1), windJux(2), length(countSpOff(:))), countSpOff(:), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth', 2);
    hold on;
    plot(linspace(windJux(1), windJux(2), length(countSpOn(:))), countSpOn(:), 'Color', [1 0.55 0.125], 'LineWidth', 2);

    ylabel('Hz');
    xlabel('sec');
    set(gca, 'box', 'off', 'color', 'none', 'TickDir', 'out', 'FontSize', fontsize);
    xlim(windJux);
    linkaxes(axis, 'x');

    % Save figure
    cd(fullfile(path, f));
    print(name + " " + f + " pb trials", '-dpdf');
end



%% Pan-behavior metrics functions

% Calculate spike train metrics
function result = computeSpikeMetrics(metricType, Motif, MO, motif, fs, varargin)
    % Generalized function to compute various spike metrics
    % metricType: 'rate', 'burst', 'cv', 'maxIFR', 'postRate'
    % Motif, MO: structures containing spike times and event trials
    % motif: motif duration
    % fs: sampling frequency
    % varargin: additional arguments (for 'postRate' metric)

    result = nan(length(Motif.center), 1);  % Initialize result array

    for k = 1:length(Motif.center)
        trial = MO.eventTrials == k;
        trialSpks = MO.eventTimesAligned(trial);
        trialIsi = diff(trialSpks);  % Compute interspike intervals (ISI)

        switch metricType
            case 'rate'
                % Compute firing rate
                result(k) = length(trialSpks) / (length(motif) / fs);

            case 'burst'
                % Compute burst rate: fraction of ISIs <= 10ms (0.01 sec)
                if isempty(trialIsi)
                    result(k) = 0;
                else
                    result(k) = nansum(trialIsi <= 0.01) / length(trialIsi);
                end

            case 'cv'
                % Compute coefficient of variation (CV) of ISIs
                if isempty(trialIsi)
                    result(k) = nan;
                else
                    result(k) = nanstd(trialIsi) / nanmean(trialIsi);
                end

            case 'maxIFR'
                % Compute max instantaneous firing rate (Hz)
                if isempty(trialIsi)
                    result(k) = nan;
                else
                    result(k) = max(1 ./ trialIsi);
                end

            case 'postRate'
                % Compute post-playback firing rate
                if length(varargin) < 3
                    error('computeSpikeMetrics: Missing arguments for postRate computation');
                end
                calcBuffer = varargin{1};  % Post-playback buffer
                clusters = varargin{2};    % Cluster structure
                buffer = varargin{3};      % Ensure buffer is passed

                % Align spike times to post-playback window
                [postStruct, ~] = alignSpikeTimesJuxta(clusters, ...
                    [length(motif) / fs + calcBuffer, buffer + length(motif) / fs + calcBuffer], ...
                    Motif);

                trial = postStruct.eventTrials == k;
                trialSpks = postStruct.eventTimesAligned(trial);

                % Compute post-rate
                result(k) = length(trialSpks) / (length(motif) / fs);

            otherwise
                error('Invalid metricType specified.');
        end
    end

    % If computing postRate, return the mean post-playback rate
    if strcmp(metricType, 'postRate')
        result = round(nanmean(result), 2);
    end
end

% Plot metrics stats
function plotMetricsStats(isPlaybackStats, fixedSegsOn, fixedSegsOff, clusters, SpontOn, SpontOff, path, f, name, fontsize)
    fStat = figure('WindowState', 'normal');
    
    txty = 0.9;
    txtx1 = 0.6;
    txtx2 = 1.4;
    xlim([0 2]);
    ylim([-0.4 1]);
    
    set(gca, 'box', 'off', 'color', 'none', 'xcolor', 'none', 'ycolor', 'none', 'FontSize', fontsize);
    
    % Title Selection
    if isPlaybackStats
        text(txtx1 + 0.1, txty + 0.2, 'PLAYBACK STATS', 'FontSize', fontsize);
    else
        text(txtx1 - 0.15, txty + 0.2, 'SPONTANEOUS STATS', 'FontSize', fontsize);
    end
    
    % Column Labels
    text(txtx1, txty, 'On', 'FontSize', fontsize);
    text(txtx2, txty, 'Off', 'FontSize', fontsize);
    
    % Row Labels
    text(0, txty - 0.2, 'Trials', 'FontSize', fontsize);
    text(0, txty - 0.4, 'Rate (Hz)', 'FontSize', fontsize);
    if isPlaybackStats
        text(0, txty - 0.6, 'Post', 'FontSize', fontsize);
        text(0, txty - 0.8, 'Burst %', 'FontSize', fontsize);
        text(0, txty - 1, 'Match Prec', 'FontSize', fontsize);
        text(0, txty - 1.2, 'Match Base', 'FontSize', fontsize);
    else
        text(0, txty - 0.6, 'Burst %', 'FontSize', fontsize);
        text(0, txty - 0.8, 'CV', 'FontSize', fontsize);
    end

    % Save Figure
    cd(fullfile(path, f));

    if isPlaybackStats
        print(name + " " + f + " pb stats", '-dpdf');
    else
        print(name + " " + f + " spont stats", '-dpdf');
    end
end

