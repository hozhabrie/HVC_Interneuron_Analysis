function [MotifFixed] = shiftMotifTimes(path, fs, Motif)
    load([path 'MotifFolder\playbackClipRange.mat']);
    MotifFixed = Motif;
    MotifFixed.start = MotifFixed.start + playbackClipRange(1)/fs;
    MotifFixed.center = MotifFixed.center + playbackClipRange(1)/fs;
    MotifFixed.stop = MotifFixed.stop + playbackClipRange(1)/fs;
end