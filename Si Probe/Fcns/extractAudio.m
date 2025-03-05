function [audio] = extractAudio(F)
% extracts audio from analog.dat and writes .lfp file for use in neuroscope
% Input: F = folder location
% Output: audio = extracted normalized vector from analogIn.dat

%% extract audio from analog signal
    chansAn = 2;
    nAmpAn = 1;
    filePath = F;

    fileinfo = dir([F,'analogin.dat']);

    num_channels = chansAn; % ADC input info from header file
    num_samples = fileinfo.bytes/(num_channels*2); % uint16 = 2 bytes
    fid = fopen([filePath,'analogin.dat'], 'r');
    v = fread(fid, [num_channels, num_samples], 'uint16');
    %v = fread(fid, 'uint16');
    fclose(fid);
    v = v * 0.000050354; % Magic number - convert to volts
    analogIn = v(nAmpAn,:);
    %clear v
    audio = analogIn; %(300*60*3e4:end);
    clear analogIn

    audio = audio - mean(audio);
    audio = audio/max(abs(audio));

%% write .lfp for use in neuroscope
    minOld = min(audio);
    maxOld = max(audio);
    minNew = -2^15+1;
    maxNew = 2^15;

    audioLin = (audio- minOld) /(maxOld - minOld)*(maxNew - minNew) + minNew;


    fidn = fopen([F,'audio.lfp'], 'w');
    fwrite(fidn, audioLin, 'int16');
    fclose(fidn);
    %% write to .wav    
    audiowrite([F,'/audio/audioFull.wav'],audio,3e4)
end