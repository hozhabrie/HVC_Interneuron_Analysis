function [an1,an2] = extractAnalog(F)
% extracts audio from analog.dat and writes .lfp file for use in neuroscope
% Input: F = folder location
% Output: audio = extracted normalized vector from analogIn.dat

%% extract audio from analog signal
    chansAn = 2;
    filePath = F;

    fileinfo = dir([F,'analogin.dat']);

    num_channels = chansAn; % ADC input info from header file
    num_samples = fileinfo.bytes/(num_channels*2); % uint16 = 2 bytes
    fid = fopen([filePath,'analogin.dat'], 'r');
    v = fread(fid, [num_channels, num_samples], 'uint16');
    %v = fread(fid, 'uint16');
    fclose(fid);
    v = v * 0.000050354; % Magic number - convert to volts
    nAudio = 1;
    an1 = v(nAudio,:);
    an1 = an1 - mean(an1);
    an1 = an1/max(abs(an1));

%% write .lfp for use in neuroscope
    minOld = min(an1);
    maxOld = max(an1);
    minNew = -2^15+1;
    maxNew = 2^15;
    audioLin = (an1 - minOld) /(maxOld - minOld)*(maxNew - minNew) + minNew;
%     fidn = fopen([F,'audio.lfp'], 'w');
%     fwrite(fidn, audioLin, 'int16');
%     fclose(fidn);

%% write to .wav    
%     audiowrite([F,'audioFull.wav'],an1,3e4);
%     clear audio;
%     
%% 
    nLight = 2;
    an2 = v(nLight,:);
   % clear v;
    an2 = an2 - mean(an2);
    an2 = an2/max(abs(an2));
    diffLight = find(round(diff(an2))==1);
    lightsOff = diffLight(1)/3e4;
    lightsOn = diffLight(end)/3e4;
    save([F,'lightSig.mat'], 'lightsOn', 'lightsOff');

%% write .lfp for use in neuroscope
    minOld = min(an2);
    maxOld = max(an2);
    minNew = -2^15+1;
    maxNew = 2^15;

    an2 = (an2 - minOld) /(maxOld - minOld)*(maxNew - minNew) + minNew;


    fidn = fopen([F,'lightSig.lfp'], 'w');
    fwrite(fidn, an2, 'int16');
    fclose(fidn);
end