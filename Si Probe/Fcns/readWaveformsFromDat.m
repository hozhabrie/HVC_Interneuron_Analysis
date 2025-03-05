function [meanWF, allWF] = readWaveformsFromDat(datFilename, chansInDat, takeChans, sampleTimes, window, nToRead)

FileInf = dir(datFilename);
nSampsInDat = (FileInf.bytes/chansInDat/2);
rawData = memmapfile(datFilename, 'Format', {'int16', [chansInDat, nSampsInDat], 'x'});

if isempty(nToRead)
    theseTimes = sampleTimes;
    nToRead = length(sampleTimes);
else
    if nToRead>length(sampleTimes)
        nToRead = length(sampleTimes);
    end
    q = randperm(length(sampleTimes));
    theseTimes = sampleTimes(sort(q(1:nToRead)));
end

theseTimes = theseTimes(theseTimes>-window(1) & theseTimes<nSampsInDat-window(2)-1);
nToRead = numel(theseTimes);

allWF = zeros(length(takeChans), diff(window)+1, nToRead);

for i=1:nToRead
    allWF(:,:,i) = double(rawData.Data.x(takeChans,theseTimes(i)+window(1):theseTimes(i)+window(2)))*.195;
end
allWF=squeeze(allWF);
%meanWF = median(allWF,length(size(allWF)));
meanWF = mean(allWF,length(size(allWF)));

