function cutDat(F,ampFile,samplingRate,nChan,times,channels)
%Input
% F: folder with the data
% ampFile:  This will be written as cut_FILE
% samplingRate: Fs
% nChan:   # of channels in your probe/amplifier
% times: matrix of N x 2, where N is the number of segments you want to
% splice together, and each column is the start and stop IN SECONDS. leave
% empty to use all
% channels: channels to cut. leave empty to use all. (i.e 1:64)
%Output: none, but there will be new files in the folder
maxPiece=10*60;%analyze X seconds at a time. here 10 minutes. use for RAM control, increase for minor time speedups
if isempty(times)&&isempty(channels)
    error('uggggghhhh')
end
if isempty(channels)
    channels=1:nChan;
end
if isempty(times)
    fileinfo = dir([F,ampFile]);
    stop = fileinfo.bytes/(nChan * 2)/samplingRate;
    times=[0,stop];
end
%%
ampFileCut = fopen([F, 'cut_', ampFile], 'w');
for t = 1:size(times,1)
    start = times(t,1);
    stop = times(t,2);
    indCut=0;
    right = start;
    fprintf(['Starting Amp Time Segment #' num2str(t) '/' num2str(size(t,1)) ': '])
    while right < stop
        left = start+indCut*maxPiece;
        right = min(stop,left+maxPiece);
        dataChunk = LoadBinary([F,ampFile],'nChannels',nChan,'channels',channels,'start',left,'duration',right-left,'frequency',samplingRate);%SWAP TO ALL SHANKS?
        fwrite(ampFileCut, dataChunk', 'int16');
        indCut=indCut+1;
        fprintf([num2str(indCut) '/' num2str(ceil((stop-start)/maxPiece)) ', '])
    end
    disp('done')
end
fclose(ampFileCut);
disp('done')