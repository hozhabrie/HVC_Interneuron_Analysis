function writeNeuroscopeEvents(starts,outFileName)
% takes event start times from 'starts' (in seconds) and writes file with name defined
%           by fileBaseName (i.e. stimID1.evt.ext) which takes start times
%           in milliseconds (conversion done in function, input must be in
%           seconds)

clusterIDSorted = [];
eventTimesSorted = [];
eventTimesA = []; 
clusterIDA = [];

for k = 1:length(starts)
    eventTimes = [];
    clustIDvec = [];
    eventTimes = starts(k)*1e3;
    eventTimes = round(eventTimes);
    clustIDvec = ones(length(eventTimes),1);
%     clustIDvec = (clustIDvec.*clusterID(k))';
    eventTimesA = [eventTimesA eventTimes]; 
    clusterIDA = [clusterIDA clustIDvec];
end
[eventTimesSorted,indEvt] = sort(eventTimesA);
clusterIDSorted = clusterIDA(indEvt); 

outFile = fopen(outFileName,'w');
for i = 1:length(eventTimesSorted)
    fprintf(outFile,'%.2f playback%i\n',[eventTimesSorted(i)', clusterIDSorted(i)]); % in ms
end
fclose(outFile);

end
