function writeNeuroscopeClass(clusters,samplingRate,fileBaseName)
% takes cluster IDs and spike times from loadKS (or variant) and writes two file types:
% fileBaseName.clu.n = list of cluster IDs in chronological order
%                       corresponding to spike times; first entry = total
%                       number of clusters on class (n = class number)
% fileBaseName.res.n = list of spike times in chronological order;
%                       (n = class number)

numClust = length(clusters);
%classID = horzcat(clusters(:).class);
clusterID = double(horzcat(clusters(:).clusterID));
numclass = unique(clusterID);
classID = clusterID;

for n = numclass
    indN = find(classID == n);
    classSpikeTimes = [];
    classClusterID = [];
    spikeTimes = [];
    clustIDvec = [];
    classClusterIDSorted = [];
    classSpikeTimesSorted = [];
    for k = indN
        spikeTimes = clusters(k).spikeTimes.*samplingRate;
        clustIDvec = ones(length(spikeTimes),1);
        clustIDvec = clustIDvec.*clusterID(k);
        classSpikeTimes = [classSpikeTimes;spikeTimes]; 
        classClusterID = [classClusterID;clustIDvec];
    end
    [classSpikeTimesSorted,indsT] = sort(classSpikeTimes);
    classClusterIDSorted = classClusterID(indsT);    
    classClusterIDSorted = [length(indN);classClusterIDSorted];
    try
        dlmwrite([fileBaseName,'.clu.',num2str(n)],classClusterIDSorted,'precision',11,'newline','pc');
        dlmwrite([fileBaseName,'.res.',num2str(n)],classSpikeTimesSorted,'precision',11,'newline','pc');
    catch
        keyboard
    end
    
end

end

