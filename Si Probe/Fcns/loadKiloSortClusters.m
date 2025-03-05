function clusters = loadKiloSortClusters(SiProbePath, SiProbeSamplingInterval)

%Load in from KiloSort
%try%regular mode
    load(fullfile(SiProbePath,'KS_Output.mat'));
    clust_id=readNPY(fullfile(SiProbePath,'spike_clusters.npy'));
    clust_group=importdata(fullfile(SiProbePath,'cluster_group.tsv'));
% catch%developer mode
%     load(fullfile(SiProbePath,'KS_Output.mat'));
%     clust_id=readNPY(fullfile(SiProbePath,'spike_clusters.npy'));
%     clust_group=importdata(fullfile(SiProbePath,'cluster_group.tsv'));
% end
ops.chanMap= 'C:\Users\User\Dropbox (NYU Langone Health)\EllieScripts\SiClustering\Extra/chanMaps/128_5_integrated.mat';
load(ops.chanMap,'connected', 'xcoords', 'ycoords','kcoords');
% cOrder=[AnatGroup{:}]+1;
% clusts=double(unique(clust_id)')';
clusterLabels = cell(length(clust_group) - 1, 2);
for i = 1:length(clust_group) - 1
    line = textscan(clust_group{i+1},'%d %s');
    clusterLabels{i, 1} = line{2};%label assigned to the cluster
    clusterLabels{i, 2} = line{1};%cluster id 
end

% ADD: load inverse whitening matrix from KiloSort output
WInv = rez.WrotInv;
unwhitened_templates = zeros(size(rez.Wraw));
for i = 1:length(rez.Wraw(1, 1, :))
    template_unwhitened = WInv * rez.Wraw(:,:,i);
    unwhitened_templates(:, :, i) = template_unwhitened(:, :);
end

goodClusterNames = strcmp([clusterLabels{:,1}],'good');%only take good
goodSpikeIndices = ismember(clust_id, [clusterLabels{goodClusterNames,2}]);
goodClusters = unique(clust_id(goodSpikeIndices));
clusters=struct('spikeTimes',[],'clusterID',[],'coordinates',[],'shank',[]);


%%
for i = [20 34 45 55]%55-1 %34-2 %20-3 %45-4 1:length(goodClusters)
    disp(i)
    clusterOfInterest = goodClusters(i);
    clusterSpikeIndices = ismember(clust_id, clusterOfInterest);
    clusterSpikeTimes = SiProbeSamplingInterval*rez.st3(clusterSpikeIndices,1);
    
    %Find (channel) location of max waveform
    tmpSpike = clust_id==clusterOfInterest;
    origC = unique(rez.st3(tmpSpike,2));
%     meanTemplate = mean(rez.Wraw(:,:,origC),3);    % CHANGE: first unwhiten template (apply inverse whitening matrix)
%     [~,KS_channel] = max(mean(abs(meanTemplate ),2));% CHANGE: max - min
    
    % consistent with phy: unwhiten template, then max - min
    cluster_templates = unwhitened_templates(:,:,origC);
    meanTemplate = mean(cluster_templates,3); 
    [~,KS_channel] = max(max(meanTemplate, [], 2) - min(meanTemplate, [], 2))
    
    %Find max waveform (not template)
    nChan = ops.NchanTOT;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ellie Changed This%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%size(cluster_templates, 1);
    wind = [-20 40];
    n_spike_samples = min(10000, length(clusterSpikeTimes));
    spike_samples_rand = randperm(length(clusterSpikeTimes), n_spike_samples);
    tmp = rez.st3(clusterSpikeIndices);
    spike_samples_load = tmp(spike_samples_rand);
    FileInf = dir([SiProbePath 'amplifier.dat']);
    nSampsInDat = (FileInf.bytes/nChan/2);
    valid_samples = (spike_samples_load > -wind(1)) & (spike_samples_load < nSampsInDat - wind(2));
    spike_samples_load = spike_samples_load(valid_samples);
    [meanWF, ~] = readWaveformsFromDat([SiProbePath 'amplifier.dat'], nChan, 1:nChan, spike_samples_load, wind, []);
    [~,WF_channel] = max(max(meanWF, [], 2) - min(meanWF, [], 2));
    shank = kcoords(KS_channel);
    addChan = find(kcoords==shank,1)-1;
    meanWF2 = meanWF(kcoords == shank,:);
    [~,WF_channel2] = max(max(meanWF2, [], 2) - min(meanWF2, [], 2));
    WF_channel2 = addChan + WF_channel2;

    clusters(i).spikeTimes = unique(clusterSpikeTimes);
    clusters(i).clusterID = clusterOfInterest;
    clusters(i).maxChannelTemplate = KS_channel - 1;
    clusters(i).maxChannelWF = WF_channel - 1;
    clusters(i).coordinates = [xcoords(KS_channel) ycoords(KS_channel)];
    clusters(i).shank= shank;
end