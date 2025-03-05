function [audioMax,audioMin,tStop]=audioBounds(audio,ds,tStop)

numCuts=ceil(length(audio)/ds);
addlen=numCuts*ds-length(audio);
audio=[audio;zeros(addlen,1)];%zero pad the end
tStop=tStop+addlen/4e4;%increase this part
audioBox=reshape(audio',ds,numCuts);
audioMax=max(audioBox,[],1);
audioMin=min(audioBox,[],1);
