function [times,trial]=timeLock(spikes,events,range)
%times=spikes relative to 
times=[];
trial=[];
for e=1:length(events)
    tN=spikes-events(e);
    rm1=tN<range(1);
    rm2=tN>range(2);
    tN(rm1|rm2)=[];
    times=[times;tN];
    trial=[trial,repmat(e,1,length(tN))];
end
[times,inds]=sort(times);
trial=trial(inds);