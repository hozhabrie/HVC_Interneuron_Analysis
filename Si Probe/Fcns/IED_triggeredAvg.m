function IED_triggeredAvg(ts,tIED)
l=-500;
r=400;
bw=20;
bins=l:bw:r;
[times,~]=timeLock(ts,tIED,[l,r]);%turn both into ms
n=histcounts(times,bins);
stairs(bins(1:end-1)+bw/2,n/bw*1e3/length(tIED))
line([-490,-200],[1,1],'color','k')
line([35,200],[1,1],'color','r')
xlabel('time (ms)')
ylabel('Hz')
axis tight
y=ylim;
ylim([0,y(2)])