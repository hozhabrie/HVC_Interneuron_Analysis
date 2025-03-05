function plotEx(file,ts,ch,waveform)
tSec=ts/2e4;%in seconds
cla;hold on;
dur=.02;
n=nan(1,length(tSec));
for i=1:length(tSec)
    n(i)=sum(tSec>tSec(i)& tSec<(tSec(i)+dur));
end
loc=find(n>2,1,'first');
if isempty(loc)
    loc=1;
end
start=tSec(loc)-.005;%move 5 ms back
v=.195*(LoadBinary(file,'nChannels',120,'channels',ch,'start',start,'duration',dur));
t=(1:length(v))/2e4+start;%in seconds
plot(t,v);
%  test others, ignore
% ts=rez.st3(:,1);
% tSec=rez.st3(:,1)/2e4;
% tSec=tSec(tSec>t(1)&tSec<t(end));
% ts=ts-t(1)*2e4+1;
% ts(ts<0|ts>length(v))=[];
% plot(tSec,v(ts),'kx')
% ts=spk_t{c};
% tSec=ts/2e4;
tSec=tSec(tSec>t(1)&tSec<t(end));
ts=round(ts-t(1)*2e4+1);
ts(ts<0|ts>length(v))=[];
plot(tSec,v(ts),'rx')
ylabel('\muV')
xlabel('time (s)')
axis tight
%%