function m=auto_corr(ts,width,binSize)
tSpike=ts/2e4;%convert from samples to seconds
if nargin<2
    width=30/1e3;%in s
    binSize=.3/1e3;%in s
end
newFs=1/binSize;
% resample at a lower Fs in tune with bins
Sinds=ceil(tSpike*newFs);
Sinds=Sinds-Sinds(1)+1;
tSpike2=zeros(1,max(Sinds));%this is more memory intensive but faster than a for loop
tSpike2(Sinds)=1;
[c,lags]=xcorr(tSpike2,newFs*width);
c(ceil(end/2))=NaN;
c=c/max(max(c),.05);
c(isnan(c))=0;
[i,d]=max(c);
m=lags(d)*binSize*1e3;
stairs(lags*binSize*1e3,c);
xlabel('time (ms)')
ylabel('corr')
title('AutoCorrelogram')
axis tight;
yy=ylim;
ylim([0,yy(2)]);