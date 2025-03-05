function plot_avgs(waveforms,pos,ch,t)
cla;hold on;
if nargin<4
    w=(size(waveforms,2)-1)/2;
    t=(-w:w)/20;
end
xO=(length(t)/20)*1.2;%offset in x
yO=20;%offset in y
% waveforms=waveforms*.195;
for wi=1:size(waveforms,1);
    w=waveforms(wi,:);
    if wi==ch
        plot(t+xO*pos(wi,1),w +yO*pos(wi,2),'r');
        plot(xO*pos(wi,1),w(t==0) +yO*pos(wi,2),'rx');
    else
        plot(t+xO*pos(wi,1),w +yO*pos(wi,2),'k');
    end
end
axis tight;
set(gca,'xtick',[])
set(gca,'ytick',[])
ylim([-yO,(max(pos(:,2))+1)*yO])
if xO<0
    xlim([(max(pos(:,1))+1)*xO,-xO])
else
    xlim(fliplr([(max(pos(:,1))+1)*xO,-xO]))
end
xLoc=12*sign(xO);xW=1;
yLoc=-yO+5;yH=10;
line([0,xW]+xLoc,yLoc+[0,0],'color','k')
line([0,0]+xLoc,yLoc+[0,yH],'color','k')