function amplitude_trace(t,k_amp,v_amp,th)

% plot(t,k_amp,'.');
% axis tight;
% y=ylim;
% ylim([0,y(2)]);
% line(xlim,[1,1]*th,'color','b')
% line(4.83*[1,1],ylim,'color','r')
plot(t,v_amp,'.');
axis tight;
% y=ylim;
% ylim([0,y(2)]);
line(xlim,[1,1]*0,'color','b')
line(4.83*[1,1],ylim,'color','r')
ylabel('\muV')
xlabel('time (m)')
% > set(get(AX(2),'Ylabel'),'string','direction') 



% [ax,h1,h2]=plotyy(t,k_amp,t,v_amp);
% h1.Marker='.';h2.Marker='.';
% h1.LineStyle='none';h2.LineStyle='none';
% axes(ax(1));
% axis tight;
% ylim(max(abs(k_amp))*[-1,1])
% line(ax(1),xlim,[1,1]*th,'color','b')
% line(ax(1),4.83*[1,1],ylim,'color','r')
% axes(ax(2));
% axis tight;
% line(xlim,[0,0])
% ylim(max(abs(v_amp))*[-1,1])
% set(gca,'ytick',-60:20:60);
% set(get(AX(2),'Ylabel'),'string','direction') 