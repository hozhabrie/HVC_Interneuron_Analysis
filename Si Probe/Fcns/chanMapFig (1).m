load('128_5_integrated.mat')

figure;clf;hold on;
p=plot(xcoords,ycoords,'o','MarkerSize',20);
xlabel('\mum')
ylabel('\mum')
set(gca,'TickDir','out')
for i=1:128
    text(xcoords(i),ycoords(i),num2str(chanMap0ind(i)),'HorizontalAlignment','center','Color', 'b')
    %text(xcoords(i)-10,ycoords(i),num2str(i),'HorizontalAlignment','right')

end
xlim(xlim+[-1,1]*30)
ylim([-20 400]); 
title('Chan Map for C22');
maxChannels = [clusters.maxChannelWF]+1;
ids = string([clusters.clusterID]);
for i = 1:128
   whichUnits = find(maxChannels == i);
   text(xcoords(i)+5, ycoords(i), ids(whichUnits))
end