clear all; close all; clc;

path = 'C:/SiProbe/C52/data/';
bird = 'C52';
load([path 'clusters.mat']);
load('128_5_integrated.mat')


figure; hold on;
p=plot(xcoords,ycoords,'o','MarkerSize',20);

set(gca,'box','off','color','none','xcolor','none','ycolor','none')

for i=1:128
    text(xcoords(i),ycoords(i),num2str(chanMap0ind(i)),'HorizontalAlignment','center','Color', 'b')
    %text(xcoords(i)-10,ycoords(i),num2str(i),'HorizontalAlignment','right')
end


title(['Chan Map for ' bird]);
maxChannels = [clusters.maxChannelWF]+1;
ids = string([clusters.clusterID]);
for i = 1:128
   whichUnits = find(maxChannels == i);
   text(xcoords(i)+5, ycoords(i), ids(whichUnits))
end

set(gcf,'PaperOrientation','landscape', 'PaperUnits', 'inches', 'PaperSize', [20 10]);

cd(path)
print([bird ' chan map'], '-dpdf','-fillpage');  
