%%

clear all; close all; clc;
path = 'C:/Users/User/Dropbox (NYU Langone Health)/Int Juxta/';
goodCells = readtable([path 'goodCellsnew.xlsx'])
birds = string(goodCells{:,1});
cells = string(goodCells{:,2});
toPrint = [(goodCells{:,5})];
identity = string(goodCells{:,9});


%%

counter = 0;
for i = 1:height(goodCells) %[1:6,12,14:15,17:19,21:25,27:45]%51,52:61,63,65:67,69]%height(goodCells)
    if toPrint(i) ==  1 & identity(i) ~= 'x'
        close all;
        counter = counter + 1;
        changedir = [path + birds(i) + '/' + cells(i) + '/'];
        im1 = imread([changedir + birds(i) + ' ' + cells(i) + ' wf.pdf']);
        im2 = imread([changedir + birds(i) + ' ' + cells(i) + ' pb trials.pdf']);
        im3 = imread([changedir + birds(i) + ' ' + cells(i) + ' pb stats.pdf']);
        im4 = imread([changedir + birds(i) + ' ' + cells(i) + ' isis.pdf']);
        im5 = imread([changedir + birds(i) + ' ' + cells(i) + ' spont.pdf']);
        % im3 = imresize(im3, [656, 875]);
        im6 = imread([changedir + birds(i) + ' ' + cells(i) + ' spont stats.pdf']);
      %  im7 = imread([changedir + birds(i) + ' ' + cells(i) + ' annotated ifr.pdf']);
        
        figure('WindowState','maximized');        
        montage({im1,im4,im2,im5,im3,im6},'Size', [3 2]);shg;
        title(string(counter),'FontSize',30)
        set(gca, 'box', 'off', 'color', 'none','xcolor','none','ycolor','none');
        set(gcf,'PaperOrientation','Landscape', 'PaperUnits', 'inches');
        x = xlim;
        y = ylim;
        text(range(x)*.75,range(y)*.95,identity(i),'FontSize',30)
        %  ha=get(gcf,'children');
        %  ha(2, 1).Position  = [0.13    0    0.5    0.4];
        mkdir([path 'donefigs/'])
        cd([path 'donefigs/'])
        print([string(counter) + ' ' + birds(i) + ' ' + cells(i) + ' comm'], '-dpdf','-fillpage');
        pause(.5)
        %     close all;
    end
end


