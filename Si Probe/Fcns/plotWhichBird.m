function plotWhichBird(allCSi,which,field)

    si = getfield(allCSi,field);
    si = si(which);
    %si = sqrt(si);
    birds = allCSi.bird(which);
    si(isinf(si)) = nan;
    edges = min(si):abs(range(si)/60):max(si);
    
    figure('Visible','OFF');
    set(gcf,'WindowState','maximized')
    ax(1) = subplot(8,1,[1:5]);
    histogram(si,edges,'DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
    title(field)
    set(gca, 'TickDir', 'out','color','none','box','off','FontSize',18)
    hold on
    high = ylim;
    yticks(0:floor(high(2)/3):high(2))
    label = {'eh22' 'eh24' 'eh50' 'eh52' 'eh57'};
    color = {'b' 'r' 'm' 'g' 'k'};
    y = 1; %.1*max(ylim);
    tickLim = [2];

    ax(2) =  subplot(8,1,[6:8]);
    for a  = 1:length(label)
       whichBird = si(birds==label(a));
       for i = 1:length(whichBird)
            %line([whichBird(i) whichBird(i)],tickLim-3*y*a,'Color', color{a})
            plot(whichBird(i),-3*y*a,'|','Color',color{a},'MarkerSize',15)
            hold on
            count = sum(whichBird == whichBird(i));
            if count > 1
                text(whichBird(i),-.5-3*y*a,string(count),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',15)
            end
            hold on
        end
    end
    ylim([-16 0])
    linkaxes(ax,'x')
    set(gca, 'TickDir', 'out','color','none','box','off','xcolor','none','ycolor','none')
    set(gcf,'PaperOrientation','Landscape')
end



