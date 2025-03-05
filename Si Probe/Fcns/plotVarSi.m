function plotVarSi(allCSi,field,which,FRs)

si = getfield(allCSi,field);
si(isinf(si)) = nan;
edges = min(si):abs(range(si)/100):max(si);

figure('Visible','OFF');
set(gcf,'WindowState','maximized')

    switch nargin
        case 2
            edges = min(si):abs(range(si)/100):max(si);
            subplot(5,1,1)
            histogram(si,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title(field)
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,2)
            siTransform = si;
            zeroed = sum(si==0);
            siTransform(si==0) = nan;
            siTransform = log(siTransform);
            edgesT = min(siTransform):abs(range(siTransform)/100):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': log...num 0s nanned: ' + string(zeroed)])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,3)
            siTransform = nthroot(si,2);
            edgesT = min(siTransform):abs(range(siTransform)/100):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': sqrt'])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,4)
            siTransform = nthroot(si,3);
            edgesT = min(siTransform):abs(range(siTransform)/100):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': cube root'])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,5)
            bins = 25;
            siTransform = discretize(si,bins);
            edgesT = min(siTransform):abs(range(siTransform)/bins):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': discretized'])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
        case 4
            si = si(which);
            edges = min(si):abs(range(si)/100):max(si);

            subplot(5,1,1)
            histogram(si,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title(field)
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,2)
            siTransform = si;
            zeroed = sum(si==0);
            siTransform(si==0) = nan;
            siTransform = log(siTransform);
            edgesT = min(siTransform):abs(range(siTransform)/100):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': log...num 0s nanned: ' + string(zeroed)])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,3)
            siTransform = nthroot(si,2);
            edgesT = min(siTransform):abs(range(siTransform)/100):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': sqrt'])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,4)
            siTransform = nthroot(si,3);
            edgesT = min(siTransform):abs(range(siTransform)/100):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': cube root'])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
            subplot(5,1,5)
            bins = 25;
            siTransform = discretize(si,bins);
            edgesT = min(siTransform):abs(range(siTransform)/bins):max(siTransform);
            histogram(siTransform,edgesT,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title([field + ': discretized'])
            set(gca, 'TickDir', 'out','color','none','box','off')
    
        case 3
            si = si(which);
            edges = min(si):abs(range(si)/100):max(si);

            subplot(5,1,1)
            histogram(si,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
            title(field)
            set(gca, 'TickDir', 'out','color','none','box','off')
            return
    end
end

%   siZoom = rmoutliers(si,'ThresholdFactor',5);
%   edgesZoom = min(siZoom):abs(range(siZoom)/100):max(siZoom);

