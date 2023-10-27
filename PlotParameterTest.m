function allData = PlotParameterTest(histData)

    nLines = length(histData);
    
    figure('Position',[100,100,1600,800]);
    
    allData = zeros(nLines,7000);
    
    hl = zeros(nLines+1,1);
    
    for iLine = 1:nLines
        hl(iLine) = plot(histData(iLine).mids,histData(iLine).rhoFull,'color',[0,0,0,0.1]);
        hold on
        if( histData(iLine).opts.addProb == 0.8 ...
                && histData(iLine).opts.c == 500 ...
                && histData(iLine).opts.structureStrength == -1.25 ...
                && histData(iLine).opts.hybridStrength == 1.25*0.48 ...
                && histData(iLine).opts.structure2consider == -11)
            iRed = iLine;
        end
            
        allData(iLine,:) = histData(iLine).rhoFull;
    end

    hl(nLines+1) = plot(histData(iRed).mids,histData(iRed).rhoFull,'color',[1,0,0]);
       
    set(gca,'fontsize', 40)
    xlabel('Position')
    ylabel('Fraction of occupancy time')
    ylim([0.9*10^(-4),4*10^(-4)])
    
    ha1 = gca;
   
    p = get(ha1, 'Position');
    ha2 = axes('Parent', gcf, 'Position', [p(1)+.15 p(2)+.48 p(3)-.52 p(4)-.54]);
    hr2 = copyobj(hl,ha2);
    axes(ha2)
    hold on
    plot(histData(iRed).mids,histData(iRed).rhoFull,'color',[1,0,0]);
    xlim([400,500])
    ylim([1*10^(-4),3.5*10^(-4)])
    set(ha2,'fontsize',30, 'box','on');

    ha3 = axes('Parent', gcf, 'Position', [p(1)+.45 p(2)+.48 p(3)-.52 p(4)-.54]);
    hr3 = copyobj(hl,ha3);
    axes(ha3)
    hold on
    plot(histData(iRed).mids,histData(iRed).rhoFull,'color',[1,0,0]);
    xlim([1150,1250])
    ylim([1*10^(-4),2.5*10^(-4)])
    set(ha3,'fontsize',30, 'box','on');

    
    
end