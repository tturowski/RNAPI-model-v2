function p = top1Prob(maxProb)
    sizeTop1 = 25+38;
    distMax = 250;
    x = [sizeTop1:distMax];
    
%     p = (distMax-x)/(distMax-sizeTop1)*maxProb; %lower distance gives higher prob to cut but usually is usual to topProb

%     p = (distMax-sizeTop1)./(distMax-x)*maxProb; %higher distance gives
%     higher prob (hyperbolic)

    p = ((x-distMax)./(distMax-sizeTop1)*maxProb)+maxProb; %higher distance gives higher prob (linearly)
    
    figure;
    plot(x,p);
end