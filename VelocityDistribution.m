function VelocityDistribution(max,vInt)
    
    sigma = 4.9/12.8; %from Adelman et al PNAS 2002
    range = [-5:0.1:max];
    normInt = normpdf(range,vInt,sigma*vInt);
    normPause = normpdf(range,0.9,1.5);
    dist = [0.922*normInt+0.078*normPause];
    cumsumdist = cumsum(dist)/10;
    
    figure;
    plot(range,dist);
    title(['Velocity Distribution for vInt = ' num2str(vInt)]);
    xlabel('Velocity of RNAP');
%     cumsumdist
end