function sigma2Torque()
    
    stdDis = 100;
    mi = 300;
    c = (mi*pi^2*1^4)/10.5;
    range = [-20:20];
    
    Ti = c*( log(stdDis./(stdDis-range)) - log(stdDis./(stdDis+range)) );
    
    figure;
    plot(range,Ti,'-o');
    grid on
    grid minor
%     title(['DNA spring constant c for velocity = ' num2str(vInt)]);
    xlabel('DNA overwind \sigma [%]');
    ylabel('Torque \tau [pN nm]');
end