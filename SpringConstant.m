function SpringConstant(vInt)
    range = [1:20];
    c = -vInt ./ ( (100./(100+range)) - (100./(100-range)) );
    figure;
    plot(range,c,'-o');
    title(['DNA spring constant c for velocity = ' num2str(vInt)]);
    xlabel('DNA overwind \sigma (%)');
end