function VelocityHistogram(l,vInt)
    
    sigma=0.4;    
    
    vRandom = rand(l,1);
    distPause = (vRandom <= 0.078);
    distInt = (vRandom > 0.078);

    nPause = nnz(distPause);
    nInt = nnz(distInt);
    
    vRandomInt = randn([nInt,1]) .* sigma.*vInt + vInt;
    vRandomPause = randn([nPause,1]) .* 1.5 + 0.9;

    
    vRandom(distPause) = vRandomPause;
    vRandom(distInt) = vRandomInt;
    
    figure;
    histogram(vRandom);
%     title(['Velocity Distribution for vInt = ' num2str(vInt)]);
    xlabel('Velocity of RNAP')
end