function vRandom = randomVelocity(N)
    
    sigma = 0.4;
    vInt = 50;
    
    vRandom = rand(N,1);
    distPause = (vRandom <= 0.078);
    distInt = (vRandom > 0.078);

    nPause = nnz(distPause);
    nInt = nnz(distInt);

    vRandomInt = randn([nInt,1]) .* sigma.*vInt + vInt;
    vRandomPause = randn([nPause,1]) .* 1.5 + 0.9;

    vRandom(distPause) = vRandomPause;
    vRandom(distInt) = vRandomInt;

end