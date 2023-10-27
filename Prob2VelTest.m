function Prob2VelTest(n)
    tMax = 1;
    dt = 0.001;
    nSteps = tMax/dt;
    
    k1 = 70;

    x = zeros(n,1); %positions

    for tStep = 1:nSteps
        prob = rand(size(x));
        
        moveMask = (prob < k1*dt);

        x(moveMask) = x(moveMask) + 1;
       
    end
    
    figure;
    histogram((x)./tMax);
end