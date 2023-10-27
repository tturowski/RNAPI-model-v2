function rhoExp = doExperiment(rho, mids, rho0, sigma,name)
 
    dx = mids(2)-mids(1);
    cutoff = 0.5*( 1 + erf( (rho - rho0)/sigma ) );
    rhoExp = cutoff.*rho ;
    nParticlesExp = sum(rhoExp*dx);    
    rhoExp = rhoExp/nParticlesExp;
    
    name = [name '_rho0Exp' num2str(rho0) '_sigmaExp' num2str(sigma)];
    
    csvDataExp = [mids; rhoExp];
    csvFileExp = ['Data' filesep 'nParticles' num2str(nParticlesExp) join(name,'Exp.csv')];
    csvwrite(csvFileExp,csvDataExp);
end
