function dataOut = analyseDataDiscrete(data,name,lenValue,doAnalysePlots)

    set(0,'DefaultLineLineWidth',2); set(0,'DefaultAxesFontSize',20);

    tFull = data(1).tFull;
    nTimes = length(tFull);
    nRuns = length(data);

    edges = 0:1:lenValue;
    mids = edges + edges(1)/2;
    mids = mids(1:end-1);

    dx = edges(2);

    tPos = 1;
    
    rhoFull = zeros(size(mids));
    dxAvgFull = zeros(size(mids));
%     dxAvg70Full = zeros(size(mids));
    structureFull = zeros(size(mids));
    hybridFull = zeros(size(mids));
    windingFull = zeros(size(mids));
    backtrackingFull = zeros(size(mids));
    collisionsFull = zeros(size(mids));
    xEndFull = [];
    xEndPrematureFull = [];
    tTakenFull = [];
    nBins = length(mids);

    rhoFullOut = [];
    dxAvgFullOut = [];
    structureFullOut = [];
    hybridFullOut = [];
    windingFullOut = [];
    backtrackingFullOut = [];
    backtrackingFullRawOut = [];
    collisionsFullOut= [];
    collisionsFullRawOut= [];
    nParticlesFull = [];


    for iRun = 1:nRuns

        xFull = data(iRun).xFull;
        dxAvg = data(iRun).dxAvgFull;
        structure = data(iRun).structureFull;
        hybrid = data(iRun).hybridFull;
        winding = data(iRun).windingFull;
        backtracking = data(iRun).backtrackingFull;    
        collisions = data(iRun).collisionsFull;  
        xEnd = data(iRun).xEnd;
        xEndPremature = data(iRun).xEndPremature;
        tTaken = data(iRun).tTaken;

        for iTime = tPos:nTimes
            
            % positional data
            x = xFull(:,iTime);
            killMask = (x>lenValue) | isnan(x);
            x(killMask) = [];
            xBoxes = floor( nBins* x / lenValue ) + 1;
            xBoxes(xBoxes==nBins+1)=nBins;
                        
            % number of particles
            rho = ( accumarray(xBoxes, ones(size(x))) ).'; %Matlab trick to make a histogram
            rhoFull(1:length(rho)) = rhoFull(1:length(rho)) + rho; % full histogram
            
            % average separation
            dxAvgT = dxAvg(:,iTime);
            dxAvgT(killMask) = [];
            dxAvgHist = accumarray(xBoxes, dxAvgT).'; 
            dxAvgFull(1:length(dxAvgHist)) = dxAvgFull(1:length(dxAvgHist)) + dxAvgHist;
            
%             % count the number with average separation less than 70
%             dxAvg70Hist = accumarray(xBoxes, (dxAvgT<70)).'; 
%             dxAvg70Full(1:length(dxAvgHist)) = dxAvg70Full(1:length(dxAvgHist)) + dxAvg70Hist;
            
            % forces
            structureT = structure(:,iTime);
            structureT(killMask) = [];
            structureHist = accumarray(xBoxes, structureT).'; 
            structureFull(1:length(structureHist)) = structureFull(1:length(structureHist)) + structureHist;

            hybridT = hybrid(:,iTime);
            hybridT(killMask) = [];
            hybridHist = accumarray(xBoxes, hybridT).'; 
            hybridFull(1:length(hybridHist)) = hybridFull(1:length(hybridHist)) + hybridHist;

            windingT = winding(:,iTime);
            windingT(killMask) = [];
            windingHist = accumarray(xBoxes, windingT).'; 
            windingFull(1:length(windingHist)) = windingFull(1:length(windingHist)) + windingHist;

            % backtracking
            backtrackingT = backtracking(:,iTime);
            backtrackingT(killMask) = [];
            backtrackingHist = accumarray(xBoxes, backtrackingT).'; 
            backtrackingFull(1:length(backtrackingHist)) = backtrackingFull(1:length(backtrackingHist)) + backtrackingHist;

            % collisions
            collisionsT = collisions(:,iTime);
            collisionsT(killMask) = [];
            collisionsHist = accumarray(xBoxes, collisionsT).'; 
            collisionsFull(1:length(collisionsHist)) = collisionsFull(1:length(collisionsHist)) + collisionsHist;

            
        end

        xEndFull = [xEndFull; xEnd];
        xEndPrematureFull = [xEndPrematureFull; xEndPremature];
        tTakenFull = concat(tTakenFull, tTaken);
%     end

        % average separation
        dxAvgFull = dxAvgFull./rhoFull;
        dxAvgFullOut = [dxAvgFullOut; dxAvgFull];
        
    %     % mean average separation over whole domain
    %     dxAvgCut = dxAvgFull(~isnan(dxAvgFull)); % strip off end
    %     dxAvgCut(isinf(dxAvgCut)) = []; % remove anywhere where the density is zero
    %     dxAvgMean = mean(dxAvgCut);
    
    %     % probability within 70
    %     dxAvg70Full = dxAvg70Full./rhoFull;
        
        % forces
        structureFull = structureFull./rhoFull;
        hybridFull = hybridFull./rhoFull;
        windingFull = windingFull./rhoFull;
        structureFullOut = [structureFullOut; structureFull];
        hybridFullOut = [hybridFullOut; hybridFull];
        windingFullOut = [windingFullOut; windingFull];

        % backtracking
        backtrackingFullRaw = backtrackingFull;
        backtrackingFull = backtrackingFull./rhoFull;
        backtrackingFullRawOut = [backtrackingFullRawOut; backtrackingFullRaw];
        backtrackingFullOut = [backtrackingFullOut; backtrackingFull];
        
        % collisions
        collisionsFullRaw = collisionsFull;
        collisionsFull = collisionsFull./rhoFull;
        collisionsFullRawOut = [collisionsFullRawOut; collisionsFullRaw];
        collisionsFullOut = [collisionsFullOut; collisionsFull];
        
        % normalise density
        rhoFull = rhoFull/(nTimes-tPos+1)/dx;
%         rhoFull = rhoFull/(nTimes-tPos+1)/dx/nRuns;
        nParticles = sum(rhoFull*dx);
        nParticlesFull = [nParticlesFull; nParticles];
        rhoFull = rhoFull/nParticles;
        rhoFullOut = [rhoFullOut; rhoFull];
    end
   
%     % output csv
%     csvData = [mids; rhoFull; dxAvgFull; structureFull; hybridFull; windingFull; backtrackingFull; backtrackingFullRaw; collisionsFull; collisionsFullRaw];
% %     csvFile = ['Data' filesep 'nParticles' num2str(nParticles) name '.csv'];
%     csvFile = ['Data' filesep 'nParticles' num2str(nParticles) '_' name '.csv'];
%     csvwrite(csvFile,csvData);
    
    % output data
    dataOut.mids = mids;
    dataOut.rhoFull = rhoFullOut; 
    dataOut.dxAvgFull = dxAvgFullOut; 
% %     dataOut.dxAvg70Full = dxAvg70Full; 
    dataOut.structureFull = structureFullOut; 
    dataOut.hybridFull = hybridFullOut; 
    dataOut.windingFull = windingFullOut; 
    dataOut.backtrackingFull = backtrackingFullOut;
    dataOut.backtrackingFullRaw = backtrackingFullRawOut;
    dataOut.collisionsFull = collisionsFullOut;
    dataOut.collisionsFullRaw = collisionsFullRawOut;
    dataOut.nParticles = nParticlesFull;
    dataOut.xEnd = xEndFull;
    dataOut.xEndPremature = xEndPrematureFull;
    dataOut.tTaken = tTakenFull;
    
    function z = concat(x,y)
        sx = size(x);
        sy = size(y);
        a = max(sx(1),sy(1));
        z = [[x;nan(abs([a 0]-sx))],[y;nan(abs([a,0]-sy))]];
    end

end