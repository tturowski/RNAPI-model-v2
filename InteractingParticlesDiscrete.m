function data = InteractingParticlesDiscrete(opts)

%--------------------------------------------------------------------------
% set up parameters
%--------------------------------------------------------------------------

vInt = opts.vInt; % intrinsic velocity
sigma = opts.sigma; % intrinsic velocity standard deviation parameter

addProb = opts.addProb; % probability to add a particle

c = opts.c;    % winding strength
d = opts.dist; % distance from when DNA works fully as a spring

lenValue = opts.lenValue;
rnt1CleavageSite = lenValue-278; %6697+25 (=7000-278) to make it exact as RNAPI CRAC peak
v5end = opts.v5end; %friction slowing down RNAP without spring

% hybrid parameters
hybridStrength = opts.hybridStrength;
temp = csvread(opts.hybridFile);
hybrid = temp(:,2);

% structure parameters
structureStrength = opts.structureStrength;
temp = csvread(opts.structureFile);
structure = temp(:,2);
structure(structure > opts.structure2consider) = 0;

if (lenValue==14000)
    hybrid = [hybrid;hybrid];
    structure = [structure;structure];
end

if (opts.rnt1Cleavage == 2)
    structure(rnt1CleavageSite:rnt1CleavageSite+50) = 0;
elseif (opts.rnt1Cleavage >= 3)
    structure(rnt1CleavageSite:lenValue) = 0;
end



% size of RNAP
minRange = opts.minRange;
RNAPbubble = opts.RNAPbubble;

% top1 parameters
if(isfield(opts,'topProb'))
    sizeTop1 = minRange + 25; % RNAPI footprint + Top1 size
    topProb = opts.topProb;
else
    topProb = 0;
end

% time stepping
tMax = opts.tMax;
dt = opts.dt;
nSteps = tMax/dt;

%premature termination
estTime = opts.outOfBound / opts.vInt; %time = distance / velocity
preTermProb = (opts.preTermProb / estTime) * dt;
preTermDist = min([opts.preTermDist, rnt1CleavageSite]);

% acceptance probability for moving
acceptProb = dt;

% length of domain
outOfBound = opts.outOfBound;

% whether to do backtracking
backtrackingVelocity = opts.backtrackingVelocity;

% start with one RNAP
nRNAP = 1;   % number of RNAP
x = 0;       % position
tIn = 0;     % entry time
tTaken = []; % time taken to leave
backtrackingState = 0; % backtracking state
xInit = inf; % distance to the next particle
xEnd = 0; % particles that finished transcription productively 
xEndPremature = 0; % particles that finished transcription non-productively

% probability of adding an RNAP in each time step
addProb = addProb*dt;

% only save at certain steps as the files become too large otherwise
nSaves = 1000;
iSave = 0;
saveFreq = ceil(nSteps/nSaves);

% preallocate for speed (doesn't set a maximum number of particles)
maxNParticles = 80;
xFull = nan(maxNParticles,nSaves);
dxAvgFull = nan(maxNParticles,nSaves);
structureFull = nan(maxNParticles,nSaves);
hybridFull = nan(maxNParticles,nSaves);
windingFull = nan(maxNParticles,nSaves);
backtrackingFull = nan(maxNParticles,nSaves);
collisionsFull = nan(maxNParticles,nSaves);
tFull = zeros(nSaves,1);

for tStep = 1:nSteps

    %----------------------------------------------------------------------
    % Addition of new RNAP
    %----------------------------------------------------------------------
 
    % random chance of adding RNAP
    toAdd = rand(1);
    
    if(toAdd < addProb) % attempt to add
        if(x(1)>minRange) % don't add if existing RNAP is too close
            xInit = [x(1);xInit];  % create new RNAP and set defaults
            x = [0;x]; 
            backtrackingState = [0;backtrackingState];
            nRNAP = nRNAP+1;
            tIn = [tStep*dt;tIn];
        end
    end
    
    %----------------------------------------------------------------------
    % Winding forces
    %----------------------------------------------------------------------
    
    % distances to RNAP in front and behind
    dxP = abs([x(2:end) - x(1:end-1) ;inf]);
    dxM = abs([inf;x(2:end) - x(1:end-1)]);
    
    % compute forces from neighbour RNAP
    vP = c*(1-(xInit-RNAPbubble)./(dxP-RNAPbubble));
    vM = c*(1-([inf;xInit(1:end-1)]-RNAPbubble)./(dxM-RNAPbubble));

    % remove  DNA spring for the very first and very last particle
    vM(1) = vP(1);
    vP(end) = vM(end);

    %----------------------------------------------------------------------
    % reduce DNA winding forces for the first 1-2 kb 
    %----------------------------------------------------------------------

    % mask for first d nt (typically 2000)
    startMask = (x < d);
    startProgress = x(startMask) ./ d;
      
    % reduce forces
    vP(startMask) = vP(startMask) .* startProgress;
    vM(startMask) = vM(startMask) .* startProgress;
    
    %reseting xInit in the first 1-2 kb - probability function
    xInitShift = xInit(startMask) - dxP(startMask);
    xInitNew = xInit(startMask) - ceil(xInitShift .* (1-startProgress)); %level of xInit reseting

    %%prob function of xInit reseting
    resetProb = 1 * dt;
    pTestxInit = rand(size(xInitNew));
    xInitResetMask = (pTestxInit < resetProb);
    xInitUpdate = xInit(startMask);
    xInitUpdate(xInitResetMask) = xInitNew(xInitResetMask);
    xInit(startMask) = xInitUpdate;

    %----------------------------------------------------------------------
    % Rnt1 cleavage
    %----------------------------------------------------------------------

    % remove DNA spring after Rnt1 cleavage
    if (opts.rnt1Cleavage >= 1) %remove spring mechanism after Rnt1 cleavage
        rnt1Mask = (x > rnt1CleavageSite);
        vM(rnt1Mask) = 0;
        vP(rnt1Mask) = 0;
        
        if any(rnt1Mask > 0 & opts.rnt1Cleavage >= 4)
            rnt1Mask2 = (~rnt1Mask);
            rnt1Mask2(end) = 1; %removing spring after the last particle
            vP(rnt1Mask2) = 0;
        end
    end
    
    %----------------------------------------------------------------------
    % draw random velocities
    %----------------------------------------------------------------------

    % new random velocity
    vRandom = randomVelocity(nRNAP,startMask,startProgress);
    
    %assign new velocity for bactracked state
    if (backtrackingVelocity) 
        bactrackedMask = (backtrackingState < 0);
        vRandomPause = randn([nnz(bactrackedMask),1]) .* 1.5 + 0.9;
        vRandom(bactrackedMask) = vRandomPause;
    end
    
    %----------------------------------------------------------------------
    % determine total velocities
    %----------------------------------------------------------------------
        
    % modification of random velocity
    % dG of structure and GC RNA:DNA hybrid
    vTilde = vRandom + structureStrength*structure(x+1) - hybridStrength*hybrid(x+1);
    % DNA winding forces
    v = vTilde + vP - vM;
    
    % fix NaN issues with only one particle
    if(nRNAP==1)
        v = vInt;
    end
    
    %----------------------------------------------------------------------
    % restrict overlapping with other RNAP
    %----------------------------------------------------------------------

    closeP = (dxP <= minRange);
    closeM = (dxM <= minRange);
    closeM(1) = (x(1)<= minRange);
    closeP(end) = false;
    collisionsState = closeP+closeM; %save information about RNAPI collisions

    v(closeP) = min(v(closeP),0); %do not move fwd if RNAP in front
    v(closeM) = max(v(closeM),0); %do not move bwd if RNAP behind
    
    closeM(1) = 0; %for collisionsState save
    
    %----------------------------------------------------------------------
    % determine if the RNAP moves or not, and which direction
    %----------------------------------------------------------------------

    p = rand(size(x))./ v;
    moveMask = (abs(p) <= acceptProb);
        
    stepMove = zeros(size(p));
    stepMove(moveMask) = sign(p(moveMask)); %recording move direction

    %----------------------------------------------------------------------
    % Rpa12 cleavage
    %----------------------------------------------------------------------

    backtrackingState = min(backtrackingState + stepMove, 0); %recording backtracking
    
    if (opts.rpa12Cleavage == 1)
        rpa12Mask = (backtrackingState < -5);
        deepBack = backtrackingState(rpa12Mask);
        deepBack(rpa12Prob(deepBack)) = -1;
        backtrackingState(rpa12Mask) = deepBack;
    end
    
    %----------------------------------------------------------------------
    % do the moves
    %----------------------------------------------------------------------
    
    x(moveMask) = x(moveMask) + sign(p(moveMask));

    %----------------------------------------------------------------------
    % kill off RPNA that run off the end
    %----------------------------------------------------------------------

    xOutOfBound = (x>outOfBound);
    if(nnz(xOutOfBound)>0)
        x(xOutOfBound) = [];
        xInit(xOutOfBound) = []; 
        tTaken = [tTaken; tStep*dt - tIn(xOutOfBound)];
        tIn(xOutOfBound) = [];
        nRNAP = length(x);
        dxP(xOutOfBound) = [];
        dxM(xOutOfBound) = [];
        backtrackingState(xOutOfBound) = [];
        collisionsState(xOutOfBound) = [];
        xEnd = xEnd+nnz(xOutOfBound);
    end

    %----------------------------------------------------------------------
    % apply premature termination
    %----------------------------------------------------------------------

    if(opts.preTermProb > 0)
        xPreTerm = false(size(x));
        preTermMask = (x<preTermDist);
        pTest = rand(size(preTermMask(preTermMask==1)));
        xPreTerm(preTermMask) = (pTest<preTermProb);
        if ((nnz(xPreTerm)>0) && nRNAP>1)
            x(xPreTerm) = [];
            xInit(xPreTerm) = [];
            tIn(xPreTerm) = [];
            nRNAP = length(x);
            dxP(xPreTerm) = [];
            dxM(xPreTerm) = [];
            backtrackingState(xPreTerm) = [];
            collisionsState(xPreTerm) = [];
            xEndPremature = xEndPremature+nnz(xPreTerm);
        end
    end

    %----------------------------------------------------------------------
    % apply Top1
    %----------------------------------------------------------------------    
    
    if(topProb>0)
        top1Mask = (dxP>=sizeTop1); %size mask

        top1MaskP = top1Mask & (dxP-xInit>=10); %overwinding
        top1MaskM = top1Mask & (dxP-xInit<=-10); %underwinding

        top1DistP = dxP(top1MaskP);
        top1DistM = dxP(top1MaskM);
        probTop1P = zeros(size(dxP));
        probTop1M = zeros(size(dxP));
        probTop1P(top1MaskP) = top1Prob(top1DistP);
        probTop1M(top1MaskM) = top1Prob(top1DistM);

        pTestP = rand(size(top1DistP));
        top1CutP = false(size(dxP));
        top1CutP(top1MaskP) = (pTestP<probTop1P(top1MaskP));
        pTestM = rand(size(top1DistM));
        top1CutM = false(size(dxP));
        top1CutM(top1MaskM) = (pTestM<probTop1M(top1MaskM));

        xInit(top1CutP) = xInit(top1CutP) + 10;
        xInit(top1CutM) = xInit(top1CutM) - 10;
    end
        
    %----------------------------------------------------------------------
    % save data
    %----------------------------------------------------------------------
    
    if(mod(tStep,saveFreq)==1)   
        iSave = iSave + 1;
        
        % position
        xFull(1:nRNAP,iSave) = x(1:nRNAP);
        xFull(nRNAP+1:end,iSave) = NaN; % pad xFull with NaNs where there are no particles
        % time
        tFull(iSave) = tStep*dt;

        
        % average distance of particle to its neighbours
        dxAvg = (dxM+dxP)/2;
        dxAvg(1) = dxP(1);
        dxAvg(end) = dxM(end);
        dxAvgFull(1:nRNAP,iSave) = dxAvg;

        % forces
        structureFull(1:nRNAP,iSave) = structureStrength*structure(x+1);
        hybridFull(1:nRNAP,iSave) = - hybridStrength*hybrid(x+1);
        winding = vP - vM;
        windingFull(1:nRNAP,iSave) = winding(1:nRNAP); % deal with ones that have left the domain
        
        % backtracking
        backtrackingFull(1:nRNAP,iSave) = backtrackingState;
        
        %collisions        
        collisionsFull(1:nRNAP,iSave) =  collisionsState;
    end
    
    
end

%--------------------------------------------------------------------------
% output data
%--------------------------------------------------------------------------

data.tMax = tMax;
data.dt = dt;
data.xFull = xFull;
data.tFull = tFull;
data.dxAvgFull = dxAvgFull;
data.structureFull = structureFull;
data.hybridFull = hybridFull;
data.windingFull = windingFull;
data.backtrackingFull = backtrackingFull;
data.collisionsFull = collisionsFull;
data.xEnd = xEnd;
data.xEndPremature = xEndPremature;
data.tTaken = tTaken;

%--------------------------------------------------------------------------
% auxiliary functions
%--------------------------------------------------------------------------

    function vRandom = randomVelocity(N,startMask,startProgress)
        
        vRandom = rand(N,1);
        distPause = (vRandom <= 0.078);
        distInt = (~startMask);
        distInt5end = (startMask);
        
        nPause = nnz(distPause);
        nInt = nnz(distInt);
        nInt5end = nnz(distInt5end);
        
        vInt5end = vInt .* (1-v5end + v5end.*startProgress);
        vRandomInt5end = randn([nInt5end,1]) .* sigma.*vInt5end + vInt5end;
        vRandomInt = randn([nInt,1]) .* sigma.*vInt + vInt;
        vRandomPause = randn([nPause,1]) .* 1.5 + 0.9;
        
        vRandom(distInt) = vRandomInt;
        vRandom(distInt5end) = vRandomInt5end;
        vRandom(distPause) = vRandomPause;

    end

    function p = rpa12Prob(x)
       maxProb = 0.2*dt;
       p = (rand(size(x)) < maxProb);
    end

    function p = top1Prob(x)
        maxProb = topProb*dt;
        distMax = 250;
        p = ((x-distMax)./(distMax-sizeTop1)*maxProb)+maxProb;
    end

end