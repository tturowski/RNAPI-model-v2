function histData = RunModels()

% experimental post-processing
Exp.doExperiment = true;
Exp.rho0 = 3e-4;
Exp.sigma = 1e-4;

% opts.vInt = 50;   % intrinsic velocity
opts.sigma = 0.4; % intrinsic velocity variance parameter

% opts.v5end = 0.2; % friction slowing down RNAP without spring 0.2
% opts.dist = 2000; % distance where DNA works fully as a spring 2000
opts.lenValue = 7000; %default 7000
opts.outOfBound = opts.lenValue-250; % 6750 outer cutoff
opts.rpa12Cleavage = 1;
% opts.preTermProb = 0.01; % probability of premature/non-productive termination 0.01
% opts.preTermDist = 6750; % distance where preTermProb is applied 6750

opts.minRange = 38;   % size of RNAP
opts.RNAPbubble = 11; % size of RNAP bubble

opts.hybridFile = 'transcriptionBubble_dG.csv'; % hybrid data file

% directory for data storage
dataDir = 'Data';
if(~exist(dataDir,'dir'))
    mkdir(dataDir);
end
addpath(genpath(pwd));    % add Data and code directories to the path
saveData = true;          % whether to save data
opts.saveHistOnly = true; % only save the histograms to save storage space
recompute = true;        % DataStorage option to force recompute

opts.nRuns = 16; % number of runs - 50% out of 150 repeats PMID: 24190922
                % to calculate for 1 dubling time multiply 16 runs by 4.7
opts.tMax = 6000; % 6000 total run time PMID: 17713537
opts.dt = 0.004; % time step 0.008

optsExtra.doParforSimple = true; % parfor loop option

%--------------------------------------------------------------------------
%      models
%--------------------------------------------------------------------------
fileName = "runParameters_cd.csv";
names = getNames(fileName)
parameters = getParams(fileName);
% names = [...
%     "discrete",         'w65_dGadd15.csv'];
% % 	c,dist,RNA,v5end,   addProb,vInt,preTermProb,preTermDist
% parameters = [...
%     0,  0,  0,  0,      0.8,  50,  0,  6750];

count = 0; % run counter

for n= 1:size(names,1)
    count = count+1;
    elems = parameters(n,:); %parameters for a model
    
%     structureFile = names(n,2);
    structureFile = 'w65_dGadd15.csv';
    opts.structureFile = structureFile; % structure data file
    
%     opts.topProb = 0.2;

    rnt1Cleavage=3; %0 none, 1 no spring, 2 and 3 no spring no dG
    backtrackingVelocity=0;
%     addProb = 0.8; % [0.7,0.8,0.9]
    % parameters
    opts.c = elems(1); % [400,500,600]
    opts.dist = elems(2); % [2000]
    Strength = elems(3); % [1,1.25,1.5]
    opts.v5end = elems(4); % friction slowing down RNAP without spring 0.2
    opts.addProb = elems(5);
    opts.vInt = elems(6);
    opts.preTermProb = elems(7);
    opts.preTermDist = elems(8);

    structure2consider = -11; % [-10,-11,-12]
    ratio = 0.48; % [0.32,0.48,0.64]
%     ratio = 0; % [0.32,0.48,0.64]

    opts.backtrackingVelocity = backtrackingVelocity;
    opts.rnt1Cleavage = rnt1Cleavage;  
    opts.structure2consider = structure2consider;
    opts.structureStrength = -1 * Strength;
    opts.hybridStrength = Strength * ratio; 

    % csv file name
%     name = ['_addProb' num2str(opts.addProb)  '_c' num2str(opts.c) ...
%             '_structStrength' num2str(opts.structureStrength) '_structConsid' num2str(opts.structure2consider) ...
%             '_hybridStrength' num2str(opts.hybridStrength) '_rnt1Cl' num2str(opts.rnt1Cleavage) ...
%             '_backVel' num2str(opts.backtrackingVelocity) nameRoot];
    name = char(names(n,1))
    opts.name = name;

    % run simulations
    if (saveData == true)
        data = DataStorage(dataDir,@generateDataDiscrete,opts,optsExtra,recompute);
    % data.xEnd
    elseif (saveData == false)
        data = generateDataDiscrete(opts,optsExtra);
    end

%     % produce histograms
%     if(opts.saveHistOnly)
%         histDataTemp = data;
%     else
%         doAnalysePlots = false;
%         histDataTemp = analyseDataDiscrete(data,name,doAnalysePlots);
%     end
% 
%     if(Exp.doExperiment)
%         histDataTemp.rhoFullExp = doExperiment(histDataTemp.rhoFull, histDataTemp.mids, Exp.rho0, Exp.sigma,opts.name);
%     end
% 
%     histDataTemp.Exp = Exp;
%     histDataTemp.opts = opts;
% 
%     % store data in a structure
%     histData(count) = histDataTemp;
end

function outArray = getNames(fileName)
    importFile = detectImportOptions(fileName);
    tableFile = readtable(fileName, importFile);
    nameArray = table2array(tableFile(:,1));
    outArray = string(nameArray);
end

function outArray = getParams(fileName)
    importFile = detectImportOptions(fileName);
    tableFile = readtable(fileName, importFile);
    outArray = table2array(tableFile(:,2:9));
end

end