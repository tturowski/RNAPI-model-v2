function histData = ParametersTest()

% setup
figure('Position',[100,100,1200,800]);
ha = axes;

% experimental post-processing
Exp.doExperiment = true;
Exp.rho0 = 3e-4;
Exp.sigma = 1e-4;

if(Exp.doExperiment)
    figure('Position',[100,100,1200,800]);
    haExp = axes;
end

opts.vInt = 50;   % intrinsic velocity
opts.sigma = 0.4; % intrinsic velocity variance parameter

opts.v5end = 0.2; % friction slowing down RNAP without spring 0.2
opts.dist = 2000; % distance where DNA works fully as a spring 2000
opts.outOfBound = 6750; % outer cutoff
opts.rpa12Cleavage = 1;

opts.minRange = 38;   % size of RNAP
opts.RNAPbubble = 11; % size of RNAP bubble

name = 'w65_dGadd15.csv';
opts.structureFile = name; % structure data file
opts.hybridFile = 'transcriptionBubble_dG.csv'; % hybrid data file

% directory for data storage
dataDir = 'Data';
if(~exist(dataDir,'dir'))
    mkdir(dataDir);
end
addpath(genpath(pwd));    % add Data and code directories to the path
saveData = true;          % whether to save data
opts.saveHistOnly = true; % only save the histograms to save storage space
recompute = false;        % DataStorage option to force recompute

opts.nRuns = 256; % number of runs
opts.tMax = 2000; % total run time
opts.dt = 0.008; % time step

optsExtra.doParforSimple = true; % parfor loop option

nameRoot = ['_dist' num2str(opts.dist) '_n' num2str(opts.nRuns) ...
            '_v5end' num2str(opts.v5end)]; 
            
count = 0; % run counter

%opts.topProb = 0.2;

for(rnt1Cleavage=[3]) %0 none, 1 no spring, 2 and 3 no spring no dG
    for(backtrackingVelocity=[0])
        for(addProb = [0.7,0.8,0.9] ) 
            for(c = [400,500,600]) 
                for(Strength = [1,1.25,1.5]) 
                    for(structure2consider = [-10,-11,-12]) 
                        for(ratio = [0.32,0.48,0.64]) 
                            
                            count = count+1

                            opts.backtrackingVelocity = backtrackingVelocity;
                            opts.rnt1Cleavage = rnt1Cleavage;

                            opts.addProb = addProb;
                            opts.c = c;

                            opts.structure2consider = structure2consider;
                            opts.structureStrength = -1 * Strength;
                            opts.hybridStrength = Strength * ratio; 

                            % csv file name
                            name = ['_addProb' num2str(opts.addProb)  '_c' num2str(opts.c) ...
                                    '_structStrength' num2str(opts.structureStrength) '_structConsid' num2str(opts.structure2consider) ...
                                    '_hybridStrength' num2str(opts.hybridStrength) '_rnt1Cl' num2str(opts.rnt1Cleavage) ...
                                    '_backVel' num2str(opts.backtrackingVelocity) nameRoot];
                            opts.name = name;

                            % run simulations
                            if (saveData == true)
                                data = DataStorage(dataDir,@generateDataDiscrete,opts,optsExtra,recompute);
%                                 data.xEnd
                            elseif (saveData == false)
                                data = generateDataDiscrete(opts,optsExtra);
                            end

                            % produce histograms
                            if(opts.saveHistOnly)
                                histDataTemp = data;
                            else
                                doAnalysePlots = false;
                                histDataTemp = analyseDataDiscrete(data,char(name),doAnalysePlots);
                            end

                            % do plotting
                            plot(ha,histDataTemp.mids,histDataTemp.rhoFull,'color',[0,0,0,0.1]);
                            hold(ha,'on');
                            shg;                                        

                            % do experimental post-processing
                            if(Exp.doExperiment)
                                histDataTemp.rhoFullExp = doExperiment(histDataTemp.rhoFull, histDataTemp.mids, Exp.rho0, Exp.sigma,char(name));
                                plot(haExp,histDataTemp.mids,histDataTemp.rhoFullExp,'color',[0,0,0,0.1]);
                                hold(haExp,'on');
                                shg;
                            end

                            histDataTemp.Exp = Exp;
                            histDataTemp.opts = opts;

                            % store data in a structure
                            histData(count) = histDataTemp;
                            
                        end
                    end
                end
            end
        end
    end
end
    
end