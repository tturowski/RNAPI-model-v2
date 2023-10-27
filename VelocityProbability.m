function VelocityProbability(n)
    tMax = 20;
    dt = 0.005;
    nSteps = tMax/dt;
    
    k1 = 80; % pre => post
    kMin1 = 5; % pre <= post
    F1 = 0.4; % back => pre
    B1 = 0.5;% back <= pre
    k3 = 40; % cat

%     k1 = 88;
%     kMin1 = 680;
%     F1 = 1.3;
%     B1 = 6.9;
%     k3 = 35;
    
    % format long
    x = zeros(n,1); %positions
    backtrackingState = zeros(n,1); %saving bactracking state
    postTranslocated = zeros(n,1); %saving post-translocated state

    saveFreq = ceil(nSteps/tMax);
    saveState = zeros(n,1);
    iSave=1;
    xSave = nan(50,tMax);
    
%     figure;
    
    for tStep = 1:nSteps
        prob = rand(size(x));
        
        %analyse current sub-states
        backMask = (backtrackingState < 0);     
        postMask = (postTranslocated > 0);

        %move post-translocated to pre
        kMin1Move = (prob < kMin1*dt);
        postTranslocated(kMin1Move) =  0;
        postMask = (postTranslocated > 0);

        %move from bactracked position
        F1MoveMask = (prob < F1*dt);
        backtrackingState(F1MoveMask) = backtrackingState(F1MoveMask)+1;
%         backtrackingState(F1MoveMask) = max(backtrackingState(F1MoveMask)+1,0);

        %translocation +1 or -1
        B1MoveMask = (prob < B1*dt);
        k1MoveMask = (prob < (B1+k1)*dt);
        
        backtrackingState(~backMask & ~postMask & B1MoveMask) = -1;
        postTranslocated(~backMask & ~postMask & k1MoveMask & ~B1MoveMask) = 1; %should be allowed to ~postMask?
        
        %catalysis
        moveMask = (prob < k3*dt);
%         postMask = (postTranslocated > 0);
        x(postMask & moveMask) = x(postMask & moveMask) + 1;
        
        
%         save every so often
        if(mod(tStep+1,saveFreq)==1 & tStep>1)
            xSave(:,iSave) = histcounts(x - saveState,50);
            saveState = x;
            iSave = iSave + 1;
        end
        
        
%         time = [0,tStep*dt];
%         
%         subplot(2,2,1);
%         plot(time,1+zeros(size(time)));
%         title('time (s)');
%         
%         subplot(2,2,2);
%         histogram((x));
%         title('position');
% 
%         subplot(2,2,3);
%         histogram(backtrackingState);
%         title('bactracking state');
% 
%         subplot(2,2,4);
%         histogram(postTranslocated);
%         title('posttranslocation');
%         drawnow;
        
    end
    
    figure;
       
    subplot(2,2,1);
    histogram((x)./tMax);
    title('velocity (position/s)');
    
    subplot(2,2,2);
    histogram(x);
    title('velocity (position/s)');
    
    subplot(2,2,3);
    histogram(backtrackingState);
    title('bactracking state');
    
    subplot(2,2,4);
    histogram(postTranslocated);
    title('posttranslocation');
    
end