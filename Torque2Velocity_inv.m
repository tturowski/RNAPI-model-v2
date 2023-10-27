% function Torque2Velocity(stalForce, vInt)
function Torque2Velocity_inv()

    TrqMax = 20;
    Trq = [0:0.01:10];
    k = 0.5; %steepness

    vInt = 20;

    figure;
    
%     t=[-10,-6,-2.5,0,5,7.5,10];
%     vv=[45,24,23,22,21,18,0];
%     plot(t,vv,'o');
    
    v = [0:40];
    Trq = TrqMax ./ (1+exp(k.*(v-vInt)));

%     Trq(1) = 0.0001;
%     Trq(end) = 0.9999;
%     v = log((Trq-1) ./ Trq);
   
%     v = log(Trq./(1-Trq));
%     v = log(1./Trq) ./ 0.1;

    plot(Trq-10,v,'-o');
    
    xlabel('Torque \tau [pN nm]');
    ylabel('Velocity v [bp/s]');
    grid on
    grid minor
end

%DNA melting Torque = -11 +- 1 pN·nm
%RNAP stall Torque 11 pN·nm or 14 pN·nm (= 25 pN)