% function Torque2Velocity(stalForce, vInt)
function Torque2Velocity()

    TrqMax = 20;
    Trq = [0:0.01:20];
    k = 0.7; %steepness
    vInt = 22;

    figure;
    hold on

    v = vInt + log(TrqMax ./ Trq - 1) ./ k;
    plot(Trq-10,v,'-o');

    t=[-10,-6,-2.5,0,5,7.5,10];
    vv=[45,24,23,22,21,18,0];
    plot(t,vv,'o');
    
    xlabel('Torque \tau [pN nm]');
    ylabel('Velocity v [bp/s]');
    grid on
    grid minor
end

%DNA melting Torque = -11 +- 1 pN·nm
%RNAP stall Torque 11 pN·nm or 14 pN·nm (= 25 pN)