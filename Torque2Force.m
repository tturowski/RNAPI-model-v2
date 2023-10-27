function Torque2Force
    % Ma et al Science 2013 Supp Data: force * d = Trq * theta
    % d represents the contour length of DNA per bp (~0.34 nm)
    % Theta represents the angular rotation of RNAP after 1 bp translocation 
    % (0.60 radian or ~34°, converted from 10.5 bp per turn).
    
    Trq = [0:30];
    force = (Trq .* 0.6) ./ 0.34;
    
    figure;
    plot(Trq,force,'-o');
    grid on
    grid minor
    title(['Force calculated from Torque \tau']);
    xlabel('Torque \tau [pN nm]');
    ylabel('Force [pN]');
end