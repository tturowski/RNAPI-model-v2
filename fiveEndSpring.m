function fiveEndSpring(dist)
    
    %set up array
    x = [1:7000];
    v = zeros(size(x))+100;
   
    %distance of the cut off and StartMask
    startMask = (x<dist);
    
    %mod DNA spring
    v_mod1 = v;
    v_mod1(startMask) = v(startMask)./(dist-x(startMask));
    
    v_mod2 = v;
    v_mod2(startMask) = v_mod2(startMask) .* x(startMask) ./ dist;
    
    v_mod3 = v;
    v_mod3(startMask) = v_mod3(startMask) .* exp((x(startMask) - dist));
    
    v_mod4 = v;
    v_mod4(startMask) = v_mod4(startMask) .* nthroot((x(startMask) ./ dist),2);
  
    v_mod5 = v;
    v_mod5(startMask) = v_mod5(startMask) .* nthroot((x(startMask) ./ dist),12);
       
    %plotting
    figure;
    hold on
    plot(x,v);
%     plot(x,v_mod1);
    plot(x,v_mod2);
%     plot(x,v_mod3);
    plot(x,v_mod4);
    plot(x,v_mod5);
end