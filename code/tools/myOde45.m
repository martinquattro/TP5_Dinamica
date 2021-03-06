function [tode, X] = myOde45(tspan, x0, odeOptions, u, b_value)
     
    [tode,X] = ode45(@odefun, tspan, x0, odeOptions, u);

    function xp = odefun(t,x,u)
        q = x(1);
        q_p = x(2);
        
        G = 2.938334*cos(q) + 0.00009068205*sin(q);
        M = 0.04598582;
        C = 0.0;
        
        q_2p = M^-1*(u - b_value*q_p - C*q_p - G);
        xp=[q_p ; q_2p];
    end
end