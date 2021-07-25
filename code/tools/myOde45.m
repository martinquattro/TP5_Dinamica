function [tode, X ] = myOde45(M, C, G, x0, tspan, u)
    
    Ts = 1E-2;
    odeOptions = odeset('RelTol',0.001,'AbsTol',0.001,'InitialStep',Ts/20,'MaxStep',Ts);

    [tode,X] = ode45(@odefun, x0,tspan,odeOptions,u);
    
    function xp = odefun(t,x,u)
        q = x(1);
        G = 2.938334137596694943805886168775*cos(q) + 0.000090682051871867628002848830959248*sin(q);
        q_p = x(2);
        q_2p = M^-1*(u - C*q_p - G);
        xp=[q_p ; q_2p];
    end
end