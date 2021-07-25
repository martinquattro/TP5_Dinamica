% -- Estimacion de los parametros dinamicos

fprintf('Estimacion de parametros...\n\n')
datos_pendulo

m =     2;
a =     0.2;
g =     9.8;
q =     datos(:,2);
q_dot =   datos(:,3);
q_2dot =  datos(:,4);
tau = datos(:,5);

phi_kn = [ q_2dot*a^2 + g*a*cos(q) ];
phi_un = [ 2*q_2dot*a + g*cos(q) , -g*sin(q) , q_2dot];

fprintf('El numero de condicion de la matriz es: %.2f \n', cond(phi_un))
p_hat = pinv(phi_un) * (tau - phi_kn * m);

fprintf('La estiamcion resulta \n')
fprintf('Xg = %.3f \n', p_hat(1)/m);
fprintf('Yg = %.3f \n', p_hat(2)/m);
fprintf('Iozz = %.3f \n', p_hat(3));