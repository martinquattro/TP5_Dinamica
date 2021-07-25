% -- Estimacion de los parametros dinamicos


datos_pendulo

m_value =     2;
a_value =     0.2;
g_value =     9.8;
q_value =     datos(:,2);
q_dot_value =   datos(:,3);
q_2dot_value =  datos(:,4);
tau_value = datos(:,5);

phi_kn = [ q_2dot_value*a_value^2 + g_value*a_value*cos(q_value) ];
phi_un = [ 2*q_2dot_value*a_value + g_value*cos(q_value) , -g_value*sin(q_value) , q_2dot_value];

fprintf('El numero de condicion de la matriz es: %.2f \n', cond(phi_un))
p_hat = pinv(phi_un) * (tau_value - phi_kn * m_value);

