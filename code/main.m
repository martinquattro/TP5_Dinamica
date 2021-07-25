clear all; close all; clc;

% 8615 - Robotica - FIUBA
% TP5 - Dinamica
% Autor: Quattrone Martin y Segura Lola

addpath('tools');
fprintf(' --- Obtencion del modelo dinamico inverso...\n\n')

% Definicion de parametros.
d = 0;
alpha = 0;

syms m a 
syms q q_dot q_2dot

% Matrices de la dinamica.
% -- Matriz de roto-traslacion.
A_01 = [cos(q) -sin(q) 0 a*cos(q);
        sin(q) cos(q) 0 a*sin(q);
        0 0 1 0;
        0 0 0 1];

% -- Centro de masa.
syms xG yG zG
rG = [xG yG zG];

% -- Matriz de pseudoinercia.
syms Ioxx Ioxy Ioxz
syms Ioyx Ioyy Ioyz
syms Iozx Iozy Iozz
J1 = [(-Ioxx+Ioyy+Iozz)/2 -Ioxy -Ioxz xG*m;
       -Ioxy (+Ioxx-Ioyy+Iozz)/2 -Ioyz yG*m;
       -Ioxz -Ioyz (+Ioxx+Ioyy-Iozz)/2 zG*m;
       xG*m yG*m zG*m m];

% -- Matriz de inercia.
dA_01 = diff(A_01, q);   
m11 = simplify(trace(dA_01 * J1 * dA_01.'));
M = m11;
fprintf('M = ')
disp(M)

% -- Matriz de componentes 
% -- de fuerzas Centripeta
% -- y de Coreolis.   
C11 = 1/2 * diff(m11,q) * q_dot;
C = C11;
fprintf('C = ')
disp(C)

% -- Matriz de fuerza o torque 
% -- sobre juntas s por peso
% -- propio.
syms g1
g = [0 -g1 0 1];   
G1 = simplify(-m * g * dA_01 * [rG.'; 1]);
G = G1;
fprintf('G = ')
disp(G)

% -- Calculo del torque.
tau = simplify(m11 * q_2dot + C11 * q_dot + G1);
fprintf('tau = ')
disp(tau)

% -- Estimacion por cuadrados minimos
fprintf('--- Estimacion de parametros...\n\n')
LSestimation
xG_hat = p_hat(1)/m_value;
yG_hat = p_hat(2)/m_value;
Iozz_hat = p_hat(3);

fprintf('La estiamcion resulta \n')
fprintf('Xg = %.3f \n', xG_hat);
fprintf('Yg = %.6f \n', yG_hat);
fprintf('Iozz = %.5f \n', Iozz_hat);


% -- Simulacion
fprintf('\n--- Simulacion...\n\n')
syms q2p qp q tau
q2p = 1/M * (tau - C -G);

fprintf('El modelo directo resulta \nq2p = ')
disp(q2p)

% -- simulador dinamico
M_evaluate = simplify(subs(M, [m a xG yG Iozz], [m_value a_value xG_hat yG_hat Iozz_hat]));
C_evaluate = simplify(subs(C, [m a xG yG Iozz], [m_value a_value xG_hat yG_hat Iozz_hat]));
G_evaluate = vpa(subs(G, [m a xG yG Iozz g1], [m_value a_value xG_hat yG_hat Iozz_hat g_value]));

u = 0;
x0 = [-pi/2 , 0];
tspan = [0 5];
Ts = 1E-2;

% [tode,X] = myOde45(M_evaluate, C_evaluate, G_evaluate, x0, tspan, u);
[tode,X] = ode45(@odefun, x0,tspan,u);
    
figure(1)
plot(tode, X(:,2))

function xp = odefun(t,x,u)
    q = x(1);
    q_p = x(2);
    G = 2.938334137596694943805886168775*cos(q) + 0.000090682051871867628002848830959248*sin(q);
    M = 1325450956090428077/28823037615171174400;
    C = 0;
    q_2p = M^-1*(u - C*q_p - G);
    xp=[q_p ; q_2p];
end






























