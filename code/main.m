clear all; close all; clc;

% 8615 - Robotica - FIUBA
% TP5 - Dinamica
% Autor: Quattrone Martin y Segura Lola

addpath('tools');

% =========================================================================
% Obtencion del modilo dinamico inverso
% =========================================================================
fprintf('--- Obtencion del modelo dinamico inverso...\n\n')

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
fprintf('M = %s\n', M)

% -- Matriz de componentes 
% -- de fuerzas Centripeta
% -- y de Coreolis.   
C11 = 1/2 * diff(m11,q) * q_dot;
C = C11;
fprintf('C = %s\n', C)

% -- Matriz de fuerza o torque 
% -- sobre juntas s por peso
% -- propio.
syms g1
g = [0 -g1 0 1];   
G1 = simplify(-m * g * dA_01 * [rG.'; 1]);
G = G1;
fprintf('G = %s\n', G)

% -- Calculo del torque.
tau = simplify(m11 * q_2dot + C11 * q_dot + G1);
fprintf('\nTAU = %s\n\n ', tau)

% =========================================================================
% Estimacion por cuadrados minimos
% =========================================================================

fprintf('--- Estimacion de parametros...\n\n')
LSestimation
xG_hat = p_hat(1)/m_value;
yG_hat = p_hat(2)/m_value;
Iozz_hat = p_hat(3);

fprintf('La estimacion de los parametros dinamicos desconocidos resulta \n')
fprintf('Xg = %.3f \n', xG_hat);
fprintf('Yg = %.6f \n', yG_hat);
fprintf('Iozz = %.5f \n', Iozz_hat);

% =========================================================================
% Simulacion
% =========================================================================

fprintf('\n--- Simulacion...\n\n')
syms q2p qp q tau b
q2p = 1/M * (tau - b*qp - C - G);

fprintf('El modelo directo resulta \nq2p = %s\n', q2p)

M_evaluate = vpa(subs(M, [m a xG yG Iozz], [m_value a_value xG_hat yG_hat Iozz_hat]),7)
C_evaluate = vpa(subs(C, [m a xG yG Iozz], [m_value a_value xG_hat yG_hat Iozz_hat]),7)
G_evaluate = vpa(subs(G, [m a xG yG Iozz g1], [m_value a_value xG_hat yG_hat Iozz_hat g_value]),7)

u = [0.1, 2, 3];                % torque aplicados a simular [Nm]
b_values = [0 , 0.1];           % coeficiente viscoso [Nm/rad/s]
q_eq = atan((a_value+xG_hat)/yG_hat);   
x0 = [ q_eq, 0];                % punto de equilibrio etable
tspan = [0 5];                  % tiempo de simulacion
Ts = 1E-3;
odeOptions = odeset('RelTol',0.001,'AbsTol',0.001,'InitialStep',Ts/20,'MaxStep',Ts);

% simulacion de trayectoria para tau = 0.1 Nm con y sin termino disipativo
[tode11,X11] = myOde45(tspan, x0, odeOptions, u(1), b_values(1));
[tode12,X12] = myOde45(tspan, x0, odeOptions, u(1), b_values(2));

% simulacion de trayectoria para tau = 2 Nm con y sin termino disipativo
[tode21,X21] = myOde45(tspan, x0, odeOptions, u(2), b_values(1));
[tode22,X22] = myOde45(tspan, x0, odeOptions, u(2), b_values(2));

% simulacion de trayectoria para tau = 3 Nm con y sin termino disipativo
[tode31,X31] = myOde45(tspan, x0, odeOptions, u(3), b_values(1));
[tode32,X32] = myOde45(tspan, x0, odeOptions, u(3), b_values(2));

% =========================================================================
% Graficos
% =========================================================================

plotTrajectory(tode11, X11(:,1), tode12, X12(:,1), 1)
plotTrajectory(tode21, X21(:,1), tode22, X22(:,1), 2)
plotTrajectory(tode31, X31(:,1), tode32, X32(:,1), 3)































