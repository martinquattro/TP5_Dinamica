clear all; close all; clc;

% 8615 - Robotica - FIUBA
% TP5 - Dinamica
% Autor: Quattrone Martin y Segura Lola

addpath('tools');
fprintf('Obtencion del modelo dinamico inverso...\n\n')

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

% -- Matriz de fuerza o torque 
% -- sobre juntas s por peso
% -- propio.
syms g1
g = [0 -g1 0 1];   
G1 = simplify(-m * g * dA_01 * [rG.'; 1]);

% -- Matriz de componentes 
% -- de fuerzas Centripeta
% -- y de Coreolis.   
C11 = 1/2 * diff(m11,q) * q_dot;

% Calculo del torque.
tau = simplify(m11 * q_2dot + C11 * q_dot + G1);
fprintf('tau = \n')
disp(tau)

% Estimacion por cuadrados minimos
datos_pendulo
LSestimation





























