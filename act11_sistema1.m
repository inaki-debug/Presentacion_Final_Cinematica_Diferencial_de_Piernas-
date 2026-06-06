%% Limpieza de entorno
clear all; close all; clc;

%% Parte 1: Modelado cinemático simbólico por sistemas
syms theta1 theta2 theta3 theta4
syms dtheta1 dtheta2 dtheta3 dtheta4
syms L1 L2 L3 L4

q_dot = [dtheta1; dtheta2; dtheta3; dtheta4];
GDL = 4;
H0 = SE3; 

% Matrices de transformación relativa de cada eslabón
h1_obj = SE3(rotz(theta1), [0 0 L1]) * SE3(rotx(pi/2), [0 0 0]);
H_local(:,:,1) = h1_obj.T;

h2_obj = SE3(rotz(theta2), [L2 0 0]);
H_local(:,:,2) = h2_obj.T;

h3_obj = SE3(rotz(theta3), [L3 0 0]);
H_local(:,:,3) = h3_obj.T;

h4_obj = SE3(rotz(theta4), [L4 0 0]);
H_local(:,:,4) = h4_obj.T;

T_global = sym(zeros(4,4,GDL));
v_sistemas = sym(zeros(3, GDL));
w_sistemas = sym(zeros(3, GDL));

z0 = [0; 0; 1];
o0 = [0; 0; 0];

%% Parte 2: Configuración de valores numéricos para evaluación
L1_num = 2; L2_num = 1; L3_num = 1; L4_num = 1;
theta1_num = 0; theta2_num = 0; theta3_num = 0; theta4_num = 0;
dtheta1_num = 1; dtheta2_num = 1; dtheta3_num = 1; dtheta4_num = 1;

vars    = [theta1, theta2, theta3, theta4, L1, L2, L3, L4, dtheta1, dtheta2, dtheta3, dtheta4];
valores = [theta1_num, theta2_num, theta3_num, theta4_num, L1_num, L2_num, L3_num, L4_num, dtheta1_num, dtheta2_num, dtheta3_num, dtheta4_num];

%% Parte 3: Cálculo y evaluación por sistema
for i = 1:GDL
    % Cálculo de la matriz de transformación global acumulada
    if i == 1
        T_global(:,:,i) = simplify(H0.T * H_local(:,:,i));
    else
        T_global(:,:,i) = simplify(T_global(:,:,i-1) * H_local(:,:,i));
    end
    
    o_actual = T_global(1:3, 4, i);
    J_v_i = sym(zeros(3, GDL));
    J_w_i = sym(zeros(3, GDL));
    
    % Construcción del Jacobiano geométrico
    for k = 1:i
        if k == 1
            z_ant = z0;
            o_ant = o0;
        else
            z_ant = T_global(1:3, 3, k-1);
            o_ant = T_global(1:3, 4, k-1);
        end
        
        J_v_i(:, k) = cross(z_ant, (o_actual - o_ant));
        J_w_i(:, k) = z_ant;
    end
    
    v_sistemas(:, i) = simplify(J_v_i * q_dot);
    w_sistemas(:, i) = simplify(J_w_i * q_dot);
    
    % Despliegue de resultados en consola
    fprintf('  Sistema / Eslabón %d\n', i)
    
    fprintf('Matriz homogénea local simbólica (H%d):\n', i)
    disp(simplify(H_local(:,:,i)))
    fprintf('Matriz homogénea local numérica (H%d) a 0°:\n', i)
    disp(double(subs(H_local(:,:,i), vars, valores)))
    
    fprintf('Matriz global acumulada simbólica (T%d):\n', i)
    disp(T_global(:,:,i))
    fprintf('Matriz global acumulada numérica (T%d) a 0°:\n', i)
    disp(double(subs(T_global(:,:,i), vars, valores)))
    
    fprintf('Vector de velocidad lineal simbólica (v%d):\n', i)
    disp(v_sistemas(:, i))
    fprintf('Vector de velocidad lineal numérica (v%d) [dq = 1 rad/s]:\n', i)
    disp(double(subs(v_sistemas(:, i), vars, valores)))
    
    fprintf('Vector de velocidad angular simbólica (w%d):\n', i)
    disp(w_sistemas(:, i))
    fprintf('Vector de velocidad angular numérica (w%d) [dq = 1 rad/s]:\n', i)
    disp(double(subs(w_sistemas(:, i), vars, valores)))
    fprintf('\n')
end

%% Parte 4: Simulación gráfica
H0_1_num = SE3(double(subs(T_global(:,:,1), vars, valores)));
H0_2_num = SE3(double(subs(T_global(:,:,2), vars, valores)));
H0_3_num = SE3(double(subs(T_global(:,:,3), vars, valores)));
H0_4_num = SE3(double(subs(T_global(:,:,4), vars, valores)));

p0 = [0; 0; 0];
p1 = H0_1_num.t;
p2 = H0_2_num.t; 
p3 = H0_3_num.t;
p4 = H0_4_num.t;

x_robot = [p0(1) p1(1) p2(1) p3(1) p4(1)];
y_robot = [p0(2) p1(2) p2(2) p3(2) p4(2)];
z_robot = [p0(3) p1(3) p2(3) p3(3) p4(3)];

figure();
plot3(x_robot, y_robot, z_robot, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 6); 
axis([-1 4 -1 4 -1 4]); grid on;
xlabel('Eje X'); ylabel('Eje Y'); zlabel('Eje Z');
title('Simulación');
hold on;
view(135, 30); 

trplot(H0, 'rgb', 'axis', [-1 4 -1 4 -1 4], 'frame', '0', 'length', 0.4);

% Ejecución de la animación paso a paso
pause; tranimate(H0, H0_1_num, 'rgb', 'axis', [-1 4 -1 4 -1 4], 'frame', '1', 'length', 0.4);
pause; tranimate(H0_1_num, H0_2_num, 'rgb', 'axis', [-1 4 -1 4 -1 4], 'frame', '2', 'length', 0.4);
pause; tranimate(H0_2_num, H0_3_num, 'rgb', 'axis', [-1 4 -1 4 -1 4], 'frame', '3', 'length', 0.4);
pause; tranimate(H0_3_num, H0_4_num, 'rgb', 'axis', [-1 4 -1 4 -1 4], 'frame', '4', 'length', 0.4);