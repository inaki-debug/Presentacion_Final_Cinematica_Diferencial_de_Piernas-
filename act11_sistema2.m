%% Limpieza de entorno
clear all; close all; clc;

%% Parte 1: Modelado cinemático simbólico por sistemas
syms theta1 theta2 theta3 theta4 theta5 real
syms dtheta1 dtheta2 dtheta3 dtheta4 dtheta5 real
syms L2 L3 real 

q_dot = [dtheta1; dtheta2; dtheta3; dtheta4; dtheta5];
GDL = 5;
H0 = SE3; 

% Matrices locales simbólicas puras
R1 = rotz(theta1);
H_local(:,:,1) = [R1, [0; 0; 0]; 0 0 0 1];

R2 = rotz(theta2) * rotx(-pi/2);
H_local(:,:,2) = [R2, [0; 0; 0]; 0 0 0 1];

R3 = rotz(theta3);
H_local(:,:,3) = [R3, [L2; 0; 0]; 0 0 0 1];

R4 = rotz(theta4 - pi/2);
H_local(:,:,4) = [R4, [L3; 0; 0]; 0 0 0 1];

R5 = rotz(theta5) * rotx(-pi/2);
H_local(:,:,5) = [R5, [0; 0; 0]; 0 0 0 1];

T_global = sym(zeros(4,4,GDL));
v_sistemas = sym(zeros(3, GDL));
w_sistemas = sym(zeros(3, GDL));

z0 = [0; 0; 1];
o0 = [0; 0; 0];

%% Parte 2: Configuración de valores numéricos para evaluación
L2_num = 2.5; 
L3_num = 2.0; 

theta1_num = 0; theta2_num = 0; theta3_num = 0; theta4_num = 0; theta5_num = 0;
dtheta1_num = 1; dtheta2_num = 1; dtheta3_num = 1; dtheta4_num = 1; dtheta5_num = 1;

vars    = [theta1, theta2, theta3, theta4, theta5, L2, L3, ...
           dtheta1, dtheta2, dtheta3, dtheta4, dtheta5];
valores = [theta1_num, theta2_num, theta3_num, theta4_num, theta5_num, L2_num, L3_num, ...
           dtheta1_num, dtheta2_num, dtheta3_num, dtheta4_num, dtheta5_num];

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
    
    fprintf('Matriz de transformación homogénea local simbólica (H%d):\n', i)
    disp(simplify(H_local(:,:,i)))
    fprintf('Matriz de transformación homogénea local numérica (H%d):\n', i)
    disp(double(subs(H_local(:,:,i), vars, valores)))
    
    fprintf('Matriz de transformación homogénea global simbólica (T%d):\n', i)
    disp(T_global(:,:,i))
    fprintf('Matriz de transformación homogénea global numérica (T%d):\n', i)
    disp(double(subs(T_global(:,:,i), vars, valores)))
    
    fprintf('Vector de velocidad lineal simbólica (v%d):\n', i)
    disp(v_sistemas(:, i))
    fprintf('Vector de velocidad lineal numérica (v%d):\n', i)
    disp(double(subs(v_sistemas(:, i), vars, valores)))
    
    fprintf('Vector de velocidad angular simbólica (w%d):\n', i)
    disp(w_sistemas(:, i))
    fprintf('Vector de velocidad angular numérica (w%d):\n', i)
    disp(double(subs(w_sistemas(:, i), vars, valores)))
    fprintf('\n')
end

%% Parte 4: Simulación gráfica
limites_cuadro = [-1 6 -3 3 -5 1];

H0_1_num = SE3(double(subs(T_global(:,:,1), vars, valores)));
H0_2_num = SE3(double(subs(T_global(:,:,2), vars, valores)));
H0_3_num = SE3(double(subs(T_global(:,:,3), vars, valores)));
H0_4_num = SE3(double(subs(T_global(:,:,4), vars, valores)));
H0_5_num = SE3(double(subs(T_global(:,:,5), vars, valores)));

p0 = [0; 0; 0];
p1 = H0_1_num.t;
p2 = H0_2_num.t; 
p3 = H0_3_num.t;
p4 = H0_4_num.t;
p5 = H0_5_num.t; 

x_robot = [p0(1) p1(1) p2(1) p3(1) p4(1) p5(1)];
y_robot = [p0(2) p1(2) p2(2) p3(2) p4(2) p5(2)];
z_robot = [p0(3) p1(3) p2(3) p3(3) p4(3) p5(3)];

figure();
plot3(x_robot, y_robot, z_robot, 'b-o', 'LineWidth', 3, 'MarkerSize', 7); 
axis(limites_cuadro); grid on;
xlabel('Eje X'); ylabel('Eje Y'); zlabel('Eje Z');
title('Simulación');
hold on;
view(120, 20); 

trplot(H0, 'rgb', 'axis', limites_cuadro, 'frame', '0', 'length', 0.5);

% Ejecución de la animación paso a paso
pause; tranimate(H0, H0_1_num, 'rgb', 'axis', limites_cuadro, 'frame', '1', 'length', 0.5);
pause; tranimate(H0_1_num, H0_2_num, 'rgb', 'axis', limites_cuadro, 'frame', '2', 'length', 0.5);
pause; tranimate(H0_2_num, H0_3_num, 'rgb', 'axis', limites_cuadro, 'frame', '3', 'length', 0.5);
pause; tranimate(H0_3_num, H0_4_num, 'rgb', 'axis', limites_cuadro, 'frame', '4', 'length', 0.5);
pause; tranimate(H0_4_num, H0_5_num, 'rgb', 'axis', limites_cuadro, 'frame', '5', 'length', 0.5);