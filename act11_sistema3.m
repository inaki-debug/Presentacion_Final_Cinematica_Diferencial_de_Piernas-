%% Limpieza de entorno
clear all; close all; clc;

%% Parte 1: Modelado cinemático simbólico por sistemas
syms q1 q2 q3 q4 q5 q6 real
syms dq1 dq2 dq3 dq4 dq5 dq6 real
syms px py pz L_femur L_tibia L_pie real 

q_dot = [dq1; dq2; dq3; dq4; dq5; dq6];
GDL = 6; 
SISTEMAS = 7; 
H0 = SE3;

% Matrices locales
H_local_obj(:,:,1) = SE3(rotz(q1), [px py pz]);
H_local_obj(:,:,2) = SE3(rotz(q2), [0 0 0]) * SE3(rotx(-pi/2)*roty(pi/2), [0 0 0]);
H_local_obj(:,:,3) = SE3(rotz(q3), [0 0 0]) * SE3(rotz(-pi/2)*rotx(pi/2), [0 0 0]);
H_local_obj(:,:,4) = SE3(rotz(q4), [0 0 L_femur]);
H_local_obj(:,:,5) = SE3(rotz(q5), [0 0 L_tibia]) * SE3(rotz(pi/2), [0 0 0]);
H_local_obj(:,:,6) = SE3(rotz(q6), [0 0 0]) * SE3(rotx(pi/2)*roty(pi/2), [0 0 0]);
H_local_obj(:,:,7) = SE3(eye(3), [0 0 L_pie]);

T_global = sym(zeros(4,4,SISTEMAS));
v_sistemas = sym(zeros(3, SISTEMAS));
w_sistemas = sym(zeros(3, SISTEMAS));

z0 = [0; 0; 1];
o0 = [0; 0; 0];

%% Parte 2: Configuración de valores numéricos para evaluación
px_n = 0.5; py_n = 0; pz_n = -0.5;  
L_femur_n = 2.0;                  
L_tibia_n = 1.8;                  
L_pie_n = 0.4;                    

q1_n = 0; q2_n = 0; q3_n = 0; q4_n = 0; q5_n = 0; q6_n = 0;
dq1_n = 1; dq2_n = 1; dq3_n = 1; dq4_n = 1; dq5_n = 1; dq6_n = 1;

vars    = [q1, q2, q3, q4, q5, q6, px, py, pz, L_femur, L_tibia, L_pie, ...
           dq1, dq2, dq3, dq4, dq5, dq6];
valores = [q1_n, q2_n, q3_n, q4_n, q5_n, q6_n, px_n, py_n, pz_n, L_femur_n, L_tibia_n, L_pie_n, ...
           dq1_n, dq2_n, dq3_n, dq4_n, dq5_n, dq6_n];

%% Parte 3: Cálculo y evaluación por sistema
for i = 1:SISTEMAS
    % Cálculo de la matriz global acumulada
    if i == 1
        T_global(:,:,i) = simplify(H0.T * H_local_obj(:,:,i).T);
    else
        T_global(:,:,i) = simplify(T_global(:,:,i-1) * H_local_obj(:,:,i).T);
    end
    
    o_actual = T_global(1:3, 4, i);
    juntas_activas = min(i, GDL);
    J_v_i = sym(zeros(3, GDL));
    J_w_i = sym(zeros(3, GDL));
    
    % Construcción del Jacobiano geométrico
    for k = 1:juntas_activas
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
limites_cuadro = [-2 4 -2 4 -5 2];

H0_1_num = SE3(double(subs(T_global(:,:,1), vars, valores)));
H0_2_num = SE3(double(subs(T_global(:,:,2), vars, valores)));
H0_3_num = SE3(double(subs(T_global(:,:,3), vars, valores)));
H0_4_num = SE3(double(subs(T_global(:,:,4), vars, valores)));
H0_5_num = SE3(double(subs(T_global(:,:,5), vars, valores)));
H0_6_num = SE3(double(subs(T_global(:,:,6), vars, valores)));
H0_7_num = SE3(double(subs(T_global(:,:,7), vars, valores)));

p0 = [0; 0; 0];
p1 = H0_1_num.t;
p2 = H0_2_num.t;
p3 = H0_3_num.t;
p4 = H0_4_num.t;
p5 = H0_5_num.t;
p6 = H0_6_num.t;
p7 = H0_7_num.t;

x_robot = [p0(1) p1(1) p4(1) p5(1) p7(1)];
y_robot = [p0(2) p1(2) p4(2) p5(2) p7(2)];
z_robot = [p0(3) p1(3) p4(3) p5(3) p7(3)];

figure();
plot3(x_robot, y_robot, z_robot, 'g-s', 'LineWidth', 3, 'MarkerSize', 6); 
axis(limites_cuadro); grid on;
xlabel('Eje X'); ylabel('Eje Y'); zlabel('Eje Z');
title('Simulación');
hold on;
view(140, 25); 

trplot(H0, 'rgb', 'axis', limites_cuadro, 'frame', '0', 'length', 0.4);

% Ejecución de la animación paso a paso
pause; tranimate(H0, H0_1_num, 'rgb', 'axis', limites_cuadro, 'frame', '1', 'length', 0.4);
pause; tranimate(H0_1_num, H0_2_num, 'rgb', 'axis', limites_cuadro, 'frame', '2', 'length', 0.4);
pause; tranimate(H0_2_num, H0_3_num, 'rgb', 'axis', limites_cuadro, 'frame', '3', 'length', 0.4);
pause; tranimate(H0_3_num, H0_4_num, 'rgb', 'axis', limites_cuadro, 'frame', '4', 'length', 0.4);
pause; tranimate(H0_4_num, H0_5_num, 'rgb', 'axis', limites_cuadro, 'frame', '5', 'length', 0.4);
pause; tranimate(H0_5_num, H0_6_num, 'rgb', 'axis', limites_cuadro, 'frame', '6', 'length', 0.4);
pause; tranimate(H0_6_num, H0_7_num, 'rgb', 'axis', limites_cuadro, 'frame', '7', 'length', 0.4);

fprintf('\nMatriz numérica de transformación final en el pie (H0_7):\n');
disp(H0_7_num.T);