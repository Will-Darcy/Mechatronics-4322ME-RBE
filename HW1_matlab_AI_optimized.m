% Homework 1: Six-Bar Linkage Optimization (Numerical Matrix Formulation)
clc;
clear;

%% ------------------------------------------------------------------
%% 1. Geometric Constants & Parameters
%% ------------------------------------------------------------------
A = [1.4,   0.485, 0];
B = [1.67,  0.99,  0];
C = [0.255, 1.035, 0];
D = [0.285, 0.055, 0];
E = [0.195, 2.54,  0];
F = [-0.98, 2.57,  0];
G = [0.05,  0.2,   0];

L_FH = 1.83;
lAB  = 0.5726;
lBC  = 1.4157;
lDE  = 2.4866;
lEF  = 1.1754;
lGF  = 2.5841;
lCD  = norm(D - C);

% Weights (N) & Applied Force
W_vecs.AB = [0, -228.865,  0];
W_vecs.BC = [0, -553.383,  0];
W_vecs.DE = [0, -954.701,  0];
W_vecs.EF = [0, -460.889,  0];
W_vecs.FG = [0, -1716.634, 0];
Fapp      = [0, -200,     0];

% Link Masses (kg)
masses.AB = 23.3377;
masses.BC = 56.42938;
masses.DE = 97.35243;
masses.EF = 46.9976;
masses.FG = 175.04798;

% Inertia calculation helper
calc_IG = @(pt1, pt2, m) (1/12) * m * ((norm(pt2 - pt1) + 0.1)^2 + 0.1^2);

% Compute initial H to find FG inertia
dGF_init = hypot(F(1) - G(1), F(2) - G(2));
H_init   = [F(1) + L_FH * (F(1) - G(1)) / dGF_init, ...
            F(2) + L_FH * (F(2) - G(2)) / dGF_init, 0];

% Mass moments of inertia (constant for rigid bodies)
IG.AB = calc_IG(A, B, masses.AB);
IG.BC = calc_IG(B, C, masses.BC);
IG.DE = calc_IG(D, E, masses.DE);
IG.EF = calc_IG(E, F, masses.EF);
IG.FG = calc_IG(H_init, G, masses.FG);

%% ------------------------------------------------------------------
%% 2. Preallocation & Storage Setup
%% ------------------------------------------------------------------
num_joints = 7;
num_links  = 5;
num_angles = 361;

positions           = NaN(2, num_joints + 1, num_angles);
Forces              = NaN(2, num_joints, num_angles);
DynamicSol          = NaN(2, num_joints, num_angles);
angularVel          = NaN(1, num_links, num_angles);
angularAcc          = NaN(1, num_links, num_angles);
InputTorque_static  = NaN(1, num_angles);
InputTorque_dynamic = NaN(1, num_angles);

% Ground pivots (A=col 1, D=col 4, G=col 7) never move
positions(:, 1, :) = repmat(A(1:2)', 1, 1, num_angles);
positions(:, 4, :) = repmat(D(1:2)', 1, 1, num_angles);
positions(:, 7, :) = repmat(G(1:2)', 1, 1, num_angles);

initial_theta = atan2(B(2) - A(2), B(1) - A(1));
if initial_theta < 0
    NewAngle = 2*pi + initial_theta;
else
    NewAngle = initial_theta;
end

% Storage for Center of Mass kinematics: [X/Y, COM 1:5, angle]
vel_COM = NaN(2, num_links, num_angles);
acc_COM = NaN(2, num_links, num_angles);

% Continuous branch tracking references
C_prev = C;
F_prev = F;

%% ------------------------------------------------------------------
%% 3. Main Numerical Sweep
%% ------------------------------------------------------------------
tic;
for theta = 1:num_angles
    % Offset by 1: theta = 1 corresponds to 0 degrees rotation
    theta_deg = theta - 1;
    theta_rad = NewAngle + deg2rad(theta_deg);

    % 1. Kinematic Loop: Joint Positions
    B_new = A + [lAB * cos(theta_rad), lAB * sin(theta_rad), 0];

    [Cx, Cy] = circcirc(B_new(1), B_new(2), lBC, D(1), D(2), lCD);
    if isnan(Cx(1)), continue; end

    C_1 = [Cx(1), Cy(1), 0];
    C_2 = [Cx(2), Cy(2), 0];
    if norm(C_1 - C_prev) < norm(C_2 - C_prev)
        C_new = C_1;
    else
        C_new = C_2;
    end
    C_prev = C_new;

    D_theta = atan2(C_new(2) - D(2), C_new(1) - D(1));
    E_new   = D + [lDE * cos(D_theta), lDE * sin(D_theta), 0];

    [Fx, Fy] = circcirc(E_new(1), E_new(2), lEF, G(1), G(2), lGF);
    if isnan(Fx(1)), continue; end

    F_1 = [Fx(1), Fy(1), 0];
    F_2 = [Fx(2), Fy(2), 0];
    if norm(F_1 - F_prev) < norm(F_2 - F_prev)
        F_new = F_1;
    else
        F_new = F_2;
    end
    F_prev = F_new;

    dx = F_new(1) - G(1);
    dy = F_new(2) - G(2);
    d_GF = hypot(dx, dy);
    H_new = [F_new(1) + L_FH * (dx / d_GF), F_new(2) + L_FH * (dy / d_GF), 0];

    % Store positions
    positions(:, 2, theta) = B_new(1:2)';
    positions(:, 3, theta) = C_new(1:2)';
    positions(:, 5, theta) = E_new(1:2)';
    positions(:, 6, theta) = F_new(1:2)';
    positions(:, 8, theta) = H_new(1:2)';

    % Centers of mass
    S.S1 = (A + B_new) / 2;
    S.S2 = (B_new + C_new) / 2;
    S.S3 = (D + E_new) / 2;
    S.S4 = (E_new + F_new) / 2;
    S.S5 = (H_new + G) / 2;

    pts.A = A; pts.B = B_new; pts.C = C_new; pts.D = D;
    pts.E = E_new; pts.F = F_new; pts.G = G;

    % 2. Angular Velocity Loops (2x2 Matrix Inversions)
    w_AB = 2.424; % rad/s
    % Loop 1: w_AB*(B-A) + w_BC*(C-B) + w_DE*(D-C) = 0
    rCB = C_new - B_new;  rDC = D - C_new;  rBA = B_new - A;
    J1 = [-rCB(2), -rDC(2); 
           rCB(1),  rDC(1)];
    b_vel1 = -w_AB * [-rBA(2); rBA(1)];
    w_sol1 = J1 \ b_vel1;
    w_BC = w_sol1(1); 
    w_DE = w_sol1(2);

    % Loop 2: w_DE*(E-D) + w_EF*(F-E) + w_FG*(G-F) = 0
    rFE = F_new - E_new;  rGF = G - F_new;  rED = E_new - D;
    J2 = [-rFE(2), -rGF(2); 
           rFE(1),  rGF(1)];
    b_vel2 = -w_DE * [-rED(2); rED(1)];
    w_sol2 = J2 \ b_vel2;
    w_EF = w_sol2(1); 
    w_FG = w_sol2(2);

    angularVel(1, :, theta) = [w_AB, w_BC, w_DE, w_EF, w_FG];

    % 3. Angular Acceleration Loops (2x2 Matrix Inversions)
    % Loop 1
    a_norm_AB = -w_AB^2 * rBA;
    a_norm_BC = -w_BC^2 * rCB;
    a_norm_DE = -w_DE^2 * rDC;
    b_acc1    = -(a_norm_AB(1:2)' + a_norm_BC(1:2)' + a_norm_DE(1:2)');
    a_sol1    = J1 \ b_acc1;
    a_BC      = a_sol1(1); 
    a_DE      = a_sol1(2);

    % Loop 2 (exact transfer to link DCE)
    rCD = C_new - D;
    a_D_E = cross([0, 0, a_DE], rCD) + cross([0, 0, w_DE], cross([0, 0, w_DE], rCD));
    a_norm_EF = -w_EF^2 * rFE;
    a_norm_FG = -w_FG^2 * rGF;
    b_acc2    = -(a_D_E(1:2)' + a_norm_EF(1:2)' + a_norm_FG(1:2)');
    a_sol2    = J2 \ b_acc2;
    a_EF      = a_sol2(1); 
    a_FG      = a_sol2(2);

    angularAcc(1, :, theta) = [0, a_BC, a_DE, a_EF, a_FG];

    % 4. Center-of-Mass Accelerations
    w_v.AB = [0, 0, w_AB]; w_v.BC = [0, 0, w_BC]; w_v.DE = [0, 0, w_DE];
    w_v.EF = [0, 0, w_EF]; w_v.FG = [0, 0, w_FG];

    a_v.AB = [0, 0, 0];    a_v.BC = [0, 0, a_BC]; a_v.DE = [0, 0, a_DE];
    a_v.EF = [0, 0, a_EF]; a_v.FG = [0, 0, a_FG];

    accCOM.S1 = cross(w_v.AB, cross(w_v.AB, S.S1 - A));
    aB_A      = cross(w_v.AB, cross(w_v.AB, B_new - A));
    aS2_B     = cross(a_v.BC, S.S2 - B_new) + cross(w_v.BC, cross(w_v.BC, S.S2 - B_new));
    accCOM.S2 = aS2_B + aB_A;
    accCOM.S3 = cross(a_v.DE, S.S3 - D) + cross(w_v.DE, cross(w_v.DE, S.S3 - D));
    aF_G      = cross(a_v.FG, F_new - G) + cross(w_v.FG, cross(w_v.FG, F_new - G));
    aS4_F     = cross(a_v.EF, S.S4 - F_new) + cross(w_v.EF, cross(w_v.EF, S.S4 - F_new));
    accCOM.S4 = aS4_F + aF_G;
    accCOM.S5 = cross(a_v.FG, S.S5 - G) + cross(w_v.FG, cross(w_v.FG, S.S5 - G));

    % Compute COM Linear Velocities (Planar 2D)
    vB_GR = cross(w_v.AB, B_new - A);
    vF_GR = cross(w_v.FG, F_new - G);

    vS1 = cross(w_v.AB, S.S1 - A);
    vS2 = cross(w_v.BC, S.S2 - B_new) + vB_GR;
    vS3 = cross(w_v.DE, S.S3 - D);
    vS4 = cross(w_v.EF, S.S4 - F_new) + vF_GR;
    vS5 = cross(w_v.FG, S.S5 - G);

    % Store [Vx; Vy] and [Ax; Ay] for S1 through S5
    vel_COM(:, 1, theta) = vS1(1:2)';
    vel_COM(:, 2, theta) = vS2(1:2)';
    vel_COM(:, 3, theta) = vS3(1:2)';
    vel_COM(:, 4, theta) = vS4(1:2)';
    vel_COM(:, 5, theta) = vS5(1:2)';

    acc_COM(:, 1, theta) = accCOM.S1(1:2)';
    acc_COM(:, 2, theta) = accCOM.S2(1:2)';
    acc_COM(:, 3, theta) = accCOM.S3(1:2)';
    acc_COM(:, 4, theta) = accCOM.S4(1:2)';
    acc_COM(:, 5, theta) = accCOM.S5(1:2)';

    % 5. Linear Static & Dynamic Force Equilibrium (15x15 System)
    [F_static,  T_static]  = solve_linkage_forces(pts, S, H_new, W_vecs, Fapp, masses, IG, accCOM, a_v, false);
    [F_dynamic, T_dynamic] = solve_linkage_forces(pts, S, H_new, W_vecs, Fapp, masses, IG, accCOM, a_v, true);

    Forces(:, :, theta)         = F_static;
    DynamicSol(:, :, theta)     = F_dynamic;
    InputTorque_static(theta)   = T_static;
    InputTorque_dynamic(theta)  = T_dynamic;
end
elapsed = toc;
fprintf('Optimized simulation completed in: %.4f seconds\n', elapsed);

%% ------------------------------------------------------------------
%% 4. Plotting (With Maximum Value Annotations)
%% ------------------------------------------------------------------
jointLabels = {'A','B','C','D','E','F','G'};
theta_axis  = 0:(num_angles - 1); % Sweeps 0 deg to 360 deg

% -------------------------------------------------------------------
% Plot 1: Coupler Curves
% -------------------------------------------------------------------
figure('Color','w','Name','Joint Traces');
hold on; axis equal; grid on;
jointCols = [2, 3, 5, 6, 8];
markerColors = lines(5);
names = {'B','C','E','F','H'};
for k = 1:5
    col = jointCols(k);
    plot(squeeze(positions(1, col, :)), squeeze(positions(2, col, :)), ...
         '.', 'Color', markerColors(k,:), 'DisplayName', names{k});
end
plot([A(1), D(1), G(1)], [A(2), D(2), G(2)], 'ks', 'MarkerFaceColor','k', 'DisplayName','Ground pivots');
legend('Location','bestoutside'); xlabel('X (m)'); ylabel('Y (m)');
title('Joint Path Traces Over Full Rotation');

% -------------------------------------------------------------------
% Plot 2: Force Magnitudes (Peak Dynamic Force Annotated per Joint)
% -------------------------------------------------------------------
figure('Color','w','Name','Force Magnitudes');
for j = 1:num_joints
    subplot(3,3,j);
    
    mag_static  = hypot(squeeze(Forces(1,j,:)), squeeze(Forces(2,j,:)));
    mag_dynamic = hypot(squeeze(DynamicSol(1,j,:)), squeeze(DynamicSol(2,j,:)));
    
    plot(theta_axis, mag_static, 'b-', 'LineWidth', 1.2); hold on;
    plot(theta_axis, mag_dynamic, 'r--', 'LineWidth', 1.2);
    grid on;
    
    % Annotate peak dynamic force for this joint
    annotate_max(theta_axis, mag_dynamic, [0.85 0 0], 'N');
    
    title(['Joint ' jointLabels{j}]); xlabel('\theta (deg)'); ylabel('|F| (N)');
    if j == 1, legend('Static','Dynamic','Location','best'); end
end
sgtitle('Joint Force Magnitudes vs Angle');

% -------------------------------------------------------------------
% Plot 3: Input Torque (Peak Static & Dynamic Annotated)
% -------------------------------------------------------------------
figure('Color','w','Name','Input Torque');
hold on; grid on;
plot(theta_axis, InputTorque_static, 'b-', 'LineWidth', 1.5);
plot(theta_axis, InputTorque_dynamic, 'r--', 'LineWidth', 1.5);

% Annotate peak magnitude of both static and dynamic torque
annotate_max(theta_axis, InputTorque_static, [0 0.3 0.8], 'N\cdotm');
annotate_max(theta_axis, InputTorque_dynamic, [0.85 0 0], 'N\cdotm');

xlabel('\theta (deg)'); ylabel('Torque (N\cdotm)');
title('Required Input Torque'); 
legend('Static','Dynamic','Location','best');

% -------------------------------------------------------------------
% Plot 4: Center of Mass Kinematics
% -------------------------------------------------------------------
comLabels = {'S_1 (Link AB)', 'S_2 (Link BC)', 'S_3 (Link DE)', 'S_4 (Link EF)', 'S_5 (Link FG)'};
comColors = lines(num_links);

figure('Color', 'w', 'Name', 'Center of Mass Kinematics');

% Velocity Magnitude Subplot
subplot(2, 1, 1);
hold on; grid on;
for i = 1:num_links
    v_mag = hypot(squeeze(vel_COM(1, i, :)), squeeze(vel_COM(2, i, :)));
    plot(theta_axis, v_mag, 'LineWidth', 1.4, 'Color', comColors(i, :), 'DisplayName', comLabels{i});
    annotate_max(theta_axis, v_mag, comColors(i, :), 'm/s');
end
xlabel('\theta (deg)'); ylabel('|v| (m/s)');
title('Center of Mass Linear Velocity Magnitudes');
legend('Location', 'bestoutside');

% Acceleration Magnitude Subplot
subplot(2, 1, 2);
hold on; grid on;
for i = 1:num_links
    a_mag = hypot(squeeze(acc_COM(1, i, :)), squeeze(acc_COM(2, i, :)));
    plot(theta_axis, a_mag, 'LineWidth', 1.4, 'Color', comColors(i, :), 'DisplayName', comLabels{i});
    annotate_max(theta_axis, a_mag, comColors(i, :), 'm/s^2');
end
xlabel('\theta (deg)'); ylabel('|a| (m/s^2)');
title('Center of Mass Linear Acceleration Magnitudes');
legend('Location', 'bestoutside');

% -------------------------------------------------------------------
% Plot 5: Angular Kinematics
% -------------------------------------------------------------------
linkNames = {'Link AB', 'Link BC', 'Link DCE', 'Link EF', 'Link FG'};
linkColors = lines(num_links);

figure('Color', 'w', 'Name', 'Link Angular Kinematics');

% Angular Velocity Subplot
subplot(2, 1, 1);
hold on; grid on;
for i = 1:num_links
    w_mag = abs(squeeze(angularVel(1, i, :)));
    plot(theta_axis, w_mag, 'LineWidth', 1.4, 'Color', linkColors(i, :), 'DisplayName', linkNames{i});
    annotate_max(theta_axis, w_mag, linkColors(i, :), 'rad/s');
end
xlabel('\theta (deg)'); ylabel('|\omega| (rad/s)');
title('Link Angular Velocity Magnitudes');
legend('Location', 'bestoutside');

% Angular Acceleration Subplot
subplot(2, 1, 2);
hold on; grid on;
for i = 1:num_links
    alpha_mag = abs(squeeze(angularAcc(1, i, :)));
    plot(theta_axis, alpha_mag, 'LineWidth', 1.4, 'Color', linkColors(i, :), 'DisplayName', linkNames{i});
    annotate_max(theta_axis, alpha_mag, linkColors(i, :), 'rad/s^2');
end
xlabel('\theta (deg)'); ylabel('|\alpha| (rad/s^2)');
title('Link Angular Acceleration Magnitudes');
legend('Location', 'bestoutside');

%% ------------------------------------------------------------------
%% Tables: Initial State (theta = 0 deg) & Cycle Summary
%% ------------------------------------------------------------------

% -------------------------------------------------------------------
% 1. Table: Initial Joint Kinematics & Reaction Forces
% -------------------------------------------------------------------
joint_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"];

% Extract X, Y positions at step 1 (theta = 0 deg)
X_pos = round(squeeze(positions(1, :, 1))', 4);
Y_pos = round(squeeze(positions(2, :, 1))', 4);

% Joint forces (A through G from solver; H carries external load)
F_stat_mag = [round(hypot(Forces(1, :, 1), Forces(2, :, 1))', 2); norm(Fapp(1:2))];
F_dyn_mag  = [round(hypot(DynamicSol(1, :, 1), DynamicSol(2, :, 1))', 2); norm(Fapp(1:2))];

T_Joints = table(joint_names, X_pos, Y_pos, F_stat_mag, F_dyn_mag, ...
    'VariableNames', {'Joint', 'X_m', 'Y_m', 'Static_Force_N', 'Dynamic_Force_N'});

% -------------------------------------------------------------------
% 2. Table: Link Properties & Center-of-Mass Kinematics
% -------------------------------------------------------------------
link_names = ["AB"; "BC"; "DCE"; "EF"; "FGH"];

m_list     = [masses.AB; masses.BC; masses.DE; masses.EF; masses.FG];
IG_list    = [IG.AB; IG.BC; IG.DE; IG.EF; IG.FG];
w_init     = round(squeeze(angularVel(1, :, 1))', 3);
a_init     = round(squeeze(angularAcc(1, :, 1))', 3);
v_com_init = round(hypot(squeeze(vel_COM(1, :, 1)), squeeze(vel_COM(2, :, 1)))', 3);
a_com_init = round(hypot(squeeze(acc_COM(1, :, 1)), squeeze(acc_COM(2, :, 1)))', 3);

T_Links = table(link_names, m_list, round(IG_list, 4), w_init, a_init, v_com_init, a_com_init, ...
    'VariableNames', {'Link', 'Mass_kg', 'Inertia_kg_m2', 'Omega_rad_s', 'Alpha_rad_s2', 'Vel_COM_m_s', 'Acc_COM_m_s2'});

% -------------------------------------------------------------------
% 3. Table: Executive Summary (Peak vs Initial Loading)
% -------------------------------------------------------------------
Metric = [
    "Input Torque - Initial Static (N*m)";
    "Input Torque - Initial Dynamic (N*m)";
    "Input Torque - Peak Absolute Dynamic (N*m)";
    "Max Reaction Force - Joint A (N)";
    "Max Reaction Force - Joint D (N)";
    "Max Reaction Force - Joint G (N)"
];

Value = [
    round(InputTorque_static(1), 2);
    round(InputTorque_dynamic(1), 2);
    round(max(abs(InputTorque_dynamic)), 2);
    round(max(hypot(DynamicSol(1, 1, :), DynamicSol(2, 1, :))), 2);
    round(max(hypot(DynamicSol(1, 4, :), DynamicSol(2, 4, :))), 2);
    round(max(hypot(DynamicSol(1, 7, :), DynamicSol(2, 7, :))), 2)
];

T_Summary = table(Metric, Value);

% Display in Command Window
disp('--- Joint Positions & Forces at Initial State (theta = 0 deg) ---');
disp(T_Joints);
disp('--- Link Kinematics & Inertial Properties (theta = 0 deg) ---');
disp(T_Links);
disp('--- Mechanism Cycle Extrema Summary ---');
disp(T_Summary);

%% ------------------------------------------------------------------
%% Render & Export High-Res Publication Figures
%% ------------------------------------------------------------------
export_table_image(T_Joints,  'Table_Joint_Data.png',    [620, 240]);
export_table_image(T_Links,   'Table_Link_Data.png',     [740, 180]);
export_table_image(T_Summary, 'Table_Cycle_Summary.png', [500, 200]);

% Optional: Export CSV files
writetable(T_Joints,  'Table_Joint_Data.csv');
writetable(T_Links,   'Table_Link_Data.csv');
writetable(T_Summary, 'Table_Cycle_Summary.csv');

%% ------------------------------------------------------------------
%% 5. Numerical Force Assembler Function
%% ------------------------------------------------------------------
function [forces, Tin] = solve_linkage_forces(pts, S, H, W, Fapp, m, IG, accCOM, a_v, isDynamic)
    % Unknown state vector x (15x1):
    % [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]'
    A_eq = zeros(15, 15);
    b_eq = zeros(15, 1);

    % Column index pointers
    iA = 1:2;  iB = 3:4;  iC = 5:6;  iD = 7:8; 
    iE = 9:10; iF = 11:12; iG = 13:14; iTin = 15;

    % Cross-product moment helper: r x F = rx*Fy - ry*Fx = [-ry, rx] * [Fx; Fy]
    x_pos = @(r) [-r(2),  r(1)];
    x_neg = @(r) [ r(2), -r(1)];

    % 1. Link AB
    A_eq(1:2, iA) = eye(2);
    A_eq(1:2, iB) = eye(2);
    b_eq(1:2)     = -W.AB(1:2)' + (isDynamic * m.AB * accCOM.S1(1:2)');

    A_eq(3, iA)   = x_pos(pts.A - S.S1);
    A_eq(3, iB)   = x_pos(pts.B - S.S1);
    A_eq(3, iTin) = 1;
    b_eq(3)       = isDynamic * IG.AB * a_v.AB(3);

    % 2. Link BC
    A_eq(4:5, iB) = -eye(2);
    A_eq(4:5, iC) = eye(2);
    b_eq(4:5)     = -W.BC(1:2)' + (isDynamic * m.BC * accCOM.S2(1:2)');

    A_eq(6, iB)   = x_neg(pts.B - S.S2);
    A_eq(6, iC)   = x_pos(pts.C - S.S2);
    b_eq(6)       = isDynamic * IG.BC * a_v.BC(3);

    % 3. Link DCE
    A_eq(7:8, iC) = -eye(2);
    A_eq(7:8, iD) = eye(2);
    A_eq(7:8, iE) = eye(2);
    b_eq(7:8)     = -W.DE(1:2)' + (isDynamic * m.DE * accCOM.S3(1:2)');

    A_eq(9, iC)   = x_neg(pts.C - S.S3);
    A_eq(9, iD)   = x_pos(pts.D - S.S3);
    A_eq(9, iE)   = x_pos(pts.E - S.S3);
    b_eq(9)       = isDynamic * IG.DE * a_v.DE(3);

    % 4. Link EF
    A_eq(10:11, iE) = -eye(2);
    A_eq(10:11, iF) = eye(2);
    b_eq(10:11)     = -W.EF(1:2)' + (isDynamic * m.EF * accCOM.S4(1:2)');

    A_eq(12, iE)    = x_neg(pts.E - S.S4);
    A_eq(12, iF)    = x_pos(pts.F - S.S4);
    b_eq(12)        = isDynamic * IG.EF * a_v.EF(3);

    % 5. Link FG
    A_eq(13:14, iF) = -eye(2);
    A_eq(13:14, iG) = eye(2);
    b_eq(13:14)     = -W.FG(1:2)' - Fapp(1:2)' + (isDynamic * m.FG * accCOM.S5(1:2)');

    A_eq(15, iF)    = x_neg(pts.F - S.S5);
    A_eq(15, iG)    = x_pos(pts.G - S.S5);
    M_app_z         = (H(1) - S.S5(1))*Fapp(2) - (H(2) - S.S5(2))*Fapp(1);
    b_eq(15)        = -M_app_z + (isDynamic * IG.FG * a_v.FG(3));

    % Direct linear solve
    x = A_eq \ b_eq;
    forces = reshape(x(1:14), 2, 7);
    Tin = x(15);
end

function annotate_max(x_data, y_data, marker_color, unit_str)
% Find peak value and corresponding angle
[max_val, max_idx] = max(y_data);
x_peak = x_data(max_idx);

% Plot peak marker (excluded from legend)
plot(x_peak, max_val, 'o', ...
    'Color', marker_color, ...
    'MarkerFaceColor', marker_color, ...
    'MarkerSize', 5, ...
    'HandleVisibility', 'off');

% Print callout text above the peak
label_str = sprintf('  %.1f %s (%d°)', max_val, unit_str, round(x_peak));
text(x_peak, max_val, label_str, ...
    'Color', marker_color, ...
    'FontSize', 7.5, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'HandleVisibility', 'off');

% Expand Y-limits slightly so the label doesn't clip against the axis border
yl = ylim;
if yl(2) <= max_val * 1.05
    ylim([yl(1), yl(1) + 1.25 * (max_val - yl(1))]);
end
end

% function export_table_image(tbl, filename, dims)
%     fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, dims(1), dims(2)], 'MenuBar', 'none', 'ToolBar', 'none');
%     uitable(fig, ...
%         'Data', tbl{:,:}, ...
%         'ColumnName', tbl.Properties.VariableNames, ...
%         'RowName', [], ...
%         'Units', 'normalized', ...
%         'Position', [0.02, 0.02, 0.96, 0.96], ...
%         'FontSize', 11);
%     drawnow;
%     exportgraphics(fig, filename, 'Resolution', 300);
%     close(fig);
% end