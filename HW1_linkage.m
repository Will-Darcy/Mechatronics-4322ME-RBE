%Homework 1 six bar linkage
clc;
clear;

A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.2 0];


%Calculating point H at the end of G-F linkage
%distance from F-H
L_FH = 1.83;
% Assume X_G, Y_G, X_F, Y_F are column vectors (N x 1)
dx = F(1) - G(1);
dy = F(2) - G(2);

% Euclidean distances for every step
d_GF = hypot(dx, dy); 

% Unit vector components
ux = dx ./ d_GF;
uy = dy ./ d_GF;
% Coordinates of H for every configuration
X_H = F(1) + L_FH .* ux;
Y_H = F(2) + L_FH .* uy;

H = [X_H Y_H 0];

%joint values
new_B_x = B(1);
new_B_y = B(2);
new_C_x = C(1);
new_C_y = C(2);
new_E_x = E(1);
new_E_y = E(2);
new_F_x = F(1);
new_F_y = F(2);

%Define the lengths of the links
lAB = 0.5726;
lBC = 1.4157;
lDE = 2.4866;
lEF = 1.1754;
lGF = 2.5841;
lCE = norm(C-E);

lCD = norm(D-C); % lengh of part of joint DE

% Weight of each link (Newtons)
WAB = [0, -228.865,  0];
WBC = [0, -553.383,  0];
WDE = [0, -954.701,  0];
WEF = [0, -460.889,  0];
WFG = [0, -1716.634, 0];

% Center of mass of each link
S1 = (A+B)/2;
S2 = (B+C)/2;
S3 = (D+E)/2;
S4 = (E+F)/2;
S5 = (H+G)/2; 

%Applied Force
AppliedForce = [0 -200 0];

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];

%Sum of forces = 0
eqn1 = ForceA + ForceB + WAB ==0;

% Sum of Moments = 0 with respect to COM of Link AB

eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;


% Equations for Link BC
% Sum of Forces = 0
eqn3 = -ForceB + ForceC + WBC ==0;

% Sum of Moments = 0
% Sum of Moments = 0 with respect to CoM of Link BC

eqn4 = cross(B-S2, -ForceB) + cross(C-S2,ForceC) == 0;


% Equations for Link CD
% Sum of Forces = 0 for link DCE
eqn5 = -ForceC + ForceD + ForceE + WDE == 0;

% Sum of Moments = 0
% Sum of Moments = 0 with respect to CoM of Link CD
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;


% Equations for Link EF

% Sum of Forces = 0
eqn7 = -ForceE + ForceF + WEF == 0;

% Sum of Moments = 0
% S4E x -FE + S4F x FF = 0;
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;


% Equations for Link FG

% Sum of Forces = 0
eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

% Sum of Moments = 0
% Sum of Moments = 0 with respect to CoM of Link FG
eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) + cross(H-S5, AppliedForce) == 0;


% Solve the 10 equations

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution = solve(eqnMatrix,[FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);


Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);
Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);
Input_Torque = double(StaticSolution.Tin);

%Display Calculated Results

syms wBC wDE
angVel_AB=[0 0 1];
omega_BC=[0 0 wBC];
omega_DE = [0 0 wDE];

eqn11 = cross(angVel_AB, B-A) + cross(omega_BC,C-B) + cross(omega_DE,D-C) == 0;

loop1Solution = solve(eqn11,[wBC wDE]);

angularVelocity_BC = double(loop1Solution.wBC);
angularVelocity_DE = double(loop1Solution.wDE);


angVel_BC=[0 0 angularVelocity_BC];
angVel_DE=[0 0 angularVelocity_DE];

syms wEF wFG

omega_EF=[0 0 wEF];
omega_FG=[0 0 wFG];

eqn12 = cross(angVel_DE,E-D) + cross(omega_EF,F-E) + cross(omega_FG,G-F) == 0;

loop2Solution = solve(eqn12,[wEF wFG]);

%extract solutions
angularVelocity_EF = double(loop2Solution.wEF);
angularVelocity_FG = double(loop2Solution.wFG);

angVel_EF = [0 0 angularVelocity_EF];
angVel_FG = [0 0 angularVelocity_FG];

%acceleration loop equations
syms aBC aDE

alpha_AB = [0 0 0]; % assuming no acceleration
alpha_BC = [0 0 aBC];
alpha_DE = [0 0 aDE];

a_A_B = cross(alpha_AB,B-A) + cross(angVel_AB,cross(angVel_AB,B-A));
a_C_B = cross(alpha_BC,C-B) + cross(angVel_BC,cross(angVel_BC,C-B));
a_D_E = cross(alpha_DE,D-C) + cross(angVel_DE,cross(angVel_DE,D-C));

eqn13 = a_A_B + a_C_B + a_D_E == 0;

loop1AccSolution = solve(eqn13,[aBC aDE]);

alphaBC = double(loop1AccSolution.aBC);
alphaDE = double(loop1AccSolution.aDE);


acc_BC = [0 0 alphaBC];
acc_DE = [0 0 alphaDE];

syms aEF aFG

alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];


a_D_E = cross(acc_DE,C-D) + cross(angVel_DE,cross(angVel_DE,C-D));
a_E_F = cross(alpha_EF,F-E) + cross(angVel_EF, cross(angVel_EF,F-E));
a_F_G = cross(alpha_FG,G-F) + cross(angVel_FG, cross(angVel_FG,G-F));

eqn14 = a_D_E + a_E_F + a_F_G == 0;

loop2AccSolution = solve(eqn14,[aEF aFG]);

alphaEF = double(loop2AccSolution.aEF);
alphaFG = double(loop2AccSolution.aFG);

acc_EF = [0 0 alphaEF];
acc_FG = [0 0 alphaFG];


%velocities of each joint relative to ground
vB_GR = cross(angVel_AB,B-A);

vC_GR = cross(angVel_DE,C-D);

vE_GR = cross(angVel_DE,E-D);

vF_GR = cross(angVel_FG,F-G);

%Calulating velocities for center of masses for each link

vS1_GR = cross(angVel_AB,S1-A); %Link AB

vS2_B = cross(angVel_BC,S2-B); 
vS2_GR = vS2_B + vB_GR; % Link BC

vS3_GR = cross(angVel_DE,S3-D); % Link DE

vS4_F = cross(angVel_EF,S4-F);
vS4_GR = vS4_F + vF_GR; % Link EF

vS5_GR = cross(angVel_FG,S5-G); % Link FGH


% Calculating the Accelerations for the center of masses
%Link AB
aS1_GR = cross(alpha_AB,S1-A) + cross(angVel_AB,cross(angVel_AB,S1-A));

aS2_B = cross(acc_BC,S2-B) + cross(angVel_BC,cross(angVel_BC,S2-B));
aB_A = cross(alpha_AB,B-A) + cross(angVel_AB,cross(angVel_AB,B-A));
%Link BC
aS2_GR = aS2_B + aB_A;

%Link BE
aS3_GR = cross(acc_DE,S3-D) + cross(angVel_DE,cross(angVel_DE,S3-D));

aS4_F = cross(acc_EF,S4-F) + cross(angVel_EF,cross(angVel_EF,S4-F));
aF_G = cross(acc_FG,F-G) + cross(angVel_FG,cross(angVel_FG,F-G));
%Link EF
aS4_GR = aS4_F + aF_G;
%Link FGH
aS5_GR = cross(acc_FG,S5-G) + cross(angVel_FG,cross(angVel_FG,S5-G));

%Masses in Kilograms
MassAB = 23.3377;
MassBC = 56.42938;
MassDE = 97.35243;
MassEF = 46.9976;
MassFG = 175.04798;

function [IG] = calc_link_inertia(firstpt, secondpt, mass)
% CALC_LINK_INERTIA Calculates planar moment of inertia for a rectangular link.
%
% Inputs:
%   firstpt       - [X, Y] start position
%   secondpt       - [X, Y] intermediate position
%   mass    - Link mass (kg)
%   
% Constants:
%   width   - Constant cross-sectional width (m)
%
% Outputs:
%   IG      - Mass moment of inertia about COM (kg*m^2)

% 1. Compute total length along the link line
d_GF = norm(secondpt - firstpt);

% 2. Planar rectangular inertia about COM
IG = (1/12) * mass * (d_GF^2 + 0.1^2);
end

%mass moment of inertia
J_AB = calc_link_inertia(A, B, MassAB);
J_BC = calc_link_inertia(B, C, MassBC);
J_DE = calc_link_inertia(D, E, MassDE);
J_EF = calc_link_inertia(E, F, MassEF);
J_FG = calc_link_inertia(H, G, MassFG);


syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

%define fores
NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NImputT = [0 0 NTin];

eqn15 = NForceA + NForceB + WAB == MassAB * aS1_GR;

eqn16 = cross(A-S1,NForceA) + cross(B-S1,NForceB) + NImputT == J_AB * alpha_AB;

eqn17 = -NForceB + NForceC + WBC == MassBC * aS2_GR;

eqn18 = cross(B-S2,-NForceB) + cross(C-S2,NForceC) == J_BC * acc_BC;

eqn19 = NForceE -NForceC + NForceD + WDE == MassDE * aS3_GR;

eqn20 = cross(E-S3,NForceE) + cross(C-S3, -NForceC) + cross(D-S3,NForceD) == J_DE * acc_DE;

eqn21 = -NForceE + NForceF + WEF == MassEF * aS4_GR;

eqn22 = cross(E-S4, -NForceE) + cross(F-S4,NForceF) == J_EF * acc_EF;

eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG * aS5_GR;

eqn24 = cross(F-S5, -NForceF) + cross(G-S5,NForceG) + cross(H-S5,AppliedForce) == J_FG * acc_FG;

%Solving equations
NeqMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
DynamicSolution = solve(NeqMatrix, [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin]);

% Extract forces from the dynamic solution
NForceAx = double(DynamicSolution.NFAx);
NForceAy = double(DynamicSolution.NFAy);
NForceBx = double(DynamicSolution.NFBx);
NForceBy = double(DynamicSolution.NFBy);
NForceCx = double(DynamicSolution.NFCx);
NForceCy = double(DynamicSolution.NFCy);
NForceDx = double(DynamicSolution.NFDx);
NForceDy = double(DynamicSolution.NFDy);
NForceEx = double(DynamicSolution.NFEx);
NForceEy = double(DynamicSolution.NFEy);
NForceFx = double(DynamicSolution.NFFx);
NForceFy = double(DynamicSolution.NFFy);
NForceGx = double(DynamicSolution.NFGx);
NForceGy = double(DynamicSolution.NFGy);
NImput_T = double(DynamicSolution.NTin);


%prealocation of memory for forces to be stored
%calulation for all the forces at each link durning a full 360 degree turn
% Dimensions: [coordinate (X=1, Y=2), joint (1:7), angle (1:360)]
%storage positions
% Joint column key used everywhere below: 1=A 2=B 3=C 4=D 5=E 6=F 7=G
num_joints = 7; %H is going to be stored he even though it isn't a joint
num_angles = 360;
positions = NaN(2, num_joints+1, num_angles);

% A, D, G are ground pivots and never move, so fill their columns (1,4,7)
% once, across every page (angle), right now.
positions(:, 1, :) = repmat(A(1:2)', 1, 1, num_angles); % A
positions(:, 4, :) = repmat(D(1:2)', 1, 1, num_angles); % D
positions(:, 7, :) = repmat(G(1:2)', 1, 1, num_angles); % G

%Will create matricies for other aspects here
%create matricies for all other aspects to solve for forces, velocities,
%acc, and dynamic solutions
Forces = NaN(2, num_joints, num_angles);

% Angular velocity/acceleration only exist per LINK (AB,BC,DE,EF,FG),
% not per joint, so they get their own smaller matrices instead of
% reusing the 7-joint size. Link column key: 1=AB 2=BEC 3=CD 4=EF 5=FG
num_links = 5;
angularVel = NaN(1, num_links, num_angles);
angularAcc = NaN(1, num_links, num_angles);
DynamicSol = NaN(2, num_joints, num_angles);



%compute intial angle of input link AB
initial_theta = atan2(B(2)-A(2),B(1)-A(1));
if (initial_theta<0)
    NewAngle = 2*pi + initial_theta;
else
    NewAngle = initial_theta;
end

for theta=1:1:360
    %starting with finding new position of B
    B_new = A + [lAB*cos(NewAngle+deg2rad(theta)) lAB*sin(NewAngle+deg2rad(theta)) 0];

    %new position of C
    %with B_new as center, BC as radius
    %with D as center and DC as radius
    [Cx, Cy] = circcirc(B_new(1),B_new(2),lBC,D(1),D(2),lCD);

    %checking if there is a NaN
    circIntersect_x_C = any(isnan(vpa(Cx)));
    circIntersect_y_C = any(isnan(vpa(Cy)));
    if circIntersect_x_C==0 && circIntersect_y_C==0
        C_1 = [Cx(1) Cy(1) 0];
        C_2 = [Cx(2) Cy(2) 0];

        %distance to determine weather C-1 or C_2 is correct

        distC1 = norm(C_1-C);
        distC2 = norm(C_2-C);

        if(distC1<distC2)
            C_new = vpa(C_1);

        else

            C_new = vpa(C_2);
        end

        % %caclulating new angle of D to find E
        % D_theta = atan2(C_new(2)-D(2),C_new(1)-D(1));
        % if (D_theta<0)
        %     D_Theta_New = 2*pi + D_theta;
        % else
        %     D_Theta_New = D_theta;
        % end
        % 
        % % Calulating new E
        % E_new = D + [lDE*cos(D_Theta_New) lDE*sin(D_Theta_New) 0];

        %with D as center and DC as radius
        [Ex, Ey] = circcirc(C_new(1),C_new(2),lCE,D(1),D(2),lDE);

        %checking if there is a NaN
        circIntersect_x_E = any(isnan(vpa(Ex)));
        circIntersect_y_E = any(isnan(vpa(Ey)));
        if circIntersect_x_E==0 && circIntersect_y_E==0
            E_1 = [Ex(1) Ey(1) 0];
            E_2 = [Ex(2) Ey(2) 0];

            %distance to determine weather C-1 or C_2 is correct

            distC1 = norm(E_1-E);
            distC2 = norm(E_2-E);

            if(distC1<distC2)
                E_new = vpa(E_1);

            else

                E_new = vpa(E_2);
            end

             %New position of Joint F using E_new and G
        % New position of Joint F using E_new and G
        [Fx, Fy] = circcirc(E_new(1), E_new(2), lEF, G(1), G(2), lGF);

        % Checking if there is a NaN
        circIntersect_x_F = any(isnan(vpa(Fx)));
        circIntersect_y_F = any(isnan(vpa(Fy)));

        if circIntersect_x_F == 0 && circIntersect_y_F == 0
            F_1 = [Fx(1) Fy(1) 0]; 
            F_2 = [Fx(2) Fy(2) 0];

            %distance to determine weather F_1 or F_2 is correct
            distF1 = norm(F_1-F);
            distF2 = norm(F_2-F);

            if (distF1 < distF2)
                F_new = vpa(F_1);
            else
                F_new = vpa(F_2);
            end

            % Assume X_G, Y_G, X_F, Y_F are column vectors (N x 1)
            dx = F_new(1) - G(1);
            dy = F_new(2) - G(2);

            d_GF = hypot(dx, dy); 

            % Unit vector components
            ux = dx ./ d_GF;
            uy = dy ./ d_GF;
            % Coordinates of H for configuration
            X_H = F(1) + L_FH .* ux;
            Y_H = F(2) + L_FH .* uy;

            H_new = [X_H Y_H 0];



            %store values in matrix
            % positions(:, jointColumn, theta) = [x; y] for this angle.
            % double() is needed because C_new/E_new/F_new are vpa
            % (symbolic) values, and positions is a plain double array.
            positions(:, 2, theta) = double(B_new(1:2))'; % B
            positions(:, 3, theta) = double(C_new(1:2))'; % C
            positions(:, 5, theta) = double(E_new(1:2))'; % E
            positions(:, 6, theta) = double(F_new(1:2))'; % F
            positions(:, 8, theta) = double(H_new(1:2))'; % H

            % Calculating forces for new positions
            % Center of mass of each link
            S1_new = (A+B_new)/2;
            S2_new = (B_new+C_new)/2;
            S3_new = (D+E_new)/2;
            S4_new = (E_new+F_new)/2;
            S5_new = (H_new+G)/2; 


            syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

            ForceA = [FAx FAy 0];
            ForceB = [FBx FBy 0];
            ForceC = [FCx FCy 0];
            ForceD = [FDx FDy 0];
            ForceE = [FEx FEy 0];
            ForceF = [FFx FFy 0];
            ForceG = [FGx FGy 0];
            InputTorque = [0 0 Tin];

            %Sum of forces = 0
            eqn1 = ForceA + ForceB + WAB ==0;

            % Sum of Moments = 0 with respect to COM of Link AB

            eqn2 = cross(A-S1_new,ForceA) + cross(B_new-S1_new,ForceB) + InputTorque == 0;


            % Equations for Link BC
            % Sum of Forces = 0
            eqn3 = -ForceB + ForceC + WBC ==0;

            % Sum of Moments = 0
            % Sum of Moments = 0 with respect to CoM of Link BC

            eqn4 = cross(B_new-S2_new, -ForceB) + cross(C_new-S2_new,ForceC) == 0;


            % Equations for Link CD
            % Sum of Forces = 0 for link DCE
            eqn5 = -ForceC + ForceD + ForceE + WDE == 0;

            % Sum of Moments = 0
            % Sum of Moments = 0 with respect to CoM of Link CD
            eqn6 = cross(C_new-S3_new, -ForceC) + cross(D-S3_new, ForceD) == 0;


            % Equations for Link EF

            % Sum of Forces = 0
            eqn7 = -ForceE + ForceF + WEF == 0;

            % Sum of Moments = 0
            % S4E x -FE + S4F x FF = 0;
            eqn8 = cross(E_new-S4_new, -ForceE) + cross(F_new-S4_new, ForceF) == 0;


            % Equations for Link FG

            % Sum of Forces = 0
            eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

            % Sum of Moments = 0
            % Sum of Moments = 0 with respect to CoM of Link FG
            eqn10 = cross(F_new-S5_new, -ForceF) + cross(G-S5_new, ForceG) + cross(H_new-S5_new, AppliedForce) == 0;


            % Solve the 10 equations

            eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];
            StaticSolution = solve(eqnMatrix,[FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);


            Force_Ax_new = double(StaticSolution.FAx);
            Force_Ay_new = double(StaticSolution.FAy);
            Force_Bx_new = double(StaticSolution.FBx);
            Force_By_new = double(StaticSolution.FBy);
            Force_Cx_new = double(StaticSolution.FCx);
            Force_Cy_new = double(StaticSolution.FCy);
            Force_Dx_new = double(StaticSolution.FDx);
            Force_Dy_new = double(StaticSolution.FDy);
            Force_Ex_new = double(StaticSolution.FEx);
            Force_Ey_new = double(StaticSolution.FEy);
            Force_Fx_new = double(StaticSolution.FFx);
            Force_Fy_new = double(StaticSolution.FFy);
            Force_Gx_new = double(StaticSolution.FGx);
            Force_Gy_new = double(StaticSolution.FGy);
            Input_Torque_new = double(StaticSolution.Tin);
            
            %store static joint forces in the Forces matrix for this angle

            Forces(:, 1, theta) = [Force_Ax_new; Force_Ay_new]; % A
            Forces(:, 2, theta) = [Force_Bx_new; Force_By_new]; % B
            Forces(:, 3, theta) = [Force_Cx_new; Force_Cy_new]; % C
            Forces(:, 4, theta) = [Force_Dx_new; Force_Dy_new]; % D
            Forces(:, 5, theta) = [Force_Ex_new; Force_Ey_new]; % E
            Forces(:, 6, theta) = [Force_Fx_new; Force_Fy_new]; % F
            Forces(:, 7, theta) = [Force_Gx_new; Force_Gy_new]; % G


            %Solving for new Angel velocities
            syms wBC wDE
            angVel_AB_new=[0 0 1];
            omega_BC=[0 0 wBC];
            omega_DE = [0 0 wDE];

            eqn11 = cross(angVel_AB, B_new-A) + cross(omega_BC,C_new-B_new) + cross(omega_DE,D-C_new) == 0;

            loop1Solution = solve(eqn11,[wBC wDE]);

            angularVelocity_BC = double(loop1Solution.wBC);
            angularVelocity_DE = double(loop1Solution.wDE);


            angVel_BC_new=[0 0 angularVelocity_BC];
            angVel_DE_new=[0 0 angularVelocity_DE];

            syms wEF wFG

            omega_EF=[0 0 wEF];
            omega_FG=[0 0 wFG];

            eqn12 = cross(angVel_DE_new,E_new-D) + cross(omega_EF,F_new-E_new) + cross(omega_FG,G-F_new) == 0;

            loop2Solution = solve(eqn12,[wEF wFG]);

            %extract solutions
            angularVelocity_EF = double(loop2Solution.wEF);
            angularVelocity_FG = double(loop2Solution.wFG);

            angVel_EF_new = [0 0 angularVelocity_EF];
            angVel_FG_new = [0 0 angularVelocity_FG];

            %store angular velocities for this angle (z-component only)
            angularVel(1, 1, theta) = 1;
            angularVel(1, 2, theta) = angVel_BC_new(3);     % BC
            angularVel(1, 3, theta) = angVel_DE_new(3);     % DE
            angularVel(1, 4, theta) = angVel_EF_new(3);     % EF
            angularVel(1, 5, theta) = angVel_FG_new(3);     % FG
            
            %acceleration loop equations
            syms aBC aDE

            alpha_AB_new = [0 0 0]; % assuming no acceleration
            alpha_BC = [0 0 aBC];
            alpha_DE = [0 0 aDE];

            a_A_B = cross(alpha_AB_new,B_new-A) + cross(angVel_AB,cross(angVel_AB,B_new-A));
            a_C_B = cross(alpha_BC,C_new-B_new) + cross(angVel_BC,cross(angVel_BC,C_new-B_new));
            a_D_E = cross(alpha_DE,D-C_new) + cross(angVel_DE,cross(angVel_DE,D-C_new));

            eqn13 = a_A_B + a_C_B + a_D_E == 0;

            loop1AccSolution = solve(eqn13,[aBC aDE]);

            alphaBC = double(loop1AccSolution.aBC);
            alphaDE = double(loop1AccSolution.aDE);


            acc_BC_new = [0 0 alphaBC];
            acc_DE_new = [0 0 alphaDE];

            syms aEF aFG

            alpha_EF = [0 0 aEF];
            alpha_FG = [0 0 aFG];


            a_D_E = cross(acc_DE,C_new-D) + cross(angVel_DE,cross(angVel_DE,C_new-D));
            a_E_F = cross(alpha_EF,F_new-E_new) + cross(angVel_EF, cross(angVel_EF,F_new-E_new));
            a_F_G = cross(alpha_FG,G-F_new) + cross(angVel_FG, cross(angVel_FG,G-F_new));

            eqn14 = a_D_E + a_E_F + a_F_G == 0;

            loop2AccSolution = solve(eqn14,[aEF aFG]);

            alphaEF = double(loop2AccSolution.aEF);
            alphaFG = double(loop2AccSolution.aFG);

            acc_EF_new = [0 0 alphaEF];
            acc_FG_new = [0 0 alphaFG];

            %store angular accelerations for this angle (z-component only)
            angularAcc(1, 1, theta) = 0;                % AB (assumed 0)
            angularAcc(1, 2, theta) = acc_BC_new(3);         % BC
            angularAcc(1, 3, theta) = acc_DE_new(3);         % DE
            angularAcc(1, 4, theta) = acc_EF_new(3);         % EF
            angularAcc(1, 5, theta) = acc_FG_new(3);         % FG

            %Calulating new velocities and acc of center of masses to GR

            %velocities of each joint relative to ground
            vB_GR_new = cross(angVel_AB_new,B_new-A);

            vC_GR_new = cross(angVel_DE_new,C_new-D);

            vE_GR_new = cross(angVel_DE_new,E_new-D);

            vF_GR_new = cross(angVel_FG_new,F_new-G);

            %Calulating velocities for center of masses for each link

            vS1_GR_new = cross(angVel_AB_new,S1_new-A); %Link AB

            vS2_B_new = cross(angVel_BC_new,S2_new-B); 
            vS2_GR_new = vS2_B_new + vB_GR_new; % Link BC

            vS3_GR_new = cross(angVel_DE_new,S3_new-D); % Link DE

            vS4_F_new = cross(angVel_EF_new,S4_new-F_new);
            vS4_GR_new = vS4_F_new + vF_GR_new; % Link EF

            vS5_GR_new = cross(angVel_FG_new,S5_new-G); % Link FGH


            % Calculating the Accelerations for the center of masses
            %Link AB
            aS1_GR_new = cross(alpha_AB_new,S1_new-A) + cross(angVel_AB_new,cross(angVel_AB_new,S1_new-A));

            aS2_B_new = cross(acc_BC_new,S2_new-B_new) + cross(angVel_BC_new,cross(angVel_BC_new,S2_new-B_new));
            aB_A_new = cross(alpha_AB_new,B_new-A) + cross(angVel_AB_new,cross(angVel_AB_new,B_new-A));
            %Link BC
            aS2_GR_new = aS2_B_new + aB_A_new;

            %Link BE
            aS3_GR_new = cross(acc_DE_new,S3_new-D) + cross(angVel_DE_new,cross(angVel_DE_new,S3_new-D));

            aS4_F_new = cross(acc_EF_new,S4_new-F_new) + cross(angVel_EF_new,cross(angVel_EF_new,S4_new-F_new));
            aF_G_new = cross(acc_FG_new,F_new-G) + cross(angVel_FG_new,cross(angVel_FG_new,F_new-G));
            %Link EF
            aS4_GR_new = aS4_F_new + aF_G_new;
            %Link FGH
            aS5_GR_new = cross(acc_FG_new,S5_new-G) + cross(angVel_FG_new,cross(angVel_FG_new,S5_new-G));

            %Calculating the applied forces in loop
            syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

            %define fores
            NForceA_new = [NFAx NFAy 0];
            NForceB_new = [NFBx NFBy 0];
            NForceC_new = [NFCx NFCy 0];
            NForceD_new = [NFDx NFDy 0];
            NForceE_new = [NFEx NFEy 0];
            NForceF_new = [NFFx NFFy 0];
            NForceG_new = [NFGx NFGy 0];
            NImputT_new = [0 0 NTin];

            eqn15 = NForceA_new + NForceB_new + WAB == MassAB * aS1_GR_new;

            eqn16 = cross(A-S1_new,NForceA_new) + cross(B_new-S1_new,NForceB_new) + NImputT_new == J_AB * alpha_AB_new;

            eqn17 = -NForceB_new + NForceC_new + WBC == MassBC * aS2_GR;

            eqn18 = cross(B_new-S2_new,-NForceB_new) + cross(C_new-S2_new,NForceC_new) == J_BC * acc_BC_new;

            eqn19 = NForceE_new -NForceC_new + NForceD_new + WDE == MassDE * aS3_GR_new;

            eqn20 = cross(E_new-S3_new,NForceE_new) + cross(C_new-S3_new, -NForceC_new) + cross(D-S3_new,NForceD_new) == J_DE * acc_DE_new;

            eqn21 = -NForceE_new + NForceF_new + WEF == MassEF * aS4_GR_new;

            eqn22 = cross(E_new-S4_new, -NForceE_new) + cross(F_new-S4_new,NForceF_new) == J_EF * acc_EF_new;

            eqn23 = -NForceF_new + NForceG_new + WFG + AppliedForce == MassFG * aS5_GR_new;

            eqn24 = cross(F_new-S5_new, -NForceF_new) + cross(G-S5_new,NForceG_new) + cross(H_new-S5_new,AppliedForce) == J_FG * acc_FG_new;

            %Solving equations
            NeqMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
            DynamicSolution = solve(NeqMatrix, [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin]);

            % Extract forces from the dynamic solution
            NForceAx_new = double(DynamicSolution.NFAx);
            NForceAy_new = double(DynamicSolution.NFAy);
            NForceBx_new = double(DynamicSolution.NFBx);
            NForceBy_new = double(DynamicSolution.NFBy);
            NForceCx_new = double(DynamicSolution.NFCx);
            NForceCy_new = double(DynamicSolution.NFCy);
            NForceDx_new = double(DynamicSolution.NFDx);
            NForceDy_new = double(DynamicSolution.NFDy);
            NForceEx_new = double(DynamicSolution.NFEx);
            NForceEy_new = double(DynamicSolution.NFEy);
            NForceFx_new = double(DynamicSolution.NFFx);
            NForceFy_new = double(DynamicSolution.NFFy);
            NForceGx_new = double(DynamicSolution.NFGx);
            NForceGy_new = double(DynamicSolution.NFGy);
            NImput_T_new = double(DynamicSolution.NTin);

            %store dynamic joint forces in the DynamicSol matrix for this angle
            DynamicSol(:, 1, theta) = [NForceAx_new; NForceAy_new]; % A
            DynamicSol(:, 2, theta) = [NForceBx_new; NForceBy_new]; % B
            DynamicSol(:, 3, theta) = [NForceCx_new; NForceCy_new]; % C
            DynamicSol(:, 4, theta) = [NForceDx_new; NForceDy_new]; % D
            DynamicSol(:, 5, theta) = [NForceEx_new; NForceEy_new]; % E
            DynamicSol(:, 6, theta) = [NForceFx_new; NForceFy_new]; % F
            DynamicSol(:, 7, theta) = [NForceGx_new; NForceGy_new]; % G
            fprintf('New Position at angle: %d degree\n',theta)

            B=B_new;
            C=C_new;
            E=E_new;
            F=F_new;

        else
            fprintf('New Position of F cannot be determined at angle: %d degree\n',theta)
        end


        else
             fprintf('New Position of E cannot be determined at angle: %d degree\n',theta)
        end

    else
        fprintf('New Position of C cannot be determined at angle: %d degree\n',theta)
    end


end

%% ------------------------------------------------------------------
%% Plot 1: joint path traces (coupler curves) for the moving joints
%% ------------------------------------------------------------------
figure('Color','w');
hold on; axis equal; grid on;

jointNames = {'B','C','E','F'};
jointCols  = [2 3 5 6];
markerColors = lines(numel(jointCols));

for k = 1:numel(jointCols)
    col = jointCols(k);
    x = squeeze(positions(1, col, :));
    y = squeeze(positions(2, col, :));
    plot(x, y, '.', 'Color', markerColors(k,:), 'MarkerSize', 6, 'DisplayName', jointNames{k});
end

plot([A(1) D(1) G(1)], [A(2) D(2) G(2)], 'ks', 'MarkerFaceColor','k', 'MarkerSize', 8, 'DisplayName','Ground pivots');
text(A(1), A(2), '  A'); text(D(1), D(2), '  D'); text(G(1), G(2), '  G');

legend('Location','bestoutside');
xlabel('X (m)'); ylabel('Y (m)');
title('Joint Path Traces Over Full Rotation');

%% ====================================================================
%% Plotting: joint forces and input torque over the full 360-deg sweep
%% ====================================================================
% Uses Forces (static, weight-only) and DynamicSol (full dynamic,
% includes inertia) populated in the loop above.
% Joint column key: 1=A 2=B 3=C 4=D 5=E 6=F 7=G
jointLabels = {'A','B','C','D','E','F','G'};
theta_axis = 1:num_angles;

%% Figure 1: force magnitude |F| at each joint, static vs dynamic
figure('Color','w','Name','Joint Force Magnitudes');
for j = 1:num_joints
    subplot(3,3,j);

    Fx_static = squeeze(Forces(1,j,:));
    Fy_static = squeeze(Forces(2,j,:));
    mag_static = hypot(Fx_static, Fy_static);

    Fx_dynamic = squeeze(DynamicSol(1,j,:));
    Fy_dynamic = squeeze(DynamicSol(2,j,:));
    mag_dynamic = hypot(Fx_dynamic, Fy_dynamic);

    plot(theta_axis, mag_static, 'b-', 'LineWidth', 1.3); hold on;
    plot(theta_axis, mag_dynamic, 'r--', 'LineWidth', 1.3);
    grid on;
    title(['Joint ' jointLabels{j}]);
    xlabel('\theta (deg)'); ylabel('|F| (N)');
    if j == 1
        legend('Static','Dynamic','Location','best');
    end
end
sgtitle('Joint Force Magnitude vs Input Angle');

%% Figure 2: Fx at each joint, static vs dynamic
figure('Color','w','Name','Joint Fx Components');
for j = 1:num_joints
    subplot(3,3,j);
    plot(theta_axis, squeeze(Forces(1,j,:)), 'b-', 'LineWidth', 1.3); hold on;
    plot(theta_axis, squeeze(DynamicSol(1,j,:)), 'r--', 'LineWidth', 1.3);
    grid on;
    title(['Joint ' jointLabels{j}]);
    xlabel('\theta (deg)'); ylabel('F_x (N)');
    if j == 1
        legend('Static','Dynamic','Location','best');
    end
end
sgtitle('Joint X-Force vs Input Angle');

%% Figure 3: Fy at each joint, static vs dynamic
figure('Color','w','Name','Joint Fy Components');
for j = 1:num_joints
    subplot(3,3,j);
    plot(theta_axis, squeeze(Forces(2,j,:)), 'b-', 'LineWidth', 1.3); hold on;
    plot(theta_axis, squeeze(DynamicSol(2,j,:)), 'r--', 'LineWidth', 1.3);
    grid on;
    title(['Joint ' jointLabels{j}]);
    xlabel('\theta (deg)'); ylabel('F_y (N)');
    if j == 1
        legend('Static','Dynamic','Location','best');
    end
end
sgtitle('Joint Y-Force vs Input Angle');

%% Figure 4: input torque, static vs dynamic
figure('Color','w','Name','Input Torque');
plot(theta_axis, InputTorque_static, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_axis, InputTorque_dynamic, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\theta (deg)'); ylabel('Input Torque (N\cdotm)');
title('Required Input Torque vs Input Angle');
legend('Static (weight only)','Dynamic (with inertia)','Location','best');