close all;
clear;
clc;

%% 
%Define all Parameters

% Parameters
mu0 = 4*pi*1e-7; % Permeability of free space
m1 = [10, 0, 0];   % Magnetic moment vector of the first dipole (aligned along z-axis)
m2 = [0, 0, 10];   % Magnetic moment vector of the second dipole (same orientation)
m3 = [0, 0, -10];   % Magnetic moment vector of the third dipole (same orientation)
m4 = [10, 0, 0];   % Magnetic moment vector of the fourth dipole (same orientation)
m5 = [10, 0, 0];   % Magnetic moment vector of the fourth dipole (same orientation)

% Position offset for the second dipole (0.15 meters above the first dipole)
p1 = [0, 0, 0];    % Position of the first dipole
p2 = [0.1, 0, 0.15]; % Position of the second dipole (change as needed)
p3 = [0, 0, 0.3]; % Position of the third dipole (change as needed)
p4 = [0, 0.10, 0.45]; % Position of the fourth dipole (change as needed)
p5 = [0, 0, 0.6]; % Position of the fifth dipole (change as needed)

% Torque Scale Factor
TSF = 10;

%%
%Calculate Magentic fields

% Grid
[x, y, z] = meshgrid(linspace(-0.3, 0.3, 25), linspace(-0.3, 0.3, 25), linspace(-0.3, 0.9, 25));

dx = (0.3-(-0.3))/(25-1);
dy = (0.3-(-0.3))/(25-1);
dz = (0.9-(-0.3))/(25-1);

% Calculate field components for the first dipole
x1 = x - p1(1);
y1 = y - p1(2);
z1 = z - p1(3);
r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;

Bx1 = mu0/(4*pi) * (3*(m1(1)*rx1 + m1(2)*ry1 + m1(3)*rz1).*rx1 - m1(1))./r1.^3;
By1 = mu0/(4*pi) * (3*(m1(1)*rx1 + m1(2)*ry1 + m1(3)*rz1).*ry1 - m1(2))./r1.^3;
Bz1 = mu0/(4*pi) * (3*(m1(1)*rx1 + m1(2)*ry1 + m1(3)*rz1).*rz1 - m1(3))./r1.^3;

% Calculate field components for the second dipole
x2 = x - p2(1);
y2 = y - p2(2);
z2 = z - p2(3);
r2 = sqrt(x2.^2 + y2.^2 + z2.^2);
rx2 = x2./r2; ry2 = y2./r2; rz2 = z2./r2;

Bx2 = mu0/(4*pi) * (3*(m2(1)*rx2 + m2(2)*ry2 + m2(3)*rz2).*rx2 - m2(1))./r2.^3;
By2 = mu0/(4*pi) * (3*(m2(1)*rx2 + m2(2)*ry2 + m2(3)*rz2).*ry2 - m2(2))./r2.^3;
Bz2 = mu0/(4*pi) * (3*(m2(1)*rx2 + m2(2)*ry2 + m2(3)*rz2).*rz2 - m2(3))./r2.^3;

% Calculate field components for the thrid dipole
x3 = x - p3(1);
y3 = y - p3(2);
z3 = z - p3(3);
r3 = sqrt(x3.^2 + y3.^2 + z3.^2);
rx3 = x3./r3; ry3 = y3./r3; rz3 = z3./r3;

Bx3 = mu0/(4*pi) * (3*(m3(1)*rx3 + m3(2)*ry3 + m3(3)*rz3).*rx3 - m3(1))./r3.^3;
By3 = mu0/(4*pi) * (3*(m3(1)*rx3 + m3(2)*ry3 + m3(3)*rz3).*ry3 - m3(2))./r3.^3;
Bz3 = mu0/(4*pi) * (3*(m3(1)*rx3 + m3(2)*ry3 + m3(3)*rz3).*rz3 - m3(3))./r3.^3;

% Calculate field components for the Fourth dipole
x4 = x - p4(1);
y4 = y - p4(2);
z4 = z - p4(3);
r4 = sqrt(x4.^2 + y4.^2 + z4.^2);
rx4 = x4./r4; ry4 = y4./r4; rz4 = z4./r4;

Bx4 = mu0/(4*pi) * (3*(m4(1)*rx4 + m4(2)*ry4 + m4(3)*rz4).*rx4 - m4(1))./r4.^3;
By4 = mu0/(4*pi) * (3*(m4(1)*rx4 + m4(2)*ry4 + m4(3)*rz4).*ry4 - m4(2))./r4.^3;
Bz4 = mu0/(4*pi) * (3*(m4(1)*rx4 + m4(2)*ry4 + m4(3)*rz4).*rz4 - m4(3))./r4.^3;

% Calculate field components for the Fourth dipole
x5 = x - p5(1);
y5 = y - p5(2);
z5 = z - p5(3);
r5 = sqrt(x5.^2 + y5.^2 + z5.^2);
rx5 = x5./r5; ry5 = y5./r5; rz5 = z5./r5;

Bx5 = mu0/(4*pi) * (3*(m5(1)*rx5 + m5(2)*ry5 + m5(3)*rz5).*rx5 - m5(1))./r5.^3;
By5 = mu0/(4*pi) * (3*(m5(1)*rx5 + m5(2)*ry5 + m5(3)*rz5).*ry5 - m5(2))./r5.^3;
Bz5 = mu0/(4*pi) * (3*(m5(1)*rx5 + m5(2)*ry5 + m5(3)*rz5).*rz5 - m5(3))./r5.^3;


% Remove singularities for all dipoles
threshold = 0.1;
Bx1(r1<threshold) = NaN; By1(r1<threshold) = NaN; Bz1(r1<threshold) = NaN;
Bx2(r2<threshold) = NaN; By2(r2<threshold) = NaN; Bz2(r2<threshold) = NaN;
Bx3(r3<threshold) = NaN; By3(r3<threshold) = NaN; Bz3(r3<threshold) = NaN;
Bx4(r4<threshold) = NaN; By4(r4<threshold) = NaN; Bz4(r4<threshold) = NaN;
Bx5(r5<threshold) = NaN; By5(r5<threshold) = NaN; Bz5(r5<threshold) = NaN;

% Sum the magnetic fields from all dipoles
Bx_total = Bx1 + Bx2 + Bx3 + Bx4 + Bx5;
By_total = By1 + By2 + By3 + By4 + By5;
Bz_total = Bz1 + Bz2 + Bz3 + Bz4 + Bz5;

% Find indices for each dipole
[idx1, idy1, idz1] = findClosestGridPoint(x, y, z, p1);
[idx2, idy2, idz2] = findClosestGridPoint(x, y, z, p2);
[idx3, idy3, idz3] = findClosestGridPoint(x, y, z, p3);
[idx4, idy4, idz4] = findClosestGridPoint(x, y, z, p4);
[idx5, idy5, idz5] = findClosestGridPoint(x, y, z, p5);

%%
% Compute Torques
% Torque = m x B (cross product)
% Calculate the magnetic field at each dipole's position due to other dipoles


% Calculate Torques Dipole1
T1_2 = f_getTorque(Bx2,By2,Bz2,idx1, idy1, idz1, m1);
T1_3 = f_getTorque(Bx3,By3,Bz3,idx1, idy1, idz1,m1);
T1_4 = f_getTorque(Bx4,By4,Bz4,idx1, idy1, idz1,m1);
T1_5 = f_getTorque(Bx5,By5,Bz5,idx1, idy1, idz1,m1);

% Sum the torques to get the total torque on dipole 1
torque1 = T1_2+ T1_3 + T1_4+ T1_5;

% Calculate Torques
T2_1 = f_getTorque(Bx1,By1,Bz1,idx2, idy2, idz2,m2);
T2_3 = f_getTorque(Bx3,By3,Bz3,idx2, idy2, idz2,m2);
T2_4 = f_getTorque(Bx4,By4,Bz4,idx2, idy2, idz2,m2);
T2_5 = f_getTorque(Bx5,By5,Bz5,idx2, idy2, idz2,m2);

% Sum the torques to get the total torque on dipole 2
torque2 = T2_1+ T2_3 + T2_4+ T2_5;


% Calculate Torques
T3_1 = f_getTorque(Bx1,By1,Bz1,idx3, idy3, idz3,m3);
T3_2 = f_getTorque(Bx2,By2,Bz2,idx3, idy3, idz3,m3);
T3_4 = f_getTorque(Bx4,By4,Bz4,idx3, idy3, idz3,m3);
T3_5 = f_getTorque(Bx5,By5,Bz5,idx3, idy3, idz3,m3);

% Sum the torques to get the total torque on dipole 3
torque3 = T3_1+ T3_2 + T3_4+ T3_5;


% Calculate Torques
T4_1 = f_getTorque(Bx1,By1,Bz1,idx4, idy4, idz4,m4);
T4_2 = f_getTorque(Bx2,By2,Bz2,idx4, idy4, idz4,m4);
T4_3 = f_getTorque(Bx3,By3,Bz3,idx4, idy4, idz4,m4);
T4_5 = f_getTorque(Bx5,By5,Bz5,idx4, idy4, idz4,m4);

% Sum the torques to get the total torque on dipole 3
torque4 = T4_1+ T4_2 + T4_3+ T4_5;


% Calculate Torques
T5_1 = f_getTorque(Bx1,By1,Bz1,idx5, idy5, idz5,m5);
T5_2 = f_getTorque(Bx2,By2,Bz2,idx5, idy5, idz5,m5);
T5_3 = f_getTorque(Bx3,By3,Bz3,idx5, idy5, idz5,m5);
T5_4 = f_getTorque(Bx4,By4,Bz4,idx5, idy5, idz5,m5);

% Sum the torques to get the total torque on dipole 5
torque5 = T5_1+ T5_2 + T5_3+ T5_4;


%%
%Scale torques
torque1 = torque1*TSF;
torque2 = torque2*TSF;
torque3 = torque3*TSF;
torque4 = torque4*TSF;
torque5 = torque5*TSF;

%% Compute forces
F1 = f_getForce(Bx_total, By_total, Bz_total, idx1, idy1, idz1, m1, dx, dy, dz)*TSF;
F2 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m2, dx, dy, dz)*TSF;
F3 = f_getForce(Bx_total, By_total, Bz_total, idx3, idy3, idz3, m3, dx, dy, dz)*TSF;
F4 = f_getForce(Bx_total, By_total, Bz_total, idx4, idy4, idz4, m4, dx, dy, dz)*TSF;
F5 = f_getForce(Bx_total, By_total, Bz_total, idx5, idy5, idz5, m5, dx, dy, dz)*TSF;

%%
%Plotting

% Plotting the total magnetic field (3D)
figure(1)
quiver3(x, y, z, Bx_total, By_total, Bz_total);
hold on;
plot3(p1(1), p1(2), p1(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % First dipole position
plot3(p2(1), p2(2), p2(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Second dipole position
plot3(p3(1), p3(2), p3(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Third dipole position
plot3(p4(1), p4(2), p4(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Fourth dipole position
plot3(p5(1), p5(2), p5(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Fifth dipole position
hold off;
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Total Magnetic Field Vector Field of Two Dipoles');


% Plotting the torques
figure(2)
quiver3(x, y, z, Bx_total, By_total, Bz_total);
hold on;
quiver3(p1(1), p1(2), p1(3), torque1(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p1(1), p1(2), p1(3), 0, torque1(2), 0, 'Color', 'r'); % Y-component
quiver3(p1(1), p1(2), p1(3), 0, 0, torque1(3), 'Color', 'r'); % Z-component
quiver3(p2(1), p2(2), p2(3), torque2(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p2(1), p2(2), p2(3), 0, torque2(2), 0, 'Color', 'r'); % Y-component
quiver3(p2(1), p2(2), p2(3), 0, 0, torque2(3), 'Color', 'r'); % Z-component
quiver3(p3(1), p3(2), p3(3), torque3(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p3(1), p3(2), p3(3), 0, torque3(2), 0, 'Color', 'r'); % Y-component
quiver3(p3(1), p3(2), p3(3), 0, 0, torque3(3), 'Color', 'r'); % Z-component
quiver3(p4(1), p4(2), p4(3), torque4(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p4(1), p4(2), p4(3), 0, torque4(2), 0, 'Color', 'r'); % Y-component
quiver3(p4(1), p4(2), p4(3), 0, 0, torque4(3), 'Color', 'r'); % Z-component
quiver3(p5(1), p5(2), p5(3), torque5(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p5(1), p5(2), p5(3), 0, torque5(2), 0, 'Color', 'r'); % Y-component
quiver3(p5(1), p5(2), p5(3), 0, 0, torque5(3), 'Color', 'r'); % Z-component
hold off;
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Total Magnetic Field Vector Field of Two Dipoles');

% Plotting the force components
figure(3)
quiver3(x, y, z, Bx_total, By_total, Bz_total);
hold on;

% Plotting the x, y, z components of the force for each dipole
% For Dipole 1
quiver3(p1(1), p1(2), p1(3), F1(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p1(1), p1(2), p1(3), 0, F1(2), 0, 'Color', 'g'); % Y-component
quiver3(p1(1), p1(2), p1(3), 0, 0, F1(3), 'Color', 'b'); % Z-component

% For Dipole 2
quiver3(p2(1), p2(2), p2(3), F2(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p2(1), p2(2), p2(3), 0, F2(2), 0, 'Color', 'g'); % Y-component
quiver3(p2(1), p2(2), p2(3), 0, 0, F2(3), 'Color', 'b'); % Z-component

% For Dipole 3
quiver3(p3(1), p3(2), p3(3), F3(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p3(1), p3(2), p3(3), 0, F3(2), 0, 'Color', 'g'); % Y-component
quiver3(p3(1), p3(2), p3(3), 0, 0, F3(3), 'Color', 'b'); % Z-component

% For Dipole 4
quiver3(p4(1), p4(2), p4(3), F4(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p4(1), p4(2), p4(3), 0, F4(2), 0, 'Color', 'g'); % Y-component
quiver3(p4(1), p4(2), p4(3), 0, 0, F4(3), 'Color', 'b'); % Z-component

% For Dipole 5
quiver3(p5(1), p5(2), p5(3), F5(1), 0, 0, 'Color', 'r'); % X-component
quiver3(p5(1), p5(2), p5(3), 0, F5(2), 0, 'Color', 'g'); % Y-component
quiver3(p5(1), p5(2), p5(3), 0, 0, F5(3), 'Color', 'b'); % Z-component

hold off;
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Components of Magnetic Forces on Dipoles');

% Function to find the closest grid point
function [idx, idy, idz] = findClosestGridPoint(x, y, z, p)
    [~, idx] = min(abs(x(1,:,1) - p(1)));
    [~, idy] = min(abs(y(:,1,1) - p(2)));
    [~, idz] = min(abs(z(1,1,:) - p(3)));
end