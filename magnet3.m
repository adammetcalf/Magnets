close all;
clear;
clc;

% Parameters
mu0 = 4*pi*1e-7; % Permeability of free space
m1 = [4, 0, 0];   % Magnetic moment vector of the first dipole (aligned along z-axis)
m2 = [0, 0, 4];   % Magnetic moment vector of the second dipole (same orientation)
m3 = [0, 4, 0];   % Magnetic moment vector of the third dipole (same orientation)
m4 = [0, 0, 4];   % Magnetic moment vector of the fourth dipole (same orientation)
m5 = [4, 0, 0];   % Magnetic moment vector of the fourth dipole (same orientation)

% Position offset for the second dipole (0.15 meters above the first dipole)
dipole_offset = [0, 0, 0.15]; 
dipole2_offset = [0, 0, 0.3]; 
dipole3_offset = [0, 0, 0.45];
dipole4_offset = [0, 0, 0.6];

% Grid
[x, y, z] = meshgrid(linspace(-0.3, 0.3, 25), linspace(-0.3, 0.3, 25), linspace(-0.3, 0.9, 25));

% Calculate field components for the first dipole
r1 = sqrt(x.^2 + y.^2 + z.^2);
rx1 = x./r1; ry1 = y./r1; rz1 = z./r1;

Bx1 = mu0/(4*pi) * (3*(m1(1)*rx1 + m1(2)*ry1 + m1(3)*rz1).*rx1 - m1(1))./r1.^3;
By1 = mu0/(4*pi) * (3*(m1(1)*rx1 + m1(2)*ry1 + m1(3)*rz1).*ry1 - m1(2))./r1.^3;
Bz1 = mu0/(4*pi) * (3*(m1(1)*rx1 + m1(2)*ry1 + m1(3)*rz1).*rz1 - m1(3))./r1.^3;

% Calculate field components for the second dipole
x2 = x - dipole_offset(1);
y2 = y - dipole_offset(2);
z2 = z - dipole_offset(3);
r2 = sqrt(x2.^2 + y2.^2 + z2.^2);
rx2 = x2./r2; ry2 = y2./r2; rz2 = z2./r2;

Bx2 = mu0/(4*pi) * (3*(m2(1)*rx2 + m2(2)*ry2 + m2(3)*rz2).*rx2 - m2(1))./r2.^3;
By2 = mu0/(4*pi) * (3*(m2(1)*rx2 + m2(2)*ry2 + m2(3)*rz2).*ry2 - m2(2))./r2.^3;
Bz2 = mu0/(4*pi) * (3*(m2(1)*rx2 + m2(2)*ry2 + m2(3)*rz2).*rz2 - m2(3))./r2.^3;

% Calculate field components for the thrid dipole
x3 = x - dipole2_offset(1);
y3 = y - dipole2_offset(2);
z3 = z - dipole2_offset(3);
r3 = sqrt(x3.^2 + y3.^2 + z3.^2);
rx3 = x3./r3; ry3 = y3./r3; rz3 = z3./r3;

Bx3 = mu0/(4*pi) * (3*(m3(1)*rx3 + m3(2)*ry3 + m3(3)*rz3).*rx3 - m3(1))./r3.^3;
By3 = mu0/(4*pi) * (3*(m3(1)*rx3 + m3(2)*ry3 + m3(3)*rz3).*ry3 - m3(2))./r3.^3;
Bz3 = mu0/(4*pi) * (3*(m3(1)*rx3 + m3(2)*ry3 + m3(3)*rz3).*rz3 - m3(3))./r3.^3;

% Calculate field components for the Fourth dipole
x4 = x - dipole3_offset(1);
y4 = y - dipole3_offset(2);
z4 = z - dipole3_offset(3);
r4 = sqrt(x4.^2 + y4.^2 + z4.^2);
rx4 = x4./r4; ry4 = y4./r4; rz4 = z4./r4;

Bx4 = mu0/(4*pi) * (3*(m4(1)*rx4 + m4(2)*ry4 + m4(3)*rz4).*rx4 - m4(1))./r4.^3;
By4 = mu0/(4*pi) * (3*(m4(1)*rx4 + m4(2)*ry4 + m4(3)*rz4).*ry4 - m4(2))./r4.^3;
Bz4 = mu0/(4*pi) * (3*(m4(1)*rx4 + m4(2)*ry4 + m4(3)*rz4).*rz4 - m4(3))./r4.^3;

% Calculate field components for the Fourth dipole
x5 = x - dipole4_offset(1);
y5 = y - dipole4_offset(2);
z5 = z - dipole4_offset(3);
r5 = sqrt(x5.^2 + y5.^2 + z5.^2);
rx5 = x5./r5; ry5 = y5./r5; rz5 = z5./r5;

Bx5 = mu0/(4*pi) * (3*(m5(1)*rx5 + m5(2)*ry5 + m5(3)*rz5).*rx5 - m5(1))./r5.^3;
By5 = mu0/(4*pi) * (3*(m5(1)*rx5 + m5(2)*ry5 + m5(3)*rz5).*ry5 - m5(2))./r5.^3;
Bz5 = mu0/(4*pi) * (3*(m5(1)*rx5 + m5(2)*ry5 + m5(3)*rz5).*rz5 - m5(3))./r5.^3;


% Remove singularities for both dipoles
threshold = 0.1;
Bx1(r1<threshold) = NaN; By1(r1<threshold) = NaN; Bz1(r1<threshold) = NaN;
Bx2(r2<threshold) = NaN; By2(r2<threshold) = NaN; Bz2(r2<threshold) = NaN;
Bx3(r3<threshold) = NaN; By3(r3<threshold) = NaN; Bz3(r3<threshold) = NaN;
Bx4(r4<threshold) = NaN; By4(r4<threshold) = NaN; Bz4(r4<threshold) = NaN;
Bx5(r5<threshold) = NaN; By5(r5<threshold) = NaN; Bz5(r5<threshold) = NaN;

% Sum the magnetic fields from both dipoles
Bx_total = Bx1 + Bx2 + Bx3 + Bx4 + Bx5;
By_total = By1 + By2 + By3 + By4 + By5;
Bz_total = Bz1 + Bz2 + Bz3 + Bz4 + Bz5;

% Plotting the total magnetic field (3D)
figure(1)
quiver3(x, y, z, Bx_total, By_total, Bz_total);
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % First dipole position
plot3(0, 0, 0.15, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Second dipole position
plot3(0, 0, 0.3, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Third dipole position
plot3(0, 0, 0.45, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Fourth dipole position
plot3(0, 0, 0.6, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Fourth dipole position
hold off;
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Total Magnetic Field Vector Field of Two Dipoles');


