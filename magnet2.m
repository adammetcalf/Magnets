close all;
clear;
clc;

% Parameters
mu0 = 4*pi*1e-7; % Permeability of free space
m = [0, 0, 1];   % Magnetic moment vector (aligned along z-axis)

% Grid
[x, y, z] = meshgrid(linspace(-0.3, 0.3, 20), linspace(-0.3, 0.3, 20), linspace(-0.3, 0.3, 20));

% Calculate field components
r = sqrt(x.^2 + y.^2 + z.^2);
rx = x./r; ry = y./r; rz = z./r;

Bx = mu0/(4*pi) * (3*(m(1)*rx + m(2)*ry + m(3)*rz).*rx - m(1))./r.^3;
By = mu0/(4*pi) * (3*(m(1)*rx + m(2)*ry + m(3)*rz).*ry - m(2))./r.^3;
Bz = mu0/(4*pi) * (3*(m(1)*rx + m(2)*ry + m(3)*rz).*rz - m(3))./r.^3;

% Remove singularities
Bx(r<0.1) = NaN;
By(r<0.1) = NaN;
Bz(r<0.1) = NaN;

% Plotting
figure(1)
quiver3(x, y, z, Bx, By, Bz);
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r'); % First dipole position
hold off
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Magnetic Field Vector Field of a Dipole');

% Plotting
figure(2)
quiver(x, y,Bx, By);
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
title('Magnetic Field Vector Field of a Dipole');