close all;
clear;
clc;


%% Parameters
mu0 = 4*pi*1e-7;                % Permeability of free space

% Vector upon which magnetic moments act
EPM_Diretion = [1, 0, 0];       % Magnetic moment vector of the first dipole (aligned along z-axis)
dipole_Direction = [0,0,1];     % Define unit vectors in local magnet frame (assuming the magnet's north pole points along the local Z-axis)

% Dipole magentic moments
mu_EPM = 970.1;                 %Magentic moment of the EPM (Magnitude)
mu_dipole = 1;              %Magentic moment of the dipole (Magnitude)

% Dipole Positions
EPM_Pos = [-1, 0, 0];            % Position of the first dipole
dipole_pos = [0, 0, 0];                   % Position of the dipole

%Force/Torque Scaling Factor
TSF = 100;

%% Calculate magnetic moment vectors
m_EPM = mu_EPM * EPM_Diretion;
m_dipole = mu_dipole * dipole_Direction;

%% Grid for visualising
% Grid
[x, y, z] = meshgrid(linspace(-0.1, 0.1, 11), linspace(-0.1, 0.1, 11), linspace(-0.1, 0.1, 11));

dx = (0.1-(-0.1))/(15-1);
dy = (0.1-(-0.1))/(15-1);
dz = (0.1-(-0.1))/(15-1);

%% Calculate Magnetic field

% Calculate field components for the EPM 
x1 = x - EPM_Pos(1);
y1 = y - EPM_Pos(2);
z1 = z - EPM_Pos(3);
r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;

Bx1 = mu0/(4*pi) * (3*(m_EPM(1)*rx1 + m_EPM(2)*ry1 + m_EPM(3)*rz1).*rx1 - m_EPM(1))./r1.^3;
By1 = mu0/(4*pi) * (3*(m_EPM(1)*rx1 + m_EPM(2)*ry1 + m_EPM(3)*rz1).*ry1 - m_EPM(2))./r1.^3;
Bz1 = mu0/(4*pi) * (3*(m_EPM(1)*rx1 + m_EPM(2)*ry1 + m_EPM(3)*rz1).*rz1 - m_EPM(3))./r1.^3;

% Calculate field components for the dipole
x2 = x - dipole_pos(1);
y2 = y - dipole_pos(2);
z2 = z - dipole_pos(3);
r2 = sqrt(x2.^2 + y2.^2 + z2.^2);
rx2 = x2./r2; ry2 = y2./r2; rz2 = z2./r2;

Bx2 = mu0/(4*pi) * (3*(m_dipole(1)*rx2 + m_dipole(2)*ry2 + m_dipole(3)*rz2).*rx2 - m_dipole(1))./r2.^3;
By2 = mu0/(4*pi) * (3*(m_dipole(1)*rx2 + m_dipole(2)*ry2 + m_dipole(3)*rz2).*ry2 - m_dipole(2))./r2.^3;
Bz2 = mu0/(4*pi) * (3*(m_dipole(1)*rx2 + m_dipole(2)*ry2 + m_dipole(3)*rz2).*rz2 - m_dipole(3))./r2.^3;

% Remove singularities for all dipoles
threshold = 0.1;
Bx1(r1<threshold) = NaN; By1(r1<threshold) = NaN; Bz1(r1<threshold) = NaN;
Bx2(r2<threshold) = NaN; By2(r2<threshold) = NaN; Bz2(r2<threshold) = NaN;

% Sum the magnetic fields from all dipoles
Bx_total = Bx1;%+ Bx2;                            
By_total = By1;%+ By2;
Bz_total = Bz1;%+ Bz2;


% Find meshgrid indices for each dipole
[idx1, idy1, idz1] = findClosestGridPoint(x, y, z, EPM_Pos);
[idx2, idy2, idz2] = findClosestGridPoint(x, y, z, dipole_pos);

%% Compute Torques on Dipole
T1 = f_getTorque(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m_dipole);

%% Compute forces
F1 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m_dipole, dx, dy, dz);

%% Create vector fields for plot

% Define the target range
lower_bound_x = -0.07;
upper_bound_x = 0.07;
lower_bound_y = -0.1;
upper_bound_y = 0.1;
lower_bound_z = -0.1;
upper_bound_z = 0.1;

%create gap in bx etc for clearer plots
[plotx, ploty, plotz] = removeVectors(x,y,z,Bx_total,By_total,Bz_total,upper_bound_x,lower_bound_x,upper_bound_y,lower_bound_y,upper_bound_z,lower_bound_z);

%Remove vecotrs from Bx etc for clearer plots.
[decX,decY,decZ] = decimateVectors(Bx_total,By_total,Bz_total,4);

%% Visualisation
% Plotting the total magnetic field (3D)
figure(1)
quiver3(x, y, z, decX, decY, decZ);
hold on;
%plot3(EPM_Pos(1), EPM_Pos(2), EPM_Pos(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'b'); % First dipole position
plot3(dipole_pos(1), dipole_pos(2),dipole_pos(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Second dipole position
axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
hold off;
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

% Plotting the torques
figure(2)
quiver3(x, y, z, plotx, ploty, plotz);
hold on
quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), T1(1)*TSF, 0, 0, 'Color', 'r'); % X-component
quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, T1(2)*TSF, 0, 'Color', 'r'); % Y-component
quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, 0, T1(3)*TSF, 'Color', 'r'); % Z-component
hold off
axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');


% Plotting the force components for the dipole
figure(3)
quiver3(x, y, z, plotx, ploty, plotz);
hold on
quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), F1(1)*TSF, 0, 0, 'Color', 'r'); % X-component
quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, F1(2)*TSF, 0, 'Color', 'g'); % Y-component
quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, 0, F1(3)*TSF, 'Color', 'b'); % Z-component
hold off
axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
