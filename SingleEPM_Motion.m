close all;
clear;
clc;

%% Visualisation
MagfieldView = true;
TorqueView = false;
Forceview = false;
VectorFieldsVisible = false;

% Initialize Video Writer
videoFilename = 'magnetic_field_simulation.avi';
writerObj = VideoWriter(videoFilename);
open(writerObj);

%% Parameters
mu0 = 4*pi*1e-7;                % Permeability of free space

% Vector upon which magnetic moments act
EPM_Diretion = [1, 0, 0];       % Magnetic moment vector of the first dipole (aligned along z-axis)
initial_dipole_Direction = [0,1,0];     % Define unit vectors in local magnet frame (assuming the magnet's north pole points along the local Z-axis)

% Dipole magentic moments
mu_EPM = 970.1;                 %Magentic moment of the EPM (Magnitude)
mu_dipole = 1;              %Magentic moment of the dipole (Magnitude)

% Dipole Positions
EPM_Pos = [-1, 0, 0];            % Position of the first dipole
dipole_pos = [0, 0, 0];                   % Position of the dipole

% Initial orientation in Euler angles (roll, pitch, yaw)
orientation_dipole = [0, 0, 0]; % Example values

%Force/Torque Scaling Factor
TSF = 100;

%For orientation visualisation of dipole
axis_length = 0.02; 

%Initialise Dipole Mass
Dipole_Mass = 0.01; %10grams

% Initialize Velocity and Acceleration
velocity_dipole = [0, 0, 0];
acceleration_dipole = [0, 0, 0];

% Initialize angular velocity
omega = [0, 0, 0]; % Angular velocity in radians per second
alpha = [0, 0, 0];

% Initialse Moment of Inertia
I = 0.01; 

% Rotational Damping coefficient
b = 0.001;

% Time step and simulation duration
dt = 0.1; % 10 milliseconds per step
total_time = 45; % total time of 1 second

%Magnetic field fthreshold
threshold = 0.1;

%% Calculate magnetic moment vectors
m_EPM = mu_EPM * EPM_Diretion;
m_dipole = mu_dipole * initial_dipole_Direction;

%% Grid for visualising
% Grid
[x, y, z] = meshgrid(linspace(-0.1, 0.1, 25), linspace(-0.1, 0.1, 25), linspace(-0.1, 0.1, 25));

dx = (0.1-(-0.1))/(25-1);
dy = (0.1-(-0.1))/(25-1);
dz = (0.1-(-0.1))/(25-1);

% Main loop for simulation
for t = 0:dt:total_time

    %% Calculate Magnetic field
    
    % Calculate field components for the EPM 
    [Bx1,By1,Bz1] = calculateField(EPM_Pos,mu0,m_EPM,x,y,z,threshold);
    
    % Calculate field components for the dipole   
    [Bx2,By2,Bz2] = calculateField(dipole_pos,mu0,m_dipole,x,y,z,threshold);
    
    % Sum the magnetic fields from all dipoles
    Bx_total = Bx1;%+ Bx2;                            
    By_total = By1;%+ By2;
    Bz_total = Bz1;%+ Bz2;
        
    % Find meshgrid indices for each dipole
    [idx1, idy1, idz1] = findClosestGridPoint(x, y, z, EPM_Pos);
    [idx2, idy2, idz2] = findClosestGridPoint(x, y, z, dipole_pos);
    
    %% Compute Torques on Dipole
    T1 = f_getTorque(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m_dipole);

    % Calculate damping torque (proportional to angular velocity)
    damping_torque = -b * omega;
    
    % Apply damping torque
    T1 = T1 + damping_torque;
    
    %% Compute forces
    F1 = f_getForce(Bx_total, By_total, Bz_total, idx2, idy2, idz2, m_dipole, dx, dy, dz);

    %% Dynamics
    % Update velocity and position using Euler's method
    acceleration_dipole = F1 / Dipole_Mass;
    velocity_dipole = velocity_dipole + acceleration_dipole * dt;

    % Calculate Angular Acceleration
    alpha = T1 / I;

    % Update Angular Velocity
    omega = omega + alpha * dt;

    % Update Orientation
    % For simplicity, this example assumes small rotations
    % For larger rotations, use rotation matrices or quaternions
    orientation_dipole = orientation_dipole + omega * dt;

    % Calculate Rotation Matrix
    roll = orientation_dipole(1);
    pitch = orientation_dipole(2);
    yaw = orientation_dipole(3);

    R_x = [1, 0, 0;
           0, cos(roll), -sin(roll);
           0, sin(roll), cos(roll)];

    R_y = [cos(pitch), 0, sin(pitch);
           0, 1, 0;
           -sin(pitch), 0, cos(pitch)];

    R_z = [cos(yaw), -sin(yaw), 0;
           sin(yaw), cos(yaw), 0;
           0, 0, 1];

    R = R_z * R_y * R_x;

    % Update Dipole Position
    dipole_pos = dipole_pos + velocity_dipole * dt;

    % Homogeneous Transformation Matrix
    T = [R, dipole_pos'; % Concatenate rotation matrix and translation vector
         0, 0, 0, 1];    % Append the row [0 0 0 1]

    % Apply Transformation to the Dipole's Direction Vector: 
    initial_direction_x = [1; 0; 0; 1]; % Initial direction of X-axis
    initial_direction_y = [0; 1; 0; 1]; % Initial direction of Y-axis
    initial_direction_z = [0; 0; 1; 1]; % Initial direction of Z-axis

    updated_direction_x = T * initial_direction_x;
    updated_direction_y = T * initial_direction_y;
    updated_direction_z = T * initial_direction_z;

    %% Recalculate magnetic moment with the updated direction
    global_dipole_Direction = R * initial_dipole_Direction'; % Transform the local direction to global frame

    m_dipole = mu_dipole * global_dipole_Direction;

    %% Create vector fields for plot
    
    % Define the target range
    lower_bound_x = -0.07; upper_bound_x = 0.07;
    lower_bound_y = -0.1;  upper_bound_y = 0.1;
    lower_bound_z = -0.1;  upper_bound_z = 0.1;
    
    if VectorFieldsVisible
        %create gap in bx etc for clearer plots
        [plotx, ploty, plotz] = removeVectors(x,y,z,Bx_total,By_total,Bz_total,upper_bound_x,lower_bound_x,upper_bound_y,lower_bound_y,upper_bound_z,lower_bound_z);
        
        %Remove vectors from Bx etc for clearer plots.
        [decX,decY,decZ] = decimateVectors(Bx_total,By_total,Bz_total,100);
    end
    %% Visualisation

    if MagfieldView
    % Plotting the total magnetic field (3D)
    figure(1)
    if VectorFieldsVisible
        quiver3(x, y, z, decX, decY, decZ);
        hold on;
    end
    plot3(EPM_Pos(1), EPM_Pos(2), EPM_Pos(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'b'); % First dipole position
    hold on
    plot3(dipole_pos(1), dipole_pos(2),dipole_pos(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Second dipole position
    hold on
    % Draw axes for dipole based on the updated direction
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), axis_length * updated_direction_x(1), axis_length * updated_direction_x(2), axis_length * updated_direction_x(3), 'r'); % X-axis
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), axis_length * updated_direction_y(1), axis_length * updated_direction_y(2), axis_length * updated_direction_y(3), 'b'); % Y-axis
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), axis_length * updated_direction_z(1), axis_length * updated_direction_z(2), axis_length * updated_direction_z(3), 'g'); % Z-axis
    axis equal
    axis([-1 0.1 -0.1 0.1 -0.1 0.1]);
    hold off
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');

    % Set the viewing angle (azimuth, elevation)
    view(-42.0000, 39); 

    % Capture the frame
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    end

    if TorqueView
    %Plotting the torques
    figure(2)
    if VectorFieldsVisible
        quiver3(x, y, z, plotx, ploty, plotz);
        hold on
    end
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), T1(1)*TSF, 0, 0, 'Color', 'r'); % X-component
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, T1(2)*TSF, 0, 'Color', 'r'); % Y-component
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, 0, T1(3)*TSF, 'Color', 'r'); % Z-component
    hold off
    axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis')  

    % Set the viewing angle (azimuth, elevation)
    view(30, 45); % Example values, adjust as needed
    end

    if Forceview
    % Plotting the force components for the dipole
    figure(3)
    if VectorFieldsVisible
        quiver3(x, y, z, plotx, ploty, plotz);
        hold on
    end
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), F1(1)*TSF, 0, 0, 'Color', 'r'); % X-component
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, F1(2)*TSF, 0, 'Color', 'g'); % Y-component
    quiver3(dipole_pos(1), dipole_pos(2), dipole_pos(3), 0, 0, F1(3)*TSF, 'Color', 'b'); % Z-component
    hold off
    axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    end

end 

% Close the Video Writer
close(writerObj);

function [Bx,By,Bz] = calculateField(Pos,mu0,m,x,y,z,threshold)
    % Calculate field components for a dipole
    x1 = x - Pos(1);
    y1 = y - Pos(2);
    z1 = z - Pos(3);
    r1 = sqrt(x1.^2 + y1.^2 + z1.^2);
    rx1 = x1./r1; ry1 = y1./r1; rz1 = z1./r1;
    
    Bx = mu0/(4*pi) * (3*(m(1)*rx1 + m(2)*ry1 + m(3)*rz1).*rx1 - m(1))./r1.^3;
    By = mu0/(4*pi) * (3*(m(1)*rx1 + m(2)*ry1 + m(3)*rz1).*ry1 - m(2))./r1.^3;
    Bz = mu0/(4*pi) * (3*(m(1)*rx1 + m(2)*ry1 + m(3)*rz1).*rz1 - m(3))./r1.^3;

    %remove singularities
    Bx(r1<threshold) = NaN; By(r1<threshold) = NaN; Bz(r1<threshold) = NaN;
end


