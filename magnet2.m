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


% Calculate the gradients
[gradBx_x, gradBx_y, gradBx_z] = gradient(Bx, x(1,:,1), y(:,1,1), z(1,1,:));
[gradBy_x, gradBy_y, gradBy_z] = gradient(By, x(1,:,1), y(:,1,1), z(1,1,:));
[gradBz_x, gradBz_y, gradBz_z] = gradient(Bz, x(1,:,1), y(:,1,1), z(1,1,:));

% Compute the magnitude of the gradient
gradB_magnitude = sqrt(gradBx_x.^2 + gradBx_y.^2 + gradBx_z.^2 + ...
                       gradBy_x.^2 + gradBy_y.^2 + gradBy_z.^2 + ...
                       gradBz_x.^2 + gradBz_y.^2 + gradBz_z.^2);

% Visualize the gradient magnitude
% For a 2D visualization on the xy-plane (z=0) choose one slice
sliceIndex = find(abs(z(1,1,:) - 0) == min(abs(z(1,1,:) - 0)), 1);
figure(3);
imagesc(linspace(-0.3, 0.3, 20), linspace(-0.3, 0.3, 20), squeeze(gradB_magnitude(:,:,sliceIndex)));
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('Magnitude of the Gradient of the Magnetic Field');

% Visualize the gradient magnitude in 3D
figure(4);
% Define slice positions (modify as needed)
sx = 0; % Slice along X at x=0
sy = 0; % Slice along Y at y=0
sz = [min(z(:)) max(z(:))]; % Slices along Z covering the full range
slice(x, y, z, gradB_magnitude, sx, sy, sz);
shading interp; % Interpolates colors across faces of slice
colorbar; % Show color scale
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Visualization of the Gradient Magnitude of the Magnetic Field');


% Set up the figure for animation
figure(5);
colormap jet; % Set the colormap

% Define fixed axis limits based on your data range
xLimits = [min(x(:)), max(x(:))];
yLimits = [min(y(:)), max(y(:))];
zLimits = [min(z(:)), max(z(:))];

% Animation parameters
zSlices = linspace(min(z(:)), max(z(:)), 50); % Define z-slices to animate through
nFrames = length(zSlices);

% Prepare VideoWriter
v = VideoWriter('magnetic_field_gradient_animation.avi');
open(v);

% Animate through z-slices
for i = 1:nFrames
    % Clear the current figure
    clf;
    
    % Slice parameters for current frame
    sz = zSlices(i);
    
    % Draw the current slice
    quiver3(x, y, z, Bx, By, Bz);
    hold on;
    hs = slice(x, y, z, gradB_magnitude, [], [], sz);
    shading interp; % Smooths the color transition
    colorbar; % Show color scale
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    % Set fixed axis limits
    xlim(xLimits);
    ylim(yLimits);
    zlim(zLimits);
    title(sprintf('Gradient Magnitude at Z = %.2f', sz));
    clim([min(gradB_magnitude(:)), max(gradB_magnitude(:))]); % Fix color axis

    % Adjust the alpha value for each slice
    for j = 1:length(hs)
        hs(j).FaceAlpha = 0.5; % Set transparency to 0.5
    end

    drawnow; % Update the figure window
    
    % Capture the current frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);


% Set up the figure for animation x axis
figure(6);
colormap jet; % Set the colormap 

% Define fixed axis limits based on your data range
xLimits = [min(x(:)), max(x(:))];
yLimits = [min(y(:)), max(y(:))];
zLimits = [min(z(:)), max(z(:))];

% Animation parameters
xSlices = linspace(min(x(:)), max(x(:)), 50); % Define x-slices to animate through
nFrames = length(xSlices);

% Prepare VideoWriter if you wish to save the animation
v = VideoWriter('magnetic_field_gradient_animation_x-axis.avi');
open(v);

% Animate through x-slices
for i = 1:nFrames
    % Clear the current figure
    clf;
    
    % Draw the current slice
    quiver3(x, y, z, Bx, By, Bz);
    hold on;
    hs = slice(x, y, z, gradB_magnitude, xSlices(i), [], []);
    shading interp; % Smooths the color transition
    colorbar; % Show color scale
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(sprintf('Gradient Magnitude at X = %.2f', xSlices(i)));
    
    % Set fixed axis limits
    xlim(xLimits);
    ylim(yLimits);
    zlim(zLimits);
    
    % fix the color axis for consistent coloring
    clim([min(gradB_magnitude(:)), max(gradB_magnitude(:))]);

    % Adjust the alpha value for each slice
    for j = 1:length(hs)
        hs(j).FaceAlpha = 0.5; % Set transparency to 0.5
    end

    drawnow; % Update the figure window
    
    % Capture the current frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);

% Set up the figure for animation x axis
figure(7);
colormap jet; % Set the colormap

% Define fixed axis limits based on your data range
xLimits = [min(x(:)), max(x(:))];
yLimits = [min(y(:)), max(y(:))];
zLimits = [min(z(:)), max(z(:))];

% Animation parameters
YSlices = linspace(min(y(:)), max(y(:)), 50); % Define y-slices to animate through
nFrames = length(YSlices);

% Prepare VideoWriter if you wish to save the animation
v = VideoWriter('magnetic_field_gradient_animation_y-axis.avi');
open(v);

% Animate through y-slices
for i = 1:nFrames
    % Clear the current figure
    clf;
    
    % Draw the current slice
    quiver3(x, y, z, Bx, By, Bz);
    hold on;
    hs = slice(x, y, z, gradB_magnitude,[], YSlices(nFrames+1-i),  []);
    shading interp; % Smooths the color transition
    colorbar; % Show color scale
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(sprintf('Gradient Magnitude at Y = %.2f', YSlices(nFrames+1-i)));
    
    % Set fixed axis limits
    xlim(xLimits);
    ylim(yLimits);
    zlim(zLimits);
    
    % fix the color axis for consistent coloring
    clim([min(gradB_magnitude(:)), max(gradB_magnitude(:))]);

    % Adjust the alpha value for each slice
    for j = 1:length(hs)
        hs(j).FaceAlpha = 0.5; % Set transparency to 0.5
    end
    
    hold off; % Release the figure for next update
    
    drawnow; % Update the figure window
    
    % Capture the current frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);

% Set up the figure for animation
figure(8);
colormap jet; % Set the colormap

% Define fixed axis limits based on your data range
xLimits = [min(x(:)), max(x(:))];
yLimits = [min(y(:)), max(y(:))];
zLimits = [min(z(:)), max(z(:))];

% Animation parameters
nFrames = 50; % Define number of frames for the animation
xSlices = linspace(min(x(:)), max(x(:)), nFrames);
ySlices = linspace(min(y(:)), max(y(:)), nFrames);
zSlices = linspace(min(z(:)), max(z(:)), nFrames);

% Prepare VideoWriter if you wish to save the animation
v = VideoWriter('magnetic_field_gradient_animation_xyz-axis_transparent.avi');
open(v);

% Animate through slices
for i = 1:nFrames
    % Clear the current figure
    clf;
    
    % Draw the magnetic field vectors
    quiver3(x, y, z, Bx, By, Bz, 'k');
    hold on;
    
    % Overlay the slices of the gradient magnitude
    hs = slice(x, y, z, gradB_magnitude, xSlices(i), ySlices(nFrames+1-i), zSlices(i));
    shading interp; % Smooths the color transition
    colorbar; % Show color scale
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(sprintf('Gradient Magnitude at Frame %d', i));
    
    % Set fixed axis limits
    xlim(xLimits);
    ylim(yLimits);
    zlim(zLimits);
    
    % fix the color axis for consistent coloring
    clim([min(gradB_magnitude(:)), max(gradB_magnitude(:))]);
    
    % Adjust the alpha value for each slice
    for j = 1:length(hs)
        hs(j).FaceAlpha = 0.5; % Set transparency to 0.5
    end
    
    hold off; % Release the figure for next update
    
    drawnow; % Update the figure window
    
    % Capture the current frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);


