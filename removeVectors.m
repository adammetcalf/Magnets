function[plotx, ploty, plotz] = removeVectors(x,y,z,Bx,By,Bz,xhigh,xlow,yhigh,ylow,zhigh,zlow)
%This fuction removes vector arrows from a vector filed for clearer
%plotting

plotx = Bx; ploty = By; plotz = Bz;

% Find indices corresponding to the range in each dimension
indices_x = x >= xlow & x <= xhigh;
indices_y = y >= ylow & y <= yhigh;
indices_z = z >= zlow & z <= zhigh;

% Combine indices for all three dimensions
combined_indices = indices_x & indices_y & indices_z;

% Set the vector field components to zero in the specified range
plotx(combined_indices) = 0;
ploty(combined_indices) = 0;
plotz(combined_indices) = 0;

end