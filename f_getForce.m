function force = f_getForce(Bx, By, Bz, x_index, y_index, z_index, m, dx, dy, dz)
%dx, dy, dz are the grid spacings in the x, y, and z directions, respectively.


    % Calculate the magnetic field at the dipole's position
    B = [Bx(x_index, y_index, z_index), By(x_index, y_index, z_index), Bz(x_index, y_index, z_index)];

    % Approximate the gradient of the magnetic field
    dBx_dx = (Bx(x_index+1, y_index, z_index) - Bx(x_index-1, y_index, z_index)) / (2*dx);
    dBy_dy = (By(x_index, y_index+1, z_index) - By(x_index, y_index-1, z_index)) / (2*dy);
    dBz_dz = (Bz(x_index, y_index, z_index+1) - Bz(x_index, y_index, z_index-1)) / (2*dz);

    % Calculate the force
    force = [m(1) * dBx_dx, m(2) * dBy_dy, m(3) * dBz_dz];
end
