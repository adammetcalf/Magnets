function torque = f_getTorque(Bx, By, Bz, x_index, y_index, z_index, m)
%This function returns the torques at a dipole wrt another dipole

% Calculate the magnetic field at dipole due to another dipole
B = [Bx(x_index,y_index,z_index), By(x_index,y_index,z_index), Bz(x_index,y_index,z_index)];
torque= cross(m, B);
end