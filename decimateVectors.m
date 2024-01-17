function [bx_out, by_out, bz_out] = decimateVectors(bx, by, bz, n)
    % Ensure bx, by, bz have the same size
    assert(isequal(size(bx), size(by), size(bz)), 'Input matrices must have the same size.');
    
    % Initialize output matrices
    bx_out = bx;
    by_out = by;
    bz_out = bz;
    
    % Get the size of the input matrices
    [nx, ny, nz] = size(bx);
    
    % Loop through each element
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                % Calculate the sum of indices
                sum_idx = ix + iy + iz;
                
                % If the sum of indices is not divisible by 3, set to zero
                if mod(sum_idx, n) ~= 0
                    bx_out(ix, iy, iz) = 0;
                    by_out(ix, iy, iz) = 0;
                    bz_out(ix, iy, iz) = 0;
                end
            end
        end
    end
end