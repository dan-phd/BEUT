function [E_s] = organizeScatteredField(operator_file, grid, in_shape )
% Convert the matrix output from C++ to an appropriate cube output
% containing the electric field at the corresponding grid points

% Load from C++
N_T=operator_file.N_T;
E = operator_file.E;

%% rearrange field to grid points
tmp = zeros(size(in_shape));
E_s = zeros([size(grid) N_T]);

for j=1:N_T
    
    tmp(~in_shape) = E(:,j);
    
    E_s(:,:,j) = reshape(tmp(1:numel(grid),:),size(grid));
end
