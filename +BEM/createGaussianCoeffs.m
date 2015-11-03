function [ Gcoeffs ] = createGaussianCoeffs(so,wo,si,wi,max_deg_test,max_deg_basis)
% Using the quadrature values and weights, create a matrix which interacts
% all basis and test quadrature points up to the maximum degree.
% Gcoeffs will be of size
% [test_points, basis_points, max_deg_test+1, max_deg_basis+1]

% Multiply quadrature points up to max polynomial degree
Gcoeffs_test =  bsxfun(@power, so, max_deg_test -(0:max_deg_test));
Gcoeffs_basis = bsxfun(@power, si, max_deg_basis-(0:max_deg_basis));

% multiply with quadrature weights
Gcoeffs_test =  bsxfun(@times, wo, Gcoeffs_test);
Gcoeffs_basis = bsxfun(@times, wi, Gcoeffs_basis);

test_points =  numel(so);
basis_points = numel(si);

% Multiply every single basis and test coefficient
Gcoeffs = bsxfun(@times,reshape(Gcoeffs_test, [  test_points, 1,max_deg_test+1]),...
                        reshape(Gcoeffs_basis,[1,basis_points,1,max_deg_basis+1]));

end

