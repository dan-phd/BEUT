% Compare operator matrices obtained by Matlab code with those obtain with the C++ code

filename = 'cyl_res21.mat';
N_T = 10;

%% Compute operator matrices using Matlab
load([BEUT.CFolder '\input\' filename]);

outer_points = 25;
inner_points = outer_points+1;
outer_points_sp = 5*outer_points;
inner_points_sp = 5*outer_points+1;

Lagrange_degree = 1;
cheat = false;

Zmatrix = BEUT.BEM.ZMatrices(N_T,dt,boundary,c);

if dual
    Zmatrix.basis_function_Z = BEUT.BEM.BasisFunction.createDualSquare(boundary,false);
    Zmatrix.basis_function_S = BEUT.BEM.BasisFunction.createDualHat(boundary,false);
    Zmatrix.test_function_Z = BEUT.BEM.BasisFunction.createDualSquare(boundary,true);
    Zmatrix.test_function_S = BEUT.BEM.BasisFunction.createDualHat(boundary,true);
else
    Zmatrix.basis_function_Z = BEUT.BEM.BasisFunction.createSquare(boundary,false);
    Zmatrix.basis_function_S = BEUT.BEM.BasisFunction.createHat(boundary,false);
    Zmatrix.test_function_Z = BEUT.BEM.BasisFunction.createSquare(boundary,true);
    Zmatrix.test_function_S = BEUT.BEM.BasisFunction.createHat(boundary,true);
end

timeBasis = BEUT.BEM.LagrangeInterpolator(dt,Lagrange_degree);
Zmatrix.timeBasis_D = timeBasis;
Zmatrix.timeBasis_Nh = int(timeBasis);
Zmatrix.timeBasis_Ns = diff(timeBasis);
Zmatrix.outer_points_sp = outer_points_sp;
Zmatrix.inner_points_sp = inner_points_sp;
Zmatrix.outer_points = outer_points;
Zmatrix.inner_points = inner_points;

tic
[S1,D1,Dp1,Nh1,Ns1] = Zmatrix.compute(cheat);
N1 = Nh1 + Ns1/c^2;
toc


%% Compare with operator matrices obtained using C++
cfile = matfile([BEUT.CFolder '\results\' filename]);
assert(N_T<=cfile.N_T,['Chosen N_T must be less than or equal to the number of timesteps'...
    'computed in the file']);
N2 = cfile.N(:,:,1:N_T);
S2 = cfile.S(:,:,1:N_T);
D2 = cfile.D(:,:,1:N_T);
Dp2 = cfile.Dp(:,:,1:N_T);
BEUT.relError(N1(:,:,1:N_T),N2(:,:,1:N_T));
BEUT.relError(S1(:,:,1:N_T),S2(:,:,1:N_T));
BEUT.relError(D1(:,:,1:N_T),D2(:,:,1:N_T));
BEUT.relError(Dp1(:,:,1:N_T),Dp2(:,:,1:N_T));
