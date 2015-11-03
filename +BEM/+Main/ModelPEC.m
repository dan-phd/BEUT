%% 2D time domain BEM for cylindrical PEC
clear all

%% Parameters
mu = 4*pi*10^-7;            % permeability of free space
eps = 8.854187817e-12;      % permittivity of free space
c = 1/sqrt(mu*eps);         % propagation speed inside vacuum
eta = sqrt(mu/eps);


%% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);
radius=max(range(vertcat(boundary.halfedges.a)))/2;
N_V = boundary.N_V;
boundary.halfedges = BEUT.BEM.Analytical.reorderFacesToMatchAnalytical(boundary.halfedges);
shadow_side = 1;
exposed_side = floor(N_V/2+1);
observation_points(1).name='Exposed side'; observation_points(1).point=exposed_side;
observation_points(2).name='Shadow side'; observation_points(2).point=shadow_side;

% Temporal parameters
N_T = 400;
omega_max = N_V*c/(10*radius);
dt = pi/(7*omega_max);               % oversampling factor = 7 here
time = 0 : dt : N_T*dt-dt;           % time vector used for graphs


%% RHS setup
T = 8 / c;            % width of pulse
t0 = 1.5;             % start time (relative to width)
direction = [1 0];

gpw = BEUT.Excitation.GaussianWave(T, t0, c, direction);
figure; plot(time,gpw.eval(time));
Fc = gpw.freq_response(time)

rhsCalc = BEUT.BEM.RHS(N_T, dt);
rhsCalc.excitation = @gpw.eval;
rhsCalc.Gaussian_points = 3;
rhsCalc.polarization = [0 -1];
rhsCalc.display_plot = false;
rhsCalc.geometry = boundary.halfedges;


%% Computation
tic

% Create basis functions, 2nd argument false means function has magnitude 1
hat_function = BEUT.BEM.BasisFunction.createHat(boundary.halfedges,false);
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,false);

% Compute the Z matrices
gauss_points_m = 4;
gauss_points_n = 5;
Zmatrix_obj = BEUT.BEM.ZMatrices(N_T,dt,boundary.halfedges,c);

% Basis and test functions applied to each segment
% S = transverse plane, Z = z-directed, d = S divergence
% TE EFIE and TM MFIE require hat functions, so lets just use
% all hat functions so we only have to compute the operators once
Zmatrix_obj.basis_function_Z = hat_function;
Zmatrix_obj.basis_function_S = hat_function;
Zmatrix_obj.test_function_Z = hat_function;
Zmatrix_obj.test_function_S = hat_function;

% Lagrange interpolators (temporal basis functions)
timeBasis = BEUT.BEM.LagrangeInterpolator(dt,1);
Zmatrix_obj.timeBasis_D = timeBasis;
Zmatrix_obj.timeBasis_Nh = int(timeBasis);
Zmatrix_obj.timeBasis_Ns = diff(timeBasis);

Zmatrix_obj.outer_points = gauss_points_m;
Zmatrix_obj.inner_points = gauss_points_n;
[S,D,Dp,Nh,Ns] = Zmatrix_obj.compute(true);
O = S*0;        % zero matrix
toc


%% TE MFIE
% RHS vector calculation
rhsCalc.test_function = square_function;
Hz_i  = rhsCalc.compute(false) / eta;

% Compute Gram matrix
Gram = BEUT.BEM.GramMatrix();
Gram.geometry = boundary.halfedges;
Gram.test_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,false);
g = compute(Gram);
G = O; G(:,:,1) = 0.5*g;

% MOT
Z = -G+D;
unknown = -BEUT.BEM.Main.MOT(Z, Hz_i);
Jxy_MFIE = -unknown * eta;

% Plot 'J' at shadow side and at exposed side of cylinder
BEUT.BEM.Main.plotCurrentDensityInTime(time, Jxy_MFIE, observation_points); title('MFIE J_x_y');

% Calculate current density in frequency domain at all points on circle
[Jfxy_MFIE,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity( time,Jxy_MFIE,gpw,c );

% Analytical solution
A_obj = BEUT.BEM.Analytical.AnalyticalPECCylinder(N_V,radius,omega);
A_obj.mu = mu; A_obj.eps = eps;
Analytical_Jz  = A_obj.calcTM_J;
Analytical_Jxy = A_obj.calcTE_J;

% Error between the analytical and numerical solutions for all stable frequencies and all theta
% (first element will be inf, so start at 2nd element)
BEUT.relError(abs(Jfxy_MFIE(:,2:limit)),abs(Analytical_Jxy(:,2:limit)) );

% Plot J(w)
BEUT.BEM.Main.plotCurrentDensityInFrequency(c, Jfxy_MFIE, omega, A, limit, observation_points, Analytical_Jxy)


%% TM MFIE
% RHS vector calculation
rhsCalc.test_function = hat_function;
Hxy_i  = rhsCalc.compute(true) / eta;

Gram.test_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,false);
Gram.basis_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,true);
g = Gram.compute;
G = O; G(:,:,1) = 0.5*g;

% MOT
Z = G+Dp;
unknown = -BEUT.BEM.Main.MOT(Z, Hxy_i);
Jz_MFIE = unknown * eta;

% Differentiate Jz so FFT computes on function tending to zero
Jz_MFIE = diff(Jz_MFIE,1,2)/dt;
Jz_MFIE(:,N_T)=zeros(1,N_V);

% Plot 'J' at shadow side and at exposed side of cylinder
BEUT.BEM.Main.plotCurrentDensityInTime(time, Jz_MFIE, observation_points); title('MFIE J_z');

% Calculate current density in frequency domain at all points on circle
[Jfz_MFIE,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity( time,Jz_MFIE,gpw,c );

% Integrate Jz to bring back to correct value
Jfz_MFIE = bsxfun(@rdivide,Jfz_MFIE , 1i*omega );

% Error between the analytical and numerical solutions for all stable frequencies and all theta
% (first element will be inf, so start at 2nd element)
BEUT.relError(abs(Jfz_MFIE(:,2:limit)),abs(Analytical_Jz(:,2:limit)) );

% Plot J(w)
BEUT.BEM.Main.plotCurrentDensityInFrequency(c, Jfz_MFIE, omega, A, limit, observation_points, Analytical_Jz)


%% TM EFIE
% RHS vector calculation
rhsCalc.test_function = square_function;
Ez_i  = rhsCalc.compute(false);

% MOT
Z = S*mu;
unknown = BEUT.BEM.Main.MOT(Z, Ez_i);
Jz_EFIE = unknown * eta;

% Differentiate Jz so FFT computes on function tending to zero
Jz_EFIE = diff(Jz_EFIE,1,2)/dt;
Jz_EFIE(:,N_T)=zeros(1,N_V);

% Plot 'J' at shadow side and at exposed side of cylinder
BEUT.BEM.Main.plotCurrentDensityInTime(time, Jz_EFIE, observation_points); title('MFIE J_z');

% Calculate current density in frequency domain at all points on circle
[Jfz_EFIE,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity( time,Jz_EFIE,gpw,c );
% Integrate Jz to bring back to correct value
Jfz_EFIE = bsxfun(@rdivide,Jfz_EFIE , 1i*omega );

% Error between the analytical and numerical solutions for all stable frequencies and all theta
% (first element will be inf, so start at 2nd element)
BEUT.relError(abs(Jfz_EFIE(:,2:limit)),abs(Analytical_Jz(:,2:limit)) );

% Plot J(w)
BEUT.BEM.Main.plotCurrentDensityInFrequency(c, Jfz_EFIE, omega, A, limit, observation_points, Analytical_Jz)


%% TE EFIE
% RHS vector calculation
rhsCalc.test_function = hat_function;
Exy_i = rhsCalc.compute(true);

% MOT
Z = (Nh*c+Ns/c)*eta;
unknown = BEUT.BEM.Main.MOT(Z, Exy_i);
Jxy_EFIE = unknown * eta;

% Plot 'J' at shadow side and at exposed side of cylinder
BEUT.BEM.Main.plotCurrentDensityInTime(time, Jxy_EFIE, observation_points); title('EFIE J_x_y');

% Calculate current density in frequency domain at all points on circle
[Jfxy_EFIE,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity( time,Jxy_EFIE,gpw,c );

% Error between the analytical and numerical solutions for all stable frequencies and all theta
% (first element will be inf, so start at 2nd element)
BEUT.relError(abs(Jfxy_EFIE(:,2:limit)),abs(Analytical_Jxy(:,2:limit)) );

% Plot J(w)
BEUT.BEM.Main.plotCurrentDensityInFrequency(c, Jfxy_EFIE, omega, A, limit, observation_points, Analytical_Jxy)
