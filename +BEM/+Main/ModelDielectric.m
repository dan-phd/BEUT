%% 2D time domain BEM for cylindrical dielectric
clear all

%% Parameters
mu0 = 4e-7*pi;               % permeability of free space
eps0 = 8.854187817e-12;      % permittivity of free space
c0 = 1/sqrt(mu0*eps0);       % propagation speed inside vacuum


%% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);
radius=max(range(vertcat(boundary.halfedges.a)))/2;
N_V = boundary.N_V;
boundary.halfedges = BEUT.BEM.Analytical.reorderFacesToMatchAnalytical(boundary.halfedges);
shadow_side = 1;
exposed_side = floor(N_V/2+1);

% Define materials
eps_r(1) = 1; mu_r(1) = 1;  % region 1 is background
eps_r(2) = 2; mu_r(2) = 1;  % region 2 is inside cylinder

for i=1:boundary.num_shapes+1
    mu(i) = mu_r(i)*mu0; eps(i) = eps_r(i)*eps0;
    eta(i) = sqrt(mu(i)/eps(i));
    c(i) = 1/sqrt(mu(i)*eps(i));
end

% Temporal parameters
N_T = 300;
f_max = c0;
dt = 1/(8*f_max);
time = 0 : dt : N_T*dt-dt;


%% RHS setup
T = 5/f_max;           % width of pulse
t0 = 1.5;              % start time (relative to width)
direction = [1 0];

% Create excitation in region 1 and note cutoff frequency < max frequency
gpw = BEUT.Excitation.GaussianWave(T, t0, c(1), direction);
figure; plot(time,gpw.eval(time));
Fc = gpw.freq_response(time)

rhsCalc = BEUT.BEM.RHS(N_T, dt);
rhsCalc.excitation = @gpw.eval;
rhsCalc.Gaussian_points = 3;
rhsCalc.polarization = [0 -1];
rhsCalc.display_plot = false;
rhsCalc.geometry = boundary.halfedges;


%% Computation of operators

% Create basis functions
hat_function = BEUT.BEM.BasisFunction.createHat(boundary.halfedges,false);
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,true);

% Gaussian quadrature points
gauss_points_m = 4;
gauss_points_n = 5;
outer_points_sp = 50;
inner_points_sp = 51;
degree = 1;

% Compute Z matrices
zm = BEUT.BEM.ZMatrices(N_T,dt,boundary.halfedges);
    
% Basis and test functions applied to each segment
% S = transverse plane, Z = z-directed
zm.basis_function_Z = square_function;
zm.basis_function_S = hat_function;
zm.test_function_Z = square_function;
zm.test_function_S = hat_function;

% Lagrange interpolators (temporal basis functions)
timeBasis = BEUT.BEM.LagrangeInterpolator(dt,degree);
zm.timeBasis_D = timeBasis;
zm.timeBasis_Nh = int(timeBasis);
zm.timeBasis_Ns = diff(timeBasis);

zm.outer_points = gauss_points_m;
zm.inner_points = gauss_points_n;
zm.outer_points_sp = outer_points_sp;
zm.inner_points_sp = inner_points_sp;

% we can use the cheat when the geometry is a cylinder with equal edge lengths
cheat = true;
    
tic
for i=1:boundary.num_shapes+1
    
    zm.c = c(i);
    
    [S{i},D{i},Dp{i},Nh,Ns] = zm.compute(cheat);
    N{i} = Nh + Ns/c(i)^2;
    
end
toc


%% TM
% RHS vector calculation
rhsCalc.test_function = square_function;
Ez_i  = rhsCalc.compute(false);
rhsCalc.test_function = hat_function;
Hxy_i = rhsCalc.compute(true) / eta(1);

% MOT
Z = [-D{1}-D{2}       (S{1}+S{2})*mu(1);...
     (N{1}+N{2})/mu(1)  Dp{1}+Dp{2}];
rhs = [ Ez_i; Hxy_i];
unknown = BEUT.BEM.Main.MOT(Z, rhs);

% split m and j
M_TM = unknown(1:N_V,:);        % E_z = M
J_TM = -unknown(N_V+1:2*N_V,:); % n x H_xy = -J

% Plot current density at shadow side and at exposed side of cylinder
points_to_plot(1).name='Exposed side'; points_to_plot(1).point=exposed_side;
points_to_plot(2).name='Shadow side'; points_to_plot(2).point=shadow_side;
BEUT.BEM.Main.plotCurrentDensityInTime(time, Ez_i, points_to_plot); title('Ez_i');
BEUT.BEM.Main.plotCurrentDensityInTime(time, M_TM, points_to_plot); title('M');
BEUT.BEM.Main.plotCurrentDensityInTime(time, Hxy_i, points_to_plot); title('Hxy_i');
BEUT.BEM.Main.plotCurrentDensityInTime(time, J_TM, points_to_plot); title('J');

% Calculate current density in frequency domain at all points on circle
[M_FFT,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity( time,M_TM,gpw,c(1) );

% Analytic solution
analytic = BEUT.BEM.Analytical.AnalyticalDielectricCylinder(N_V,radius,omega,mu0,eps0);
analytic.eps_r=eps_r(2); analytic.mu_r=mu_r(2);
[analytic_M_FFT,~] = analytic.calcSurfaceCurrents;

% Error between the analytical and numerical solutions
BEUT.relError(abs(M_FFT(:,2:limit)),abs(analytic_M_FFT(:,2:limit)) );

% Frequency domain plots
probe_freq_idx = 5;
BEUT.BEM.Main.plotCurrentDensityInFrequency(c(1), M_FFT, omega, A, limit, points_to_plot, analytic_M_FFT)
BEUT.BEM.Main.plotCurrentDensityAtOneFrequency( N_V, M_FFT, omega, probe_freq_idx, analytic_M_FFT )


%% TE
%{
% RHS vector calculation
rhsCalc.test_function = square_function;
Hz_i = rhsCalc.compute(false) / eta(1);
rhsCalc.test_function = hat_function;
Exy_i  = rhsCalc.compute(true);

% MOT
Z = [-D{1}-D{2}         -(S{1}+S{2})*eps(1);...
     -(N{1}+N{2})/eps(1)  Dp{1}+Dp{2}];
rhs_vec = [ Hz_i; Exy_i ];
unknown = BEUT.BEM.Main.MOT(Z, rhs_vec);

% split m and j
J_TE = unknown(1:N_V,:);        % H_z
M_TE = -unknown(N_V+1:2*N_V,:); % n x E_xy

%% Plot current density at shadow side and at exposed side of cylinder
BEUT.BEM.Main.plotCurrentDensityInTime(time, Exy_i, points_to_plot); title('Exy_i');
BEUT.BEM.Main.plotCurrentDensityInTime(time, M_TE, points_to_plot); title('M');
BEUT.BEM.Main.plotCurrentDensityInTime(time, Hz_i, points_to_plot); title('Hz_i');
BEUT.BEM.Main.plotCurrentDensityInTime(time, J_TE, points_to_plot); title('J');

% Calculate current density in frequency domain at all points on circle
[J_FFT,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity( time,J_TE,gpw,c(1) );

% Error between the analytical and numerical solutions
analytic = BEUT.BEM.Analytical.AnalyticalDielectricCylinder(N_V,radius,omega,mu0,eps0);
analytic.eps_r=eps_r(2); analytic.mu_r=mu_r(2);
[~,analytic_J_f] = analytic.calcSurfaceCurrents;
analytic_J_f = analytic_J_f/eta(1);
BEUT.relError(abs(J_FFT(:,2:limit)),abs(analytic_J_f(:,2:limit)) );

% Frequency domain plots
probe_freq_idx = 2;
BEUT.BEM.Main.plotCurrentDensityInFrequency(c(1), J_FFT, omega, A, limit, points_to_plot, analytic_J_f)
BEUT.BEM.Main.plotCurrentDensityAtOneFrequency( N_V, J_FFT, omega, probe_freq_idx, analytic_J_f )

%}
