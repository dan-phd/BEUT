%% 2D time domain BEM for cylindrical dielectric
clear all;

%% Load the geometry
filename = 'cyl_res21';

load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename '.mat']);
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
N_T = 400;
f_max = c(1);
dt = 1/(8*f_max);
time = 0 : dt : N_T*dt-dt;
mesh.dt = dt;


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
Z = [-D{1}-D{2}         (S{1}+S{2})*mu(1);...
     (N{1}+N{2})/mu(1)  Dp{1}+Dp{2}];
rhs = [ Ez_i; Hxy_i];
unknown = BEUT.BEM.Main.MOT(Z, rhs);

% split m and j
M_TM = unknown(1:N_V,:);        % E_z = M
J_TM = -unknown(N_V+1:2*N_V,:); % n x H_xy = -J

%% Plot current density at shadow side and at exposed side of cylinder
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


%% Output file to compute scattered field in C++
% then run the C++ code using flags: -fcyl_res21_scattered -t400 -S
domain_X = [-3:0.1:3];
domain_Y = [-3:0.1:3];
[X,Y] = meshgrid(domain_X,domain_Y);
x_coords = X(:); y_coords = Y(:);
dual = false;
c_file = [BEUT.CFolder filesep 'input' filesep filename '_scattered.mat'];
in = BEUT.BEM.Main.saveScatteredFieldPoints(mesh,x_coords,y_coords, M_TM, J_TM,dual,c_file);

%% Once C++ has computed the scattered field operators, run this section to animate
%{
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename '_scattered.mat']);
E_s = BEUT.BEM.Main.organizeScatteredField(operator_file, X, in );

% animate
material_vertices = vertcat(boundary.halfedges.a);
BEUT.animate_fields(2,'domain',X,Y,...
    'animation',E_s,...
    'overlay',material_vertices,'dimensions',2);
%}
