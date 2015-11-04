%% 2D time domain BEM for 2 separated dielectrics
clear all

%% Parameters
mu0 = 4e-7*pi;               % permeability of free space
eps0 = 8.854187817e-12;      % permittivity of free space
c0 = 1/sqrt(mu0*eps0);       % propagation speed inside vacuum


%% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21rec_res42_combined.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);
N_V = boundary.N_V;

% Define materials
eps_r(1) = 1; mu_r(1) = 1;  % region 1 is 1st shape
eps_r(2) = 1; mu_r(2) = 1;  % region 2 is 2nd shape
eps_r(3) = 1; mu_r(3) = 1;  % region 3 is background (free space)

for i=1:boundary.num_shapes+1
    mu(i) = mu_r(i)*mu0; eps(i) = eps_r(i)*eps0;
    eta(i) = sqrt(mu(i)/eps(i));
    c(i) = 1/sqrt(mu(i)*eps(i));
end

% Define regions (third region is the combination of regions 1 and 2)
region{1} = vertcat(boundary.halfedges.shape)==1;
region{2} = vertcat(boundary.halfedges.shape)==2;
region{3} = region{1} | region{2};
cheat(1) = true;    % only the 1st shape is a cylinder (so we can use the cheat)
cheat(2) = false;
cheat(3) = false;

% Temporal parameters
N_T = 300;
f_max = c0;
dt = 1/(8*f_max);
time = 0 : dt : N_T*dt-dt;


%% RHS setup
%
T = 5/f_max;           % width of pulse
t0 = 1.5;              % start time (relative to width)
direction = [1 0];

% Create excitation in each region
for i=1:boundary.num_shapes+1
    gpw(i) = BEUT.Excitation.GaussianWave(T, t0, c(i), direction);
%     figure; plot(time,gpw(i).eval(time));
%     Fc = gpw(i).freq_response(time)
end

rhsCalc(1) = BEUT.BEM.RHS(N_T, dt);
rhsCalc(1).Gaussian_points = 3;
rhsCalc(1).polarization = [0 -1];
rhsCalc(1).display_plot = false;

for i=1:boundary.num_shapes+1
    rhsCalc(i) = rhsCalc(1);
    gpw_=gpw(i);
    rhsCalc(i).excitation = @gpw_.eval;
    rhsCalc(i).geometry = boundary.halfedges(region{i});
end


%% Computation of operators

% Create basis functions
for i=1:boundary.num_shapes+1
    hat_function(i) = BEUT.BEM.BasisFunction.createHat(boundary.halfedges(region{i}),false);
    square_function(i) = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges(region{i}),true);
end

% Gaussian quadrature points and temporal degree
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
    
tic
for i=1:boundary.num_shapes+1
    
    zm.c = c(i);
    
    [S{i},D{i},Dp{i},Nh,Ns] = zm.compute(cheat);
    N{i} = Nh + Ns/c(i)^2;
    
end
toc


%% Compute Gram matrices
gram = BEUT.BEM.GramMatrix();
for i=1:boundary.num_shapes+1
    gram.geometry = boundary.halfedges(region{i});
    
    gram.test_function = square_function(i);
    gram.basis_function = hat_function(i);
    G1{i} = zeros(size(D{i})); G1{i}(:,:,1) = compute(gram);
    
    gram.test_function = hat_function(i);
    gram.basis_function = square_function(i);
    G2{i} = zeros(size(D{i})); G2{i}(:,:,1) = compute(gram);
end


%% TM
% RHS vector calculation
for i=1:boundary.num_shapes+1
    
    rhsCalc(i).test_function = square_function(i);
    Ez_i{i}  = rhsCalc(i).compute(false);
    
    rhsCalc(i).test_function = hat_function(i);
    Hxy_i{i} = rhsCalc(i).compute(true) / eta(i);
    
end

% no incident field in regions 1 and 2
for i=[1 2]
    Ez_i{i}=zeros(size(Ez_i{i}));
    Hxy_i{i}=zeros(size(Hxy_i{i}));
end

% MOT
O{1} = zeros(size(D{2},1),size(D{1},2),N_T);
O{2} = zeros(size(D{1},1),size(D{2},2),N_T);

Z3 = [G1{3}+D{3}   -S{3}*mu(3);...
      -N{3}/mu(3)  G2{3}-Dp{3}];

Z12 = [ G1{1}-D{1} O{2}       S{1}*mu(1)  O{2}       ;...
        O{1}       G1{2}-D{2} O{1}        S{2}*mu(2)  ;...
        N{1}/mu(1) O{2}       G2{1}+Dp{1} O{2}       ;...
        O{1}       N{2}/mu(2) O{1}        G2{2}+Dp{2} ];

Z = Z12-Z3;

rhs = [ Ez_i{3}(region{1},:) - Ez_i{1};...
        Ez_i{3}(region{2},:) - Ez_i{2};...
        Hxy_i{3}(region{1},:) - Hxy_i{1};...
        Hxy_i{3}(region{2},:) - Hxy_i{2}];
unknown = BEUT.BEM.Main.MOT(Z, rhs);

% split m and j
M = unknown(1:N_V,:);        % E_z = M
J = -unknown(N_V+1:2*N_V,:); % n x H_xy = -J


%% Plot current density at various points in the mesh
positions = 1:10:N_V;
boundary.plot(positions)
for i=1:numel(positions), points_to_plot(i).point=positions(i); end
time_to_plot=time(1:100);
BEUT.BEM.Main.plotCurrentDensityInTime(time_to_plot, Ez_i{3}, points_to_plot); title('Ez_i');
BEUT.BEM.Main.plotCurrentDensityInTime(time_to_plot, M, points_to_plot); title('M');
BEUT.BEM.Main.plotCurrentDensityInTime(time_to_plot, Hxy_i{3}, points_to_plot); title('Hxy_i');
BEUT.BEM.Main.plotCurrentDensityInTime(time_to_plot, J, points_to_plot); title('J');


