% BEUT - simple test 1
% Region 1 = UTLM cylinder mesh of free space
% Region 2 = free space BEM
% Excite with a plane wave and see the wave pass straight through
clear all;

% NOTE: this is the only test in +BEUT/+Main that allows you to choose whether to use Matlab to compute the
% BEM operators, the others will require using the C++ program because Matlab will be too slow!
use_Matlab_to_compute_operators = false;

global mu0 eps0;


%% Load the geometry (along with dt, mu0 and eps0)
filename = 'cyl_res21.mat';
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename]);
boundary=BEUT.Meshing.MeshBoundary(mesh);
N_V = boundary.N_V;

% Define material parameters
eps_r(1) = 1; mu_r(1) = 1;  % region 1 is background
eps_r(2) = 1; mu_r(2) = 1;  % region 2 is inside cylinder

for i=1:boundary.num_shapes+1
    mu(i) = mu_r(i)*mu0; eps(i) = eps_r(i)*eps0;
    eta(i) = sqrt(mu(i)/eps(i));
    c(i) = 1/sqrt(mu(i)*eps(i));
end
mesh.setMaterial(eps_r(2),mu_r(2));
mesh.calcAdmittances;

% Temporal parameters
N_T = 1000;
time = 0:dt:(N_T-1)*dt;


%% Set up excitation
% Gaussian pulse properties
edge_lengths = vertcat(boundary.halfedges.l);
width = min(edge_lengths)/c(1)*30;
delay = 1.5;
inc_wave = BEUT.Excitation.GaussianWave(width,delay,c(1));
inc_wave.direction = [1 0];
V_source = inc_wave.eval(time);
figure; plot(time,V_source)
title('Incident wave in the time domain'); xlabel('time');

% check stability
min_wavelength = c(1)/inc_wave.freq_response(time,true);
if min(edge_lengths)>min_wavelength/10
    warning(['Minimum edge length (' num2str(min(edge_lengths)) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end

% Plane wave from external BEM region
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,true);
dual_hat_function = BEUT.BEM.BasisFunction.createDualHat(boundary.dual,true);

rhsCalc = BEUT.BEM.RHS(N_T, dt);
rhsCalc.excitation = @inc_wave.eval;
rhsCalc.Gaussian_points = 3;
rhsCalc.polarization = [0 -1];
rhsCalc.display_plot = true;

rhsCalc.geometry = boundary.halfedges;
rhsCalc.test_function = square_function;
Ez_i  = rhsCalc.compute(false);
Hxy_i = rhsCalc.compute(true) / eta(1);


%% Probe positions
observation_edges = [32 253 70];
boundary_points = [1 3];
mesh.plot_halfedge(observation_edges)


%% Compute or load the operators
if use_Matlab_to_compute_operators
    
    % Compute operators
    outer_points = 25;
    inner_points = 26;
    outer_points_sp = 5*outer_points;
    inner_points_sp = 5*outer_points+1;
    
    Lagrange_degree = 1;
    cheat = true;
    
    Zmatrix = BEUT.BEM.ZMatrices(N_T,dt,boundary.dual,c(1));
    
    Zmatrix.basis_function_Z = BEUT.BEM.BasisFunction.createDualSquare(boundary.dual,false);
    Zmatrix.basis_function_S = BEUT.BEM.BasisFunction.createDualHat(boundary.dual,false);
    Zmatrix.test_function_Z = BEUT.BEM.BasisFunction.createDualSquare(boundary.dual,true);
    Zmatrix.test_function_S = BEUT.BEM.BasisFunction.createDualHat(boundary.dual,true);
    
    timeBasis = BEUT.BEM.LagrangeInterpolator(dt,Lagrange_degree);
    Zmatrix.timeBasis_D = timeBasis;
    Zmatrix.timeBasis_Nh = int(timeBasis);
    Zmatrix.timeBasis_Ns = diff(timeBasis);
    Zmatrix.outer_points_sp = outer_points_sp;
    Zmatrix.inner_points_sp = inner_points_sp;
    Zmatrix.outer_points = outer_points;
    Zmatrix.inner_points = inner_points;
    
    tic
    operator_file.N_T = N_T;
    [operator_file.S,operator_file.D,operator_file.Dp,Nh,Ns] = Zmatrix.compute(cheat);
    operator_file.N = Nh + Ns/c(1)^2;
    toc
    
else
    operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename]);
end

mesh = BEUT.Main.MOT(mesh, boundary, operator_file, observation_edges,...
        time, mu(1), Ez_i, Hxy_i);


%% Plot results against incident waves (analytic result)
tstop=size(mesh.fields.E_z,2);
figure; plot(time(1:tstop),mesh.fields.E_z(observation_edges,1:tstop))
BEM_observation(1) = find(mesh.mesh_boundary==observation_edges(boundary_points(1)));
BEM_observation(2) = find(mesh.mesh_boundary==observation_edges(boundary_points(2)));
hold all; plot(time(1:tstop),Ez_i(BEM_observation,1:tstop),':')
entries = cell(1,numel(observation_edges));
for i=1:numel(observation_edges)
    entries(i) = {sprintf('E_z at halfedge %i',observation_edges(i))};
end
entries(i+1) = {sprintf('E_z^i at halfedge %i',observation_edges(boundary_points(1)))};
entries(i+2) = {sprintf('E_z^i at halfedge %i',observation_edges(boundary_points(2)))};
legend('String',entries);
xlabel('time (s)'); ylabel('E_z')

