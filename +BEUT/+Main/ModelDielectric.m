% BEUT - model dielectric cylinder
% Region 1 = meshed free space background
% Region 2 = dielectric cylinder
% Excite with a plane wave using BEUT and UTLM, and compare numerical results with analytical solution
clear all;
global mu0 eps0;


%% Load the geometry (along with dt, mu0 and eps0)
% choose 'squ_cyl_res40.mat' for a coarse mesh
% or 'squ_cyl_res80.mat' for a finer mesh
filename = 'squ_cyl_res40.mat';

load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename]);
boundary=BEUT.Meshing.MeshBoundary(mesh);

% Define material layers
mu_r(1)=1; eps_r(1)=1;      % region 1 is meshed free space background
mu_r(2)=1; eps_r(2)=3;      % inner cylinder

for i=1:boundary.num_shapes+1
    mu(i) = mu_r(i)*mu0; eps(i) = eps_r(i)*eps0;
    eta(i) = sqrt(mu(i)/eps(i));
    c(i) = 1/sqrt(mu(i)*eps(i));
end
mesh.setMaterial(eps_r,mu_r);
mesh.calcAdmittances;

% Temporal parameters
N_T = 8000;
time = 0:dt:(N_T-1)*dt;


%% Set up excitation
% Gaussian pulse properties
min_edge_length = min(vertcat(boundary.halfedges.l));
width = min_edge_length/c(1)*30;
delay = 1.5;
inc_wave = BEUT.Excitation.GaussianWave(width,delay,c(1));
inc_wave.direction = [1 0];
V_source = inc_wave.eval(time);
figure; plot(time,V_source);
title('Incident wave in the time domain');

% check stability
min_wavelength = c(1)/inc_wave.freq_response(time,false);
if min_edge_length>min_wavelength/10
    warning(['Minimum edge length (' num2str(min_edge_length) ...
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


%% set boundary condition: 1(=open), -1(=short), 0(=absorbing)
if strcmp(filename,'squ_cyl_res40.mat')
    bottom_side = 1:8; right_side = 9:16; top_side = 17:24; left_side = 25:32;
elseif strcmp(filename,'squ_cyl_res80.mat')
    bottom_side = 1:20; right_side = 21:42; top_side = 43:62; left_side = 63:84;
end

condition([bottom_side,top_side]) = 1;
condition([left_side,right_side]) = 0;
mesh.setBoundary(condition);


%% Probe positions (exposed and shadow side)
if strcmp(filename,'squ_cyl_res40.mat')
    observation_edges = [316 493];
elseif strcmp(filename,'squ_cyl_res80.mat')
    observation_edges = [1811 713];
end


%% BEUT MOT
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename]);
mesh = BEUT.Main.MOT(mesh, boundary, operator_file, observation_edges,...
    time, mu(1), Ez_i, Hxy_i);
Ez_BEUT = mesh.fields.E_z;


%% Use UTLM to run the same tests
source_edges=mesh.mesh_boundary(left_side);
mesh = BEUT.UTLM.Main.run(mesh, N_T, V_source, source_edges, 750, 0);
tstop=size(mesh.fields.E_z,2);
window=ones(1,tstop);
w = 0.5*(1-cos(2*pi*(1:tstop)/(tstop/2)));
window(ceil(3*tstop/4):tstop)=w(ceil(3*tstop/4):tstop);
Ez_UTLM = bsxfun(@times,mesh.fields.E_z,window);


%% Time domain plots
max_BEUT = max(max(Ez_BEUT));
Ez_UTLM_resized = Ez_UTLM*max_BEUT/max(max(Ez_UTLM));
figure; hold on;
plot(time(1:tstop),Ez_BEUT(observation_edges,1:tstop))
plot(time(1:tstop),Ez_UTLM_resized(observation_edges,1:tstop),'--')
entries = cell(1,numel(observation_edges));
for i=1:numel(observation_edges)
    entries(i) = {sprintf('BEUT - E_z at halfedge %i',observation_edges(i))};
    entries(numel(observation_edges)+i) = {sprintf('UTLM - E_z at halfedge %i',observation_edges(i))};
end
legend('String',entries);


%% Compare analytic and numerical J
% Calculate J in frequency domain at all points on the internal circle
cylinder_halfedges = mesh.material_boundaries{2};
radius = 0.35;
N_V_cylinder = numel(cylinder_halfedges);

[M_FFT_UTLM,omega,A,limit] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity...
    ( time,Ez_UTLM_resized(cylinder_halfedges,:),inc_wave,c(1) );
M_FFT_UTLM = M_FFT_UTLM/M_FFT_UTLM(1);
[M_FFT_BEUT] = BEUT.BEM.Main.findFrequencyDomainCurrentDensity...
    ( time,Ez_BEUT(cylinder_halfedges,:),inc_wave,c(1) );

analytic = BEUT.BEM.Analytical.AnalyticalDielectricCylinder(N_V_cylinder,radius,omega,mu0,eps0);
analytic.eps_r=eps_r(2); analytic.mu_r=mu_r(2);
[analytic_M_FFT,~] = analytic.calcSurfaceCurrents;

exposed_side = find(cylinder_halfedges==observation_edges(1));
shadow_side = find(cylinder_halfedges==observation_edges(2));
analytic_M_FFT_reordered([exposed_side shadow_side],:) = ...
    analytic_M_FFT([floor(N_V_cylinder/2+1) 1],:);
BEUT.BEM.Main.plotCurrentDensityInFrequency...
    (c(1), M_FFT_BEUT, omega, A, limit, exposed_side,shadow_side, analytic_M_FFT_reordered)
BEUT.BEM.Main.plotCurrentDensityInFrequency...
    (c(1), M_FFT_UTLM, omega, A, limit, exposed_side,shadow_side, analytic_M_FFT_reordered)
title('Current density at exposed and shadow side of shape (using UTLM)');

