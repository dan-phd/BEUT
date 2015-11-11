% BEUT - test 5
% Region 1 = Luneburg lens 1
% Region 2 = Luneburg lens 2
% Region 3 = free space
% Excite with a plane wave and watch animation!
global mu0 eps0;


%% Load the geometry (along with dt, mu0 and eps0)
filename = 'cyl_res69.mat';
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename]);
boundary=BEUT.Meshing.MeshBoundary(mesh);
radius=max(range(vertcat(boundary.halfedges.a)))/2;

% Define material layers
% Luneburg lens has sn eps_r that is a function of the distance from the center
CC = circumcenter(mesh.TR);
for i=1:mesh.nF
    r = norm(CC(i,:));
    eps_r(i) = (2-(r/radius).^2);
end
mu_r(1)=1;
% mu_r(2)=1;

eta = sqrt(mu0/eps0);
c = 1/sqrt(mu0*eps0);

mesh.setMaterial(eps_r,mu_r);
mesh.calcAdmittances;
mesh.plot_materials('eps_r')

% Temporal parameters
N_T = 3000;
time = 0:dt:(N_T-1)*dt;


%% Set up excitation
% Sinewave signal properties
min_edge_length = min(vertcat(boundary.halfedges.l));
f_width = 0.3 * c;
f_mod = 0.9 * c;
pulse_width = 4/f_width;    % this is the corresponding width of the time domain pulse
t_end = 2*pulse_width;      % this is the time the pulse will sufficiently subside
direction = [1 0];
inc_wave = BEUT.Excitation.SineWave(f_width, f_mod, c, direction, 0.8);
V_source = inc_wave.eval(time);
figure; plot(time,V_source)

% check stability
min_wavelength = c/inc_wave.freq_response(time,false);
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
Hxy_i = rhsCalc.compute(true) / eta;


%% Probe positions (exposed and shadow side)
observation_edges = [117 1048];
BEM_points = [1];


%% BEUT MOT
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename]);
mesh = BEUT.Main.MOT(mesh, boundary, operator_file, observation_edges,...
    time, mu0(1), Ez_i, Hxy_i);


%% Time domain plots
tstop = size(mesh.fields.E_z,2);
BEM_observation(1) = find(mesh.mesh_boundary==observation_edges(BEM_points(1)));
figure; plot(time(1:tstop),mesh.fields.E_z(observation_edges,1:tstop))
hold all; plot(time(1:tstop),Ez_i(BEM_observation,1:tstop),':')
entries = cell(1,numel(observation_edges));
for i=1:numel(observation_edges)
    entries(i) = {sprintf('E_z at halfedge %i',observation_edges(i))};
end
entries(i+1) = {sprintf('E_z^i at halfedge %i',observation_edges(BEM_points(1)))};
legend('String',entries);

