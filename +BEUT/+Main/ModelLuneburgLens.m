% BEUT - Luneburg lens simulation
% Region 1 = Luneburg lens
% Excite with a plane wave or point source and watch animation!
clear all;
global mu0 eps0;


%% Load the geometry (along with dt, mu0 and eps0)
filename = 'cyl_res69';
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename '.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);
radius=max(range(vertcat(boundary.halfedges.a)))/2;

% Define material layers
% Luneburg lens has an eps_r that is a function of the distance from the center
CC = circumcenter(mesh.TR);
for i=1:mesh.nF
    r = norm(CC(i,:));
    eps_r(i) = (2-(r/radius).^2);
end
mu_r=1;

eta0 = sqrt(mu0/eps0);
c0 = 1/sqrt(mu0*eps0);

mesh.setMaterial(eps_r,mu_r);
mesh.calcAdmittances;
mesh.plot_materials('eps_r')

% Temporal parameters
N_T = 3000;
time = 0:dt:(N_T-1)*dt;


%% Set up excitation
% Sinewave signal properties
min_edge_length = min(vertcat(boundary.halfedges.l));
f_width = 0.3 * c0;
f_mod = 0.9 * c0;
direction = [1 0];
inc_wave = BEUT.Excitation.SineWave(f_width, f_mod, c0, direction, 0.8);
V_source = inc_wave.eval(time);
figure; plot(time,V_source)
title('Incident wave in the time domain'); xlabel('time');

% check stability
min_wavelength = c0/inc_wave.freq_response(time,false);
if min_edge_length>min_wavelength/10
    warning(['Minimum edge length (' num2str(min_edge_length) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end


%% Probe positions (exposed and shadow side)
observation_edges = [117 1048];
mesh.plot_halfedge(observation_edges);


%% BEUT MOT
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename '.mat']);
[mesh, M_TM, J_TM] = BEUT.Main.MOT(mesh, boundary, operator_file, observation_edges,...
    time, mu0, 0, 0, V_source, observation_edges(1));


%% Time domain plots
tstop = size(mesh.fields.E_z,2);
figure; plot(time(1:tstop),mesh.fields.E_z(observation_edges,1:tstop))
entries = cell(1,numel(observation_edges));
for i=1:numel(observation_edges)
    entries(i) = {sprintf('E_z at halfedge %i',observation_edges(i))};
end
legend('String',entries);


%% Animate field inside the UTLM region
mesh.animate('E')


%% Output file to compute scattered field in C++
% then run the C++ code using flags: -fcyl_res69_scattered -t3000 -S
[X,Y] = meshgrid([-1.5:0.2:2],[-1.5:0.2:1.5]);
x_coords = X(:); y_coords = Y(:);
M = mesh.fields.E_z(mesh.mesh_boundary,:);
J = -mesh.fields.H_xy(mesh.mesh_boundary,:);
dual = true;
c_file = [BEUT.CFolder filesep 'input' filesep filename '_scattered.mat'];
in_scatterer = BEUT.BEM.Main.saveScatteredFieldPoints(mesh,x_coords,y_coords,M,J,dual,c_file);


%% Run this AFTER scattered fields have been computed in C++ to plot all fields at a timestep
BEUT.Main.plotFields( filename,mesh,X,Y,in_scatterer,1000,true )


%% Run this AFTER scattered fields have been computed in C++ to animate the BEM scattered fields
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename '_scattered.mat']);
E_s = BEUT.BEM.Main.organizeScatteredField(operator_file, X, in_scatterer );
material_vertices = vertcat(boundary.halfedges.a);
BEUT.animate_fields(2,'domain',X,Y,...
    'animation',E_s/max(max(max(E_s))),...
    'overlay',material_vertices,'dimensions',2,...
    'skipTimesteps',10,...
    'max_amplitude',1,'min_amplitude',-1);

