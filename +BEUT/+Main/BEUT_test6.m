% BEUT - test 6
% Region 1 = Luneburg lens 1
% Region 2 = Luneburg lens 2
% Region 3 = free space
% Excite left lens with a point source, and observe right lens for transmitted signal
global mu0 eps0;


%% Load the geometry (along with dt, mu0 and eps0)
filename = 'cyl_res138_x2.mat';
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename]);
boundary=BEUT.Meshing.MeshBoundary(mesh);
radius=max(range(vertcat(boundary.halfedges(vertcat(boundary.halfedges.shape)==1).a)))/2;

% Define material layers
% Luneburg lens has an eps_r that is a function of the distance from the center
CC = circumcenter(mesh.TR);
eps_r = zeros(mesh.nF,1);
for i=1:mesh.nF
    if mesh.faces(i).fnum==1
        r = norm(CC(i,:));
        eps_r(i) = (2-(r/radius).^2);
    elseif mesh.faces(i).fnum==2
        r = norm(CC(i,:)-[5,0]);
        eps_r(i) = (2-(r/radius).^2);
    end
end
mu_r=1;

eta = sqrt(mu0/eps0);
c = 1/sqrt(mu0*eps0);

mesh.setMaterial(eps_r,mu_r);
mesh.calcAdmittances;
mesh.plot_materials('eps_r')

% Temporal parameters
N_T = 5000;
time = 0:dt:(N_T-1)*dt;


%% Set up excitation
% Sinewave signal properties
min_edge_length = min(vertcat(boundary.halfedges.l));
f_mod = 2 * c;          % 2 wavelengths per radius
f_width = f_mod/8;
direction = [1 0];
inc_wave = BEUT.Excitation.SineWave(f_width, f_mod, c, direction, 0.65);
inc_wave.envelope = @cosine;
V_source = inc_wave.eval(time);
figure; plot(time,V_source)

% check stability
min_wavelength = c/inc_wave.freq_response(time,false);
if min_edge_length>min_wavelength/10
    warning(['Minimum edge length (' num2str(min_edge_length) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end


%% Probe positions (exposed and shadow side)
observation_edges = [2419 1316 21328];


%% BEUT MOT
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename]);
mesh = BEUT.Main.MOT(mesh, boundary, operator_file, observation_edges,...
    time, mu0(1), 0, 0, V_source, observation_edges(1));


%% Time domain plots
tstop = size(mesh.fields.E_z,2);
figure; plot(time(1:tstop),mesh.fields.E_z(observation_edges,1:tstop))
entries = cell(1,numel(observation_edges));
for i=1:numel(observation_edges)
    entries(i) = {sprintf('E_z at halfedge %i',observation_edges(i))};
end
legend('String',entries);

