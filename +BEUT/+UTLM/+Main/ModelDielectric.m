% Model scattering by a cylindrical dielectric using UTLM
global mu0 eps0;


%% Paramaters
NT = 4000;
mu0 = 4*pi*10^-7;           % permeability of free space
eps0 = 8.854187817e-12;     % permittivity of free space
c0 = 1/sqrt(eps0*mu0);


%% Read mesh
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_dielectric_cav20.mat']);
radius=max(range(mesh.vertices))/2;


%% Setup TLM
% Set the material parameters
eps_r(1) = 1; mu_r(1) = 1;
eps_r(2) = 3; mu_r(2) = 1;
mesh.setMaterial(eps_r,mu_r);
mesh.plot_materials('eps_r');

% Find dt
dt = mesh.shortestLinkLength*sqrt(2*eps0*mu0) /4;       % dt < minLinklength*sqrt(2*eps0*mu0)
mesh.dt=dt;
time = 0:dt:(NT-1)*dt;

% Calculate admittances using dt
mesh.calcAdmittances;

% set boundary condition:
% "absorbing" = 0
% "open circuit" = 1
% "short circuit" = -1
mesh.setBoundary(0);


%% Set up excitation
% Find face closest to centre [0 0] for source node
CC=vertcat(mesh.halfedges.circumcenter);
for i=1:length(CC), normalised_CC(i)=norm(CC(i,:)); end
centerEdges = find(normalised_CC==min(normalised_CC));
sourceEdges = centerEdges(2);

min_edge_length = min(vertcat(mesh.halfedges.edgeLength));
width = min_edge_length/c0*10;
delay = 1;

inc_wave = BEUT.Excitation.GaussianWave(width, delay, c0);
inc_wave.A=0.5;
V_source = inc_wave.eval(time);
figure; plot(time,V_source)
min_wavelength = c0/inc_wave.freq_response(time,true);
if min_edge_length>min_wavelength/10
    warning(['Minimum edge length (' num2str(min_edge_length) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end


%% Run TLM
mesh = BEUT.UTLM.Main.run(mesh, NT, V_source, sourceEdges);


%% Plot at various points around the boundary
time=0:dt:(NT-1)*dt;
observation_point = mesh.mesh_boundary(1:5:end);

% mesh.plot_halfedge(observation_point)
figure, plot(time,mesh.fields.E_z(observation_point,:));
for i=1:numel(observation_point)
    [~,~,~,names] = legend;
    h = legend([names {sprintf('Halfedge %i',observation_point(i))}]);
end


%% Animate
mesh.animate('E');

