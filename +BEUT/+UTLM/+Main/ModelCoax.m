% Model PEC co-ax
clear all;
global mu0 eps0;


%% Read mesh
load([fileparts(which...
    ('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res20.mat']);
radius=max(range(mesh.vertices))/2;
mesh.plot_mesh;
c0 = 1/sqrt(eps0*mu0);
NT = 2000;
time = 0:dt:(NT-1)*dt;


%% Setup TLM
% Set the material parameters
eps_r(1) = 1; mu_r(1) = 1;      % outside sylinder is free space
eps_r(2) = inf; mu_r(2) = 1;    % inside cylinder is PEC
mesh.setMaterial(eps_r,mu_r);

% Calculate admittances using dt
mesh.calcAdmittances;

% set boundary condition:
% "absorbing" = 0
% "open circuit" = 1
% "short circuit" = -1
mesh.setBoundary(-1);


%% Set up excitation
min_edge_length = min(vertcat(mesh.halfedges.edgeLength));

width = 3*min_edge_length/c0*10;
delay = 1;
inc_wave = BEUT.Excitation.GaussianWave(width, delay, c0);
min_wavelength = c0/inc_wave.freq_response(time,true);

V_source = inc_wave.eval(time);
figure; plot(time,V_source)
title('Incident wave in the time domain'); xlabel('time');
if min_edge_length>min_wavelength/10
    warning(['Minimum edge length (' num2str(min_edge_length) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end


%% Run TLM
sourceEdges = 20;
mesh = BEUT.UTLM.Main.run(mesh, NT, V_source, sourceEdges);


%% Animate
mesh.animate('E',50);

