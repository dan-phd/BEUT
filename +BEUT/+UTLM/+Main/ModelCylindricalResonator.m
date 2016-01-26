% Demonstrate scattering within a cylindrical resonator using UTLM
% Compare results with analytical
% NOTE: to hit more analytical resonant frequencies, increase the requency
% spectrum of the excitation. However, increasing it too far will lead to 
% artificial dispersion, only using a finer mesh will solve this.
clear all;
global mu0 eps0;


%% Read mesh
load([fileparts(which...
    ('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
radius=max(range(mesh.vertices))/2;
c0 = 1/sqrt(eps0*mu0);
NT = 10000;
time = 0:dt:(NT-1)*dt;


%% Setup TLM
% Set the material parameters
eps_r = 1; mu_r = 1;
mesh.setMaterial(eps_r,mu_r);

% Calculate admittances using dt
mesh.calcAdmittances;

% set boundary condition:
% "absorbing" = 0
% "open circuit" = 1
% "short circuit" = -1
mesh.setBoundary(-1);


%% Set up excitation
sourceEdges = 1;

min_edge_length = min(vertcat(mesh.halfedges.edgeLength));
width = 2.5*min_edge_length/c0*10;
delay = 1;

inc_wave = BEUT.Excitation.GaussianWave(width,delay,c0);
inc_wave.A=0.5;
V_source = inc_wave.eval(time);
figure; plot(time,V_source)
title('Incident wave in the time domain'); xlabel('time');
min_wavelength = c0/inc_wave.freq_response(time,true);
if min_edge_length>min_wavelength/10
    warning(['Minimum edge length (' num2str(min_edge_length) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end


%% Run TLM
mesh = BEUT.UTLM.Main.run(mesh, NT, V_source, sourceEdges);


%% Check results vs analytical solution
BEUT.UTLM.Main.plotAgainstAnalyticalCylinder(mesh.fields.E_z,time,c0,radius);
