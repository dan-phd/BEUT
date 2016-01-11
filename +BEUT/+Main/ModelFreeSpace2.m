% BEUT - simple test 2
% Region 1 = UTLM cylinder mesh of free space
% Region 2 = free space BEM
% Excite with a point source at various source edges and compare against UTLM
clear all;
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
mesh.setBoundary(0);

% Temporal parameters
N_T = 1000;
time = 0:dt:(N_T-1)*dt;


%% Set up excitation
% Gaussian pulse properties
edge_lengths = vertcat(boundary.halfedges.l);
width = min(edge_lengths)/c(1)*30;
delay = 1.2;
inc_wave = BEUT.Excitation.GaussianWave(width,delay,c(1));
inc_wave.direction = [1 0];
V_source = inc_wave.eval(time);
figure; plot(time,V_source)
title('Incident wave in the time domain'); xlabel('time');

% check stability
min_wavelength = c(1)/inc_wave.freq_response(time,false);
if min(edge_lengths)>min_wavelength/10
    warning(['Minimum edge length (' num2str(min(edge_lengths)) ...
        ') should be less than a tenth of the minimum wavelength ('...
        num2str(min_wavelength) ')'])
end


%% Probe positions
source_edges = [263 64 58];
observation_edges = [277 106 206];
mesh.plot_halfedge([source_edges observation_edges])


%% MOT, exciting from 3 different positions
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename]);

for i=1:numel(source_edges)
    mesh = BEUT.Main.MOT(mesh, boundary, operator_file, 0,...
        time, mu(1), 0, 0, V_source, source_edges(i));
    Ez_BEUT{i} = mesh.fields.E_z;
end


%% Use UTLM to run the same tests
for i=1:numel(source_edges)
    mesh = BEUT.UTLM.Main.run(mesh, N_T, V_source, source_edges(i));
    Ez_TLM{i} = mesh.fields.E_z;
end


%% Plot results
tstop=size(Ez_TLM{1},2);
figure1 = figure; axes1 = axes('Parent',figure1,...
    'FontSize',14,'YTickLabel',{});
hold(axes1,'on');
title('Comparison between BEUT and UTLM at different locations in free space')
xlabel('Time (s)','FontSize',20);
ylabel('Normalized electric field','FontSize',20);
interval=8; marker = {'k' 'ok' 'sk'};
for i=1:numel(observation_edges)
    % Decimate plots to use markers
    h(i) = plot(time(1:interval:tstop),Ez_BEUT{i}(observation_edges(i),1:interval:tstop),marker{i});
    set(h(i),'DisplayName',['BEUT - position ' num2str(observation_edges(i))]);
end
for i=1:numel(observation_edges)
    h(i+numel(observation_edges)) = plot(time(1:tstop),Ez_TLM{i}(observation_edges(i),1:tstop),'LineWidth',2);
    set(h(i+numel(observation_edges)),'DisplayName',['UTLM - position ' num2str(observation_edges(i))]);
end
legend(axes1,h);

