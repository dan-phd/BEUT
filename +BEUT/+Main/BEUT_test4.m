% BEUT - test 3 (roughly 2.5 hours)
% Region 1 = background
% Region 2 and 3 = spatially distinct meshes
% Compare BEUT with unmeshed background against UTLM with meshed background,
% with different sizes of domain
global mu0 eps0;


%% Load the geometry (along with dt, mu0 and eps0)
filename_BEUT = 'cyl_res24_x2.mat';
mesh_file = matfile([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename_BEUT]);
dt_BEUT = mesh_file.dt;
mu0 = mesh_file.mu0;
eps0 = mesh_file.eps0;
mesh_BEUT = mesh_file.mesh;
boundary=BEUT.Meshing.MeshBoundary(mesh_BEUT);

filename_UTLM{1} = 'test5_2000.mat';
filename_UTLM{2} = 'test5_4000.mat';
filename_UTLM{3} = 'test5_8000.mat';
filename_UTLM{4} = 'test5_12000.mat';
for i=1:numel(filename_UTLM)
    mesh_file = matfile([fileparts(which('BEUT.Meshing.load'))...
        filesep 'meshes' filesep filename_UTLM{i}]);
    dt_UTLM(i) = mesh_file.dt;
    mesh_UTLM{i} = mesh_file.mesh;
end

% Define material layers
mu_r(1)=1; eps_r(1)=1;
mu_r(2)=1; eps_r(2)=2;      % left circle
mu_r(3)=1; eps_r(3)=3;      % right circle

for i=1:boundary.num_shapes+1
    mu(i) = mu_r(i)*mu0; eps(i) = eps_r(i)*eps0;
    eta(i) = sqrt(mu(i)/eps(i));
    c(i) = 1/sqrt(mu(i)*eps(i));
end
mesh_BEUT.setMaterial(eps_r(2:3),mu_r(2:3));
mesh_BEUT.calcAdmittances;
for i=1:numel(filename_UTLM)
    mesh_UTLM{i}.setMaterial(eps_r,mu_r);
    mesh_UTLM{i}.calcAdmittances;
    mesh_UTLM{i}.setBoundary(0);
end

% Temporal parameters
N_T_BEUT = 3500;
time_BEUT = 0:dt_BEUT:(N_T_BEUT-1)*dt_BEUT;
N_T_UTLM(1) = 8600;
N_T_UTLM(2) = 7450;
N_T_UTLM(3) = 8600;
N_T_UTLM(4) = 5750;


%% Set up point source excitation
min_edge_length = min(vertcat(boundary.halfedges.l));
width = 8/c(1);
delay = 1.2;
inc_wave = BEUT.Excitation.GaussianWave(width,delay,c(1));
inc_wave.direction = [1 0];


%% Probe positions
observation_edges_BEUT = [168 58 673 613];
source_edges_BEUT = observation_edges_BEUT(1);
for i=1:numel(filename_UTLM)-1
    observation_edges_UTLM{i} = [168 58 673 613];
    source_edges_UTLM{i} = observation_edges_UTLM{i}(1);
end
observation_edges_UTLM{4} = [162 58 391 598];
source_edges_UTLM{4} = observation_edges_UTLM{4}(1);


%% BEUT MOT
V_source = inc_wave.eval(time_BEUT);
operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename_BEUT]);
mesh_BEUT = BEUT.Main.MOT(mesh_BEUT, boundary, operator_file, observation_edges_BEUT,...
    time_BEUT, mu(1), 0, 0, V_source, source_edges_BEUT);


%% Use UTLM to run the same tests with different sized domains
for i=1:numel(filename_UTLM)
    time_UTLM{i} = 0:dt_UTLM(i):(N_T_UTLM(i)-1)*dt_UTLM(i);
    V_source = inc_wave.eval(time_UTLM{i});
    mesh_UTLM{i} = BEUT.UTLM.Main.run(mesh_UTLM{i}, N_T_UTLM(i), V_source, source_edges_UTLM{i});
end


%% Interpolate graphs (since they all have different dt) and plot comparison
observation_idx = 4;
time = linspace(0,3e-7,10000);
E_z_BEUT = interp1(time_BEUT,mesh_BEUT.fields.E_z(observation_edges_BEUT(observation_idx),:),time);
for i=1:numel(time_UTLM)
    E_z_UTLM{i} = interp1(time_UTLM{i},mesh_UTLM{i}.fields.E_z(observation_edges_UTLM{i}(observation_idx),:),time);
end

figure; hold on;
plot(time,E_z_BEUT,'LineWidth',4);
legend('BEUT');
for i=1:numel(time_UTLM)
    plot(time,E_z_UTLM{i})
    names = legend;
    h = legend([names.String {sprintf('UTLM with %i triangles',mesh_UTLM{i}.nF)}]);
end
set(gca,'yticklabel',{})
xlabel('time'); ylabel('E_z');
title('UTLM convergence to BEUT solution');


%% Find error between each plot
BEUT.relError(E_z_BEUT,E_z_UTLM{1},E_z_UTLM{2},E_z_UTLM{3},E_z_UTLM{4})

