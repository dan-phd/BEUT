function [mesh, M_TM, J_TM] = MOT( mesh, boundary, operator_file, observation_points, time,...
    material_param, Ez_i, Hxy_i, V_source, source_he  )
%BEUT MOT


% Parameters
N_V = boundary.N_V;
N_T = length(time);

% Check if there should be a point source excitation
if nargin < 9 || all(V_source==0)
    UTLM_source = false;
else
    UTLM_source = true;
end

% Check if there should be a plane wave excitation
if Ez_i == 0
    Ez_i = zeros(N_V, N_T);
end
if Hxy_i == 0
    Hxy_i = zeros(N_V, N_T);
end

% Check if the user wants a real-time plot
realTimePlot = true;
if observation_points == 0
    realTimePlot = false;
end

% Get operators from file output by C++ program
assert(N_T<=operator_file.N_T,...
    'Chosen N_T must be less than or equal to the number of timesteps computed in the operators');
N = operator_file.N(:,:,1:N_T)/material_param;
S = operator_file.S(:,:,1:N_T)*material_param;
D = operator_file.D(:,:,1:N_T);
Dp = operator_file.Dp(:,:,1:N_T);


% Hybridisation
Y_TLM_vec(1:N_V) = vertcat(mesh.halfedges(mesh.mesh_boundary(1:N_V)).Y_link) +...
                   vertcat(mesh.halfedges(mesh.mesh_boundary(1:N_V)).Y_stub);
Z_TLM = diag(1./Y_TLM_vec);
Y_TLM = diag(Y_TLM_vec);
edge_lengths = vertcat(boundary.halfedges.l);
l = diag(edge_lengths);
inv_l = inv(l);

Z = zeros(2*size(S,1),2*size(S,2),N_T);
Z(:,:,1) = [-D(:,:,1)           Z_TLM/2+S(:,:,1)*inv_l ;...
            Y_TLM/2+N(:,:,1)*l  Dp(:,:,1) ];

Z0_inv = inv(Z(:,:,1));


% Initialisation
mesh.reset;
M_TM = zeros(boundary.N_V,N_T);
J_TM = zeros(boundary.N_V,N_T);
unknown = zeros(2*boundary.N_V,N_T);
V_open = zeros(boundary.N_V,N_T);
I_closed = zeros(boundary.N_V,N_T);

if realTimePlot
    figure;
    title('This figure can be closed safely if not required')
    for i=1:numel(observation_points)
        h(i) = animatedline;
    end;
    h(numel(observation_points)+1) = animatedline('LineStyle',':');
end

bar = waitbar(0,'Timestepping...');


% Timestepping procedure
tic
for j=1:N_T
    waitbar(j/N_T,bar);
    
    % TLM scatter and store V_open and I_closed
    mesh.scatter(j);
    V_open(:,j) = mesh.V_open;
    I_closed(:,j) = mesh.I_open;
    
    
    % MOT
    rhs_vec = [ Ez_i(:,j) - V_open(:,j)/2; l*Hxy_i(:,j) + I_closed(:,j)/2 ];
    
    Z(:,:,j) = [-D(:,:,j)   S(:,:,j)*inv_l ;...
                N(:,:,j)*l  Dp(:,:,j) ];

    for k=1:j-2
        rhs_vec = rhs_vec - Z(:,:,k+1)*unknown(:,j-k+1);
    end
    
    unknown(:,j) = Z0_inv * rhs_vec;
    M_TM(:,j) = unknown(1:N_V,j);           % E_z = M = V
    J_TM(:,j) = -unknown(N_V+1:2*N_V,j);    % n x H_xy = -J = I
    
    
    % Plug back into TLM
    mesh.connect(j,M_TM(:,j));
    

    % Point source excitation from TLM (if specified)
    if UTLM_source
        mesh.excite_E(source_he, V_source(j));
    end
    
    
    % Display plot in real time (if required)
    if realTimePlot
        for p=1:numel(observation_points)
            addpoints(h(p),time(j),mesh.fields.E_z(observation_points(p),j))
        end
        drawnow update
    end
    
end
toc
close(bar);


end

