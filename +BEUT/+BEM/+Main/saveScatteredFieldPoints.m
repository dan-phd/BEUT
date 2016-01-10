function in_shape = saveScatteredFieldPoints( mesh, x_coords,y_coords, M,J, dual, c_file, material_parameter, c )
%saves a file to be used by C++ program to calculate the scattered field
% where rho is a list of co-ordinates to resolve the field at

assert(isa(mesh,'BEUT.UTLM.UTLMClass'),'mesh input must be a UTLMClass object');
boundary_mesh=BEUT.Meshing.MeshBoundary(mesh);

% decide whether to use dual basis functions
if dual
    boundary = boundary_mesh.dual; dual = 1;
else
    boundary = boundary_mesh.halfedges; dual = 0;
end

% Remove points that are inside mesh
material_vertices = vertcat(boundary.a);
in_shape = inpolygon(x_coords,y_coords,material_vertices(:,1),material_vertices(:,2));
X_out = x_coords(~in_shape); Y_out = y_coords(~in_shape);
rho = [X_out Y_out];
figure; plot(material_vertices(:,1),material_vertices(:,2)); axis equal; hold on
plot(rho(:,1),rho(:,2),'k.') % points strictly outside
hold off
title('Observation points to find scattered field');

% Convert rho from matrix to struct
assert(size(rho,2)==2,'rho must be a list of coordinates (with 2 columns for x and y)');
rho_tbl = table;
rho_tbl.x = rho(:,1);
rho_tbl.y = rho(:,2);
rho = table2struct(rho_tbl);
disp(['Number of grid points: ' num2str(size(rho,1))])

% Speed of propagation
eps0 = 8.854187817e-12;     % permittivity of free space
mu0 = 4*pi*10^-7;           % permeability of free space
if nargin<8
    material_parameter = mu0;
end
if nargin<9
    c = 1/sqrt(mu0*eps0);
end

% timestep
dt=mesh.dt;

% Find number of shapes in mesh
num_shapes = boundary_mesh.num_shapes;


%% output for C++ if required
if(nargin>3)
    if(~isempty(c_file))
        save(c_file, 'boundary', 'dt', 'c', 'num_shapes', 'dual', 'rho', 'M', 'J', 'material_parameter', '-v7');
        disp(['C++ file output to: ' c_file])
    end
end


end

