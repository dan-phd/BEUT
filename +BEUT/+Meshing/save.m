function save( mesh, dual, mat_file, c_file, dt, mu0, eps0 )
%save(mesh) converts the UTLM mesh for input to other parts of the system
%   mesh is a UTLMClass object, dual (boolean) 
%   c_file is the path and filename to output the file that C++ will use
%   mat_file is the path and filename to output the file that Matlab will use
%   dt, mu0 and eps0 are optional and will be calculated automatically if not given
%
% example:
% load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat'])
% mat_file = [fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat'];
% c_file = '<path to 2DTDBEM project>\input\cyl_res21.mat';
% BEUT.Meshing.save( mesh, true, mat_file, c_file );

assert(isa(mesh,'BEUT.UTLM.UTLMClass'),'mesh input must be a UTLMClass object');
boundary_mesh=BEUT.Meshing.MeshBoundary(mesh);

% decide whether to use dual basis functions
if dual
    boundary = boundary_mesh.dual; dual = 1;
else
    boundary = boundary_mesh.halfedges; dual = 0;
end

% Speed of propagation
if nargin<7
    eps0 = 8.854187817e-12;     % permittivity of free space
    if nargin<6
        mu0 = 4*pi*10^-7;           % permeability of free space
    end
end
c = 1/sqrt(mu0*eps0);

% Check dt (if given) < minLinklength*sqrt(2*eps0*mu0)
dt_ = mesh.shortestLinkLength*sqrt(2*eps0*mu0) /1.6;
if nargin<5
    dt = dt_;
else
    assert(dt<=dt_,'Chosen timestep is unstable, dt < minLinklength*sqrt(2*eps0*mu0)');
end
mesh.dt=dt;

% Find number of shapes in mesh
num_shapes = boundary_mesh.num_shapes;


%% output for Matlab
if(~isempty(mat_file))
    save(mat_file, 'mesh', 'dt', 'mu0', 'eps0');
    disp(['Matlab file output to: ' mat_file])
end


%% output for C++ if required
if(nargin>2)
    if(~isempty(c_file))
        save(c_file, 'boundary', 'dt', 'c', 'num_shapes', 'dual', '-v7');
        disp(['C++ file output to: ' c_file])
    end
end



end

