% Load a mesh from file, extract materials from it, plot it, and then save it as a halfedge mesh we can use
clear all


%% Load mesh
mesh_name = 'test5_2000'; extension = 'mphtxt';
[v,f,fnum]=BEUT.Meshing.load(...
    [fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep...
    'unconverted' filesep mesh_name '.' extension]);


%% Extract materials 2 and 3, make them materials 1 and 2
% [f,fnum] = BEUT.Meshing.Main.extractMaterials(f,fnum,[2 3],[1 2]);


%% Extract inner circle and set its triangles to a different face number
% fnum = BEUT.Meshing.Main.setInnerCircleAsDifferentMaterial(v,f,0.35);


%% Convert to halfedge mesh
mesh = BEUT.UTLM.UTLMClass(v,f,fnum);
mesh.check_mesh;                    % look for a ratio below 10, ideally below 5
mesh.plot_mesh;

filename = [mesh_name '.mat'];
c_file = [BEUT.CFolder filesep 'input' filesep filename];
mat_file = [fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep filename];

BEUT.Meshing.save( mesh, mat_file, '' );
