% Demonstrate computation of Gram matrix

% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);

Gram = BEUT.BEM.GramMatrix();

% Compute the Gram matrix on the dual mesh
applied_geometry = boundary.dual;
Gram.geometry = applied_geometry;
Gram.test_function = BEUT.BEM.BasisFunction.createDualHat(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createDualSquare(applied_geometry,true);
G1 = compute(Gram);

% Compute the Gram matrix on the normal mesh
applied_geometry = boundary.halfedges;
Gram.geometry = applied_geometry;
Gram.test_function = BEUT.BEM.BasisFunction.createHat(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createSquare(applied_geometry,true);
G2 = compute(Gram);
