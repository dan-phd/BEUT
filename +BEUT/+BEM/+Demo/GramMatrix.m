% Demonstrate computation of Gram matrix and perform eigenanalysis of resulting matrices

% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);

Gram = BEUT.BEM.GramMatrix();

% Compute the Gram matrix on the normal mesh
applied_geometry = boundary.halfedges;
Gram.geometry = applied_geometry;
Gram.test_function = BEUT.BEM.BasisFunction.createSquare(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createSquare(applied_geometry,true);
G_ss = compute(Gram);
Gram.test_function = BEUT.BEM.BasisFunction.createSquare(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createHat(applied_geometry,true);
G_sh = compute(Gram);
Gram.test_function = BEUT.BEM.BasisFunction.createHat(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createSquare(applied_geometry,true);
G_hs = compute(Gram);
Gram.test_function = BEUT.BEM.BasisFunction.createHat(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createHat(applied_geometry,true);
G_hh = compute(Gram);

% Plot eigenanalysis
figure;
subplot(2,2,1), plot(real(svd(G_ss))); title('<squ,squ>');
subplot(2,2,2), plot(real(svd(G_sh))); title('<squ,hat>');
subplot(2,2,3), plot(real(svd(G_hs))); title('<hat,squ>');
subplot(2,2,4), plot(real(svd(G_hh))); title('<hat,hat>');

% Compute the Gram matrix on the dual mesh
applied_geometry = boundary.dual;
Gram.geometry = applied_geometry;
Gram.test_function = BEUT.BEM.BasisFunction.createDualSquare(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createDualSquare(applied_geometry,true);
G_ss = compute(Gram);
Gram.test_function = BEUT.BEM.BasisFunction.createDualSquare(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createDualHat(applied_geometry,true);
G_sh = compute(Gram);
Gram.test_function = BEUT.BEM.BasisFunction.createDualHat(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createDualSquare(applied_geometry,true);
G_hs = compute(Gram);
Gram.test_function = BEUT.BEM.BasisFunction.createDualHat(applied_geometry,true);
Gram.basis_function = BEUT.BEM.BasisFunction.createDualHat(applied_geometry,true);
G_hh = compute(Gram);

% Plot eigenanalysis
figure;
subplot(2,2,1), plot(real(svd(G_ss))); title('dual <squ,squ>');
subplot(2,2,2), plot(real(svd(G_sh))); title('dual <squ,hat>');
subplot(2,2,3), plot(real(svd(G_hs))); title('dual <hat,squ>');
subplot(2,2,4), plot(real(svd(G_hh))); title('dual <hat,hat>');
