% Demonstrate computation of Gram matrix and perform eigenanalysis of resulting matrices

% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res69.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);

Gram = BEUT.BEM.GramMatrix();

%% Compute the Gram matrix on the normal mesh
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
figure; hold on;
plot(real(svd(G_ss)),'LineWidth',2);
plot(real(svd(G_sh)),'LineWidth',2);
plot(real(svd(G_hs)),'LineWidth',2);
plot(real(svd(G_hh)),'LineWidth',2);
L=legend('$\left\langle \sqcap,\sqcap \right\rangle$','$\left\langle \sqcap,\wedge \right\rangle$',...
    '$\left\langle \wedge,\sqcap \right\rangle$','$\left\langle \wedge,\wedge \right\rangle$');
set(L,'Interpreter','latex', 'FontSize', 15);
ylabel('Singular values', 'FontSize', 18)
xlabel('Number of elements off diagonal', 'FontSize', 18)


%% Compute the Gram matrix on the dual mesh
%{
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
figure; hold on;
plot(real(svd(G_ss)),'LineWidth',2);
plot(real(svd(G_sh)),'LineWidth',2);
plot(real(svd(G_hs)),'LineWidth',2);
plot(real(svd(G_hh)),'LineWidth',2);
L=legend('Dual $\left\langle \sqcap,\sqcap \right\rangle$','Dual $\left\langle \sqcap,\wedge \right\rangle$',...
    'Dual $\left\langle \wedge,\sqcap \right\rangle$','Dual $\left\langle \wedge,\wedge \right\rangle$');
set(L,'Interpreter','latex', 'FontSize', 15);
ylabel('Singular values', 'FontSize', 18)
xlabel('Number of elements off diagonal', 'FontSize', 18)
%}
