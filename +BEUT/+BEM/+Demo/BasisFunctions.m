% Demonstrate basis functions
clear all

%% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);


%% Create basis functions
% Create hat functions that have a unit height:
hat_function = BEUT.BEM.BasisFunction.createHat(boundary.halfedges,false);
% Create square functions that have a height equal to the edge length:
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,true);
% Create functions that are the divergence of the hat functions
% (these should be equal to the previously created square functions):
div_hat_function = BEUT.BEM.BasisFunction.divergence(hat_function,boundary.halfedges);

%% Plot
BEUT.BEM.BasisFunction.plot_basis(boundary.halfedges,[hat_function;square_function;div_hat_function])


%% Create basis functions on dual mesh
hat_function = BEUT.BEM.BasisFunction.createDualHat(boundary.dual,true);
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.dual,true);

%% Plot
BEUT.BEM.BasisFunction.plot_basis(boundary.dual,[hat_function;square_function])