% Demonstrate basis functions
clear all;

%% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);


%% Create basis functions
% Create square functions that have a unit height:
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,false);
% Create hat functions that have a height equal to 1/edge length:
hat_function = BEUT.BEM.BasisFunction.createHat(boundary.halfedges,true);
% Create functions that are the divergence of the hat functions:
div_hat_function = BEUT.BEM.BasisFunction.divergence(hat_function,boundary.halfedges);

%% Plot
BEUT.BEM.BasisFunction.plot_basis(boundary.halfedges,[hat_function;square_function;div_hat_function])
title('Square (with unit height), hat and div hat basis functions')

%% Now create basis functions on the dual mesh
hat_function = BEUT.BEM.BasisFunction.createDualHat(boundary.dual,true);
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.dual,true);

%% Plot
BEUT.BEM.BasisFunction.plot_basis(boundary.dual,[hat_function;square_function])
title('Dual square and dual hat basis functions')