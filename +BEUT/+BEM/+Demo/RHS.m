% Demonstrate computation of RHS

c = 3e8;
N_T = 1000;
dt = 1 / c;

%% Load the geometry
load([fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep 'cyl_res21.mat']);
boundary=BEUT.Meshing.MeshBoundary(mesh);
N_V = boundary.N_V;

%% Create excitation
% Lets make 2 different kinds of waves, we can choose which one to use
% later on
pulseWidth = 5*dt;
startTime = pulseWidth;
polarization = [0 1];
direction = [1 0];
gpw = BEUT.Excitation.GaussianWave(pulseWidth, startTime, c, direction);
desiredFreqWidth = 1e8;
desiredModulatedFreq = 2e8;
sine = BEUT.Excitation.SineWave(desiredFreqWidth, desiredModulatedFreq,...
    c, direction, 1);

%% Create basis functions
square_function = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,false);
hat_function = BEUT.BEM.BasisFunction.createHat(boundary.halfedges,false);
dual_square_function = BEUT.BEM.BasisFunction.createDualSquare(boundary.dual,true);
dual_hat_function = BEUT.BEM.BasisFunction.createDualHat(boundary.dual,false);

%% Make RHS object
rhsCalc = BEUT.BEM.RHS(N_T, dt);
rhsCalc.excitation = @sine.eval;    % this is where we specify which wave to use
rhsCalc.Gaussian_points = 3;        % how many Gaussian quadrature points to use
rhsCalc.polarization = polarization;
rhsCalc.display_plot = true;
tangent = true;                     % this boolen (used in the compute argument)
                                    % enables/disables taking the geometry 
                                    % tangent into account (usually used
                                    % for fields transverse to the plane)

%% Compute RHS using square testing functions
rhsCalc.geometry = boundary.halfedges;
rhsCalc.test_function = square_function;
tic
V1 = rhsCalc.compute(tangent); title('RHS tested with square functions')
toc

%% Compute RHS using square testing functions
rhsCalc.geometry = boundary.halfedges;
rhsCalc.test_function = hat_function;
tic
V2 = rhsCalc.compute(tangent); title('RHS tested with hat functions')
toc

%% Compute RHS using dual hat testing functions
rhsCalc.geometry = boundary.dual;
rhsCalc.test_function = dual_hat_function;
tic
V3 = rhsCalc.compute(tangent); title('RHS tested with dual hat functions')
toc

