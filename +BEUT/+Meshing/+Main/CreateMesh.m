% Create 2D meshes using distmesh2d function and save it ready for simulations
% see: http://persson.berkeley.edu/distmesh/
% Note this is for simple geometries only, don't expect to mesh anything too complex!
% Use a dedicated mesher and import the file for more complex geometries.


%% cylinder resonator (aka circle mesh)
%
radius = 1;

xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
BoundingBox = [xmin,ymin;xmax,ymax];

% Define discretisation by number of wavelengths
num_wavelengths_per_radius = 2;
largest_wavelength = 1/radius/num_wavelengths_per_radius;
initialEdgeLength = largest_wavelength/11;    % <wavelength/10
N_V = round((2*pi)/acos(1-0.5*(initialEdgeLength/radius)^2));

% Or uncomment these lines to define discretisation by number of boundary edges
% N_V=21;
% initialEdgeLength = radius*sqrt(2-2*cos(2*pi/N_V));

distance_function = @(p) (1/radius)*sqrt(sum(p.^2,2))-1; % same as dcircle
theta = 2 * pi / N_V;
alpha = (0 : theta : 2*pi - theta)';
vertices = radius * [ cos(alpha) sin(alpha) ];
fixed_node_positions = vertices;
filename = 'cyl_res';
%}


%% rectangular resonator (aka rectangular mesh)
%{
N_V=40;
xmin = -2; xmax = 2;
ymin = -1; ymax = 1;

initialEdgeLength = 2*(xmax-xmin+ymax-ymin)/N_V;
distance_function=@(p) BEUT.Meshing.distmesh.drectangle(p,xmin,xmax,ymin,ymax);
fixed_node_positions = [xmin,ymin; xmax,ymin; xmin,ymax; xmax,ymax];
BoundingBox = [xmin,ymin;xmax,ymax];
filename = 'rec_res';
%}


%% cylinder PEC cavity (aka rectangular mesh with subtracted circle)
%{
N_V=60;
xmin = -2; xmax = 2;
ymin = -1; ymax = 1;
inner_radius = 0.5;
inner_circle_centre = [-1,0];

initialEdgeLength = 2*(xmax-xmin+ymax-ymin)/N_V;
distance_function=@(p) BEUT.Meshing.distmesh.ddiff(...
    BEUT.Meshing.distmesh.drectangle(p,xmin,xmax,ymin,ymax),...
    BEUT.Meshing.distmesh.dcircle(p,inner_circle_centre(1),inner_circle_centre(2),inner_radius) );
fixed_node_positions = [xmin,ymin; xmax,ymin; xmin,ymax; xmax,ymax];
BoundingBox = [xmin,ymin;xmax,ymax];
filename = 'PEC_cav';
%}


%% cylinder dielectric cavity mesh (aka circle mesh with meshed circle inside)
%{
inner_radius = 1;
outer_radius = 4;
inner_circle_centre = [-1,1];
outer_circle_centre = [0,0];

xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
BoundingBox = [xmin,ymin;xmax,ymax];

N_V_outer=45;
initialEdgeLength = outer_radius*sqrt(2-2*cos(2*pi/N_V_outer));

distance_function=@(p) BEUT.Meshing.distmesh.dcircle(p,outer_circle_centre(1),outer_circle_centre(2),outer_radius);

N_V_inner = round((2*pi)/acos(1-0.5*(initialEdgeLength/inner_radius)^2));
theta = 2 * pi / N_V_inner;
alpha = (0 : theta : 2*pi - theta)';
vertices = inner_radius * [ inner_circle_centre(1)+cos(alpha) inner_circle_centre(2)+sin(alpha) ];
fixed_node_positions = [vertices];
filename = 'cyl_dielectric_cav';
%}


%% rectangular dielectric cavity mesh (aka rectangular mesh with meshed circle inside)
%{
N_V = 20;
inner_radius = 1;
inner_circle_centre = [1,0];

xmin = -1; xmax = 6;
ymin = -3; ymax = 3;
BoundingBox = [xmin,ymin;xmax,ymax];

initialEdgeLength = 0.8*inner_radius*sqrt(2-2*cos(2*pi/N_V));
distance_function=@(p) BEUT.Meshing.distmesh.drectangle(p,xmin,xmax,ymin,ymax);

theta = 2 * pi / N_V; alpha = (0 : theta : 2*pi - theta)';
vertices = inner_radius * [ inner_circle_centre(1)+cos(alpha) inner_circle_centre(2)+sin(alpha) ];

% vertices = [vertices; vertices + repmat([10 0],size(vertices,1),1)];      % 2 circles

fixed_node_positions = [vertices];
filename = 'rec_dielectric_cav';
%}


%% elliptical resonator (aka ellipse mesh)
%{
N_V=20;
radius_x = 1;
radius_y = 0.6;
slopeAngle = pi/4;
co=cos(slopeAngle);
si=sin(slopeAngle);
theta = 2 * pi / N_V;
alpha = (0 : theta : 2*pi - theta)';
vertices = [ radius_x*cos(alpha)*co - radius_y*sin(alpha)*si,...
    radius_y*cos(alpha)*si + radius_y*sin(alpha)*co ];

xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
BoundingBox = [xmin,ymin;xmax,ymax];

initialEdgeLength = radius_y*sqrt(2-2*cos(2*pi/N_V));
distance_function=@(p) BEUT.Meshing.distmesh.dpoly(p,vertices);
fixed_node_positions=[];

filename = 'ellipse_res';
%}


%% airfoil mesh (NACA0025 airfoil)
%{
hlead=0.005; htrail=0.05;
t = 0.25; c = 1;
a=5*t*c*[0.2969/sqrt(c),-0.1260/c,-0.3516/c^2,0.2843/c^3,-0.1015/c^4];

distance_function=@(p) ...
    (abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1);

fixx=1-htrail*cumsum(1.3.^(0:4)');
fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
fixed_node_positions=[fixx,fixy; fixx,-fixy];
xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
BoundingBox = [xmin,ymin;xmax,ymax];
initialEdgeLength=min(hlead,htrail);

filename = 'airfoil';
%}


%% Create mesh
scaled_edge_function = @BEUT.Meshing.distmesh.huniform;       % advancing wavefront mesh
[v,f] = BEUT.Meshing.distmesh.distmesh2d(distance_function,...
        scaled_edge_function,...
        initialEdgeLength,...
        BoundingBox,...
        fixed_node_positions);

numberFaces = size(f,1);
fnum = ones(1,numberFaces);     % all traingle faces belong to the same material


%% Extract inner circle and set its triangles to a different face number
% fnum = BEUT.Meshing.Main.setInnerCircleAsDifferentMaterial(v,f,inner_radius,inner_circle_centre);


%% Save for Matlab for C++ code compatibility
mesh = BEUT.UTLM.UTLMClass(v,f,fnum);
mesh.checkMesh;                    % look for a ratio below 10, ideally below 5
mesh.plot_mesh;

file_name = [filename num2str(numel(mesh.mesh_boundary)) '.mat'];
c_file = [BEUT.CFolder filesep 'input' filesep file_name];
mat_file = [fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep file_name];

use_dual = true;
BEUT.Meshing.save( mesh, use_dual, mat_file, c_file );
