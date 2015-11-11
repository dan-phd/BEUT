%% Animate scattered fields (test demo - TODO: need to fix)
domain_X = -3:0.3:3;
domain_Y = -3:0.3:3;
[X,Y] = meshgrid(domain_X,domain_Y);
fieldN_T = 80:83;

% Vectorise all points in grid
rho = [X(:) Y(:)];

material_vertices = vertcat(boundary.halfedges.a);

% Compute Z matrices
zm2 = BEUT.BEM.ZMatrices(N_T,dt,boundary.halfedges,c(1));
zm2.basis_function_Z = BEUT.BEM.BasisFunction.createSquare(boundary.halfedges,true);
zm2.basis_function_S = BEUT.BEM.BasisFunction.createHat(boundary.halfedges,false);
timeBasis = BEUT.BEM.LagrangeInterpolator(dt,1);
zm2.timeBasis_D = timeBasis;
zm2.timeBasis_Nh = int(timeBasis);
zm2.timeBasis_Ns = diff(timeBasis);
zm2.inner_points = 4;
tic
[S,D,Dp,Nh,Ns] = zm2.computeField(rho, fieldN_T);
toc
N = Nh+Ns/c(1)^2;

Z = [D         -S*mu(1);...
     -N/mu(1)  -Dp     ];
x = [ M_TM; J_TM];


%% Plot
%{
rhs = Z(:,:,fieldN_T) * x(:,fieldN_T);
E_vec = rhs(1:numel(X));
E_s = reshape(E_vec,size(X));
surf(X,Y,E_s); view(2); colorbar;
hold on;
above_animation = max(max(max(E_s)))*ones(length(material_vertices),1);
scatter3(material_vertices(:,1),material_vertices(:,2),above_animation,10,...
    'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
%}

%% animate
%
E_s = zeros([size(X) numel(fieldN_T)]);
for i=fieldN_T
    rhs = Z(:,:,i)*x(:,i);
    E_vec = rhs(1:numel(X),:);
    E_s(:,:,i) = reshape(E_vec,size(X));
end

% [E_s] = BEUT.BEM.Main.findScatteredField(boundary.halfedges, gpw, c(1), time, M_TM, X,Y );

BEUT.animate_fields(2,'domain',X,Y,...
    'animation',E_s,...
    'overlay',material_vertices);
%}