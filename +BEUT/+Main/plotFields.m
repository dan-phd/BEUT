function plotFields( filename,mesh,X,Y,in_scattered,plot_timestep,interpolate )
% Once C++ has computed the scattered field operators, run this
% to plot fields inside and outside scatterer

operator_file = matfile([BEUT.CFolder filesep 'results' filesep filename '_scattered.mat']);
E_s = BEUT.BEM.Main.organizeScatteredField(operator_file, X, in_scattered );


%% plot fields outside scatterer
max_E = max(max(max(abs(E_s))));
E_s_ = 2*E_s(:,:,plot_timestep)/max_E;
% E_s_ = 20*log10(abs(E_s)+eps);       % in dB

Parent2 = figure;
colormap('jet');
axes2 = axes('Parent',Parent2,...
    'FontSize',14,...
    'TickDir','out',...
    'CLim',[-1 1]);
hold(axes2,'on');
background = surf2patch(X,Y,zeros(size(E_s_)),E_s_);

background.FaceColor = 'flat'; shading faceted;       % for block colors
if nargin>6
    if interpolate==true
        background.FaceColor = 'interp';  shading interp;   % for interpolated fields
    end
end

bg = patch(background,'Parent',axes2,'EdgeAlpha',0.1);

colorbar('peer',axes2);
xlabel('x (m)','FontSize',18);
ylabel('y (m)','FontSize',18);


%% plot fields inside scatterer
S.Vertices = mesh.vertices;
S.Faces = vertcat(mesh.faces.vertices);
S.LineStyle = '-';
S.FaceColor = 'flat';
S.FaceVertexCData = 2*mesh.V0(:,plot_timestep)/max_E;
scatterer = patch(S,'Parent',axes2,'EdgeAlpha',0.1);
axis(axes2,'tight');
uistack(scatterer, 'top')

%% plot scatterer outline
for m=1:mesh.num_materials
    
    halfedges = mesh.material_boundaries{m};
    
    for he_ind = 1:numel(halfedges)
        % Plot edge line
        vec = mesh.vertices(mesh.halfedges(halfedges(he_ind)).vertices,:);
        plot(vec(:,1),vec(:,2),'LineWidth',3,'Color','black')
        
    end
    
end

axis equal;

end

