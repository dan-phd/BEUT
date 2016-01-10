function plot_interfaces(obj)
% Plot mesh interfaces (between different materials)

figure;
axis equal;
hold on

for m=1:obj.num_materials
    
    halfedges = obj.material_boundaries{m};
    
    for he_ind = 1:numel(halfedges)
        % Plot edge line
        vec = obj.vertices(obj.halfedges(halfedges(he_ind)).vertices,:);
        plot(vec(:,1),vec(:,2),'LineWidth',4,'Color',obj.color(1,:))
        
    end
    
end

end
