function checkMesh(obj)
% Plot Voronoi diagram and check mesh has sufficiently large TLM link lengths

% Plot triangulation
figure, triplot(obj.TR,'k'); axis equal; hold on;

min_length = inf;
for he_ind = 1:obj.nH
    
    % Face circumcenter
    CC = obj.halfedges(he_ind).circumcenter;
    
    % Opposite halfedge index
    opp_he = obj.halfedges(he_ind).flip;
    
    % set the color number
    colnum = obj.faces(obj.halfedges(he_ind).face).fnum;
    colnum = mod(colnum-1,size(obj.color,1))+1;
    
    % If we're not on a boundary edge...
    if  opp_he~=0
        
        % If the material is different on both sides of the edge
        if ~(obj.faces(obj.halfedges(he_ind).face).fnum...
                ==obj.faces(obj.halfedges(opp_he).face).fnum)
            
            plot_voronoi(CC,obj.halfedges(he_ind).midpoint,obj.color(colnum,:));
            
        else
            
            % The link length is simply half the distance from circumcenter to circumcenter
            NCC = obj.halfedges(opp_he).circumcenter;     % neighbour circumcenter
            
            plot_voronoi(CC,NCC,obj.color(colnum,:));
            
        end
        
        % save minimum link length
        if obj.halfedges(he_ind).linkLength<min_length
            min_length=obj.halfedges(he_ind).linkLength;
            min_link_start = CC;
        end
        
    else        % If we are linking to the boundary...
        
        plot_voronoi(CC,obj.halfedges(he_ind).midpoint,obj.color(colnum,:));
        
        % save minimum link length
        if obj.halfedges(he_ind).linkLength<min_length
            min_length=obj.halfedges(he_ind).linkLength;
            min_link_start = CC;
        end
        
    end
    
    plot_voronoi(CC,obj.halfedges(he_ind).midpoint,obj.color(colnum,:));
    
    % save minimum link length
    if obj.halfedges(he_ind).linkLength<min_length
        min_length=obj.halfedges(he_ind).linkLength;
        min_link_start = CC;
    end
    
end

% plot where the minimum link length is on voronoi diagram
text(min_link_start(:,1),min_link_start(:,2),'O',...
    'HorizontalAlignment','center','FontWeight','bold','FontSize',24,'Color',[0 0.8 1])
hold off

fprintf('\nShortest link length = \t%f\n', min_length)
fprintf('Average link length = \t%f\n', mean(vertcat(obj.halfedges.linkLength)))
fprintf('Ratio = \t\t\t\t%f\n\n', mean(vertcat(obj.halfedges.linkLength))/min_length)


    % plot Voronoi diagram
    function plot_voronoi(v1,v2,col)
        
        plot([v1(:,1) v2(:,1)],[v1(:,2) v2(:,2)],'LineWidth',3,...
            'Color',col);
        
    end

end