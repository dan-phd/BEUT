classdef HalfedgeMesh < handle
    %Create a halfedge structured mesh (suitable for UTLM)
    
    properties
        vertices;
        material_boundaries;        % halfedge indices for material boundaries
        mesh_boundary;              % halfedge indices for domain boundary
        mesh_body;                  % halfedge indices for domain body
        
        faces = struct('vertices',{},'fnum',{});
        
        halfedges = struct('face',{},'vertices',{},'flip',{},...
            'circumcenter',{},'midpoint',{},'edgeLength',{},...
            'linkLength',{});
        
        TR;  % triangulation object which stores Points and ConnectivityList
        nF;         % number of faces
        nV_face;    % number of vertices per face
        nH;         % number of halfedges
        num_materials = 1;     % number of different materials used
        shortestLinkLength;
        
        % color array for cycling through material colors
        color = [   0,      0.447,	0.741;...
                    0.85,   0.325,  0.098;...
                    0.929,  0.694,  0.125;...
                    0.494,  0.184,  0.556;...
                    0.466,  0.674,  0.188;...
                    0.301,  0.745,  0.933;...
                    0.635,  0.078,  0.184 ];
    end
    
    methods
        
        % Input vertices and faces of geometry
        function obj = HalfedgeMesh(vertices, faces, fnum)
            
            % convert to 2D if 3D is given
            if (size(vertices,2)==3)
                vertices(:,3) = [];
            end
            
            % Useful parameters
            obj.nF = size(faces,1);       % number of faces
            obj.nV_face = size(faces,2);  % number of vertices per face
            assert(obj.nV_face==3,'Faces are degenerate or not 2D (do not have 3 vertices)')
            
            face_tbl = table;
            face_tbl.vertices = faces;
            
            % if there are multiple materials
            face_tbl.fnum = ones(obj.nF,1);
            if exist('fnum','var')==1
                assert(length(fnum)==obj.nF,'fnum and faces must have the same lengths')
                if size(fnum,1)~=size(faces,1), fnum=fnum'; end
                face_tbl.fnum = fnum;
                obj.num_materials = max(fnum);
                for ii=1:obj.num_materials
                    assert(any(fnum==ii),['fnum must have materials starting from 1 up to number of materials,'...
                        'fnum is missing material ' num2str(ii)])
                end
            end
            
            obj.vertices = vertices;
            obj.faces = table2struct(face_tbl);
            
            % Convert to Matlab triangulation object (for easy triangulation calculations later on)
            obj.TR=triangulation(faces,vertices);
            
            % Create halfedge structure and set lengths
            create_halfedges(obj);
            setLengths(obj);
            
            % Find boundary vertices between different materials
            find_boundary_edges(obj);
            
            
            % Create halfedge structure
            function create_halfedges(obj)
                bar = waitbar(0,'Creating halfedge structure...');
                for f = 1:obj.nF
                    waitbar(f/obj.nF,bar);
                    for v1 = 1:obj.nV_face
                        v2 = mod(v1, obj.nV_face)+1;
                        
                        % start and end of the edge
                        a = obj.faces(f).vertices(v1);
                        b = obj.faces(f).vertices(v2);
                        
                        % Store halfedge data
                        he_ind = (f-1)*obj.nV_face+v1;        % halfedge index
                        obj.halfedges(he_ind).vertices = [a b];
                        obj.halfedges(he_ind).face = f;
                        
                        % Is there a flip edge to the current halfedge? If not, flip=0 (boundary)
                        [~,flip]=ismember([b a],vertcat(obj.halfedges.vertices),'rows');
                        obj.halfedges(he_ind).flip = flip;
                        if flip~=0
                            obj.halfedges(flip).flip = he_ind;          % set the flip's flip to current halfedge
                        end
                        
                    end
                end
                close(bar);
                
                obj.nH = numel(obj.halfedges);
                
            end % create_halfedges
            
            % Create cells containing boundary halfedges between different materials and different shapes
            function find_boundary_edges(obj)
                
                boundary_index=zeros(obj.num_materials,1);
                for he_ind = 1:obj.nH
                    
                    % Current face index
                    cur_face = obj.halfedges(he_ind).face;
                    
                    % Opposite halfedge index
                    flip_he = obj.halfedges(he_ind).flip;
                    
                    % If we're on  a boundary edge...
                    if  flip_he==0
                        
                        % each cell contains the boundary halfedges of each material
                        boundary_index(obj.faces(cur_face).fnum) = boundary_index(obj.faces(cur_face).fnum)+1;
                        obj.material_boundaries{obj.faces(cur_face).fnum}(boundary_index(obj.faces(cur_face).fnum)) = he_ind;
                        
                    else
                        
                        % Opposite face index
                        flip_face = obj.halfedges(flip_he).face;
                        
                        % If the material is different on both sides of the edge
                        if ~(obj.faces(cur_face).fnum==obj.faces(flip_face).fnum)
                            
                            % each cell contains the boundary halfedges of each material
                            boundary_index(obj.faces(cur_face).fnum) = boundary_index(obj.faces(cur_face).fnum)+1;
                            obj.material_boundaries{obj.faces(cur_face).fnum}(boundary_index(obj.faces(cur_face).fnum)) = he_ind;
                            
                        end
                        
                    end
                end
                
                
                % Find boundary of entire mesh domain
                boundary_vertices = freeBoundary(obj.TR);       % obtains vertex connections
                
                % order the vertices so they are connected (this will allow basis functions to easily act on
                % neighboring boundary edges)
                ordered_boundary_vertices = boundary_vertices(1,:);
                boundary_vertices(1,:) = [];
                for i=2:size(boundary_vertices,1)+1
                    
                    % link the current head vertex to the previous tail vertex 
                    current_he = find(boundary_vertices(:,1) == ordered_boundary_vertices(i-1,2));
                    
                    % check if a new shape has to be mapped
                    if isempty(current_he)
                        current_he = 1;
                    end
                        
                    ordered_boundary_vertices(i,:) = boundary_vertices(current_he,:);
                    
                    % delete the already used vertices in case more than 1 shape 
                    boundary_vertices(current_he,:) = [];
                    
                end
                
                % Find the halfedges associated to each connection
                hes = vertcat(obj.halfedges.vertices);
                for i=1:length(ordered_boundary_vertices)
                    
                    obj.mesh_boundary(i) = find(hes(:,1)==ordered_boundary_vertices(i,1)...
                                              & hes(:,2)==ordered_boundary_vertices(i,2));
                end
                
                
                % Store remaining halfedge indices as body vertices
                obj.mesh_body = 1:obj.nH;
                obj.mesh_body(obj.mesh_boundary)=[];
                
            end % find_boundary_edges
            
            % Calculate edge lengths, midpoints, and link lengths
            function setLengths(obj)
                
                bar = waitbar(0,'Calculating edge and link lengths...');
                for he_ind = 1:obj.nH
                    waitbar(he_ind/obj.nH,bar);
                    
                    % Vertices of both ends of edge
                    v=obj.vertices(obj.halfedges(he_ind).vertices,:);
                    
                    % Edge length (from vertex to vertex)
                    obj.halfedges(he_ind).edgeLength = norm(v(1,:)-v(2,:));
                    
                    % Midpoint
                    obj.halfedges(he_ind).midpoint = (v(1,:)+v(2,:))/2;
                    
                    % Face circumcenter
                    cur_face = obj.halfedges(he_ind).face;
                    CC = circumcenter(obj.TR,cur_face);
                    obj.halfedges(he_ind).circumcenter = CC;
                    
                    % Opposite halfedge index
                    flip_he = obj.halfedges(he_ind).flip;
                    
                    % If we're not on  a boundary edge...
                    if  flip_he~=0
                        
                        % Opposite face index
                        flip_face = obj.halfedges(flip_he).face;
                        
                        % If the material is the same on both sides of the edge
                        if obj.faces(cur_face).fnum==obj.faces(flip_face).fnum
                            
                            % The link length is simply half the distance from circumcenter to circumcenter
                            NCC = circumcenter(obj.TR,flip_face);     % neighbour circumcenter
                            obj.halfedges(he_ind).linkLength = 0.5*norm(CC - NCC);
                            
                            
                        else    % If we are linking to a different material...
                            
                            % The link length is the distance from circumcenter to halfedge midpoint
                            obj.halfedges(he_ind).linkLength = norm(CC - obj.halfedges(he_ind).midpoint);
                            
                        end
                        
                    else        % If we are linking to the boundary...
                        
                        % The link length is the distance from circumcenter to halfedge midpoint
                        obj.halfedges(he_ind).linkLength = norm(CC - obj.halfedges(he_ind).midpoint);
                        
                    end
                    
                end
                close(bar);
                
                % Find shortest link length
                obj.shortestLinkLength = min(vertcat(obj.halfedges.linkLength));
                
                
                % Find triangle areas
                num_edges_per_face = obj.nV_face;
                l = zeros(obj.nF,num_edges_per_face);
                for i=1:num_edges_per_face
                    l(:,i) = vertcat(obj.halfedges(i:num_edges_per_face:obj.nH).edgeLength);
                end
                p = sum(l,2)/2;     % half perimeter
                tri_area = sqrt( p .* (p-l(:,1)) .* (p-l(:,2)) .* (p-l(:,3)) );
                for f=1:obj.nF
                    obj.faces(f).area = tri_area(f);
                end
                
            end % setLengths
            
        end % constructor
        
        
        % Find the average triangle area of the mesh
        function av_area = average_area(obj)
            
            av_area = mean(vertcat(obj.faces.area));
            
        end
        
        % Find average mesh quality factor, Q
        function Q = Q_factor(obj)
            
            QF = zeros(obj.nF,1);
            
            for f = 1:obj.nF
                
                % The halfedges associated to the triangle
                he_ind = f*obj.nV_face - obj.nV_face + 1 : f*obj.nV_face;
                
                % Any vertex associated to the triangle
                v = obj.faces(f).vertices(1);
                
                % Face circumcenter
                CC = obj.halfedges(he_ind(1)).circumcenter;
                
                % Circumradius
                CR = norm(CC - obj.vertices(v,:));
                
                % Triangle shortest edge
                SE = min(vertcat(obj.halfedges(he_ind).edgeLength));
                
                QF(f) = CR/SE;
            end
            
            Q = mean(QF);
        end
        
        function plot_mesh(obj)
            
            figure, axis equal; camlight
            
            if obj.num_materials>1
                hold on
                for f=1:obj.nF
                    material_num = vertcat(obj.faces(f).fnum);
                    hb(material_num) = patch('Vertices',obj.vertices,'Faces',vertcat(obj.faces(f).vertices),...
                        'FaceColor',obj.color(mod(material_num-1,size(obj.color,1))+1,:),'EdgeColor','black');
                end
                hold off
                
                % legend
                for i=1:numel(hb)
                    set(hb(i),'DisplayName',sprintf('Material %i',i))
                end
                legend(hb)
                
            else
                patch('Vertices',obj.vertices,'Faces',vertcat(obj.faces.vertices),...
                'FaceColor',obj.color(1,:),'EdgeColor','black')
            end
            
            % Alternatively:
            
            % plot(obj.vertices(1,:),obj.vertices(2,:),'LineWidth',3)
            
            % triplot(obj.TR);
            
            
        end % plot_mesh
        
        
        function plot_face(obj, faces)
            
            if nargin<2, faces=1:obj.nF; end
            figure; triplot(obj.TR); axis equal;
            hold on
            for f = 1:numel(faces)
                he_ind = (faces(f)-1)*obj.nV_face+1:faces(f)*obj.nV_face;
                v = vertcat(obj.halfedges(he_ind).vertices);
                patch(obj.vertices(v(:,1),1),obj.vertices(v(:,1),2), 'r')
                
                % label face
                CC = obj.halfedges(he_ind).circumcenter;
                h = text(CC(1),CC(2), [num2str(faces(f))]);
                set(h,'BackgroundColor',[1 1 .6],'VerticalAlignment','bottom')
            end
            hold off
            
        end % plot_node_location
        
        
        function plot_vertex(obj, vertices)
            
            figure;
            triplot(obj.TR); axis equal;
            hold on
            for v = 1:numel(vertices)
                plot(obj.vertices(vertices(v),1),obj.vertices(vertices(v),2),...
                    'o','MarkerSize',12,'MarkerEdgeColor','r','MarkerFaceColor','none')
                
                % label vertex
                h = text(obj.vertices(vertices(v),1),obj.vertices(vertices(v),2),...
                    [num2str(vertices(v))]);
                set(h,'BackgroundColor',[1 1 .6],'VerticalAlignment','bottom')
            end
            hold off
            
        end % plot_vertex_location
        
        
        function plot_halfedge(obj, halfedges, col)
            
            figure;
            triplot(obj.TR);
            axis equal;
            hold on
            
            if nargin<2, halfedges=1:obj.nH; end
            if nargin<3, col=1;
            else title(sprintf('Material %i',col));
            end
            
            for he_ind = 1:numel(halfedges)
                % Plot edge line
                vec = obj.vertices(obj.halfedges(halfedges(he_ind)).vertices,:);
                plot(vec(:,1),vec(:,2),'LineWidth',4, 'Color',obj.color(col,:))
                
                % Plot line for midpoint to face circumcenter
                mid = obj.halfedges(halfedges(he_ind)).midpoint;
                CC = obj.halfedges(halfedges(he_ind)).circumcenter;
                plot([mid(1) CC(1)],[mid(2) CC(2)],'-ok','MarkerSize',5)
                
                % label halfedge
                h = text((mid(1)+CC(1))/2,(mid(2)+CC(2))/2,[num2str(halfedges(he_ind))]);
                set(h,'BackgroundColor',[1 1 1])
            end
            hold off
            
            
        end % plot_vertex_location
        
        
        function plot_boundary(obj)
            
            for m=1:obj.num_materials

                obj.plot_halfedge(obj.material_boundaries{m},m);

            end
            
        end % plot_boundary
        
        
    end
    
end

