classdef MeshBoundary < handle
    % Extract boundary edges of a halfedge mesh
%     
%     For example:
% 
%       [v,f]=BEUT.Meshing.loadMesh([fileparts(which('BEUT.Meshing.loadMesh'))...
%             filesep 'meshes' filesep 'cyl_res21.mat']);
%       geometry=BEUT.Meshing.MeshBoundary(BEUT.Meshing.HalfedgeMesh(v,f));
%       geometry.plot
    
    properties
        
        % a is the vertex at the head of the edge
        % b is the vertex at the tail of the edge
        % l is the length of the edge
        % t is the vector tangent to the edge
        % n is the vector normal to the edge
        % shape is the shape number that this edge belongs to
        halfedges = struct('a',{},'b',{},'l',{},'t',{},'n',{},'shape',{});
        dual = struct('a',{},'b',{},'l',{},'t',{},'n',{},'shape',{});
        
        N_V;         % number of vertices (equal to number of halfedges for closed surface)
        num_shapes;
    end
    
    methods
        
        % Constructor to define vertices
        function obj = MeshBoundary(halfedge_mesh)
            
            boundary_halfedges = halfedge_mesh.mesh_boundary;
            obj.N_V = size(boundary_halfedges,2);
            boundary_vertices = vertcat(halfedge_mesh.halfedges(boundary_halfedges).vertices);
            
            % Vertex head and tail
            halfedge_tbl = table;
            halfedge_tbl.he_idx = boundary_halfedges';
            halfedge_tbl.a = halfedge_mesh.vertices(boundary_vertices(:,1),:);
            halfedge_tbl.b = halfedge_mesh.vertices(boundary_vertices(:,2),:);
            obj.halfedges = table2struct(halfedge_tbl);
            
            % Edge lengths
            [obj.halfedges.l]=deal(halfedge_mesh.halfedges(boundary_halfedges).edgeLength);
            
            % Tangent and normal unit vectors
            for he=1:obj.N_V
                
                % Tangent unit vector
                obj.halfedges(he).t = (obj.halfedges(he).b - obj.halfedges(he).a) / obj.halfedges(he).l;
                
                % normal = [dy, -dx] = [y2-y1, -x2+x1]
                obj.halfedges(he).n = [obj.halfedges(he).b(2)-obj.halfedges(he).a(2),...
                                      -obj.halfedges(he).b(1)+obj.halfedges(he).a(1)] / obj.halfedges(he).l;
            end
            
            % Determine spatially distint shape boundaries
            shape_count = 1; he_found_for_this_material = 0;
            for i=1:halfedge_mesh.num_materials
                for he=1:obj.N_V
                    
                    if ~isempty(find(halfedge_mesh.material_boundaries{i}==boundary_halfedges(he), 1))
                        he_found_for_this_material = 1;
                        obj.halfedges(he).shape = shape_count;
                    end
                end
                
                % increment number of shapes found so far
                if he_found_for_this_material
                    shape_count = shape_count+1;
                    he_found_for_this_material = 0;
                end
            end
            
            obj.num_shapes = shape_count-1;
            
            
            % Create dual mesh by splitting each edge in 2
            dual_mesh = obj.halfedges;
            for i=1:obj.N_V
                
                % split this edge in 2
                a1 = obj.halfedges(i).a;
                b2 = obj.halfedges(i).b;
                b1 = 0.5 * (a1 + b2);
                a2 = b1;
                l = 0.5 * obj.halfedges(i).l;
                
                % first half
                idx = i*2-1;
                dual_mesh(idx).he_idx = obj.halfedges(i).he_idx;
                dual_mesh(idx).a = a1;
                dual_mesh(idx).b = b1;
                dual_mesh(idx).l = l;
                dual_mesh(idx).t = obj.halfedges(i).t;
                dual_mesh(idx).n = obj.halfedges(i).n;
                dual_mesh(idx).shape = obj.halfedges(i).shape;
                
                % second half
                idx = i*2;
                dual_mesh(idx).he_idx = obj.halfedges(i).he_idx;
                dual_mesh(idx).a = a2;
                dual_mesh(idx).b = b2;
                dual_mesh(idx).l = l;
                dual_mesh(idx).t = obj.halfedges(i).t;
                dual_mesh(idx).n = obj.halfedges(i).n;
                dual_mesh(idx).shape = obj.halfedges(i).shape;
                
                
                obj.dual = dual_mesh;
                
            end
            
        end % constructor
        
        
        % Plot shape and label halfedges
        function plot(obj, edges_to_label)
            
            if nargin<2, edges_to_label=1:obj.N_V; end
            
            vertices = vertcat(obj.halfedges);
            a = vertcat(vertices.a);
            b = vertcat(vertices.b);
                        
            figure; axis equal;
            hold on
            % Plot edge lines
            for he_ind = 1:obj.N_V
                
                vecx = [a(he_ind,1), b(he_ind,1)];
                vecy = [a(he_ind,2), b(he_ind,2)];
                plot(vecx,vecy,'LineWidth',4)
                
            end
            % label halfedges at midpoints
            for he_ind = edges_to_label
                
                mid = ( a(he_ind,:)+b(he_ind,:) )/2;
                h = text(mid(1),mid(2),num2str(he_ind));
                set(h,'BackgroundColor',[1 1 1])
                
            end
            hold off
            
            
        end % plot
        
        
        % Stability check: largest link length < lambda/10
        function [] = checkStability(obj, lambda)
            
            longest_l = max(vertcat(obj.halfedges.l));
            stable = (lambda/10)/longest_l;        % will be decimal if unstable
            if stable < 1
                error('The frequency and/or number of vertices chosen will lead to instability')
            end
            
        end
        
    end
    
end

