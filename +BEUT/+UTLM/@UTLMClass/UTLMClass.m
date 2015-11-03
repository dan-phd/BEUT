classdef UTLMClass < BEUT.Meshing.HalfedgeMesh & handle
    %Read and store unstructured meshes
    
    properties
        % vertices, faces and halfedge properties defined
        % in halfedge_mesh superclass
        
        % fields at the edges
        fields = struct('E_z',{},'H_xy',{},'H_x',{},'H_y',{});
        V0; I0;
        
        dt;         % timestep
        
        reflection_coeff;
        Y_boundary;
        
        V_open; I_open;
    end
    
    methods
        
        % Input vertices and faces of geometry
        function obj = UTLMClass(vertices, faces, varargin)
            
            obj = obj@BEUT.Meshing.HalfedgeMesh(vertices, faces, varargin{:});
            
            obj.fields(1).E_z = zeros(obj.nH,1);
            obj.fields(1).H_xy = zeros(obj.nH,1);
            obj.fields(1).H_x = zeros(obj.nH,1);
            obj.fields(1).H_y = zeros(obj.nH,1);
            obj.I0 = zeros(obj.nF,1);
            obj.V0 = zeros(obj.nF,1);
            
            % Initialise voltages and connect flag
            for i=1:obj.nH
                
                obj.halfedges(i).V_linki = 0;
                obj.halfedges(i).V_linkr = 0;
                obj.halfedges(i).V_stub = 0;
                obj.halfedges(i).doConnect = 1;
            
            end
            
        end % constructor
        
        
        
        % Excite the source halfedge as an electric field source
        function excite_E(obj, sourceEdges, V_source)
            
            for edge = 1:numel(sourceEdges)
                
                % current and neighbour halfedge index
                he_ind = sourceEdges(edge);
                flip_ind = obj.halfedges(he_ind).flip;
                
                if flip_ind==0
                    total_admittance = obj.halfedges(he_ind).Y_link + obj.halfedges(he_ind).Y_stub;
                else
                    total_admittance = obj.halfedges(he_ind).Y_link + obj.halfedges(flip_ind).Y_link +...
                        obj.halfedges(he_ind).Y_stub + obj.halfedges(flip_ind).Y_stub;
                end
                
                % Add source excitation
                obj.halfedges(he_ind).V_linki = obj.halfedges(he_ind).V_linki + V_source/total_admittance;
                obj.halfedges(he_ind).V_stub = obj.halfedges(he_ind).V_stub + V_source/total_admittance;
                
            end
            
        end % excite_E
        
        
        % Excite the source halfedge as a magentic field source
        function excite_H(obj, sourceEdges, V_source)
            
            for edge = 1:numel(sourceEdges)
                
                % current halfedge index
                he_ind = sourceEdges(edge);
                
                % Add source excitation
                obj.halfedges(he_ind).V_linki = obj.halfedges(he_ind).V_linki + V_source*obj.halfedges(he_ind).Y_link;
                obj.halfedges(he_ind).V_stub = obj.halfedges(he_ind).V_stub + V_source*obj.halfedges(he_ind).Y_stub;
                
            end
            
        end % excite_H
        
        
        % Excite  the source halfedge and scale by edge_length/Y (J source)
        function excite_J(obj, sourceEdges, V_source)
            
            for edge = 1:numel(sourceEdges)
                
                % current halfedge index
                he_ind = sourceEdges(edge);
                
                % Add source excitation
                obj.halfedges(he_ind).V_linki = obj.halfedges(he_ind).V_linki + V_source * ...
                    obj.halfedges(he_ind).edgeLength/obj.halfedges(he_ind).Y_link;
                
            end
            
        end % excite_J
        
        
    end % methods
    
end
