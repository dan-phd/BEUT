function calcAdmittances(obj)
% Calculate admittances and set connect flags


global mu0 eps0;

for f = 1:obj.nF
    
    Y_link_total = 0;
    
    for edge = 1:obj.nV_face
        
        % halfedge index
        he_ind = (f-1)*obj.nV_face+edge;
        
        % Admittance calculations
        Y_L_link = (obj.halfedges(he_ind).edgeLength* obj.dt)/...
            (obj.halfedges(he_ind).linkLength * mu0);
        Y_C_link = (obj.halfedges(he_ind).edgeLength *...
            obj.halfedges(he_ind).linkLength * eps0)/(2*obj.dt);
        Y_L_stub = (obj.halfedges(he_ind).edgeLength* obj.dt)/...
            (2*obj.halfedges(he_ind).linkLength * mu0);
        Y_C_stub = (obj.halfedges(he_ind).edgeLength *...
            obj.halfedges(he_ind).linkLength * eps0 * obj.faces(f).eps_r)/obj.dt;
        
        
        obj.halfedges(he_ind).Y_link = Y_L_link/2;
        obj.halfedges(he_ind).Y_stub = Y_C_stub - Y_L_stub;
        
        
        % Total admittance for face, f
        Y_link_total = Y_link_total + obj.halfedges(he_ind).Y_link;
        
        
        % set connect flags (so we don't loop over and rewrite V_link twice)
        flip_ind = obj.halfedges(he_ind).flip;        % neighbour halfedge index
        if obj.halfedges(he_ind).doConnect && flip_ind~=0
            obj.halfedges(flip_ind).doConnect = 0;
            obj.mesh_body(obj.mesh_body==flip_ind) = [];
        end
        
    end
    
    % Total short circuit impedance for current face
    obj.faces(f).Z0 = 1/Y_link_total;
    
end

end % set_admittances
