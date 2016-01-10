function scatter(obj, k)

% Copy object variables to local workspace (for speed benefit)
num_edges_per_face = obj.nV_face;
num_hes = numel(obj.halfedges);
Z_0 = vertcat(obj.faces.Z0);

% Work out voltage at centre of node
I_0 = zeros(obj.nF,num_edges_per_face);
for i = 1:num_edges_per_face
    I_0(:,i) = 2*vertcat(obj.halfedges(i:num_edges_per_face:num_hes).V_linki) .* ...
        vertcat(obj.halfedges(i:num_edges_per_face:num_hes).Y_link);
end
I_0 = sum(I_0,2);

V_0 = I_0 .* Z_0;


% Work out reflected link voltages (Vr=V-Vi)
for he_ind = 1:num_hes
    
    f = obj.halfedges(he_ind).face;         % face index
    
    % If this is a PEC node, simply reflect the signal
    if isinf(obj.faces(f).eps_r)
        
        V_0(f) = 0; I_0(f) = 0;
        obj.halfedges(he_ind).V_linkr = - obj.halfedges(he_ind).V_linki;
        
    else
    
        obj.halfedges(he_ind).V_linkr = V_0(f) - obj.halfedges(he_ind).V_linki;
        
    end
    
end


% Store voltage and currents in object
obj.I0(:,k) = I_0;
obj.V0(:,k) = V_0;


% Find open circuit voltage and closed circuit current
obj.V_open = zeros(numel(obj.mesh_boundary),1);
for i=1:numel(obj.mesh_boundary)
    
    he_ind = obj.mesh_boundary(i);
    
    % If this is a PEC node, simply reflect the signal
    if isinf(obj.halfedges(he_ind).Y_stub)
        
        obj.I_closed(i) =  0;
        obj.V_open(i) = obj.halfedges(he_ind).V_linkr;
        
    else
        
        Vl = obj.halfedges(he_ind).V_linkr;
        Vs = obj.halfedges(he_ind).V_stub;
        Yl = obj.halfedges(he_ind).Y_link;
        Ys = obj.halfedges(he_ind).Y_stub;
        obj.I_closed(i) =  2*Vl*Yl + 2*Vs*Ys;
        obj.V_open(i) = ( (2*Vl*Yl + 2*Vs*Ys)/(Ys+Yl) );
        
    end
    
end

end % scatter
