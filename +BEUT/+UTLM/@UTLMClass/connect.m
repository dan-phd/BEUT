function connect(obj, k, E)
% run the connect process for timestep k. If the electric field E is given, then we use this as the incident
% voltage at the boundary halfedges, as is stipulated by the BEUT method. If E is not given, compute the
% connect process as normal.

num_hes = numel(obj.halfedges);
V_i = zeros(num_hes,1);
E_i = zeros(num_hes,1);
H_i = zeros(num_hes,1);
Hx = zeros(num_hes,1);
Hy = zeros(num_hes,1);

% Work out incident link voltages for boundary nodes
for i=1:numel(obj.mesh_boundary)
    
    he_ind = obj.mesh_boundary(i);
    
    
    % Use the electric field calculated by BEM (3rd arg)
    if nargin>2
        
        % Output Vi and normalised current as fields at each edge
        E_i(he_ind) = E(i);
        V_i(he_ind) = E(i);
        
        
        % If this is a PEC boundary, simply reflect the signal
        if isinf(obj.halfedges(he_ind).Y_stub)
            
            H_i(he_ind) = 0;
            
            obj.halfedges(he_ind).V_linki = V_i(he_ind);
            
        else
            
            I_sources = 2 * (obj.halfedges(he_ind).V_linkr * obj.halfedges(he_ind).Y_link +...
                obj.halfedges(he_ind).V_stub * obj.halfedges(he_ind).Y_stub);
            H_i(he_ind) = I_sources;
            
            % Get new (incident) link and stub voltages - using boundary admittance
            obj.halfedges(he_ind).V_linki = V_i(he_ind) - obj.halfedges(he_ind).V_linkr;
            obj.halfedges(he_ind).V_stub = V_i(he_ind) - obj.halfedges(he_ind).V_stub;
            
        end
        
        
    % if it is an aborbing boundary condition, use matched boundary admittance
    elseif obj.reflection_coeff(he_ind)==0
        
        I_sources = 2 * (obj.halfedges(he_ind).V_linkr * obj.halfedges(he_ind).Y_link +...
            obj.halfedges(he_ind).V_stub * obj.halfedges(he_ind).Y_stub);
        
        % Use the boundary admittance to find total admittance
        total_admittance = obj.halfedges(he_ind).Y_link + obj.halfedges(he_ind).Y_stub + ...
            obj.Y_boundary(he_ind);
        
        V_i(he_ind) = I_sources / total_admittance;
        E_i(he_ind) = V_i(he_ind);
        H_i(he_ind) = I_sources / obj.halfedges(he_ind).edgeLength;
        
        % Get new (incident) link and stub voltages - using boundary admittance
        obj.halfedges(he_ind).V_linki = V_i(he_ind) - obj.halfedges(he_ind).V_linkr;
        obj.halfedges(he_ind).V_stub = V_i(he_ind) - obj.halfedges(he_ind).V_stub;
        
        
    % else just use the boundary coefficient
    else
        
        I_sources = 2 * (obj.halfedges(he_ind).V_linkr * obj.halfedges(he_ind).Y_link +...
            obj.halfedges(he_ind).V_stub * obj.halfedges(he_ind).Y_stub);
        
        total_admittance = obj.halfedges(he_ind).Y_link + obj.halfedges(he_ind).Y_stub;
        
        % Output Vi and normalised current as fields at each edge
        V_i(he_ind) = I_sources / total_admittance;
        E_i(he_ind) = V_i(he_ind);
        H_i(he_ind) = I_sources / obj.halfedges(he_ind).edgeLength;
        
        % Get new (incident) link and stub voltages - using reflection coefficient
        obj.halfedges(he_ind).V_linki = obj.halfedges(he_ind).V_linkr * obj.reflection_coeff(he_ind);
        obj.halfedges(he_ind).V_stub = obj.halfedges(he_ind).V_stub * obj.reflection_coeff(he_ind);
        
    end
    
end

% Work out incident link voltages for remaining PEC boundaries (not on domain boundary)
for i=1:numel(obj.PEC_boundary)
    
    he_ind = obj.PEC_boundary(i);
    
    I_sources = 2 * (obj.halfedges(he_ind).V_linkr * obj.halfedges(he_ind).Y_link +...
        obj.halfedges(he_ind).V_stub * obj.halfedges(he_ind).Y_stub);
    
    total_admittance = obj.halfedges(he_ind).Y_link + obj.halfedges(he_ind).Y_stub;
    
    % Output Vi and normalised current as fields at each edge
    V_i(he_ind) = I_sources / total_admittance;
    E_i(he_ind) = V_i(he_ind);
    H_i(he_ind) = I_sources / obj.halfedges(he_ind).edgeLength;
    
    % Get new (incident) link and stub voltages - using reflection coefficient
    obj.halfedges(he_ind).V_linki = -obj.halfedges(he_ind).V_linkr;
    obj.halfedges(he_ind).V_stub = -obj.halfedges(he_ind).V_stub;
        
    
end

% Work out incident link voltages for internal nodes
for i=1:numel(obj.mesh_body)
    
    he_ind = obj.mesh_body(i);
    flip_ind = obj.halfedges(he_ind).flip;        % neighbour halfedge index
    
    % If this node is part of a PEC, do not scatter
    if isinf(obj.halfedges(he_ind).Y_stub) || isinf(obj.halfedges(flip_ind).Y_stub)
        
        % Output Vi and normalised current as fields at each edge
        E_i(he_ind) = 0;
        H_i(he_ind) = 0;
        Hx(he_ind) = 0;
        Hy(he_ind) = 0;
        
        % Get new (incident) link and stub voltages
        obj.halfedges(he_ind).V_linki = 0;
        obj.halfedges(he_ind).V_stub = 0;
        obj.halfedges(flip_ind).V_linki = 0;
        obj.halfedges(flip_ind).V_stub = 0;
        
    else
        
        % Remember that the voltage doubles when a signal encounters an open cct on the TL
        I_sources = 2 * (obj.halfedges(he_ind).V_linkr * obj.halfedges(he_ind).Y_link +...
            obj.halfedges(flip_ind).V_linkr * obj.halfedges(flip_ind).Y_link +...
            obj.halfedges(flip_ind).V_stub * obj.halfedges(flip_ind).Y_stub +...
            obj.halfedges(he_ind).V_stub * obj.halfedges(he_ind).Y_stub);
        
        total_admittance = obj.halfedges(he_ind).Y_link + obj.halfedges(flip_ind).Y_link +...
            obj.halfedges(he_ind).Y_stub + obj.halfedges(flip_ind).Y_stub;
        
        
        % Output Vi and normalised current as fields at each edge
        V_i(he_ind) = I_sources / total_admittance;
        E_i(he_ind) = V_i(he_ind);% / obj.halfedges(he_ind).linkLength;
        H_i(he_ind) = I_sources / obj.halfedges(he_ind).edgeLength;
        
        % Get new (incident) link and stub voltages
        obj.halfedges(he_ind).V_linki = V_i(he_ind) - obj.halfedges(he_ind).V_linkr;
        obj.halfedges(he_ind).V_stub = V_i(he_ind) - obj.halfedges(he_ind).V_stub;
        obj.halfedges(flip_ind).V_linki = V_i(he_ind) - obj.halfedges(flip_ind).V_linkr;
        obj.halfedges(flip_ind).V_stub = V_i(he_ind) - obj.halfedges(flip_ind).V_stub;
        
    end
    
    % To find x and y components of H_xy, we must find the angle between the edge vector and
    % the x-axis:
    v=obj.vertices(obj.halfedges(he_ind).vertices,:);
    vec = v(2,:)-v(1,:);
    % vectarrow([0 0],vec)        % plot vector
    angle = atan(vec(2)/vec(1));
    Hx(he_ind) = H_i(he_ind)*sin(angle);
    Hy(he_ind) = H_i(he_ind)*cos(angle);
    
    % We need to remember to output the flip halfedge values too
    E_i(flip_ind) = E_i(he_ind);
    H_i(flip_ind) = H_i(he_ind);
    Hx(flip_ind) = Hx(he_ind);
    Hy(flip_ind) = Hy(he_ind);
    
end

% Store fields
obj.fields(1).E_z(:,k) = E_i;
obj.fields(1).H_xy(:,k) = H_i;
obj.fields(1).H_x(:,k) = Hx;
obj.fields(1).H_y(:,k) = Hy;

end % connect
