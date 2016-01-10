function setBoundary(obj, condition)
% Set boundary condition, condition can be 1(=open), -1(=short), 0(=absorbing),
% or any other number (to equal the desired reflection coefficient)


% Number of halfedges on boundary
num_boundary_hes = numel(obj.mesh_boundary);

if (numel(condition)==1)
    condition = ones(num_boundary_hes,1) * condition;
else
    assert(numel(condition)==num_boundary_hes,...
        'Each boundary halfedge must have a boundary condition')
end

num_hes = numel(obj.halfedges);
obj.reflection_coeff = zeros(num_hes,1);
obj.Y_boundary = zeros(num_hes,1);

% Find boundary halfedge indices
boundary_ind = obj.mesh_boundary;

% Set reflection coefficient for each boundary edge
for i=1:num_boundary_hes
    
    switch condition(i)
        
        % open
        case 1
            obj.reflection_coeff(boundary_ind(i)) = 1;
            obj.Y_boundary(boundary_ind(i)) = 0;
            
        % short
        case -1
            obj.reflection_coeff(boundary_ind(i)) = -1;
            obj.Y_boundary(boundary_ind(i)) = 0;
            
        % absorbing
        case 0
            obj.reflection_coeff(boundary_ind(i)) = 0;
            
            global mu0 eps0;
            Y0 = sqrt(eps0/mu0);
            
            he_ind = boundary_ind(i);
            obj.Y_boundary(he_ind) = Y0*sqrt(obj.faces(obj.halfedges(he_ind).face).eps_r) *...
                obj.halfedges(he_ind).edgeLength;% / sqrt(2);
            
        otherwise
            obj.reflection_coeff(boundary_ind(i)) = condition(i);
            
    end
end

end
