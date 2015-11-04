function setMaterial(obj, relative_eps, relative_mu)
% Set material of UTLM mesh

% find number of materials
num_materials_eps = numel(relative_eps);
num_materials_mu = numel(relative_mu);
num_materials = max(num_materials_eps,num_materials_mu);
assert(num_materials_eps+num_materials_mu==2*num_materials ||...
    num_materials_eps+num_materials_mu==num_materials+1,...
    'number of material parameters must equal each other or equal 1');
assert(num_materials==obj.num_materials || num_materials==obj.nF,...
    'number of materials must be equal to the number of faces or materials specified in the mesh');

% retreive material parameters
if num_materials_eps==1
    if num_materials_mu==1
        for i=1:num_materials
            eps_r(i) = relative_eps;
            mu_r(i) = relative_mu;
        end
    else
        for i=1:num_materials
            eps_r(i) = relative_eps;
            mu_r(i) = relative_mu(i);
        end
    end
else
    if num_materials_mu==1
        for i=1:num_materials
            eps_r(i) = relative_eps(i);
            mu_r(i) = relative_mu;
        end
    else
        for i=1:num_materials
            eps_r(i) = relative_eps(i);
            mu_r(i) = relative_mu(i);
        end
    end
end


% set each face's material parameters
if num_materials>obj.num_materials
    
    % set each face individually
    for f = 1:obj.nF
        
        obj.faces(f).eps_r = eps_r(f);
        obj.faces(f).mu_r = mu_r(f);
        
    end
    
else
    
    % set each material
    for f = 1:obj.nF
        
        obj.faces(f).eps_r = eps_r(obj.faces(f).fnum);
        obj.faces(f).mu_r = mu_r(obj.faces(f).fnum);
        
    end
    
end

end % setMaterial
