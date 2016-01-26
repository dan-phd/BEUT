function plot_materials(obj, material_parameter)
% Plot the mesh with the color of each triangle equal to the values of permittivity (using material_parameter
% = 'eps_r') or permeability (using material_parameter = 'mu_r').

if strcmp(material_parameter,'eps_r')
    material_param = 1;
    material_name = '$\varepsilon_r$';
elseif strcmp(material_parameter,'mu_r')
    material_param = 2;
    material_name = '$\mu_r$';
else
    error('material_parameter must be either ''eps_r'' or ''mu_r''');
end

% Set plot properties
S.Vertices = obj.vertices;
S.Faces = vertcat(obj.faces.vertices);
S.FaceColor = 'flat';
S.LineStyle = 'None';
figure; axis equal; hold on;
colorbar
title([material_name ' of mesh'],'Interpreter','LaTex')

% Plot each face with color equal to material parameter
if material_param == 1
    S.FaceVertexCData = vertcat(obj.faces.eps_r);
elseif material_param == 2
    S.FaceVertexCData = vertcat(obj.faces.mu_r);
end

patch(S)
axis tight;

xlabel('x (m)');
ylabel('y (m)');

end
