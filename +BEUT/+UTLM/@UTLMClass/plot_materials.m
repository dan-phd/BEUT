function plot_materials(obj, material_parameter)
% Plot mesh with materials

if strcmp(material_parameter,'eps_r')
    material_param = 1;
elseif strcmp(material_parameter,'mu_r')
    material_param = 2;
else
    error('material_parameter must be either ''eps_r'' or'' mu_r''');
end

% Set plot properties
S.Vertices = obj.vertices;
S.Faces = vertcat(obj.faces.vertices);
S.FaceColor = 'flat';
S.LineStyle = 'None';
figure; axis equal; hold on;
colorbar
title([material_parameter ' of mesh'])

% Plot each face with color equal to material parameter
if material_param == 1
    S.FaceVertexCData = vertcat(obj.faces.eps_r);
elseif material_param == 2
    S.FaceVertexCData = vertcat(obj.faces.mu_r);
end

patch(S)

end
