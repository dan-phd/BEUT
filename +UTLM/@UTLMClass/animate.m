function animate(obj, field)
% Animate E or H field of the mesh

% Set plot properties
S.Vertices = obj.vertices;
S.Faces = vertcat(obj.faces.vertices);
S.FaceColor = 'flat';
S.LineStyle = 'None';

if nargin<2, field='E'; end
if strcmp(field,'E')
    S.FaceVertexCData = obj.V0;
elseif strcmp(field,'H')
    S.FaceVertexCData = obj.I0;
else
    error('field must be either ''E'' or'' H''');
end

% find vertices around materials
material_vertices = unique(...
    obj.vertices(...
    vertcat(...
    obj.halfedges(...
    horzcat(obj.material_boundaries{:})...
    ).vertices...
    ),:...
    ), 'rows');

BEUT.animate_fields(2,'skipTimesteps',10,...
    [field ' field animation'],S,...
    'overlay',material_vertices);

end % animate_mesh
