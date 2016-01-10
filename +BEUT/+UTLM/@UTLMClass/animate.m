function animate(obj, field)
% Animate E or H field of the mesh

% Set plot properties
S.Vertices = obj.vertices;
S.Faces = vertcat(obj.faces.vertices);
S.LineStyle = 'None';

if nargin<2, field='E'; end
if strcmp(field,'E')
    
    % interpolate across faces by averaging halfedge data across vertices:
%     S.FaceColor = 'interp';
%     vData = zeros(length(obj.vertices),size(obj.V0,2));
%     numel_vData = zeros(length(obj.vertices),1);
%     for he=1:obj.nH
%         verts = obj.halfedges(he).vertices;
%         for i=1:2
%             vData(verts(i),:)=vData(verts(i),:)+obj.fields.E_z(he,:);
%             numel_vData(verts(i)) = numel_vData(verts(i))+1;
%         end
%     end
%     % averaging
%     for v=1:length(obj.vertices)
%         vData(v,:)=vData(v,:)/numel_vData(v);
%     end
%     S.FaceVertexCData = vData;
    
    % Or just use face values:
    S.FaceColor = 'flat';
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

max_amplitude = max(max(max(S.FaceVertexCData)));
min_amplitude = min(min(min(S.FaceVertexCData)));
% max_amplitude = 1;
% min_amplitude = -max_amplitude;

S.FaceVertexCData = 6*S.FaceVertexCData/max(max(max(abs(S.FaceVertexCData))));

BEUT.animate_fields(2,'skipTimesteps',10,...
    [field ' field animation'],S,...
    'overlay',material_vertices,...
    'max_amplitude',max_amplitude,...
    'min_amplitude',min_amplitude);

end
