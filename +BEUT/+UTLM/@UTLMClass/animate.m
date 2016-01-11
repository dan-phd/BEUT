function animate(obj, field, scale)
% Animate E or H field of the mesh

% Scale the colors so that the fields are more obvious, this is required because
% UTLM usually has a much larger field at the excitation than everywhere else
if nargin<3
    scale = 5;
end

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


% Normalize
S.FaceVertexCData = scale*S.FaceVertexCData/max(max(max(abs(S.FaceVertexCData))));
max_amplitude = 1;
min_amplitude = -max_amplitude;

BEUT.animate_fields(2,'skipTimesteps',10,...
    [field ' field animation'],S,...
    'overlay',material_vertices,...
    'max_amplitude',max_amplitude,...
    'min_amplitude',min_amplitude);

end
