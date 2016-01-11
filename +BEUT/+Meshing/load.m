function [ vertices, faces, fnum ] = load(filename, type)
% Find vertices, faces, and face numbers (if simulating different
% materials) using input file
%
%  [ vertices, faces, fnum ] = load_mesh(filename, type)
%  filename: char array with file to be read.
%  type:     optional, indicates file type.
%
% Supported types are:
%  .in
%  .gid
%  .gmsh (or .msh)
%  .mphtxt (Comsol output)
%  .obj
%  .node with associated .ele
%  .mat (Matlab triangulation object)
%  .poly
%
% example:
%
%   mesh_name = 'cyl_res21.mat';
%   [v,f]=BEUT.Meshing.loadMesh([fileparts(which('BEUT.Meshing.loadMesh')) filesep 'meshes' filesep mesh_name]);
%   patch('Vertices',v,'Faces',f,'FaceColor','red','EdgeColor','black'); axis equal

% Last edited 04/08/2014 - Daniel Simmons - dansphd.com

fid = fopen(filename,'r');
assert(fid~=0&&fid~=-1,'Could not open %s for reading',filename);

if nargin < 2
    pos = strfind(filename,'.');
    if isempty(pos);
        error('Could not determine file type');
    end
    pos = pos(end);
    type = filename(pos+1:end);
end

switch type
    case 'in'
        readInMesh(fid);
    case 'gid'
        readGIDMesh(fid);
    case {'gmsh', 'msh'}
        readGmshMesh(fid);
    case 'mphtxt'
        read_mphtxt(fid);
    case 'obj'
        OBJ_out=read_wobj(filename);
        vertices=OBJ_out.vertices;
        faces=OBJ_out.objects(end).data.vertices;
    case {'node', 'ele'}
        file_without_ext = filename(1:pos-1);
        readNodeMesh(file_without_ext);
    case {'poly'}
        readPolyMesh(fid);
    case 'mat'
        obj_out = load(filename);
        fnum = vertcat(obj_out.mesh.faces.fnum);
        vertices = obj_out.mesh.TR.Points;
        faces = obj_out.mesh.TR.ConnectivityList;
    otherwise
        error('type argument invalid.');
end
fclose(fid);    % close the file

% If the material is the same throughout, fnum is all ones
if ~exist('fnum','var')
    fnum = ones(length(faces),1);
end

    function readPolyMesh(fid)
        
        % nodes for p (points by coordinates)
        N_V = fscanf(fid,'%d',1);               % number of vertices
        num_dim = fscanf(fid,'%d',1);
        mat_param = fscanf(fid,'%d',1);
        
        % Read in data
        data1=zeros(N_V,num_dim+1);
        for i=1:N_V
            data1(i,1)  = fscanf(fid, '%d',1);
            for j=1:num_dim
                temp  = fscanf(fid, '%s',1);        % must read in like this because of the "/" used in file
                data1(i,j+1) = str2num(temp);
            end
        end
        data1=data1(:,2:4);           % ignore first column
        
        % convert to 2D by extracting just the base layer
        if num_dim==3
            base_layer = find(data1(:,3)==0);
            data1 = data1(:,1:2);           % ignore 3rd dimension column
%             data1 = data1(base_layer,:);    % ignore vertices on other layers
        end
        
        vertices = data1*1e-10;
        
        
        % elements for t (tetras/triangles by points)
        N_P = fscanf(fid,'%d',1);                   % total number of polygons
        num_holes = fscanf(fid,'%d',1);
        
        % Read in data
        data2=zeros(N_P,num_dim+1);
        i=0;
        while i<N_P
            num_poly = fscanf(fid,'%d',1);      % number of polygons to read on next line
            fscanf(fid, '%d',1);              % ignore rest of line
            layer(i+1) = fscanf(fid, '%d',1);
            for j=1:num_poly
                num_vertices = fscanf(fid,'%d',1);      % number of vertices on this line
                for k=1:num_vertices
                    data2(j+i,k)=fscanf(fid,'%d',1);
                end
            end
            i=i+num_poly;
        end
        
        % indices start from 1 in Matlab (instead of 0)
        data2=data2+1;
        
        % convert to 2D
        if num_dim==3
            
            % At this point, there will be unused vertices, use matlab delaunay mesher to sort this out
            TRI = delaunayTriangulation(data1*1e-10);
            vertices = TRI.Points;
            faces = TRI.ConnectivityList;
        else
            vertices = data1*1e-10;     % convert from long int to double
            faces = data2;
        end
        
        
    end % readPolyMesh


    function readInMesh(fid)
        % Read .in mesh
        %
        % 1st line:
        %   1st variable: number of shapes
        % 2nd line:
        %   1st variable: number of vertices, NV
        %   2nd variable: number of edges, NE
        % Line 3 to NV+2
        %   x y z-coordinates
        % Line NV+3 to NV+NE+2
        %   vertex_index connected_vertex_index
        
        N_shapes = fscanf(fid,'%d',1);
        
        % most commonly number of shapes is 1
        if N_shapes==1
            NV = fscanf(fid,'%d',1);
            NE = fscanf(fid,'%d',1);
            
            v = fscanf(fid,'%f\t%f\t%f\n',[3,NV]);  % 3 coordinates per vertex
            % f = fscanf(fid,'%d\t%d\n',[3,NE]);      % 3 vertices per face
            e = fscanf(fid,'%d\t%d\n',[2,NE]);      % 2 vertices per edge
            
            j=1;
            for i=1:size(e,2)-1
                edge = e(:,i)';
                next_edge = e(:,i+1)';
                f(j,:) = [edge next_edge(2)];
                j=j+1;
            end
            vertices=v';
            faces=f;
        else
            
            for k=1:N_shapes
                NV = fscanf(fid,'%d',1);
                NE = fscanf(fid,'%d',1);
                
                v = fscanf(fid,'%f\t%f\t%f\n',[3,NV]);
                f = fscanf(fid,'%d\t%d\n',[2,NE]);
                
                obj.vertices(k).coords=v';
                obj.faces(k).vertices=f';
            end
        end
        
    end % readInMesh


    function readGIDMesh(fid)
        
        NV = 0; % number of vertices
        NF = 0; % number of faces
        
        previousLine = ''; % starting value
        nextLine = strtrim(fgetl(fid));
        
        goOn = 1;
        while goOn==1
            if strcmp(nextLine,'end coordinates')
                
                NV = sscanf(previousLine,'%d',1);
            end
            
            if strcmp(nextLine,'end elements')
                NF = sscanf(previousLine,'%d',1);
                break;
            end
            
            previousLine = nextLine;
            nextLine = strtrim(fgetl(fid));
        end
        
        %back up!
        fseek(fid,0,'bof');
        
        % read until Coordinates
        while goOn==1
            if strcmp(nextLine,'Coordinates')
                break
            end
            nextLine = strtrim(fgetl(fid));
        end
        
        % fill up the vertices
        tempV = fscanf(fid,'%g',[4 NV]);
        vertices = transpose(tempV(2:4,:));
        
        % read until Elements
        while goOn==1
            if strcmp(nextLine,'Elements')
                break
            end
            nextLine = strtrim(fgetl(fid));
        end
        
        % fill up the vertices
        tempV = fscanf(fid,'%g',[4 NF]);
        faces = transpose(tempV(2:4,:));
        
    end % readGIDMesh


    function readGmshMesh(fid)
        
        %previousLine = '';
        thisLine = strtrim(fgetl(fid));
        
        % find beginning of vertex section
        while ~strcmp(thisLine,'$Nodes')
            thisLine = strtrim(fgetl(fid));
        end
        
        % read the number of vertices
        thisLine = strtrim(fgetl(fid));
        NV = str2double(thisLine);
        
        % read vertices
        V = transpose( fscanf(fid, '%f', [4, NV]) );
        V = V(:,2:end);
        
        % find the beginning of the faces section
        while ~strcmp(thisLine,'$Elements')
            thisLine = strtrim(fgetl(fid));
        end
        
        % read the number of faces
        thisLine = strtrim(fgetl(fid));
        NF = str2double(thisLine);
        max_cols = 10;
        F = zeros(NF,10);
        
        % read the faces
        for i=1:NF
            thisLine = strtrim(fgetl(fid));
            elements = strread(thisLine);
            F(i,1:numel(elements)) = elements;
        end
        
        % trim unnecessary columns, we only need triangles (type 2)
        ele_type = F(:,2);
        triangle_idx = ele_type==2;
        F = F(triangle_idx,:);
        z = all(F==0);
        i=max_cols+1;
        endcol=0;
        while endcol==0
            i=i-1;
            if z(i)==0
                endcol=i;
            end
        end
        F = F( :, endcol-2:endcol);
        
        vertices = V;
        faces = F;
        
    end % readGmshMesh

    function read_mphtxt(fid)
        
        thisLine = strtrim(fgetl(fid));
        
        % find number of vertices
        while isempty(strfind(thisLine,'# number of mesh points'))
            thisLine = strtrim(fgetl(fid));
        end
        splitLine = strsplit(thisLine,' # ');
        NV = str2double(splitLine{1});
        
        % find starting index
        while isempty(strfind(thisLine,'# lowest mesh point index'))
            thisLine = strtrim(fgetl(fid));
        end
        splitLine = strsplit(thisLine,' # ');
        start_idx = str2double(splitLine{1});
        
        % read vertices
        while isempty(strfind(thisLine,'# Mesh point coordinates'))
            thisLine = strtrim(fgetl(fid));
        end
        V = transpose( fscanf(fid, '%f', [2, NV]) );
        
        % find number of faces
        while isempty(strfind(thisLine,'3 # number of nodes per element'))
            thisLine = strtrim(fgetl(fid));
        end
        thisLine = strtrim(fgetl(fid));
        splitLine = strsplit(thisLine,' # ');
        NF = str2double(splitLine{1});
        
        % read faces
        while isempty(strfind(thisLine,'# Elements'))
            thisLine = strtrim(fgetl(fid));
        end
        F = transpose( fscanf(fid, '%f', [3, NF]) );
        if start_idx~=1
            F = F + 1 - start_idx;
        end
        
        % find shape indices
        while isempty(strfind(thisLine,strcat(num2str(NF),' # number of geometric entity indices')))
            thisLine = strtrim(fgetl(fid));
        end
        
        % read faces
        while isempty(strfind(thisLine,'# Geometric entity indices'))
            thisLine = strtrim(fgetl(fid));
        end
        
        fnum = transpose( fscanf(fid, '%f', [1, NF]) );
        vertices = V;
        faces = F;
        
    end % readGmshMesh


    function readNodeMesh(file)
        
        % nodes for p (points by coordinates)
        node_file=fopen([file,'.node'], 'rt');         %open the file
        title = fgetl(node_file);                      %read in the header
        title_num = str2num(title);
        cols  = title_num(1);                          % 1st value encountered is number of vertices
        rows  = title_num(2) + title_num(4) + 1;       % 2nd value is dimension
        % 3rd value is # of attributes, 4th value is # of boundary markers (0 or 1)
        data  = fscanf(node_file, '%e',[rows,cols]);   %read in data
        p     = data(2:rows,1:cols);                   %ignore the first row
        fclose(node_file);                             %close the file
        
        % elements for t (tetras by points)
        ele_file=fopen([file,'.ele'], 'rt');          %open the file
        title = fgetl(ele_file);                      %read in the header
        title_num = str2num(title);
        cols  = title_num(1);
        rows  = title_num(2) + title_num(3) + 1;
        data  = fscanf(ele_file, '%i',[rows,cols]);   %read in data
        t     = data(2:rows,1:cols);                  %ignore the first row
        fclose(ele_file);                             %close the file
        
        % convert to object
        vertices = p(1:end-1,:)';    % end row is boundary points
        faces = t(1:end-1,:)';
        
    end % readNodeMesh


end

