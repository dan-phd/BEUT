% Create a single mesh which combines 2 different meshes

path = [fileparts(which('BEUT.Meshing.load')) filesep 'meshes' filesep];

filename1 = [path 'cyl_res21.mat'];
filename2 = [path 'cyl_res21.mat'];
[v1,f1]=BEUT.Meshing.load(filename1);
[v2,f2]=BEUT.Meshing.load(filename2);


%% offset the second geometry
offset = [5 0];

offset_ = repmat(offset,size(v2,1),1);
v2_ = v2 + offset_;
f2_ = f2 + max(max(f1));        % vertices for 2nd shape will begin after all vetices for 1st shape
v = [v1; v2_];
f = [f1; f2_];
TR=triangulation(f,v);

% Plot
figure; triplot(TR); axis equal;


%% Save mesh with fnum (1 for 1st shape, 2 for 2nd)
fnum = 2*ones(1,size(f,1));
for i=1:size(f1,1)
    fnum(i) = 1;
end

mesh = BEUT.UTLM.UTLMClass(v,f,fnum);
mesh.checkMesh;                    % look for a ratio below 10, ideally below 5

[~,name1,~] = fileparts(filename1);
[~,name2,~] = fileparts(filename2);
if strcmp(name1,name2)
    filename = [name1 '_x2'];
else
    filename = [name1 name2];
end

c_file = [BEUT.CFolder filesep 'input' filesep filename];
mat_file = [path filename];

use_dual = true;
BEUT.Meshing.save( mesh, use_dual, mat_file, c_file );
