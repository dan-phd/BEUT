function [ f_new, fnum_new ] = extractMaterials( f, fnum, wanted_materials, change_materials_to )
%Extract materials and rename the material numbers if required.
% wanted_materials and change_materials_to are arrays of the same size
% e.g. to extract materials 2 and 3, and change them to 1 and 2,
% wanted_materials = [2,3], change_materials_to = [1,2]

numMaterials = max(fnum);
for i=1:numMaterials
    materialFaces{i} = find(fnum==i);
end

assert(all(numMaterials>=wanted_materials),'Mesh must include your desired materials');

wantedFaces = materialFaces{wanted_materials(1)};
for i=2:numel(wanted_materials)
    wantedFaces = [ wantedFaces; materialFaces{wanted_materials(i)} ];
end

fnum_new = fnum(wantedFaces);
f_new = f(wantedFaces,:);

if nargin>3
    for i=1:numel(change_materials_to)
        fnum_new(fnum_new==wanted_materials(i)) = change_materials_to(i);
    end
end

end

