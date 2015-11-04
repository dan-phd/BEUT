function [ fnum ] = setInnerCircleAsDifferentMaterial( v,f,inner_radius,inner_circle_centre )
%Extract inner circle of a mesh and set its triangles to a different face number.
% Assumes all face numbers begin as one material

if nargin<4
    inner_circle_centre = [0 0];
end

TR=triangulation(f,v);
circumcentres = circumcenter(TR);

for i=1:size(circumcentres,1)
    dist_from_center(i)=norm(circumcentres(i,:)-[inner_circle_centre(1),inner_circle_centre(2)]);
end
centerFaces = find(dist_from_center<inner_radius);

fnum = ones(1,length(f));
for i=1:numel(centerFaces)
    fnum(centerFaces(i)) = 2;
end

end

