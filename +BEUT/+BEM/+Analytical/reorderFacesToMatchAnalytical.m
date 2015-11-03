function [ halfedges, row ] = reorderFacesToMatchAnalytical( halfedges )
%Re-order boundary edges of the cylinder to match the analytical
%solution(i.e. rotate the cylinder)

a=vertcat(halfedges.a);

% analytical vertices start from far right, so find the maximum
% x-coordinate
[row,~]=find(a(:,1)==max(a(:,1),[],1));

N_V = size(a,1);
halfedges = halfedges([row:N_V 1:row-1],:);


end

