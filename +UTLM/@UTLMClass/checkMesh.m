function checkMesh(obj)
% Check mesh for short link lengths

l=vertcat(obj.halfedges.linkLength);

fprintf('\nShortest link length = \t%f\n', min(l) )
fprintf('Median link length = \t%f\n', median(l) )
fprintf('Ratio = \t\t\t\t%f\n\n', median(l)/min(l) )

end
