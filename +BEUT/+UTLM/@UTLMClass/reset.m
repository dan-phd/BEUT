function reset(obj)
% Re-initialise all voltages

for i=1:obj.nH
    
    obj.halfedges(i).V_linki = 0;
    obj.halfedges(i).V_linkr = 0;
    obj.halfedges(i).V_stub = 0;
    
end

% reset all fields
obj.fields.E_z = zeros(obj.nH,1);
obj.fields.H_xy = zeros(obj.nH,1);
obj.fields.H_x = zeros(obj.nH,1);
obj.fields.H_y = zeros(obj.nH,1);
obj.I0 = zeros(obj.nF,1);
obj.V0 = zeros(obj.nF,1);

end
