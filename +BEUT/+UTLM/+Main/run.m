function mesh = run( mesh, N_T, V_source, sourceEdges, ...
    changeBoundaryTimestep, boundaryConditionToChangeTo )
%Run the UTLM algorithm by timestepping
% optional parameters are to change the boundary condition at timestep
% "changeBoundaryTimestep" to "boundaryConditionToChangeTo"

if nargin < 6
    changeBoundaryTimestep = N_T+1;
end

mesh.reset;
tic; bar = waitbar(0,'TLM in time...');
for k=1:N_T
    waitbar(k/N_T,bar);
    
    if k==changeBoundaryTimestep
        mesh.setBoundary(boundaryConditionToChangeTo);
    end
    
    mesh.scatter(k);
    
    mesh.connect(k);
    
    mesh.excite_E(sourceEdges, V_source(k));
    
end
toc; close(bar);

end

