function animate2D(wave, time, X, Y, polarization)

N_T = length(time);

% Use quiver for vectors when polarization is a vector, otherwise use surf
if numel(polarization)>1
    
    % x, y and z field component appears separately
    field = zeros(size(X,1)*2,size(X,2),N_T);
    
    for j=1:N_T
        for x=1:size(X,1)
            rho = [X(x,:); Y(x,:)]';
            field(2*x-1:2*x,:,j) = (wave.eval(time(j),rho) * polarization)';
        end
    end
    
    value.u = field(1:2:size(X,2)*2-1,:,:);
    value.v = field(2:2:size(X,2)*2,:,:);
    
    BEUT.animate_fields(2,'domain',Y(:,1),X(1,:), 'Plane wave',value);
    
else
    
    % if polarization is in z-axis, no need for vector plots
    field = zeros(size(X,1),size(X,2),N_T);
    
    for j=1:N_T
        for x=1:size(X,1)
            rho = [X(x,:); Y(x,:)]';
            field(x,:,j) = wave.eval(time(j),rho) * polarization;
        end
    end
    
    BEUT.animate_fields(2,'domain',Y(:,1),X(1,:), 'Plane wave',field);
end

end