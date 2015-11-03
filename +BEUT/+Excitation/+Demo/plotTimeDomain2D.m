function plotTimeDomain2D(wave, polarization, rho, time, domain)

N_T = length(time);

% Use quiver to plot vectors when polarization is a vector, otherwise use surf
if numel(polarization)>1
    
    % Field components appear separately
    field = zeros(numel(domain),1*2,N_T);
    for j=1:N_T
        field(:,2*1-1:2*1,j) = wave.eval(time(j),rho) * polarization;
    end
    
    u = squeeze(field(:,1,:));
    v = squeeze(field(:,2,:));
    % Matlab scales quiver arrows weirdly when using very small dt, so not using time here
    figure; quiver(0:N_T-1,domain,u,v)
    xlabel('timestep');
    
else
    
    for j=1:N_T
        V_source(:,j) = wave.eval(time(j),rho)*polarization;
    end
    
    figure; surf(time,domain,V_source,'LineStyle','None')
    xlabel('time');
    
end
ylabel('y');

end