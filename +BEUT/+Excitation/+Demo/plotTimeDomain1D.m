function plotTimeDomain1D(wave, polarization, rho, time)

N_T = length(time);
% Use quiver to plot vectors when polarization is a vector, otherwise use plot
if numel(polarization)>1
    for j=1:N_T
        field(:,j) = wave.eval(time(j),rho) * polarization;
    end
    figure; quiver(time,0,field(1,:),field(2,:));
else
    figure; plot(time,wave.eval(time,rho) * polarization);
end
xlabel('time');

end