function plotCurrentDensityAtOneFrequency( N_V, J_f, omega, probe_freq_idx, Analytical_J )
% Plot J across geometry at a select frequency

theta = 0 : 2*pi/N_V : 2*pi-2*pi/N_V;

figure; hold all;
title({sprintf('J with respect to theta at omega = %.3g Hz', omega(probe_freq_idx)/2/pi)})
set(gca,'XTick',0:pi/2:2*pi)
set(gca,'XTickLabel',{'0','pi/2','pi','3/2 pi','2pi'})
xlabel('theta'); ylabel('J(\rho,\omega)');
plot(theta,abs(J_f(:,probe_freq_idx)));                 % numerical results
if nargin>4
    plot(theta,abs(Analytical_J(:,probe_freq_idx)),'.k');   % analytical results
    legend('Numerical','Analytical')
end


end

