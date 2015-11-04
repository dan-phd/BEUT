%% Resonant frequecies of dielectric cylinder

% Constants
eps0 = 8.854187817e-12;
mu0 = 4e-7*pi;
c = 1/sqrt(eps0*mu0);

% Parameters of cylinder
mu_int = 1 * mu0;
mu_r = mu_int/mu0;
eps_int = 1.5 * eps0;
eps_r = eps_int/eps0;

radius = 1;
numberOfVertices = 120;

% Frequency range
f = (0.001:0.001:2)*c;
omega = 2*pi*f;

analytic = BEUT.BEM.Analytical.AnalyticalDielectricCylinder(numberOfVertices,radius,omega);
analytic.eps_r=eps_r;
analytic.mu_r=mu_r;
[TM_J,TE_J] = analytic.calcSurfaceCurrents;



%% TM plots
% Plot 'J' at shadow side and at exposed side of cylinder with respect to frequency
exposed = numberOfVertices/2; shadow = numberOfVertices;
figure
plot(f/c,abs(TM_J(exposed,:)), f/c,abs(TM_J(shadow,:)))
title('TM J with respect to frequency')
legend('Exposed side','Shadow side')
xlabel('f/c'); ylabel('J(\rho,\omega)');

% Plot J at a particular frequency with respect to theta
figure
freq = 100;          % Pick frequency
theta = 2*pi/numberOfVertices : 2*pi/numberOfVertices : 2*pi;
plot(theta,abs(TM_J(:,freq)),'-r.')
title(sprintf('TM J with respect to theta at %g Hz (radius = %.1flambda)', f(freq), radius*f(freq)/c))
set(gca,'XTick',0:pi/2:2*pi); set(gca,'XTickLabel',{'0','pi/2','pi','3/2 pi','2pi'})
xlabel('theta'); ylabel('J(\rho,\omega)');


%% TE plots
% Plot 'J' at shadow side and at exposed side of cylinder with respect to frequency
exposed = numberOfVertices/2; shadow = numberOfVertices;
figure
plot(f/c,abs(TE_J(exposed,:)), f/c,abs(TE_J(shadow,:)))
title('TE J with respect to frequency')
legend('Exposed side','Shadow side')
xlabel('f/c'); ylabel('J(\rho,\omega)');

% Plot J at a particular frequency with respect to theta
figure
freq = 100;          % Pick frequency
theta = 2*pi/numberOfVertices : 2*pi/numberOfVertices : 2*pi;
plot(theta,abs(TE_J(:,freq)),'-r.')
title(sprintf('TE J with respect to theta at %g Hz (radius = %.1flambda)', f(freq), radius*f(freq)/c))
set(gca,'XTick',0:pi/2:2*pi); set(gca,'XTickLabel',{'0','pi/2','pi','3/2 pi','2pi'})
xlabel('theta'); ylabel('J(\rho,\omega)');


%% Plot fields
domain_width = 10*radius;
resolution = 0.1;
chosen_omega = 2*pi * c;
plot_additional_fields = false;
[E,H] = analytic.calcFields(chosen_omega,domain_width,resolution,false);
