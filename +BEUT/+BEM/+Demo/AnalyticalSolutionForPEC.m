% BEM analytical solutions for a PEC cylinder

% Constants
eps = 8.854187817e-12;
mu = 4e-7*pi;
c = 1/sqrt(eps*mu);

% Variables
N = 80;             % Number of segments
radius = 1;
f = (0.01:0.01:3)*c;
omega = 2*pi*f;             % angular frequency = 2pi* phase velocity
eta = sqrt(mu/eps);         % characteristic impedance of free space


%% Find analytical solutions
analytic = BEUT.BEM.Analytical.AnalyticalPECCylinder(N,radius,omega);
analytic.mu=mu;
analytic.eps=eps;


%% TM J
J = analytic.calcTM_J;

J1 = abs(J(N/2,:));       % exposed side
J2 = abs(J(N,:));         % shadow side

%% Plot 'J' at shadow side and at exposed side of cylinder with respect to frequency
figure
plot(f/c,J1, f/c,J2 )
title('J with respect to frequency')
legend('Exposed side','Shadow side')
xlabel('f/c'); ylabel('J(\rho,\omega)');

%% Plot J at a particular frequency with respect to theta
figure
freq = 100;          % Pick frequency
theta = 2*pi/N : 2*pi/N : 2*pi;
plot(theta,abs(J(:,freq)),'-r.')
title(sprintf('J with respect to theta at %g Hz (radius = %.1flambda)', f(freq), radius*f(freq)/c))
set(gca,'XTick',0:pi/2:2*pi)
set(gca,'XTickLabel',{'0','pi/2','pi','3/2 pi','2pi'})
xlabel('theta'); ylabel('J(\rho,\omega)');


%% TE J
J = analytic.calcTE_J;

J1 = abs(J(N/2,:));       % exposed side
J2 = abs(J(N,:));         % shadow side

%% Plot 'J' at shadow side and at exposed side of cylinder with respect to frequency
figure
plot(f/c,J1, f/c,J2 )
title('J with respect to frequency')
legend('Exposed side','Shadow side')
xlabel('f/c'); ylabel('J(\rho,\omega)');

%% Plot J at a particular frequency with respect to theta
figure
freq = 100;          % Pick frequency
theta = 2*pi/N : 2*pi/N : 2*pi;
plot(theta,abs(J(:,freq)),'-r.')
title(sprintf('J with respect to theta at %g Hz (radius = %.1flambda)', f(freq), radius*f(freq)/c))
set(gca,'XTick',0:pi/2:2*pi)
set(gca,'XTickLabel',{'0','pi/2','pi','3/2 pi','2pi'})
xlabel('theta'); ylabel('J(\rho,\omega)');



%% Scattered E field (TM)
domain_width = 10*radius;
resolution = 0.1;
chosen_omega = 2*pi * c;
actual_E_s = analytic.calcTM_E(domain_width,resolution,chosen_omega);

%% Scattered H field (TE)
actual_H_s = analytic.calcTE_H(domain_width,resolution,chosen_omega);

