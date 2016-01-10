% Demonstrate the use of the GaussianWave class

% Parameters
mu = 4*pi*10^-7;            % permeability of free space
eps = 8.854187817e-12;      % permittivity of free space
c = 1/sqrt(mu*eps);         % propagation speed inside vacuum

% plane wave direction and polarization
direction = [1 0];            % [1 0] = x-direction, [0 -1] = negative y-direction
polarization = [0 -1];               % 1 = normal to plane

% Gaussian wave parameters
pulseWidth = 40 / c;
startTimeRatio = 1.5;        % this will be scaled by pulseWidth


%% The following parameters are automatically set but can be modified if needed
% 2D domain (scaled to width of wave for appropriate plots)
domain_x = linspace(0,(1+5*abs(direction(1)))*pulseWidth*c,100);
domain_y = linspace(0,(1+5*abs(direction(2)))*pulseWidth*c,100);
[X,Y] = meshgrid(domain_x,domain_y);

% Temporal discretization
% The maximum frequency from Gaussian pulse is roughly:
maxFreq = 2.5/pulseWidth;
dt = 1/14/maxFreq;      % time for each timestep
% The maximum time where the pulse will sufficiently subside:
maxTime = 2*(startTimeRatio+1)*pulseWidth;
N_T = round(maxTime / dt);
time = (0:N_T-1)*dt;

% Create Gaussian pulse wave
gpw = BEUT.Excitation.GaussianWave(pulseWidth, startTimeRatio, c, direction);
gpw.A = 4;


%% Plot in frequency domain and output cutoff frequency
Fc = gpw.freq_response(time,true)
title('Frequency domain plot of Gaussian pulse signal')


%% Plot 1D wave at particular value of x and y
x_idx = 1; y_idx = 2;
rho = [X(y_idx,x_idx) Y(y_idx,x_idx)];
BEUT.Excitation.Demo.plotTimeDomain1D(gpw, polarization, rho, time);
title( sprintf('Gaussian wave at x=%g, y=%g',domain_x(x_idx),domain_y(y_idx)) )


%% Plot and test differential and integral
V = zeros(N_T,1); V_diff = zeros(N_T,1); V_int = zeros(N_T,1);
for k=1:N_T
    V(k) = gpw.eval((k-1)*dt);
    V_diff(k) = gpw.evalDifferential((k-1)*dt);
    V_int(k) = gpw.evalIntegral((k-1)*dt,dt); 
end
figure; plot(time,V_diff,time(2:end),diff(V)/dt,'.');
legend('Derivative using class','Derivative using Matlab function')
figure; plot(time,V_int,time,cumtrapz(time,V),'.')
legend('Integral using class','Integral using Matlab function')


%% Animate in 1D
BEUT.Excitation.Demo.animate1D(gpw, time, X(1,:))


%% Plot 2D wave at particular value of x
x_idx = 1;
rho = [X(:,x_idx) Y(:,x_idx)];
BEUT.Excitation.Demo.plotTimeDomain2D(gpw, polarization, rho, time, domain_y);
title(['Gaussian pulse plane wave at x=' num2str(domain_x(x_idx))]); view(2)


%% Animate in 2D
BEUT.Excitation.Demo.animate2D(gpw, time, X, Y, polarization)

