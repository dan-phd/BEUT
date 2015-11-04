% Demonstrate the use of the SineWave class

% Parameters
mu = 4*pi*10^-7;            % permeability of free space
eps = 8.854187817e-12;      % permittivity of free space
c = 1/sqrt(mu*eps);         % propagation speed inside vacuum

% plane wave direction and polarization
direction = [0 1];            % [1 0] = x-direction, [0 -1] = negative y-direction
polarization = [1];               % 1 = normal to plane

% Sine wave parameters
startTimeRatio = 0.9;
desiredFreqWidth = 0.5e8;
desiredModulatedFreq = 1.2e8;


%% The following parameters are automatically set but can be modified if needed
% 2D domain (scaled to width of wave for appropriate plots)
domain_x = linspace(0,(1+5*abs(direction(1)))*(2/desiredFreqWidth)*c,100);
domain_y = linspace(0,(1+5*abs(direction(2)))*(2/desiredFreqWidth)*c,100);
[X,Y] = meshgrid(domain_x,domain_y);

% Temporal discretization
% The maximum frequency from sinewave pulse is roughly:
maxFreq = desiredModulatedFreq+desiredFreqWidth;
dt = 1/30/maxFreq;      % time for each timestep
% The maximum time where the pulse will sufficiently subside:
maxTime = (1+startTimeRatio)*8/desiredFreqWidth;
N_T = round(maxTime / dt);
time = (0:N_T-1)*dt;

% Create Sinewave signal
sinewave = BEUT.Excitation.SineWave(desiredFreqWidth, desiredModulatedFreq,...
    c, direction, startTimeRatio);
% sinewave.envelope = @cosine;    % envelope can be "Gaussian" (default) or "cosine" function


%% Plot in frequency domain
sinewave.freq_response(time,true);
title('Frequency domain plot of sinusoidal signal')


%% Plot 1D wave at particular value of x and y
x_idx = 1; y_idx = 2;
rho = [X(y_idx,x_idx) Y(y_idx,x_idx)];
BEUT.Excitation.Demo.plotTimeDomain1D(sinewave,polarization, rho, time);
title( sprintf('Sine wave at x=%g, y=%g',domain_x(x_idx),domain_y(y_idx)) )


%% Animate in 1D
BEUT.Excitation.Demo.animate1D(sinewave, time, X(1,:))


%% Plot 2D wave at particular value of x
x_idx = 1;
rho = [X(:,x_idx) Y(:,x_idx)];
BEUT.Excitation.Demo.plotTimeDomain2D(sinewave, polarization, rho, time, domain_y);
title(['Sinusoidal plane wave at x=' num2str(domain_x(x_idx))]); view(2)


%% Animate
BEUT.Excitation.Demo.animate2D(sinewave, time, X, Y, polarization)
