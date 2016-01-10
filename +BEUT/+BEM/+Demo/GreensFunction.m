%% Plot a variety of 2D Greens function in the time domain

c = 299792458;
dt = 1/c;
t=(-2:0.01:10)*dt;
P=[0 1.5 4 6];

% 2D Greens function
g = @(t,P) heaviside(t-P/c)./(2*pi*sqrt(t.^2-(P/c).^2));

figure;
for p = 1:length(P)
    plot(t/dt,g(t,P(p)),'LineWidth',2);
    hold all
end
axis([0 max(t)/dt 0 c/3])
set(gca,'yticklabel',{})
xlabel('t / \Deltat'); ylabel('g(P,t)');
title('Greens function (\Deltat=1/c)');

% Plot legend
for i=1:length(P)
    names = legend;
    h = legend([names.String {sprintf('P=%g',P(i))}]);
end


%% Plot a variety of Lagrange interpolators in the time domain
x=(-2:0.01:10)*dt;
k=[0 2 5 6];

L = BEUT.BEM.LagrangeInterpolator(dt,2);

figure;
for K = 1:length(k)
    plot(x/dt,L.translate(k(K)*dt,-1).eval(x), 'LineWidth',2);
    hold all
end
axis([-2 max(x)/dt -0.2 1.2])
xlabel('P / \Deltat'); ylabel('T(k\Deltat - t)');
title('Shifted Lagrange interpolator');

% Plot legend
for i=1:length(P)
    names = legend;
    h = legend([names.String {sprintf('k=%g',k(i))}]);
end
