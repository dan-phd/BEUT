% Demonstrate Lagrange interpolator class
clear all

dt = 0.1;
x=(-2:0.01:6)*dt;
degree = 2;

L = BEUT.BEM.LagrangeInterpolator(dt,degree);


%% Plot derivatives and integrals
dL = L.diff;                    % differentiated
d2L = diff(diff(L));            % double differentiated
d3L = L.diff.diff.diff;         % triple differentiated
intL = (int(L));                % integrated
int2L = int(intL);              % double integrated
subplot(4,2,[1 2]); h(1)=plot(x,L.eval(x)); title(sprintf('Lagrange interpolator function of degree %i', degree));
subplot(4,2,3); h(3)=plot(x,dL.eval(x)); title('Differentiated');
subplot(4,2,4); h(2)=plot(x,intL.eval(x)); title('Integrated');
subplot(4,2,5); h(5)=plot(x,d2L.eval(x)); title('Double differentiated');
subplot(4,2,6); h(4)=plot(x,int2L.eval(x)); title('Double integrated');
subplot(4,2,7); h(6)=plot(x,d3L.eval(x)); title('Triple differentiated');


%% Plot translated function
shiftFactor = 3*dt;
scaleFactor = -1;
shiftedL = L.translate(shiftFactor,scaleFactor); % shift and scale
subplot(4,2,8); h(7)=plot(x,shiftedL.eval(x)); title(sprintf('Shifted by %g, scaled by %g', shiftFactor,scaleFactor));

set(h(:),'LineWidth',2)