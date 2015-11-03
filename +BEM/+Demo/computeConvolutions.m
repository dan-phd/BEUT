% Demonstrate ComputeConvolutions function
c = 1;
dt = 0.1/c;

P=(1e-6:1e-6:1)';
k=[0 1 2 6 9];
degree = 2;

figure;
tF = tic;

% Pad the interpolators so the sizes match (Ns and D to match Nh)
timeBasis = BEUT.BEM.LagrangeInterpolator(dt,degree);
timeBasis_D = timeBasis;
timeBasis_S = diff(timeBasis);
timeBasis_Nh = int(timeBasis);
timeBasis_Ns = diff(timeBasis);
[timeBasis_Ns, timeBasis_D] = padCoeffs(timeBasis_Nh, timeBasis_Ns, timeBasis_D);

Fh = zeros(length(P),length(k));
Fs = zeros(length(P),length(k));
dF = zeros(length(P),length(k));
for K=1:length(k)
    
    % Shifted time basis - shift and flip the time basis functions to get T(k*dt-t)
    shiftedTB_D  = timeBasis_D. translate(k(K)*dt,-1);
    shiftedTB_S  = timeBasis_S. translate(k(K)*dt,-1);
    shiftedTB_Nh = timeBasis_Nh.translate(k(K)*dt,-1);
    shiftedTB_Ns = timeBasis_Ns.translate(k(K)*dt,-1);
    
    % Compute the temporal convolutions
    [Fh(:,K), Fs(:,K), dF(:,K)] = BEUT.BEM.computeConvolutions(P/c, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D);
    dF=dF/c;
    
    h1=subplot(2,2,1); plot(P,Fh(:,K), 'LineWidth',2); title('Fh'); hold on;
    h2=subplot(2,2,2); plot(P,Fs(:,K), 'LineWidth',2); title('F'); hold on;
    h3=subplot(2,2,3); plot(P,dF(:,K), 'LineWidth',2); title('dF');  hold on;
    
end
F_time = toc(tF)

% Set plot properties
axis(h1,[0 P(end) -dt/2 dt])
axis(h2,[0 P(end) -1/dt/2 1/dt])
axis(h3,[0 P(end) -5 3])
xlabel(h1,'P / c'); ylabel(h1,'g \ast \intT dt')
xlabel(h2,'P / c'); ylabel(h2,'g \ast \deltaT/\deltat')
xlabel(h3,'\rho_m - \rho`'); ylabel(h3,'\delta_P(g \ast T)');

% Insert plot legend
for K=1:length(k)
    names = legend;
    legend1 = legend([names.String {sprintf('k=%i',k(K))}]);
end
set(legend1,'Position',[0.6 0.2 0.1 0.25]);



%% Compare with C++ (when running ./bin/2DTDBEM --test computeConvolutions)
CFile = matfile([BEUT.CFolder '\results\computeConvolutions']);
TempConvs = CFile.TempConvs;

BEUT.relError(dF,TempConvs.dF,'NaN',true);
BEUT.relError(Fh,TempConvs.Fh,'NaN',true);
BEUT.relError(Fs,TempConvs.Fs,'NaN',true);
