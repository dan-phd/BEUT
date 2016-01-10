function plotAgainstAnalyticalCylinder(V,time,c,radius)
% Plot comparison of numerical results, V, with analytical resonances
% of a cylindrical cavity

num_nodes = size(V,1);
NT = size(V,2);
num_overlays = 100;     % how many frequency plots to overlay

% Create window to filter time signal
w = 0.5*(1-cos(2*pi*(1:NT)/(NT/2)));
window=ones(1,NT);
window(1,round(3*NT/4):NT)=w(1,round(3*NT/4):NT);

% Plot time and frequency domain
max_height = 0;
figure; xlabel('frequency'); hold on;
for n=1:round(num_nodes/num_overlays):num_nodes
    
    % Apply window for better FFT results
    V(n,:) = V(n,:).*window;
    
    % FFT
    [V_FFT,f,~,~,pks,locs]=BEUT.customFFT(time,V(n,:),'threshold',50,'limit',20);
%     plot(f,abs(V_FFT))
    
    % Display just the dominant frequencies
    stem(f(locs),pks,'bx')
    
    % Find the maximum height of stem plot
    if max(pks)>max_height, max_height=max(pks); end
end


% Calculate theoretical resonances for cylindrical cavity
[~,fr] = BEUT.UTLM.Analytical.findCylinderResonantFrequencies(radius,c,12,8);


% Show where the resonant frequencies should be on frequency graph
for i=1:numel(fr)
    text(fr(i),max_height/2,'\downarrow','HorizontalAlignment','center','FontWeight','bold')
end
hold off

% Limit axis to the frequencies that UTLM is accurate for (low frequencies)
axis([0 fr(12) -inf inf])
title('Comparison with analytical resonant frequencies (arrows)');

end

