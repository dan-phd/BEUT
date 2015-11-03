function [ J_FFT,omega,A,limit ] = findFrequencyDomainCurrentDensity( time,J_time,excitation,c )
% Calculate current density in frequency domain at all points on circle

[J_FFT,f_axis,~] = BEUT.customFFT(time,J_time,'dimension',2);
omega = 2*pi*f_axis;

% Scale the response to a plane wave of unit amplitude
A = excitation.evalAmplitudeResponse(omega);
J_FFT = bsxfun(@rdivide,J_FFT,A);

% Frequency domain limit
limit = find(abs(A)<0.001/c, 1 );

end

