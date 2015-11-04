function [FFT,varargout] = customFFT(t,x,varargin)
% FFT = customFFT(t,x) performs FFT on incoming time domain signal
% [x] and outputs FFT amplitudes array the same size as [x].
%
% [FFT, f] = customFFT(...) outputs frequency array [f].
%
% [FFT, f, zero_point, limit] = customFFT(...) outputs the zero-point
% index of [f] which corresponds to the start of the right quadrant 
% and a limit for the x-axis.
%
% [FFT, f, zero_point, limit, pks, locs] = customFFT(...)
% outputs peak amplitudes [pks] and their locations [locs] only if a 
% threshold to minimum peak height is given.
%
% [...] = customFFT(t,x,'parameter',parameter_value) where the 
% 'parameter' can be:
%   dimension, where one can perform an FFT on multiple time domain 
%   signals [x] along chosen dimension and output FFT amplitudes of
%   size(x).
%   
%   threshold, which can be used to specify the minimum peak height
%   (relative to the maximum peak height) for [pks] and [locs] to be
%   accounted for.
%
%   limit, which can be used to specify a limit for the x-axis [10]
%
% Example:
%     frequencies = [15 35 56];   % Signal frequencies
%     total_time = 10;            % A larger total time gives a more accurate result
%     dt = min(frequencies)/1e4;  % A shorter dt allows a larger frequency range
%     threshold = 8;              % A larger threshold will find more frequency peaks
%     time = (0:dt:total_time-dt);
%     x = sin(2*pi*frequencies(1)*time)+...  % Create time signal
%         0.5*sin(2*pi*frequencies(2)*time)+...
%         0.4*sin(2*pi*frequencies(3)*time);
%     [FFT, f, ~, ~, pks, locs] = BEUT.customFFT(time, x, 'threshold',threshold);
%     figure; hold on; plot(f,abs(FFT));
%     stem(f(locs),pks,'x');
%     for i = 1:length(locs)
%         text(f(locs(i)),pks(i),...
%             sprintf('Frequency = %1.2f \nAmplitude = %1.2f',f(locs(i)),pks(i)),...
%             'HorizontalAlignment','center','VerticalAlignment','top');
%     end
%

%   Author: Daniel Simmons - DansPhD.com
%   Edited: 23/10/2015


% Find values for threshold or dimension
D = 0; threshold = 0; user_limit=5; num_signals=2;
nArgs = length(varargin); a=1;
while a < nArgs+1
    if ischar(varargin{a})
        switch varargin{a}
            
            case 'threshold'
                threshold=varargin{a+1};
                a=a+2;
                
            case 'dimension'
                D=varargin{a+1};
                a=a+2;
                
            case 'limit'
                user_limit=varargin{a+1};
                a=a+2;
                
        end
    end
end

% If D is not given, assume x is an array and set it to 2 and make sure x is a row vector
if D==0
    num_signals=1;
    if isrow(x)==0, x=x'; end
    D=2;
end

N = size(x,D);                      % Length of FFT
dt = t(2)-t(1);

%% Create frequency array
k = ceil(-N/2):floor((N-1)/2);
f = k/(N*dt);                       % scale
varargout{1} = f;
omega = 2 * pi * f;                 % angular frequency


%% Perform FFT
% Matlab "FFT" operation: $$ \hat{f_{k}} = \sum_{j=0}^{N-1} f(j\Delta t)
% e^{-i\frac{2\pi}{N}jk} \qquad k=0,1,...,N-1 $$
% Fourier transform equation: $$  \mathcal{F}(f(t)) = \int\limits_0^T f(t)
% e^{-i\omega t}dt \approx \Delta t \sum_{j=0}^{N-1}f(j\Delta t) e^{-i\omega_{k} j\Delta t}  $$
% Compare frequencies to get scaling factor: $$ \omega_{k}= \frac{2\pi}{N\Delta t}k $$.
% Similarly, comparing amplitudes gives a scaling factor of $$ \Delta t $$.
FFT = fftshift(fft(x,N,D),D);    % shift the negative part of the spectrum to the left quadrant
% FFT = 2*FFT/length(FFT);                       % this scale displays physically correct amplitudes
FFT = dt*bsxfun(@times,FFT,exp(-1i*omega*t(1))); % this scale displays analytically correct amplitudes


%% variables that can be used when plotting frequency spectrum
if nargout>2
    [~,zero_point]=min(abs(f));                 % find where the right quadrant of frequency starts
    % zero_point=round(FFT_length/2+1);         % right quadrant usually starts at FFT_length/2
    
    % Resolve the outputted arrays
    f = f(zero_point:end);
    % FFT can consist of multiple signals
    if num_signals==1
        FFT = FFT(zero_point:end);
    elseif D==2
        FFT = FFT(:,zero_point:end);
    else
        FFT = FFT(zero_point:end,:);
    end
    
    varargout{1} = f;
    varargout{2} = zero_point;
end
if nargout>3
    limit=ceil(length(f)/user_limit);       % specify a limit for the x-axis
    if (limit>length(f))
        limit=length(f);
    end
    
    % Resolve the outputted arrays
    f = f(1:limit);
    % FFT can consist of multiple signals
    if num_signals==1
        FFT = FFT(1:limit);
    elseif D==2
        FFT = FFT(:,1:limit);
    else
        FFT = FFT(1:limit,:);
    end
    
    varargout{1} = f;
    varargout{3} = limit;
end


%% Find peaks for a cleaner graph
if nargout>4
    assert(threshold~=0,'Must provide a threshold if location of peaks are to be output');
    [pks,locs]=findpeaks(abs(FFT),'MINPEAKHEIGHT',max(abs(FFT))/threshold);
    
    varargout{4} = pks;
    varargout{5} = locs;
end

end
