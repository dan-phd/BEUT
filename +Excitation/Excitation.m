classdef Excitation
    %Abstract superclass for excitation classes
    
    properties
        c;              % speed of wave propagation
        direction;      % row vector specifying direction as a 2D unit vector
                        % e.g. [1 0] = x-direction, [0 -1] = negative y-direction
                        % 1 = normal to plane (z-direction)
        T;              % width of pulse
        t0;             % time of arrival (scale of T)
        A;              % Amplitude to scale wave by
    end
    
    methods
        
        
        % Plot frequency response
        function Fc = freq_response(obj,time_array,plot_fig)
            
            % Example:
            %
%             c = 10;
%             dt = 0.01/c;
%             time = 0:dt:1;
%             desiredFreqWidth = 10;
%             desiredModulatedFreq = 20;
%             y = BEUT.BEM.Excitation.SineWave(desiredFreqWidth,desiredModulatedFreq,[1 0],c);
%             plot(time,y.eval(0,time'))
%             y.freq_response(time)
%             width = 200*dt;
%             t0 = 2.5*width;
%             G = BEUT.BEM.Excitation.GaussianWave(t0,width,[1 0],c);
%             plot(time,G.eval(0,time'))
%             cutoff_freq = G.freq_response(time)
            
            if isrow(time_array), time_array=time_array'; end
            
            [FFT,f_axis,~,~] = BEUT.customFFT(time_array,obj.eval(time_array));
            
            % output frequency range
            FFT = abs(FFT);
            
            if nargin<3
                plot_fig=true;
            end
            if plot_fig
                figure; plot(f_axis,abs(FFT));
                xlabel('frequency');
            end
            
            below_threshold = FFT(FFT<(max(FFT)/100));
            if isempty(below_threshold)
                error('frequency plot does not reach threshold')
            end
            Fc = f_axis(FFT==max(below_threshold));
            
        end
    end
    
    methods (Abstract)
        
        eval(self, t, rho);
        
    end
    
    methods (Access = protected)
        
        % Get the x vector
        function x = get_x(self, rho)
            
            assert(size(rho,2)<3,'3D co-ordinates are not supported by this function');
            
            % x = dot(rho,self.direction)
            x = sum(bsxfun(@times, self.direction, rho),2);
            
        end
        
    end
    
end

