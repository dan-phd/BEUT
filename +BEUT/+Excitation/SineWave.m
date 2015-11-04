classdef SineWave < BEUT.Excitation.Excitation
% 2D Sinusoidal wave
    
    properties
        f;
        envelope=@Gaussian; % envelope can be "Gaussian" or "cosine" function
        G;  % Gaussian pulse (for envelope)
    end
    
    methods
        
        % Constructor
        function obj = SineWave(freq_width, modulated_frequency, c, direction, timeOfArrivalRratio)
            
            if nargin<5
                timeOfArrivalRratio=1.5;
                if nargin<4
                    direction = [1 0];
                    if nargin<3
                        c = 1;
                    end
                end
            end
            
            obj.direction = direction;
            obj.c = c;
            obj.A = 1;
            
            obj.f = modulated_frequency;
            obj.T = 4/freq_width;
            obj.t0 = timeOfArrivalRratio;
            
            % Create Gaussian signal
            obj.G = BEUT.Excitation.GaussianWave(obj.T, obj.t0, obj.c, obj.direction);
        end
        
        
        % Evaluate E_i (Gaussian pulse envelope)
        function E = eval(self, t, rho)
            
            if nargin<3, rho=0; end;
            x = self.get_x(rho);
            
            y = self.A * sin( 2*pi*self.f/self.c* (x - t*self.c) );
            
            env = self.envelope(self, t, rho);
            
            E = reshape(env,size(y)) .* y;
            
        end
        
        % Cosine function for window
        function window = cosine(~,t,rho)
            L = max(length(rho),length(t));
            w = 0.5*(1-cos(2*pi*(1:L)/round(L/2)));  % Window derived from cosine wave
            window=ones(1,L);
            quarter_L=round(L/4);
            window(1,1:quarter_L)=w(1,1:quarter_L);         % Rising part of window
            window(1,3*quarter_L:L)=w(1,3*quarter_L:L);     % Falling part of window
        end
        
    end
    
    methods (Access = protected)
        
        function x = get_x(self, rho)
            x = get_x@BEUT.Excitation.Excitation(self, rho);
        end
        
    end
    
    methods (Access = private)
        
        % Gaussian function for window
        function g = Gaussian(self,t,rho)
            g = self.G.eval(t, rho);
        end
        
        
    end
    
end

