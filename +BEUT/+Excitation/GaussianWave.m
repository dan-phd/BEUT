classdef GaussianWave < BEUT.Excitation.Excitation
% Gaussian pulse wave
    
    properties
    end
    
    methods
        
        % Constructor
        function obj = GaussianWave(width, timeOfArrivalRatio, c, direction)
            
            if nargin<4
                direction = [1 0];
                if nargin<3
                    c = 1;
                    if nargin<2
                        timeOfArrivalRatio=1.5;
                    end
                end
            end
            
            obj.direction = direction;
            obj.c = c;
            obj.A = 1;
            obj.t0 = timeOfArrivalRatio * width;
            obj.T = width;
        end
        
        
        % Evaluate E_i
        function E = eval(self, t, rho)
            
            if nargin<3, rho=0; end;
            x = self.get_x(rho);
            
            gamma = (4/self.T/self.c)*(self.c*(t-self.t0)-x);
            
            E = self.A * (4/((self.T*self.c)*sqrt(pi)))*exp(-gamma.^2);
        end
        
        
        % Evaluate differential of E_i
        function d_E = evalDifferential(self, t, rho)
            
            if nargin<3, rho=0; end;
            x = self.get_x(rho);
            
            gamma = (4/self.T/self.c)*(self.c*(t-self.t0)-x);
            
            d_E = -self.A * (32/((self.T*self.c).^2*sqrt(pi)))*gamma*self.c.*exp(-gamma.^2);
        end
        
        
        % Evaluate integral of E_i
        function E = evalIntegral(self, t, dt, rho)
            
            if nargin<4, rho=0; end;
            x = self.get_x(rho);
            
            time_so_far = (0:dt:t);
            
            V = self.eval(time_so_far, x);
            
            if numel(time_so_far)>1
                E = cumtrapz(time_so_far,V);
                E = E(end);
            else
                E = V;
            end
            
        end
        
        
        % Evaluate amplitude response at a certain frequency
        function A = evalAmplitudeResponse(obj,omega)
            A = exp(-1i*omega*obj.t0 - (obj.T*omega/(8)).^2)...
                /obj.c;
        end
        
        
    end
    
    methods (Access = protected)
        
        function x = get_x(self, rho)
            x = get_x@BEUT.Excitation.Excitation(self, rho);
        end
        
    end
    
end

