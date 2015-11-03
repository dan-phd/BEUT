classdef LagrangeInterpolator < BEUT.BEM.PiecewisePolynomial
    % The Lagrange interpolator function
    
    properties
        dt;
    end
    
    methods
        % Constructor
        function obj = LagrangeInterpolator(dt, degree)

            % Set the partitions for which each polynomial will act upon
            partition = -dt : dt : degree*dt;
            
            % Get the coefficients for each polynomial
            coefficients = zeros(degree+1);
            for i = 0 : degree
                
                f = 1;
                for phi = 1 : i
                     f = conv(f, [-1/phi/dt 1]);
                end
                
                g = 1;
                for phi = 1 : degree - i
                    g = conv(g, [1/phi/dt 1]);
                end
                
                coefficients(i+1,:) = conv(f, g);
            
            end
            
            % Make this function a piecewise polynomial
            obj = obj@BEUT.BEM.PiecewisePolynomial(partition, coefficients, degree);
            obj.dt = dt;
        end
        
        % pad the coefficients of any number of Lagrange interpolators so that they match obj1
        function varargout = padCoeffs(obj1, varargin)
            
            varargout = varargin;
            for i=1:length(varargin)
                varargout{i} = varargin{i};
                
                % pad columns
                varargout{i}.coeffs = padarray(varargout{i}.coeffs, ...
                    [0, size(obj1.coeffs,2)-size(varargin{i}.coeffs,2)],'pre');
                
                % pad rows
                varargout{i}.coeffs = padarray(varargout{i}.coeffs, ...
                    [size(obj1.coeffs,1)-size(varargin{i}.coeffs,1), 0],'post');
                
                % pad partition
                varargout{i}.partition = padarray(varargout{i}.partition, ...
                    [0, size(obj1.coeffs,1)-size(varargin{i}.coeffs,1)],'post');
            end
            
        end
    end
    
end

