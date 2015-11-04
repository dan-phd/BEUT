classdef PiecewisePolynomial
    % The Lagrange interpolator function
    
    properties
        degree;
        coeffs;
        partition;
    end
    
    methods
        % Constructor
        function obj = PiecewisePolynomial(partition, coeffs, degree)
            obj.partition = partition;
            obj.coeffs = coeffs;
            obj.degree = degree;
        end
        
        % Evaluate the function using time axis, t
        function y = eval(self, t)
            
            y = zeros(size(t));
            for interval = 1:length(self.partition)-1
                idxs = (t >= self.partition(interval)) & (t <= self.partition(interval+1));
                y(idxs) = polyval(self.coeffs(interval,:), t(idxs));
            end
        end
        
        % Shift function using k (includes dt) and scale using p, to obtain T(k+pt)
        function obj = translate(self, k, p)
            
            obj = self;     % Copy the current object to create a new one
            eff_degree = size(self.coeffs,2);       % effective number of degrees
            
            % Calculate new coefficients from transform and shift procedure
            new_coeffs = zeros(size(self.coeffs,1),eff_degree);
            for i = 1 : size(self.coeffs,1)
                q = self.coeffs(i,1);
                for j = 2 : eff_degree
                    q = conv([p k], q);
                    q(end) = q(end) + self.coeffs(i,j);
                end
                new_coeffs(i,:) = q;
            end
            
            % Sort out partitioning for new object
            new_partition = (self.partition-k)/p;
            % Mirror the function if needed
            if p < 0
                new_partition = flipdim(new_partition,2);
                new_coeffs = flipdim(new_coeffs, 1);
            end
            
            obj.partition = new_partition;
            obj.coeffs = new_coeffs;

        end
        
        % Differentiate the incoming object
        function obj = diff(self)
            
            obj = self;     % Copy the current object to create a new one
            
            % Resize the coefficient matrix so there is 1 less coefficient in each partition
            obj.coeffs=obj.coeffs(:,1:end-1);
            
            % Derivative for each polynomial row of the coefficient matrix
            for i=1:size(obj.coeffs,1)
                obj.coeffs(i,:) = polyder(self.coeffs(i,:));
                
            end
            
            obj.degree = self.degree - 1;
            
        end
        
        % Integrate the incoming object
        function obj = int(self)
            
            obj = self;     % Copy the current object to create a new one
            
            % Divide each coefficient by its degree+1
            obj.coeffs = bsxfun(@rdivide,obj.coeffs,self.degree+1:-1:1);
            
            % Make sure the end partition goes to infinity
            if obj.partition(end) ~= inf
                obj.partition(end+1) = inf;
                obj.coeffs(self.degree+2,self.degree+2) = 0;
            else
                obj.coeffs(self.degree+1,self.degree+2) = 0;
            end
            
            % Since this function will be smoother, compute the extra
            % (constant) coefficient for each polynomial using initial condition
            % that the very beginning will always be 0
            f_x = 0;
            for i=1:length(obj.partition)-1
                obj.coeffs(i,end) = f_x - polyval(obj.coeffs(i,:),obj.partition(i));
                f_x = polyval(obj.coeffs(i,:), obj.partition(i+1));
            end
            
            obj.degree = self.degree + 1;
            
        end
    end
end

