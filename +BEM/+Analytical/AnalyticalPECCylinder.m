classdef AnalyticalPECCylinder
    %Analytical solutions for current density and field in TE and TM modes
    %for a 2D perfectly conducting cylinder
    
    properties
        N_V; radius; omega;
        mu = 4*pi*10^-7;            % permeability of free space
        eps = 8.854187817e-12;      % permittivity of free space
    end
    
    properties (Dependent)
        k; x; numMultipoles; eta;
    end
    
    methods
        function obj = AnalyticalPECCylinder(N_V, radius, omega)
            obj.N_V = N_V;
            obj.radius = radius;
            obj.omega = omega;
        end
        
        %% Getter functions
        function k = get.k(obj)             % Wavenumber
            k = obj.omega*sqrt(obj.mu*obj.eps);
        end
        function x = get.x(obj)             % Hankel/Bessel function argument
            x = obj.k * obj.radius;
        end
        function numMultipoles = get.numMultipoles(obj) % summation limits
            numMultipoles = ceil(nthroot(abs(obj.x),3)*6 +abs(obj.x));
        end
        function eta = get.eta(obj)         % Characteristic impedance of medium
            eta = sqrt(obj.mu/obj.eps);
        end
        
        %% Normal functions
        
        % Output J at each vertex for each frequency (TM case)
        function J = calcTM_J(obj)
            J = zeros(obj.N_V, numel(obj.omega));
            
            for i = 1 : numel(obj.omega)
                
                n = -obj.numMultipoles(i) : obj.numMultipoles(i);   % orders to compute Hankel over
                theta = 2 * pi / obj.N_V;
                phi = (0 : theta : 2*pi - theta);             % compute for all points on circle
                
                hankels = besselh(n,2,obj.x(i));
                
                [Phi,N] = ndgrid(phi,n);
                
                exponentials = exp(1i * N .* Phi);
                
                currentCoeffs = ( -2 / (pi*obj.eta*obj.x(i)) ) * 1i.^(-n) ./ hankels;
                
                % Sum each row to sum over all orders
                J(:,i) =  sum(bsxfun(@times, exponentials, currentCoeffs),2);
            end
            
            J = J * obj.eta;

        end
        
        % Output J at each vertex for each frequency (TE case)
        function J = calcTE_J(obj)
            J = zeros(obj.N_V, numel(obj.omega));
            
            for i = 1 : numel(obj.omega)
                
                orders = -obj.numMultipoles(i)-1 : obj.numMultipoles(i)+1;    % orders to differentiate Hankel over
                theta = 2 * pi / obj.N_V;
                phi = (0 : theta : 2*pi - theta);             % compute for all points on circle
                
                hankels = besselh(orders,2,obj.x(i));
                dhankels = 0.5 * (hankels(1:end-2) - hankels(3:end));
                
                n = orders(2:end-1);
                
                [Phi,N] = ndgrid(phi,n);
                
                exponentials = exp(1i * N .* Phi);
                
                currentCoeffs = ( 2*1i / (pi*obj.x(i)) ) * 1i.^(-n) ./ dhankels;
                
                % Sum each row to sum over all orders
                J(:,i) =  sum(bsxfun(@times, exponentials, currentCoeffs),2);
                
            end
        end
        
        % Output E_z at whole domain for 1 frequency (TM case)
        function E = calcTM_E(obj, domain_width, resolution, chosen_omega)
            
            % Find where the chosen frequency is
            idx = find(obj.omega==chosen_omega);
            
            % Calculate the summation limits
            orders = -obj.numMultipoles(idx) : obj.numMultipoles(idx);
            
            % Find the magnetic coefficients using Hankel and Bessel functions
            % and place them in the 3rd dimension
            bessels = besselj(orders,obj.x(idx));
            hankels = besselh(orders,2,obj.x(idx));
            MagneticCoeffs(1,1,:) = -1i.^(-orders) .* (bessels ./ hankels);
            
            % Angle, phi and distance from centre, rho
            theta = pi / 90;
            phi = 0 : theta : 2*pi;
            rho = 1:resolution:domain_width;
            
            % Calculate values of orders and place them in the 3rd dimension
            [Rho,Phi,N] = ndgrid(rho,phi,orders);
            hankels2 = besselh(N, 2, obj.k(idx) * Rho);
            exponentials = exp(1i * N .* Phi);
            
            % Sum each 3rd dim to sum over all orders and end up with a 2D array of
            % field values at phi and rho
            E_s =  sum(bsxfun(@times, MagneticCoeffs, hankels2.*exponentials),3);
            
            % There are no fields inside scatterer
%             E_s(rho<obj.radius) = 0;
            
            % Convert E_s(rho,phi) to E(x,y)
            [E(:,:,1),E(:,:,2),E(:,:,3)] = pol2cart(Phi(:,:,1),Rho(:,:,1),abs(E_s));
            
            figure('color','white');
            surf(E(:,:,1),E(:,:,2),E(:,:,3),'Linestyle','none')
            colorbar; view(2);
        end
        
        % Output H_z at whole domain for 1 frequency (TE case)
        function H = calcTE_H(obj, domain_width, resolution, chosen_omega)
            
            % Find where the chosen frequency is
            idx = find(obj.omega==chosen_omega);
            
            % Differential of Hankel and Bessel functions
            orders = -obj.numMultipoles(idx)-1 : obj.numMultipoles(idx)+1;
            bessels = besselj(orders,  obj.x(idx));
            hankels = besselh(orders,2,obj.x(idx));
            dbessels = 0.5 * (bessels(1:end-2) - bessels(3:end));
            dhankels = 0.5 * (hankels(1:end-2) - hankels(3:end));
            
            % Magnetic coefficients
            n = orders(2:end-1);
            MagneticCoeffs(1,1,:) = -1i.^(-n) .* (dbessels ./ dhankels);
            
            % Angle, phi and distance from centre, rho
            theta = pi / 90;
            phi=0 : theta : 2*pi;
            rho = 1:resolution:domain_width;
            
            % Calculate values of orders and place them in the 3rd dimension
            [Rho,Phi,N] = ndgrid(rho,phi,n);
            hankels2 = besselh(N, 2, obj.k(idx) * Rho);
            exponentials = exp(1i * N .* Phi);
            
            % Sum each 3rd dim to sum over all orders and end up with a 2D array of
            % field values at phi and rho
            H_s =  sum(bsxfun(@times, MagneticCoeffs, hankels2.*exponentials),3);
            
            % There are no fields inside scatterer
            H_s(rho<obj.radius) = 0;
            
            % Convert H_s(rho,phi) to H(x,y)
            [H(:,:,1),H(:,:,2),H(:,:,3)] = pol2cart(Phi(:,:,1),Rho(:,:,1),abs(H_s));
            
            figure('color','white');
            surf(H(:,:,1),H(:,:,2),H(:,:,3),'Linestyle','none')
            colorbar; view(2);
        end
    end
    
end

