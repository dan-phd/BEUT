classdef AnalyticalDielectricCylinder
    %Analytical solutions for current density and field in TE and TM modes
    %for a 2D dielectric cylinder
    % ref: Theory and Computation of Electromagnetic Fields - Jian-Ming Jin p. 227
    
    properties
        N_V; R; omega;
        mu0 = 4e-7*pi;               % permeability of free space
        eps0 = 8.854187817e-12;      % permittivity of free space
        
        % Internal material parameters relative to background medium
        mu_r = 1; eps_r = 1;
    end
    
    methods
        function obj = AnalyticalDielectricCylinder(N_V, radius, omega, mu, eps)
            obj.N_V = N_V;
            obj.R = radius;
            obj.omega = omega;
            
            if nargin>3
                obj.mu0=mu;
                obj.eps0=eps;
            end
        end
        
        
        % Output J at each vertex for each frequency (TM and TE case)
        function [TM_J, TE_J] = calcSurfaceCurrents(obj)
            
            k = obj.omega*sqrt(obj.mu0*obj.eps0);         % wavenumber
            k_int = obj.omega*sqrt(obj.mu_r*obj.mu0*obj.eps_r*obj.eps0);   % wavenumber inside cylinder
            numMultipoles = ceil(nthroot(abs(k*obj.R),3)*6 +abs(k*obj.R));  % Summation limits
            
            TE_J_int = zeros(obj.N_V, numel(obj.omega));
            TE_J_ext = zeros(obj.N_V, numel(obj.omega));
            TM_J_int = zeros(obj.N_V, numel(obj.omega));
            TM_J_ext = zeros(obj.N_V, numel(obj.omega));
            inc = zeros(obj.N_V, numel(obj.omega));
            
            % Points to compute on the circle
            theta = 2 * pi / obj.N_V;
            phi = (0 : theta : 2*pi - theta);
            
            for i = 1 : numel(obj.omega)
                
                orders = -numMultipoles(i)-1 : numMultipoles(i)+1;    % orders to compute Hankel over
                
                %  Bessel and Hankel functions outside the dielectric
                bessels   = besselj(orders,  k(i)*obj.R);
                hankels   = besselh(orders,2,k(i)*obj.R);
                dbessels  = 0.5 * (bessels(1:end-2) - bessels(3:end));
                dhankels  = 0.5 * (hankels(1:end-2) - hankels(3:end));
                bessels   = bessels(2:end-1);
                hankels   = hankels(2:end-1);
                
                % Bessel and Hankel functions inside the dielectric
                bessels_int = besselj(orders,k_int(i)*obj.R);
                dbessels_int = 0.5 * (bessels_int(1:end-2) - bessels_int(3:end));
                bessels_int = bessels_int(2:end-1);
                n = orders(2:end-1);
                
                % Coefficents
                a_fraction_top =     sqrt(obj.mu_r)* dbessels.*bessels_int - sqrt(obj.eps_r)*bessels.*dbessels_int;
                b_fraction_top =     sqrt(obj.eps_r)*dbessels.*bessels_int - sqrt(obj.mu_r)* bessels.*dbessels_int;
                TM_fraction_bottom = sqrt(obj.mu_r)* dhankels.*bessels_int - sqrt(obj.eps_r)*hankels.*dbessels_int;
                TE_fraction_bottom = sqrt(obj.eps_r)*dhankels.*bessels_int - sqrt(obj.mu_r)* hankels.*dbessels_int;
                
                a_n = -1i.^(-n) .* a_fraction_top./TM_fraction_bottom;
                b_n = -1i.^(-n) .* b_fraction_top./TE_fraction_bottom;
                c_n = ( 1i.^(-n+1) / (pi*k(i)*obj.R) ) .* 2*sqrt(obj.mu_r) ./TM_fraction_bottom;
                d_n = ( 1i.^(-n+1) / (pi*k(i)*obj.R) ) .* 2*sqrt(obj.eps_r)./TE_fraction_bottom;
                
                [Phi,N] = ndgrid(phi,n);
                exponentials = exp(1i * N .* Phi);
                
                % Sum each row to sum over all orders
                inc(:,i) =    sum(bsxfun(@times,bsxfun(@times, exponentials, bessels), 1i.^(-n) ),2);
                TM_J_ext(:,i) =  sum(bsxfun(@times,bsxfun(@times, exponentials, a_n),hankels  ),2);
                TE_J_ext(:,i) =  sum(bsxfun(@times,bsxfun(@times, exponentials, b_n),hankels  ),2);
                TM_J_int(:,i) = sum(bsxfun(@times,bsxfun(@times, exponentials, c_n),bessels_int),2);
                TE_J_int(:,i) = sum(bsxfun(@times,bsxfun(@times, exponentials, d_n),bessels_int),2);
            end
            
            % The field at the surface is continuous
            err1 = BEUT.relError(abs(inc+TM_J_ext),abs(TM_J_int),'display',false);
            err2 = BEUT.relError(abs(inc+TE_J_ext),abs(TE_J_int),'display',false);
            if err1+err2<1e-10
                TM_J = TM_J_int;
                TE_J = TE_J_int;
            else
                error('Continuity is not conserved.')
            end
            
        end
        
        
        % Output E and H fields at whole domain for 1 frequency
        function [E,H] = calcFields(obj, chosen_omega, domain_width, resolution, plot)
            
            if nargin<5
                plot=false;
                if nargin<4
                    if nargin<3
                        domain_width = 5*obj.R;
                        if nargin<2
                            chosen_omega = obj.omega(end);
                        end
                    end
                    resolution = domain_width/500;
                end
            end
            
            % Find where the chosen frequency is
            idx = find(obj.omega==chosen_omega);
            k = obj.omega(idx)*sqrt(obj.mu0*obj.eps0);         % wavenumber
            k_int = obj.omega(idx)*sqrt(obj.mu_r*obj.mu0*obj.eps_r*obj.eps0);   % wavenumber inside cylinder
            numMultipoles = ceil(nthroot(abs(k*obj.R),3)*6 +abs(k*obj.R));  % Summation limits
            orders = -numMultipoles-1 : numMultipoles+1;    % orders to compute Hankel over
            
            %  Bessel and Hankel functions outside the dielectric
            bessels   = besselj(orders,k  *obj.R);
            hankels   = besselh(orders,2,k*obj.R);
            dbessels  = 0.5 * (bessels(1:end-2) - bessels(3:end));
            dhankels  = 0.5 * (hankels(1:end-2) - hankels(3:end));
            bessels   = bessels(2:end-1);
            hankels   = hankels(2:end-1);
            
            % Bessel and Hankel functions inside the dielectric
            bessels_int = besselj(orders,k_int*obj.R);
            dbessels_int = 0.5 * (bessels_int(1:end-2) - bessels_int(3:end));
            bessels_int = bessels_int(2:end-1);
            
            n = orders(2:end-1);
            
            % Coefficents
            a_fraction_top =     sqrt(obj.mu_r)* dbessels.*bessels_int - sqrt(obj.eps_r)*bessels.*dbessels_int;
            b_fraction_top =     sqrt(obj.eps_r)*dbessels.*bessels_int - sqrt(obj.mu_r)* bessels.*dbessels_int;
            TM_fraction_bottom = sqrt(obj.mu_r)* dhankels.*bessels_int - sqrt(obj.eps_r)*hankels.*dbessels_int;
            TE_fraction_bottom = sqrt(obj.eps_r)*dhankels.*bessels_int - sqrt(obj.mu_r)* hankels.*dbessels_int;
            
            % Place coefficients in the 3rd dimension
            a_n(1,1,:) = -1i.^(-n) .* a_fraction_top./TM_fraction_bottom;
            b_n(1,1,:) = -1i.^(-n) .* b_fraction_top./TE_fraction_bottom;
            c_n(1,1,:) = ( 1i.^(-n+1) / (pi*k*obj.R) ) .* 2*sqrt(obj.mu_r) ./TM_fraction_bottom;
            d_n(1,1,:) = ( 1i.^(-n+1) / (pi*k*obj.R) ) .* 2*sqrt(obj.eps_r)./TE_fraction_bottom;
            
            
            % Angle, phi and distance from centre, rho
            theta = pi / 90;
            phi = 0 : theta : 2*pi;
            rho_int = 0:resolution:obj.R-resolution;
            rho_ext = obj.R:resolution:domain_width;
            
            % Calculate values of orders and place them in the 3rd dimension
            [Rho_int,Phi_int,N_int] = ndgrid(rho_int,phi,n);
            [Rho_ext,Phi_ext,N_ext] = ndgrid(rho_ext,phi,n);
            bessels_int = besselj(N_int, k_int * Rho_int);
            bessels = besselj(N_ext,    k * Rho_ext);
            hankels = besselh(N_ext, 2, k * Rho_ext);
            exponentials_int = exp(1i * N_int .* Phi_int);
            exponentials_ext = exp(1i * N_ext .* Phi_ext);
            
            % Sum each 3rd dim to sum over all orders and end up with a 2D array of
            % field values at phi and rho
            % Outside cylinder
            E_ext =  sum(bsxfun(@times, a_n, hankels.*exponentials_ext),3);
            H_ext =  sum(bsxfun(@times, b_n, hankels.*exponentials_ext),3);
            
            % Inside
            E_int = sum(bsxfun(@times, c_n, bessels_int.*exponentials_int),3);
            H_int = sum(bsxfun(@times, d_n, bessels_int.*exponentials_int),3);
            
            % Incident
            im(1,1,:) = 1i.^(-n);
            inc = sum(bsxfun(@times, im, bessels.*exponentials_ext),3);
            E_tot = inc + E_ext;
            H_tot = inc + H_ext;
            
            % Convert E_int(rho,phi)+E_ext(rho,phi) to total field in cartesian coords, E(x,y)
            [E(:,:,1),E(:,:,2),E(:,:,3)] = pol2cart([Phi_int(:,:,1);Phi_ext(:,:,1)],...
                [Rho_int(:,:,1);Rho_ext(:,:,1)],[abs(E_int);abs(E_tot)]);
            [H(:,:,1),H(:,:,2),H(:,:,3)] = pol2cart([Phi_int(:,:,1);Phi_ext(:,:,1)],...
                [Rho_int(:,:,1);Rho_ext(:,:,1)],[abs(H_int);abs(H_tot)]);
            
            figure('color','white');
            surf(E(:,:,1),E(:,:,2),E(:,:,3),'Linestyle','none')
            colorbar; view(2); axis equal; title('Total E field');
            
            figure('color','white');
            surf(H(:,:,1),H(:,:,2),H(:,:,3),'Linestyle','none')
            colorbar; view(2); axis equal; title('Total H field');
            
            
            % Plot internal and external fields
            if plot
                [E_ext_car(:,:,1),E_ext_car(:,:,2),E_ext_car(:,:,3)] = pol2cart(Phi_ext(:,:,1),Rho_ext(:,:,1),abs(E_ext));
                [E_int_car(:,:,1),E_int_car(:,:,2),E_int_car(:,:,3)] = pol2cart(Phi_int(:,:,1),Rho_int(:,:,1),abs(E_int));
                [E_tot_car(:,:,1),E_tot_car(:,:,2),E_tot_car(:,:,3)] = pol2cart(Phi_ext(:,:,1),Rho_ext(:,:,1),abs(E_tot));
                [H_ext_car(:,:,1),H_ext_car(:,:,2),H_ext_car(:,:,3)] = pol2cart(Phi_ext(:,:,1),Rho_ext(:,:,1),abs(H_ext));
                [H_int_car(:,:,1),H_int_car(:,:,2),H_int_car(:,:,3)] = pol2cart(Phi_int(:,:,1),Rho_int(:,:,1),abs(H_int));
                [H_tot_car(:,:,1),H_tot_car(:,:,2),H_tot_car(:,:,3)] = pol2cart(Phi_ext(:,:,1),Rho_ext(:,:,1),abs(H_tot));
                
                figure('color','white');
                surf(E_ext_car(:,:,1),E_ext_car(:,:,2),E_ext_car(:,:,3),'Linestyle','none')
                colorbar; view(2); axis equal; title('Scattered E field outside dielectric');
                figure('color','white');
                surf(E_int_car(:,:,1),E_int_car(:,:,2),E_int_car(:,:,3),'Linestyle','none')
                colorbar; view(2); axis equal; title('Scattered E field inside dielectric');
                figure('color','white');
                surf(E_tot_car(:,:,1),E_tot_car(:,:,2),E_tot_car(:,:,3),'Linestyle','none')
                colorbar; view(2); axis equal; title('E_s + E_i');
                
                figure('color','white');
                surf(H_ext_car(:,:,1),H_ext_car(:,:,2),H_ext_car(:,:,3),'Linestyle','none')
                colorbar; view(2); axis equal; title('Scattered H field outside dielectric');
                figure('color','white');
                surf(H_int_car(:,:,1),H_int_car(:,:,2),H_int_car(:,:,3),'Linestyle','none')
                colorbar; view(2); axis equal; title('Scattered H field inside dielectric');
                figure('color','white');
                surf(H_tot_car(:,:,1),H_tot_car(:,:,2),H_tot_car(:,:,3),'Linestyle','none')
                colorbar; view(2); axis equal; title('H_s + H_i');
            end
            
        end
    end
    
end

