classdef RHS
    % Compute RHS (V vector) for the TDBEM which includes the incoming wave
    
    properties
        N_T;
        dt;
        geometry;
        test_function;
        display_plot=false;
        excitation;
        polarization;
        Gaussian_points = 3;    % quadrature points
    end
    
    methods
        % Constructor
        function obj = RHS( N_T, dt )
            obj.N_T = N_T;
            obj.dt = dt;
        end
        
        
        function V = compute(obj, tangent)
            
            N_E = numel(obj.geometry);
            
            if nargin < 2, tangent=false; end;
            
            [s,w] = BEUT.BEM.lgquad(obj.Gaussian_points,[0 1]);
            
            pol = obj.polarization;
            if ~tangent, pol = 1; end;
            
            % Concatenate the test function coefficients
            coeffs = vertcat(obj.test_function.pol{:});
            max_deg = max(vertcat(coeffs.degree));
            
            % make ready to multiply with spatial basis functions
            G_coeffs = BEUT.BEM.createGaussianCoeffs(s,w,1,1,max_deg,0);
            
            % Get index table to determine which edges interact
            idx = obj.test_function.idx_table;
            N_F = size(idx,1);
            assert(N_E==max(max(idx)),'The geometry does not support the chosen test functions.');
            
            % Containers
            V = zeros(N_F,obj.N_T);                % V is an array of size N_F for all timesteps
            V_temp = num2cell(zeros(N_F,1));
            
            wait=waitbar(0,'V vector calculation...');   % status bar
            
            for j = 1:obj.N_T
                
                time = (j)*obj.dt;
                
                % First get the quadrature values
                for m=1:N_E
                    
                    a = obj.geometry(m).a;
                    b = obj.geometry(m).b;
                    l = obj.geometry(m).l;
                    t = obj.geometry(m).t;
                    
                    % Location of the points along the current segment
                    rho = bsxfun(@plus, a, s*(b-a));
                    
                    % p f( c*t_j - dot(k,rho) )
                    field = obj.excitation(time,rho) * pol;
                    
                    % Find the field tangent to current segment
                    if tangent
                        field = bsxfun(@times, t, field);
                    end
                    
                    % Integrate over Gaussian quadrature points to determine
                    % values for each test function coefficient
                    v_temp=bsxfun(@times, G_coeffs, field);
                    V_temp{m} = squeeze(sum(sum(v_temp))) * l;
                    
                end
                
                % Now apply test function polynomial coefficients from all functions
                for p=1:N_F
                    
                    % Current test segments
                    index_test = idx(p,:);
                    
                    % functions on current test edges
                    TF = obj.test_function.pol(:,index_test);
                    
                    % create array of function coeffs
                    TC = zeros(size(TF,2),numel(TF{1,1}.coeffs));
                    for i = 1:size(TF,2)
                        TC(i,:) = TF{i,i}.coeffs;
                    end
                    
                    % loop through all contributing edges and multiply test coefficients with quadrature
                    % coefficients
                    for alpha = 1:numel(index_test)
                        
                        active_coeffs = TC(alpha,:).';
                        
                        % Integral of dot( S_m(rho), E_i(rho,t_J) ) using index table
                        V(p,j) = V(p,j) + sum( active_coeffs .* V_temp{index_test(alpha)} );
                        
                    end
                    
                end
                
                waitbar(j/obj.N_T)
                
            end
            close(wait)
            
            if obj.display_plot
                obj.plot_wave(V);
            end
            
        end  % compute
        
        
        % TODO: compute incident field in whole domain (same as excitationDemo?)
        function e_i = compute_domain(obj, domain_X,domain_Y, num_time_steps)
            
            % Number of time steps to sample
            if nargin<4
                num_time_steps = obj.N_T;
            end
            
            % Create <NxNxNT> container for field values
            X = repmat(domain_X,[1,1,num_time_steps]);
            Y = repmat(domain_Y,[1,1,num_time_steps]);
            e_i = zeros(size(X,1),size(Y,2),num_time_steps);
            
            wait=waitbar(0,'Calculating incident field...');   % status bar
            for k = 1:num_time_steps
                for x = 1:size(X,1)
                    
                    % Observed points inside grid
                    r = [X(x,:,k);Y(x,:,k)]';
                    t = k*obj.dt;
                    
                    e_i(x,:,k) = obj.excitation(r,t);
                    
                end
                waitbar(k/num_time_steps)
            end
            close(wait);
            
        end
        
        
        % Plot incident waves in the time domain shown at various boundary edges
        function [] = plot_wave(obj, V, t_lim)
            
            if nargin<3, t_lim=obj.N_T; end
            
            N_E = size(V,1);
            
            time = (0:t_lim-1)*obj.dt;
            
            figure;
            position = 1:4:N_E;
            plot(time,V(position,1:t_lim))
                
            % legend
            for i=1:length(position)
                entries(i) = {sprintf('%i',position(i))};
            end
            legend('String',entries);

        end
        
        end % public methods
        
end

