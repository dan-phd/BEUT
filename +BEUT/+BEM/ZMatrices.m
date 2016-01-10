classdef ZMatrices
    % Compute Z matrix for the TDBEM
    
    properties
        N_T;
        dt;
        c;
        
        timeBasis_D;   % for convolution
        timeBasis_Nh;  % for integral convolution
        timeBasis_Ns;  % for differentiated convolution
        
        basis_function_Z; % Z-directed
        basis_function_S; % transverse to plane
        test_function_Z;
        test_function_S;
        
        outer_points_sp = 50;
        inner_points_sp = 51;
        outer_points = 3;
        inner_points = 4;
    end
    
    properties (Access = private)
        
        a_m; b_m; l_m; t_m; n_m;
        a_n; b_n; l_n; t_n; n_n;
        
        shiftedTB_D;
        shiftedTB_Nh;
        shiftedTB_Ns;
        
        N_E; geom_obj;
        N_F;        % number basis functions
        
    end
    
    methods
        
        % Constructor
        function obj = ZMatrices(N_T,dt,geom_obj, c)
            
            if nargin<4
                c = 1;
            end
            
            obj.c = c;
            obj.N_T = N_T;
            obj.dt = dt;
            
            % Object to decide the current segement/vertex during loops
            obj.geom_obj = geom_obj;
            obj.N_E = numel(geom_obj);
        end
        
        
        % Loop over all time and spatial points in the domain and find the Gaussian quadrature tables
        % associated with each edge interaction. Then loop over all spatial points in the domain and, using
        % the Basis functions, determine how to combine the elements from the Gaussian quadrature coefficient
        % table
        function [S,D,Dp,Nh,Ns] = compute( obj, cheat )
            
            if nargin<2, cheat=false; end
            
            % parameterise each segment for Gaussian quadrature integration
            [so,wo] = BEUT.BEM.lgquad(obj.outer_points,[0 1]);
            [si,wi] = BEUT.BEM.lgquad(obj.inner_points,[0 1]);
            
            % Gaussian quadrature integration for singularity points
            [so_sp,wo_sp] = BEUT.BEM.lgquad(obj.outer_points_sp,[0 1]);
            [si_sp,wi_sp] = BEUT.BEM.lgquad(obj.inner_points_sp,[0 1]);
                
            assert(isempty(intersect(so,si)),...
                ['Self-path integrals contain singularities, use a different number of '...
                'Gaussian quadrature points for inner and outer integrals.'])
            
            % Find divergence of transverse basis and test functions
            basis_function_d = BEUT.BEM.BasisFunction.divergence(obj.basis_function_S, obj.geom_obj);
            test_function_d = BEUT.BEM.BasisFunction.divergence(obj.test_function_S, obj.geom_obj);
            
            % Concatenate the basis and test coefficients and find the maximum degree
            basis_coeffs_Z = vertcat(obj.basis_function_Z.pol{:});
            max_deg_basis_Z = max(vertcat(basis_coeffs_Z.degree));
            test_coeffs_Z = vertcat(obj.test_function_Z.pol{:});
            max_deg_test_Z = max(vertcat(test_coeffs_Z.degree));
            basis_coeffs_S = vertcat(obj.basis_function_S.pol{:});
            max_deg_basis_S = max(vertcat(basis_coeffs_S.degree));
            test_coeffs_S = vertcat(obj.test_function_S.pol{:});
            max_deg_test_S = max(vertcat(test_coeffs_S.degree));
            basis_coeffs_d = vertcat(basis_function_d.pol{:});
            max_deg_basis_d = max(vertcat(basis_coeffs_d.degree));
            test_coeffs_d = vertcat(test_function_d.pol{:});
            max_deg_test_d = max(vertcat(test_coeffs_d.degree));
            
            % Make the gaussian quadrature tables, ready to multiply with spatial basis functions
            Gcoeffs_ZZ = BEUT.BEM.createGaussianCoeffs(so,wo,si,wi,max_deg_test_Z,max_deg_basis_Z);
            Gcoeffs_ZS = BEUT.BEM.createGaussianCoeffs(so,wo,si,wi,max_deg_test_Z,max_deg_basis_S);
            Gcoeffs_SZ = BEUT.BEM.createGaussianCoeffs(so,wo,si,wi,max_deg_test_S,max_deg_basis_Z);
            Gcoeffs_SS = BEUT.BEM.createGaussianCoeffs(so,wo,si,wi,max_deg_test_S,max_deg_basis_S);
            Gcoeffs_dd = BEUT.BEM.createGaussianCoeffs(so,wo,si,wi,max_deg_test_d,max_deg_basis_d);
            Gcoeffs_ZZ_sp = BEUT.BEM.createGaussianCoeffs(so_sp,wo_sp,si_sp,wi_sp,max_deg_test_Z,max_deg_basis_Z);
            Gcoeffs_ZS_sp = BEUT.BEM.createGaussianCoeffs(so_sp,wo_sp,si_sp,wi_sp,max_deg_test_Z,max_deg_basis_S);
            Gcoeffs_SZ_sp = BEUT.BEM.createGaussianCoeffs(so_sp,wo_sp,si_sp,wi_sp,max_deg_test_S,max_deg_basis_Z);
            Gcoeffs_SS_sp = BEUT.BEM.createGaussianCoeffs(so_sp,wo_sp,si_sp,wi_sp,max_deg_test_S,max_deg_basis_S);
            Gcoeffs_dd_sp = BEUT.BEM.createGaussianCoeffs(so_sp,wo_sp,si_sp,wi_sp,max_deg_test_d,max_deg_basis_d);
            
            % Get basis function index table to determine which edges interact
            idx_basis_Z = obj.basis_function_Z.idx_table;
            idx_basis_S = obj.basis_function_S.idx_table;
            idx_basis_d = basis_function_d.idx_table;
            idx_test_Z = obj.test_function_Z.idx_table;
            idx_test_S = obj.test_function_S.idx_table;
            idx_test_d = test_function_d.idx_table;
            assert(all([size(idx_basis_S,1) size(idx_basis_d,1) size(idx_test_Z,1) size(idx_test_S,1) size(idx_test_d,1)]...
                == size(idx_basis_Z,1)),'basis functions must act over the same geometry');
            obj.N_F = size(idx_basis_Z,1);
            
            % Containers for operators (2D regions with 3rd dimension varying with time)
            S  = zeros(obj.N_F, obj.N_F, obj.N_T);
            D  = zeros(obj.N_F, obj.N_F, obj.N_T);
            Dp = zeros(obj.N_F, obj.N_F, obj.N_T);
            Ns = zeros(obj.N_F, obj.N_F, obj.N_T);
            Nh = zeros(obj.N_F, obj.N_F, obj.N_T);
            
            % Pad the interpolators so the sizes match Nh (hypersingular contribution of N operator)
            [obj.timeBasis_Ns, obj.timeBasis_D] = padCoeffs(obj.timeBasis_Nh, obj.timeBasis_Ns, obj.timeBasis_D);
            
            % Loop over all segments so they all act as observation and source (m and n) for all time steps
            wait=waitbar(0,'Operator calculations...');   % status bar
            max_n = obj.N_E; max_q = obj.N_F;
            for k = 1:obj.N_T
                
                % Shifted time basis - shift and flip the time basis functions to get T(k*dt-t)
                obj.shiftedTB_D  = obj.timeBasis_D. translate((k-1)*obj.dt,-1);
                obj.shiftedTB_Nh = obj.timeBasis_Nh.translate((k-1)*obj.dt,-1);
                obj.shiftedTB_Ns = obj.timeBasis_Ns.translate((k-1)*obj.dt,-1);
                
                % Containers for operator coefficients
                coeffs_nh = num2cell(zeros(obj.N_E, obj.N_E));
                coeffs_ns = num2cell(zeros(obj.N_E, obj.N_E));
                coeffs_d = num2cell(zeros(obj.N_E, obj.N_E));
                coeffs_s = num2cell(zeros(obj.N_E, obj.N_E));
                coeffs_dp = num2cell(zeros(obj.N_E, obj.N_E));
                
                for m = 1:obj.N_E
                    
                    % Find geometry around observation point
                    obj.a_m = obj.geom_obj(m).a;
                    obj.b_m = obj.geom_obj(m).b;
                    obj.l_m = obj.geom_obj(m).l;
                    obj.t_m = obj.geom_obj(m).t;
                    obj.n_m = obj.geom_obj(m).n;
                    
                    % Determine what source element to go up to
                    if cheat==1
                        if any(m==1:obj.N_E/obj.N_F)
                            max_n=obj.N_E;
                        else max_n=obj.N_E/obj.N_F;
                        end
                    end
                    
                    for n = 1:max_n
                            
                            % Find geometry around source point
                            obj.a_n = obj.geom_obj(n).a;
                            obj.b_n = obj.geom_obj(n).b;
                            obj.l_n = obj.geom_obj(n).l;
                            obj.t_n = obj.geom_obj(n).t;
                            obj.n_n = obj.geom_obj(n).n;
                            
                            % When dealing with singularities at self patch and neighbouring
                            % edges, increase number of quadrature points
                            % TODO: this will change for dual basis or functions spanning more than 2 edges?
                            if n==mod(m,obj.N_E)+1 || n==m || n==mod(m-2,obj.N_E)+1
                                
                                [coeffs_nh{m,n},coeffs_ns{m,n},coeffs_d{m,n},coeffs_s{m,n},coeffs_dp{m,n}] = ...
                                    Z_calc( obj,...
                                    si_sp, so_sp,...
                                    Gcoeffs_dd_sp,...
                                    Gcoeffs_SZ_sp,...
                                    Gcoeffs_ZS_sp,...
                                    Gcoeffs_SS_sp,...
                                    Gcoeffs_ZZ_sp);
                                
                            else
                                
                                % Normal calculation when there are no singularities
                                [coeffs_nh{m,n},coeffs_ns{m,n},coeffs_d{m,n},coeffs_s{m,n},coeffs_dp{m,n}] =...
                                    Z_calc( obj,...
                                    si, so,...
                                    Gcoeffs_dd,...
                                    Gcoeffs_SZ,...
                                    Gcoeffs_ZS,...
                                    Gcoeffs_SS,...
                                    Gcoeffs_ZZ);
                            end
                            
                            
                    end % n
                    
                end % m
                
                if cheat==1
                    number_to_skip = obj.N_E/obj.N_F;
                    coeffs_s  = SPD_cheat(coeffs_s,number_to_skip);
                    coeffs_nh = SPD_cheat(coeffs_nh,number_to_skip);
                    coeffs_ns = SPD_cheat(coeffs_ns,number_to_skip);
                    coeffs_d  = SPD_cheat(coeffs_d,number_to_skip);
                    coeffs_dp = SPD_cheat(coeffs_dp,number_to_skip);
                end % cheat
                
                
                % Use integrated convolution values along with basis/test function coefficients
                % and index table to compute final matrix entries
                for p=1:obj.N_F
                    
                    % Current test edges
                    index_test_Z = idx_test_Z(p,:);
                    index_test_S = idx_test_S(p,:);
                    index_test_d = idx_test_d(p,:);
                    
                    % Functions on current test edges
                    TFZ = obj.test_function_Z.pol(:,index_test_Z);
                    TFS = obj.test_function_S.pol(:,index_test_S);
                    TFd = test_function_d.pol(:,index_test_d);
                    
                    % Create array of function coeffs
                    TCZ = create_function_coeffs_tbl(TFZ);
                    TCS = create_function_coeffs_tbl(TFS);
                    TCd = create_function_coeffs_tbl(TFd);
                    
                    if cheat==1, if p==1, max_q=obj.N_F; else max_q=1; end; end;
                    for q = 1:max_q
                        
                        % current basis edges
                        index_basis_Z = idx_basis_Z(q,:);
                        index_basis_S = idx_basis_S(q,:);
                        index_basis_d = idx_basis_d(q,:);
                        
                        % functions on current basis edges
                        BFZ = obj.basis_function_Z.pol(:,index_basis_Z);
                        BFS = obj.basis_function_S.pol(:,index_basis_S);
                        BFd = basis_function_d.pol(:,index_basis_d);

                        % create array of function coeffs
                        BCZ = create_function_coeffs_tbl(BFZ);
                        BCS = create_function_coeffs_tbl(BFS);
                        BCd = create_function_coeffs_tbl(BFd);
                        
                        % Combine contributions
                        Ns(p,q,k) = combine_contributions(coeffs_ns,index_test_S,index_basis_S,TCS,BCS);
                        Nh(p,q,k) = combine_contributions(coeffs_nh,index_test_d,index_basis_d,TCd,BCd);
                        S(p,q,k)  = combine_contributions(coeffs_s, index_test_Z,index_basis_Z,TCZ,BCZ);
                        Dp(p,q,k) = combine_contributions(coeffs_dp,index_test_S,index_basis_Z,TCS,BCZ);
                        D(p,q,k)  = combine_contributions(coeffs_d, index_test_Z,index_basis_S,TCZ,BCS);
                        
                    end
                end
                
                if cheat==1
                    S(:,:,k)  = SPD_cheat(S(:,:,k));
                    Nh(:,:,k) = SPD_cheat(Nh(:,:,k));
                    Ns(:,:,k) = SPD_cheat(Ns(:,:,k));
                    D(:,:,k)  = SPD_cheat(D(:,:,k));
                    Dp(:,:,k) = SPD_cheat(Dp(:,:,k));
                end % cheat
                
                waitbar(k/obj.N_T)
            end % k
            close(wait)
            
            % create array of function coeffs
            function coeffs = create_function_coeffs_tbl(functions)
                
                coeffs = zeros(size(functions,2),numel(functions{1,1}.coeffs));
                for i = 1:size(functions,2)
                    coeffs(i,:) = functions{i,i}.coeffs;
                end
                    
            end
            
            function Z = combine_contributions(operator_coeffs,test_index,basis_index,test_coeffs,basis_coeffs)
                
                % Loop through all interacting edges and multiply the basis and test polynomial coefficients with
                % the integrated convolution (operator) coefficients
                Z = 0;
                for alpha=1:numel(test_index)
                    for beta=1:numel(basis_index)
                        
                        operator_coeff = operator_coeffs{test_index(alpha),basis_index(beta)};
                        polynomial_coeffs = test_coeffs(alpha,:).' * basis_coeffs(beta,:);
                        Z = Z + sum(sum( polynomial_coeffs .* operator_coeff ));
                        
                    end
                end
                
            end
            
            % Cheat to act on operator coeffs
            function [Z] = SPD_cheat(Z,num_to_skip)
                
                if nargin < 2
                    num_to_skip=1;
                end
                
                % ONLY APPLICABLE FOR PEC CYLINDER (since Z matrix is SPD)
                
                % Using just the first n rows and columns, copy elements diagonally down and right
                for mm = num_to_skip+1:size(Z,1)
                    for nn = num_to_skip+1:size(Z,2)
                        Z(mm,nn) = Z(mm-num_to_skip,nn-num_to_skip);
                    end
                end
                
            end
            
        end  % compute
        
        
        % Find scattered field in a grid of points, rho
        function [S,D] = computeField( obj, rho )
            
            N_points = size(rho,1);
            
            % parameterise each segment for Gaussian quadrature integration
            [si,wi] = BEUT.BEM.lgquad(obj.inner_points,[0 1]);
            
            % Concatenate the basis coefficients and find the maximum degree
            basis_coeffs_Z = vertcat(obj.basis_function_Z.pol{:});
            max_deg_basis_Z = max(vertcat(basis_coeffs_Z.degree));
            basis_coeffs_S = vertcat(obj.basis_function_S.pol{:});
            max_deg_basis_S = max(vertcat(basis_coeffs_S.degree));
            
            % Make the gaussian quadrature tables, ready to multiply with spatial basis functions
            Gcoeffs_ZZ = BEUT.BEM.createGaussianCoeffs(0.5,1,si,wi,0,max_deg_basis_Z);
            Gcoeffs_ZS = BEUT.BEM.createGaussianCoeffs(0.5,1,si,wi,0,max_deg_basis_S);
            
            % Get basis function index table to determine which edges interact
            idx_basis_Z = obj.basis_function_Z.idx_table;
            idx_basis_S = obj.basis_function_S.idx_table;
            assert(size(idx_basis_S,1)== size(idx_basis_Z,1),...
                'basis functions must act over the same geometry');
            obj.N_F = size(idx_basis_Z,1);
            
            % Containers for operators (2D regions with 3rd dimension varying with time)
            S  = zeros(N_points, obj.N_F, obj.N_T);
            D  = zeros(N_points, obj.N_F, obj.N_T);
            
            % Pad the interpolators so the sizes match
            [obj.timeBasis_Ns, obj.timeBasis_D] = padCoeffs(obj.timeBasis_Ns, obj.timeBasis_Ns, obj.timeBasis_D);
            
            % Loop over all segments so they all act as source for all time steps
            wait=waitbar(0,'Operator calculations...');   % status bar
            for k = 1:obj.N_T
                
                % Shifted time basis - shift and flip the time basis functions to get T(k*dt-t)
                obj.shiftedTB_D  = obj.timeBasis_D.translate((k-1)*obj.dt,-1);
                obj.shiftedTB_Ns = obj.timeBasis_Ns.translate((k-1)*obj.dt,-1);
                
                % Containers for operator coefficients
                coeffs_d = num2cell(zeros(N_points, obj.N_E));
                coeffs_s = num2cell(zeros(N_points, obj.N_E));
                
                for m = 1:N_points
                    
                    % Find geometry around observation point
                    rho_m = rho(m,:);
                    
                    for n = 1:obj.N_E
                        
                        % Find geometry around source point
                        obj.a_n = obj.geom_obj(n).a;
                        obj.b_n = obj.geom_obj(n).b;
                        obj.l_n = obj.geom_obj(n).l;
                        obj.n_n = obj.geom_obj(n).n;
                        
                        [coeffs_d{m,n},coeffs_s{m,n}] =...
                            Z_calc_field( obj,...
                            rho_m, si,...
                            Gcoeffs_ZS,...
                            Gcoeffs_ZZ);
                        
                    end % n
                    
                end % m
                
                
                % Use integrated convolution values along with basis/test function coefficients
                % and index table to compute final matrix entries
                for p=1:N_points
                    
                    for q = 1:obj.N_F
                        
                        % current basis edges
                        index_basis_Z = idx_basis_Z(q,:);
                        index_basis_S = idx_basis_S(q,:);
                        
                        % functions on current basis edges
                        BFZ = obj.basis_function_Z.pol(:,index_basis_Z);
                        BFS = obj.basis_function_S.pol(:,index_basis_S);
                        
                        % create array of function coeffs
                        BCZ = create_function_coeffs_tbl(BFZ);
                        BCS = create_function_coeffs_tbl(BFS);
                        
                        % Combine contributions
                        S(p,q,k)  = combine_contributions(coeffs_s, p,index_basis_Z,1,BCZ);
                        D(p,q,k)  = combine_contributions(coeffs_d, p,index_basis_S,1,BCS);
                        
                    end
                end
                
                waitbar(k/obj.N_T)
            end % k
            close(wait)
            
            % create array of function coeffs
            function coeffs = create_function_coeffs_tbl(functions)
                
                coeffs = zeros(size(functions,2),numel(functions{1,1}.coeffs));
                for i = 1:size(functions,2)
                    coeffs(i,:) = functions{i,i}.coeffs;
                end
                
            end
            
            function Z = combine_contributions(operator_coeffs,test_index,basis_index,test_coeffs,basis_coeffs)
                
                % Loop through all interacting edges and multiply the basis and test polynomial coefficients with
                % the integrated convolution (operator) coefficients
                Z = 0;
                for alpha=1:numel(test_index)
                    for beta=1:numel(basis_index)
                        
                        operator_coeff = operator_coeffs{test_index(alpha),basis_index(beta)};
                        polynomial_coeffs = test_coeffs(alpha,:).' * basis_coeffs(beta,:);
                        Z = Z + sum(sum( polynomial_coeffs .* operator_coeff ));
                        
                    end
                end
            end
            
            
        end  % computeField
        

    end % public methods
    
    methods (Access = private)
        
        % Computation of single elements of all operator matrices
        function [coeffs_nh,coeffs_ns,coeffs_d,coeffs_s,coeffs_dp] = Z_calc( obj,...
                s_i,s_o,G_dd,G_SZ,G_ZS,G_SS,G_ZZ )
            
            % Find inner and outer Gaussian quadrature points
            inner_quad_points = numel(s_i);
            outer_quad_points = numel(s_o);
            
            % Make rho_m a vector of outer Gaussian quadrature points along observation edge
            rho_m = bsxfun(@plus, obj.a_m, s_o*(obj.b_m - obj.a_m));
            rho_m = reshape(rho_m,[outer_quad_points, 1, 2]);
            
            % Make rho_n a vector of inner Gaussian quadrature points along source edge
            rho_n = bsxfun(@plus, obj.a_n, s_i*(obj.b_n-obj.a_n));
            rho_n = reshape(rho_n,[1, inner_quad_points, 2]);
            
            % P is an array of distances from rho_n to observation point, rho_m
            rho_mn = bsxfun(@minus,rho_m,rho_n);
            P = sqrt(sum(rho_mn.^2,3));
            
            % Compute the temporal convolutions
            [Fh, F, dF] = BEUT.BEM.computeConvolutions(P/obj.c, obj.shiftedTB_Nh, obj.shiftedTB_Ns, obj.shiftedTB_D);
            dF=dF/obj.c;
            
            % dF/dn = n dot P/|P| * dF
            unit_P = bsxfun(@rdivide, rho_mn, P);                                       % P/|P|
            dF_dnp = sum(bsxfun(@times, reshape(-obj.n_n,[1,1,2]), unit_P),3) .* dF;    % dF/dn'
            dF_dn = sum(bsxfun(@times, reshape(obj.n_m,[1,1,2]), unit_P),3) .* dF;      % dF/dn
            
            % Perform quadrature
            nh_intF = sum(sum(bsxfun(@times,G_dd,Fh),1),2);
            ns_intF = sum(sum(bsxfun(@times,G_SS,F),1),2);
            dp_intF = sum(sum(bsxfun(@times,G_SZ,dF_dn),1),2);
            d_intF = sum(sum(bsxfun(@times,G_ZS,dF_dnp),1),2);
            s_intF = sum(sum(bsxfun(@times,G_ZZ,F),1),2);
            
            % Scale using edge lengths
            coeffs_nh = shiftdim(nh_intF,2) * obj.l_n * obj.l_m;
            coeffs_ns = shiftdim(ns_intF,2) * obj.l_n * obj.l_m * dot(obj.t_m,obj.t_n);
            coeffs_d  = shiftdim(d_intF,2)  * obj.l_n * obj.l_m;
            coeffs_s  = shiftdim(s_intF,2)  * obj.l_n * obj.l_m;
            coeffs_dp = shiftdim(dp_intF,2) * obj.l_n * obj.l_m;
            
        end % Z_calc
        
        
        % Computation of single elements of all operator matrices
        function [coeffs_d,coeffs_s] = Z_calc_field( obj,...
                rho_m, s_i,G_ZS,G_ZZ )
            
            % Find inner and outer Gaussian quadrature points
            inner_quad_points = numel(s_i);
            
            % Make rho_n a vector of inner Gaussian quadrature points along source edge
            rho_n = bsxfun(@plus, obj.a_n, s_i*(obj.b_n-obj.a_n));
            rho_n = reshape(rho_n,[1, inner_quad_points, 2]);
            rho_m = reshape(rho_m,[1, 1, 2]);
            
            % P is an array of distances from rho_n to observation point, rho_m
            rho_mn = bsxfun(@minus,rho_m,rho_n);
            P = sqrt(sum(rho_mn.^2,3));
            
            % Compute the temporal convolutions
            [~, F, dF] = BEUT.BEM.computeConvolutions(P/obj.c, obj.shiftedTB_Ns, obj.shiftedTB_Ns, obj.shiftedTB_D);
            
            % dF/dn = n dot P/|P| * dF
            unit_P = bsxfun(@rdivide, rho_mn, P);                                       % P/|P|
            dF_dnp = sum(bsxfun(@times, reshape(-obj.n_n,[1,1,2]), unit_P),3) .* dF;    % dF/dn'
            
            % Perform quadrature
            d_intF = sum(sum(bsxfun(@times,G_ZS,dF_dnp),1),2);
            s_intF = sum(sum(bsxfun(@times,G_ZZ,F),1),2);
            
            % Scale using edge lengths
            coeffs_d  = shiftdim(d_intF,2)  * obj.l_n;
            coeffs_s  = shiftdim(s_intF,2)  * obj.l_n;
            
        end % Z_calc_field
        
        
    end % private methods
    
end % class


