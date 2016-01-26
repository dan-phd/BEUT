classdef GramMatrix
% compute the Gram matrix and output a matrix of size(number of edges)^2.
    
    properties
        basis_function;     % the function used for sampling, in the form of a BasisFunction
        test_function;      % the function used for testing, in the form of a BasisFunction
        geometry;           % a list of halfedges in the form of MeshBoundary.halfedges or 
                            % MeshBoundary.dual
        test_points = 3;    % how many Gaussian quadrature points to use per edge
    end
    
    methods
        
        function obj = GramMatrix()
            
        end % constructor
        
        function G = compute(obj)
            
            % Concatenate the test function coefficients and find max degree
            test_coeffs = vertcat(obj.test_function.pol{:});
            max_deg_test = max(vertcat(test_coeffs.degree));
            basis_coeffs = vertcat(obj.basis_function.pol{:});
            max_deg_basis = max(vertcat(basis_coeffs.degree));
            
            % Gaussian quadrature integration
            [s,w] = BEUT.BEM.lgquad(obj.test_points,[0 1]);
            
            % make ready to multiply with spatial basis functions
            G_coeffs = BEUT.BEM.createGaussianCoeffs(s,w,s,w.^0,max_deg_test,max_deg_basis);
            G_coeffs = bsxfun(@times,G_coeffs,eye(obj.test_points));        % Gram matrix operates when m=n
            G_coeffs = shiftdim(sum(sum(G_coeffs)),2);
            
            % Index table to determine which field coefficient to use
            test_idx = obj.test_function.idx_table;
            basis_idx = obj.basis_function.idx_table;
            assert(size(test_idx,1)==size(basis_idx,1));
            N_F = size(test_idx,1);
            N_E = numel(obj.geometry);
            
            % Interact test and basis quadrature points
            C = num2cell(zeros(N_E, N_E));
            for m = 1:N_E
                
                l_m = obj.geometry(m).l;
                
                C{m,m} = G_coeffs * l_m * l_m;
                    
            end
            
            G = zeros(N_F,N_F);
            for p = 1:N_F
                
                % Current test segments
                index_test = test_idx(p,:);
                
                % Current test function coefficients
                TF = obj.test_function.pol(:,index_test);
                TC = zeros(size(TF,2),numel(TF{1,1}.coeffs));
                for i = 1:size(TF,2)
                    TC(i,:) = TF{i,i}.coeffs;
                end
                
                % Gram matrix only operates on self-patch elements
                for q = 1:N_F;
                    
                    % Current basis segments
                    index_basis = basis_idx(q,:);
                    
                    % Current basis function coefficients
                    BF = obj.basis_function.pol(:,index_basis);
                    BC = zeros(size(BF,2),numel(BF{1,1}.coeffs));
                    for i = 1:size(BF,2)
                        BC(i,:) = BF{i,i}.coeffs;
                    end
                    
                    
                    % Loop through all interacting edges and multiply the basis and test
                    % polynomial coefficients with the Gram coefficients
                    for alpha=1:numel(index_test)
                        for beta=1:numel(index_basis)
                            
                            G_coeffs = C{index_test(alpha),index_basis(beta)};
                            polynomial_coeffs = TC(alpha,:).' * BC(beta,:);
                            G(p,q) = G(p,q) + sum(sum( polynomial_coeffs .* G_coeffs ));
                            
                        end % alpha
                    end % beta
                    
                end % q
            end % p
            
        end % compute
        
        
    end
end

