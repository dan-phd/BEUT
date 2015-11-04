function [ Fh, Fs, dF ] = ComputeConvolutions( distances, intTB, TB, dTB )
% Compute temporal convolutions for a temporal basis function and its integrated & derivative form
% P is made up of a number of distances (from source to observation)
% intTB, TB, dTB are the temporal basis functions for Fh, Fs, dF respectively
% Fh, Fs, dF are matrices of the same size as distances

P = reshape(distances,numel(distances),1);

num_partitions = length(intTB.partition)-1;
num_degree = size(intTB.coeffs,2)-1;

Ia = zeros(length(P),num_degree+1);
Ib = zeros(length(P),num_degree+1);
Fh = zeros(length(P),1);
Fs = zeros(length(P),1);
dF = zeros(length(P),1);

dIa=Ia; dIb=Ia;
for i = 1:num_partitions
    
    % Current partition properties
    partition_start = intTB.partition(i);     % (k-i)dt
    partition_end   = intTB.partition(i+1);   % (k-i+1)dt
    
    % Integration limits
    a = max(partition_start, P);
    b = max(partition_end, P);
    
    % Differentiated limits wrt P
    da = heaviside(P-partition_start);
    db = heaviside(P-partition_end);
    
    % Compute dI_1 term and set this term to 0 after threshold point
    dI1a = (a.*da - P)./sqrt(a.^2-P.^2); dI1a(P>partition_start)=0;
    dI1b = (b.*db - P)./sqrt(b.^2-P.^2); dI1b(P>partition_end)=0;
    
    % First integral equation, I_0 and dI_0
    Ia(:,1) = log(a + sqrt(a.^2 - P.^2));
    Ib(:,1) = log(b + sqrt(b.^2 - P.^2));
    dIa(:,1) = (da + dI1a)./(a+sqrt(a.^2-P.^2));
    dIb(:,1) = (db + dI1b)./(b+sqrt(b.^2-P.^2));
    
    % Convolve polynomial coefficients with solved integral
    Fh = Fh + ( Ib(:,1) -  Ia(:,1)) .* intTB.coeffs(i, num_degree+1);
    Fs = Fs + ( Ib(:,1) -  Ia(:,1)) .* TB.coeffs   (i, num_degree+1);
    dF = dF + (dIb(:,1) - dIa(:,1)) .* dTB.coeffs  (i, num_degree+1);
    
    % Second integral equation, I_1
    if num_degree>0
        
        % I_1
        Ia(:,2) = sqrt(a.^2 - P.^2);
        Ib(:,2) = sqrt(b.^2 - P.^2);
        
        % dI_1
        dIa(:,2) = dI1a;
        dIb(:,2) = dI1b;
        
        % Convolve
        Fh = Fh + ( Ib(:,2) -  Ia(:,2)) .* intTB.coeffs(i, num_degree);
        Fs = Fs + ( Ib(:,2) -  Ia(:,2)) .* TB.coeffs   (i, num_degree);
        dF = dF + (dIb(:,2) - dIa(:,2)) .* dTB.coeffs  (i, num_degree);
        
        
        % All other integral equations, I_2...I_p
        for n = 2:num_degree
            
            % I_2...I_p
            Ia(:,n+1) = (a.^(n-1).*Ia(:,2) + (n-1)*P.^2.*Ia(:,n-1))/n;
            Ib(:,n+1) = (b.^(n-1).*Ib(:,2) + (n-1)*P.^2.*Ib(:,n-1))/n;
            
            % dI_2...dI_p
            dIa(:,n+1) = (da.^(n-1).*Ia(:,2) + a.^(n-1).*dIa(:,2) + (n-1)*(2*P.*Ia(:,n-1) + P.^2.*dIa(:,n-1)))/n;
            dIb(:,n+1) = (db.^(n-1).*Ib(:,2) + b.^(n-1).*dIb(:,2) + (n-1)*(2*P.*Ib(:,n-1) + P.^2.*dIb(:,n-1)))/n;
            
            % Convolve
            Fh = Fh + ( Ib(:,n+1) -  Ia(:,n+1)) .* intTB.coeffs(i, num_degree-n+1);
            Fs = Fs + ( Ib(:,n+1) -  Ia(:,n+1)) .* TB.coeffs   (i, num_degree-n+1);
            dF = dF + (dIb(:,n+1) - dIa(:,n+1)) .* dTB.coeffs  (i, num_degree-n+1);
        end
    end
    
end

Fh = Fh/(2*pi); Fs = Fs/(2*pi);
dF = dF/(2*pi);
% dF(1:end-1) = diff(dF)./diff(P);        % Use Matlab diff

Fh = reshape(Fh,size(distances));
Fs = reshape(Fs,size(distances));
dF = reshape(dF,size(distances));

end
