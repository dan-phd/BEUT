function [ x ] = MOT( Z, V )
%Marching-on-in-time function to compute x given Z and V, where {V}={Z}[x]

Z0_inv = inv(Z(:,:,1));

N_T = size(Z,3);

x = zeros(size(V,1),N_T);
bar = waitbar(0,'Marching on in time...');
for j=0:N_T-1
    waitbar(j/N_T,bar);
    
    rhs_vec = V(:,j+1);
    
    for k=1:j-1
        rhs_vec = rhs_vec - Z(:,:,k+1)*x(:,j-k+1);
    end
    
    x(:,j+1) = Z0_inv * rhs_vec;
    
end
close(bar)

end

