function [ f, sorted_frequencies ] = findCylinderResonantFrequencies(radius,c,M,N)
% Calculate theoretical resonances for a cylindrical cavity
% of radius 'a' and propagation 'c'
% up to M orders and N Bessel roots/zeros

%% Find roots of Bessel function of first kind of order m
X = zeros(N,M+1);
for m = 0:M;  % order
    J0 = @(x) besselj(m,x);
    r = fzeros(J0);
    r = r(r>1);     % don't count the negligible roots
    X(:,m+1) = r(1:N);
end

%% Find resonant frequencies of PEC cylinder using Balanis - Advanced Engineering Electromagnetics eq. 9-27
f = (X * c)/( 2*pi*radius );

% print results in ascending order
sorted_frequencies = reshape(f,1,[]);
sorted_frequencies = sort(sorted_frequencies)';


    function found_roots = fzeros(F)
        
        interval = [0 50];
        number_points = 100;
        incremental_intervals = linspace(interval(1),interval(2),number_points);
        found_roots = [];
        for i=1:number_points-1
            try
                found_roots(end+1) = fzero(F,[incremental_intervals(i),incremental_intervals(i+1)]);
            end
        end
        
    end

end
