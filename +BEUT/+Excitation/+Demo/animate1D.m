function animate1D(wave, time, domain)

N_T = length(time);

field = zeros(numel(domain),N_T);

for j=1:N_T
    
    field(:,j) = wave.eval(time(j),domain');
    
end

BEUT.animate_fields(1,'domain',domain, '1D wave',field);


end