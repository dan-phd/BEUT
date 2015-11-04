function plotCurrentDensityInFrequency(c, J_f, omega, spectrum, limit, ...
    struct_of_points, Analytical_J)
% Plot J(w) at points around the cylinder and compare with analytical solution if used.
% struct_of_points is a struct array, with each struct requiring
% a 'point' parameter corresponding to the observation point and 
% an optional 'name' parameter corresponding to a legend name

omega_axis = omega(1:limit);
Amplitude_response = abs(spectrum(1:limit) * c);

figure; hold on;
plot(omega_axis,abs(J_f(vertcat(struct_of_points.point),1:limit)));
title('Frequency domain current density');
xlabel('\omega'); ylabel('J(\rho,\omega)');
for i=1:numel(struct_of_points)
    names = legend;
    
    if isfield(struct_of_points(i),'name')
        if ~isempty(struct_of_points(i).name)
            legend([names.String strcat(struct_of_points(i).name,' - Numerical')]);
        else
            legend([names.String strcat(num2str(struct_of_points(i).point),' - Numerical')]);
        end
    else
        legend([names.String strcat(num2str(struct_of_points(i).point),' - Numerical')]);
    end
end

if nargin>6
    for i=1:numel(struct_of_points)
        Analytical_to_plot{i} = abs(Analytical_J(struct_of_points(i).point,1:limit));
        plot(omega_axis,Analytical_to_plot{i},'.');
        names = legend;
        if isfield(struct_of_points(i),'name')
            if ~isempty(struct_of_points(i).name)
                legend([names.String strcat(struct_of_points(i).name,' - Analytical')]);
            else
                legend([names.String strcat(num2str(struct_of_points(i).point),' - Analytical')]);
            end
        else
            legend([names.String strcat(num2str(struct_of_points(i).point),' - Analytical')]);
        end
    end
    
    
    % scale the amplitude response if needed
    plot_max = max(max(vertcat(Analytical_to_plot{:})));
    if plot_max > 2 || plot_max < 0.5
        Amplitude_response = Amplitude_response*1.1*plot_max;
    end
    
end

plot(omega_axis,Amplitude_response)
names = legend;
legend([names.String 'Excitation spectrum (\omega)'])

end
