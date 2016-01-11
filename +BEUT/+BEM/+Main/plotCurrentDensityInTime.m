function plotCurrentDensityInTime(time, J, struct_of_points)
% Plot current density at points around the cylinder
% struct_of_points is a struct array, with each struct requiring
% a 'point' parameter corresponding to the observation point and 
% an optional 'name' parameter corresponding to a legend name

figure;

t_limit = numel(time);
plot( time(1:t_limit), J(vertcat(struct_of_points.point),1:t_limit) )

for i=1:numel(struct_of_points)
    names = legend;
    
    if isfield(struct_of_points(i),'name')
        if ~isempty(struct_of_points(i).name)
            legend([names.String struct_of_points(i).name]);
        else
            legend([names.String num2str(struct_of_points(i).point)]);
        end
    else
        legend([names.String num2str(struct_of_points(i).point)]);
    end
end

xlabel('time'); ylabel('J(\rho,t)');

end
