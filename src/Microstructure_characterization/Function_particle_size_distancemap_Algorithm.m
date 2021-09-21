function [distance_transform, fitted_diameter, dist_fct, analytical, error_radius] = Function_particle_size_distancemap_Algorithm(binary_phase, voxel_size)
% Fit distance map cumulative function with distance map cumulative
% function of a sphere of diameter D

%% NOTATION
% Binary phase
% =1 phase investigated
% =0 complementary phase
Domain_size=size(binary_phase); % Domain size of the microstructure

%% EUCLIDEAN DISTANCE MAP
distance_transform=double(bwdist(~binary_phase,'euclidean')); % Euclidean distance map
distance_transform=distance_transform-0.5; % From voxel center to surface
distance_transform(distance_transform<0)=0; % Complementary volume
% Unit conversion
distance_transform = distance_transform*voxel_size/1000; % nm->um

distance_transform_vector=reshape(distance_transform,[prod(Domain_size),1]); % Convert in a vector
distance_transform_vector(distance_transform_vector<=0)=[]; % % Remove all 0

%% CUMULATIVE AND DISTRIBUTION FUNCTIONS
% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;
[dist_fct, ~] = Function_probability_density(distance_transform_vector,[],density_fct_parameters);
distance = dist_fct.cumulative_fct(:,1);
cumulative= dist_fct.cumulative_fct(:,2);

%% FIT DIAMETER FOR THE ANALYTICAL CUMULATIVE FUNCTIONS FOR A SPHERE OF DIAMETER F
number_tested_radius=400;
number_of_point = 400; % Numerical resolution
min_radius=0.1*voxel_size/1000;
max_radius=2*max(distance);
error_radius=zeros(number_tested_radius,2);
tested_radius = linspace(min_radius,max_radius,number_tested_radius);

if length(distance)>=2
    for current_iteration=1:1:number_tested_radius
        % Sphere
        current_radius=tested_radius(current_iteration);
        volume_sphere=4/3*pi*(current_radius^3);
        
        % All r values between 0 and the sphere radius
        allradius=linspace(0,current_radius,number_of_point);
        
        % Maximum distance to boundary for each radius
        max_distance_to_boundary=current_radius-allradius;
        % Volume with a distance to boundary superior to max_distance_to_boundary
        volume_with_distancetoboundary_superior = (4/3*pi*(allradius.^3))/volume_sphere;
        
        % The cumulative function is max_distance_to_boundary=f(volume_with_distancetoboundary_superior)
        
        % To be compared the two functions need to share the same x-axis,
        common_xaxis=linspace(0,max(distance),2*number_tested_radius);
        
        % The analytical function: above its max, value is equal to 0
        analytical_ = interp1(max_distance_to_boundary,volume_with_distancetoboundary_superior,common_xaxis,'linear',0);
        % The numerical function: below its min, value is equal to 1
        numerical_ = interp1(distance,cumulative,common_xaxis,'linear',1);
        
        difference_=analytical_-numerical_; % Difference between analytical and numerical
        integral_=trapz(common_xaxis,difference_); % integral
        error_radius(current_iteration,2)=abs(integral_); % Absolute error...
        error_radius(current_iteration,1)=current_radius; % ...obtained with this radius
    end
    % Find the mimimum
    index_=find(error_radius(:,2)==min(error_radius(:,2)));
    fitted_radius=error_radius(index_,1);
    fitted_diameter=2*fitted_radius;
    % Recalculate the fitted cumulative fct.
    allradius=linspace(0,fitted_radius,number_of_point);
    volume_sphere=4/3*pi*(fitted_radius^3);
    max_distance_to_boundary=fitted_radius-allradius;
    volume_with_distancetoboundary_superior = (4/3*pi*(allradius.^3))/volume_sphere;
    
    analytical=zeros(length(max_distance_to_boundary),2);
    analytical(:,1)=max_distance_to_boundary;
    analytical(:,2)=volume_with_distancetoboundary_superior;
    
else
    fitted_diameter = NaN;
    analytical=NaN;
    error_radius=NaN;
end

end