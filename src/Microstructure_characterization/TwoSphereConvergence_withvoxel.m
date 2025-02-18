function [] = TwoSphereConvergence_withvoxel(maxdiameter,distance_btw_center)
    arguments
        maxdiameter {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(maxdiameter,9)}
    end


%% DIAMETER LIST
if (rem(maxdiameter, 2) == 0) % Even
    list_odd_voxelperdiameter = 7:2:maxdiameter+1; % n+1+n
    list_even_voxelperdiameter = 6:2:maxdiameter; % n+n
else
    list_odd_voxelperdiameter = 7:2:maxdiameter; % n+1+n
    list_even_voxelperdiameter = 6:2:maxdiameter-1; % n+n
end

%% INITIALIZATION
n_odd = length(list_odd_voxelperdiameter);
n_even = length(list_even_voxelperdiameter);

% Volume (area)
numerical_odd_spherevolume = zeros(n_odd,2);
numerical_even_spherevolume = zeros(n_even,2); 
numerical_odd_discarea = zeros(n_odd,2);
numerical_even_discarea = zeros(n_even,2); 


end