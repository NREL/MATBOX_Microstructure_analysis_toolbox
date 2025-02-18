function [Label_lake] = Function_particle_size_Watershed_algorithm(BW,details_convergence)

%% REFERENCE
% This algorithm is inspired by this publication:
% Watersheds in  Digital  Spaces: An  Efficient  Algorithm  Based on  Immersion  Simulations 
% Luc  Vincent  and  Pierre  Soille
% IEEE  TRANSACTIONS  ON  PATTERN  ANALYSIS  AND  MACHINE  INTELLIGENCE, VOL.  13,  NO.  6,  JUNE  1991

%% TOPOGRAPHY ANALOGY
dmap=bwdist(~BW,'euclidean'); % Euclidean distance map
dmap = -dmap; % Inverse the distance sign: 0 is the reference altitude, minimum is the lower altitude point
dmap(~BW) = +Inf; % A negative infinite value is attributed to the complementary

%% SORTING THE DISTANCE MAP VALUES
all_altitude = unique(dmap); % Get all altitudes
all_altitude = sortrows(all_altitude,1); % Matlab has already sorted this array from -inf to +inf. Let's be sure of this
all_altitude(end)=[]; % Remove last value value (= +Inf)
number_altitude = length(all_altitude); % Number of altitude level

%% LABELING LAKE ITERATIVELY
Domain_size=size(BW);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end

Label_lake = zeros(Domain_size(1), Domain_size(2), Domain_size(3)); % Initialisation
if details_convergence
    disp 'Flooding simulation in progress...(%)';
end
next_step_display=1;
for current_altitude_iteration = 1:1:number_altitude % We start from the lower position and go up
    % Progress of the calculation in percents
    progress_force_calculation = 100*current_altitude_iteration/number_altitude;
    if progress_force_calculation>=next_step_display && details_convergence
        fprintf (' %4.1f',progress_force_calculation);
        next_step_display=next_step_display+1;
    end
    current_altitude = all_altitude(current_altitude_iteration); % Current altitude value
    
    % Binary image of this level of altitude
    binary_altitude_map = zeros(Domain_size(1), Domain_size(2), Domain_size(3)); % Initialisation
    binary_altitude_map(dmap<=current_altitude)=1; % Mark 1 when voxel are located at the current altitude OR LOWER: LAKES ARE BEING FILLING FROM THEIR DEEPEST POINT
    
    % Label the current level image
    Label_current_altitude = bwlabeln(binary_altitude_map,6);
    unique_cluster = unique(Label_current_altitude); % Number of different cluster identified
    unique_cluster(1)=[]; % Remove 0 (complementary phase)
    number_cluster = length(unique_cluster); % Number of different cluster

    % We label the different lake
    if current_altitude_iteration>1
        for current_iteration = 1:1:number_cluster % Loop over all the cluster identified
            Current_label = unique_cluster(current_iteration); % Current label id
            binary_cluster = zeros(Domain_size(1), Domain_size(2), Domain_size(3)); % Create binary cluster
            binary_cluster(Label_current_altitude==Current_label)=1;
            index_cluster = find(binary_cluster==1); % Index voxels of the cluster
            intersection = binary_cluster.*Label_lake; % Calculate intersection with the already identified lake
            index_intersection = find(intersection~=0); % Index of the intersected voxels
            number_intersection = length(index_intersection); % Number of voxels intersected
            
            % % % Four cases are possible
            
            % % Case 1 and 2: intersection is empty AND...
            if number_intersection==0
                % Binary image of the existing lake
                binary_existing_lake = zeros(Domain_size(1), Domain_size(2), Domain_size(3));
                binary_existing_lake(Label_lake~=0)=1;
                % Distance bewteen all non attributes voxels and the existing lakes
                [Dist_existing_lake,Index_nearest_lake] = bwdist(binary_existing_lake,'chessboard');
                % Restrict result to the current new cluster
                Dist_existing_lake=Dist_existing_lake.*binary_cluster;
                % All distance
                unique_distance = unique(Dist_existing_lake);
                % Remove 0
                unique_distance(1)=[];
                % Minimum distance
                min_distance = unique_distance(1);
                if min_distance>2
                    % Case 1: AND... far from existing lake: this cluster is a new lake
                    % Set new id for this new lake
                    max_lake_id = max(unique(Label_lake))+1;
                    % Attribute this new id to the cluster
                    Label_lake(index_cluster)=max_lake_id;
                else
                    % Case 2: AND... 'adjacent' with one or more existing lake
                    % Index of voxels of the cluster close to the existing lake
                    index_near = find(Dist_existing_lake==min_distance);
                    % Id of these lake
                    index_adjacent_lake = Index_nearest_lake(index_near);
                    id_adjacent_lake = Label_lake(index_adjacent_lake);
                    % Number of adjacent lake
                    number_adjacent_lake = length(id_adjacent_lake);
                    if number_adjacent_lake==1
                        % Simple case: Attribute all voxels of the label to this unique lake
                        Label_lake(index_cluster)=id_adjacent_lake(1);
                    else
                        % Voxel will be attributed to the closest lake
                        for current_voxel=1:1:length(index_cluster)
                            choose_lake = zeros(number_adjacent_lake,8);
                            [x_tmp,y_tmp,z_tmp] = ind2sub(Domain_size,index_cluster(current_voxel));
                            choose_lake(:,1)=x_tmp;
                            choose_lake(:,2)=y_tmp;
                            choose_lake(:,3)=z_tmp;
                            for current_lake=1:1:number_adjacent_lake
                                [x_tmp,y_tmp,z_tmp] = ind2sub(Domain_size,index_adjacent_lake(current_lake));
                                choose_lake(current_lake,4)=x_tmp;
                                choose_lake(current_lake,5)=y_tmp;
                                choose_lake(current_lake,6)=z_tmp;
                                choose_lake(current_lake,7)=Label_lake(index_adjacent_lake(current_lake));
                            end
                            % Calculate distance
                            choose_lake(:,8)= (choose_lake(:,1)-choose_lake(:,4)).^2 + (choose_lake(:,2)-choose_lake(:,5)).^2 + (choose_lake(:,3)-choose_lake(:,5)).^2;
                            % Sort by increasing distance
                            choose_lake = sortrows(choose_lake,8);
                            % Attribute the current voxel
                            Label_lake(index_cluster(current_voxel))=choose_lake(1,7);
                        end
                    end
                end
            end
            
            % % Case 3 and 4: intersection is not empty AND...
            if number_intersection>0
                % Does intersected voxels belong to different lake?
                different_lake = unique(Label_lake(index_intersection));
                % Number different lake
                number_different_lake = length(different_lake);
                % % Case 3: AND... intersection with only one existing lake
                if number_different_lake==1
                    % Attribute all voxels of the label to this unique lake
                    Label_lake(index_cluster)=different_lake;
                end
                % % Case 4: AND... intersection with more than 1 existing lake
                if number_different_lake>1
                    % Each voxel of the cluster will be attributed to the closest lake
                    % % Step 1: Get subdomain that contains the cluster
                    % Get all coordinates
                    [C1,C2,C3] = ind2sub(Domain_size,index_cluster);
                    % Min, max coordinates
                    x_min = min(C1); y_min = min(C2); z_min = min(C3);
                    x_max = max(C1); y_max = max(C2); z_max = max(C3);
                    % Extract subdomain for the binary cluster
                    subdomain_cluster = zeros(x_max-x_min+1,y_max-y_min+1,z_max-z_min+1);
                    subdomain_cluster = binary_cluster(x_min:x_max,y_min:y_max,z_min:z_max);
                    % Extract subdomain for the lake
                    subdomain_lake = zeros(x_max-x_min+1,y_max-y_min+1,z_max-z_min+1);
                    subdomain_lake = Label_lake(x_min:x_max,y_min:y_max,z_min:z_max);
                    % Extract subdomain for the intersection
                    subdomain_intersection = zeros(x_max-x_min+1,y_max-y_min+1,z_max-z_min+1);
                    subdomain_intersection = intersection(x_min:x_max,y_min:y_max,z_min:z_max);
                    % % Step 2: Prepare the matrix that will calculate the shortest distance
                    subdomain_distance = subdomain_intersection;
                    subdomain_distance (subdomain_distance~=0)=1;
                    % Step 3: Calculate shortest distance from voxel to lake
                    [~,Distance_map_nearest_lake] = bwdist(subdomain_distance,'euclidean');
                    % Step 4: Mark the voxels that does not belong to the current cluster
                    Distance_map_nearest_lake=double(Distance_map_nearest_lake);
                    Distance_map_nearest_lake(subdomain_cluster==0)=-1;
                    % Step 5: Attribute to each voxel of the cluster the nearest lake
                    attraction_point = unique(Distance_map_nearest_lake);
                    number_attraction_point = length(attraction_point);
                    for current_=1:1:number_attraction_point
                        current_attraction_point = attraction_point(current_);
                        if current_attraction_point~= -1
                            % Find all voxels attracted by this point
                            index_attracted = find(Distance_map_nearest_lake==current_attraction_point);
                            % Attribute these voxel to the corresponding lake
                            subdomain_lake(index_attracted)=subdomain_lake(current_attraction_point);
                        end
                    end
                    % Step 6: Re-introduced the subomain to the whole domain
                    Label_lake(x_min:x_max,y_min:y_max,z_min:z_max)=subdomain_lake;
                end
            end
        end
    else
        % Firt iteration?
        % Then all the different cluster at the minimum voxel identifies
        % different lake (their lower altitude)
        Label_lake = Label_current_altitude;
    end
end

% Get all discrete particle id
unique_particle = unique(Label_lake);
unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
% Get number of particle
number_particle = length(unique_particle);
if details_convergence
    fprintf ('Initial number of discrete particle identified: %i\n',number_particle);
end
% 
% %% FINAL DISCRETE PARTICLE SIZE
% % Initialisation
% D_PSD_particle_size = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% % Unique id
% unique_particle = unique(Label_lake);
% unique_particle(1)=[];
% % Number of different particle
% number_particle = length(unique_particle);
% % loop over all particles
% for current_=1:1:number_particle
%     % Get id_
%     particle_id = unique_particle(current_);
%     % Find all voxels of the particle
%     index_particle=find(Label_lake==particle_id);
%     % Number of voxel
%     number_voxel_particle = length(index_particle);
%     volume_ = number_voxel_particle;
%     area_ = number_voxel_particle;
%     % Equivalent diameter size
%     if Domain_size(3)>1
%         % 3D case
%         equivalent_diameter_size= 2 * ((3*volume_  /(4*pi))^(1/3));
%     else
%         % 2D case
%         equivalent_diameter_size= 2 * ((area_/pi)^(1/2));
%     end
%     D_PSD_particle_size(index_particle)= equivalent_diameter_size;
% end

end