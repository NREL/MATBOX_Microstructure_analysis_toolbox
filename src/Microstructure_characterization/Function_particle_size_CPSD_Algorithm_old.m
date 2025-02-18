function [Particle_size, distance_transform, IDX] = Function_particle_size_CPSD_Algorithm_old(binary_phase,smooth_surface)
% Return largest diameter in voxel length of the largest sphere that contain each voxel noted 1

%% NOTATION
% Binary phase
% =1 phase investigated
% =0 complementary phase


%% STEP 1: CALCULATE LARGEST SPHERE CENTERED IN EACH VOXEL
% We use the matlab built-in function bwdist to calculate the distance transform (euclidean)
% bwdist(BW) : For each pixel in BW, the distance transform assigns a number that is the distance between that pixel and the nearest nonzero pixel of BW
% then, bwdist(~binary_phase,'euclidean') :
[distance_transform, IDX] = bwdist(~binary_phase,'euclidean');

%distance_transform=double(bwdist(~binary_phase,'euclidean'));

%% STEP 2: ASSIGN VOXELS TO THE LARGEST SPHERE THAT CONTAINS THEM
% Definition of the continum spherical PSD
% In order to reduce CPU time, we will (i) check the condition "voxel belong to
% sphere" from largest sphere to smallest sphere and (ii) remove from
% further condition check all voxels which have been already assigned

tolerance=1e-6; % There is no need to go below: the value increment is around 1

% Initialize distance map
Domain_size=size(binary_phase); % Domain size of the microstructure
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end
Particle_size = zeros(Domain_size(1),Domain_size(2),Domain_size(3));

% Current maximum distance
current_distance = max(max(max(distance_transform)));
current_distance=double(current_distance);

% The iteration continues while the distance transform contains values superior to sqrt(3)
% = -1 voxels have been already assigned to largest possible sphere
% = 1, sqrt(2), sqrt(3) voxels belongs to sphere of size 1, sqrt(2) or sqrt(3)
% > sqrt(3) algorithm will check these voxels
while current_distance>sqrt(3)
    % Create distance map
    current_distance_map = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    % Voxels which are the center of the largest sphere are set equal to 1
    current_distance_map( abs(distance_transform-current_distance)<tolerance ) = 1;
    % Calculate distance in the whole domain
    % (inefficient if very few sheres have been selected with current_distance, efficient if a lot of spheres have been selected)
    current_distance_map=double(bwdist(current_distance_map,'euclidean'));
    % Assign all voxels that does not belong to a largest sphere
    cond1 = current_distance_map<(current_distance+tolerance);
    Particle_size ( cond1 & Particle_size<(current_distance+tolerance) ) = current_distance;
    % Update distance transform
    distance_transform( cond1 & Particle_size<(current_distance+tolerance) )=-1;
    % Update current maximum distance
    current_distance = double( max(max(max(distance_transform))) );
end

%%

% Missingvoxels = ~double(binary_phase).*double(Particle_size);
% Particle_size(Missingvoxels~=0)=10;
% Fig = figure; imagesc(binary_phase); axis equal; axis tight; colormap gray;
% Fig = figure; imagesc(Particle_size); axis equal; axis tight; colormap turbo;
% 
% 
% keyboard
% 
% 
% cmap = turbo(length(unique(Particle_size)));
% cmap(1,:) = [0.5 0.5 0.5];
% Fig = figure; imagesc(Particle_size); axis equal; axis tight; colormap(cmap);
% ddd

%% ALLOCATING THE VOXELS WITH A DISTANCE < SQRT(3)
% The unassigned voxels of the investigated phase are contained within a sphere a diameter 1, sqrt(2) or sqrt(3)
% We will apply a commom number for particle of size sqrt(3)
Particle_size( abs(distance_transform-sqrt(3))<tolerance ) = sqrt(3);
% We will apply a commom number for particle of size sqrt(2)
Particle_size( abs(distance_transform-sqrt(2))<tolerance ) = sqrt(2);
% We will apply a commom number for particle of size sqrt(1)
Particle_size( abs(distance_transform-1)<tolerance ) = 1;
% Them the three last particles id corresponds to the spheres of diameter
% sqrt(3),sqrt(2) and 1
% The 9e9 distance value for the complemetray phase is removed
Particle_size(binary_phase==0)=0;

%% CORRECTION
do_correction = false; % Correction provides an underestimation. Do not use for now.

if do_correction
    % 1 0 0 0 1 % Binary phase
    % ->
    % 0 1 2 1 0 % Particle size (~radius)
    % -> diameter = 2*2 = 4 ? not accurate
    %    Correct diameter is 2*2 -0.5-0.5 = 3

    distance_values = unique(Particle_size);
    distance_values=[distance_values zeros(length(distance_values),1)];

    % For distance <=sqrt(3): no change diameter=2*distance-distance
    % No error
    linear_index_0 =find( distance_values(:,1)<(sqrt(3)+tolerance));

    % For distance=integer: diameter=2*integer-1=2*distance-1
    % No error on the overlapped diameter
    linear_index_1 =find( abs(distance_values(:,1)-round(distance_values(:,1)))<tolerance & distance_values(:,1)>(sqrt(3)+tolerance));

    % For distance=integer*sqrt(2): diameter=2*integer*sqrt(2)-sqrt(2) = =2*distance-sqrt(2)
    % No error on the overlapped diameter
    linear_index_2 =find( abs(distance_values(:,1)/sqrt(2)-round(distance_values(:,1)/sqrt(2)))<tolerance & distance_values(:,1)>(sqrt(3)+tolerance));

    % For distance=integer*sqrt(3): diameter=2*integer*sqrt(3)-sqrt(3) = =2*distance-sqrt(3)
    % No error on the overlapped diameter
    linear_index_3 =find( abs(distance_values(:,1)/sqrt(3)-round(distance_values(:,1)/sqrt(3)))<tolerance & distance_values(:,1)>(sqrt(3)+tolerance));

    % New distances are:
    distance_values(linear_index_0,2)=distance_values(linear_index_0,1);
    distance_values(linear_index_1,2)=2*distance_values(linear_index_1,1)-1;
    distance_values(linear_index_2,2)=2*distance_values(linear_index_2,1)-sqrt(2);
    distance_values(linear_index_3,2)=2*distance_values(linear_index_3,1)-sqrt(3);
    % For others distance, diameter=2*distance-(1+sqrt(3))/2
    % Actual overlapped diameter range from 1 to sqrt(3), we choose the mean value
    linear_index_4 = find( distance_values(:,2)==0 & distance_values(:,1)>(sqrt(3)+tolerance));
    distance_values(linear_index_4,2)=2*distance_values(linear_index_4,1)-((1+sqrt(3))/2);

    % Apply corrections
    % We first grab all the index, and then apply the changes, so that we assure
    % each voxel get only 1 change as desired
    for current_correction = 1:1:length(distance_values(:,1))
        index.correction(current_correction).value = find( Particle_size==distance_values(current_correction,1));
    end
    for current_correction = 1:1:length(distance_values(:,1))
        Particle_size(index.correction(current_correction).value)=distance_values(current_correction,2);
    end

else
    Particle_size = 2*Particle_size; % Diameter
end

%% SMOOTHING
% Pixels at 1-2 from the interface have numerical errors
if smooth_surface
    dmap_chessboard=double(bwdist(~binary_phase,'chessboard'));
    cond1 = dmap_chessboard<=2;
    idx_toreplace = find(cond1.*binary_phase==1);
    tmp = zeros(size(binary_phase));
    tmp(idx_toreplace)=1;
    tmp2=tmp+~binary_phase;

    % Re-assign voxels to nearest correct particle diameter values
    [~, IDX] = bwdist(~tmp2);
    Particle_size(idx_toreplace) = Particle_size(IDX(idx_toreplace));
end

% cmap = turbo(length(unique(Particle_size)));
% cmap(1,:) = [0.5 0.5 0.5];
% Fig = figure; imagesc(Particle_size); axis equal; axis tight; colormap(cmap);

end