function [shortest_distances_opposite_from_source_to_target,shortest_distances_opposite_from_target_to_source,slice_start,slice_end]=Function_shortest_distance_from_skeleton(G,node_coordinate,Domain_size,voxel_size,direction,Particle_label);

%% FIND ALL PARTICLE ID AT THE EDGES

% Get slice at the edges, normal to the direction
if direction==1
    slice_start = zeros(Domain_size(2),Domain_size(3));
    slice_end = zeros(Domain_size(2),Domain_size(3));
    slice_start(:,:) = Particle_label(1,:,:);
    slice_end(:,:)   = Particle_label(end,:,:);
elseif direction==2
    slice_start = zeros(Domain_size(1),Domain_size(3));
    slice_end = zeros(Domain_size(1),Domain_size(3));    
    slice_start(:,:) = Particle_label(:,1,:);
    slice_end(:,:)   = Particle_label(:,end,:);    
elseif direction==3
    slice_start = zeros(Domain_size(1),Domain_size(2));
    slice_end = zeros(Domain_size(1),Domain_size(2));     
    slice_start(:,:) = Particle_label(:,:,1);
    slice_end(:,:)   = Particle_label(:,:,end);    
end

% Get unique particle id at the edges, normal to the direction
id_start = unique(slice_start);
id_end = unique(slice_end);

% Remove 0
id_start(id_start==0)=[];
id_end(id_end==0)=[];

% Number of particles
n_start = length(id_start);
n_end = length(id_end);

%% CALCULATE DISTANCES

% % From source to target

% Dijkstra algorithm 
all_distances_opposite_faces = distances(G,id_start,id_end);
% Note, we may have some +Inf values. It occurs when two nodes does not
% belong to the same graph component (i.e. not connected).

% We need to add the distance from each particle to the edge
for current_particle = 1:1:n_start
    current_node = id_start(current_particle);
    distance_from_start_edge =  node_coordinate(current_node,direction+1);
    all_distances_opposite_faces(current_particle,:)=all_distances_opposite_faces(current_particle,:)+distance_from_start_edge;
end
for current_particle = 1:1:n_end
    current_node = id_end(current_particle);
    distance_from_start_edge =  node_coordinate(current_node,direction+1);
    distance_from_end_edge= (Domain_size(direction)*voxel_size/1000)-distance_from_start_edge;
    all_distances_opposite_faces(:,current_particle)=all_distances_opposite_faces(:,current_particle)+distance_from_end_edge;
end

% Keep only the shortest distance from each source node to the target face
shortest_distances_opposite_from_source_to_target = min(all_distances_opposite_faces,[],2) ;

% Normalize with the distance
shortest_distances_opposite_from_source_to_target=shortest_distances_opposite_from_source_to_target/(Domain_size(direction)*voxel_size/1000);

% Get coordinates (voxel) as well
for current_particle = 1:1:n_start
    current_node = id_start(current_particle);
    slice_start(slice_start==current_node)=shortest_distances_opposite_from_source_to_target(current_particle);
end


% % From target to source

% Dijkstra algorithm 
all_distances_opposite_faces_opposite = distances(G,id_end,id_start);
% Note, we may have some +Inf values. It occurs when two nodes does not
% belong to the same graph component (i.e. not connected).

% We need to add the distance from each particle to the edge
for current_particle = 1:1:n_start
    current_node = id_start(current_particle);
    distance_from_start_edge =  node_coordinate(current_node,direction+1);
    all_distances_opposite_faces_opposite(:,current_particle)=all_distances_opposite_faces_opposite(:,current_particle)+distance_from_start_edge;
end
for current_particle = 1:1:n_end
    current_node = id_end(current_particle);
    distance_from_start_edge =  node_coordinate(current_node,direction+1);
    distance_from_end_edge= (Domain_size(direction)*voxel_size/1000)-distance_from_start_edge;
    all_distances_opposite_faces_opposite(current_particle,:)=all_distances_opposite_faces_opposite(current_particle,:)+distance_from_end_edge;
end

% Keep only the shortest distance from each source node to the target face
shortest_distances_opposite_from_target_to_source = min(all_distances_opposite_faces_opposite,[],2) ;

% Normalize with the distance
shortest_distances_opposite_from_target_to_source=shortest_distances_opposite_from_target_to_source/(Domain_size(direction)*voxel_size/1000);

% Get coordinates (voxel) as well
for current_particle = 1:1:n_end
    current_node = id_end(current_particle);
    slice_end(slice_end==current_node)=shortest_distances_opposite_from_target_to_source(current_particle);
end


end

