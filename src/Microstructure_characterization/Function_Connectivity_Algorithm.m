function [Connectivity_structure] = Function_Connectivity_Algorithm(binary_array, voxel_size, connected_id, unknown_id, isolated_id)
 
Domain_size = size(binary_array); % Domain_size
Numbervoxel_domain = numel(binary_array); % Number voxel domain
Numbervoxel_phase =  sum(sum(sum(binary_array==1))); % Number voxel phase

%% CLUSTERS DETECTION
% L = bwlabeln(BW,6)
% L is containing labels for the connected components in BW
% L is the same size as BW
% The elements of L are integer values greater than or equal to 0.
% The pixels labeled 0 are the background. The pixels labeled 1 make up one connected cluster. The pixels labeled 2 make up a second connected cluster, and so on
Clusters_tmp = bwlabeln(binary_array,6);

%% CALCULATE CLUSTER SIZE
unique_cluster = unique(Clusters_tmp); % Unique cluster
unique_cluster(unique_cluster==0)=[]; % Remove background (0) if any      
number_cluster = length(unique_cluster); % Number of cluster
cluster_size = zeros(number_cluster,17); % Initialise
% Column
% 1 : cluster number of voxel
% 2 : cluster volume in um3
% 3 : cluster volume fraction relative to domain volume, in percents
% 4 : cluster volume fraction relative to phase volume, in percents
% 5 : centroid position x in physical length
% 6 : centroid position y in physical length
% 7 : centroid position z in physical length
% 8 : 1= largest cluster, 2=unknown cluster, 3=isolated cluster
% 9 : 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face2face connection direction 1)
% 10: 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face2face connection direction 2)
% 11: 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face2face connection direction 3)
% 12 : 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face connection 1)
% 13: 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face connection 1)
% 14: 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face connection 1)
% 15 : 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face2face connection direction 1)
% 16: 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face2face connection direction 2)
% 17: 1= connected cluster, 2=unknown cluster, 3=isolated cluster (face2face connection direction 3)

Connectivity_structure.connected_id=connected_id;
Connectivity_structure.unknown_id=unknown_id;
Connectivity_structure.isolated_id=isolated_id;

for k=1:1:number_cluster % Loop over all clusters
    cluster_size(k,1) = unique_cluster(k); % Cluster id
    cluster_size(k,2) = sum(sum(sum(Clusters_tmp==cluster_size(k,1)))); % Number of voxel for current cluster
end
cluster_size(:,3) = cluster_size(:,2) * (voxel_size^3); % Volume of cluster in physical length
cluster_size(:,4) = 100 * cluster_size(:,2) ./ Numbervoxel_domain; % Volume fraction of cluster in regards with the domain, in percents
cluster_size(:,5) = 100 * cluster_size(:,2) ./ Numbervoxel_phase; % Volume fraction of cluster in regards with the phase, in percents
        
%% SORT CLUSTER PER SIZE
% 0: complementary phase
% 1: largest connected cluster
% 2: second largest connected cluster
% etc.
Clusters_ = zeros(Domain_size,'uint16'); % Initialise
cluster_size = sortrows(cluster_size,-2); % Sort, by decreasing order, the number of voxel
for k=1:1:number_cluster
    Clusters_(Clusters_tmp==cluster_size(k,1))=k;
end
cluster_size(:,1)=[]; % Remove old numerotation
clear Clusters_tmp % Delete old array
% Save results in structure
Connectivity_structure.Clusters_sortedpersize_3Darray = Clusters_;

%% HOW MANY ILL-DETAILED CLUSTERS ?

illdetailed_cluster_number = sum(cluster_size(:,1)==1); % Number of cluster with a size of 1 voxel
illdetailed_cluster_domainvolumefraction = 100*illdetailed_cluster_number/Numbervoxel_domain; % percent
illdetailed_cluster_phasevolumefraction = 100*illdetailed_cluster_number/Numbervoxel_phase; % percent

% Save results in structure
Connectivity_structure.illdetailed_cluster_number = illdetailed_cluster_number;
Connectivity_structure.illdetailed_cluster_domainvolumefraction = illdetailed_cluster_domainvolumefraction;
Connectivity_structure.illdetailed_cluster_phasevolumefraction = illdetailed_cluster_phasevolumefraction;

%% CLUSTERS IDENTIFICATION
% Main cluster is the largest cluster: =1
% Unknown clusters are not connected with the largest clusters but touch the domain's boundaries: =2
% Isolated clusters are not connected with the largest clusters and do not touch the domain's boundaries: =3
Clusters_LargestIsolatedUnknown.array = zeros(Domain_size,'uint8'); % Initialise

% All clusters that connect the two opposite faces normal to the investigated direction
% Cluster connect the two faces: =1
% Cluster does not connect the faces, but touch the domain's boundaryies: =2
% Cluster does not connect the faces, and do not touch the domain's boundaryies: =3
for current_direction=1:1:3 % Loop on every direction
    Clusters_TransportFace2Face.direction(current_direction).array = zeros(Domain_size,'uint8'); % Initialise
    Clusters_TransportFromFace1.direction(current_direction).array = zeros(Domain_size,'uint8'); % Initialise
    Clusters_TransportFromFace2.direction(current_direction).array = zeros(Domain_size,'uint8'); % Initialise
end

Clusters_LargestIsolatedUnknown.array(Clusters_==1)=connected_id; % Main (largest) cluster
cluster_size(1,8)=1;
for k=1:1:number_cluster % Loop over all clusters
    I=find(Clusters_== k); % Linear index of the voxels which belong to the current cluster
    [II,JJ,KK]=ind2sub([Domain_size(1),Domain_size(2),Domain_size(3)],I); % Coordinnates of these voxels
    % Finding extremums
    x_min = min(II); x_max = max(II);
    y_min = min(JJ); y_max = max(JJ);
    z_min = min(KK); z_max = max(KK);
    % center of mass (centroid) in physical length
    cluster_size(k,5) = mean(II)*voxel_size;
    cluster_size(k,6) = mean(JJ)*voxel_size;
    cluster_size(k,7) = mean(KK)*voxel_size;
    
    % Check if the cluster is in contact with the domain boundary
    if k>1 % Exclude largest cluster
        if (x_min==1 || y_min==1 || z_min==1 || x_max==Domain_size(1) || y_max==Domain_size(2) || z_max==Domain_size(3))
            connectivity_id = unknown_id; % The current cluster connectivity is unknown
        else
            connectivity_id = isolated_id; % The current cluster is not connected to the main cluster
        end
        Clusters_LargestIsolatedUnknown.array(I)=connectivity_id;
        cluster_size(k,8)=connectivity_id;
    end
   
    % % Check if cluster connect opposite faces
    % Along direction 1
    if x_min==1 && x_max==Domain_size(1)
        connectivity_id_1 = connected_id; % Connect
    elseif y_min==1 || z_min==1 || y_max==Domain_size(2) || z_max==Domain_size(3)
        connectivity_id_1 = unknown_id; % Does not connect, but connection actually unknown due to limited field of view
    else
        connectivity_id_1 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFace2Face.direction(1).array(I) = connectivity_id_1;
    cluster_size(k,9)=connectivity_id_1;
    
    % Along direction 2
    if y_min==1 && y_max==Domain_size(2)
        connectivity_id_2 = 1; % Connect
    elseif x_min==1 || z_min==1 || x_max==Domain_size(1) || z_max==Domain_size(3)
        connectivity_id_2 = unknown_id; % Does not connect, but connection actually unknown due to limited field of view
    else
        connectivity_id_2 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFace2Face.direction(2).array(I) = connectivity_id_2;
    cluster_size(k,10)=connectivity_id_2;
    
    % Along direction 3
    if z_min==1 && z_max==Domain_size(3)
        connectivity_id_3 = 1; % Does not connect (true isolated cluster)
    elseif y_min==1 || x_min==1 || y_max==Domain_size(2) || x_max==Domain_size(1)
        connectivity_id_3 = unknown_id; % Does not connect (true isolated cluster)
    else
        connectivity_id_3 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFace2Face.direction(3).array(I) = connectivity_id_3;
    cluster_size(k,11)=connectivity_id_3;


    % % Check if cluster connect to one face
    % Along direction 1
    if x_min==1
        connectivity_id_1 = connected_id; % Connect
    elseif y_min==1 || z_min==1 || y_max==Domain_size(2) || z_max==Domain_size(3)
        connectivity_id_1 = unknown_id; % Does not connect, but connection actually unknown due to limited field of view
    else
        connectivity_id_1 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFromFace1.direction(1).array(I) = connectivity_id_1;
    cluster_size(k,12)=connectivity_id_1;
    
    % Along direction 2
    if y_min==1
        connectivity_id_2 = 1; % Connect
    elseif x_min==1 || z_min==1 || x_max==Domain_size(1) || z_max==Domain_size(3)
        connectivity_id_2 = unknown_id; % Does not connect, but connection actually unknown due to limited field of view
    else
        connectivity_id_2 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFromFace1.direction(2).array(I) = connectivity_id_2;
    cluster_size(k,13)=connectivity_id_2;
    
    % Along direction 3
    if z_min==1
        connectivity_id_3 = 1; % Does not connect (true isolated cluster)
    elseif y_min==1 || x_min==1 || y_max==Domain_size(2) || x_max==Domain_size(1)
        connectivity_id_3 = unknown_id; % Does not connect (true isolated cluster)
    else
        connectivity_id_3 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFromFace1.direction(3).array(I) = connectivity_id_3;
    cluster_size(k,14)=connectivity_id_3;

    % % Check if cluster connect to one face
    % Along direction 1
    if x_max==Domain_size(1)
        connectivity_id_1 = connected_id; % Connect
    elseif y_min==1 || z_min==1 || y_max==Domain_size(2) || z_max==Domain_size(3)
        connectivity_id_1 = unknown_id; % Does not connect, but connection actually unknown due to limited field of view
    else
        connectivity_id_1 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFromFace2.direction(1).array(I) = connectivity_id_1;
    cluster_size(k,15)=connectivity_id_1;
    
    % Along direction 2
    if y_max==Domain_size(2)
        connectivity_id_2 = 1; % Connect
    elseif x_min==1 || z_min==1 || x_max==Domain_size(1) || z_max==Domain_size(3)
        connectivity_id_2 = unknown_id; % Does not connect, but connection actually unknown due to limited field of view
    else
        connectivity_id_2 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFromFace2.direction(2).array(I) = connectivity_id_2;
    cluster_size(k,16)=connectivity_id_2;
    
    % Along direction 3
    if z_max==Domain_size(3)
        connectivity_id_3 = 1; % Does not connect (true isolated cluster)
    elseif y_min==1 || x_min==1 || y_max==Domain_size(2) || x_max==Domain_size(1)
        connectivity_id_3 = unknown_id; % Does not connect (true isolated cluster)
    else
        connectivity_id_3 = isolated_id; % Does not connect (true isolated cluster)
    end
    Clusters_TransportFromFace2.direction(3).array(I) = connectivity_id_3;
    cluster_size(k,17)=connectivity_id_3;


end

% Save results in structure
Connectivity_structure.Clusters_sortedpersize = cluster_size;

%% ANALYSE RESULT
% values in percents

% Main, unknown, isolated
Clusters_LargestIsolatedUnknown.main_cluster_phasefraction     = 100 * sum(sum(sum( Clusters_LargestIsolatedUnknown.array==connected_id))) / Numbervoxel_phase; 
Clusters_LargestIsolatedUnknown.unknown_cluster_phasefraction  = 100 * sum(sum(sum( Clusters_LargestIsolatedUnknown.array==unknown_id))) / Numbervoxel_phase; 
Clusters_LargestIsolatedUnknown.isolated_cluster_phasefraction = 100 * sum(sum(sum( Clusters_LargestIsolatedUnknown.array==isolated_id))) / Numbervoxel_phase; 
% Connection along direction
for current_direction=1:1:3
    Clusters_TransportFace2Face.direction(current_direction).connected_cluster_phasefraction = 100 * sum(sum(sum( Clusters_TransportFace2Face.direction(current_direction).array==connected_id))) / Numbervoxel_phase; 
    Clusters_TransportFace2Face.direction(current_direction).unknown_cluster_phasefraction   = 100 * sum(sum(sum( Clusters_TransportFace2Face.direction(current_direction).array==unknown_id))) / Numbervoxel_phase; 
    Clusters_TransportFace2Face.direction(current_direction).isolated_cluster_phasefraction  = 100 * sum(sum(sum( Clusters_TransportFace2Face.direction(current_direction).array==isolated_id))) / Numbervoxel_phase; 

    Clusters_TransportFromFace1.direction(current_direction).connected_cluster_phasefraction = 100 * sum(sum(sum( Clusters_TransportFromFace1.direction(current_direction).array==connected_id))) / Numbervoxel_phase; 
    Clusters_TransportFromFace1.direction(current_direction).unknown_cluster_phasefraction   = 100 * sum(sum(sum( Clusters_TransportFromFace1.direction(current_direction).array==unknown_id))) / Numbervoxel_phase; 
    Clusters_TransportFromFace1.direction(current_direction).isolated_cluster_phasefraction  = 100 * sum(sum(sum( Clusters_TransportFromFace1.direction(current_direction).array==isolated_id))) / Numbervoxel_phase; 

    Clusters_TransportFromFace2.direction(current_direction).connected_cluster_phasefraction = 100 * sum(sum(sum( Clusters_TransportFromFace2.direction(current_direction).array==connected_id))) / Numbervoxel_phase; 
    Clusters_TransportFromFace2.direction(current_direction).unknown_cluster_phasefraction   = 100 * sum(sum(sum( Clusters_TransportFromFace2.direction(current_direction).array==unknown_id))) / Numbervoxel_phase; 
    Clusters_TransportFromFace2.direction(current_direction).isolated_cluster_phasefraction  = 100 * sum(sum(sum( Clusters_TransportFromFace2.direction(current_direction).array==isolated_id))) / Numbervoxel_phase; 

end
% Save results in structure
Connectivity_structure.Clusters_LargestIsolatedUnknown = Clusters_LargestIsolatedUnknown;
Connectivity_structure.Clusters_TransportFace2Face = Clusters_TransportFace2Face;
Connectivity_structure.Clusters_TransportFromFace1 = Clusters_TransportFromFace1;
Connectivity_structure.Clusters_TransportFromFace2 = Clusters_TransportFromFace2;

end