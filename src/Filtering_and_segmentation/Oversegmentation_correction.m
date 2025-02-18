function [Label_lake] = Oversegmentation_correction(binary_phase, Label_lake)

% STEP 1: search for non-continguous discrete particle and correct them
% STEP 2: treat 1 voxel size particle
% STEP 3: calculate D-PSD size
% STEP 4: calculate C-PSD size from another algorithm (will be done only for the first iteration)
% STEP 5: treat voxel whom D-PSD size < C-PSD size

% Initialize iteration
iteration_=0;
number_change=1e9;

% Domain size of the microstructure
Domain_size=size(binary_phase);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end

% Initialise Discrete_particle checksum and list
% Checksum_state0_Label_lake = DataHash(Label_lake);
% That may induce an OUT OF MEMORY error if Label_lake is too large
% Instead, we are calculating the hash slice by slice
Checksum_state0_Label_lake=[];
for current_=1:1:Domain_size(3)
    current_checksum = DataHash(Label_lake(:,:,current_));
    Checksum_state0_Label_lake=[Checksum_state0_Label_lake current_checksum];
end
Checksum_list.state(1).hash=Checksum_state0_Label_lake;
% Checksum (hash) on Label_lake is used to verify that the iterative process does not bring us back to a previous state and run indefinitely
% The function used to check the Hask is DataHash, from Jan Simon, downloaded from Matlab File Exchange (Copyright (c) 2016, Jan Simon, All rights reserved)
% Copyright notice, list of conditions, and disclaimer are in the third party licences folder of the repo

while number_change~=0
    iteration_=iteration_+1;
    fprintf ('Current iteration: %i\n',iteration_);

    %% STEP 1: SEARCH FOR NON-CONTINGUOUS PARTICLE

    % We use the Matlab built-in function for detetecting and labeling non continguous particle/cluster
    % L = bwlabeln(BW,6)
    % L is containing labels for the connected components in BW
    % L is the same size as BW
    % The elements of L are integer values greater than or equal to 0.
    % The pixels labeled 0 are the background. The pixels labeled 1 make up one connected cluster. The pixels labeled 2 make up a second connected cluster, and so on

    % Get all unique id particle
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove id 0 (complementary phase)
    % Get number of different id
    number_id = length(unique_particle);
    % Loop over id
    for current_=1:1:number_id
        % Get id
        current_id = unique_particle(current_);
        % Create binary matrix for this particle id
        binary_particle = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
        binary_particle(Label_lake==current_id)=1;
        % Label particle
        L_particle = bwlabeln(binary_particle,6);
        % Unique label
        unique_cluster = unique(L_particle);
        unique_cluster(1)=[]; % remobe background cluster id
        % Number of distinct cluster
        cluster_number=length(unique_cluster);
        % Treat the case where there are more than 1 cluster
        if cluster_number>1
            % The largest cluster will be unchanged
            % All the other cluster will be treated as follow:
            % Case a) If they are surrounded by ONLY voxels of the complementary phase, then this cluster is considered as a new particle
            % Case b) otherwise, voxels of this cluster will be attributed to the particle id of the nearest voxel that belong to the phase and
            % that does not belong to this particular cluster
            cluster_size=zeros(cluster_number,2);
            % Column: cluster id, size
            for current_cluster=1:1:cluster_number
                cluster_id = unique_cluster(current_cluster);
                cluster_size(current_cluster,1)=cluster_id;
                cluster_size(current_cluster,2)=sum(sum(sum(L_particle==cluster_id)));
            end
            % Sort per cluster size (decreasing order)
            cluster_size = sortrows(cluster_size,-2);
            % We will deal only with the second, third etc. cluster
            for analysed_cluster=2:1:cluster_number
                % Get cluster id
                cluster_id = cluster_size(analysed_cluster,1);
                % Create cluster binary matrix
                binary_cluster = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
                binary_cluster(L_particle==cluster_id)=1;
                % Find all voxels of this cluster
                index_cluster=find(binary_cluster==1);
                % Get all their coordinate
                [Cl1,Cl2,Cl3] = ind2sub(Domain_size,index_cluster);
                % Number of voxel
                number_voxel_cluster = length(Cl1);
                % % Identify voxel of the complemtary that belong at the border
                % They will be marked 2 in the binary phase
                surrounding_values=[];
                pos_surrounding=zeros(1,5);
                for current_voxel=1:1:number_voxel_cluster
                    % Get back voxel coordinate
                    x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                    % Check x-
                    if x_>1
                        if Label_lake(x_-1,y_,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_-1,y_,z_)];
                            if Label_lake(x_-1,y_,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_-1,y_,z_,0,Label_lake(x_-1,y_,z_)]];
                            end
                        end
                    end
                    % Check y-
                    if y_>1
                        if Label_lake(x_,y_-1,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_-1,z_)];
                            if Label_lake(x_,y_-1,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_-1,z_,0,Label_lake(x_,y_-1,z_)]];
                            end
                        end
                    end
                    % Check z-
                    if z_>1
                        if Label_lake(x_,y_,z_-1)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_,z_-1)];
                            if Label_lake(x_,y_,z_-1)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_,z_-1,0,Label_lake(x_,y_,z_-1)]];
                            end
                        end
                    end
                    % Check x+
                    if x_<Domain_size(1)
                        if Label_lake(x_+1,y_,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_+1,y_,z_)];
                            if Label_lake(x_+1,y_,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_+1,y_,z_,0,Label_lake(x_+1,y_,z_)]];
                            end
                        end
                    end
                    % Check y+
                    if y_<Domain_size(2)
                        if Label_lake(x_,y_+1,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_+1,z_)];
                            if Label_lake(x_,y_+1,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_+1,z_,0,Label_lake(x_,y_+1,z_)]];
                            end
                        end
                    end
                    % Check z+
                    if z_<Domain_size(3)
                        if Label_lake(x_,y_,z_+1)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_,z_+1)];
                            if Label_lake(x_,y_,z_+1)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_,z_+1,0,Label_lake(x_,y_,z_+1)]];
                            end
                        end
                    end
                end
                % Remove fist 0 line
                pos_surrounding(1,:)=[];
                [n_surrounding,~]=size(pos_surrounding);
                % Replace accoring to case a or case b
                if length(unique(surrounding_values))==1 && surrounding_values(1)==0
                    % Case a:
                    unique_particle = unique(Label_lake);
                    max_id = max(unique_particle);
                    new_id = max_id+1;
                    Label_lake(index_cluster)=new_id;
                else
                    % Case b:
                    % Calculate minimum distance to nearest voxel
                    for current_voxel=1:1:number_voxel_cluster
                        % Get back voxel coordinate
                        x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                        % Calculate all distance betwenen the voxel and the surrounding voxels
                        for voxel_surround=1:1:n_surrounding
                            x_s=double(pos_surrounding(voxel_surround,1));
                            y_s=double(pos_surrounding(voxel_surround,2));
                            z_s=double(pos_surrounding(voxel_surround,3));
                            dist_surround = sqrt((x_s-x_)^2 + (y_s-y_)^2 + (z_s-z_)^2);
                            pos_surrounding(voxel_surround,4)=dist_surround;
                        end
                        % Sort by increasing order
                        pos_surrounding = sortrows(pos_surrounding,4);
                        % Get minimum distance
                        min_dist = pos_surrounding(1,4);
                        % Get all index with this distance
                        index_min = find(pos_surrounding(:,4)==min_dist);
                        % Get all id
                        list_id=[];
                        for kkk=1:1:length(index_min)
                            list_id = [list_id pos_surrounding(kkk,5)];
                        end
                        % Get unique id
                        unique_surrounding = unique(list_id);
                        % Get the value that is preominemt
                        sum_surrounding=zeros(length(unique_surrounding),1);
                        for n_=1:1:length(unique_surrounding)
                            sum_surrounding(n_)=sum(list_id==unique_surrounding(n_));
                        end
                        max_index=find(sum_surrounding==max(sum_surrounding));
                        max_index=max_index(1); % In case there is equality
                        preominent_value = unique_surrounding(max_index);
                        % Replace value
                        Label_lake(x_,y_,z_)=preominent_value;
                    end
                end

            end
        end
    end
    % Get all discrete particle id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
    % Get number of particle
    number_particle = length(unique_particle);
    fprintf ('   - Search for non connected particle: number of discrete particle: %i\n',number_particle);

    %% STEP 2: CASE OF THE X VOXEL SIZE PARTICLE

    % X voxel-size particle will be removed
    % Since I do NOT want parameter, i will consider only the case X=1
    critical_minimal_size = 1;

    % If you really want to use it, you should consider iterate from smaller particle to
    % larger particle (modify the algorithm in consequence)

    % Get all unique id particle
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove id 0 (complementary phase)
    % Get number of different id
    number_id = length(unique_particle);
    % Loop over id
    for current_=1:1:number_id
        % Get id
        current_id = unique_particle(current_);
        % Size
        D_PSD_particle_size = sum(sum(sum(Label_lake==current_id)));
        if D_PSD_particle_size<=critical_minimal_size
            index_particle=find(Label_lake==current_id);
            % Get all their coordinate
            [Cl1,Cl2,Cl3] = ind2sub(Domain_size,index_particle);
            % Number of voxel
            number_voxel_particle = length(Cl1);
            % % Identify voxel of the complemtary that belong at the border
            surrounding_values=[];
            pos_surrounding=zeros(1,5);
            % Loop over all voxel of the particle
            for current_voxel=1:1:number_voxel_particle
                % Get back voxel coordinate
                x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                % Check x-
                if x_>1
                    if Label_lake(x_-1,y_,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_-1,y_,z_)];
                        if Label_lake(x_-1,y_,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_-1,y_,z_,0,Label_lake(x_-1,y_,z_)]];
                        end
                    end
                end
                % Check y-
                if y_>1
                    if Label_lake(x_,y_-1,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_-1,z_)];
                        if Label_lake(x_,y_-1,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_-1,z_,0,Label_lake(x_,y_-1,z_)]];
                        end
                    end
                end
                % Check z-
                if z_>1
                    if Label_lake(x_,y_,z_-1)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_,z_-1)];
                        if Label_lake(x_,y_,z_-1)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_,z_-1,0,Label_lake(x_,y_,z_-1)]];
                        end
                    end
                end
                % Check x+
                if x_<Domain_size(1)
                    if Label_lake(x_+1,y_,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_+1,y_,z_)];
                        if Label_lake(x_+1,y_,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_+1,y_,z_,0,Label_lake(x_+1,y_,z_)]];
                        end
                    end
                end
                % Check y+
                if y_<Domain_size(2)
                    if Label_lake(x_,y_+1,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_+1,z_)];
                        if Label_lake(x_,y_+1,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_+1,z_,0,Label_lake(x_,y_+1,z_)]];
                        end
                    end
                end
                % Check z+
                if z_<Domain_size(3)
                    if Label_lake(x_,y_,z_+1)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_,z_+1)];
                        if Label_lake(x_,y_,z_+1)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_,z_+1,0,Label_lake(x_,y_,z_+1)]];
                        end
                    end
                end
            end

            % Remove fist 0 line
            pos_surrounding(1,:)=[];
            [n_surrounding,~]=size(pos_surrounding);

            % Replace accoring to case a or case b
            if length(unique(surrounding_values))==1 && surrounding_values(1)==0
                % Case a: is removed
                Label_lake(index_particle)=0;
            else
                % Case b:
                % Calculate minimum distance to nearest voxel
                for current_voxel=1:1:number_voxel_particle
                    % Get back voxel coordinate
                    x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                    % Calculate all distance betwenen the voxel and the surrounding voxels
                    for voxel_surround=1:1:n_surrounding
                        x_s=double(pos_surrounding(voxel_surround,1));
                        y_s=double(pos_surrounding(voxel_surround,2));
                        z_s=double(pos_surrounding(voxel_surround,3));
                        dist_surround = sqrt((x_s-x_)^2 + (y_s-y_)^2 + (z_s-z_)^2);
                        pos_surrounding(voxel_surround,4)=dist_surround;
                    end
                    % Sort by increasing order
                    pos_surrounding = sortrows(pos_surrounding,4);
                    % Get minimum distance
                    min_dist = pos_surrounding(1,4);
                    % Get all index with this distance
                    index_min = find(pos_surrounding(:,4)==min_dist);
                    % Get all id
                    list_id=[];
                    for kkk=1:1:length(index_min)
                        list_id = [list_id pos_surrounding(kkk,5)];
                    end
                    % Get unique id
                    unique_surrounding = unique(list_id);
                    % Get the value that is preominemt
                    sum_surrounding=zeros(length(unique_surrounding),1);
                    for n_=1:1:length(unique_surrounding)
                        sum_surrounding(n_)=sum(list_id==unique_surrounding(n_));
                    end
                    max_index=find(sum_surrounding==max(sum_surrounding));
                    max_index=max_index(1); % In case there is equality
                    preominent_value = unique_surrounding(max_index);
                    % Replace value
                    Label_lake(x_,y_,z_)=preominent_value;
                end
            end

        end
    end
    % Re-order seed id
    unique_id = unique(Label_lake);
    number_id = length(unique_id);
    for current_id=2:1:number_id % 0 is the complementary phase
        id_to_be_replaced = unique_id(current_id);
        index_ = find(Label_lake==id_to_be_replaced);
        Label_lake(index_)=current_id-1;
    end
    % Get all discrete particle id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
    % Get number of particle
    number_particle = length(unique_particle);
    fprintf ('   - One-voxel size particle treated:   number of discrete particle: %i\n',number_particle);

    %% STEP 3: Discrete-PSD
    % Initialisation
    D_PSD_particle_size = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    % Unique id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[];
    % Number of different particle
    number_particle = length(unique_particle);
    % loop over all particles
    for current_=1:1:number_particle
        % Get id_
        particle_id = unique_particle(current_);
        % Find all voxels of the particle
        index_particle=find(Label_lake==particle_id);
        % Number of voxel
        number_voxel_particle = length(index_particle);
        volume_ = number_voxel_particle;
        area_ = number_voxel_particle;
        % Equivalent diameter size
        if Domain_size(3)>1
            % 3D case
            equivalent_diameter_size= 2 * ((3*volume_  /(4*pi))^(1/3));
        else
            % 2D case
            equivalent_diameter_size= 2 * ((area_/pi)^(1/2));
        end
        D_PSD_particle_size(index_particle)= equivalent_diameter_size;
    end

    %% STEP 4: Continuum-PSD
    % It will be calculated only one time
    if iteration_==1
        % Run the C-PSD algorithm
        [C_PSD_particle_size] = Function_particle_size_CPSD_Algorithm(binary_phase);
    end

    %% STEP 5: COMPARE SIZE BETWEEN C-PSD AND D-PSD
    % Initialise
    smaller_than_cpsd = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    % Set =1 when D-PSD is smaller than C-PSD
    smaller_than_cpsd(D_PSD_particle_size<C_PSD_particle_size)=1;
    % % Two cases can explaines D-PSD < C-PSD
    % Case a) The particle is at the domain's edge: C-PSD algorithm does not exibit
    % edge effect compare to the D-PSD that will decrease the particle size at
    % the border. Then, there is nothing to do for these particles
    % Case b) these voxels are likely stuck between two large particles. We will replace them

    % L = bwlabeln(BW,6)
    % L is containing labels for the connected components in BW
    % L is the same size as BW
    % The elements of L are integer values greater than or equal to 0.
    % The pixels labeled 0 are the background. The pixels labeled 1 make up one connected cluster. The pixels labeled 2 make up a second connected cluster, and so on
    % Label zone
    L_zone = bwlabeln(smaller_than_cpsd,6);

    % % Step1: remove from the analysis zone located at the domain's edge
    % Unique zone
    unique_zone = unique(Label_lake);
    unique_zone(1)=[]; % remove background
    % Number of distinct zone
    zone_number=length(unique_zone);
    % Loop over all zones
    for current_=1:1:zone_number
        % Get id
        current_zone = unique_zone(current_);
        % Get all voxels
        index_zone = find(Label_lake==current_zone);
        % Get all coordinates
        [Z1,Z2,Z3] = ind2sub(Domain_size,index_zone);
        % Get min max
        x_min=min(Z1); y_min=min(Z2); z_min=min(Z3);
        x_max=max(Z1); y_max=max(Z2); z_max=max(Z3);
        % Check if the zone is located at the border
        if Domain_size(3)>1
            % 3D case
            at_the_border = (x_min==1 || y_min==1 || z_min==1 || x_max==Domain_size(1) || y_max==Domain_size(2) || z_max==Domain_size(3));
        else
            % 2D case
            at_the_border = (x_min==1 || y_min==1 || x_max==Domain_size(1) || y_max==Domain_size(2));
        end
        % Remove the zone from the analysis if true
        if at_the_border==1
            L_zone(index_zone)=0;
            smaller_than_cpsd(index_zone)=0;
        end
    end
    
    % % Step 2: replace the remaining voxel with nearest voxel
    % Unique zone
    unique_zone = unique(L_zone);
    unique_zone(1)=[]; % remove background
    % Number of distinct zone
    zone_number=length(unique_zone);
    % Loop over all zones
    for current_=1:1:zone_number
        % Get id
        current_zone = unique_zone(current_);
        % Get all voxels
        index_zone = find(L_zone==current_zone);
        % Get all coordinates
        [Z1,Z2,Z3] = ind2sub(Domain_size,index_zone);
        % Number of voxel
        number_voxel_zone  = length(Z1);
        % % Identify voxel of the complemtary that belong at the border
        pos_surrounding=zeros(1,5);
        for current_voxel=1:1:number_voxel_zone
            % Get back voxel coordinate
            x_=Z1(current_voxel); y_=Z2(current_voxel); z_=Z3(current_voxel);
            % Check x-
            if x_>1
                if L_zone(x_-1,y_,z_)~=current_zone
                    if Label_lake(x_-1,y_,z_)~=0
                        pos_surrounding=[pos_surrounding; [x_-1,y_,z_,0,Label_lake(x_-1,y_,z_)]];
                    end
                end
            end
            % Check y-
            if y_>1
                if L_zone(x_,y_-1,z_)~=current_zone
                    if Label_lake(x_,y_-1,z_)~=0
                        pos_surrounding=[pos_surrounding; [x_,y_-1,z_,0,Label_lake(x_,y_-1,z_)]];
                    end
                end
            end
            % Check z-
            if z_>1
                if L_zone(x_,y_,z_-1)~=current_zone
                    if Label_lake(x_,y_,z_-1)~=0
                        pos_surrounding=[pos_surrounding; [x_,y_,z_-1,0,Label_lake(x_,y_,z_-1)]];
                    end
                end
            end
            % Check x+
            if x_<Domain_size(1)
                if L_zone(x_+1,y_,z_)~=current_zone
                    if Label_lake(x_+1,y_,z_)~=0
                        pos_surrounding=[pos_surrounding; [x_+1,y_,z_,0,Label_lake(x_+1,y_,z_)]];
                    end
                end
            end
            % Check y+
            if y_<Domain_size(2)
                if L_zone(x_,y_+1,z_)~=current_zone
                    if Label_lake(x_,y_+1,z_)~=0
                        pos_surrounding=[pos_surrounding; [x_,y_+1,z_,0,Label_lake(x_,y_+1,z_)]];
                    end
                end
            end
            % Check z+
            if z_<Domain_size(3)
                if L_zone(x_,y_,z_+1)~=current_zone
                    if Label_lake(x_,y_,z_+1)~=0
                        pos_surrounding=[pos_surrounding; [x_,y_,z_+1,0,Label_lake(x_,y_,z_+1)]];
                    end
                end
            end
        end
        % Remove fist 0 line
        pos_surrounding(1,:)=[];
        [n_surrounding,~]=size(pos_surrounding);

        if n_surrounding>0
            % Calculate minimum distance to nearest voxel
            for current_voxel=1:1:number_voxel_zone
                % Get back voxel coordinate
                x_=Z1(current_voxel); y_=Z2(current_voxel); z_=Z3(current_voxel);
                % Calculate all distance betwenen the voxel and the surrounding voxels
                for voxel_surround=1:1:n_surrounding
                    x_s=double(pos_surrounding(voxel_surround,1));
                    y_s=double(pos_surrounding(voxel_surround,2));
                    z_s=double(pos_surrounding(voxel_surround,3));
                    dist_surround = sqrt((x_s-x_)^2 + (y_s-y_)^2 + (z_s-z_)^2);
                    pos_surrounding(voxel_surround,4)=dist_surround;
                end
                % Sort by increasing order
                pos_surrounding = sortrows(pos_surrounding,4);
                % Get minimum distance
                min_dist = pos_surrounding(1,4);
                % Get all index with this distance
                index_min = find(pos_surrounding(:,4)==min_dist);
                % Get all id
                list_id=[];
                for kkk=1:1:length(index_min)
                    list_id = [list_id pos_surrounding(kkk,5)];
                end
                % Get unique id
                unique_surrounding = unique(list_id);
                % Get the value that is preominemt
                sum_surrounding=zeros(length(unique_surrounding),1);
                for n_=1:1:length(unique_surrounding)
                    sum_surrounding(n_)=sum(list_id==unique_surrounding(n_));
                end
                max_index=find(sum_surrounding==max(sum_surrounding));
                max_index=max_index(1); % In case there is equality
                preominent_value = unique_surrounding(max_index);
                % Replace value
                Label_lake(x_,y_,z_)=preominent_value;
            end
        end
    end
    % Re-order seed id
    unique_id = unique(Label_lake);
    number_id = length(unique_id);
    for current_id=2:1:number_id % 0 is the complementary phase
        id_to_be_replaced = unique_id(current_id);
        index_ = find(Label_lake==id_to_be_replaced);
        Label_lake(index_)=current_id-1;
    end
    % Get all discrete particle id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
    % Get number of particle
    number_particle = length(unique_particle);
    fprintf ('   - C-PSD & D-PSD analysis:            number of discrete particle: %i\n',number_particle);
   

    %% CHANGE FROM PREVIOUS ITERATION: EXIT WHILE LOOP TEST
    % Compare with previous iteration resutt
    if iteration_>1
        number_change=sum(sum(sum(previous_discrete_particle~=Label_lake)));
        fprintf ('   CONVERGENCE CHECK: Number of voxels changed compared with previous iteration: %i\n',number_change);
    end
    % For the next iteration
    previous_discrete_particle = Label_lake;

    %% CHECKSUM

    % Current Discrete_particle checksum
    % Checksum_currentstate_Label_lake = DataHash(Label_lake);
    % That may induce an OUT OF MEMORY error if Label_lake is too large
    % Instead, we are calculating the hash slice by slice
    Checksum_currentstate_Label_lake=[];
    for current_=1:1:Domain_size(3)
        current_checksum = DataHash(Label_lake(:,:,current_));
        Checksum_currentstate_Label_lake=[Checksum_currentstate_Label_lake current_checksum];
    end

    % Compare with previous checksum
    exit_loop_checksum=0;
    for previous_state=1:1:iteration_
        Checksum_previousstate_Label_lake = Checksum_list.state(previous_state).hash;
        check_checksum = isequal(Checksum_previousstate_Label_lake,Checksum_currentstate_Label_lake);
        if check_checksum==1
            exit_loop_checksum=1;
        end
    end
    % Save current checksum in the list
    Checksum_list.state(iteration_+1).hash=Checksum_currentstate_Label_lake;
    % Exit loop due to checksum equality
    if (exit_loop_checksum==1 && number_change>0)
        disp 'Iterative process has been founded to back to a previous state. Iterations are stopped.'
        break
    end

end

end

