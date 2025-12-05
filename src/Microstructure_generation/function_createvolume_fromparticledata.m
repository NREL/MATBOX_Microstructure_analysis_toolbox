function [PhaseLabel, ParticleId, overlapping_stat] = function_createvolume_fromparticledata(particle_data,domain_size,cropparameters,scale_diameterratio,force_separation,distance_separation,porosity_target, label_original_vf, shuffleparticles, keep_relativesolidvf, keep_solidvf)

[n_label,~] = size(label_original_vf);
original_solidvolumefraction = label_original_vf(:,2);

PhaseLabel = zeros(domain_size); % Initialize volume
ParticleId = zeros(domain_size); % Initialize volume

[number_particle,~] = size(particle_data);
overlappingids = [];

if porosity_target>0 && shuffleparticles
    particle_data = particle_data(randperm(size(particle_data, 1)), :); % Shuffle the particle order
end

relative_solidvolumefraction = original_solidvolumefraction/sum(original_solidvolumefraction);
label_upscaled_vf = label_original_vf;
label_upscaled_vf(:,2) = (1-porosity_target)*relative_solidvolumefraction;

current_label_vf = zeros(n_label,1);

number_voxel = numel(PhaseLabel);
number_voxelpore = number_voxel;
for k_particle = 1:1:number_particle % Loop over all particles
    if porosity_target>0
        current_porosity = number_voxelpore / number_voxel;
        %current_porosity = sum(sum(sum(PhaseLabel==0)))/number_voxel;
        if current_porosity<=porosity_target
            break
        end
    end

    % Particle metric
    id = particle_data(k_particle,1);
    label = particle_data(k_particle,2);
    x_center = particle_data(k_particle,3);
    y_center = particle_data(k_particle,4);
    z_center = particle_data(k_particle,5);
    dx = round(particle_data(k_particle,6)*scale_diameterratio);
    dy = round(particle_data(k_particle,7)*scale_diameterratio);
    dz = round(particle_data(k_particle,8)*scale_diameterratio);
    angle_x_deg = particle_data(k_particle,9);
    angle_y_deg = particle_data(k_particle,10);
    angle_z_deg = particle_data(k_particle,11);

    if x_center<1 || y_center<1 || z_center<1 || x_center>domain_size(1) || y_center>domain_size(2) || z_center>domain_size(3)
        continue
    end

    idx_label = 1;
    if porosity_target>0
        if keep_relativesolidvf
            idx_label = find(label_upscaled_vf(:,1)==label);
            if current_label_vf(idx_label)/number_voxel >= label_upscaled_vf(idx_label,2)
                continue
            end
        elseif keep_solidvf
            idx_label = find(label_original_vf(:,1)==label);
            if current_label_vf(idx_label)/number_voxel >= label_original_vf(idx_label,2)
                continue
            end
        end
    end

    % Create particle
    binary_ellipsoid = create_ellipsoid(dx,dy,dz); % Create ellipsoid
    if (angle_x_deg~=0 && angle_x_deg~=180) || (angle_y_deg~=0 && angle_y_deg~=180) || (angle_z_deg~=0 && angle_z_deg~=180)  % Only if one or more rotation
        binary_ellipsoid = rotate_domain(binary_ellipsoid,angle_x_deg, angle_y_deg, angle_z_deg); % Apply rotation matrix
        binary_ellipsoid(binary_ellipsoid>=0.5)=1;
        binary_ellipsoid(binary_ellipsoid~=1)=0;
    end
    size_ellipsoid = size(binary_ellipsoid); % Dimension of the ellipsoid

    % Subdomain and domain bounds that contain the ellipsoid
    [x_sub_min, x_sub_max, x_min, x_max] = bounds_subdomain_domain(size_ellipsoid(1),x_center,1,domain_size(1));
    [y_sub_min, y_sub_max, y_min, y_max] = bounds_subdomain_domain(size_ellipsoid(2),y_center,1,domain_size(2));
    [z_sub_min, z_sub_max, z_min, z_max] = bounds_subdomain_domain(size_ellipsoid(3),z_center,1,domain_size(3));
    subdomain_id = ParticleId(x_min:x_max, y_min:y_max, z_min:z_max); % Subdomain
    subdomain_label = PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max);
    subdomain_ellipsoid = binary_ellipsoid(x_sub_min:x_sub_max, y_sub_min:y_sub_max, z_sub_min:z_sub_max); % Ellipsoid subdomain

    % Ellipsoids insertion
    unique_subdomain = unique(subdomain_id);
    if length(unique_subdomain)==1 && unique_subdomain==0 % Subdomain is empty
        PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_ellipsoid * label; % Insert ellipsoid within the domain.
        ParticleId(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_ellipsoid * id;
        if porosity_target>0
            volparticle = sum(sum(sum(subdomain_ellipsoid)));
            number_voxelpore = number_voxelpore - volparticle;
            current_label_vf(idx_label) = current_label_vf(idx_label) + volparticle; 
        end
    else % Subdomain contains other particles
        subdomain_binary = subdomain_label;
        subdomain_binary(subdomain_binary~=0)=1;
        union_subdomain_ellipsoid = subdomain_binary + subdomain_ellipsoid;
        idx_overlapping = find(union_subdomain_ellipsoid==2);
        if isempty(idx_overlapping) % There is no overlapping
            PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_label + subdomain_ellipsoid * label; % Insert ellipsoid within the domain, in addition to existing particles
            ParticleId(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_id + subdomain_ellipsoid * id;
            if porosity_target>0
                volparticle = sum(sum(sum(subdomain_ellipsoid)));
                number_voxelpore = number_voxelpore - volparticle;
                current_label_vf(idx_label) = current_label_vf(idx_label) + volparticle; 
            end
        else % Overlapping, we need to distribute conflicting voxels based on nearest particle
            subdomain_union_minus_intersection = abs(subdomain_ellipsoid - subdomain_binary); % (particle union others) - (particle intersection others)
            if sum(sum(sum(subdomain_union_minus_intersection))) > 0
                [Euclidean_distance_map,Idx_distance_map] = bwdist(subdomain_union_minus_intersection);
                index_nearest_particle = Idx_distance_map(idx_overlapping);

                Sub_label = subdomain_label + subdomain_ellipsoid * label;
                Sub_Id = subdomain_id + subdomain_ellipsoid * id;
                Sub_label(idx_overlapping) = Sub_label(index_nearest_particle);
                Sub_Id(idx_overlapping) = Sub_Id(index_nearest_particle);

                % Pverlapping stat
                [idx_interpenetration_domain] = Index_from_subdomain_to_domain(union_subdomain_ellipsoid, idx_overlapping, domain_size, [x_min y_min z_min]);
                overlappingids = [overlappingids; idx_interpenetration_domain];

                if force_separation
                    bw = zeros(size(Sub_Id));
                    bw(Sub_Id==id)=1;
                    dmap = bwdist(bw,'Chessboard');
                    cond1 = dmap<=distance_separation;
                    cond2 = dmap>0;
                    id_cond = find(cond1 .* cond2);
                    if ~isempty(id_cond)
                        Sub_label(id_cond)=0;
                        Sub_Id(id_cond)=0;
                    end

                end

                if porosity_target>0
                    npore_before = sum(sum(sum( PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max)==0 )));
                    npore_after = sum(sum(sum( Sub_label==0 )));
                    number_voxelpore = number_voxelpore - (npore_before-npore_after);
                    current_label_vf(idx_label) = current_label_vf(idx_label) + (npore_before-npore_after);
                end

                PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max) = Sub_label; % Insert ellipsoid within the domain, in addition to existing particles
                ParticleId(x_min:x_max, y_min:y_max, z_min:z_max) = Sub_Id;
            end

        end
    end
end

%% OVERLAPPING STAT AND CROPPING
overlappingids = unique(overlappingids); % Remove duplicate
cs = [cropparameters(1,1) cropparameters(2,1) cropparameters(3,1)];
sz = size(PhaseLabel);
if sum(cropparameters(:,1))>0
    [IX,IY,IZ]=ind2sub(sz,overlappingids);
    ids_notcropped = ones(length(IX),1);
    ids_notcropped(IX < cs(1)+1) = 0;
    ids_notcropped(IY < cs(2)+1) = 0;
    ids_notcropped(IZ < cs(3)+1) = 0;
    ids_notcropped(IX > sz(1) - cs(1))=0;
    ids_notcropped(IY > sz(2) - cs(2))=0;
    ids_notcropped(IZ > sz(3) - cs(3))=0;
    PhaseLabel = PhaseLabel(cs(1)+1:end-cs(1), cs(2)+1:end-cs(2), cs(3)+1:end-cs(3));
    ParticleId =ParticleId(cs(1)+1:end-cs(1), cs(2)+1:end-cs(2), cs(3)+1:end-cs(3));
    sz = size(PhaseLabel);
    overlapping_stat = sum(ids_notcropped) / sum(sum(sum( PhaseLabel~=0 ))); % Ratio of overlapping over solid phases volume
else
    overlapping_stat = length(overlappingids) / sum(sum(sum( PhaseLabel~=0 ))); % Ratio of overlapping over solid phases volume
end

%% REMOVE TRUNCATED PARTICLES
remove_ids = [];
for dir=1:1:3
    if sum(cropparameters(dir,2:4))>0
        if dir==1
            Idedges_all = [ParticleId(1,:,:);ParticleId(end,:,:)];
        elseif dir==2
            Idedges_all = [ParticleId(:,1,:);ParticleId(:,end,:)];
        else
            Idedges_all = [ParticleId(:,:,1);ParticleId(:,:,end)];
        end
        Idedges_all = reshape(Idedges_all,[numel(Idedges_all),1]);
        Idedges_all(Idedges_all==0)=[];
        Idedges = unique(Idedges_all);
        n_Idedge = length(Idedges);
        if ~isempty(Idedges)
            if cropparameters(dir,2) % Remove particles with more than 4 voxels at the edges
                for kid=1:1:n_Idedge
                    current_particle_data = particle_data(Idedges(kid),:);
                    particle_center = current_particle_data(2+dir);
                    new_particle_center = particle_center-cs(dir);                    
                    if new_particle_center>=1 && new_particle_center<=sz(dir)
                        n_voxel = sum(Idedges_all==Idedges(kid));
                        if n_voxel >4
                            remove_ids = [remove_ids; Idedges(kid)];
                        end
                    else
                        remove_ids = [remove_ids; Idedges(kid)];
                    end
                end
            end
            if cropparameters(dir,3) % Remove particles with more than x% truncated volume
                for kid=1:1:n_Idedge
                    n_voxel_FOV = sum(sum(sum(ParticleId==Idedges(kid) )));
                    current_particle_data = particle_data(Idedges(kid),:);
                    [binary_ellipsoid] = create_ellipsoid(current_particle_data(6),current_particle_data(7),current_particle_data(8));
                    n_voxel_particle = sum(sum(sum( binary_ellipsoid )));
                    truncated_percentage = max([0, 100*(n_voxel_particle-n_voxel_FOV)/n_voxel_particle]);
                    if truncated_percentage>=cropparameters(dir,3)
                        remove_ids = [remove_ids; Idedges(kid)];
                    end
                end
            end
            if cropparameters(dir,4) % Remove particles with center outside FOV (can happen if domain has been cropped)
                for kid=1:1:n_Idedge
                    current_particle_data = particle_data(Idedges(kid),:);
                    % Id / phase label / x / y / z / dx / dy / dz / Rx / Ry / Rz
                    particle_center = current_particle_data(2+dir);
                    new_particle_center = particle_center-cs(dir);
                    if new_particle_center<1 ||  new_particle_center>sz(dir)
                        remove_ids = [remove_ids; Idedges(kid)];
                    end
                end
            end
        end
    end
end
if ~isempty(remove_ids)
    remove_ids = unique(remove_ids); % Remove duplicates
    for kid = 1:1:length(remove_ids)
        remove_idx = find(ParticleId == remove_ids(kid) );
        ParticleId(remove_idx) = 0;
        PhaseLabel(remove_idx) = 0;
    end
end

%% FUNCTIONS
    function [binary_ellipsoid] = create_ellipsoid(dx,dy,dz)
        if ~mod(dx,2) && ~mod(dy,2) && ~mod(dz,2) && dx==dy && dx==dz
            % Sphere, even diameter
            binary_ellipsoid = zeros(dx,dx,dx);
            center1 = dx/2;
            center2 = center1+1;
            binary_ellipsoid(center1,center1,center1) = 1;
            binary_ellipsoid(center2,center1,center1) = 1;
            binary_ellipsoid(center1,center2,center1) = 1;
            binary_ellipsoid(center2,center2,center1) = 1;
            binary_ellipsoid(center1,center1,center2) = 1;
            binary_ellipsoid(center2,center1,center2) = 1;
            binary_ellipsoid(center1,center2,center2) = 1;
            binary_ellipsoid(center2,center2,center2) = 1;
            dmap_sphere = bwdist(binary_ellipsoid);
            radius = (dx-2)/2;
            binary_ellipsoid(dmap_sphere<=radius) = 1;

        else
            % Exact for sphere and ellipsoid with odd diameters. Approximation if one diameter is even
            % Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix
            Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
            [X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
            R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
            binary_ellipsoid = zeros(size(X));
            %binary_ellipsoid(R <= 1 ) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
            idr = find(R<=1);
            if isempty(idr)
                idr = find(R==min(min(min(R))));
            end
            binary_ellipsoid(idr) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
        end

    end

    function [binary_ellipsoid] = rotate_domain(binary_ellipsoid,angle_x_deg, angle_y_deg, angle_z_deg)
        % Angle rotation in radians
        angle_x_rad=deg2rad(angle_x_deg);
        angle_y_rad=deg2rad(angle_y_deg);
        angle_z_rad=deg2rad(angle_z_deg);

        % Rotation matrix
        % https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
        % Transformation matrix that rotates the image around the x-axis
        tx = [1 0 0 0
            0 cos(angle_x_rad)  sin(angle_x_rad) 0
            0 -sin(angle_x_rad) cos(angle_x_rad) 0
            0             0              0       1];
        % Transformation matrix that rotates the image around the y-axis
        ty = [cos(angle_y_rad)  0      -sin(angle_y_rad)   0
            0             1              0     0
            sin(angle_y_rad)    0       cos(angle_y_rad)   0
            0             0              0     1];
        % Transformation matrix that rotates the image around the z-axis
        tz = [cos(angle_z_rad) sin(angle_z_rad) 0 0
            -sin(angle_z_rad)  cos(angle_z_rad) 0 0
            0 0 1 0
            0 0 0 1];
        % Then pass the matrix to the affine3d object constructor.
        tx_form = affine3d(tx);
        ty_form = affine3d(ty);
        tz_form = affine3d(tz);

        % Apply the transformation to the image
        if angle_x_deg~=0
            binary_ellipsoid = imwarp(binary_ellipsoid,tx_form);
        end
        if angle_y_deg~=0
            binary_ellipsoid = imwarp(binary_ellipsoid,ty_form);
        end
        if angle_z_deg~=0
            binary_ellipsoid = imwarp(binary_ellipsoid,tz_form);
        end

        % Crop
        idx = find(binary_ellipsoid~=0);
        [I,J,K]=ind2sub(size(binary_ellipsoid),idx);
        binary_ellipsoid = binary_ellipsoid(min(I):max(I),min(J):max(J),min(K):max(K));
    end

    function [sub_min, sub_max, min_, max_] = bounds_subdomain_domain(dim,center,boundmin,boundmax)
        is_odd = mod(dim,2);
        if is_odd % Odd case
            min_ = center - floor(dim/2);
            max_ = center + floor(dim/2);
        else % Even case
            min_ = center - (dim/2 - 1);
            max_ = center + dim/2;
        end
        if min_<boundmin
            sub_min = boundmin-min_+1;
            min_ = boundmin;
        else
            sub_min = 1;
        end
        if max_>boundmax
            sub_max = dim - (max_-boundmax);
            max_ = boundmax;
        else
            sub_max = dim;
        end
    end

    function [linear_indices_domain] = Index_from_subdomain_to_domain(binary_subdomain, idx, domain_size, subdomain_min_location)
        if isempty(idx)
            idx = find(binary_subdomain==1); % Index of all voxels within the subdomain
        end
        [Isub,Jsub, Ksub] = ind2sub(size(binary_subdomain),idx); % Coordinates of all voxels within the subdomain
        Idomain = Isub + subdomain_min_location(1)-1; % Coordinates of all voxels within the full domain
        Jdomain = Jsub + subdomain_min_location(2)-1;
        Kdomain = Ksub + subdomain_min_location(3)-1;
        linear_indices_domain = sub2ind(domain_size, Idomain, Jdomain, Kdomain);
    end


end