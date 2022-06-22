function [PhaseLabel, ParticleId] = function_createvolume_fromparticledata(particle_data,domain_size,scale_diameterratio,force_separation,distance_separation)

PhaseLabel = zeros(domain_size); % Initialize volume
ParticleId = zeros(domain_size); % Initialize volume

[number_particle,~] = size(particle_data);

for k_particle = 1:1:number_particle % Loop over all particles

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

    % Create particle
    binary_ellipsoid = create_ellipsoid(dx,dy,dz); % Create ellipsoid
    if angle_x_deg~=0 || angle_y_deg~=0 || angle_z_deg~=0 % Only if one or more rotation
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

    else % Subdomain contains other particles
        subdomain_binary = subdomain_label;
        subdomain_binary(subdomain_binary~=0)=1;
        union_subdomain_ellipsoid = subdomain_binary + subdomain_ellipsoid;
        idx_overlapping = find(union_subdomain_ellipsoid==2);
        if isempty(idx_overlapping) % There is no overlapping
            PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_label + subdomain_ellipsoid * label; % Insert ellipsoid within the domain, in addition to existing particles
            ParticleId(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_id + subdomain_ellipsoid * id;
        else % Overlapping, we need to distribute conflicting voxels based on nearest particle
            subdomain_union_minus_intersection = abs(subdomain_ellipsoid - subdomain_binary); % (particle union others) - (particle intersection others)
            [Euclidean_distance_map,Idx_distance_map] = bwdist(subdomain_union_minus_intersection);
            index_nearest_particle = Idx_distance_map(idx_overlapping);

            Sub_label = subdomain_label + subdomain_ellipsoid * label;
            Sub_Id = subdomain_id + subdomain_ellipsoid * id;
            Sub_label(idx_overlapping) = Sub_label(index_nearest_particle);
            Sub_Id(idx_overlapping) = Sub_Id(index_nearest_particle);

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

            PhaseLabel(x_min:x_max, y_min:y_max, z_min:z_max) = Sub_label; % Insert ellipsoid within the domain, in addition to existing particles
            ParticleId(x_min:x_max, y_min:y_max, z_min:z_max) = Sub_Id;

        end
    end    
end


    function [binary_ellipsoid] = create_ellipsoid(dx,dy,dz)
        % Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix

        Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
        [X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
        R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
        binary_ellipsoid = zeros(size(X));
        binary_ellipsoid(R <= 1 ) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
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


end