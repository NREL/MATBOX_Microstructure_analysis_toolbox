function [microstructure] = generate_coil(coil_diameter,coil_thickness,coil_length,rotation_speed,additional_distance_throughplane,additional_distance_inplane,coil_is_background)


%% SET DISTANCE AND DOMAIN SIZE
coil_radius = coil_diameter/2;
additional_distance_throughplane=round(additional_distance_throughplane*coil_length);
additional_distance_inplane=round(additional_distance_inplane*coil_diameter);

% Domain size
n_1 = round(coil_length + additional_distance_throughplane);
n_2 = round(coil_diameter + 2*additional_distance_inplane);
n_3=n_2;

% Domain size
Domain_size=[n_1 n_2 n_3];

%% ASSIGN VOXELS
% Initialise microstructure
microstructure = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% Coil center in the plane
a=(n_2/2);
b=(n_3/2);
theta=0;
incrementpervoxel_coil_rotation = rotation_speed * 2*pi/coil_diameter;
for x=Domain_size(1)-additional_distance_throughplane:-1:1
    x_ = max(a + coil_radius*cos(theta),0.51);
    y_ = max(b + coil_radius*sin(theta),0.51);    
    microstructure(round(x),round(x_),round(y_))=1;
    theta=theta+incrementpervoxel_coil_rotation;
    if theta>=2*pi
        theta=0;
    end
end

%% GENERATE PARTICLES

% Calculate Euclidean distance map
Distance_map = bwdist(microstructure);
% Distance lower than the radius means the voxel belong to one particle
microstructure( Distance_map<=round(coil_thickness/2)-1) = 1;

%% COIL IS BACKGROUN
if coil_is_background
    id_coil = microstructure==1;
    microstructure = ones(size(microstructure));
    microstructure(id_coil)=0;
end

microstructure=uint8(microstructure); % Convert in 8 bits
    
end

