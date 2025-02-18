function [Particle_size] = Calculate_particlesize_fromlabel(Labelmap,voxel_size)


[C,ia,ic] = unique(Labelmap);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
id_background = find(C(:,1)==0);
if ~isempty(id_background)
    value_counts(id_background,:)=[];
end
[number_particle,~]=size(value_counts);

% Initialization
Particle_size = zeros(size(Labelmap));
for k=1:1:number_particle
    Particle_size(Particle_size==value_counts(k,1))=value_counts(k,2);
end



% Initialisation
sz = size(Labelmap);
Particle_size = zeros(sz);
% Unique id
unique_particle = unique(Labelmap);
unique_particle(unique_particle==0)=[]; % Remove background
% Number of different particle
number_particle = length(unique_particle);
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    % Find all voxels of the particle
    index_particle=find(Labelmap==particle_id);
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
    Particle_size(index_particle)= equivalent_diameter_size;
end


end