function [microstructure] = generate_spheres(p)

%% PARAMETERS
number_particles_through_plane = p.number_particles_through_plane;
number_particles_in_plane = p.number_particles_in_plane;
cut_particle_inplane_extremities = p.cut_particle_inplane_extremities;
particle_diameter = p.particle_diameter;
distance_between_center = p.distance_between_center;
additional_distance_throughplane = p.additional_distance_throughplane;
additional_distance_inplane = p.additional_distance_inplane;
clear p;

%% Dimension in voxel

%particle_diameter = 10; % micrometers
%voxel_size = 0.2; % micrometers
%distance_between_center = 0.9; % Distance between particles center (expressed in particle diameter)
%additional_distance_throughplane=0.5; % Distance betwenn particle extremity to domain boundary (expressed in particle diameter)
%additional_distance_inplane=0.5; % Distance betwenn particle extremity to domain boundary (expressed in particle diameter)

%number_particles_through_plane=4; % Number of particles along through plane
%number_particles_in_plane=4; % Number of particles in-plane
%cut_particle_inplane_extremities=1; % Do you want to cut particles at their middle (in-plane)
% 
% % Saving options
% sav_folder='C:\Users\fussegli\Desktop\Matlab_data\Electrode\'; % Save folder
% name_electrode='Particle_3pTP_3pIP_contactSEP_200um'; % Filename
% 

if number_particles_in_plane==1
    cut_particle_inplane_extremities=false; % Overwritte
end

%% SET DISTANCE AND DOMAIN SIZE
% Number of voxel for 1 particle diameter
%particle_diameter=particle_diameter/voxel_size;
% Convert distance in number of voxel
distance_between_center=round(distance_between_center*particle_diameter);
additional_distance_throughplane=round(additional_distance_throughplane*particle_diameter);
additional_distance_inplane=round(additional_distance_inplane*particle_diameter);
% Domain size
n_1 = round(distance_between_center*(number_particles_through_plane-1) + particle_diameter/2 + additional_distance_throughplane);
n_2 = round(distance_between_center*(number_particles_in_plane-1) + particle_diameter + 2*additional_distance_inplane);
n_3=n_2;
% Domain size
Domain_size=[n_1 n_2 n_3];

%% ASSIGN PARTICLE CENTERS
% Initialise microstructure
microstructure = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
all_centerparticle=[];
% Assign center of the particle
for x_position=1:1:number_particles_through_plane
    % along thickness coordinate
    x=round(((x_position-1)*distance_between_center)+1);
    for y_position=1:1:number_particles_in_plane
        % First in-plane coordinate
        y= round( additional_distance_inplane + particle_diameter/2 +  (y_position-1)*distance_between_center );
        for z_position=1:1:number_particles_in_plane
            % Second in-plane coordinate
            z= round( additional_distance_inplane + particle_diameter/2 +  (z_position-1)*distance_between_center );
            % Assign particle center
            center_particle = [x y z];
            all_centerparticle=[all_centerparticle; center_particle];
            microstructure(center_particle(1),center_particle(2),center_particle(3))=1;
        end
    end
end

%% CUT
if cut_particle_inplane_extremities==1
    min_y=min(all_centerparticle(:,2));
    max_y=max(all_centerparticle(:,2));
    min_z=min(all_centerparticle(:,3));
    max_z=max(all_centerparticle(:,3));    
    % Crop
    microstructure=microstructure(:,min_y:max_y,min_z:max_z);
end

%% GENERATE PARTICLES
% Calculate Euclidean distance map
Distance_map = bwdist(microstructure);
% Distance lower than the radius means the voxel belong to one particle
microstructure( Distance_map<=round(particle_diameter/2)-1 ) = 1;

microstructure=uint8(microstructure); % Convert in 8 bits


%%
% keyboard
% anode_id =1;
% cathode_id =2;
% separator_id = 4; % 4
% electrolyte_id = 3;
% 
% sz = size(microstructure);
% % a=zeros(5,sz(2),sz(3))+1;
% % microstructure = [a; microstructure];
% 
% domain_size = size(microstructure);
% tmp=zeros(domain_size,'uint8');
% for k=1:1:domain_size(1)
%     slice=microstructure(k,:,:);
%     tmp(domain_size(1)-k+1,:,:)=slice;
% end
% 
% microstructure(microstructure==1)=cathode_id; % Cathode
% tmp(tmp==1)=anode_id; % anode
% sep = zeros(20,sz(2),sz(3))+separator_id;
% 
% microstructure = [microstructure;sep;tmp];
% microstructure(microstructure==0)=electrolyte_id; % Electrolyte
% 
% figure
% microstructure(1:13,:,:) = [];
% %microstructure(205:end,:,:) = [];
% microstructure(122:end,:,:) = [];
% 
% sz = size(microstructure);
% imagesc(microstructure(:,:,round(sz(3)/2)));
% % 
% % a=microstructure(:,1:9,:);
% % a(a==electrolyte_id)=4;
% % microstructure(:,1:9,:)=a;
% %  
% function_save_tif(microstructure,'Periodicspheres_beforeCBD.tif');
% % function_save_tif(microstructure,'Periodic_spheres_withpad_withsep_smallrev.tif');
% % function_save_tif(microstructure,'Periodic_spheres_withpad_withsep_smallpar.tif');


end

