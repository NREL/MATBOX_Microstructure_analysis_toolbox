function [Connectivity_particle,interface_location] = Function_particleconnectivity(label_lake)

Domain_size=size(label_lake); % Domain size of the microstructure
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end

% Unique id
unique_particle = unique(label_lake);
unique_particle(unique_particle==0)=[]; % Remove background
% Number of different particle
number_particle = length(unique_particle);

% Interface initialisation
Interface_particle = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
interface_particleparticle=3; % interface between 2 particle of the phase
interface_particlecomplemetary=2; % interface between particle of the phase and the complementary phase

% Id interface is calculated with a pairing function (Cantor)
interface_location = zeros(Domain_size(1)*Domain_size(2)*Domain_size(3),4,'uint16'); % Pairing number (unique), X,Y,Z

% Connectivity initialisation
Connectivity_particle = zeros(number_particle+1,number_particle+1);
% Set table axis
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    % Write axis
    Connectivity_particle(current_+1,1)=particle_id;
    Connectivity_particle(1,current_+1)=particle_id;
    % Write identity
    Connectivity_particle(current_+1,current_+1)=1;
end
% Connectivity_particle(i,j) is the number of voxels identified at the interface (if any) between the two particles
% with a face-to-face identification for the interface, then the number is equal to the number of voxel face at the interface

number_voxel_atinterface=0; % Initialize
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    % Find all voxels of the particle
    index_particle=find(label_lake==particle_id);
    % Get all their coordinate
    [Pa1,Pa2,Pa3] = ind2sub(Domain_size,index_particle);
    % Number of voxel
    number_voxel_particle = length(Pa1);
    % Loop over all voxels of the phase
    for current_voxel=1:1:number_voxel_particle
        % Get back voxel coordinate
        x_=Pa1(current_voxel); y_=Pa2(current_voxel); z_=Pa3(current_voxel);
        
        % Check x-
        if x_>1
            % Tested voxel
            x_tested = x_-1; y_tested = y_; z_tested = z_;
            % Value founded there
            Id_tested = label_lake(x_tested,y_tested,z_tested);
            if Id_tested~=particle_id
                % Interface founded
                if Id_tested==0
                    % Interface with the complementary phase
                    Interface_particle(x_,y_,z_)=interface_particlecomplemetary;
                else
                    % Interface with the another particle of the phase
                    Interface_particle(x_,y_,z_)=interface_particleparticle;
                    % Update Connectivity_particle as consequence
                    row_ = find( Connectivity_particle(:,1)==particle_id);
                    column_ = find( Connectivity_particle(1,:)==Id_tested);
                    Connectivity_particle(row_,column_)=Connectivity_particle(row_,column_)+1;
                    
                    % Interface id
                    k1=min(Id_tested,particle_id);
                    k2=max(Id_tested,particle_id);
                    pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2;
                    number_voxel_atinterface=number_voxel_atinterface+1; % Increment
                    interface_location(number_voxel_atinterface,1) = pairing_value;
                    interface_location(number_voxel_atinterface,2) = x_;
                    interface_location(number_voxel_atinterface,3) = y_;
                    interface_location(number_voxel_atinterface,4) = z_;
                    
                end
            end
        end
        
        % Check y-
        if y_>1
            % Tested voxel
            x_tested = x_; y_tested = y_-1; z_tested = z_;
            % Value founded there
            Id_tested = label_lake(x_tested,y_tested,z_tested);
            if Id_tested~=particle_id
                % Interface founded
                if Id_tested==0
                    % Interface with the complementary phase
                    Interface_particle(x_,y_,z_)=interface_particlecomplemetary;
                else
                    % Interface with the another particle of the phase
                    Interface_particle(x_,y_,z_)=interface_particleparticle;
                    % Update Connectivity_particle as consequence
                    row_ = find( Connectivity_particle(:,1)==particle_id);
                    column_ = find( Connectivity_particle(1,:)==Id_tested);
                    Connectivity_particle(row_,column_)=Connectivity_particle(row_,column_)+1;
                    
                    % Interface id
                    k1=min(Id_tested,particle_id);
                    k2=max(Id_tested,particle_id);
                    pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2;
                    number_voxel_atinterface=number_voxel_atinterface+1; % Increment
                    interface_location(number_voxel_atinterface,1) = pairing_value;
                    interface_location(number_voxel_atinterface,2) = x_;
                    interface_location(number_voxel_atinterface,3) = y_;
                    interface_location(number_voxel_atinterface,4) = z_;                   
                    
                end
            end
        end
        
        % Check z-
        if z_>1
            % Tested voxel
            x_tested = x_; y_tested = y_; z_tested = z_-1;
            % Value founded there
            Id_tested = label_lake(x_tested,y_tested,z_tested);
            if Id_tested~=particle_id
                % Interface founded
                if Id_tested==0
                    % Interface with the complementary phase
                    Interface_particle(x_,y_,z_)=interface_particlecomplemetary;
                else
                    % Interface with the another particle of the phase
                    Interface_particle(x_,y_,z_)=interface_particleparticle;
                    % Update Connectivity_particle as consequence
                    row_ = find( Connectivity_particle(:,1)==particle_id);
                    column_ = find( Connectivity_particle(1,:)==Id_tested);
                    Connectivity_particle(row_,column_)=Connectivity_particle(row_,column_)+1;
                    
                    % Interface id
                    k1=min(Id_tested,particle_id);
                    k2=max(Id_tested,particle_id);
                    pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2;
                    number_voxel_atinterface=number_voxel_atinterface+1; % Increment
                    interface_location(number_voxel_atinterface,1) = pairing_value;
                    interface_location(number_voxel_atinterface,2) = x_;
                    interface_location(number_voxel_atinterface,3) = y_;
                    interface_location(number_voxel_atinterface,4) = z_;               
                    
                end
            end
        end
        
        % Check x+
        if x_<Domain_size(1)
            % Tested voxel
            x_tested = x_+1; y_tested = y_; z_tested = z_;
            % Value founded there
            Id_tested = label_lake(x_tested,y_tested,z_tested);
            if Id_tested~=particle_id
                % Interface founded
                if Id_tested==0
                    % Interface with the complementary phase
                    Interface_particle(x_,y_,z_)=interface_particlecomplemetary;
                else
                    % Interface with the another particle of the phase
                    Interface_particle(x_,y_,z_)=interface_particleparticle;
                    % Update Connectivity_particle as consequence
                    row_ = find( Connectivity_particle(:,1)==particle_id);
                    column_ = find( Connectivity_particle(1,:)==Id_tested);
                    Connectivity_particle(row_,column_)=Connectivity_particle(row_,column_)+1;
                    
                    % Interface id
                    k1=min(Id_tested,particle_id);
                    k2=max(Id_tested,particle_id);
                    pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2;
                    number_voxel_atinterface=number_voxel_atinterface+1; % Increment
                    interface_location(number_voxel_atinterface,1) = pairing_value;
                    interface_location(number_voxel_atinterface,2) = x_;
                    interface_location(number_voxel_atinterface,3) = y_;
                    interface_location(number_voxel_atinterface,4) = z_;                 
                    
                end
            end
        end
        
        % Check y+
        if y_<Domain_size(2)
            % Tested voxel
            x_tested = x_; y_tested = y_+1; z_tested = z_;
            % Value founded there
            Id_tested = label_lake(x_tested,y_tested,z_tested);
            if Id_tested~=particle_id
                % Interface founded
                if Id_tested==0
                    % Interface with the complementary phase
                    Interface_particle(x_,y_,z_)=interface_particlecomplemetary;
                else
                    % Interface with the another particle of the phase
                    Interface_particle(x_,y_,z_)=interface_particleparticle;
                    % Update Connectivity_particle as consequence
                    row_ = find( Connectivity_particle(:,1)==particle_id);
                    column_ = find( Connectivity_particle(1,:)==Id_tested);
                    Connectivity_particle(row_,column_)=Connectivity_particle(row_,column_)+1;
                    
                    % Interface id
                    k1=min(Id_tested,particle_id);
                    k2=max(Id_tested,particle_id);
                    pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2;
                    number_voxel_atinterface=number_voxel_atinterface+1; % Increment
                    interface_location(number_voxel_atinterface,1) = pairing_value;
                    interface_location(number_voxel_atinterface,2) = x_;
                    interface_location(number_voxel_atinterface,3) = y_;
                    interface_location(number_voxel_atinterface,4) = z_;                  
                    
                end
            end
        end
        
        % Check z+
        if z_<Domain_size(3)
            % Tested voxel
            x_tested = x_; y_tested = y_; z_tested = z_+1;
            % Value founded there
            Id_tested = label_lake(x_tested,y_tested,z_tested);
            if Id_tested~=particle_id
                % Interface founded
                if Id_tested==0
                    % Interface with the complementary phase
                    Interface_particle(x_,y_,z_)=interface_particlecomplemetary;
                else
                    % Interface with the another particle of the phase
                    Interface_particle(x_,y_,z_)=interface_particleparticle;
                    % Update Connectivity_particle as consequence
                    row_ = find( Connectivity_particle(:,1)==particle_id);
                    column_ = find( Connectivity_particle(1,:)==Id_tested);
                    Connectivity_particle(row_,column_)=Connectivity_particle(row_,column_)+1;
                    
                    % Interface id
                    k1=min(Id_tested,particle_id);
                    k2=max(Id_tested,particle_id);
                    pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2;
                    number_voxel_atinterface=number_voxel_atinterface+1; % Increment
                    interface_location(number_voxel_atinterface,1) = pairing_value;
                    interface_location(number_voxel_atinterface,2) = x_;
                    interface_location(number_voxel_atinterface,3) = y_;
                    interface_location(number_voxel_atinterface,4) = z_;                   
                    
                end
            end
        end
        
    end
end
% Check the matrix is symmetric
Connectivity_particle_symmetry = issymmetric(Connectivity_particle);
if Connectivity_particle_symmetry==1
    disp 'Connectivity between particle is symmetrical (Normal)'
else
    disp 'Connectivity between particle is non-symmetrical (Error)'
end


end

