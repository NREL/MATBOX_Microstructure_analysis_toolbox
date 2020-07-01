function [results_correlation, results_visualization, Table_discrete_particleinterface_analysis, Connectivity_particle] = Function_analyse_discrete_particle_interface(Particle_size, Particle_label, unique_particle, mass_center_particle, results_correlation, results_visualization, current_phase, Current_folder, INFO, OPTIONS, dpsdname)

voxel_size = INFO.voxel_size;
% Number of different particle
number_particle = length(unique_particle);
% Domain size
Domain_size = size(Particle_size);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end


% Interface initialisation
Interface_particle = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
interface_particleparticle=3; % interface between 2 particle of the phase
interface_particlecomplemetary=2; % interface between particle of the phase and the complementary phase
interface_bulk=1; % phase
interface_complementary=0; % complementary phase

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
    index_particle=find(Particle_label==particle_id);
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
            Id_tested = Particle_label(x_tested,y_tested,z_tested);
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
            Id_tested = Particle_label(x_tested,y_tested,z_tested);
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
            Id_tested = Particle_label(x_tested,y_tested,z_tested);
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
            Id_tested = Particle_label(x_tested,y_tested,z_tested);
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
            Id_tested = Particle_label(x_tested,y_tested,z_tested);
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
            Id_tested = Particle_label(x_tested,y_tested,z_tested);
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

results_visualization(current_phase).(['Interface_particle_' dpsdname]) = Interface_particle;
%if strcmp(dpsdname,'watershed')
%    results_visualization(current_phase).Interface_particle = Interface_particle;
%end

% Number of connection
Number_connection=zeros(number_particle,3);
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    Number_connection(current_,1)=particle_id;
    % Get size
    Number_connection(current_,2)=mass_center_particle(current_,6);
    % Get number of connection
    row_ = find( Connectivity_particle(:,1)==particle_id);
    n_connection = sum(Connectivity_particle(row_,:)~=0)-1;
    Number_connection(current_,3)=n_connection;
end

Variable_name_table={'Particle_Id' 'Equivalent_diameter' 'Number_connection'};
Table_connection = array2table(Number_connection,...
    'VariableNames',Variable_name_table);

% Save
Table_discrete_particleinterface_analysis.Interface_particle=Interface_particle;
Table_discrete_particleinterface_analysis.code_interface_particleparticle=interface_particleparticle;
Table_discrete_particleinterface_analysis.code_interface_particlecomplemetary=interface_particlecomplemetary;
Table_discrete_particleinterface_analysis.code_bulk=interface_bulk;
Table_discrete_particleinterface_analysis.code_complementary=interface_complementary;
Table_discrete_particleinterface_analysis.Connectivity_particle=Connectivity_particle;
Table_discrete_particleinterface_analysis.Table_connection=Table_connection;

% % Text Results are saved
if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Particle_connection_' dpsdname '_' INFO.phase(current_phase).name];
    full_path=[Current_folder filename];
    % Prepare the data
    clear DATA_writetable;
    % Data : Size ratio
    DATA_writetable.sheet(1).name='Connection';
    DATA_writetable.sheet(1).table=Table_connection;    
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

% Connectivity as a function of particle size
% Particle with only two connections are likely to be considered as
% connectors.
% Then, this plot will tell you if you have captured these connections, or
% if they have been distributed between the larger particles

% - Create figure
Fig = figure;
Fig.Name= sprintf('Number of particle connection, %s',INFO.phase(current_phase).name);
Fig.Color='white'; % Background colour
% - Create axes
axes_ = axes('Parent',Fig);
hold(axes_,'on');
% - Title
t_=title (' ','FontName','Times New Roman','FontSize',16);
t_.String= sprintf('Number of particle connection, %s',INFO.phase(current_phase).name);
h_pointcloud = scatter(Number_connection(:,2),Number_connection(:,3));
% Colors, thickness, markers
% - Axis label
t_ = xlabel(' ');
t_1 = sprintf('Equivalent particle diameter');
t_2 = ' (\mum)';
t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
t_ = ylabel('Number of connection');

% - Grid
if strcmp(OPTIONS.grid,'on')
    grid(axes_,'on'); % Display grid
    set(axes_,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
end
% - Fontname and fontsize
set(axes_,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
t_.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
% - Figure has been done
hold(axes_,'off');
axis tight
if OPTIONS.save_fig == true % Save figure
    filename= sprintf('Particle_connection_diameter_%s_%s',dpsdname,INFO.phase(current_phase).name);
    function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end


%% INTERFACES DIMENSION

% Remove 0
interface_location(number_voxel_atinterface+1:end,:)=[];
% Unique pairing values
unique_paring_values = unique(interface_location(:,1));
% Number of interfaces
number_interface=length(unique_paring_values);

% Allocate result array
% pairing id, meanx , meany, meanz, normalx, normaly, normalz, number of face, area um2, equivalent disc diamter um, area um2 x 2/3, corrected equivalent disc diameter um
Interface_results = zeros(number_interface,12);

% Pairing id
Interface_results(:,1) = unique_paring_values;

% Mean position and normal
for current_interface=1:1:number_interface % Loop over all interfaces
    current_pairing_value = unique_paring_values(current_interface);
    index_=find(interface_location(:,1)==current_pairing_value); % Get all linear index
    all_x = interface_location(index_,2); % Get all voxels coordinates
    all_y = interface_location(index_,3);
    all_z = interface_location(index_,4);
    
    mean_x = mean(all_x); % Get mean position
    mean_y = mean(all_y);
    mean_z = mean(all_z);
    Interface_results(current_interface,2) = mean_x; % And save it
    Interface_results(current_interface,3) = mean_y;
    Interface_results(current_interface,4) = mean_z;
    
    if length(all_x)>=3
        normal = function_fitplane3dpoints(double([all_x all_y all_z]), 'Orthogonalregression_PCA', false); % Plane coefficients
        Interface_results(current_interface,5) = normal(1); % And save it
        Interface_results(current_interface,6) = normal(2);
        Interface_results(current_interface,7) = normal(3);
    else
        Interface_results(current_interface,5) = 1/(3^0.5); % Not enough point to perform a fit
        Interface_results(current_interface,6) = 1/(3^0.5);
        Interface_results(current_interface,7) = 1/(3^0.5);
    end
end

% Number of face
for row=2:number_particle
    for column=row+1:number_particle+1
        number_of_facets_interface = Connectivity_particle(row,column);
        if number_of_facets_interface>0 % Interface exist
            % Get unique interface id
            particle_id_1 = Connectivity_particle(row,1);
            particle_id_2 = Connectivity_particle(1,column);
            k1=particle_id_1;
            k2=particle_id_2;
            pairing_value = 0.5*(k1+k2)*(k1+k2+1)+k2; % k1<k2
            index_=find(Interface_results(:,1)==pairing_value);
            Interface_results(index_,8) = number_of_facets_interface; % And save it
        end
    end
end

% Area
Interface_results(:,9) = Interface_results(:,8)* (voxel_size/1000) * (voxel_size/1000) ; % um2

% Equivalent disc diameter
Interface_results(:,10) = 2 * ((Interface_results(:,9) ./ pi).^(0.5)); % um

% Corrected area (spherical assumption)
Interface_results(:,11) = Interface_results(:,9) * 2/3; %um2

% Corrected equivalent disc diameter
Interface_results(:,12) = 2 * ((Interface_results(:,11) ./ pi).^(0.5)); % um

% Deduce mean diameter along each direction
Interface_meandiameter_alongdirection = zeros(3,3); % 3rd column is with corrected diameter
Interface_meanarea_alongdirection = zeros(3,3); % 3rd column is with corrected diameter
for direction=1:1:3
    Interface_meandiameter_alongdirection(direction,1)=direction;
    Interface_meanarea_alongdirection(direction,1)=direction;
end


area_ = Interface_results(:,9);
corrected_area_ = Interface_results(:,11);
diameter_ = Interface_results(:,10);
corrected_diameter_ = Interface_results(:,12);
weight_dir1 = abs(Interface_results(:,5));
weight_dir2 = abs(Interface_results(:,6));
weight_dir3 = abs(Interface_results(:,7));

Interface_meandiameter_alongdirection(1,2) = sum( diameter_.*weight_dir1 ) / sum( weight_dir1 );
Interface_meandiameter_alongdirection(2,2) = sum( diameter_.*weight_dir2 ) / sum( weight_dir2 );
Interface_meandiameter_alongdirection(3,2) = sum( diameter_.*weight_dir3 ) / sum( weight_dir3 );
Interface_meandiameter_alongdirection(1,3) = sum( corrected_diameter_.*weight_dir1 ) / sum( weight_dir1 );
Interface_meandiameter_alongdirection(2,3) = sum( corrected_diameter_.*weight_dir2 ) / sum( weight_dir2 );
Interface_meandiameter_alongdirection(3,3) = sum( corrected_diameter_.*weight_dir3 ) / sum( weight_dir3 );

Interface_meanarea_alongdirection(1,2) = sum( area_.*weight_dir1 ) / sum( weight_dir1 );
Interface_meanarea_alongdirection(2,2) = sum( area_.*weight_dir2 ) / sum( weight_dir2 );
Interface_meanarea_alongdirection(3,2) = sum( area_.*weight_dir3 ) / sum( weight_dir3 );
Interface_meanarea_alongdirection(1,3) = sum( corrected_area_.*weight_dir1 ) / sum( weight_dir1 );
Interface_meanarea_alongdirection(2,3) = sum( corrected_area_.*weight_dir2 ) / sum( weight_dir2 );
Interface_meanarea_alongdirection(3,3) = sum( corrected_area_.*weight_dir3 ) / sum( weight_dir3 );

results_correlation(current_phase).(['ParticleInterfaceArea_along_dir1_' dpsdname]) = Interface_meanarea_alongdirection(1,3);
results_correlation(current_phase).(['ParticleInterfaceArea_along_dir2_' dpsdname]) = Interface_meanarea_alongdirection(2,3);
results_correlation(current_phase).(['ParticleInterfaceArea_along_dir3_' dpsdname]) = Interface_meanarea_alongdirection(3,3);
%if strcmp(dpsdname,'watershed')
%    results_correlation(current_phase).ParticleInterfaceArea_along_dir1_Watershed = Interface_meanarea_alongdirection(1,3);
%    results_correlation(current_phase).ParticleInterfaceArea_along_dir2_Watershed = Interface_meanarea_alongdirection(2,3);
%    results_correlation(current_phase).ParticleInterfaceArea_along_dir3_Watershed = Interface_meanarea_alongdirection(3,3);
%end

% % Save results
clear Variable_name_table;
Variable_name_table={'Direction' 'Weighted_diameter_um' 'Weighted_corrected_diameter_um'};
Table_interface_diameter = array2table(Interface_meandiameter_alongdirection,...
    'VariableNames',Variable_name_table);

clear Variable_name_table;
Variable_name_table={'Direction' 'Weighted_area_um2' 'Weighted_corrected_area_um2'};
Table_interface_area = array2table(Interface_meanarea_alongdirection,...
    'VariableNames',Variable_name_table);

Table_discrete_particleinterface_analysis.Interface_meandiameter_alongdirection = Table_interface_diameter;
Table_discrete_particleinterface_analysis.Interface_meanarea_alongdirection = Table_interface_area;

if OPTIONS.displaytext==true
    fprintf('PARTICLE-to-PARTICLE CONNECTION DIMENSION (area)\n');
    fprintf('   sum(equivalent disc area * interface plane normal direction component) / sum(interface plane normal direction component)\n');
    fprintf('   i.e. boundaries for which fitted plane normal is not aligned with the direction have less weight\n');
    fprintf('   3rd column: area multiplied by 2/3 (spherical assumption)\n\n');
    disp(Table_interface_area)
    fprintf('PARTICLE-to-PARTICLE CONNECTION DIMENSION (diameter)\n');
    fprintf('   sum(equivalent disc diameter * interface plane normal direction component) / sum(interface plane normal direction component)\n');
    fprintf('   Disc diameter deduced from the area analysis\n\n');
    disp(Table_interface_diameter)
end

if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Particle_Interface_dimension_' dpsdname '_' INFO.phase(current_phase).name];
    full_path=[Current_folder filename];
    % Prepare the data
    clear DATA_writetable;
    DATA_writetable.sheet(1).name='Diameter';
    DATA_writetable.sheet(1).table=Table_interface_diameter;
    DATA_writetable.sheet(2).name='Area';
    DATA_writetable.sheet(2).table=Table_interface_area;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


end
