function [results_correlation, results_visualization, Table_discrete_particle_graphnetowrk] = Function_analyse_discrete_particle_graphnetwork(Particle_size, Particle_label, unique_particle, mass_center_particle, Connectivity_particle, results_correlation, results_visualization, current_phase, Current_folder, INFO, OPTIONS, dpsdname)

%%
%% WARNING NEED REFACTORING !
%%

voxel_size = INFO.voxel_size;
% Number of different particle
number_particle = length(unique_particle);
% Domain size
Domain_size = size(Particle_size);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end

%%
%% SKELETON
%%

%% 2D-3D SKELETON: PREPARE DATA

if Domain_size(3)==1
    % % 2D Preparation: Prepare background image
    % Initialisation
    skeleton_ = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    % Phase
    skeleton_(binary_phase==1)=1;
    % Interface
    skeleton_(Interface_particle==interface_particleparticle)=2;
    skeleton_(Interface_particle==interface_particlecomplemetary)=3;
    % % 2D Preparation: equivalent circle
    % loop over all particles
    circle_data = zeros(number_particle,3);
    for current_=1:1:number_particle
        circle_data(current_,1)=mass_center_particle(current_,3); % x coordinate of circle center
        circle_data(current_,2)=mass_center_particle(current_,2); % y coordinate of circle center
        circle_data(current_,3)=mass_center_particle(current_,6)/2; % Radius (micrometer)
    end
end

% % 2D-3D Preparation: Scatter point for center of mass
% 3D scatter
x_scatter3 = zeros(number_particle,1);
y_scatter3 = zeros(number_particle,1);
z_scatter3 = zeros(number_particle,1);
% loop over all particles
for current_=1:1:number_particle
    % 3D scatter
    x_scatter3(current_,1)=mass_center_particle(current_,2);
    y_scatter3(current_,1)=mass_center_particle(current_,3);
    z_scatter3(current_,1)=mass_center_particle(current_,4);
end

% % 2D-3D Preparation: Point coordinate for skeleton line defined simply as line joining mass center of connected particles
number_line=0;
skeleton_line=zeros(3,2,number_particle*number_particle);
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id_1 = unique_particle(current_);
    % Get row
    row_ = find( Connectivity_particle(:,1)==particle_id_1);
    % loop over all other particles
    for column=row_+1:1:number_particle
        connection_ = Connectivity_particle(row_,column);
        if connection_>0
            number_line=number_line+1;
            % Particle id 1 and id 2 are connected
            particle_id_2 = Connectivity_particle(1,column);
            % Find associate coordinate
            index_1 = find( mass_center_particle(:,1)==particle_id_1);
            index_2 = find( mass_center_particle(:,1)==particle_id_2);
            x_1 = mass_center_particle(index_1,2); x_2 = mass_center_particle(index_2,2);
            y_1 = mass_center_particle(index_1,3); y_2 = mass_center_particle(index_2,3);
            z_1 = mass_center_particle(index_1,4); z_2 = mass_center_particle(index_2,4);
            % Create vector
            skeleton_line(1,1,number_line)=x_1; skeleton_line(1,2,number_line)=x_2;
            skeleton_line(2,1,number_line)=y_1; skeleton_line(2,2,number_line)=y_2;
            skeleton_line(3,1,number_line)=z_1; skeleton_line(3,2,number_line)=z_2;
        end
    end
end


%%
%% SKELETON IMAGE
%%

%% 2D CASE
if Domain_size(3)==1
    % - Create figure
    Fig = figure;
    Fig.Name= sprintf('Skeleton of phase: %s',INFO.phase(current_phase).name);
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= sprintf('Skeleton, %s',INFO.phase(current_phase).name);
    % Background
    image(skeleton_(:,:,1),'CDataMapping','scaled')
    % Skeleton line defined simply as line joining mass center of connected particles
    for current_line=1:1:number_line
        y_line = [skeleton_line(1,1,current_line) skeleton_line(1,2,current_line)];
        x_line = [skeleton_line(2,1,current_line) skeleton_line(2,2,current_line)];
        z_line = [skeleton_line(3,1,current_line) skeleton_line(3,2,current_line)];
        h_skeleton_line = plot3(x_line,y_line,z_line);
        h_skeleton_line.LineStyle='-';
        h_skeleton_line.Color='w';
        h_skeleton_line.LineWidth=2;
    end
    % Mass center
    h_center = scatter3(y_scatter3,x_scatter3,z_scatter3,'+','MarkerEdgeColor','r','LineWidth',2);
    % Circle
    for current_=1:1:number_particle
        h_circle = viscircles([circle_data(current_,1),circle_data(current_,2)],circle_data(current_,3),'Color','r','LineWidth',2);
    end
    % - Axis label
    t_ = xlabel(' ');
    t_1 = ['Position along ' INFO.direction(1).name];
    t_2 = ' (\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' ');
    t_1 = ['Position along ' INFO.direction(2).name];
    t_2 = ' (\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    % - Fontname and fontsize
    set(axes_,'FontName','Times New Roman','FontSize',18);
    % - Figure has been done
    hold(axes_,'off');
    axis equal
    axis tight
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('Skeleton_%s_%s',dpsdname,INFO.phase(current_phase).name);
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end

%% 3D CASE
if Domain_size(3)>1
    
    % - Create figure
    Fig = figure;
    Fig.Name= sprintf('Skeleton, %s',INFO.phase(current_phase).name);
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= sprintf('Skeleton, %s',INFO.phase(current_phase).name);
    % Skeleton line defined simply as line joining mass center of connected particles
    for current_line=1:1:number_line
        y_line = [skeleton_line(1,1,current_line) skeleton_line(1,2,current_line)];
        x_line = [skeleton_line(2,1,current_line) skeleton_line(2,2,current_line)];
        z_line = [skeleton_line(3,1,current_line) skeleton_line(3,2,current_line)];
        h_skeleton_line = plot3(x_line,y_line,z_line);
        h_skeleton_line.LineStyle='-';
        h_skeleton_line.Color='b';
        h_skeleton_line.LineWidth=2;
    end
    % Mass center
    h_center = scatter3(y_scatter3,x_scatter3,z_scatter3,'+','MarkerEdgeColor','r','LineWidth',2);
    % Circle
    %for current_=1:1:number_particle
    %    h_circle = viscircles([circle_data(current_,1),circle_data(current_,2)],circle_data(current_,3),'Color','r','LineWidth',2);
    %end
    % - Axis label
    t_ = xlabel(' ');
    t_1 = ['Position along ' INFO.direction(1).name];
    t_2 = ' (\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' ');
    t_1 = ['Position along ' INFO.direction(2).name];
    t_2 = ' (\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = zlabel(' ');
    t_1 = ['Position along ' INFO.direction(3).name];
    t_2 = ' (\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    % - Fontname and fontsize
    set(axes_,'FontName','Times New Roman','FontSize',18);
    % - Grid
    grid(axes_,'on'); % Display grid
    set(axes_,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    % View
    view(3)
    % - Figure has been done
    hold(axes_,'off');
    axis equal
    axis tight
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('Skeleton_%s_%s',dpsdname,INFO.phase(current_phase).name);
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
    
end


%%
%% GRAPH
%%

% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;

% Skeleton graph
% Geometric tortuosity (from skeleton graph)
% Extract only the required data
node_coordinate = mass_center_particle;
node_connectivity = Connectivity_particle;
% Call the function to create the graph
[G] = Function_create_graph_from_skeleton(node_coordinate,node_connectivity);

% Characteristic length
%[all_distances] = Function_charateristic_length_from_anygraph(G_bis);

% Loop over direction
for direction=1:1:3
    
    % Calculate shortest path (normalizd with the domain's dimension)
    [shortest_distances_opposite_from_source_to_target,shortest_distances_opposite_from_target_to_source, slice_source, slice_target] = Function_shortest_distance_from_skeleton(G,node_coordinate,Domain_size,voxel_size,direction,Particle_label);
    
    % Save results
    Table_discrete_particle_graphnetowrk.direction(direction).shortest_distances_opposite_from_source_to_target=shortest_distances_opposite_from_source_to_target;
    Table_discrete_particle_graphnetowrk.direction(direction).shortest_distances_opposite_from_target_to_source=shortest_distances_opposite_from_target_to_source;
    Table_discrete_particle_graphnetowrk.direction(direction).slice_source=slice_source;
    Table_discrete_particle_graphnetowrk.direction(direction).slice_target=slice_target;
    
    % Asymetry between tortusosity calculated form source or target plane reveals the investigated volume is below the RVE
    
    % % Step 4: Analyse from plane (x=0 Source) to plane end (y=end target)
    
    % calculate the shortest distance distribution
    array=shortest_distances_opposite_from_source_to_target;
    array(array==Inf)=[];
    
    [tmp, ~] = Function_probability_density(array,[],density_fct_parameters);
    if ~isempty(tmp.smoothed_cumulative_fct)
        cumulative_function = tmp.smoothed_cumulative_fct;
        geometric_tortuosity_distribtuion = tmp.smoothed_probability_density_fct;
        geo_tor_integral = tmp.integral_probability_density_fct; if isempty(geo_tor_integral); geo_tor_integral= NaN; end
    else
        cumulative_function = tmp.cumulative_fct;
        geometric_tortuosity_distribtuion = tmp.probability_density_fct;
        geo_tor_integral = tmp.integral_smoothed_probability_density_fct; if isempty(geo_tor_integral); geo_tor_integral= NaN; end
    end
    %[cumulative_function, geometric_tortuosity_distribtuion, geo_tor_integral] = Function_size_distribution(array,length(array),10);
    %[cumulative_function, geometric_tortuosity_distribtuion, geo_tor_integral] = Function_size_distribution_updated(array,true,250,5);
    
    % Save it
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.cumulative_function=cumulative_function;
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.geometric_tortuosity_distribtuion=geometric_tortuosity_distribtuion;
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.geo_tor_integral=geo_tor_integral;
    
    % calculate statistics
    % t_50: the mean geometric tortuosity
    % 1D Interpolation  to get the d_50 from the cumulative function
    x = cumulative_function(:,1); y = cumulative_function(:,2);
    t_50 = interp1(y,x,0.5);
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.t_50 = t_50;
    % Min and max
    min_geotortuosity = min(geometric_tortuosity_distribtuion(:,1));
    max_geotortuosity = max(geometric_tortuosity_distribtuion(:,1));
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.min_geotortuosity = min_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.max_geotortuosity = max_geotortuosity;
    % The weighted starndard deviation formula is
    % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
    % With wi the weight of the xi, and <x> the weighted mean (mean_size)
    mean_geotortuosity = sum( geometric_tortuosity_distribtuion(:,1).*geometric_tortuosity_distribtuion(:,2) ) / sum (geometric_tortuosity_distribtuion(:,2));
    wi = geometric_tortuosity_distribtuion(:,2);
    xi = geometric_tortuosity_distribtuion(:,1);
    n = length(xi);
    std_geotortuosity = sqrt( sum( wi.*((xi-mean_geotortuosity).^2)) / ( (n-1)/n * sum(wi)  ));
    std_geotortuosity_percent = std_geotortuosity*100/mean_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.mean_geotortuosity = mean_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.std_geotortuosity = std_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).Plane0.std_geotortuosity_percent = std_geotortuosity_percent;
    
    % % Step 5: Analyse from plane end (target) to plane 0 (source)
    
    % calculate the shortest distance distribution
    array=shortest_distances_opposite_from_target_to_source;
    array(array==Inf)=[];
    
    [tmp, ~] = Function_probability_density(array,[],density_fct_parameters);
    if ~isempty(tmp.smoothed_cumulative_fct)
        cumulative_function = tmp.smoothed_cumulative_fct;
        geometric_tortuosity_distribtuion = tmp.smoothed_probability_density_fct;
        geo_tor_integral = tmp.integral_probability_density_fct; if isempty(geo_tor_integral); geo_tor_integral= NaN; end
    else
        cumulative_function = tmp.cumulative_fct;
        geometric_tortuosity_distribtuion = tmp.probability_density_fct;
        geo_tor_integral = tmp.integral_smoothed_probability_density_fct; if isempty(geo_tor_integral); geo_tor_integral= NaN; end
    end
    % [cumulative_function, geometric_tortuosity_distribtuion, geo_tor_integral] = Function_size_distribution(array,length(array),10);
    % Save it
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.cumulative_function=cumulative_function;
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.geometric_tortuosity_distribtuion=geometric_tortuosity_distribtuion;
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.geo_tor_integral=geo_tor_integral;
    
    % calculate statistics
    % t_50: the mean geometric tortuosity
    % 1D Interpolation  to get the d_50 from the cumulative function
    x = cumulative_function(:,1); y = cumulative_function(:,2);
    t_50 = interp1(y,x,0.5);
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.t_50 = t_50;
    % Min and max
    min_geotortuosity = min(geometric_tortuosity_distribtuion(:,1));
    max_geotortuosity = max(geometric_tortuosity_distribtuion(:,1));
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.min_geotortuosity = min_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.max_geotortuosity = max_geotortuosity;
    % The weighted starndard deviation formula is
    % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
    % With wi the weight of the xi, and <x> the weighted mean (mean_size)
    mean_geotortuosity = sum( geometric_tortuosity_distribtuion(:,1).*geometric_tortuosity_distribtuion(:,2) ) / sum (geometric_tortuosity_distribtuion(:,2));
    wi = geometric_tortuosity_distribtuion(:,2);
    xi = geometric_tortuosity_distribtuion(:,1);
    n = length(xi);
    std_geotortuosity = sqrt( sum( wi.*((xi-mean_geotortuosity).^2)) / ( (n-1)/n * sum(wi)  ));
    std_geotortuosity_percent = std_geotortuosity*100/mean_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.mean_geotortuosity = mean_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.std_geotortuosity = std_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.std_geotortuosity_percent = std_geotortuosity_percent;
    
    
    % % Step 6: Analyse all distances
    
    % Merge distances
    shortest_distances_opposite_faces=[shortest_distances_opposite_from_source_to_target; shortest_distances_opposite_from_target_to_source];
    array=shortest_distances_opposite_faces;
    array(array==Inf)=[];
    
    % calculate the shortest distance distribution
    [tmp, ~] = Function_probability_density(array,[],density_fct_parameters);
    if ~isempty(tmp.smoothed_cumulative_fct)
        cumulative_function = tmp.smoothed_cumulative_fct;
        geometric_tortuosity_distribtuion = tmp.smoothed_probability_density_fct;
        geo_tor_integral = tmp.integral_probability_density_fct; if isempty(geo_tor_integral); geo_tor_integral= NaN; end
    else
        cumulative_function = tmp.cumulative_fct;
        geometric_tortuosity_distribtuion = tmp.probability_density_fct;
        geo_tor_integral = tmp.integral_smoothed_probability_density_fct; if isempty(geo_tor_integral); geo_tor_integral= NaN; end
    end
    % [cumulative_function, geometric_tortuosity_distribtuion, geo_tor_integral] = Function_size_distribution(array,length(array),10);
    
    % Save it
    Table_discrete_particle_graphnetowrk.direction(direction).All.cumulative_function=cumulative_function;
    Table_discrete_particle_graphnetowrk.direction(direction).All.geometric_tortuosity_distribtuion=geometric_tortuosity_distribtuion;
    Table_discrete_particle_graphnetowrk.direction(direction).All.geo_tor_integral=geo_tor_integral;
    
    % calculate statistics
    % t_50: the mean geometric tortuosity
    % 1D Interpolation  to get the d_50 from the cumulative function
    x = cumulative_function(:,1); y = cumulative_function(:,2);
    t_50 = interp1(y,x,0.5);
    Table_discrete_particle_graphnetowrk.direction(direction).All.t_50 = t_50;
    % Min and max
    min_geotortuosity = min(geometric_tortuosity_distribtuion(:,1));
    max_geotortuosity = max(geometric_tortuosity_distribtuion(:,1));
    Table_discrete_particle_graphnetowrk.direction(direction).All.min_geotortuosity = min_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).All.max_geotortuosity = max_geotortuosity;
    % The weighted starndard deviation formula is
    % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
    % With wi the weight of the xi, and <x> the weighted mean (mean_size)
    mean_geotortuosity = sum( geometric_tortuosity_distribtuion(:,1).*geometric_tortuosity_distribtuion(:,2) ) / sum (geometric_tortuosity_distribtuion(:,2));
    wi = geometric_tortuosity_distribtuion(:,2);
    xi = geometric_tortuosity_distribtuion(:,1);
    n = length(xi);
    std_geotortuosity = sqrt( sum( wi.*((xi-mean_geotortuosity).^2)) / ( (n-1)/n * sum(wi)  ));
    std_geotortuosity_percent = std_geotortuosity*100/mean_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).All.mean_geotortuosity = mean_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).All.std_geotortuosity = std_geotortuosity;
    Table_discrete_particle_graphnetowrk.direction(direction).All.std_geotortuosity_percent = std_geotortuosity_percent;
    
    [str_direction] = function_remove_emptyandspecialcharacter_string(INFO.direction(direction).name);
    results_correlation(current_phase).(['Geometric_Tortuosity_mean_' dpsdname '_' num2str(direction) '_' str_direction]) = mean_geotortuosity;
    results_correlation(current_phase).(['Geometric_Tortuosity_std_' dpsdname '_' num2str(direction) '_' str_direction]) = std_geotortuosity;
    results_correlation(current_phase).(['Geometric_Tortuosity_relstd_' dpsdname '_' num2str(direction) '_' str_direction]) = std_geotortuosity_percent;
    
end


%% MANAGING RESULTS (GEO TORTUOSITY)
% Results are saved in a table

clear column
column(1).name='Sens_plus';
column(2).name='Sens_minus';
column(3).name='Both_sens';
for direction = 1:1:3
    Tau_50 = [Table_discrete_particle_graphnetowrk.direction(direction).Plane0.t_50; Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.t_50; Table_discrete_particle_graphnetowrk.direction(direction).All.t_50];
    Tau_min = [Table_discrete_particle_graphnetowrk.direction(direction).Plane0.min_geotortuosity; Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.min_geotortuosity; Table_discrete_particle_graphnetowrk.direction(direction).All.min_geotortuosity];
    Tau_max = [Table_discrete_particle_graphnetowrk.direction(direction).Plane0.max_geotortuosity; Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.max_geotortuosity; Table_discrete_particle_graphnetowrk.direction(direction).All.max_geotortuosity];
    Tau_std = [Table_discrete_particle_graphnetowrk.direction(direction).Plane0.std_geotortuosity; Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.std_geotortuosity; Table_discrete_particle_graphnetowrk.direction(direction).All.std_geotortuosity];
    Tau_stdpercent = [Table_discrete_particle_graphnetowrk.direction(direction).Plane0.std_geotortuosity_percent; Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.std_geotortuosity_percent; Table_discrete_particle_graphnetowrk.direction(direction).All.std_geotortuosity_percent];
    Tau_integral = [Table_discrete_particle_graphnetowrk.direction(direction).Plane0.geo_tor_integral; Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.geo_tor_integral; Table_discrete_particle_graphnetowrk.direction(direction).All.geo_tor_integral];
    Table_geometric_tortuosity_skeleton.direction(direction).table = table(char(column(:).name),Tau_50,Tau_min,Tau_max,Tau_std,Tau_stdpercent,Tau_integral,...
        'VariableNames',{'Sens' 'Tau_50' 'Min' 'Max' 'Std' 'Std_percents' 'Integral'});
    % Note: 'char' vertically concatenates character arrays, padding each input array as needed so that each row contains the same number of characters.
end
Table_discrete_particle_graphnetowrk.Table_geometric_tortuosity_skeleton = Table_geometric_tortuosity_skeleton;

%% DISPLAY TEXT RESULTS (GEO TORTUOSITY)
if OPTIONS.displaytext==true
    fprintf(['GEOMETRIC TORTUOSITY, skeleton (from ' dpsdname ') method:\n\n']);
    for direction = 1:1:3
        fprintf ('      For the %s\n',INFO.direction(direction).name);
        disp(Table_geometric_tortuosity_skeleton.direction(direction).table)
    end
end


%% SAVE RESULTS (GEO TORTUOSITY)
% % Text Results are saved
if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Geometric_tortuosity_skeleton_' dpsdname '_' INFO.phase(current_phase).name];
    % Prepare the data
    clear DATA_writetable;
    for direction = 1:1:3
        DATA_writetable.sheet(direction).name=['Geo_tort_dir_' num2str(direction)];
        DATA_writetable.sheet(direction).table=Table_geometric_tortuosity_skeleton.direction(direction).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


%% FIGURES (GEO TORTUOSITY)

% We will display the geometric tortuosity calculated for each voxel on the
% plane x=0 and at the plane x=end, depending on which particle they are
% belonging
% With the geometric tortuosity distribution

for direction=1:1:3
    % Direction name
    if direction ==1
        dir_name ='In plane dir. 1';
        dir_name_bis ='In_plane_direction_1';
        dir_Y = 'In plane dir. 2';
        dir_X = 'Through-plane dir';
    elseif direction ==2
        dir_name ='In plane dir. 2';
        dir_name_bis ='In_plane_direction_2';
        dir_Y = 'In plane dir. 1';
        dir_X = 'Through-plane dir';
    elseif direction ==3
        dir_name ='Through-plane dir.';
        dir_name_bis ='Through_plane_direction';
        dir_Y = 'In plane dir. 1';
        dir_X = 'In plane dir. 2';
    end    
    
    % % Create figure
    Fig_ = figure;
    Fig_.Name= ['Geometric_tortuosity_skeleton_dir' num2str(direction) '_' dpsdname '_' INFO.phase(current_phase).name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % Scale
    max_sensA = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.max_geotortuosity;
    max_sensB = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.max_geotortuosity;
    max_geo_tortuosity=max(max_sensA,max_sensB);
    
    % % GEO. TORTUOSITY MAP
    
    % Create axes as a subplot
    for id_axe=1:1:2
        
        % Sens
        if id_axe ==1
            sens_name ='Sens +';
            matrix_geometric_tortuosity = Table_discrete_particle_graphnetowrk.direction(direction).slice_source;
        elseif id_axe ==2
            sens_name ='Sens -';
            matrix_geometric_tortuosity = Table_discrete_particle_graphnetowrk.direction(direction).slice_target;
        end
        
        if direction==1
            if id_axe==1
                t_name = '\tau_{geo}^{In-plane dir.1 +}';
            elseif id_axe==2
                t_name = '\tau_{geo}^{In-plane dir.1 -}';
            end
        elseif direction==2
            if id_axe==1
                t_name = '\tau_{geo}^{In-plane dir.2 +}';
            elseif id_axe==2
                t_name = '\tau_{geo}^{In-plane dir.2 -}';
            end
        elseif direction==3
            if id_axe==1
                t_name = '\tau_{geo}^{Through-plane +}';
            elseif id_axe==2
                t_name = '\tau_{geo}^{Through-plane -}';
            end
        end
        
        % Subplot axe
        sub_axes{id_axe}=subplot(2,2,id_axe,'Parent',Fig_);
        % Active subplot
        hold(sub_axes{id_axe},'on');
        % - Title
        t_=title (' ','FontName','Times New Roman','FontSize',16);
        %t_.String= sprintf('Geo. Tortuosity map, %s, %s, %s',INFO.phase(current_phase).name,dir_name,sens_name);
        t_.String= sprintf('Geo. Tortuosity map, %s, %s',INFO.phase(current_phase).name,t_name);
        
        % Set NaN for voxels of the complementary phase: they will appear in white
        matrix_geometric_tortuosity(matrix_geometric_tortuosity==0)=NaN;
        matrix_geometric_tortuosity(matrix_geometric_tortuosity==Inf)=NaN;
        % Save the matrix for further plots
        matrix_geo_tortuosity.phase(current_phase).direction(direction).sens(id_axe).matrix = matrix_geometric_tortuosity;
        
        % Display geometric tortuosity map
        h_=pcolor(matrix_geometric_tortuosity);
        
        % Set properties
        set(h_,'LineStyle','none'); % No voxel line
        colormap(jet); % Jet color map
        set(sub_axes{id_axe},'CLim',[1 max_geo_tortuosity]); % color map scale
        colorbar('peer',sub_axes{id_axe}); % Display color bar
        
        % Axis label
        xlabel([sprintf('%s ',dir_X) '(\mum)']);
        ylabel([sprintf('%s ',dir_Y) '(\mum)']);
        
        % Fit the axes box
        %axis(sub_axes{id_axe},'tight');
        % Aspect ratio is 1:1
        axis(sub_axes{id_axe},'equal');
        
        % - Fontname and fontsize
        set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
        
        % Axis limit
        if direction==1
            ylim([0 Domain_size(2)]);
            xlim([0 Domain_size(3)]);
        elseif direction==2
            ylim([0 Domain_size(1)]);
            xlim([0 Domain_size(3)]);
        elseif direction==3
            ylim([0 Domain_size(1)]);
            xlim([0 Domain_size(2)]);
        end
        
        % Get Axis thick
        axis_xtick = sub_axes{id_axe}.XTick;
        axis_ytick = sub_axes{id_axe}.YTick;
        % Update Axis thick (micrometers)
        x_tick_newlabel=axis_xtick*voxel_size/1000;
        y_tick_newlabel=axis_ytick*voxel_size/1000;
        set(sub_axes{id_axe},'XTickLabel',x_tick_newlabel);
        set(sub_axes{id_axe},'YTickLabel',y_tick_newlabel);
        
        % Release current subplot
        hold(sub_axes{id_axe},'off');
        
    end
    
    
    % % CUMULATIVE FUNCTION
    
    % Create axes as a subplot
    id_axe=3;
    % Subplot axe
    sub_axes{id_axe}=subplot(2,2,id_axe,'Parent',Fig_);
    % Active subplot
    hold(sub_axes{id_axe},'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= sprintf('Geo. Tortuosity cumulative function, %s, %s',INFO.phase(current_phase).name,dir_name);
    % - Plot graphs
    % Curves
    cumulative_function = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.cumulative_function;
    h_sens_plus=plot(cumulative_function(:,1),cumulative_function(:,2));
    cumulative_function = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.cumulative_function;
    h_sens_minus=plot(cumulative_function(:,1),cumulative_function(:,2));
    cumulative_function = Table_discrete_particle_graphnetowrk.direction(direction).All.cumulative_function;
    h_sens_all=plot(cumulative_function(:,1),cumulative_function(:,2));
    % Colors
    set(h_sens_plus, 'Color', 'r');
    set(h_sens_minus, 'Color', 'b');
    set(h_sens_all, 'Color', 'k');
    % Axis limit
    x_lim = sub_axes{id_axe}.XLim;
    x_min = max(1,x_lim(1));
    % Add tau 50
    t_50 = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.t_50;
    l1 = plot([x_min t_50],[0.5 0.5]);
    l2 = plot([t_50 t_50],[0 0.5]);
    set(l1, 'Color', 'r','MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    set(l2, 'Color', 'r','MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    t_50 = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.t_50;
    l1 = plot([x_min t_50],[0.5 0.5]);
    l2 = plot([t_50 t_50],[0 0.5]);
    set(l1, 'Color', 'b','MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    set(l2, 'Color', 'b','MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    t_50 = Table_discrete_particle_graphnetowrk.direction(direction).All.t_50;
    l1 = plot([x_min t_50],[0.5 0.5]);
    l2 = plot([t_50 t_50],[0 0.5]);
    set(l1, 'Color', 'k','MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    set(l2, 'Color', 'k','MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    % - Legend
    if direction==1
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.1 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.1}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif direction==2
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.2 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif direction==3
        legend(sub_axes{id_axe},['\tau_{geo}^{Through-plane dir. +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir. -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir.}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    end
    %legend(sub_axes{id_axe},'show','Location','best','Sens +','Sens -','Both');
    % - Axis label
    %if direction==1
    %    xlabel('Geometric tortuosity \tau_{geo}^{In-plane dir.1}');
    %elseif direction==2
    %    xlabel('Geometric tortuosity \tau_{geo}^{In-plane dir.2}');
    %elseif direction==3
    %    xlabel('Geometric tortuosity \tau_{geo}^{Through-plane}');
    %end
    xlabel('Geometric tortuosity \tau_{geo}');
    ylabel('Normalized cumulative function');
    % Axis min
    sub_axes{id_axe}.XLim = [x_min x_lim(2)];
    sub_axes{id_axe}.YLim = [0 1];
    % - Grid
    grid(sub_axes{id_axe},'on'); % Display grid
    set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
    % Release current subplot
    hold(sub_axes{id_axe},'off');
    
    % % GEOMETRIC TORTUOSITY DISTRIBUTION
    
    % Create axes as a subplot
    id_axe=4;
    % Subplot axe
    sub_axes{id_axe}=subplot(2,2,id_axe,'Parent',Fig_);
    % Active subplot
    hold(sub_axes{id_axe},'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= sprintf('Geo. Tortuosity distribution, %s, %s',INFO.phase(current_phase).name,dir_name);
    % - Plot graphs
    % Curves
    geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.geometric_tortuosity_distribtuion;
    h_sens_plus=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
    geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.geometric_tortuosity_distribtuion;
    h_sens_minus=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
    geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(direction).All.geometric_tortuosity_distribtuion;
    h_sens_all=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
    % Colors
    set(h_sens_plus, 'Color', 'r');
    set(h_sens_minus, 'Color', 'b');
    set(h_sens_all, 'Color', 'k');
    % Axis limit
    x_lim = sub_axes{id_axe}.XLim;
    x_min = max(1,x_lim(1));
    % - Legend
    if direction==1
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.1 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.1}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif direction==2
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.2 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif direction==3
        legend(sub_axes{id_axe},['\tau_{geo}^{Through-plane dir. +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir. -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir.}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    end
    %legend(sub_axes{id_axe},'show','Location','best','Sens +','Sens -','Both');
    % - Axis label
    % if direction==1
    %    xlabel('Geometric tortuosity \tau_{geo}^{In-plane dir.1}');
    % elseif direction==2
    %    xlabel('Geometric tortuosity \tau_{geo}^{In-plane dir.2}');
    % elseif direction==3
    %    xlabel('Geometric tortuosity \tau_{geo}^{Through-plane}');
    %end
    xlabel('Geometric tortuosity \tau_{geo}');
    ylabel('Distribution');
    % Axis min
    sub_axes{id_axe}.XLim = [x_min x_lim(2)];
    % - Grid
    grid(sub_axes{id_axe},'on'); % Display grid
    set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
    % Release current subplot
    hold(sub_axes{id_axe},'off');
    
    if OPTIONS.save_fig == true % Save figure
        filename= Fig_.Name;
        function_savefig(Fig_, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig_); % Do not keep open figures
    end
end




%% FIGURES (GEO TORTUOSITY)

% The cumulative function and geometric tortuosity distribution

% % Create figure
Fig_ = figure;
Fig_.Name= ['Geometric_tortuosity_skeleton_' dpsdname '_' INFO.phase(current_phase).name];
Fig_.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig_,'position',scrsz); % Full screen figure

% % CUMULATIVE FUNCTION
for id_axe=1:1:3
    % Create axes as a subplot
    % Subplot axe
    sub_axes{id_axe}=subplot(2,3,id_axe,'Parent',Fig_);
    % Active subplot
    hold(sub_axes{id_axe},'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    if id_axe==1
        t_.String= sprintf('Geo. Tortuosity cumulative function, %s, Sens +',INFO.phase(current_phase).name);
    elseif id_axe==2
        t_.String= sprintf('Geo. Tortuosity cumulative function, %s, Sens -',INFO.phase(current_phase).name);
    elseif id_axe==3
        t_.String= sprintf('Geo. Tortuosity cumulative function, %s, Both',INFO.phase(current_phase).name);
    end
    % - Plot graphs
    % Curves
    if id_axe==1
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(1).Plane0.cumulative_function;
        h_direction1=plot(cumulative_function(:,1),cumulative_function(:,2));
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(2).Plane0.cumulative_function;
        h_direction2=plot(cumulative_function(:,1),cumulative_function(:,2));
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(3).Plane0.cumulative_function;
        h_direction3=plot(cumulative_function(:,1),cumulative_function(:,2));
    elseif id_axe==2
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(1).PlaneEnd.cumulative_function;
        h_direction1=plot(cumulative_function(:,1),cumulative_function(:,2));
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(2).PlaneEnd.cumulative_function;
        h_direction2=plot(cumulative_function(:,1),cumulative_function(:,2));
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(3).PlaneEnd.cumulative_function;
        h_direction3=plot(cumulative_function(:,1),cumulative_function(:,2));
    elseif id_axe==3
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(1).All.cumulative_function;
        h_direction1=plot(cumulative_function(:,1),cumulative_function(:,2));
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(2).All.cumulative_function;
        h_direction2=plot(cumulative_function(:,1),cumulative_function(:,2));
        cumulative_function = Table_discrete_particle_graphnetowrk.direction(3).All.cumulative_function;
        h_direction3=plot(cumulative_function(:,1),cumulative_function(:,2));
    end
    % Colors
    set(h_direction1, 'Color', [0 114 189]/255);
    set(h_direction2, 'Color', [217 83 25]/255);
    set(h_direction3, 'Color', [237 177 32]/255);
    % Axis limit
    x_lim = sub_axes{id_axe}.XLim;
    x_min = max(1,x_lim(1));
    % Add tau 50
    if id_axe==1
        t_50 = Table_discrete_particle_graphnetowrk.direction(1).Plane0.t_50;
    elseif id_axe==2
        t_50 = Table_discrete_particle_graphnetowrk.direction(1).PlaneEnd.t_50;
    elseif id_axe==3
        t_50 = Table_discrete_particle_graphnetowrk.direction(1).All.t_50;
    end
    l1 = plot([x_min t_50],[0.5 0.5]);
    l2 = plot([t_50 t_50],[0 0.5]);
    set(l1, 'Color', [0 114 189]/255,'MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    set(l2, 'Color', [0 114 189]/255,'MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    if id_axe==1
        t_50 = Table_discrete_particle_graphnetowrk.direction(2).Plane0.t_50;
    elseif id_axe==2
        t_50 = Table_discrete_particle_graphnetowrk.direction(2).PlaneEnd.t_50;
    elseif id_axe==3
        t_50 = Table_discrete_particle_graphnetowrk.direction(2).All.t_50;
    end
    l1 = plot([x_min t_50],[0.5 0.5]);
    l2 = plot([t_50 t_50],[0 0.5]);
    set(l1, 'Color', [217 83 25]/255,'MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    set(l2, 'Color', [217 83 25]/255,'MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    if id_axe==1
        t_50 = Table_discrete_particle_graphnetowrk.direction(3).Plane0.t_50;
    elseif id_axe==2
        t_50 = Table_discrete_particle_graphnetowrk.direction(3).PlaneEnd.t_50;
    elseif id_axe==3
        t_50 = Table_discrete_particle_graphnetowrk.direction(3).All.t_50;
    end
    t_50 = Table_discrete_particle_graphnetowrk.direction(3).All.t_50;
    l1 = plot([x_min t_50],[0.5 0.5]);
    l2 = plot([t_50 t_50],[0 0.5]);
    set(l1, 'Color', [237 177 32]/255,'MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    set(l2, 'Color', [237 177 32]/255,'MarkerSize',12,'Marker','none','LineWidth',1,'LineStyle','--');
    % - Legend
    if id_axe==1
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir. +}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif id_axe==2
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir. -}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif id_axe==3
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir.}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    end
    % legend(sub_axes{id_axe},'show','Location','best','In plane dir. 1','In plane dir. 2','Through-plane dir.');
    % - Axis label
    if id_axe==1
        xlabel(['Geometric tortuosity ' '\tau_{geo}^{+}']);
    elseif id_axe==2
        xlabel(['Geometric tortuosity ' '\tau_{geo}^{-}']);
    elseif id_axe==3
        xlabel(['Geometric tortuosity ' '\tau_{geo}']);
    end
    ylabel('Normalized cumulative function');
    % Axis min
    sub_axes{id_axe}.XLim = [x_min x_lim(2)];
    sub_axes{id_axe}.YLim = [0 1];
    % - Grid
    grid(sub_axes{id_axe},'on'); % Display grid
    set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
    % Release current subplot
    hold(sub_axes{id_axe},'off');
end

% % GEOMETRIC TORTUOSITY DISTRIBUTION
for id_axe=4:1:6
    % Create axes as a subplot
    % Subplot axe
    sub_axes{id_axe}=subplot(2,3,id_axe,'Parent',Fig_);
    % Active subplot
    hold(sub_axes{id_axe},'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    if id_axe==4
        t_.String= sprintf('Geo. Tortuosity distribution, %s, Sens +',INFO.phase(current_phase).name);
    elseif id_axe==5
        t_.String= sprintf('Geo. Tortuosity distribution, %s, Sens -',INFO.phase(current_phase).name);
    elseif id_axe==6
        t_.String= sprintf('Geo. Tortuosity distribution, %s, Both',INFO.phase(current_phase).name);
    end
    % - Plot graphs
    % Curves
    if id_axe==4
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(1).Plane0.geometric_tortuosity_distribtuion;
        h_direction1=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(2).Plane0.geometric_tortuosity_distribtuion;
        h_direction2=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(3).Plane0.geometric_tortuosity_distribtuion;
        h_direction3=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
    elseif id_axe==5
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(1).PlaneEnd.geometric_tortuosity_distribtuion;
        h_direction1=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(2).PlaneEnd.geometric_tortuosity_distribtuion;
        h_direction2=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(3).PlaneEnd.geometric_tortuosity_distribtuion;
        h_direction3=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
    elseif id_axe==6
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(1).All.geometric_tortuosity_distribtuion;
        h_direction1=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(2).All.geometric_tortuosity_distribtuion;
        h_direction2=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
        geometric_tortuosity_distribtuion = Table_discrete_particle_graphnetowrk.direction(3).All.geometric_tortuosity_distribtuion;
        h_direction3=plot(geometric_tortuosity_distribtuion(:,1),geometric_tortuosity_distribtuion(:,2));
    end
    % Colors
    set(h_direction1, 'Color', [0 114 189]/255);
    set(h_direction2, 'Color', [217 83 25]/255);
    set(h_direction3, 'Color', [237 177 32]/255);
    % Axis limit
    x_lim = sub_axes{id_axe}.XLim;
    x_min = max(1,x_lim(1));
    % - Legend
    if id_axe==4
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2 +}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir. +}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif id_axe==5
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2 -}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir. -}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    elseif id_axe==6
        legend(sub_axes{id_axe},['\tau_{geo}^{In-plane dir.1}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{In-plane dir.2}' ' (' INFO.phase(current_phase).name ')'],['\tau_{geo}^{Through-plane dir.}' ' (' INFO.phase(current_phase).name ')'],'Location','best');
    end
    %legend(sub_axes{id_axe},'show','Location','best','In plane dir. 1','In plane dir. 2','Through-plane dir.');
    % - Axis label
    if id_axe==4
        xlabel(['Geometric tortuosity ' '\tau_{geo}^{+}']);
    elseif id_axe==5
        xlabel(['Geometric tortuosity ' '\tau_{geo}^{-}']);
    elseif id_axe==6
        xlabel(['Geometric tortuosity ' '\tau_{geo}']);
    end
    ylabel('Distribution');
    % Axis min
    sub_axes{id_axe}.XLim = [x_min x_lim(2)];
    % - Grid
    grid(sub_axes{id_axe},'on'); % Display grid
    set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
    % Release current subplot
    hold(sub_axes{id_axe},'off');
end

if OPTIONS.save_fig == true % Save figure
    filename= ['Geometric_tortuosity_skeleton2_' dpsdname '_' INFO.phase(current_phase).name];
    function_savefig(Fig_, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig_); % Do not keep open figures
end


%% EVOLUTION OF THE GEOMETRIC TORTUOSITY

% WARNING: this is a not a classic 'plot along direction' as done many
% times for the others properties, such as volume fraction.
% Here for the geo. tortuosity calculated for the direction 1, for the sens
% plus, we have a variation along the two other directions: this variation
% will be ploted.

% The geometric map is displayed again, for a better understanding

for direction=1:1:3
    
    % Direction name
    if direction ==1
        dir_name ='In plane dir. 1';
        dir_name_bis ='In_plane_direction_1';
        dir_Y = 'In plane dir. 2';
        dir_X = 'Through-plane dir';
    elseif direction ==2
        dir_name ='In plane dir. 2';
        dir_name_bis ='In_plane_direction_2';
        dir_Y = 'In plane dir. 1';
        dir_X = 'Through-plane dir';
    elseif direction ==3
        dir_name ='Through-plane dir.';
        dir_name_bis ='Through_plane_direction';
        dir_Y = 'In plane dir. 1';
        dir_X = 'In plane dir. 2';
    end
    
    % Scale
    max_sensA = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.max_geotortuosity;
    max_sensB = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.max_geotortuosity;
    max_geo_tortuosity=max(max_sensA,max_sensB);
    
    % % Create figure
    Fig_ = figure;
    Fig_.Name= ['Geometric_tortuosity_skeletonWatershed_inplane_evolution_' INFO.phase(current_phase).name '_' dir_name_bis];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    
    % % GEO. TORTUOSITY MAP
    
    % Create axes as a subplot
    for id_axe=1:3:4
        
        % Sens
        if id_axe ==1
            sens_name ='Sens +';
            map_ = matrix_geo_tortuosity.phase(current_phase).direction(direction).sens(1).matrix;
        elseif id_axe ==4
            sens_name ='Sens -';
            map_ = matrix_geo_tortuosity.phase(current_phase).direction(direction).sens(2).matrix;
        end
        
        if direction==1
            if id_axe==1
                t_name = '\tau_{geo}^{In-plane dir.1 +}';
            elseif id_axe==4
                t_name = '\tau_{geo}^{In-plane dir.1 -}';
            end
        elseif direction==2
            if id_axe==1
                t_name = '\tau_{geo}^{In-plane dir.2 +}';
            elseif id_axe==4
                t_name = '\tau_{geo}^{In-plane dir.2 -}';
            end
        elseif direction==3
            if id_axe==1
                t_name = '\tau_{geo}^{Through-plane +}';
            elseif id_axe==4
                t_name = '\tau_{geo}^{Through-plane -}';
            end
        end
        
        % Subplot axe
        sub_axes{id_axe}=subplot(2,3,id_axe,'Parent',Fig_);
        % Active subplot
        hold(sub_axes{id_axe},'on');
        % - Title
        t_=title (' ','FontName','Times New Roman','FontSize',16);
        %t_.String= sprintf('Geo. Tortuosity map, %s, %s, %s',INFO.phase(current_phase).name,dir_name,sens_name);
        t_.String= sprintf('Geo. Tortuosity map, %s, %s',INFO.phase(current_phase).name,t_name);
        
        % Display geometric tortuosity map
        h_=pcolor(map_);
        
        % Set properties
        set(h_,'LineStyle','none'); % No voxel line
        colormap(jet); % Jet color map
        set(sub_axes{id_axe},'CLim',[1 max_geo_tortuosity]); % color map scale
        colorbar('peer',sub_axes{id_axe}); % Display color bar
        
        % Axis label
        xlabel([sprintf('%s ',dir_X) '(\mum)']);
        ylabel([sprintf('%s ',dir_Y) '(\mum)']);
        
        % Fit the axes box
        %axis(sub_axes{id_axe},'tight');
        % Aspect ratio is 1:1
        axis(sub_axes{id_axe},'equal');
        
        % - Fontname and fontsize
        set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
        
        % Axis limit
        if direction==1
            ylim([0 Domain_size(2)]);
            xlim([0 Domain_size(3)]);
        elseif direction==2
            ylim([0 Domain_size(1)]);
            xlim([0 Domain_size(3)]);
        elseif direction==3
            ylim([0 Domain_size(1)]);
            xlim([0 Domain_size(2)]);
        end
        
        % Get Axis thick
        axis_xtick = sub_axes{id_axe}.XTick;
        axis_ytick = sub_axes{id_axe}.YTick;
        % Update Axis thick (micrometers)
        x_tick_newlabel=axis_xtick*voxel_size/1000;
        y_tick_newlabel=axis_ytick*voxel_size/1000;
        set(sub_axes{id_axe},'XTickLabel',x_tick_newlabel);
        set(sub_axes{id_axe},'YTickLabel',y_tick_newlabel);
        
        % Release current subplot
        hold(sub_axes{id_axe},'off');
        
    end
    
    % % Evolution along direction perpendicular to the direction
    
    % Create axes as a subplot
    for id_axe=1:1:6
        if (id_axe==2 || id_axe==3 || id_axe==5 || id_axe==6)
            if (id_axe==2 || id_axe==3)
                % Map
                Map_ = matrix_geo_tortuosity.phase(current_phase).direction(direction).sens(1).matrix;
            elseif (id_axe==5 || id_axe==6)
                % Map
                Map_ = matrix_geo_tortuosity.phase(current_phase).direction(direction).sens(2).matrix;
            end
            [n_dir_X, n_dir_Y]=size(Map_);
            
            % Set identical min-max
            min_1 = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.min_geotortuosity;
            min_2 = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.min_geotortuosity;
            max_1 = Table_discrete_particle_graphnetowrk.direction(direction).Plane0.max_geotortuosity;
            max_2 = Table_discrete_particle_graphnetowrk.direction(direction).PlaneEnd.max_geotortuosity;
            y_minaxis=min(min_1,min_2);
            y_maxaxis=max(max_1,max_2);
            % Set a margin, for a better visualisation
            margin = ((y_maxaxis-y_minaxis)*0.20)/2;
            y_minaxis=y_minaxis-margin;
            y_maxaxis=y_maxaxis+margin;
            
            % Initialization
            % :,1 position in micrometers
            % :,2 mean diameter
            % :,3 max diameter
            % :,4 std diameter
            % :,5 min diameter
            size_direction_X=zeros(n_dir_X,5);
            size_direction_Y=zeros(n_dir_Y,5);
            % Position axis in micrometers
            size_direction_X(:,1)=(1:1:n_dir_X)*voxel_size/1000;
            size_direction_Y(:,1)=(1:1:n_dir_Y)*voxel_size/1000;
            
            if (id_axe==2 || id_axe==5)
                if direction==1
                    if id_axe==2
                        t_name = '\tau_{geo}^{In-plane dir.1 +}';
                    elseif id_axe==5
                        t_name = '\tau_{geo}^{In-plane dir.1 -}';
                    end
                elseif direction==2
                    if id_axe==2
                        t_name = '\tau_{geo}^{In-plane dir.2 +}';
                    elseif id_axe==5
                        t_name = '\tau_{geo}^{In-plane dir.2 -}';
                    end
                elseif direction==3
                    if id_axe==2
                        t_name = '\tau_{geo}^{Through-plane +}';
                    elseif id_axe==5
                        t_name = '\tau_{geo}^{Through-plane -}';
                    end
                end
                % Direction X
                for position=1:1:n_dir_X
                    % Current slice of data
                    slice_ = Map_(position,:);
                    % Get all values
                    different_values =  unique(slice_);
                    % Remove NaN
                    different_values(isnan(different_values))=[];
                    if ~isempty(different_values)
                        % Maximum and minimum value
                        size_direction_X(position,3) = max(different_values);
                        size_direction_X(position,5) = min(different_values);
                        % Weighed values
                        % Initialisation
                        values=zeros(length(different_values),2);
                        for current_size=1:1:length(different_values)
                            values(current_size,1)=different_values(current_size);
                            values(current_size,2)=sum(sum( slice_== different_values(current_size)));
                        end
                        % Mean diameter
                        size_direction_X(position,2)=sum(values(:,1).*values(:,2))/sum(values(:,2));
                        % Standard deviation
                        % The weighted starndard deviation formula is
                        % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
                        % With wi the weight of the xi, and <x> the weighted mean (mean_size)
                        wi = values(:,2);
                        xi = values(:,1);
                        n = length(xi);
                        mean_size = size_direction_X(position,2);
                        size_direction_X(position,4) = sqrt( sum( wi.*((xi-mean_size).^2)) / ( (n-1)/n * sum(wi)  ));
                    else
                        size_direction_X(position,:)=NaN;
                    end
                end
                x_=size_direction_X(:,1);
                y_mean = size_direction_X(:,2);
                y_max = size_direction_X(:,3);
                y_std = size_direction_X(:,4);
                y_min = size_direction_X(:,5);
                % Subplot axe
                sub_axes{id_axe}=subplot(2,3,id_axe,'Parent',Fig_);
                % Active subplot
                hold(sub_axes{id_axe},'on');
                % - Title
                t_=title (' ','FontName','Times New Roman','FontSize',16);
                t_.String= sprintf('%s, %s, along %s',t_name,INFO.phase(current_phase).name,dir_Y);
                % Mean
                h_mean=plot(x_,y_mean); % For the legend order
                % Extremums
                h_max=plot(x_,y_max);
                h_min=plot(x_,y_min);
                % Colors, thickness, markers
                set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
                set(h_max, 'Color', INFO.phase(current_phase).color,'LineStyle','--','MarkerSize',12,'Marker','none');
                set(h_min, 'Color', INFO.phase(current_phase).color,'LineStyle','--','MarkerSize',12,'Marker','none');
                % Mean with error bar (+- standard deviation)
                h_mean_witherrorbar = errorbar(x_,y_mean,y_std);
                set(h_mean_witherrorbar, 'Color', INFO.phase(current_phase).color,'LineWidth',1,'MarkerSize',12,'Marker','none');
                h_mean=plot(x_,y_mean); % Plot over the other
                set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
                % - Axis label
                t_ = xlabel(' ');
                t_1 = sprintf('Position along %s ',dir_Y);
                t_2 = '(\mum)';
                t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
                t_ = ylabel(' ');
                t_.String= [t_name]; % Sprintf does not accept greek characters
                % - Legend
                legend(sub_axes{id_axe},'Mean \tau_{geo} (with std)','Min. & Max. \tau_{geo}','Location','best');
                % - Grid
                grid(sub_axes{id_axe},'on'); % Display grid
                set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                % - Fontname and fontsize
                set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
                % Axis limit
                if direction==1
                    xlim([0 Domain_size(2)*voxel_size/1000]);
                elseif direction==2
                    xlim([0 Domain_size(1)*voxel_size/1000]);
                elseif direction==3
                    xlim([0 Domain_size(1)*voxel_size/1000]);
                end
                ylim([y_minaxis y_maxaxis]);
                % - Figure has been done
                hold(sub_axes{id_axe},'off');
                
            end
            
            if (id_axe==3 || id_axe==6)
                if direction==1
                    if id_axe==3
                        t_name = '\tau_{geo}^{In-plane dir.1 +}';
                    elseif id_axe==6
                        t_name = '\tau_{geo}^{In-plane dir.1 -}';
                    end
                elseif direction==2
                    if id_axe==3
                        t_name = '\tau_{geo}^{In-plane dir.2 +}';
                    elseif id_axe==6
                        t_name = '\tau_{geo}^{In-plane dir.2 -}';
                    end
                elseif direction==3
                    if id_axe==3
                        t_name = '\tau_{geo}^{Through-plane +}';
                    elseif id_axe==6
                        t_name = '\tau_{geo}^{Through-plane -}';
                    end
                end
                % Direction Y
                for position=1:1:n_dir_Y
                    % Current slice of data
                    slice_ = Map_(:,position);
                    % Get all values
                    different_values =  unique(slice_);
                    % Remove NaN
                    different_values(isnan(different_values))=[];
                    if ~isempty(different_values)
                        % Maximum and minimum value
                        size_direction_Y(position,3) = max(different_values);
                        size_direction_Y(position,5) = min(different_values);
                        % Weighed values
                        % Initialisation
                        values=zeros(length(different_values),2);
                        for current_size=1:1:length(different_values)
                            values(current_size,1)=different_values(current_size);
                            values(current_size,2)=sum(sum( slice_== different_values(current_size)));
                        end
                        % Mean diameter
                        size_direction_Y(position,2)=sum(values(:,1).*values(:,2))/sum(values(:,2));
                        % Standard deviation
                        % The weighted starndard deviation formula is
                        % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
                        % With wi the weight of the xi, and <x> the weighted mean (mean_size)
                        wi = values(:,2);
                        xi = values(:,1);
                        n = length(xi);
                        mean_size = size_direction_Y(position,2);
                        size_direction_Y(position,4) = sqrt( sum( wi.*((xi-mean_size).^2)) / ( (n-1)/n * sum(wi)  ));
                    else
                        size_direction_Y(position,:)=NaN;
                    end
                end
                x_=size_direction_Y(:,1);
                y_mean = size_direction_Y(:,2);
                y_max = size_direction_Y(:,3);
                y_std = size_direction_Y(:,4);
                y_min = size_direction_Y(:,5);
                % Subplot axe
                sub_axes{id_axe}=subplot(2,3,id_axe,'Parent',Fig_);
                % Active subplot
                hold(sub_axes{id_axe},'on');
                % - Title
                t_=title (' ','FontName','Times New Roman','FontSize',16);
                t_.String= sprintf('%s, %s, along %s',t_name,INFO.phase(current_phase).name,dir_X);
                % Mean
                h_mean=plot(x_,y_mean); % For the legend order
                % Extremums
                h_max=plot(x_,y_max);
                h_min=plot(x_,y_min);
                % Colors, thickness, markers
                set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
                set(h_max, 'Color', INFO.phase(current_phase).color,'LineStyle','--','MarkerSize',12,'Marker','none');
                set(h_min, 'Color', INFO.phase(current_phase).color,'LineStyle','--','MarkerSize',12,'Marker','none');
                % Mean with error bar (+- standard deviation)
                h_mean_witherrorbar = errorbar(x_,y_mean,y_std);
                set(h_mean_witherrorbar, 'Color', INFO.phase(current_phase).color,'LineWidth',1,'MarkerSize',12,'Marker','none');
                h_mean=plot(x_,y_mean); % Plot over the other
                set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
                % - Axis label
                t_ = xlabel(' ');
                t_1 = sprintf('Position along %s ',dir_X);
                t_2 = '(\mum)';
                t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
                t_ = ylabel(' ');
                t_.String= [t_name]; % Sprintf does not accept greek characters
                % - Legend
                legend(sub_axes{id_axe},'Mean \tau_{geo}','Min. & Max. \tau_{geo}','Location','best');
                % - Grid
                grid(sub_axes{id_axe},'on'); % Display grid
                set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                % - Fontname and fontsize
                set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
                % Axis limit
                if direction==1
                    xlim([0 Domain_size(3)*voxel_size/1000]);
                elseif direction==2
                    xlim([0 Domain_size(3)*voxel_size/1000]);
                elseif direction==3
                    xlim([0 Domain_size(2)*voxel_size/1000]);
                end
                ylim([y_minaxis y_maxaxis]);
                % - Figure has been done
                hold(sub_axes{id_axe},'off');
            end
            
            
        end
    end
    
    if OPTIONS.save_fig == true % Save figure
        filename= ['Geometric_tortuosity_skeleton_' dpsdname '_' INFO.phase(current_phase).name 'inplane_evolution_' dir_name_bis];
        function_savefig(Fig_, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig_); % Do not keep open figures
    end
    
end





end

