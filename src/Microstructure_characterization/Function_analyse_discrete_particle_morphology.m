function [results_correlation, results_visualization, Table_discrete_particlemorpholgy_analysis, unique_particle, mass_center_particle] = Function_analyse_discrete_particle_morphology(Particle_size, Particle_label, results_correlation, results_visualization, current_phase, Current_folder, INFO, OPTIONS, dpsdname)

voxel_size = INFO.voxel_size;
% Unique id
unique_particle = unique(Particle_label);
unique_particle(unique_particle==0)=[]; % Remove background
% Number of different particle
number_particle = length(unique_particle);
% Domain size
Domain_size = size(Particle_size);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end
if OPTIONS.displaytext==true
    fprintf ('- Number of particle: %i\n',number_particle)
end

% Save
Table_discrete_particlemorpholgy_analysis.number_particle=number_particle;
Variable_name_table={'Number_particle'};
Table_tmp = array2table(number_particle,...
    'VariableNames',Variable_name_table);

%%
%% MASS CENTER
%%
%%

% Mass center initialisation
% Id, position x,y,z, number of voxel, equivalent diameter in micrometer
mass_center_particle = zeros(number_particle,6);
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    mass_center_particle(current_,1)=particle_id;
    % Get index
    index_particle = find(Particle_label==particle_id);
    % Get all coordinates
    [M1,M2,M3] = ind2sub(Domain_size,index_particle);
    % Center of voxel position
    M1=M1-0.5; M2=M2-0.5; M3=M3-0.5;
    % Calculate mass center position
    x_m = mean(M1); y_m = mean(M2); z_m = mean(M3);
    mass_center_particle(current_,2)=x_m*voxel_size/1000;
    mass_center_particle(current_,3)=y_m*voxel_size/1000;
    mass_center_particle(current_,4)=z_m*voxel_size/1000;
    % Number of voxel
    volume_=length(M1); % 3D case
    area_ = length(M1); % 2D case
    mass_center_particle(current_,5)=volume_;
    % The equivalent diameter is
    if Domain_size(3)>1
        % 3D case
        mass_center_particle(current_,6)= 2 * ((3*volume_  /(4*pi))^(1/3));
    else
        % 2D case
        mass_center_particle(current_,6)= 2 * ((area_/pi)^(1/2));
    end
    % Convert in micrometer
    mass_center_particle(current_,6) = mass_center_particle(current_,6)*voxel_size/1000;
end

% Save
Variable_name_table={'Particle_Id' 'Position_x' 'Position_y' 'Position_z' 'Number_voxel' 'Equivalent_diameter'  };
Table_mass_center_particle = array2table(mass_center_particle,...
    'VariableNames',Variable_name_table);
Table_discrete_particlemorpholgy_analysis.Table_mass_center_particle = Table_mass_center_particle;
% % Text Results are saved
if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Particle_Centroid_' dpsdname '_' INFO.phase(current_phase).name];
    full_path=[Current_folder filename];
    % Prepare the data
    clear DATA_writetable;
    % Data : Number of particle
    DATA_writetable.sheet(1).name='Number_particle';
    DATA_writetable.sheet(1).table=Table_tmp;
    % Data : Mass center
    DATA_writetable.sheet(2).name='Mass_center';
    DATA_writetable.sheet(2).table=Table_mass_center_particle;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


%%
%% PARTICLE SIZE RATIO
%%
%%

%% CALCULATION
% Particle ratio initialisation
% Id, Equivalent diameter, ratio xy, ratio xz, ratio yz
particle_ratio = zeros(number_particle,5);
% Voxel per voxel
particle_ratio_xy_voxel = zeros(Domain_size);
particle_ratio_xz_voxel = zeros(Domain_size);
particle_ratio_yz_voxel = zeros(Domain_size);
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    particle_ratio(current_,1)=particle_id;
    % Equivalent diameter
    particle_ratio(current_,2)=mass_center_particle(current_,6);
    % Get index
    index_particle = find(Particle_label==particle_id);
    % Get all coordinates
    [M1,M2,M3] = ind2sub(Domain_size,index_particle);
    % Get minimum and maximum
    x_min = min(M1); y_min = min(M2); z_min = min(M3); 
    x_max = max(M1); y_max = max(M2); z_max = max(M3); 
    % Delta in voxel length
    delta_x = x_max-x_min+1;
    delta_y = y_max-y_min+1;
    delta_z = z_max-z_min+1;
    % Ratio xy
    current_ratio_xy = delta_x/delta_y;
    particle_ratio(current_,3) = current_ratio_xy;
    particle_ratio_xy_voxel(index_particle) = current_ratio_xy;
    % Ratio xz
    current_ratio_xz = delta_x/delta_z;    
    particle_ratio(current_,4) = current_ratio_xz;
    particle_ratio_xz_voxel(index_particle) = current_ratio_xz;
    % Ratio yz
    current_ratio_yz = delta_y/delta_z;
    particle_ratio(current_,5) = current_ratio_yz;
    particle_ratio_yz_voxel(index_particle) = current_ratio_yz;
end
% Save
Variable_name_table={'Particle_Id' 'Equivalent_diameter' 'Ratio_xy' 'Ratio_xz' 'Ratio_yz'};
Table_particle_ratio = array2table(particle_ratio,...
    'VariableNames',Variable_name_table);
Table_discrete_particlemorpholgy_analysis.Table_particle_ratio = Table_particle_ratio;

results_visualization(current_phase).(['particle_ratio_xy_' dpsdname]) = particle_ratio_xy_voxel;
results_visualization(current_phase).(['particle_ratio_xz_' dpsdname]) = particle_ratio_xz_voxel;
results_visualization(current_phase).(['particle_ratio_yz_' dpsdname]) = particle_ratio_yz_voxel;
% if strcmp(dpsdname,'watershed')
%     results_visualization(current_phase).particle_ratio_xy_watershed = particle_ratio_xy_voxel;
%     results_visualization(current_phase).particle_ratio_xz_watershed = particle_ratio_xz_voxel;
%     results_visualization(current_phase).particle_ratio_yz_watershed = particle_ratio_yz_voxel;
% elseif strcmp(dpsdname,'watershed')
%     results_visualization(current_phase).particle_ratio_xy_PCRF = particle_ratio_xy_voxel;
%     results_visualization(current_phase).particle_ratio_xz_PCRF = particle_ratio_xz_voxel;
%     results_visualization(current_phase).particle_ratio_yz_PCRF = particle_ratio_yz_voxel;    
% end
    

%% STATISTICS
% Statistics on particle ratio
% Particle ratio statistics initialisation
% Row 1: min
% Row 2: mean, weighted by the particle volume
% Row 3: max
% Row 4: std
% Row 5: std % of the mean
% Column 1: Ratio xy
% Column 2: Ratio xz
% Column 3: Ratio yz
particle_ratio_statistics = zeros(5,3);
for column=1:1:3
    % Min
    particle_ratio_statistics(1,column)=min(particle_ratio(:,column+2));
    % Weighted mean
    x_ = particle_ratio(:,column+2); % value
    w_ = mass_center_particle(:,5); % Weigth
    particle_ratio_statistics(2,column)=sum(x_.*w_)/sum(w_);
    % Max
    particle_ratio_statistics(3,column)=max(particle_ratio(:,column+2));
    % Weighted standard deviation
    mw_ = particle_ratio_statistics(2,column);
    n_ = length(x_);
    particle_ratio_statistics(4,column)=sqrt(sum(w_.*(x_-mw_).^2) / ((n_-1)*sum(w_)/n_));
    % Weighted standard deviation, in percent of the weigthed mean
    particle_ratio_statistics(5,column)=particle_ratio_statistics(4,column)*100/particle_ratio_statistics(2,column);
end

table_.line(1).name = 'Minimum';
table_.line(2).name = 'Weigt_Mean';
table_.line(3).name = 'Maximum';
table_.line(4).name = 'Weigh_Std';
table_.line(5).name = 'Std%';
Table_particle_ratio_statistics = table(char(table_.line(:).name),particle_ratio_statistics(:,1),particle_ratio_statistics(:,2),particle_ratio_statistics(:,3),...
    'VariableNames',{'Property' 'Ratio_x_over_y' 'Ratio_x_over_z' 'Ratio_y_over_z'});
% Note: 'char' vertically concatenates character arrays, padding each input array as needed so that each row contains the same number of characters.
if OPTIONS.displaytext==true
    fprintf('PARTICLE SIZE RATIO:\n\n');
    disp(Table_particle_ratio_statistics)
end

results_correlation(current_phase).(['Particle_aspect_ratio_x_over_y_mean_' dpsdname]) = particle_ratio_statistics(2,1);
results_correlation(current_phase).(['Particle_aspect_ratio_x_over_z_mean_' dpsdname]) = particle_ratio_statistics(2,2);
results_correlation(current_phase).(['Particle_aspect_ratio_y_over_z_mean_' dpsdname]) = particle_ratio_statistics(2,3);
results_correlation(current_phase).(['Particle_aspect_ratio_x_over_y_std_' dpsdname]) = particle_ratio_statistics(4,1);
results_correlation(current_phase).(['Particle_aspect_ratio_x_over_z_std_' dpsdname]) = particle_ratio_statistics(4,2);
results_correlation(current_phase).(['Particle_aspect_ratio_y_over_z_std_' dpsdname]) = particle_ratio_statistics(4,3);
results_correlation(current_phase).(['Particle_aspect_ratio_x_over_y_relstd_' dpsdname]) = particle_ratio_statistics(5,1);
results_correlation(current_phase).(['Particle_aspect_ratio_x_over_z_relstd_' dpsdname]) = particle_ratio_statistics(5,2);
results_correlation(current_phase).(['Particle_aspect_ratio_y_over_z_relstd_' dpsdname]) = particle_ratio_statistics(5,3);  
% if strcmp(dpsdname,'watershed')
%     results_correlation(current_phase).Particle_aspect_ratio_x_over_y_mean_Watershed = particle_ratio_statistics(2,1);
%     results_correlation(current_phase).Particle_aspect_ratio_x_over_z_mean_Watershed = particle_ratio_statistics(2,2);
%     results_correlation(current_phase).Particle_aspect_ratio_y_over_z_mean_Watershed = particle_ratio_statistics(2,3);
%     results_correlation(current_phase).Particle_aspect_ratio_x_over_y_std_Watershed = particle_ratio_statistics(4,1);
%     results_correlation(current_phase).Particle_aspect_ratio_x_over_z_std_Watershed = particle_ratio_statistics(4,2);
%     results_correlation(current_phase).Particle_aspect_ratio_y_over_z_std_Watershed = particle_ratio_statistics(4,3);    
%     results_correlation(current_phase).Particle_aspect_ratio_x_over_y_relstd_Watershed = particle_ratio_statistics(5,1);
%     results_correlation(current_phase).Particle_aspect_ratio_x_over_z_relstd_Watershed = particle_ratio_statistics(5,2);
%     results_correlation(current_phase).Particle_aspect_ratio_y_over_z_relstd_Watershed = particle_ratio_statistics(5,3);      
% end
% Save
Table_discrete_particlemorpholgy_analysis.Table_particle_ratio_statistics = Table_particle_ratio_statistics;

%% CUMULATIVE AND PROBABILITY DENSITY DISTRIBUTION FUNCTIONS
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;
[Aspect_ratio_xy_dist, ~] = Function_probability_density(particle_ratio(:,3),mass_center_particle(:,5),density_fct_parameters);
[Aspect_ratio_xz_dist, ~] = Function_probability_density(particle_ratio(:,4),mass_center_particle(:,5),density_fct_parameters);
[Aspect_ratio_yz_dist, ~] = Function_probability_density(particle_ratio(:,5),mass_center_particle(:,5),density_fct_parameters);

array_tmp(:,1) = Aspect_ratio_xy_dist.cumulative_fct(:,1);
array_tmp(:,2) = Aspect_ratio_xy_dist.cumulative_fct(:,2);
array_tmp(:,3) = Aspect_ratio_xy_dist.probability_density_fct(:,2);
if ~isempty(Aspect_ratio_xy_dist.smoothed_cumulative_fct)
    array_tmp(:,4) = Aspect_ratio_xy_dist.smoothed_cumulative_fct(:,1);
    array_tmp(:,5) = Aspect_ratio_xy_dist.smoothed_cumulative_fct(:,2);
    array_tmp(:,6) = Aspect_ratio_xy_dist.smoothed_probability_density_fct(:,2);
    Variable_name_table={'Ratio_x_over_y' 'Cumulative_function' 'Probability_density_distribution_function' 'Micrometers_smoothed' 'Cumulative_function_smoothed' 'Probability_density_function_smoothed'};
else
    Variable_name_table={'Ratio_x_over_y' 'Cumulative_function' 'Probability_density_distribution_function'};
end
Table_Aspect_ratio_distribution.xy = array2table(array_tmp,'VariableNames',Variable_name_table);
clear array_tmp

array_tmp(:,1) = Aspect_ratio_xz_dist.cumulative_fct(:,1);
array_tmp(:,2) = Aspect_ratio_xz_dist.cumulative_fct(:,2);
array_tmp(:,3) = Aspect_ratio_xz_dist.probability_density_fct(:,2);
if ~isempty(Aspect_ratio_xz_dist.smoothed_cumulative_fct)
    array_tmp(:,4) = Aspect_ratio_xz_dist.smoothed_cumulative_fct(:,1);
    array_tmp(:,5) = Aspect_ratio_xz_dist.smoothed_cumulative_fct(:,2);
    array_tmp(:,6) = Aspect_ratio_xz_dist.smoothed_probability_density_fct(:,2);
    Variable_name_table={'Ratio_x_over_z' 'Cumulative_function' 'Probability_density_distribution_function' 'Micrometers_smoothed' 'Cumulative_function_smoothed' 'Probability_density_function_smoothed'};
else
    Variable_name_table={'Ratio_x_over_z' 'Cumulative_function' 'Probability_density_distribution_function'};
end
Table_Aspect_ratio_distribution.xz = array2table(array_tmp,'VariableNames',Variable_name_table);
clear array_tmp

array_tmp(:,1) = Aspect_ratio_yz_dist.cumulative_fct(:,1);
array_tmp(:,2) = Aspect_ratio_yz_dist.cumulative_fct(:,2);
array_tmp(:,3) = Aspect_ratio_yz_dist.probability_density_fct(:,2);
if ~isempty(Aspect_ratio_yz_dist.smoothed_cumulative_fct)
    array_tmp(:,4) = Aspect_ratio_yz_dist.smoothed_cumulative_fct(:,1);
    array_tmp(:,5) = Aspect_ratio_yz_dist.smoothed_cumulative_fct(:,2);
    array_tmp(:,6) = Aspect_ratio_yz_dist.smoothed_probability_density_fct(:,2);
    Variable_name_table={'Ratio_y_over_z' 'Cumulative_function' 'Probability_density_distribution_function' 'Micrometers_smoothed' 'Cumulative_function_smoothed' 'Probability_density_function_smoothed'};
else
    Variable_name_table={'Ratio_y_over_z' 'Cumulative_function' 'Probability_density_distribution_function'};
end
Table_Aspect_ratio_distribution.yz = array2table(array_tmp,'VariableNames',Variable_name_table);
clear array_tmp

Table_discrete_particlemorpholgy_analysis.Aspect_ratio_distribution = Table_Aspect_ratio_distribution; % Save in main table result

% % Text Results are saved
if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Particle_AspectRatio_' dpsdname '_' INFO.phase(current_phase).name];
    full_path=[Current_folder filename];
    % Prepare the data
    clear DATA_writetable;
    % Data : Size ratio
    DATA_writetable.sheet(1).name='Ratio_size';
    DATA_writetable.sheet(1).table=Table_particle_ratio;    
    % Data : Size ratio
    DATA_writetable.sheet(2).name='Ratio_size_statistics';
    DATA_writetable.sheet(2).table=Table_particle_ratio_statistics;
    % Distribution
    DATA_writetable.sheet(3).name='Ratio_x_over_y';
    DATA_writetable.sheet(3).table=Table_Aspect_ratio_distribution.xy;
    DATA_writetable.sheet(4).name= 'Ratio_x_over_z';
    DATA_writetable.sheet(4).table=Table_Aspect_ratio_distribution.xz;        
    DATA_writetable.sheet(5).name= 'Ratio_y_over_z';
    DATA_writetable.sheet(5).table=Table_Aspect_ratio_distribution.yz;        
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

scrsz = get(0,'ScreenSize'); % Screen resolution
% CUMULATIVE AND DISTRIBUTION FUNCTIONS PLOT
parameters_distributionfigure.figureposition = [100 100 1500 800];
parameters_distributionfigure.fontname = OPTIONS.fontname;
parameters_distributionfigure.grid = OPTIONS.grid;
parameters_distributionfigure.minorgrid = OPTIONS.minorgrid;
parameters_distributionfigure.fullpath = Current_folder;
parameters_distributionfigure.save = true;

parameters_distributionfigure.figurename =  ['Particle aspect ratio x over y, ' INFO.phase(current_phase).name];
parameters_distributionfigure.filename = ['Particle_aspect_ratio_xy_' dpsdname '_' INFO.phase(current_phase).name];
parameters_distributionfigure.subaxe1_title = 'Cumulative function';
parameters_distributionfigure.subaxe2_title = 'Distribution function';
parameters_distributionfigure.title = {['Particle aspect ratio ' INFO.direction(1).name '/' INFO.direction(2).name],[INFO.phase(current_phase).name ', ' dpsdname ' method']};
parameters_distributionfigure.xlabel = ['Aspect ratio ' INFO.direction(1).name '/' INFO.direction(2).name];
Aspect_ratio_xy_dist.unit = [];
function_probability_distribution_figure(Aspect_ratio_xy_dist,parameters_distributionfigure);

parameters_distributionfigure.figurename =  ['Particle aspect ratio x over z, ' INFO.phase(current_phase).name];
parameters_distributionfigure.filename = ['Particle_aspect_ratio_xz_' dpsdname '_' INFO.phase(current_phase).name];
parameters_distributionfigure.title = {['Particle aspect ratio ' INFO.direction(1).name '/' INFO.direction(3).name],[INFO.phase(current_phase).name ', ' dpsdname ' method']};
parameters_distributionfigure.xlabel = ['Aspect ratio ' INFO.direction(1).name '/' INFO.direction(3).name];
Aspect_ratio_xz_dist.unit = [];
function_probability_distribution_figure(Aspect_ratio_xz_dist,parameters_distributionfigure);

parameters_distributionfigure.figurename =  ['Particle aspect ratio y over z, ' INFO.phase(current_phase).name];
parameters_distributionfigure.filename = ['Particle_aspect_ratio_yz_' dpsdname '_' INFO.phase(current_phase).name];
parameters_distributionfigure.title = {['Particle aspect ratio ' INFO.direction(2).name '/' INFO.direction(3).name],[INFO.phase(current_phase).name ', ' dpsdname ' method']};
parameters_distributionfigure.xlabel = ['Aspect ratio ' INFO.direction(2).name '/' INFO.direction(3).name];
Aspect_ratio_yz_dist.unit = [];
function_probability_distribution_figure(Aspect_ratio_yz_dist,parameters_distributionfigure);

%% ASPECT RATIO AS FUNCTION OF PARTICLE DIAMETER
% Display aspect ratio (scatter point) as a function of the particle diameter
% - Create figure
Fig_ratio_withdiameter_cloud = figure;
Fig_ratio_withdiameter_cloud.Name= sprintf('Particle aspect ratio of phase: %s',INFO.phase(current_phase).name);
Fig_ratio_withdiameter_cloud.Color='white'; % Background colour
% - Create axes
axes_ = axes('Parent',Fig_ratio_withdiameter_cloud);
hold(axes_,'on');
% - Title
t_=title (' ','FontName','Times New Roman','FontSize',16);
t_.String= {'Particle aspect ratio',[INFO.phase(current_phase).name ', ' dpsdname ' method']};
x_= particle_ratio(:,2);
ratio_xy = particle_ratio(:,3);
ratio_xz = particle_ratio(:,4);
ratio_yz = particle_ratio(:,5);
h_xy_pointcloud = scatter(x_,ratio_xy);
h_xz_pointcloud = scatter(x_,ratio_xz);
h_yz_pointcloud = scatter(x_,ratio_yz);
% Colors, thickness, markers
set(h_xy_pointcloud, 'MarkerEdgeColor', [0 114 189]/255);
set(h_xz_pointcloud, 'MarkerEdgeColor', [217 83 25]/255);
set(h_yz_pointcloud, 'MarkerEdgeColor', [237 177 32]/255);
% - Axis label
t_ = xlabel(' ');
t_1 = sprintf('Equivalent diameter');
t_2 = ' (\mum)';
t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
t_ = ylabel('Particle aspect size ratio ');
% - Legend
h_legend = legend(axes_,{['Ratio ' INFO.direction(1).name '/' INFO.direction(2).name]',['Ratio ' INFO.direction(1).name '/' INFO.direction(3).name]',['Ratio ' INFO.direction(2).name '/' INFO.direction(3).name]'},'Location','best');
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
    filename= sprintf('Particle_aspect_ratio_diameter_%s_%s',dpsdname,INFO.phase(current_phase).name);
    function_savefig(Fig_ratio_withdiameter_cloud, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig_ratio_withdiameter_cloud); % Do not keep open figures
end


% % % EVOLUTION OF PARTICLE SIZE RATIO ALONG DIRECTIONS (GRADIENTS)
   
% % ALGORITHM

% Initialization
% :,1,: position in micrometers
% :,2,: min ratio
% :,3,: mean ratio
% :,4,: max ratio
% :,5,: std
% Position axis in micrometers
for direction = 1:1:3
    Alongdir(direction).ratio = zeros(Domain_size(direction),5,3);
    for ratio_=1:1:3
        Alongdir(direction).ratio(:,1,ratio_) = (1:1:Domain_size(direction))*voxel_size/1000;
    end
end

% Calculation
for direction = 1:1:3
    for ratio_=1:1:3
        if ratio_==1
            data_ratio = particle_ratio_xy_voxel;
        elseif ratio_==2
            data_ratio = particle_ratio_xz_voxel;
        else
            data_ratio = particle_ratio_yz_voxel;
        end
        
        for position=1:1:Domain_size(direction)
            % Current slice of data
            if direction==1
                slice_ = data_ratio(position,:,:);
            elseif direction==2
                slice_ = data_ratio(:,position,:);
            else
                slice_ = data_ratio(:,:,position);
            end
            % Get all ratio
            different_values =  unique(slice_);
            % Remove 0
            different_values(1)=[];
            % Minimum ratio
            Alongdir(direction).ratio(position,2,ratio_) = min(different_values);
            % Maximum ratio
            Alongdir(direction).ratio(position,4,ratio_) = max(different_values);
            % Weighed values
            % Initialisation
            values=zeros(length(different_values),2);
            for current_size=1:1:length(different_values)
                values(current_size,1)=different_values(current_size);
                values(current_size,2)=sum(sum( slice_== different_values(current_size)));
            end
            % Mean ratio
            Alongdir(direction).ratio(position,3,ratio_)=sum(values(:,1).*values(:,2))/sum(values(:,2));
            % Standard deviation
            % The weighted starndard deviation formula is
            % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
            % With wi the weight of the xi, and <x> the weighted mean (mean_size)
            wi = values(:,2);
            xi = values(:,1);
            n = length(xi);
            mean_size = Alongdir(direction).ratio(position,3,ratio_);
            Alongdir(direction).ratio(position,5,ratio_) = sqrt( sum( wi.*((xi-mean_size).^2)) / ( (n-1)/n * sum(wi)  ));
        end
    end
end

% % MANAGING RESULTS
% Results are saved in a table
Variable_name_table={'Position_um' 'min_ratio' 'mean_ratio' 'max_ratio' 'std'};
for direction = 1:1:3
    Table_discrete_particlemorpholgy_analysis.Direction(direction).ratio_xy = array2table(Alongdir(direction).ratio(:,:,1),...
        'VariableNames',Variable_name_table);
    Table_discrete_particlemorpholgy_analysis.Direction(direction).ratio_xz = array2table(Alongdir(direction).ratio(:,:,2),...
        'VariableNames',Variable_name_table);
    Table_discrete_particlemorpholgy_analysis.Direction(direction).ratio_yz = array2table(Alongdir(direction).ratio(:,:,3),...
        'VariableNames',Variable_name_table);
end
clear Variable_name_table;

% % Text Results are saved
if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Particle_AspectRatio_' dpsdname '_along_direction_' INFO.phase(current_phase).name];
    full_path=[Current_folder filename];
    % Prepare the data
    clear DATA_writetable;
    for direction=1:1:3
        s=3*(direction-1);
        DATA_writetable.sheet(s+1).name= ['Dir_' num2str(direction) '_ratio_1_over_2'];
        DATA_writetable.sheet(s+1).table=Table_discrete_particlemorpholgy_analysis.Direction(direction).ratio_xy;
        DATA_writetable.sheet(s+2).name= ['Dir_' num2str(direction) 'ratio_1_over_3'];
        DATA_writetable.sheet(s+2).table=Table_discrete_particlemorpholgy_analysis.Direction(direction).ratio_xz;
        DATA_writetable.sheet(s+3).name= ['Dir_' num2str(direction) 'ratio_2_over_3'];
        DATA_writetable.sheet(s+3).table=Table_discrete_particlemorpholgy_analysis.Direction(direction).ratio_yz;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


% % DISPLAY FIGURES
% % Evolution of particle size along direction
Fig_ = figure;
Fig_.Name= ['Particle size ratio, ' INFO.phase(current_phase).name];
Fig_.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig_,'position',scrsz); % Full screen figure
% Create axes as a subplot
id_axe=0;
for ratio_ = 1:1:3
    if ratio_==1
        ratio_name=['Ratio ' INFO.direction(1).name '/' INFO.direction(2).name];
    elseif ratio_==2
        ratio_name=['Ratio ' INFO.direction(1).name '/' INFO.direction(3).name];
    else
        ratio_name=['Ratio ' INFO.direction(2).name '/' INFO.direction(3).name];
    end
    for direction =1:1:3
        id_axe=id_axe+1;
        sub_axes{id_axe}=subplot(3,3,id_axe,'Parent',Fig_);        
        % Active subplot
        hold(sub_axes{id_axe},'on');         
        % - Title
        t_=title (' ','FontName','Times New Roman','FontSize',16);
        t_.String= ratio_name;
        % Curves data
        x_=Alongdir(direction).ratio(:,1,ratio_);
        y_min = Alongdir(direction).ratio(:,2,ratio_);
        y_mean = Alongdir(direction).ratio(:,3,ratio_);
        y_max = Alongdir(direction).ratio(:,4,ratio_);
        y_std = Alongdir(direction).ratio(:,5,ratio_);
        % Mean
        h_mean=plot(x_,y_mean); % For the legend order
        % Extremums
        %h_min=plot(x_,y_min);
        %h_max=plot(x_,y_max);
        % Colors, thickness, markers
        set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
        %set(h_min, 'Color', OPTIONS.color(current_phase,:),'LineStyle','--','MarkerSize',12,'Marker','none');
        %set(h_max, 'Color', OPTIONS.color(current_phase,:),'LineStyle','--','MarkerSize',12,'Marker','none');
        % Mean with error bar (+- standard deviation)
        h_mean_witherrorbar = errorbar(x_,y_mean,y_std);
        set(h_mean_witherrorbar, 'Color', INFO.phase(current_phase).color,'LineWidth',1,'MarkerSize',12,'Marker','none');
        h_mean=plot(x_,y_mean); % Plot over the other
        set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
        % - Axis label
        t_ = xlabel(' ');
        t_1 = ['Position along ' INFO.direction(direction).name];
        t_2 = ' (\mum)';
        t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
        t_ = ylabel(' ');
        t_1 = 'Particle size aspect ratio';        
        t_.String= [t_1]; % Sprintf does not accept greek characters
        % - Legend
        % legend(axes_direction,'show','Location','best','Minimum ratio','Mean ratio','Maximun ratio');
        % Axis limit
        xlim([0 Domain_size(direction)*voxel_size/1000])
        % - Grid
        grid(sub_axes{id_axe},'on'); % Display grid
        set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
        % - Fontname and fontsize
        set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
        % Release current subplot
        hold(sub_axes{id_axe},'off');
    end
end
% Save figures
filename= sprintf('Particle_size_aspectratio_%s_alongdirection_%s',dpsdname,INFO.phase(current_phase).name);
fullpath=[Current_folder filename];
if OPTIONS.save_fig == true % Save figure
    function_savefig(Fig_, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig_); % Do not keep open figures
end


%%
%% SPHERICITY
%%
%%

% Sphericity = surface area of sphere with same volume / surface area of particle
% Modifief sphericity: the numerator is modified to reflect the surface
% area of a sphere constituted of cubic elements (then *3/2)

% Particle sphericity initialisation
% Id, Equivalent diameter, volume, surface (direct counting), sphericity (direct counting), modified sphericity (direct mounting).
particle_sphericity = zeros(number_particle,6);
% Voxel per voxel (modified sphericity and geo. covariogram)
particle_modified_sphericity_voxel = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% particle_geocova_sphericity_voxel = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    particle_sphericity(current_,1)=particle_id;
    % Equivalent diameter
    particle_sphericity(current_,2)=mass_center_particle(current_,6);
    % Get index
    index_particle = find(Particle_label==particle_id);
    % Get all coordinates
    [M1,M2,M3] = ind2sub(Domain_size,index_particle);
    % Get volume (in number of voxel)
    particle_sphericity(current_,3)=length(M1);
    % Get minimum and maximum
    x_min = min(M1); y_min = min(M2); z_min = min(M3); 
    x_max = max(M1); y_max = max(M2); z_max = max(M3);
    % %Get surface
    % Step 1: create subdomain
    subdomain_size = [x_max-x_min+1+2 y_max-y_min+1+2 z_max-z_min+1+2];
    subdomain = zeros(subdomain_size(1),subdomain_size(2),subdomain_size(3));
    % Step 2: Insert particle id
    subdomain(2:end-1,2:end-1,2:end-1)=Particle_label(x_min:x_max, y_min:y_max, z_min:z_max);
    % Step 3 :binary subdomain
    subdomain(subdomain~=particle_id)=0;
    subdomain(subdomain==particle_id)=1;
    % Step 4: get surface
    
    % Method a: direct counting
    [~,particle_sphericity(current_,4)] = Function_Specificsurface_direct_Algorithm(subdomain);
    
    % Method b: geometric covariogram (see the related function for details)
    % [covariogram_direction1,~] = Function_Covariogram_covariance_Algorithm(subdomain,1);
    % [covariogram_direction2,~] = Function_Covariogram_covariance_Algorithm(subdomain,2);
    % [covariogram_direction3,~] = Function_Covariogram_covariance_Algorithm(subdomain,3);
    % Deduce the specific surface area
    % specific_surface_area = -4*(d(K(V,h,alpha))/dh for h=0)/mes(V)
    % Calculation of D_ = d(K(V,h,alpha))/dh for h=0
    % D_1= (covariogram_direction1(2,1)-covariogram_direction1(1,1))/1;
    % D_2= (covariogram_direction2(2,1)-covariogram_direction2(1,1))/1;
    % D_3= (covariogram_direction3(2,1)-covariogram_direction3(1,1))/1;
    % Surface specific = -4*D_/mes(V)
    % surface_area_1=-4*D_1; surface_area_2=-4*D_2; surface_area_3=-4*D_3;
    % Mean specific surface area
    % particle_sphericity(current_,5) = (surface_area_1+surface_area_2+surface_area_3)/3;
    % % Sphericity
    V_  = particle_sphericity(current_,3);
    %V_modified  = particle_sphericity(current_,3)*3/2;
    Ad_ = particle_sphericity(current_,4);
    % Ac_ = particle_sphericity(current_,5);
    % Direct counting method
    particle_sphericity(current_,5)= min([ (pi^(1/3) * (6*V_)^(2/3))/Ad_, 1 ]);
    % Geo. covariogram method
    % particle_sphericity(current_,7)= (pi^(1/3) * (6*V_)^(2/3))/Ac_;
    % Modified sphericity
    particle_sphericity(current_,6)= min([ (9/(2*Ad_)) * V_^(2/3) * (4*pi/3)^(1/3) , 1 ]);
        
    % Save location
    particle_modified_sphericity_voxel(index_particle) = particle_sphericity(current_,6);
    
    % For very small particles, the geometric covariogram may induce sphericity higher than 1
    % if particle_sphericity(current_,7)>1
    %     particle_sphericity(current_,7)=nan;
    %     particle_geocova_sphericity_voxel(index_particle) = nan;
    % else
    %     particle_geocova_sphericity_voxel(index_particle) = particle_sphericity(current_,7);
    % end
    % if particle_sphericity(current_,8)>1
    %     particle_sphericity(current_,8)=nan;
    %     particle_modified_sphericity_voxel(index_particle) = nan;
    % else
    %     particle_modified_sphericity_voxel(index_particle) = particle_sphericity(current_,8);
    % end    
end
% Check if nan is significant or not
% n_voxel_binaryphase = sum(sum(sum(Particle_label~=0)));
% [row, ~] = find(isnan(particle_geocova_sphericity_voxel));
% n_nan_geocova = length(row);
% ratio_nan_geocova = n_nan_geocova*100/n_voxel_binaryphase;
% [row, ~] = find(isnan(particle_modified_sphericity_voxel));
% n_nan_modifiedsphericity = length(row);
% ratio_nan_modifiedsphericity = n_nan_modifiedsphericity*100/n_voxel_binaryphase;


results_visualization(current_phase).(['particle_sphericity_' dpsdname]) = particle_modified_sphericity_voxel;
%if strcmp(dpsdname,'watershed')
%    results_visualization(current_phase).particle_sphericity_watershed = particle_modified_sphericity_voxel;
%end

% Save
%Variable_name_table={'Particle_Id' 'Equivalent_diameter' 'Volume' 'Surface_direct' 'Surface_geocova' 'Sphericity_direct' 'Sphericity_geocova' 'Sphericity_modified'};
Variable_name_table={'Particle_Id' 'Equivalent_diameter' 'Volume' 'Surface_direct' 'Sphericity_direct' 'Sphericity_modified'};
Table_particle_sphericity = array2table(particle_sphericity,'VariableNames',Variable_name_table);
Table_discrete_particlemorpholgy_analysis.particle_sphericity = Table_particle_sphericity;

% Statistics on particle sphericity
% Particle sphericity statistics initialisation
% Row 1: min
% Row 2: mean, weighted by the particle volume
% Row 3: max
% Row 4: std
% Row 5: std % of the mean
% Column 1: Sphericity (direct counting)
% Column 2: Modified sphericity (direct counting)
particle_sphericity_statistics = zeros(5,2);
for column=1:1:2
    % Min
    particle_sphericity_statistics(1,column)=min(particle_sphericity(:,column+4));
    % value and weight
    x_ = particle_sphericity(:,column+4); % value
    w_ = particle_sphericity(:,3); % Weigth
    % Remove nan
    [row, ~] = find(isnan(x_));
    w_(row)=nan;
    x_(isnan(x_)) = [];
    w_(isnan(w_)) = [];
    % Weighted mean
    particle_sphericity_statistics(2,column)=sum(x_.*w_)/sum(w_);
    % Max
    particle_sphericity_statistics(3,column)=max(particle_sphericity(:,column+4));
    % Weighted standard deviation
    mw_ = particle_sphericity_statistics(2,column);
    n_ = length(x_);
    particle_sphericity_statistics(4,column)=sqrt(sum(w_.*(x_-mw_).^2) / ((n_-1)*sum(w_)/n_));
    % Weighted standard deviation, in percent of the weigthed mean
    particle_sphericity_statistics(5,column)=particle_sphericity_statistics(4,column)*100/particle_sphericity_statistics(2,column);
end

table_.line(1).name = 'Minimum';
table_.line(2).name = 'Weigt_Mean';
table_.line(3).name = 'Maximum';
table_.line(4).name = 'Weigh_Std';
table_.line(5).name = 'Std%';
Table_particle_sphericity_statistics = table(char(table_.line(:).name),particle_sphericity_statistics(:,1),particle_sphericity_statistics(:,2),...
    'VariableNames',{'Property' 'Sphericity_direct_counting' 'Modified_Sphericity_direct_counting'});
% Note: 'char' vertically concatenates character arrays, padding each input array as needed so that each row contains the same number of characters.
if OPTIONS.displaytext==true
    fprintf('PARTICLE SIZE SPHERICITY:\n\n');
    disp(Table_particle_sphericity_statistics)
end


results_correlation(current_phase).(['Particle_sphericity_mean_' dpsdname]) = particle_sphericity_statistics(2,2);
results_correlation(current_phase).(['Particle_sphericity_std_' dpsdname]) = particle_ratio_statistics(4,2);
results_correlation(current_phase).(['Particle_sphericity_relstd_' dpsdname]) = particle_ratio_statistics(5,2);
%if strcmp(dpsdname,'watershed')
%    results_correlation(current_phase).Particle_sphericity_mean_Watershed = particle_sphericity_statistics(2,2);
%    results_correlation(current_phase).Particle_sphericity_std_Watershed = particle_ratio_statistics(4,2);
%    results_correlation(current_phase).Particle_sphericity_relstd_Watershed = particle_ratio_statistics(5,2);
%end

% Save
Table_discrete_particlemorpholgy_analysis.particle_sphericity_statistics = Table_particle_sphericity_statistics;

% % Text Results are saved
if OPTIONS.save_xls==true
    filename = ['Sphericity_' dpsdname '_' INFO.phase(current_phase).name]; % Filename without extension
    % Prepare the data
    clear DATA_writetable;
    % Data : Size ratio
    DATA_writetable.sheet(1).name='Sphericity';
    DATA_writetable.sheet(1).table=Table_particle_sphericity;    
    % Data : Size ratio
    DATA_writetable.sheet(2).name='Sphericity_statistics';
    DATA_writetable.sheet(2).table=Table_particle_sphericity_statistics;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

% Sphericity as a function of diameter

% Display sphericity (scatter point) as a function of the particle diameter
% - Create figure
Fig = figure;
Fig.Name= ['Particle sphericity, ' INFO.phase(current_phase).name];
Fig.Color='white'; % Background colour
% - Create axes
axes_ = axes('Parent',Fig);
hold(axes_,'on');
% - Title
t_=title (' ','FontName','Times New Roman','FontSize',16);
t_.String= ['Particle sphericity, ' INFO.phase(current_phase).name];
x_= particle_sphericity(:,2);
% Modified sphericity
modified_sphericity = particle_sphericity(:,6);
h_modified_sphericity_pointcloud = scatter(x_,modified_sphericity);
% Colors, thickness, markers
set(h_modified_sphericity_pointcloud, 'MarkerEdgeColor',INFO.phase(current_phase).color);
% - Axis label
t_ = xlabel(' ');
t_1 = sprintf('Equivalent diameter');
t_2 = ' (\mum)';
t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
t_ = ylabel('Modified sphericity ');
% yaxis
y_min = floor(10*min(modified_sphericity))/10;
ylim(axes_,[y_min 1]);
% - Grid
grid(axes_,'on'); % Display grid
set(axes_,'XMinorGrid','on','YMinorGrid','on','YMinorTick','on','YTick',[y_min:0.1:1]); % Display grid for minor thicks also
% - Fontname and fontsize
set(axes_,'FontName','Times New Roman','FontSize',14);
axis tight
% - Figure has been done
ylim(axes_,[y_min 1]);
hold(axes_,'off');
% Save figures
filename= sprintf('Sphericity_diameter_%s_%s',dpsdname,INFO.phase(current_phase).name);
if OPTIONS.save_fig == true % Save figure
    function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end

% % % EVOLUTION OF PARTICLE SPHERICITY ALONG DIRECTIONS
    
% % ALGORITHM

% Initialization
% :,1 position in micrometers
% :,2 min sphericity
% :,3 mean sphericity
% :,4 max sphericity
% :,5 std
for direction = 1:1:3
    Alongdir(direction).sphericity = zeros(Domain_size(direction),5);
    Alongdir(direction).sphericity(:,1) = (1:1:Domain_size(direction))*voxel_size/1000;
end

% Calculation
data_ = particle_modified_sphericity_voxel;
for direction = 1:1:3
    for position=1:1:Domain_size(direction)
        % Current slice of data
        if direction==1
            slice_ = data_(position,:,:);
        elseif direction==2
            slice_ = data_(:,position,:);
        else
            slice_ = data_(:,:,position);
        end
        % Get all ratio
        different_values =  unique(slice_);
        % Remove 0
        different_values(1)=[];
        % Remove nan
        different_values(isnan(different_values)) = [];
        % Minimum ratio
        Alongdir(direction).sphericity(position,2) = min(different_values);
        % Maximum ratio
        Alongdir(direction).sphericity(position,4) = max(different_values);
        % Weighed values
        % Initialisation
        values=zeros(length(different_values),2);
        for current_size=1:1:length(different_values)
            values(current_size,1)=different_values(current_size);
            values(current_size,2)=sum(sum( slice_== different_values(current_size)));
        end
        % Mean ratio
        Alongdir(direction).sphericity(position,3)=sum(values(:,1).*values(:,2))/sum(values(:,2));
        % Standard deviation
        % The weighted starndard deviation formula is
        % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
        % With wi the weight of the xi, and <x> the weighted mean (mean_size)
        wi = values(:,2);
        xi = values(:,1);
        n = length(xi);
        mean_size = Alongdir(direction).sphericity(position,3);
        Alongdir(direction).sphericity(position,5) = sqrt( sum( wi.*((xi-mean_size).^2)) / ( (n-1)/n * sum(wi)  ));
    end
end


% % MANAGING RESULTS
% Results are saved in a table
Variable_name_table={'Position_um' 'min_ratio' 'mean_ratio' 'max_ratio' 'std'};
for direction = 1:1:3
    Table_discrete_particlemorpholgy_analysis.Direction(direction).sphericity = array2table(Alongdir(direction).sphericity(:,:),...
        'VariableNames',Variable_name_table);
end
clear Variable_name_table;

% % Text Results are saved
if OPTIONS.save_xls==true
    % Filename without extension
    filename = ['Particle_Sphericity_' dpsdname '_along_direction_' INFO.phase(current_phase).name];
    full_path=[Current_folder filename];
    % Prepare the data
    clear DATA_writetable;
    for direction=1:1:3
        DATA_writetable.sheet(direction).name= ['Dir_' num2str(direction) '_ratio_1_over_2'];
        DATA_writetable.sheet(direction).table=Table_discrete_particlemorpholgy_analysis.Direction(direction).sphericity;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

% % DISPLAY FIGURES
% % Evolution of particle sphericity along direction
Fig_ = figure;
Fig_.Name= ['Particle sphericity, ' INFO.phase(current_phase).name];
Fig_.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig_,'position',scrsz); % Full screen figure
% Create axes as a subplot
id_axe=0;
for direction =1:1:3
    id_axe=id_axe+1;
    sub_axes{id_axe}=subplot(1,3,id_axe,'Parent',Fig_);
    % Active subplot
    hold(sub_axes{id_axe},'on');
    % - Title
    t_=title ('Sphericity','FontName','Times New Roman','FontSize',16);
    % Curves data
    x_=Alongdir(direction).sphericity(:,1);
    y_min = Alongdir(direction).sphericity(:,2);
    y_mean = Alongdir(direction).sphericity(:,3);
    y_max = Alongdir(direction).sphericity(:,4);
    y_std = Alongdir(direction).sphericity(:,5);
    % Mean
    h_mean=plot(x_,y_mean); % For the legend order
    % Extremums
    %h_min=plot(x_,y_min);
    %h_max=plot(x_,y_max);
    % Colors, thickness, markers
    set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
    %set(h_min, 'Color', OPTIONS.color(current_phase,:),'LineStyle','--','MarkerSize',12,'Marker','none');
    %set(h_max, 'Color', OPTIONS.color(current_phase,:),'LineStyle','--','MarkerSize',12,'Marker','none');
    % Mean with error bar (+- standard deviation)
    h_mean_witherrorbar = errorbar(x_,y_mean,y_std);
    set(h_mean_witherrorbar, 'Color', INFO.phase(current_phase).color,'LineWidth',1,'MarkerSize',12,'Marker','none');
    h_mean=plot(x_,y_mean); % Plot over the other
    set(h_mean, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
    % - Axis label
    t_ = xlabel(' ');
    t_1 = ['Position along ' INFO.direction(direction).name];
    t_2 = ' (\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' ');
    t_1 = 'Particle sphericity';
    t_.String= [t_1]; % Sprintf does not accept greek characters
    % - Legend
    % legend(axes_direction,'show','Location','best','Minimum ratio','Mean ratio','Maximun ratio');
    % Axis limit
    xlim([0 Domain_size(direction)*voxel_size/1000])
    % - Grid
    grid(sub_axes{id_axe},'on'); % Display grid
    set(sub_axes{id_axe},'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(sub_axes{id_axe},'FontName','Times New Roman','FontSize',14);
    % Release current subplot
    hold(sub_axes{id_axe},'off');
end
% Save figures
filename= sprintf('Particle_sphericity_%s_alongdirection_%s',dpsdname,INFO.phase(current_phase).name);
fullpath=[Current_folder filename];
if OPTIONS.save_fig == true % Save figure
    function_savefig(Fig_, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig_); % Do not keep open figures
end


end

