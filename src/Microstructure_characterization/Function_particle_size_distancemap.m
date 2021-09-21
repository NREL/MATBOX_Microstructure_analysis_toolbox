function [] = Function_particle_size_distancemap(Phase_microstructure, PROPERTY, OPTIONS, INFO)
%Calculate Particle size fitting the cumulative function of the euclidean
%distance map
% Function_particle_size_distancemap(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_particle_size_distancemap(array, voxelsize) - when use as a standalone function

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 2 % Case for function called as: Function_particle_size_CPSD(Phase_microstructure, voxel_size)
        voxel_size = PROPERTY;
        clear PROPERTY OPTIONS
        
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        INFO.resultfoldername = ['Volume_characterization_' char(t)];
        desktop=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
        OPTIONS.mainsavefolder = [desktop '\' INFO.resultfoldername];           
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        
        % Set default phase information
        unique_code = unique(Phase_microstructure);
        INFO.number_phase = length(unique_code);
        INFO.voxel_number = numel(Phase_microstructure);
        for k=1:1:INFO.number_phase % Loop over all unique values
            INFO.phase(k).code = unique_code(k); % Assing code
            INFO.phasename(k,1) = {['Phase ' num2str(INFO.phase(k).code)]}; % Assign phase name based on code value
            INFO.phase(k).filename = ['Phase_' num2str(INFO.phase(k).code)]; % Assign phase name based on code value
            INFO.phase(k).name = INFO.phase(k).filename;
            INFO.phase(k).color = allcolor(k,:);
        end
        INFO.voxel_size = voxel_size; % nanometers
        
        % Set default direction information
        INFO.direction(1).name = 'in-plane direction 1';
        INFO.direction(2).name = 'in-plane direction 2';
        INFO.direction(3).name = 'through-plane direction';
        INFO.direction(1).filename = 'inplane_direction_1';
        INFO.direction(2).filename = 'inplane_direction_2';
        INFO.direction(3).filename = 'throughplane_direction';
        
        % Set default options
        OPTIONS.save_resultsmat = true;
        OPTIONS.save_xls = true;
        OPTIONS.save_fig = true;
        OPTIONS.savefig_infig = true;
        OPTIONS.savefig_informat = {'png'};
        
        % Set display options
        OPTIONS.fontname = 'Times New Roman';
        OPTIONS.displaytext = true;
        OPTIONS.closefigureaftercreation = false;
        OPTIONS.grid = 'on'; OPTIONS.minorgrid = 'on';
        OPTIONS.Linewidth = 2;
        OPTIONS.Fontsize_axe =  12;
        OPTIONS.Fontsize_legend =  12;
        OPTIONS.Fontsize_title =  14;
        
        % No Voxel size dependence analysis and RVE analysis
        PROPERTY.particlesize_dmap.voxel_size_dependence.todo = false;
        PROPERTY.particlesize_dmap.number_RVE = 0;
        
    else % Incorrect number of argument
        disp 'Error calling Function_particle_size_distancemap. Wrong number of argument.'
        help Function_particle_size_distancemap
    end
    
end


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Particle_size_dmap\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Particle_size_dmap/'; % Result are saved in this subfolder
end
Current_folder = [main_folder Sub_folder]; 
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% INFO
Domain_size = size(Phase_microstructure);
if length(Domain_size)==2 % Add third column with unit value
    Domain_size = [Domain_size 1];
end
if min(Domain_size) == 1 % 2D case
    number_dimension = 2;
else
    number_dimension =3; % 3D case
end
number_phase = INFO.number_phase; % Number of phase
voxel_number = INFO.voxel_number; % Number of voxel
voxel_size = INFO.voxel_size; % Voxel size nm

%% INITIALIZE RESULTS (USE FOR CORRELATION and VISUALIZATION)
for current_phase=1:1:number_phase
    results_correlation(current_phase).name = INFO.phase(current_phase).name;
    results_visualization(current_phase).name = INFO.phase(current_phase).name;
end

%% DISPLAY
if OPTIONS.displaytext==true
    disp '    PARTICLE SIZE - DISTANCE MAP CUMULATIVE FUNCTION FITTING METHOD';
    disp '    ---------------------------------------------------------------';
    disp ' ';
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% CALCULATION
time_cpu_start = cputime; % CPU start
tic; % Stopwatch start
% Initialization
dmap_results = zeros(number_phase,5); % min, mean, max, std, std%
d50_results = zeros(number_phase,1); % mean
for current_phase=1:1:number_phase % Loop over all phases
    code_tmp = INFO.phase(current_phase).code;
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3)); % Initialization
    binary_phase(Phase_microstructure == code_tmp) = 1;
    [distance_transform, fitted_diameter, numericalpsd(current_phase).psd, analyticalpsd(current_phase).psd, error_radius(current_phase).val] = Function_particle_size_distancemap_Algorithm(binary_phase, voxel_size);
    
    % Statistics
    all_distance = distance_transform(binary_phase==1);
    dmap_results(current_phase,1) = min(all_distance); % Minimum
    dmap_results(current_phase,2) = mean(all_distance);% Mean
    dmap_results(current_phase,3) = max(all_distance); % Maximum   
    dmap_results(current_phase,4) = std(all_distance); % Standard deviation    
    dmap_results(current_phase,5) = 100*dmap_results(current_phase,4)/dmap_results(current_phase,2); % Relative standard deviation      
    
    d50_results(current_phase,1) = fitted_diameter;
    
    PSD_results(current_phase,1) = numericalpsd(current_phase).psd.x50; % d50
    PSD_results(current_phase,2) = numericalpsd(current_phase).psd.smoothed_x50; % smoothed d50   
    PSD_results(current_phase,3) = numericalpsd(current_phase).psd.integral_probability_density_fct; 
    PSD_results(current_phase,4) = numericalpsd(current_phase).psd.integral_smoothed_probability_density_fct;    
    
    % Save (correlation, visualization)
    results_correlation(current_phase).Distance_to_surface_mean = dmap_results(current_phase,2);
    results_correlation(current_phase).Distance_to_surface_std = dmap_results(current_phase,4);
    results_correlation(current_phase).Distance_to_surface_relstd = dmap_results(current_phase,5);    
    results_correlation(current_phase).Particle_diameter_mean_dmap = d50_results(current_phase,1);
    
    results_visualization(current_phase).distance_to_boundary = distance_transform;    
end
% CPU and stopwatch time - end
time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc; % Stopwatch elapsed time

%% TABLES
% Time
Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
Table_time = table(Time_measure(1)*1e-6,Time_measure(2),Time_measure(3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_dmap.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_dmap = table(INFO.phasename,dmap_results(:,1),dmap_results(:,2),dmap_results(:,3),dmap_results(:,4),dmap_results(:,5),...
    'VariableNames',{'Name' 'Min' 'Mean' 'Max' 'Std' 'Std_percents'});
Table_dmap_cumulative = table(INFO.phasename,PSD_results(:,1),PSD_results(:,2),PSD_results(:,3),PSD_results(:,4),...
    'VariableNames',{'Name' 'd50' 'Smoothed_d50' 'Distribution_integral' 'Smoothed_distribution_integral'});
Table_dmap_d50 = table(INFO.phasename,d50_results(:,1),...
    'VariableNames',{'Name' 'Fitted_diameter_um'});
Results_dmap.Table_dmap = Table_dmap; % Save in main table result
Results_dmap.Table_dmap_cumulative = Table_dmap_cumulative;
Results_dmap.Table_dmap_d50 = Table_dmap_d50;

for current_phase=1:1:number_phase
    array_tmp(:,1) = numericalpsd(current_phase).psd.cumulative_fct(:,1);
    array_tmp(:,2) = numericalpsd(current_phase).psd.cumulative_fct(:,2);
    array_tmp(:,3) = numericalpsd(current_phase).psd.probability_density_fct(:,2);
    if ~isempty(numericalpsd(current_phase).psd.smoothed_cumulative_fct)
        array_tmp(:,4) = numericalpsd(current_phase).psd.smoothed_cumulative_fct(:,1);
        array_tmp(:,5) = numericalpsd(current_phase).psd.smoothed_cumulative_fct(:,2);
        array_tmp(:,6) = numericalpsd(current_phase).psd.smoothed_probability_density_fct(:,2);
        Variable_name_table={'Micrometers' 'Cumulative_function' 'Probability_density_distribution_function' 'Micrometers_smoothed' 'Cumulative_function_smoothed' 'Probability_density_function_smoothed'};
    else
        Variable_name_table={'Micrometers' 'Cumulative_function' 'Probability_density_distribution_function'};
    end
    Table_cumulative_sizedistribution.phase(current_phase).table = array2table(array_tmp,...
        'VariableNames',Variable_name_table);
    clear array_tmp
end
Results_dmap.Table_cumulative_sizedistribution = Table_cumulative_sizedistribution; % Save in main table result

%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Euclidean_distance'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Distance_from_boundary';
    DATA_writetable.sheet(1).table=Table_dmap;
    DATA_writetable.sheet(2).name='From_cumulative_fct';
    DATA_writetable.sheet(2).table=Table_dmap_cumulative;        
    DATA_writetable.sheet(3).name='Fitted_diameter';
    DATA_writetable.sheet(3).table=Table_dmap_d50;       
    % Cumulative and size distribution
    sheet_=3;
    for current_phase = 1:1:number_phase
        sheet_=sheet_+1;
        DATA_writetable.sheet(sheet_).name=[INFO.phase(current_phase).filename '_DD'];
        DATA_writetable.sheet(sheet_).table=Table_cumulative_sizedistribution.phase(current_phase).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
if OPTIONS.displaytext==true
    fprintf('> Calculated on the whole domain:\n\n');
    disp 'Distance map, i.e distance from boundary';
    disp(Table_dmap)
    disp 'Distance map, i.e distance from boundary (from cumulative function)';
    disp(Table_dmap_cumulative)
    disp 'Fitted diameter';
    disp(Table_dmap_d50)    
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_time)
end

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%
scrsz = get(0,'ScreenSize'); % Screen resolution

%% CUMULATIVE FUNCTION DIAMTER FITTING
for current_phase = 1:1:number_phase
    Fig= figure;
    Fig.Name= sprintf('Distance map: fitted diameter of phase %s',INFO.phase(current_phase).name);
    Fig.Color='white'; % Background colour
    set(Fig, 'Position', [scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*1/2]);
    sub_axes=subplot(1,2,1,'Parent',Fig); % Create axes
    hold(sub_axes,'on'); % Active subplot
    h_title=title (' ','FontName',OPTIONS.fontname);
    h_title.String= {'Analytical - numerical integral difference of','the distance to surface map cumulative fct'};
    h_error=plot(error_radius(current_phase).val(:,1) ,error_radius(current_phase).val(:,2)); % Curves
    set(h_error,'LineStyle','-','Color','k','LineWidth',OPTIONS.Linewidth); % Colors
    xlabel('Sphere radius (\mum)'); % Axis label
    ylabel('Error');
    if strcmp(OPTIONS.grid,'on') % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
    end
    set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % - Fontname and fontsize
    h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
    hold(sub_axes,'off'); % - Figure has been done
    
    sub_axes=subplot(1,2,2,'Parent',Fig); % Create axes
    hold(sub_axes,'on'); % Active subplot
    h_title=title (' ','FontName',OPTIONS.fontname);
    h_title.String= {'Distance to the surface cumulative fct','calculated for the phase and a unique sphere with a fitted diameter'};
    h_phase = plot(numericalpsd(current_phase).psd.cumulative_fct(:,1), numericalpsd(current_phase).psd.cumulative_fct(:,2));
    h_sphere=plot(analyticalpsd(current_phase).psd(:,1),analyticalpsd(current_phase).psd(:,2));
    set(h_phase,'LineStyle','none','Marker','o','Color',INFO.phase(current_phase).color,'LineWidth',OPTIONS.Linewidth); % Colors
    set(h_sphere,'LineStyle','-','Marker','none','Color','k','LineWidth',OPTIONS.Linewidth);
    xlabel('Distance to the surface (\mum)'); % - Axis label
    ylabel('Cumulative function');
    str1_=['Numerical cumulative funtion of ' INFO.phase(current_phase).name]; % Legend
    str2_=['Analytical cumulative function of a sphere of diameter ' num2str(d50_results(current_phase,1),'%4.2f') ' \mum'];
    legend(sub_axes,str1_,str2_,'Location','best');
    if strcmp(OPTIONS.grid,'on') % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
    end
    set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % - Fontname and fontsize
    h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
    hold(sub_axes,'off'); % - Figure has been done
    sgtitle(Fig,['Euclidean distance map cumulative function diameter fitting, ' INFO.phase(current_phase).name] ,'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('Distance_to_surface_fitted_%s',INFO.phase(current_phase).filename);
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end

%% CUMULATIVE AND DISTRIBUTION FUNCTIONS PLOT
parameters_distributionfigure.figureposition = [100 100 1500 800];
parameters_distributionfigure.fontname = OPTIONS.fontname;
parameters_distributionfigure.grid = OPTIONS.grid;
parameters_distributionfigure.minorgrid = OPTIONS.minorgrid;
parameters_distributionfigure.fullpath = Current_folder;
for current_phase=1:1:number_phase % Loop over all phases
    parameters_distributionfigure.figurename =  ['Distance to boundary, ' INFO.phase(current_phase).name];
    parameters_distributionfigure.filename = ['Distance_to_boundary_' INFO.phase(current_phase).name];
    parameters_distributionfigure.subaxe1_title = 'Cumulative function';
    parameters_distributionfigure.subaxe2_title = 'Distribution function';
    parameters_distributionfigure.title = ['Distance to surface, ' INFO.phase(current_phase).name];
    parameters_distributionfigure.xlabel = 'Distance to surface (\mum)';
    numericalpsd(current_phase).psd.unit = '\mum';
    function_probability_distribution_figure(numericalpsd(current_phase).psd,parameters_distributionfigure);
end

%% ALONG DIRECTIONS
alongdirection_parameters.number_phase = number_phase;
alongdirection_parameters.data = results_visualization;
alongdirection_parameters.field = 'distance_to_boundary';
alongdirection_parameters.number_dimension = number_dimension;
alongdirection_parameters.Domain_size = Domain_size;
alongdirection_parameters.voxel_size = voxel_size/1000;
alongdirection_parameters.ignore_value = 0;
alongdirection_parameters.ignore_min = true;
alongdirection_parameters.Variable_name_table = {'Position_um' 'mean_distance_um' 'max_distance_um' 'std_um'};
alongdirection_parameters.OPTIONS = OPTIONS;
alongdirection_parameters.INFO = INFO;
alongdirection_parameters.Table_filename = 'Distance_to_surface';
alongdirection_parameters.figure_name = 'Distance_to_surface';
alongdirection_parameters.axe_title = 'Distance to surface';
alongdirection_parameters.figure_title = 'Distance to surface along directions';
alongdirection_parameters.figure_filename = 'Distance_to_surface_along_directions';
alongdirection_parameters.ylabel = 'Distance to surface (\mum)';
alongdirection_parameters.ylabel_unit = '\mum';
alongdirection_parameters.legendname = 'distance';
alongdirection_parameters.mean_val = PSD_results(:,2);
alongdirection_parameters.Current_folder = Current_folder;
[Table_evolution] = Function_along_direction(alongdirection_parameters); % Call function
Results_dmap.Table_evolution = Table_evolution; % Save in main table result

%% PARTICLE DISTANCE TO SURFACE MAP
for current_phase=1:1:number_phase % Loop over phases
    data_distance = results_visualization(current_phase).distance_to_boundary;
    Fig = figure; % Create figure
    Fig.Name= ['Distance to surface, phase ' INFO.phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
    for current_direction=1:1:number_dimension % Iterate over axe
        sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        h_title=title ({'Slice in the middle',['View normal to ' INFO.direction(current_direction).name]}); % Set title font
        % - Plot graphs
        if current_direction==1
            tmp=squeeze(data_distance(round(Domain_size(1)/2),:,:));
            image(tmp,'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',INFO.direction(3).name);
            t_1y = sprintf('Position along %s ',INFO.direction(2).name);
        elseif current_direction==2
            tmp=squeeze(data_distance(:,round(Domain_size(2)/2),:));
            image(tmp,'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',INFO.direction(3).name);
            t_1y = sprintf('Position along %s ',INFO.direction(1).name);
        elseif current_direction==3
            image(data_distance(:,:,round(Domain_size(3)/2)),'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',INFO.direction(1).name);
            t_1y = sprintf('Position along %s ',INFO.direction(2).name);
        end
        axis equal; axis tight;
        x_value = get(sub_axes,'XTick');
        set(sub_axes,'XtickLabel',x_value*voxel_size/1000);
        y_value = get(sub_axes,'YTick');
        set(sub_axes,'YtickLabel',y_value*voxel_size/1000);
        % - Axis label
        t_ = xlabel(' ');
        t_2 = '(\mum)';
        t_.String= [t_1x t_2]; % Sprintf does not accept greek characters
        t_ = ylabel(' ');
        t_2 = '(\mum)';
        t_.String= [t_1y t_2]; % Sprintf does not accept greek characters
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
        % Color map
        myColorMap = jet(256);
        myColorMap(1,:) = 1; % background color is white
        colormap(myColorMap);
        % Create colorbar
        %colorbar('peer',sub_axes);
        h=colorbar(sub_axes);
        ylabel(h, 'Distance to surface (\mum)');
        set(h,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe);        
        h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,['Distance to surface, ' INFO.phase(current_phase).name],'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('View_distance_%s', INFO.phase(current_phase).filename);
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end


%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%
interpolation_voxelsize_order=1;

if PROPERTY.particlesize_dmap.voxel_size_dependence.todo % Check if voxel size analysis is asked
    size_choice = PROPERTY.particlesize_dmap.voxel_size_dependence.voxel;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed

    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_phase+1,2);
    property_voxelsizedependence(1,1,:)=voxel_size;
    property_voxelsizedependence(1,2:end,1)=d50_results(:,1)'; % Fitted diameter
    property_voxelsizedependence(1,2:end,2)=dmap_results(:,2)'; % Mean distance to surface
    % Loop on each voxel size
    for current_phase=1:1:number_phase
        PSDresized(current_phase).iteration(1).voxelsize = voxel_size/1000;
        PSDresized(current_phase).iteration(1).psd = numericalpsd(current_phase).psd;
    end
    for current_iteration=1:1:number_resize
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1,:)=current_voxel_size;
        % Microstructure resized
        [Phase_microstructure_resized] = function_scale_array(Phase_microstructure, voxel_size, current_voxel_size, INFO.phaseinfo);        
        Current_domain_size=size(Phase_microstructure_resized);
        % CPU and stopwatch time - start
        time_cpu_start = cputime;
        tic;
        for current_phase=1:1:number_phase
            code_tmp = INFO.phase(current_phase).code; % Code of the phase
            binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Binary phase
            binary_phase(Phase_microstructure_resized == code_tmp) = 1;
            
            [distance_transform_resized, fitted_diameter_resized, numericalpsd_resized(current_phase).psd, ~, ~] = Function_particle_size_distancemap_Algorithm(binary_phase, current_voxel_size);
            all_distance_resized = distance_transform_resized(binary_phase==1);
            property_voxelsizedependence(current_iteration+1,current_phase+1,1) = fitted_diameter_resized;
            property_voxelsizedependence(current_iteration+1,current_phase+1,2) = mean(all_distance_resized);

            PSDresized(current_phase).iteration(current_iteration+1).voxelsize = current_voxel_size/1000;
            PSDresized(current_phase).iteration(current_iteration+1).psd = numericalpsd_resized(current_phase).psd;
        end        
        % Number of voxel of the current resized microstructure
        voxel_number_tmp=numel(Phase_microstructure_resized);       
        % CPU and stopwatch time - end
        time_cpu_elapsed = cputime-time_cpu_start;
        time_stopwatch_elapsed = toc;
        Time_tmp = [voxel_number_tmp time_cpu_elapsed time_stopwatch_elapsed];
        Time_measure = [Time_measure;Time_tmp];
    end
    clear Phase_microstructure_resized;
    
    
    %% EXTRAPOLATION TO 0 nm
    str_correlation(1).name = 'Particle_diameter_mean_dmap';
    str_correlation(2).name = 'Distance_to_surface_mean';
    tmp = zeros(number_resize+2,number_phase+1,2);
    for property_resized = 1:1:2
        x=property_voxelsizedependence(:,1,property_resized);
        for current_phase=1:1:number_phase
            y=property_voxelsizedependence(:,current_phase+1,property_resized);
            p = polyfit(x,y,interpolation_voxelsize_order);
            vq = polyval(p,0);
            tmp(1,current_phase+1,property_resized)=vq;
            interpolation_voxelsize(current_phase,property_resized).p=p;
            % For correlation
            results_correlation(current_phase).([str_correlation(property_resized).name '_extrapolated']) = vq;
        end
        tmp(2:end,:,property_resized) = property_voxelsizedependence(:,:,property_resized);
    end
    property_voxelsizedependence = tmp; clear tmp;
        
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={'Voxel_size_nm'}; % Columns name
    for current_phase=1:1:number_phase
        Variable_name_table(1+current_phase)={INFO.phase(current_phase).filename};
    end
    % Table
    Table_d50_Dmap_voxelsizedependence = array2table(property_voxelsizedependence(:,:,1),...
        'VariableNames',Variable_name_table);
    Table_distance2surface_voxelsizedependence = array2table(property_voxelsizedependence(:,:,2),...
        'VariableNames',Variable_name_table);    
    
    %% DISPLAY TEXT RESULTS
    if (OPTIONS.displaytext==1)
        fprintf('> Particle diameter dependence with the voxel size:\n\n');
        disp(Table_d50_Dmap_voxelsizedependence)
        fprintf('> Distance to surface dependence with the voxel size:\n\n');
        disp(Table_distance2surface_voxelsizedependence)        
    end
       
    %% SAVE RESULTS
    Results_dmap.voxelsizedependence_d50 = Table_d50_Dmap_voxelsizedependence; % Save in main table result
    Results_dmap.voxelsizedependence_d2surface = Table_distance2surface_voxelsizedependence; % Save in main table resu    
    if OPTIONS.save_xls==true
        filename = 'Distancemap_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Particle_diameter';
        DATA_writetable.sheet(1).table=Table_d50_Dmap_voxelsizedependence;
        DATA_writetable.sheet(2).name='Distance_to_surface';
        DATA_writetable.sheet(2).table=Table_distance2surface_voxelsizedependence;        
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end        
        
    %% FIGURES
    parameters_figure.number_phase = number_phase;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.INFO = INFO;
    parameters_figure.OPTIONS = OPTIONS;
    
    parameters_figure.propertyname = 'Mean diameter';
    parameters_figure.method = 'Distance map';
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,1);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,1);
    parameters_figure.str_ylabel = 'D_{50} (\mum)';
    parameters_figure.propertynameunit = '\mum';
    parameters_figure.filename = 'D50_dmap_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures

    parameters_figure.propertyname = 'Mean distance to surface';
    parameters_figure.method = [];
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,2);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,2);
    parameters_figure.str_ylabel = 'Distance to surface (\mum)';
    parameters_figure.propertynameunit = '\mum';
    parameters_figure.filename = 'Distance_to_surface_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    
    
    %% FIGURES 2 (DISTRIBUTION)
    for current_phase=1:1:number_phase % Loop over all phases
        parameters_distributionfigure.figurename =  ['Distance to surface, ' INFO.phase(current_phase).name];
        parameters_distributionfigure.filename = ['Distance_to_surface_voxelsize_' INFO.phase(current_phase).name];
        parameters_distributionfigure.subaxe1_title = 'Cumulative function';
        parameters_distributionfigure.subaxe2_title = 'Distribution function';
        parameters_distributionfigure.title = ['Distance to surface voxel size dependence, ' INFO.phase(current_phase).name];
        parameters_distributionfigure.xlabel = 'Distance to surface (\mum)';
        PSDresized(current_phase).unit = '\mum';
        function_probability_distribution_size_figure(PSDresized(current_phase),parameters_distributionfigure);
    end
    
end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) ANALYSIS
%%

if PROPERTY.particlesize_dmap.number_RVE>0
    for n_RVE=1:1:PROPERTY.particlesize_dmap.number_RVE % Loop over all RVE asked
        RVEparameters.name = PROPERTY.particlesize_dmap.RVE(n_RVE).name;
        RVEparameters.savename = PROPERTY.particlesize_dmap.RVE(n_RVE).savename;
        RVEparameters.type = PROPERTY.particlesize_dmap.RVE(n_RVE).type;
        RVEparameters.divisions = PROPERTY.particlesize_dmap.RVE(n_RVE).divisions;
        RVEparameters.subs2 = PROPERTY.particlesize_dmap.RVE(n_RVE).subs2;
        RVEparameters.subs4 = PROPERTY.particlesize_dmap.RVE(n_RVE).subs4;
        RVEparameters.Aspectratio = PROPERTY.particlesize_dmap.RVE(n_RVE).Aspectratio;
        if  strcmp(PROPERTY.particlesize_dmap.RVE(n_RVE).type,'A')
            RVEparameters.Aspectratio_name = [num2str(Domain_size(1)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(2)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(3)/Domain_size(3),'%1.3f\t')];
        elseif strcmp(PROPERTY.particlesize_dmap.RVE(n_RVE).type,'B') || strcmp(PROPERTY.particlesize_dmap.RVE(n_RVE).type,'D')
            RVEparameters.Aspectratio_name = [num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
        end
        RVEparameters.Constantdirection = PROPERTY.particlesize_dmap.RVE(n_RVE).Constantdirection;
        RVEparameters.Growthdirection = PROPERTY.particlesize_dmap.RVE(n_RVE).Growthdirection;
        RVEparameters.Growthperstep = PROPERTY.particlesize_dmap.RVE(n_RVE).Growthperstep;
        RVEparameters.Growthrelativeto = PROPERTY.particlesize_dmap.RVE(n_RVE).Growthrelativeto;
        RVEparameters.threshold_std = PROPERTY.particlesize_dmap.RVE(n_RVE).threshold_std;
        RVEparameters.threshold_numbersubvolumes = PROPERTY.particlesize_dmap.RVE(n_RVE).threshold_numbersubvolumes;
        RVEparameters.firstuniquevolume_size = PROPERTY.particlesize_dmap.RVE(n_RVE).firstuniquevolume_size;
        RVEparameters.firstuniquevolume_unit = PROPERTY.particlesize_dmap.RVE(n_RVE).firstuniquevolume_unit;
        res(1).RVE(n_RVE).RVEparameters = RVEparameters;
        res(2).RVE(n_RVE).RVEparameters = RVEparameters;
        
        if OPTIONS.save_xls==true || OPTIONS.save_fig == true
            if ispc
                Sub_folder_RVE = [Current_folder RVEparameters.savename '\'];
            else
                Sub_folder_RVE = [Current_folder RVEparameters.savename '/'];
            end
            if exist(Sub_folder_RVE,'dir')==0 % Folder existence is checked, and created if necessary
                mkdir(Sub_folder_RVE);
            end
        end
        
        [All_subdomain,GROUP_SUBDOMAIN, Wholevolume_size] = Function_get_Subdomains(RVEparameters, Domain_size, voxel_size); % Location of all subdomains
        Wholevolume_size(1:2) = Wholevolume_size(1:2)*voxel_size/1000; % Size are set in micrometer
        
        % Information about subdomains
        tmp = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),...
            'VariableNames',{'Subdomain_Id' 'Group_Id' 'Number_subdomain' 'Equivalent_cubic_length' 'Section_length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
        res(1).RVE(n_RVE).info = tmp;
        res(2).RVE(n_RVE).info = tmp;     
        
        [number_subdomain,~] = size(All_subdomain); % The number of subdomain
        number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size
        
        %% ALGORITHM
        % Initialisation
        % Colunm 1 is the subdomain id
        % Colunm 2 and 3 are the sizes of the subdomain.
        Property_eachsubdomain = zeros(number_subdomain,number_phase+3,2); %1 fitted diameter, 2 mean distance to surface
        % Property calculated for each subdomain
        for subdomain_id = 1:1:number_subdomain
            % Boundary of the subdomain
            x0 = All_subdomain(subdomain_id,6); x1 = All_subdomain(subdomain_id,7);
            y0 = All_subdomain(subdomain_id,8); y1 = All_subdomain(subdomain_id,9);
            z0 = All_subdomain(subdomain_id,10); z1 = All_subdomain(subdomain_id,11);
            clear current_subdomain;
            current_subdomain = Phase_microstructure(x0:x1,y0:y1,z0:z1);
            Current_domain_size = size(current_subdomain);
            current_voxel_number = numel(current_subdomain);
            Property_eachsubdomain(subdomain_id,1,:)=subdomain_id;
            % Equivalent size of the subdomain
            Property_eachsubdomain(subdomain_id,2,:)=All_subdomain(subdomain_id,4)*voxel_size/1000; % Size are set in micrometer
            Property_eachsubdomain(subdomain_id,3,:)=All_subdomain(subdomain_id,5)*voxel_size/1000;            
            % CPU and stopwatch time - start
            time_cpu_start = cputime;
            tic;
            for current_phase=1:1:number_phase % Loop over all phases
                code_tmp = INFO.phase(current_phase).code; % Code of the phase
                binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Initialization
                binary_phase(current_subdomain == code_tmp) = 1;
                [distance_transform_subdomain, fitted_diameter_subdomain, ~, ~, ~] = Function_particle_size_distancemap_Algorithm(binary_phase, voxel_size);
                all_distance_subdomain = distance_transform_subdomain(binary_phase==1);
                Property_eachsubdomain(subdomain_id,current_phase+3,1)=fitted_diameter_subdomain;
                Property_eachsubdomain(subdomain_id,current_phase+3,2)=mean(all_distance_subdomain);
            end
            % CPU and stopwatch time - end
            time_cpu_elapsed = cputime-time_cpu_start;
            time_stopwatch_elapsed = toc;
            Time_tmp = [current_voxel_number time_cpu_elapsed time_stopwatch_elapsed];
            Time_measure = [Time_measure;Time_tmp];
        end
        
        %% STATISTICAL ANALYSIS and RVE SIZE
        str_RVE(1).name = 'D50_dmap_RVE_';
        str_RVE(2).name = 'Distance_to_surface_RVE_';
        str_RVE(1).propertyname = 'Mean diameter';
        str_RVE(2).propertyname = 'Distance to surface';      
        str_RVE(1).propertynameunit = '\mum';
        str_RVE(2).propertynameunit = '\mum';      
        res(1).Wholevolume_results = d50_results(:,1);
        res(2).Wholevolume_results = dmap_results(:,2);        
        
        for property_RVE=1:1:2
            [res(property_RVE).Property_subdomains_statistics, res(property_RVE).Size_RVE] = Function_subdomains_statistical_analysis(number_group_size,number_phase,GROUP_SUBDOMAIN,Property_eachsubdomain(:,:,property_RVE),voxel_size,RVEparameters);
            % For correlation
            for current_phase=1:1:number_phase
                if max(res(property_RVE).Size_RVE(2,:,1))~=0
                    results_correlation(current_phase).([str_RVE(property_RVE).name RVEparameters.savename]) = res(property_RVE).Size_RVE(2,current_phase,1);
                end
                if strcmp(RVEparameters.type,'C')
                    if max(res(property_RVE).Size_RVE(2,:,2))~=0
                        results_correlation(current_phase).([str_RVE(property_RVE).name 'length_' RVEparameters.savename]) = res(property_RVE).Size_RVE(2,current_phase,2);
                    end
                end
            end
            
            %% MANAGING RESULTS
            [res(property_RVE).RVE] = Function_subdomains_manage_results(Property_eachsubdomain(:,:,property_RVE), res(property_RVE).Property_subdomains_statistics, RVEparameters,res(property_RVE).RVE,res(property_RVE).Size_RVE,n_RVE,number_phase,INFO);
            
            %% TEXT DISPLAY AND SAVE RESULTS
            propertyname=str_RVE(property_RVE).propertyname;
            if property_RVE>1
                INFO.showrveparameters=false;
            else
                INFO.showrveparameters=true;
            end
            Function_subdomains_display_and_save(OPTIONS,INFO,res(property_RVE).RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)
            
            %% FIGURES
            parameters_figure.propertyname = propertyname;
            parameters_figure.propertynameunit = str_RVE(property_RVE).propertynameunit;
            parameters_figure.RVE = RVEparameters;
            parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
            parameters_figure.savefolder = Sub_folder_RVE;
            parameters_figure.OPTIONS = OPTIONS;
            parameters_figure.INFO = INFO;
            parameters_figure.number_phase = number_phase;
            parameters_figure.Property_subdomains_statistics = res(property_RVE).Property_subdomains_statistics;
            parameters_figure.Property_eachsubdomain = Property_eachsubdomain(:,:,property_RVE);
            parameters_figure.Size_RVE = res(property_RVE).Size_RVE;
            parameters_figure.Wholevolume_results = res(property_RVE).Wholevolume_results;
            parameters_figure.Wholevolume_size =  Wholevolume_size;
            Function_create_figures_RVE(parameters_figure) % Figures
        end
                
    end
    Results_dmap.RVE.d50 = res(1).RVE; % Save in main table result
    Results_dmap.RVE.dist2surface = res(2).RVE; % Save in main table result
end


%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Diameter_Dmap_calculation_time';
[Results_dmap] = function_time_figure(Time_measure,date_start, date_end, Results_dmap, Current_folder, filename, 'Particle size (distance map)', OPTIONS);
 
%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_particlediameter_Dmap'],'Results_dmap')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Correlation_particlediameter_Dmap'],'results_correlation')
    Sub_folder = 'Visualization\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Visualization_particlediameter_Dmap'],'results_visualization','-v7.3');    
end


end