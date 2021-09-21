function [] = Function_particle_size_CPSD(Phase_microstructure, PROPERTY, OPTIONS, INFO)
%Calculate Particle size with a spherical assumption (C-PSD)
% Function_particle_size_CPSD(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_particle_size_CPSD(array, voxelsize) - when use as a standalone function

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
        PROPERTY.particlesize_cpsd.voxel_size_dependence.todo = false;
        PROPERTY.particlesize_cpsd.number_RVE = 0;
        
    else % Incorrect number of argument
        disp 'Error calling Function_particle_size_CPSD. Wrong number of argument.'
        help Function_particle_size_CPSD
    end
    
end


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Particle_size_Cpsd\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Particle_size_Cpsd/'; % Result are saved in this subfolder
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
    disp '    PARTICLE SIZE - CONTINUUM PARTICLE SIZE DISTRIBUTION (C-PSD) METHOD';
    disp '    -------------------------------------------------------------------';
    disp ' ';
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% CALCULATION
time_cpu_start = cputime; % CPU start
tic; % Stopwatch start
% Initialization
PSD_results = zeros(number_phase,5); % min, mean, max, std, std%
Particle_lvl_description = zeros(number_phase,1);

% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;

for current_phase=1:1:number_phase % Loop over all phases
    code_tmp = INFO.phase(current_phase).code;
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3)); % Initialization
    binary_phase(Phase_microstructure == code_tmp) = 1;
    
    [Particle_size] = Function_particle_size_CPSD_Algorithm(binary_phase); % Call algorithm
    
    % Particle level description 1- (number of voxels that belong to particles
    % of size <= sqrt(3) normalized with the number of voxel of the phase).
    % The higher the better.
    Particle_lvl_description(current_phase,1) = 1 - (sum(sum(sum( logical(double(binary_phase==1) .* double(Particle_size<=sqrt(3))) )))) / (sum(sum(sum(binary_phase==1))));
    
    % Statistics
    Particle_size = Particle_size*voxel_size/1000; % nm->um
    all_diameters = Particle_size(binary_phase==1);
    PSD_results(current_phase,1) = min(all_diameters); % Minimum
    PSD_results(current_phase,2) = mean(all_diameters);% Mean
    PSD_results(current_phase,3) = max(all_diameters); % Maximum   
    PSD_results(current_phase,4) = std(all_diameters); % Standard deviation    
    PSD_results(current_phase,5) = 100*PSD_results(current_phase,4)/PSD_results(current_phase,2); % Relative standard deviation       
   
    % Cumulative and probability density distribution functions
    [PSD(current_phase).psd, ~] = Function_probability_density(all_diameters,[],density_fct_parameters);
    PSD_results(current_phase,6) = PSD(current_phase).psd.x50; % d50
    PSD_results(current_phase,7) = PSD(current_phase).psd.smoothed_x50; % smoothed d50   
    PSD_results(current_phase,8) = PSD(current_phase).psd.integral_probability_density_fct; 
    PSD_results(current_phase,9) = PSD(current_phase).psd.integral_smoothed_probability_density_fct;

    % Save (correlation, visualization)
    results_correlation(current_phase).Particle_diameter_mean_CPSD = PSD_results(current_phase,2);
    results_correlation(current_phase).Particle_diameter_std_CPSD = PSD_results(current_phase,4);
    results_correlation(current_phase).Particle_diameter_relstd_CPSD = PSD_results(current_phase,5);
    %results_correlation(current_phase).Particle_radius_mean_CPSD = PSD_results(current_phase,2)/2; % Example on how to add new correlation
    results_visualization(current_phase).Particle_diameter_CPSD = Particle_size;
    %results_visualization(current_phase).Particle_radius_CPSD = Particle_size/2; % Example on how to add new visualization
end
% CPU and stopwatch time - end
time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc; % Stopwatch elapsed time

%% TABLES
% Time
Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
Table_time = table(Time_measure(1)*1e-6,Time_measure(2),Time_measure(3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_cpsd.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_cpsd = table(INFO.phasename,PSD_results(:,1),PSD_results(:,2),PSD_results(:,3),PSD_results(:,4),PSD_results(:,5),...
    'VariableNames',{'Name' 'Min' 'Mean' 'Max' 'Std' 'Std_percents'});
%Table_radius = table(INFO.phasename,PSD_results(:,1)/2,PSD_results(:,2)/2,PSD_results(:,3)/2,PSD_results(:,4)/2,PSD_results(:,5),...
%   'VariableNames',{'Name' 'Min' 'Mean' 'Max' 'Std' 'Std_percents'}); % Example
Table_cpsd_d50 = table(INFO.phasename,PSD_results(:,6),PSD_results(:,7),PSD_results(:,8),PSD_results(:,9),...
    'VariableNames',{'Name' 'd50' 'Smoothed_d50' 'Distribution_integral' 'Smoothed_distribution_integral'});
Table_cpsd_lvldetail = table(INFO.phasename,Particle_lvl_description(:,1),...
    'VariableNames',{'Name' 'Particle_level_detail'});
Results_cpsd.Table_cpsd = Table_cpsd; % Save in main table result
%Results_cpsd.Table_radius = Table_radius; % Save in main table result (example)
Results_cpsd.Table_cpsd_d50 = Table_cpsd_d50;
Results_cpsd.Table_cpsd_lvldetail = Table_cpsd_lvldetail;

for current_phase=1:1:number_phase
    array_tmp(:,1) = PSD(current_phase).psd.cumulative_fct(:,1);
    array_tmp(:,2) = PSD(current_phase).psd.cumulative_fct(:,2);
    array_tmp(:,3) = PSD(current_phase).psd.probability_density_fct(:,2);
    if ~isempty(PSD(current_phase).psd.smoothed_cumulative_fct)
        array_tmp(:,4) = PSD(current_phase).psd.smoothed_cumulative_fct(:,1);
        array_tmp(:,5) = PSD(current_phase).psd.smoothed_cumulative_fct(:,2);
        array_tmp(:,6) = PSD(current_phase).psd.smoothed_probability_density_fct(:,2);
        Variable_name_table={'Micrometers' 'Cumulative_function' 'Probability_density_distribution_function' 'Micrometers_smoothed' 'Cumulative_function_smoothed' 'Probability_density_function_smoothed'};
    else
        Variable_name_table={'Micrometers' 'Cumulative_function' 'Probability_density_distribution_function'};
    end
    Table_cumulative_sizedistribution.phase(current_phase).table = array2table(array_tmp,...
        'VariableNames',Variable_name_table);
    clear array_tmp
end
Results_cpsd.Table_cumulative_sizedistribution = Table_cumulative_sizedistribution; % Save in main table result

%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Continuum_PSD'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Results';
    DATA_writetable.sheet(1).table=Table_cpsd;
    DATA_writetable.sheet(2).name='Level_details';
    DATA_writetable.sheet(2).table=Table_cpsd_lvldetail;        
    DATA_writetable.sheet(3).name='D50_from_cumulativefct';
    DATA_writetable.sheet(3).table=Table_cpsd_d50;       
    % Cumulative and size distribution
    sheet_=3;
    for current_phase = 1:1:number_phase
        sheet_=sheet_+1;
        DATA_writetable.sheet(sheet_).name=[INFO.phase(current_phase).filename '_PSD'];
        DATA_writetable.sheet(sheet_).table=Table_cumulative_sizedistribution.phase(current_phase).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
if OPTIONS.displaytext==true
    fprintf('> Calculated on the whole domain:\n\n');
    disp 'Particle diameter';
    disp(Table_cpsd)
    disp 'Particle diameter (from cumulative function)';
    disp(Table_cpsd_d50)
    disp 'Particle level of details: 1 - number voxel phase with particle size <sqrt(3)*voxel length / number voxel phase';
    disp(Table_cpsd_lvldetail)    
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_time)
end


%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%
scrsz = get(0,'ScreenSize'); % Screen resolution

%% CUMULATIVE AND DISTRIBUTION FUNCTIONS PLOT
parameters_distributionfigure.figureposition = [100 100 1500 800];
parameters_distributionfigure.fontname = OPTIONS.fontname;
parameters_distributionfigure.grid = OPTIONS.grid;
parameters_distributionfigure.minorgrid = OPTIONS.minorgrid;
parameters_distributionfigure.fullpath = Current_folder;
parameters_distributionfigure.save = true;
for current_phase=1:1:number_phase % Loop over all phases
    parameters_distributionfigure.figurename =  ['Particle diameter, ' INFO.phase(current_phase).name];
    parameters_distributionfigure.filename = ['Diameter_CPSD_' INFO.phase(current_phase).name];
    parameters_distributionfigure.subaxe1_title = 'Cumulative function';
    parameters_distributionfigure.subaxe2_title = 'Distribution function';
    parameters_distributionfigure.title = ['Particle diameter (C-PSD), ' INFO.phase(current_phase).name];
    parameters_distributionfigure.xlabel = 'Particle diameter (\mum)';
    PSD(current_phase).psd.unit = '\mum';
    function_probability_distribution_figure(PSD(current_phase).psd,parameters_distributionfigure);
end
    
%% ALONG DIRECTIONS
% Calculate the function Diameter(x) defined as 1/L*int(Diameter(x)*dx,0,L)=d_50
% For each direction and each phase

alongdirection_parameters.number_phase = number_phase;
alongdirection_parameters.data = results_visualization;
alongdirection_parameters.field = 'Particle_diameter_CPSD';
alongdirection_parameters.number_dimension = number_dimension;
alongdirection_parameters.Domain_size = Domain_size;
alongdirection_parameters.voxel_size = voxel_size/1000;
alongdirection_parameters.ignore_value = 0;
alongdirection_parameters.ignore_min = true;
alongdirection_parameters.Variable_name_table = {'Position_um' 'mean_diameter_um' 'max_diameter_um' 'std_um'};
alongdirection_parameters.OPTIONS = OPTIONS;
alongdirection_parameters.INFO = INFO;
alongdirection_parameters.Table_filename = 'Diameter_CPSD';
alongdirection_parameters.figure_name = 'Diameter_CPSD';
alongdirection_parameters.axe_title = 'Particle diameter (C-PSD)';
alongdirection_parameters.figure_title = 'Particle diameter along directions (C-PSD)';
alongdirection_parameters.figure_filename = 'Diameter_CPSD_along_directions';
alongdirection_parameters.ylabel = 'Diameter (\mum)';
alongdirection_parameters.ylabel_unit = '\mum';
alongdirection_parameters.legendname = 'diameter';
alongdirection_parameters.mean_val = PSD_results(:,2);
alongdirection_parameters.Current_folder = Current_folder;
[Table_evolution] = Function_along_direction(alongdirection_parameters); % Call function
Results_cpsd.Table_evolution = Table_evolution; % Save in main table result

%% PARTICLE DIAMETER MAP
for current_phase=1:1:number_phase % Loop over phases
    data_Particle_size = results_visualization(current_phase).Particle_diameter_CPSD;
    Fig = figure; % Create figure
    Fig.Name= ['Particle diameter, phase ' INFO.phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
    for current_direction=1:1:number_dimension % Iterate over axe
        sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        h_title=title ({'Slice in the middle',['View normal to ' INFO.direction(current_direction).name]}); % Set title font
        % - Plot graphs
        if current_direction==1
            tmp=squeeze(data_Particle_size(round(Domain_size(1)/2),:,:));
            h=image(tmp,'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',INFO.direction(3).name);
            t_1y = sprintf('Position along %s ',INFO.direction(2).name);
            set(h, 'XData', [0, Domain_size(3)*voxel_size/1000]);
            set(h, 'YData', [0, Domain_size(2)*voxel_size/1000]);
        elseif current_direction==2
            tmp=squeeze(data_Particle_size(:,round(Domain_size(2)/2),:));
            h=image(tmp,'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',INFO.direction(3).name);
            t_1y = sprintf('Position along %s ',INFO.direction(1).name);
            set(h, 'XData', [0, Domain_size(3)*voxel_size/1000]);
            set(h, 'YData', [0, Domain_size(1)*voxel_size/1000]);
        elseif current_direction==3
            h=image(data_Particle_size(:,:,round(Domain_size(3)/2)),'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',INFO.direction(1).name);
            t_1y = sprintf('Position along %s ',INFO.direction(2).name);
            set(h, 'XData', [0, Domain_size(1)*voxel_size/1000]);
            set(h, 'YData', [0, Domain_size(2)*voxel_size/1000]);
        end
        axis equal; axis tight;
%         x_value = get(sub_axes,'XTick');
%         set(sub_axes,'XtickLabel',x_value*voxel_size/1000);
%         y_value = get(sub_axes,'YTick');
%         set(sub_axes,'YtickLabel',y_value*voxel_size/1000);
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
        ylabel(h, 'Diameter (\mum)');
        set(h,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe);        
        h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,['Particle diameter (C-PSD), ' INFO.phase(current_phase).name],'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('View_CPSD_%s', INFO.phase(current_phase).filename);
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

if PROPERTY.particlesize_cpsd.voxel_size_dependence.todo % Check if voxel size analysis is asked
    size_choice = PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed

    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_phase+1,2);
    property_voxelsizedependence(1,1,:)=voxel_size;
    property_voxelsizedependence(1,2:end,1)=PSD_results(:,2)';
    property_voxelsizedependence(1,2:end,2)=Particle_lvl_description(:,1)';
    % Loop on each voxel size
    for current_phase=1:1:number_phase
        PSDresized(current_phase).iteration(1).voxelsize = voxel_size/1000;
        PSDresized(current_phase).iteration(1).psd = PSD(current_phase).psd;
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
            [Particle_size_resized] = Function_particle_size_CPSD_Algorithm(binary_phase); % Call algorithm
            
            property_voxelsizedependence(current_iteration+1,current_phase+1,2) = 1 - (sum(sum(sum( logical(double(binary_phase==1) .* double(Particle_size_resized<=sqrt(3))) )))) / (sum(sum(sum(binary_phase==1))));
            
            Particle_size_resized = Particle_size_resized*current_voxel_size/1000; % nm->um
            all_diameters_resized = Particle_size_resized(binary_phase==1);            
            property_voxelsizedependence(current_iteration+1,current_phase+1,1) = mean(all_diameters_resized);
            
            PSDresized(current_phase).iteration(current_iteration+1).voxelsize = current_voxel_size/1000;
            [PSDresized(current_phase).iteration(current_iteration+1).psd, ~] = Function_probability_density(all_diameters_resized,[],density_fct_parameters);
        end        
        % Number of voxel of the current resized microstructure
        voxel_number_tmp=numel(Phase_microstructure_resized);       
        % CPU and stopwatch time - end
        time_cpu_elapsed = cputime-time_cpu_start;
        time_stopwatch_elapsed = toc;
        Time_tmp = [voxel_number_tmp time_cpu_elapsed time_stopwatch_elapsed];
        Time_measure = [Time_measure;Time_tmp];
    end
    clear Phase_microstructure_resized Particle_size_resized;
    
    
    %% EXTRAPOLATION TO 0 nm
    str_correlation(1).name = 'Particle_diameter_mean_CPSD';
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
            if property_resized==1
                results_correlation(current_phase).([str_correlation(property_resized).name '_extrapolated']) = vq;
            end
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
    Table_d50_Cpsd_voxelsizedependence = array2table(property_voxelsizedependence(:,:,1),...
        'VariableNames',Variable_name_table);
    Table_particle_lvl_detail_Cpsd_voxelsizedependence = array2table(property_voxelsizedependence(:,:,2),...
        'VariableNames',Variable_name_table);    
    
    %% DISPLAY TEXT RESULTS
    if (OPTIONS.displaytext==1)
        fprintf('> Mean particle diameter dependence with the voxel size:\n\n');
        disp(Table_d50_Cpsd_voxelsizedependence)
        fprintf('> Particle level of details dependence with the voxel size:\n\n');
        disp(Table_particle_lvl_detail_Cpsd_voxelsizedependence)        
    end
       
    %% SAVE RESULTS
    Results_cpsd.voxelsizedependence_d50 = Table_d50_Cpsd_voxelsizedependence; % Save in main table result
    Results_cpsd.voxelsizedependence_lvldetail = Table_particle_lvl_detail_Cpsd_voxelsizedependence; % Save in main table resu    
    if OPTIONS.save_xls==true
        filename = 'D50_Cpsd_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='D50_Cpsd';
        DATA_writetable.sheet(1).table=Table_d50_Cpsd_voxelsizedependence;
        DATA_writetable.sheet(2).name='Particle_lvl_detail';
        DATA_writetable.sheet(2).table=Table_particle_lvl_detail_Cpsd_voxelsizedependence;        
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end        
        
    %% FIGURES
    parameters_figure.number_phase = number_phase;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.INFO = INFO;
    parameters_figure.OPTIONS = OPTIONS;
    
    parameters_figure.propertyname = 'Mean diameter';
    parameters_figure.method = 'C-PSD';
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,1);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,1);
    parameters_figure.str_ylabel = 'D_{50} (\mum)';
    parameters_figure.propertynameunit = '\mum';
    parameters_figure.filename = 'D50_cpsd_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures

    parameters_figure.propertyname = 'Particle level of details';
    parameters_figure.method = [];
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,2);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,2);
    parameters_figure.str_ylabel = 'Particle level of details';
    parameters_figure.propertynameunit = [];
    parameters_figure.filename = 'Particle_lvl_details_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    
    
    %% FIGURES 2 (DISTRIBUTION)
    for current_phase=1:1:number_phase % Loop over all phases
        parameters_distributionfigure.figurename =  ['Particle diameter, ' INFO.phase(current_phase).name];
        parameters_distributionfigure.filename = ['Diameter_CPSD_voxelsize_' INFO.phase(current_phase).name];
        parameters_distributionfigure.subaxe1_title = 'Cumulative function';
        parameters_distributionfigure.subaxe2_title = 'Distribution function';
        parameters_distributionfigure.title = ['Particle diameter voxel size dependence (C-PSD), ' INFO.phase(current_phase).name];
        parameters_distributionfigure.xlabel = 'Particle diameter (\mum)';
        PSDresized(current_phase).unit = '\mum';
        function_probability_distribution_size_figure(PSDresized(current_phase),parameters_distributionfigure);
    end
    
end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) ANALYSIS
%%

if PROPERTY.particlesize_cpsd.number_RVE>0
    for n_RVE=1:1:PROPERTY.particlesize_cpsd.number_RVE % Loop over all RVE asked
        RVEparameters.name = PROPERTY.particlesize_cpsd.RVE(n_RVE).name;
        RVEparameters.savename = PROPERTY.particlesize_cpsd.RVE(n_RVE).savename;
        RVEparameters.type = PROPERTY.particlesize_cpsd.RVE(n_RVE).type;
        RVEparameters.divisions = PROPERTY.particlesize_cpsd.RVE(n_RVE).divisions;
        RVEparameters.subs2 = PROPERTY.particlesize_cpsd.RVE(n_RVE).subs2;
        RVEparameters.subs4 = PROPERTY.particlesize_cpsd.RVE(n_RVE).subs4;
        RVEparameters.Aspectratio = PROPERTY.particlesize_cpsd.RVE(n_RVE).Aspectratio;
        if  strcmp(PROPERTY.particlesize_cpsd.RVE(n_RVE).type,'A')
            RVEparameters.Aspectratio_name = [num2str(Domain_size(1)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(2)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(3)/Domain_size(3),'%1.3f\t')];
        elseif strcmp(PROPERTY.particlesize_cpsd.RVE(n_RVE).type,'B') || strcmp(PROPERTY.particlesize_cpsd.RVE(n_RVE).type,'D')
            RVEparameters.Aspectratio_name = [num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
        end
        RVEparameters.Constantdirection = PROPERTY.particlesize_cpsd.RVE(n_RVE).Constantdirection;
        RVEparameters.Growthdirection = PROPERTY.particlesize_cpsd.RVE(n_RVE).Growthdirection;
        RVEparameters.Growthperstep = PROPERTY.particlesize_cpsd.RVE(n_RVE).Growthperstep;
        RVEparameters.Growthrelativeto = PROPERTY.particlesize_cpsd.RVE(n_RVE).Growthrelativeto;
        RVEparameters.threshold_std = PROPERTY.particlesize_cpsd.RVE(n_RVE).threshold_std;
        RVEparameters.threshold_numbersubvolumes = PROPERTY.particlesize_cpsd.RVE(n_RVE).threshold_numbersubvolumes;
        RVEparameters.firstuniquevolume_size = PROPERTY.particlesize_cpsd.RVE(n_RVE).firstuniquevolume_size;
        RVEparameters.firstuniquevolume_unit = PROPERTY.particlesize_cpsd.RVE(n_RVE).firstuniquevolume_unit;
        RVE(n_RVE).RVEparameters = RVEparameters;
        
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
        RVE(n_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),...
            'VariableNames',{'Subdomain_Id' 'Group_Id' 'Number_subdomain' 'Equivalent_cubic_length' 'Section_length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
        
        [number_subdomain,~] = size(All_subdomain); % The number of subdomain
        number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size
        
        %% ALGORITHM
        % Initialisation
        % Colunm 1 is the subdomain id
        % Colunm 2 and 3 are the sizes of the subdomain.
        Property_eachsubdomain = zeros(number_subdomain,number_phase+3);
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
            Property_eachsubdomain(subdomain_id,1)=subdomain_id;
            % Equivalent size of the subdomain
            Property_eachsubdomain(subdomain_id,2)=All_subdomain(subdomain_id,4)*voxel_size/1000; % Size are set in micrometer
            Property_eachsubdomain(subdomain_id,3)=All_subdomain(subdomain_id,5)*voxel_size/1000;            
            % CPU and stopwatch time - start
            time_cpu_start = cputime;
            tic;
            for current_phase=1:1:number_phase % Loop over all phases
                code_tmp = INFO.phase(current_phase).code; % Code of the phase
                binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Initialization
                binary_phase(current_subdomain == code_tmp) = 1;
                [Particle_size_subdomain] = Function_particle_size_CPSD_Algorithm(binary_phase); % Call algorithm
                Particle_size_subdomain = Particle_size_subdomain*voxel_size/1000; % nm->um
                all_diameters_subdomain = Particle_size_subdomain(binary_phase==1);
                Property_eachsubdomain(subdomain_id,current_phase+3)=mean(all_diameters_subdomain);
            end
            % CPU and stopwatch time - end
            time_cpu_elapsed = cputime-time_cpu_start;
            time_stopwatch_elapsed = toc;
            Time_tmp = [current_voxel_number time_cpu_elapsed time_stopwatch_elapsed];
            Time_measure = [Time_measure;Time_tmp];
        end
        
        %% STATISTICAL ANALYSIS and RVE SIZE
        [Property_subdomains_statistics, Size_RVE] = Function_subdomains_statistical_analysis(number_group_size,number_phase,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters);
                
        % For correlation
        for current_phase=1:1:number_phase
            if max(Size_RVE(2,:,1))~=0
                results_correlation(current_phase).(['D50_CPSD_RVE_' RVEparameters.savename]) = Size_RVE(2,current_phase,1);
            end
            if strcmp(RVEparameters.type,'C')
                if max(Size_RVE(2,:,2))~=0
                    results_correlation(current_phase).(['D50_CPSD_RVE_length_' RVEparameters.savename]) = Size_RVE(2,current_phase,2);
                end
            end
        end 
        
        %% MANAGING RESULTS
        [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, RVEparameters,RVE,Size_RVE,n_RVE,number_phase,INFO);        
        
        %% TEXT DISPLAY AND SAVE RESULTS
        propertyname='Mean diameter';
        Function_subdomains_display_and_save(OPTIONS,INFO,RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)
        
        %% FIGURES
        parameters_figure.propertyname = propertyname;
        parameters_figure.propertynameunit = '\mum';
        parameters_figure.RVE = RVEparameters;
        parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
        parameters_figure.savefolder = Sub_folder_RVE;
        parameters_figure.OPTIONS = OPTIONS;
        parameters_figure.INFO = INFO;
        parameters_figure.number_phase = number_phase;
        parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
        parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
        parameters_figure.Size_RVE = Size_RVE;
        parameters_figure.Wholevolume_results = PSD_results(:,2);
        parameters_figure.Wholevolume_size = Wholevolume_size;
        Function_create_figures_RVE(parameters_figure) % Figures        
                
    end
    Results_cpsd.RVE.d50 = RVE; % Save in main table result
end


%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Diameter_CPSD_calculation_time';
[Results_cpsd] = function_time_figure(Time_measure,date_start, date_end, Results_cpsd, Current_folder, filename, 'Particle size (C-PSD)', OPTIONS);
 
%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_particlediameter_CPSD'],'Results_cpsd')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Correlation_particlediameter_CPSD'],'results_correlation')
    Sub_folder = 'Visualization\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Visualization_particlediameter_CPSD'],'results_visualization','-v7.3');    
end

end