function [] = Function_Discrete_particle_size_PCRF(Phase_microstructure, PROPERTY, OPTIONS, INFO, foo)
%Identify particles with a PCRF approach and calculate their morphology and size 
% Function_Discrete_particle_size_PCRF(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_Discrete_particle_size_PCRF(array, voxelsize, customfunction, cpsd_refining, details_convergence, visualize_2D) - when use as a standalone function

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 5 % Case for function called as: Function_Discrete_particle_size_PCRF(Phase_microstructure, voxel_size, cpsd_refining, details_convergence, visualize_2D)
        voxel_size = PROPERTY;
        cpsd_refining = OPTIONS;
        details_convergence = INFO;
        visualize_2D = foo;
        clear PROPERTY OPTIONS INFO
        
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
        
        % Parameters
        PROPERTY.particlesize_PCRF.cpsd_refining = cpsd_refining;
        PROPERTY.particlesize_PCRF.details_convergence = details_convergence;        
        PROPERTY.particlesize_PCRF.visualize_2D = visualize_2D; 
        
    else % Incorrect number of argument
        disp 'Error calling Function_Discrete_particle_size_PCRF. Wrong number of argument.'
        help Function_Discrete_particle_size_PCRF
    end
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Particle_size_PCRF\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Particle_size_PCRF/'; % Result are saved in this subfolder
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
    PROPERTY.particlesize_PCRF.visualize_2D = false;
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
    disp '    PARTICLE SIZE - PCRF (D-PSD) METHOD';
    disp '    -------------------------------------------------------------------';
    disp ' ';
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% PARAMETERS
cpsd_refining = PROPERTY.particlesize_PCRF.cpsd_refining;
details_convergence = PROPERTY.particlesize_PCRF.details_convergence;
visualize_2D = PROPERTY.particlesize_PCRF.visualize_2D;

% Hard coded parameter value
exponent_law = 3; % Parameter k of article. Choose 2 or 3.
max_dist = 40; % PCRF local value is calculated on a restricted field of view 2*max_dist;2*max_dist;2*max_dist centered on the local position to reduce CPU time


%% CALCULATION
time_cpu_start = cputime; % CPU start
tic; % Stopwatch start
% Initialization
PSD_results = zeros(number_phase,5); % min, mean, max, std, std%
% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;

% Edge detection
background = 0;
edgewithbackground = false;

for current_phase=1:1:number_phase % Loop over all phases
    fprintf ('Current phase calculated: %s\n',INFO.phase(current_phase).name)
    code_tmp = INFO.phase(current_phase).code;
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3)); % Initialization
    binary_phase(Phase_microstructure == code_tmp) = 1;    
    
    % Particle size and particle identification
    [Particle_size,Particle_label_tmp, Repulsive_matrix, Repulsive_matrix_sign] = Function_Discrete_particle_size_PCRF_algorithm(binary_phase,cpsd_refining,details_convergence,visualize_2D,exponent_law,max_dist);
    
    % Edge, watershed
    [index_border_phase,~,~,~] = Function_identify_edges(binary_phase);
    [index_border_label,~,~,~] = Function_identify_labelsedges(Particle_label_tmp, background, edgewithbackground); 
    index_(current_phase).border_phase = index_border_phase;
    index_(current_phase).border_label = index_border_label;
        
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
    
    % Randomize label
    Particle_label = zeros(size(Particle_label_tmp));
    currentvalues = unique(Particle_label_tmp); tmp = currentvalues; tmp(tmp==0)=[]; tmp=1:1:length(tmp);
    for k=1:1:length(currentvalues)
        if currentvalues(k)~=0
            idx = randi(length(tmp));
            Particle_label( Particle_label_tmp==currentvalues(k) ) = tmp(idx);
            tmp(idx)=[];
        end
    end    
    
    % Save (correlation, visualization)
    results_correlation(current_phase).Particle_diameter_mean_PCRF = PSD_results(current_phase,2);
    results_correlation(current_phase).Particle_diameter_std_PCRF = PSD_results(current_phase,4);
    results_correlation(current_phase).Particle_diameter_relstd_PCRF = PSD_results(current_phase,5);
    results_visualization(current_phase).Particle_diameter_PCRF = Particle_size;
    results_visualization(current_phase).Particle_ParticleLabel_PCRF = Particle_label;
    results_visualization(current_phase).Repulsive_matrix_x = Repulsive_matrix(:,:,:,1);
    results_visualization(current_phase).Repulsive_matrix_y = Repulsive_matrix(:,:,:,2);
    results_visualization(current_phase).Repulsive_matrix_z = Repulsive_matrix(:,:,:,3);
    results_visualization(current_phase).Repulsive_matrix_sign_x = Repulsive_matrix_sign(:,:,:,1);
    results_visualization(current_phase).Repulsive_matrix_sign_y = Repulsive_matrix_sign(:,:,:,2);
    results_visualization(current_phase).Repulsive_matrix_sign_z = Repulsive_matrix_sign(:,:,:,3);
end
% CPU and stopwatch time - end
time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc; % Stopwatch elapsed time


%% TABLES
% Time
Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
Table_time = table(Time_measure(1)*1e-6,Time_measure(2),Time_measure(3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_PCRF.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_PCRF_size = table(INFO.phasename,PSD_results(:,1),PSD_results(:,2),PSD_results(:,3),PSD_results(:,4),PSD_results(:,5),...
    'VariableNames',{'Name' 'Min' 'Mean' 'Max' 'Std' 'Std_percents'});
Table_PCRF_d50 = table(INFO.phasename,PSD_results(:,6),PSD_results(:,7),PSD_results(:,8),PSD_results(:,9),...
    'VariableNames',{'Name' 'd50' 'Smoothed_d50' 'Distribution_integral' 'Smoothed_distribution_integral'});
Results_PCRF.Table_PCRF_size = Table_PCRF_size; % Save in main table result
Results_PCRF.Table_PCRF_d50 = Table_PCRF_d50;

for current_phase=1:1:number_phase
    array_tmp(:,1) = PSD(current_phase).psd.cumulative_fct(:,1);
    array_tmp(:,2) = PSD(current_phase).psd.cumulative_fct(:,2);
    array_tmp(:,3) = PSD(current_phase).psd.probability_density_fct(:,2);
    array_tmp(:,4) = PSD(current_phase).psd.smoothed_cumulative_fct(:,1);
    array_tmp(:,5) = PSD(current_phase).psd.smoothed_cumulative_fct(:,2);
    array_tmp(:,6) = PSD(current_phase).psd.smoothed_probability_density_fct(:,2);
    Variable_name_table={'Micrometers' 'Cumulative_function' 'Probability_density_distribution_function' 'Micrometers_smoothed' 'Cumulative_function_smoothed' 'Probability_density_function_smoothed'};
    Table_cumulative_sizedistribution.phase(current_phase).table = array2table(array_tmp,...
        'VariableNames',Variable_name_table);
    clear array_tmp
end
Results_PCRF.Table_cumulative_sizedistribution = Table_cumulative_sizedistribution; % Save in main table result

%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'PCRF_PSD'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Results';
    DATA_writetable.sheet(1).table=Table_PCRF_size;
    DATA_writetable.sheet(2).name='D50_from_cumulativefct';
    DATA_writetable.sheet(2).table=Table_PCRF_d50;       
    % Cumulative and size distribution
    sheet_=2;
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
    disp(Table_PCRF_size)
    disp 'Particle diameter (from cumulative function)';
    disp(Table_PCRF_d50)
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
    parameters_distributionfigure.filename = ['Diameter_PCRF_' INFO.phase(current_phase).name];
    parameters_distributionfigure.subaxe1_title = 'Cumulative function';
    parameters_distributionfigure.subaxe2_title = 'Distribution function';
    parameters_distributionfigure.title = ['Particle diameter (PCRF), ' INFO.phase(current_phase).name];
    parameters_distributionfigure.xlabel = 'Particle diameter (\mum)';
    PSD(current_phase).psd.unit = '\mum';
    function_probability_distribution_figure(PSD(current_phase).psd,parameters_distributionfigure);
end
    
%% ALONG DIRECTIONS
% Calculate the function Diameter(x) defined as 1/L*int(Diameter(x)*dx,0,L)=d_50
% For each direction and each phase

alongdirection_parameters.number_phase = number_phase;
alongdirection_parameters.data = results_visualization;
alongdirection_parameters.field = 'Particle_diameter_PCRF';
alongdirection_parameters.number_dimension = number_dimension;
alongdirection_parameters.Domain_size = Domain_size;
alongdirection_parameters.voxel_size = voxel_size/1000;
alongdirection_parameters.ignore_value = 0;
alongdirection_parameters.ignore_min = false;
alongdirection_parameters.Variable_name_table = {'Position_um' 'min_diameter_um' 'mean_diameter_um' 'max_diameter_um' 'std_um'};
alongdirection_parameters.OPTIONS = OPTIONS;
alongdirection_parameters.INFO = INFO;
alongdirection_parameters.Table_filename = 'Diameter_PCRF';
alongdirection_parameters.figure_name = 'Diameter_PCRF';
alongdirection_parameters.axe_title = 'Particle diameter (PCRF)';
alongdirection_parameters.figure_title = 'Particle diameter along directions (PCRF)';
alongdirection_parameters.figure_filename = 'Diameter_PCRF_along_directions';
alongdirection_parameters.ylabel = 'Diameter (\mum)';
alongdirection_parameters.ylabel_unit = '\mum';
alongdirection_parameters.legendname = 'diameter';
alongdirection_parameters.mean_val = PSD_results(:,2);
alongdirection_parameters.Current_folder = Current_folder;
[Table_evolution] = Function_along_direction(alongdirection_parameters); % Call function
Results_PCRF.Table_evolution = Table_evolution; % Save in main table result


%% PARTICLE DIAMETER/LABEL MAP

% Color
lake_id_RGB_color = 0.1 + (0.9-0.1).*rand(1e6,3);
lake_id_grey_color = sum(lake_id_RGB_color,2)/3;
for current_phase=1:1:number_phase % Loop over phases
    data_Particle_size = results_visualization(current_phase).Particle_diameter_PCRF;
    data_Particle_label = results_visualization(current_phase).Particle_ParticleLabel_PCRF;
    volume_color = zeros(Domain_size(1),Domain_size(2),Domain_size(3),3); % RGB color map
    volume_grey = zeros(Domain_size); % Blue color map
    % Grey
    currentvalues = unique(data_Particle_label);
    for k=1:1:length(currentvalues)
        idx=find(data_Particle_label==currentvalues(k));
        if currentvalues(k)==0
            volume_grey(idx) = 1; % Background
        else
            volume_grey(idx) = lake_id_grey_color(k);
        end
    end
    % Edge
    volume_grey(index_(current_phase).border_phase) = 0;
    volume_color(:,:,:,1) = volume_grey;
    volume_color(:,:,:,2) = volume_grey;
    volume_color(:,:,:,3) = volume_grey;
    % Watershed lines
    tmp = volume_grey; tmp(index_(current_phase).border_label) = 1;
    volume_color(:,:,:,1)=tmp; % Attribute RGB color
    tmp = volume_grey; tmp(index_(current_phase).border_label) = 0;
    volume_color(:,:,:,2)=tmp;
    tmp = volume_grey; tmp(index_(current_phase).border_label) = 0;
    volume_color(:,:,:,3)=tmp;

    Fig = figure; % Create figure
    Fig.Name= ['Particle label and diameter, phase ' INFO.phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',scrsz); % Full screen figure
    kk=0;
    for row=1:1:3
        if row==1
            str_ = 'Particle size'; tmpdata = data_Particle_size; myColorMap = jet(256); myColorMap(1,:) = 1; % background color is white
        elseif row==2
            str_ = 'Particle label'; myColorMap = parula(256); myColorMap(1,:) = 1; % background color is white
%             % Randomize label
%             tmpdata = data_Particle_label;
%             currentvalues = unique(data_Particle_label); tmp = currentvalues; tmp(tmp==0)=[]; tmp=1:1:length(tmp);
%             for k=1:1:length(currentvalues)
%                 if currentvalues(k)~=0
%                     idx = randi(length(tmp));
%                     tmpdata( data_Particle_label==currentvalues(k) ) = tmp(idx);
%                     tmp(idx)=[];
%                 end
%             end
        else
            str_ = 'Boundary lines';
        end
        for current_direction=1:1:number_dimension % Iterate over axe
            kk=kk+1;
            sub_axes=subplot(3,number_dimension,kk,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            h_title=title ({[str_ ', slice in the middle'],['View normal to ' INFO.direction(current_direction).name]}); % Set title font
            % - Plot graphs
            if current_direction==1
                if row==3
                    tmpdata = volume_color(:,:,:,1);
                    r_ = squeeze(tmpdata(round(Domain_size(1)/2),:,:));
                    tmpdata = volume_color(:,:,:,2);
                    g_ = squeeze(tmpdata(round(Domain_size(1)/2),:,:));
                    tmpdata = volume_color(:,:,:,3);
                    b_ = squeeze(tmpdata(round(Domain_size(1)/2),:,:));
                    tmp = r_; tmp(:,:,2) = g_; tmp(:,:,3) = b_; 
                else
                    tmp=squeeze(tmpdata(round(Domain_size(1)/2),:,:));
                end
                h=image(tmp,'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',INFO.direction(3).name);
                t_1y = sprintf('Position along %s ',INFO.direction(2).name);
                set(h, 'XData', [0, Domain_size(3)*voxel_size/1000]);
                set(h, 'YData', [0, Domain_size(2)*voxel_size/1000]);
            elseif current_direction==2
                if row==3
                    tmpdata = volume_color(:,:,:,1);
                    r_ = squeeze(tmpdata(:,round(Domain_size(2)/2),:));
                    tmpdata = volume_color(:,:,:,2);
                    g_ = squeeze(tmpdata(:,round(Domain_size(2)/2),:));
                    tmpdata = volume_color(:,:,:,3);
                    b_ = squeeze(tmpdata(:,round(Domain_size(2)/2),:));
                    tmp = r_; tmp(:,:,2) = g_; tmp(:,:,3) = b_;                     
                else
                    tmp=squeeze(tmpdata(:,round(Domain_size(2)/2),:));
                end
                h=image(tmp,'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',INFO.direction(3).name);
                t_1y = sprintf('Position along %s ',INFO.direction(1).name);
                set(h, 'XData', [0, Domain_size(3)*voxel_size/1000]);
                set(h, 'YData', [0, Domain_size(1)*voxel_size/1000]);
            elseif current_direction==3
                if row==3
                    r_ = volume_color(:,:,round(Domain_size(3)/2),1);
                    g_ = volume_color(:,:,round(Domain_size(3)/2),2);
                    b_ = volume_color(:,:,round(Domain_size(3)/2),3);
                    tmp = r_; tmp(:,:,2) = g_; tmp(:,:,3) = b_;  
                    h=image(tmp,'CDataMapping','scaled');
                else                
                    h=image(tmpdata(:,:,round(Domain_size(3)/2)),'CDataMapping','scaled');
                end
                t_1x = sprintf('Position along %s ',INFO.direction(1).name);
                t_1y = sprintf('Position along %s ',INFO.direction(2).name);
                set(h, 'XData', [0, Domain_size(1)*voxel_size/1000]);
                set(h, 'YData', [0, Domain_size(2)*voxel_size/1000]);
            end
            axis equal; axis tight;
            % - Axis label
            t_ = xlabel(' ');
            t_2 = '(\mum)';
            t_.String= [t_1x t_2]; % Sprintf does not accept greek characters
            t_ = ylabel(' ');
            t_2 = '(\mum)';
            t_.String= [t_1y t_2]; % Sprintf does not accept greek characters
            set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe-2); % Fontname and fontsize
            % Color map
            colormap(sub_axes,myColorMap);
            if row==1 % Create colorbar
                h=colorbar(sub_axes);
                ylabel(h, 'Diameter (\mum)');
                set(h,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe-2);
            end
            h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
    end
    sgtitle(Fig,['Particle diameter, label and boundary line, ' INFO.phase(current_phase).name],'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('View_PCRF_%s', INFO.phase(current_phase).filename);
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end

%%
%% PARTICLE ANALYSIS
%%
%% WARNING: THIS SECTION NEED REFACTORING !

% Require particle identification
for current_phase=1:1:number_phase

    disp ' '
    fprintf ('Deep analysis of the discrete particle size, for the phase: %s\n',INFO.phase(current_phase).name)    
    Particle_size = results_visualization(current_phase).Particle_diameter_PCRF; % Size
    Particle_label = results_visualization(current_phase).Particle_ParticleLabel_PCRF; % Label
    
    % Number of particle
    % Mass center
    % Particle ratio
    % Connectivity
    % Skeleton (particle mass-center to particle mass-center)
    % particle sphericity
    [results_correlation, results_visualization, Table_discrete_particlemorpholgy_analysis, unique_particle, mass_center_particle] = Function_analyse_discrete_particle_morphology(Particle_size, Particle_label, results_correlation, results_visualization, current_phase, Current_folder, INFO, OPTIONS, 'PCRF');
    Results_PCRF.Phase(current_phase).Table_discrete_particlemorpholgy_analysis = Table_discrete_particlemorpholgy_analysis;
        
    % Particle-Particle connectivity
    % Particle-Particle interface
    [results_correlation, results_visualization, Table_discrete_particleinterface_analysis, Connectivity_particle] = Function_analyse_discrete_particle_interface(Particle_size, Particle_label, unique_particle, mass_center_particle, results_correlation, results_visualization, current_phase, Current_folder, INFO, OPTIONS, 'PCRF');
    Results_PCRF.Phase(current_phase).Table_discrete_particleinterface_analysis = Table_discrete_particleinterface_analysis;

    % Particle network
    [results_correlation, results_visualization, Table_discrete_particle_graphnetowrk] = Function_analyse_discrete_particle_graphnetwork(Particle_size, Particle_label, unique_particle, mass_center_particle, Connectivity_particle, results_correlation, results_visualization, current_phase, Current_folder, INFO, OPTIONS, 'PCRF');
    Results_PCRF.Phase(current_phase).Table_discrete_particle_graphnetowrk = Table_discrete_particle_graphnetowrk;
    
    % Constriction factor
            
    
end


%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Diameter_PCRF_calculation_time';
[Results_PCRF] = function_time_figure(Time_measure,date_start, date_end, Results_PCRF, Current_folder, filename, 'Particle size (PCRF)', OPTIONS);
 
%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_particlediameter_PCRF'],'Results_PCRF')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Correlation_particlediameter_PCRF'],'results_correlation')
    Sub_folder = 'Visualization\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Visualization_particlediameter_PCRF'],'results_visualization','-v7.3');    
end

end
