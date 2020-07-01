function [] = Function_Volume_fractions(Phase_microstructure, PROPERTY, OPTIONS, INFO )
%Calculate volume fractions
% Function_Volume_fractions(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_Volume_fractions(array, voxelsize) - when use as a standalone function

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 2 % Case for function called as: Function_Volume_fractions(Phase_microstructure, voxel_size)
        voxel_size = PROPERTY;
        clear PROPERTY
        
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
        PROPERTY.volumefractions.voxel_size_dependence.todo = false;
        PROPERTY.volumefractions.number_RVE = 0;
        
    else % Incorrect number of argument
        disp 'Error calling Function_Volume_fractions. Wrong number of argument.'
        help Function_Volume_fractions
    end
    
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Volume_fractions\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Volume_fractions/'; % Result are saved in this subfolder
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

%% INITIALIZE RESULTS (USE FOR CORRELATION)
for current_phase=1:1:number_phase
    results_correlation(current_phase).name = INFO.phase(current_phase).name;
end

% Domain length is saved (only for the volume fraction file)
for current_phase=1:1:number_phase
    results_correlation(current_phase).domain_length = prod(Domain_size)^(1/number_dimension) * voxel_size/1000;
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% CALCULATION
time_cpu_start = cputime; % CPU start
tic; % Stopwatch start
Numbervoxel_phase=zeros(number_phase,1); % Initialization
Volumefraction_phase=zeros(number_phase,1);
for current_phase=1:1:number_phase % Loop over all phases
    code_tmp = INFO.phase(current_phase).code;
    Numbervoxel_phase(current_phase,1)= sum(sum(sum(Phase_microstructure==code_tmp)));
    Volumefraction_phase(current_phase,1) = Numbervoxel_phase(current_phase,1)/voxel_number; % Defintion of volume fraction
    results_correlation(current_phase).volume_fraction = Volumefraction_phase(current_phase,1);
end
% CPU and stopwatch time - end
time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc; % Stopwatch elapsed time


%% TABLES
% Time
Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
Table_time = table(Time_measure(1)*1e-6,Time_measure(2),Time_measure(3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_Volumefraction.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_Volumefraction = table(INFO.phasename,[INFO.phase(:).code]',Volumefraction_phase,Numbervoxel_phase,...
    'VariableNames',{'Name' 'Code' 'Volume_fraction' 'Number_of_voxel'});%
Results_Volumefraction.Table_Volumefraction = Table_Volumefraction; % Save in main table result


%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Volume_fraction'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Volume_fraction';
    DATA_writetable.sheet(1).table=Table_Volumefraction;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


%% DISPLAY
if OPTIONS.displaytext==true
    disp '    VOLUME FRACTIONS';
    disp '    ----------------';
    disp ' ';
    fprintf('> Calculated on the whole domain:\n\n');
    disp(Table_Volumefraction)
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_time)
end


%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME
%%

%% ALONG DIRECTIONS
% Initialize
for current_direction=1:1:number_dimension % Loop over all directions
    tmp = Domain_size; tmp(current_direction)=1; Results.direction(current_direction).numbervoxelslice = prod(tmp); clear tmp; % Number of voxel within each slice normal to the direction
    Results.direction(current_direction).volumefraction = zeros(Domain_size(current_direction),number_phase+1); % Initialization
    Results.direction(current_direction).volumefraction(:,1) = (1:1:Domain_size(current_direction))*voxel_size/1000; % Axis (position along direction)
end
% Calculate results
for current_phase = 1:1:number_phase % Loop over all phases
    code_tmp = INFO.phase(current_phase).code; % Phase code
    % Direction 1
    for current_position=1:1:Domain_size(1) % Loop over postion
        Numbervoxel_position= sum(sum(Phase_microstructure(current_position,:,:)==code_tmp));
        Results.direction(1).volumefraction(current_position,current_phase+1) = Numbervoxel_position / Results.direction(1).numbervoxelslice;
    end
    % Direction 2
    for current_position=1:1:Domain_size(2) % Loop over postion
        Numbervoxel_position= sum(sum(Phase_microstructure(:,current_position,:)==code_tmp));
        Results.direction(2).volumefraction(current_position,current_phase+1) = Numbervoxel_position / Results.direction(2).numbervoxelslice;
    end
    if number_dimension==3 % 3D case
        % Direction 3
        for current_position=1:1:Domain_size(3) % Loop over postion
            Numbervoxel_position= sum(sum(Phase_microstructure(:,:,current_position)==code_tmp));
            Results.direction(3).volumefraction(current_position,current_phase+1) = Numbervoxel_position / Results.direction(3).numbervoxelslice;
        end
    end
end

%% TABLES

% Prepare header name
clear Variable_name_table;
Variable_name_table(1)={'Position_um'};
for current_phase=1:1:number_phase
    Variable_name_table(current_phase+1)={INFO.phase(current_phase).filename};
end
% Create table
for current_direction=1:1:number_dimension % Loop over all directions
    EvolutionVolumefraction.direction(current_direction).table = array2table(Results.direction(current_direction).volumefraction,'VariableNames',Variable_name_table);
end
Results_Volumefraction.EvolutionVolumefraction = EvolutionVolumefraction; % Save in main table result

%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Volume_fraction_along_directions'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    for current_direction=1:1:number_dimension % Loop over all directions
        DATA_writetable.sheet(current_direction).name = INFO.direction(current_direction).filename;
        DATA_writetable.sheet(current_direction).table = EvolutionVolumefraction.direction(current_direction).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% FIGURES
scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Volume fractions'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
for current_direction=1:1:number_dimension % Iterate over axe
    sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    h_title=title (' '); % Set title font
    h_title.String = sprintf('Volume fractions along %s',INFO.direction(current_direction).name); % Set title font
    % Plot graphs
    for current_phase=1:1:number_phase % Loop over phases
        x_ = Results.direction(current_direction).volumefraction(:,1);
        y_ = Results.direction(current_direction).volumefraction(:,current_phase+1);
        h=plot(x_,y_);
        set(h, 'Color', INFO.phase(current_phase).color, 'LineWidth',OPTIONS.Linewidth); % Colors
    end
    % Axis label
    t_ = xlabel(' ');
    t_1 = sprintf('Position along %s ',INFO.direction(current_direction).name);
    t_2 = '(\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    ylabel('Volume fractions');
    % Legend
    for current_phase=1:1:number_phase
        str_legend(current_phase).name = [INFO.phase(current_phase).name ', ' num2str(Volumefraction_phase(current_phase),'%1.3f')];
    end
    h_legend = legend(sub_axes,str_legend.name,'Location','best');
    % - Grid
    if strcmp(OPTIONS.grid,'on')
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
    end    
    set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
    h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
    h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
    hold(sub_axes,'off'); % Relase figure    
end
sgtitle(Fig,'Volume fractions along directions','FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
if OPTIONS.save_fig == true % Save figure
    filename= 'Volume_fractions_along_direction';
    function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end


%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%
interpolation_voxelsize_order=1;

if PROPERTY.volumefractions.voxel_size_dependence.todo % Check if voxel size analysis is asked
    size_choice = PROPERTY.volumefractions.voxel_size_dependence.voxel;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed
    
    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_phase+1);
    property_voxelsizedependence(1,1)=voxel_size;
    property_voxelsizedependence(1,2:end)=Volumefraction_phase(:,1)';
    % Loop on each voxel size
    for current_iteration=1:1:number_resize
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1)=current_voxel_size;
        % Microstructure resized
        [Phase_microstructure_resized] = function_scale_array(Phase_microstructure, voxel_size, current_voxel_size, INFO.phaseinfo);
        % CPU and stopwatch time - start
        time_cpu_start = cputime;
        tic;
        % Number of voxel of the current resized microstructure
        voxel_number_tmp=numel(Phase_microstructure_resized);
        % Volume fraction algorithm
        for current_phase=1:1:number_phase
            code_tmp = INFO.phase(current_phase).code;
            Numbervoxel_phase_tmp= sum(sum(sum(Phase_microstructure_resized==code_tmp)));
            property_voxelsizedependence(current_iteration+1,current_phase+1)=Numbervoxel_phase_tmp/voxel_number_tmp;
        end
        % CPU and stopwatch time - end
        time_cpu_elapsed = cputime-time_cpu_start;
        time_stopwatch_elapsed = toc;
        Time_tmp = [voxel_number_tmp time_cpu_elapsed time_stopwatch_elapsed];
        Time_measure = [Time_measure;Time_tmp];
    end
    clear Phase_microstructure_resized;
    
    %% EXTRAPOLATION TO 0 nm
    tmp = zeros(number_resize+2,number_phase+1);
    x=property_voxelsizedependence(:,1);
    for current_phase=1:1:number_phase
        y=property_voxelsizedependence(:,current_phase+1);
        p = polyfit(x,y,interpolation_voxelsize_order);
        vq = polyval(p,0);
        tmp(1,current_phase+1)=vq;
        interpolation_voxelsize(current_phase).p=p;
    end
    tmp(2:end,:) = property_voxelsizedependence;
    property_voxelsizedependence = tmp; clear tmp;
    
    
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={'Voxel_size_nm'}; % Columns name
    for current_phase=1:1:number_phase
        Variable_name_table(1+current_phase)={INFO.phase(current_phase).filename};
    end
    % Table
    Table_Volumefraction_voxelsizedependence = array2table(property_voxelsizedependence,...
        'VariableNames',Variable_name_table);
    
    %% DISPLAY TEXT RESULTS
    if (OPTIONS.displaytext==1)
        fprintf('> Volume fractions dependence with the voxel size (value at 0 nm has been linearly extrapolated):\n\n');
        disp(Table_Volumefraction_voxelsizedependence)
    end
    
    %% SAVE RESULTS
    Results_Volumefraction.voxelsizedependence = Table_Volumefraction_voxelsizedependence; % Save in main table result
    if OPTIONS.save_xls==true
        filename = 'Volume_fraction_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Volume_fractions';
        DATA_writetable.sheet(1).table=Table_Volumefraction_voxelsizedependence;
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    %% FIGURES
    parameters_figure.propertyname = 'Volume fractions';
    parameters_figure.method = [];
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence;
    parameters_figure.number_phase = number_phase;
    parameters_figure.str_ylabel = 'Volume fractions';
    parameters_figure.propertynameunit = [];
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'Volume_fractions_voxel_size_dependence';
    parameters_figure.INFO = INFO;
    parameters_figure.OPTIONS = OPTIONS;
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    
end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) ANALYSIS
%%

if PROPERTY.volumefractions.number_RVE>0
    for n_RVE=1:1:PROPERTY.volumefractions.number_RVE % Loop over all RVE asked
        RVEparameters.name = PROPERTY.volumefractions.RVE(n_RVE).name;
        RVEparameters.savename = PROPERTY.volumefractions.RVE(n_RVE).savename;
        RVEparameters.type = PROPERTY.volumefractions.RVE(n_RVE).type;
        RVEparameters.divisions = PROPERTY.volumefractions.RVE(n_RVE).divisions;
        RVEparameters.subs2 = PROPERTY.volumefractions.RVE(n_RVE).subs2;
        RVEparameters.subs4 = PROPERTY.volumefractions.RVE(n_RVE).subs4;
        RVEparameters.Aspectratio = PROPERTY.volumefractions.RVE(n_RVE).Aspectratio;
        if  strcmp(PROPERTY.volumefractions.RVE(n_RVE).type,'A')
            RVEparameters.Aspectratio_name = [num2str(Domain_size(1)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(2)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(3)/Domain_size(3),'%1.3f\t')];
        elseif strcmp(PROPERTY.volumefractions.RVE(n_RVE).type,'B') || strcmp(PROPERTY.volumefractions.RVE(n_RVE).type,'D')
            RVEparameters.Aspectratio_name = [num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
        end
        RVEparameters.Constantdirection = PROPERTY.volumefractions.RVE(n_RVE).Constantdirection;
        RVEparameters.Growthdirection = PROPERTY.volumefractions.RVE(n_RVE).Growthdirection;
        RVEparameters.Growthperstep = PROPERTY.volumefractions.RVE(n_RVE).Growthperstep;
        RVEparameters.Growthrelativeto = PROPERTY.volumefractions.RVE(n_RVE).Growthrelativeto;
        RVEparameters.threshold_std = PROPERTY.volumefractions.RVE(n_RVE).threshold_std;
        RVEparameters.threshold_numbersubvolumes = PROPERTY.volumefractions.RVE(n_RVE).threshold_numbersubvolumes;
        RVEparameters.firstuniquevolume_size = PROPERTY.volumefractions.RVE(n_RVE).firstuniquevolume_size;
        RVEparameters.firstuniquevolume_unit = PROPERTY.volumefractions.RVE(n_RVE).firstuniquevolume_unit;
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
            current_voxel_number = numel(current_subdomain);
            Property_eachsubdomain(subdomain_id,1)=subdomain_id;
            % Equivalent size of the subdomain
            Property_eachsubdomain(subdomain_id,2)=All_subdomain(subdomain_id,4)*voxel_size/1000; % Size are set in micrometer
            Property_eachsubdomain(subdomain_id,3)=All_subdomain(subdomain_id,5)*voxel_size/1000;
            % CPU and stopwatch time - start
            time_cpu_start = cputime;
            tic;
            for current_phase=1:1:number_phase
                code_tmp = INFO.phase(current_phase).code;
                Numbervoxel_current_subdomain_current_phase = sum(sum(sum(current_subdomain==code_tmp)));
                Volumefraction_current_subdomain_current_phase = Numbervoxel_current_subdomain_current_phase/current_voxel_number;
                Property_eachsubdomain(subdomain_id,current_phase+3)=Volumefraction_current_subdomain_current_phase;
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
                str_ = ['volume_fractions_RVE_' RVEparameters.savename];
                results_correlation(current_phase).(str_) = Size_RVE(2,current_phase,1);
            end
            if strcmp(RVEparameters.type,'C')
                if max(Size_RVE(2,:,2))~=0
                    str_ = ['volume_fractions_RVE_length_' RVEparameters.savename];
                    results_correlation(current_phase).(str_) = Size_RVE(2,current_phase,2);
                end
            end
        end        
        
        %% MANAGING RESULTS
        [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, RVEparameters,RVE,Size_RVE,n_RVE,number_phase,INFO);        
        
        %% TEXT DISPLAY AND SAVE RESULTS
        propertyname='Volume fractions';
        Function_subdomains_display_and_save(OPTIONS,INFO,RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)

        %% FIGURES
        parameters_figure.propertyname = propertyname;
        parameters_figure.propertynameunit = [];
        parameters_figure.RVE = RVEparameters;
        parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
        parameters_figure.savefolder = Sub_folder_RVE;
        parameters_figure.OPTIONS = OPTIONS;
        parameters_figure.INFO = INFO;
        parameters_figure.number_phase = number_phase;
        parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
        parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
        parameters_figure.Size_RVE = Size_RVE;
        parameters_figure.Wholevolume_results = Volumefraction_phase;
        parameters_figure.Wholevolume_size = Wholevolume_size;
        Function_create_figures_RVE(parameters_figure) % Figures
        
    end
    Results_Volumefraction.RVE.volumefractions = RVE; % Save in main table result
    
end

%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Volume_fraction_calculation_time';
[Results_Volumefraction] = function_time_figure(Time_measure,date_start, date_end, Results_Volumefraction, Current_folder, filename, 'Volume fractions', OPTIONS);


%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_volume_fraction'],'Results_Volumefraction')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end    
    save([Current_folder 'Correlation_volume_fraction'],'results_correlation')
end

end