function [] = Function_Specificsurface_direct(Phase_microstructure, PROPERTY, OPTIONS, INFO)
%Calculate specific surface area
% Function_Specificsurface_direct(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_Specificsurface_direct(array, voxelsize, corrective_factor) - when use as a standalone function

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 3 % Case for function called as: Function_Specificsurface_direct(Phase_microstructure, voxel_size, corrective_factor)
        voxel_size = PROPERTY;
        corrective_factor_specificsurface = OPTIONS;
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
        PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.todo = false;
        PROPERTY.specificsurfacearea_directmethod.number_RVE = 0;
        
        % Parameter
        PROPERTY.specificsurfacearea_directmethod.correctivefactor = corrective_factor_specificsurface;
        
    else % Incorrect number of argument
        disp 'Error calling Function_Specificsurface_direct. Wrong number of argument.'
        help Function_Specificsurface_direct
    end
    
end


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Specific_surface_direct\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Specific_surface_direct/'; % Result are saved in this subfolder
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

%% DISPLAY
if OPTIONS.displaytext==true
    disp '    SPECIFIC SURFACE AREA - DIRECT METHOD';
    disp '    -------------------------------------';
    disp ' ';
end

%% PARAMETERS
% Correct voxel cubic representation
corrective_factor_specificsurface = PROPERTY.specificsurfacearea_directmethod.correctivefactor;
% Phase surface / phase volume or Phase surface / domain volume ?
definition_sp = PROPERTY.specificsurfacearea_directmethod.definition;
if strcmp(definition_sp,'Phase surface divided by phase volume')
    refvolume = 2;
    refvalue = 4;
    str_figure = 'Phase_{surface} / Phase_{volume}';
    str_correlation = 'Sp_phase_surface_over_phase_volume_directcounting';
elseif strcmp(definition_sp,'Phase surface divided by domain volume')
    refvolume = 5;
    refvalue = 7;    
    str_figure = 'Phase_{surface} / Domain_{volume}';
    str_correlation = 'Sp_phase_surface_over_domain_volume_directcounting';
end


%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% CALCULATION
time_cpu_start = cputime; % CPU start
tic; % Stopwatch start
% Colunm 1: Surface area in square micrometer
% Column 2: Volume of the phase in cubic micrometer
% Column 3: Specific surface area (surface phase/volume phase) in micrometer-1
% Column 4: Specific surface area (surface phase/volume phase) * corrective_factor_specificsurface in micrometer-1
% Column 5: Volume of the domain in cubic micrometer
% Column 6: Specific surface area (surface phase/volume domain) in micrometer-1
% Column 7: Specific surface area (surface phase/volume domain) * corrective_factor_specificsurface in micrometer-1
Specificsurface_phase=zeros(number_phase,7); % Initialization
for current_phase=1:1:number_phase % Loop over all phases
    code_tmp = INFO.phase(current_phase).code;
    % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3)); % Initialization
    binary_phase(Phase_microstructure == code_tmp) = 1; % Binary phase
    n_=sum(sum(sum(binary_phase==1))); % Number of voxel of the phase
    Specificsurface_phase(current_phase,2)=n_*(voxel_size/1000)*(voxel_size/1000)*(voxel_size/1000); % Volume in cubic micrometer
    [Position_surface.phase(current_phase),surface_] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
    n_=surface_; clear surface_; % Number of voxel surface at the boundary of the phase
    Specificsurface_phase(current_phase,1)=n_*(voxel_size/1000)*(voxel_size/1000); % Surface area in square micrometer
    Specificsurface_phase(current_phase,3)=Specificsurface_phase(current_phase,1)/Specificsurface_phase(current_phase,2); % Specific surface area in micrometer-1
    Specificsurface_phase(current_phase,4)=Specificsurface_phase(current_phase,3)*corrective_factor_specificsurface; % Specific surface area * corrective factor in micrometer-1
    Specificsurface_phase(current_phase,5)=voxel_number*(voxel_size/1000)*(voxel_size/1000)*(voxel_size/1000); % Volume in cubic micrometer
    Specificsurface_phase(current_phase,6)=Specificsurface_phase(current_phase,1)/Specificsurface_phase(current_phase,5); % Specific surface area in micrometer-1
    Specificsurface_phase(current_phase,7)=Specificsurface_phase(current_phase,6)*corrective_factor_specificsurface; % Specific surface area * corrective factor in micrometer-1
    results_correlation(current_phase).Sp_phase_surface_over_phase_volume_directcounting = Specificsurface_phase(current_phase,4);
    results_correlation(current_phase).Sp_phase_surface_over_domain_volume_directcounting = Specificsurface_phase(current_phase,7);
end
% CPU and stopwatch time - end
time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc; % Stopwatch elapsed time

%% TABLES
% Time
Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
Table_time = table(Time_measure(1)*1e-6,Time_measure(2),Time_measure(3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_specificsurfacearea.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_Specificsurface_phasevolume_direct = table(INFO.phasename,[INFO.phase(:).code]',Specificsurface_phase(:,1),Specificsurface_phase(:,2),Specificsurface_phase(:,3),ones(number_phase,1)*corrective_factor_specificsurface,Specificsurface_phase(:,4),...
    'VariableNames',{'Name' 'Code' 'Surface_micrometer2' 'Volume_phase_micrometer3' 'Specific_surface_area_micrometer_m1' 'Corrective_factor' 'Corrected_Specific_surface_area_micrometer_m1'});
Results_specificsurfacearea.Table_Specificsurface_phasevolume_direct = Table_Specificsurface_phasevolume_direct; % Save in main table result
Table_Specificsurface_domainvolume_direct = table(INFO.phasename,[INFO.phase(:).code]',Specificsurface_phase(:,1),Specificsurface_phase(:,5),Specificsurface_phase(:,6),ones(number_phase,1)*corrective_factor_specificsurface,Specificsurface_phase(:,7),...
    'VariableNames',{'Name' 'Code' 'Surface_micrometer2' 'Volume_domain_micrometer3' 'Specific_surface_area_micrometer_m1' 'Corrective_factor' 'Corrected_Specific_surface_area_micrometer_m1'});
Results_specificsurfacearea.Table_Specificsurface_domainvolume_direct = Table_Specificsurface_domainvolume_direct; % Save in main table result


%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Specific_surface_area_direct'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Phasesurface_over_Phasevolume';
    DATA_writetable.sheet(1).table=Table_Specificsurface_phasevolume_direct;
    DATA_writetable.sheet(2).name='Phasesurface_over_Domainvolume';
    DATA_writetable.sheet(2).table=Table_Specificsurface_domainvolume_direct;    
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


%% DISPLAY
if OPTIONS.displaytext==true
    fprintf('> Calculated on the whole domain:\n\n');
    disp 'Specific surface area = phase surface / phase volume';
    disp(Table_Specificsurface_phasevolume_direct)
    disp 'Specific surface area = phase surface / domain volume';
    disp(Table_Specificsurface_domainvolume_direct)
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_time)
end


%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%

%% ALONG DIRECTIONS
% README

% Specific surface area Sp =  sum(g(x)), with g(x)=n(x)*face area/phase volume
% and n(x) the number of face at the position x along a given direction
% This discrete sum defintion will be converted in a continous (integral)
% sum definition: Sp = 1/D*integral(f(x),dx,0,D) = sum(g(x))
% To get f(x) [um-1], sum(g(x)) is derived.

% This is analogous with calculating a probability density function.
% (Although, the integral calculated here is not a probability but [um-1])

% Pros: obtain correct shape of the curve Sp=f(x).
% For a randomly distributed phase, f(x)=<Sp>

%% CALCULATION OF g(x), with g(x)=n(x)*face area/phase volume
for current_phase=1:1:number_phase
    for direction=1:1:3
        % Set position in micrometer
        Position_surface.phase(current_phase).direction(direction).value(:,1)=Position_surface.phase(current_phase).direction(direction).value(:,1)*voxel_size/1000;
        % Set face surface in squared micrometer
        Position_surface.phase(current_phase).direction(direction).value(:,3)=Position_surface.phase(current_phase).direction(direction).value(:,2)*(voxel_size/1000)*(voxel_size/1000);
        %  Set specific surface area
        Position_surface.phase(current_phase).direction(direction).value(:,3)=Position_surface.phase(current_phase).direction(direction).value(:,3)/Specificsurface_phase(current_phase,refvolume);
        % Corrective factor
        Position_surface.phase(current_phase).direction(direction).value(:,3)=Position_surface.phase(current_phase).direction(direction).value(:,3)*corrective_factor_specificsurface;
    end
end

%% DISTRIBUTION
% Cumulative function sum(g(x))
for current_phase=1:1:number_phase
    for direction=1:1:3
        Cumulative_surface.phase(current_phase).direction(direction).value=zeros(1,2*Domain_size(direction)-1);
        for i=1:1:2*Domain_size(direction)-1
            for j=1:1:i
                Cumulative_surface.phase(current_phase).direction(direction).value(1,j)=Cumulative_surface.phase(current_phase).direction(direction).value(1,j)+Position_surface.phase(current_phase).direction(direction).value(i,3);
            end
        end
    end
end
% Derivate function (i.e. 1/D*f(x))
for current_phase=1:1:number_phase
    for direction=1:1:3
        c_length = 2*Domain_size(direction)-1;
        for i=2:1:(c_length-1)
            Derivate_Cumulative_surface.phase(current_phase).direction(direction).value(1,i)=-( Cumulative_surface.phase(current_phase).direction(direction).value(1,i+1)-Cumulative_surface.phase(current_phase).direction(direction).value(1,i-1) ) / (2*0.5*voxel_size/1000);
        end
        Derivate_Cumulative_surface.phase(current_phase).direction(direction).value(1,1)=-(Cumulative_surface.phase(current_phase).direction(direction).value(1,2)-Cumulative_surface.phase(current_phase).direction(direction).value(1,1))/(0.5*voxel_size/1000);
        Derivate_Cumulative_surface.phase(current_phase).direction(direction).value(1,c_length)=-(Cumulative_surface.phase(current_phase).direction(direction).value(1,c_length)-Cumulative_surface.phase(current_phase).direction(direction).value(1,c_length-1))/(0.5*voxel_size/1000);
    end
end
% The specific surface area f(x)
for current_phase=1:1:number_phase
    for direction=1:1:3
        Specific_surface.phase(current_phase).direction(direction).value = zeros(1,2*Domain_size(direction)-1);
        Specific_surface.phase(current_phase).direction(direction).value = Derivate_Cumulative_surface.phase(current_phase).direction(direction).value*Domain_size(direction)*(voxel_size/1000);
    end
end

%% VERIFICATION
% the relation Sp = 1/D*integral(f(x),dx,0,D) is checked
Specificsurface_integralvalue=zeros(number_phase,7);
for current_phase=1:1:number_phase
    for direction=1:1:3
        X_= Position_surface.phase(current_phase).direction(direction).value(:,1);
        Y_= Specific_surface.phase(current_phase).direction(direction).value(1,:);
        Specificsurface_integralvalue(current_phase,1)=Specificsurface_phase(current_phase,refvalue); % Reference value
        Specificsurface_integralvalue(current_phase,direction+1)=(1/(Domain_size(direction)*(voxel_size/1000)))*trapz(X_,Y_); % Integral value
        % Relative error in percent
        Specificsurface_integralvalue(current_phase,direction+1+3)=100*(Specificsurface_integralvalue(current_phase,direction+1)-Specificsurface_integralvalue(current_phase,1))/Specificsurface_integralvalue(current_phase,1);
    end
end


%% MANAGING RESULTS
% Results are saved in a table

clear Variable_name_table;
Variable_name_table(1)={'Position_um'};
for current_phase=1:1:number_phase
    Variable_name_table(current_phase+1)={INFO.phase(current_phase).filename};
end

% g(x)
% Direction 1
array_g_dir1 = zeros(2*Domain_size(1)-1,number_phase+1);
array_g_dir1(:,1)= Position_surface.phase(current_phase).direction(1).value(:,1);
for current_phase=1:1:number_phase
    array_g_dir1(:,current_phase+1)= Position_surface.phase(current_phase).direction(1).value(:,3);
end
Table_evolution_g_direction_1 = array2table(array_g_dir1,...
    'VariableNames',Variable_name_table);
% Direction 2
array_g_dir2 = zeros(2*Domain_size(2)-1,number_phase+1);
array_g_dir2(:,1)= Position_surface.phase(current_phase).direction(2).value(:,1);
for current_phase=1:1:number_phase
    array_g_dir2(:,current_phase+1)= Position_surface.phase(current_phase).direction(2).value(:,3);
end
Table_evolution_g_direction_2 = array2table(array_g_dir2,...
    'VariableNames',Variable_name_table);
% Direction 3
array_g_dir3 = zeros(2*Domain_size(3)-1,number_phase+1);
array_g_dir3(:,1)= Position_surface.phase(current_phase).direction(3).value(:,1);
for current_phase=1:1:number_phase
    array_g_dir3(:,current_phase+1)= Position_surface.phase(current_phase).direction(3).value(:,3);
end
Table_evolution_g_direction_3 = array2table(array_g_dir3,...
    'VariableNames',Variable_name_table);

% f(x)
% Direction 1
array_f_dir1 = zeros(2*Domain_size(1)-1,number_phase+1);
array_f_dir1(:,1)= Position_surface.phase(current_phase).direction(1).value(:,1);
for current_phase=1:1:number_phase
    array_f_dir1(:,current_phase+1)= Specific_surface.phase(current_phase).direction(1).value(1,:);
end
Table_evolution_f_direction_1 = array2table(array_f_dir1,...
    'VariableNames',Variable_name_table);
% Direction 2
array_f_dir2 = zeros(2*Domain_size(2)-1,number_phase+1);
array_f_dir2(:,1)= Position_surface.phase(current_phase).direction(2).value(:,1);
for current_phase=1:1:number_phase
    array_f_dir2(:,current_phase+1)= Specific_surface.phase(current_phase).direction(2).value(1,:);
end
Table_evolution_f_direction_2 = array2table(array_f_dir2,...
    'VariableNames',Variable_name_table);
% Direction 3
array_f_dir3 = zeros(2*Domain_size(3)-1,number_phase+1);
array_f_dir3(:,1)= Position_surface.phase(current_phase).direction(3).value(:,1);
for current_phase=1:1:number_phase
    array_f_dir3(:,current_phase+1)= Specific_surface.phase(current_phase).direction(3).value(1,:);
end
Table_evolution_f_direction_3 = array2table(array_f_dir3,...
    'VariableNames',Variable_name_table);

% Verification
Table_verification_integral = table(INFO.phasename,Specificsurface_integralvalue(:,1),Specificsurface_integralvalue(:,2),Specificsurface_integralvalue(:,3),Specificsurface_integralvalue(:,4),Specificsurface_integralvalue(:,5),Specificsurface_integralvalue(:,6),Specificsurface_integralvalue(:,7),...
    'VariableNames',{'Name' 'Reference' 'Integration_direction_1' 'Integration_direction_2' 'Integration_direction_3' 'Relative_error_percent_direction_1' 'Relative_error_percent_direction_2' 'Relative_error_percent_direction_3'});
if OPTIONS.displaytext==true
    fprintf('> Calculated along direction:\n\n');
    fprintf('  Integral verification: relative error should be close to 0 percent:\n\n');
    disp(Table_verification_integral)
end

Results_specificsurfacearea.Table_evolution_g_direction_1 = Table_evolution_g_direction_1; % Save in main table result
Results_specificsurfacearea.Table_evolution_g_direction_2=Table_evolution_g_direction_2;
Results_specificsurfacearea.Table_evolution_g_direction_3=Table_evolution_g_direction_3;
Results_specificsurfacearea.Table_evolution_f_direction_1=Table_evolution_f_direction_1;
Results_specificsurfacearea.Table_evolution_f_direction_2=Table_evolution_f_direction_2;
Results_specificsurfacearea.Table_evolution_f_direction_3=Table_evolution_f_direction_3;
Results_specificsurfacearea.Table_verification_integral=Table_verification_integral;


%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Sp_direct_along_directions'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    % Data : Evolution of volume fraction along direction 1
    DATA_writetable.sheet(1).name='Along_direction_1_g';
    DATA_writetable.sheet(1).table=Table_evolution_g_direction_1;
    % Data : Evolution of volume fraction along direction 2
    DATA_writetable.sheet(2).name='Along_direction_2_g';
    DATA_writetable.sheet(2).table=Table_evolution_g_direction_2;
    % Data : Evolution of volume fraction along direction 3
    DATA_writetable.sheet(3).name='Along_direction_3_g';
    DATA_writetable.sheet(3).table=Table_evolution_g_direction_3;
    % Data : Evolution of volume fraction along direction 1
    DATA_writetable.sheet(4).name='Along_direction_1_f';
    DATA_writetable.sheet(4).table=Table_evolution_f_direction_1;
    % Data : Evolution of volume fraction along direction 2
    DATA_writetable.sheet(5).name='Along_direction_2_f';
    DATA_writetable.sheet(5).table=Table_evolution_f_direction_2;
    % Data : Evolution of volume fraction along direction 3
    DATA_writetable.sheet(6).name='Along_direction_3_f';
    DATA_writetable.sheet(6).table=Table_evolution_f_direction_3;    
    % Data : Integral summation
    DATA_writetable.sheet(7).name='Integration_f';
    DATA_writetable.sheet(7).table=Table_verification_integral;       
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end


%% FIGURES
scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Specific surface area, direct method (integral definition)'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
for current_direction=1:1:number_dimension % Iterate over axe
    sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    h_title=title ({'Specific surface area, direct method',sprintf('S_{p} = %s ',str_figure)}); % Set title font
    % Plot graphs
    for current_phase=1:1:number_phase % Loop over phases
        x_ = Position_surface.phase(current_phase).direction(current_direction).value(:,1);
        y_ = Specific_surface.phase(current_phase).direction(current_direction).value(1,:);
        h=plot(x_,y_);
        set(h, 'Color', INFO.phase(current_phase).color, 'LineWidth',OPTIONS.Linewidth); % Colors
    end
    % Axis label
    t_ = xlabel(' ');
    t_1 = sprintf('Position along %s ',INFO.direction(current_direction).name);
    t_2 = '(\mum)';
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' '); 
    str1 = 'Specific surface distribution function f(x), such as';
    str2 = 'S_{p} = 1/L \times \int_{0}^{L} f(x) dx (\mum^{-1})';
    t_.String= {str1,str2};
    % Legend
    for current_phase=1:1:number_phase
        str_legend(current_phase).name = [INFO.phase(current_phase).name ', ' num2str(Specificsurface_phase(current_phase,refvalue),'%1.3f') '\mum^{-1}'];
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
sgtitle(Fig,'Specific surface area distribution function','FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
if OPTIONS.save_fig == true % Save figure
    filename= 'Specific_surface_direct_along_direction';
    function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end


%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%
interpolation_voxelsize_order=1;

if PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.todo % Check if voxel size analysis is asked
    size_choice = PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed

    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_phase+1);
    property_voxelsizedependence(1,1)=voxel_size;
    property_voxelsizedependence(1,2:end)=Specificsurface_phase(:,refvalue)';
    % Loop on each voxel size
    for current_iteration=1:1:number_resize
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1)=current_voxel_size;
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
            if refvolume==2
                n_=sum(sum(sum(binary_phase==1))); % Number of voxel of the phase
            else
                n_=prod(Current_domain_size); % Number of voxel of the domain
            end
            Volume_=n_*(current_voxel_size/1000)*(current_voxel_size/1000)*(current_voxel_size/1000); % Volume in cubic micrometer
            [~,surface_] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
            surface_=surface_*(current_voxel_size/1000)*(current_voxel_size/1000); % Surface area in square micrometer
            Specificsurface_=surface_/Volume_; % Specific surface area in micrometer-1
            property_voxelsizedependence(current_iteration+1,current_phase+1)=Specificsurface_*corrective_factor_specificsurface; % Specific surface area * corrective factor in micrometer-1
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
    tmp = zeros(number_resize+2,number_phase+1);
    x=property_voxelsizedependence(:,1);
    for current_phase=1:1:number_phase
        y=property_voxelsizedependence(:,current_phase+1);
        p = polyfit(x,y,interpolation_voxelsize_order);
        vq = polyval(p,0);
        tmp(1,current_phase+1)=vq;
        interpolation_voxelsize(current_phase).p=p;
        % For correlation
        results_correlation(current_phase).([str_correlation '_extrapolated']) = vq;
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
    Table_Specficsurfacearea_direct_voxelsizedependence = array2table(property_voxelsizedependence,...
        'VariableNames',Variable_name_table);
    
    %% DISPLAY TEXT RESULTS
    if (OPTIONS.displaytext==1)
        fprintf('> Specific surface area dependence with the voxel size (with corrective factor):\n\n');
        disp(Table_Specficsurfacearea_direct_voxelsizedependence)
    end
       
    %% SAVE RESULTS
    Results_specificsurfacearea.voxelsizedependence = Table_Specficsurfacearea_direct_voxelsizedependence; % Save in main table result
    if OPTIONS.save_xls==true
        filename = 'Sp_direct_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Specific_surface_area_direct';
        DATA_writetable.sheet(1).table=Table_Specficsurfacearea_direct_voxelsizedependence;
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end        
        
    %% FIGURES
    parameters_figure.propertyname = 'Specific surface area';
    parameters_figure.method = 'direct counting';
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence;
    parameters_figure.number_phase = number_phase;
    parameters_figure.str_ylabel = [sprintf('S_{p} = %s ',str_figure) '(\mum^{-1})'];
    parameters_figure.propertynameunit = '\mum^{-1}';
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'Sp_direct_voxel_size_dependence';
    parameters_figure.INFO = INFO;
    parameters_figure.OPTIONS = OPTIONS;
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures
end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) ANALYSIS
%%

if PROPERTY.specificsurfacearea_directmethod.number_RVE>0
    for n_RVE=1:1:PROPERTY.specificsurfacearea_directmethod.number_RVE % Loop over all RVE asked
        RVEparameters.name = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).name;
        RVEparameters.savename = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).savename;
        RVEparameters.type = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).type;
        RVEparameters.divisions = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).divisions;
        RVEparameters.subs2 = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).subs2;
        RVEparameters.subs4 = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).subs4;
        RVEparameters.Aspectratio = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).Aspectratio;
        if  strcmp(PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).type,'A')
            RVEparameters.Aspectratio_name = [num2str(Domain_size(1)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(2)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(3)/Domain_size(3),'%1.3f\t')];
        elseif strcmp(PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).type,'B') || strcmp(PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).type,'D')
            RVEparameters.Aspectratio_name = [num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
        end
        RVEparameters.Constantdirection = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).Constantdirection;
        RVEparameters.Growthdirection = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).Growthdirection;
        RVEparameters.Growthperstep = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).Growthperstep;
        RVEparameters.Growthrelativeto = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).Growthrelativeto;
        RVEparameters.threshold_std = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).threshold_std;
        RVEparameters.threshold_numbersubvolumes = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).threshold_numbersubvolumes;
        RVEparameters.firstuniquevolume_size = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).firstuniquevolume_size;
        RVEparameters.firstuniquevolume_unit = PROPERTY.specificsurfacearea_directmethod.RVE(n_RVE).firstuniquevolume_unit;
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
                if refvolume==2
                    n_=sum(sum(sum(binary_phase==1))); % Number of voxel of the phase
                else
                    n_=prod(Current_domain_size); % Number of voxel of the domain
                end
                Volume_=n_*(voxel_size/1000)*(voxel_size/1000)*(voxel_size/1000); % Volume in cubic micrometer
                [~,surface_] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
                surface_=surface_*(voxel_size/1000)*(voxel_size/1000); % Surface area in square micrometer
                Specificsurface_=surface_/Volume_; % Specific surface area in micrometer-1
                Property_eachsubdomain(subdomain_id,current_phase+3)=Specificsurface_*corrective_factor_specificsurface; % Specific surface area * corrective factor in micrometer-1
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
                results_correlation(current_phase).(['Sp_direct_RVE_' RVEparameters.savename]) = Size_RVE(2,current_phase,1);
            end
            if strcmp(RVEparameters.type,'C')
                if max(Size_RVE(2,:,2))~=0
                    results_correlation(current_phase).(['Sp_direct_RVE_length_' RVEparameters.savename]) = Size_RVE(2,current_phase,2);
                end
            end
        end 
        
        %% MANAGING RESULTS
        [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, RVEparameters,RVE,Size_RVE,n_RVE,number_phase,INFO);        
        
        %% TEXT DISPLAY AND SAVE RESULTS
        propertyname='Specific surface area';
        Function_subdomains_display_and_save(OPTIONS,INFO,RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)
        
        %% FIGURES
        parameters_figure.propertyname = propertyname;
        parameters_figure.propertynameunit = '\mum^{-1}';
        parameters_figure.RVE = RVEparameters;
        parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
        parameters_figure.savefolder = Sub_folder_RVE;
        parameters_figure.OPTIONS = OPTIONS;
        parameters_figure.INFO = INFO;
        parameters_figure.number_phase = number_phase;
        parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
        parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
        parameters_figure.Size_RVE = Size_RVE;
        parameters_figure.Wholevolume_results = Specificsurface_phase(:,refvalue);
        parameters_figure.Wholevolume_size = Wholevolume_size;
        Function_create_figures_RVE(parameters_figure) % Figures        
                
    end
    Results_specificsurfacearea.RVE.sp = RVE; % Save in main table result
end


%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Sp_direct_calculation_time';
[Results_specificsurfacearea] = function_time_figure(Time_measure,date_start, date_end, Results_specificsurfacearea, Current_folder, filename, 'Specific surface area', OPTIONS);
 
%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_specific_surface_direct'],'Results_specificsurfacearea')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Correlation_specific_surface_area'],'results_correlation')
end


end