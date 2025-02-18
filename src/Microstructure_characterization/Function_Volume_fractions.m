function [] = Function_Volume_fractions(Phase_microstructure,infovol,opt,p)
% Calculate volume fractions
% Function_Volume_fractions(Phase_microstructure, infovol, opt, p) - when used with the MATBOX toolbox
% or
% Function_Volume_fractions(Phase_microstructure, voxelsize, unit) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       e.g.: Function_Volume_fractions(<your_3d_array>, 0.4, 'um');

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 3 % Case for function called as: Function_Volume_fractions(Phase_microstructure, voxelsize, unit). Standalone use.
        voxelsize = infovol; clear infovol;
        unit = opt; clear opt;
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['Volumefractions_' char(t)];
        if ispc
            infovol.mainfolder=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
            separator = '\';
        else
            infovol.mainfolder = pwd;
            separator = '/';
        end
        infovol.volpath = [infovol.mainfolder separator infovol.volumesubfolder separator];    

        % Set default phase information
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        infovol.phaselabel = unique(Phase_microstructure);
        for k=1:1:length(infovol.phaselabel) % Loop over all unique values
            infovol.phasename(k,1) = {['Phase ' num2str(infovol.phaselabel(k))]}; % Assign phase name based on code value
            infovol.phasecolor(k,:) = allcolor(k,:);
            p.todo(k)=1;
        end
        infovol.voxelsize = voxelsize;
        infovol.unit = unit;
        
        % Set default direction information
        infovol.directionname = {'in-plane direction 1';'in-plane direction 2';'through-plane direction'};
        
        % Set default options
        opt.save.mat = true;
        opt.save.xls = true;
        opt.save.savefig = true;
        opt.save.fig_infig = true;
        opt.save.fig_format = {'png'};
        
        % Set display options
        opt.format.fontname = 'Times New Roman';
        opt.format.axefontsize =  12;
        opt.format.legendfontsize =  10;
        opt.format.titlefontsize =  14;
        opt.format.sgtitlefontsize =  16;
        opt.format.linewidth = 2;
        opt.format.grid = 'on'; opt.format.minorgrid = 'on';
        opt.format.autoclosefig = false;
        
        % No Voxel size dependence analysis
        p.scaling = 1;
        % No Representative Volume Element analysis
        p.RVE.number_RVE = 0;
        % Fractal analysis
        p.fractal_bertei.todo = false;     
        p.fractal_boxcounting.todo = false;  

    else % Incorrect number of argument
        disp 'Error calling Function_Volume_fractions. Wrong number of argument.'
        help Function_Volume_fractions
        return
    end
    
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Volume_fractions' separator];
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% VOLUME INFORMATION
Domain_size = size(Phase_microstructure);
if length(Domain_size)==2 % Add third column with unit value
    Domain_size = [Domain_size 1];
    number_dimension = 2;
else
    number_dimension =3; % 3D case
end
number_phase = length(infovol.phaselabel); % Number of phase
voxel_number = prod(Domain_size); % Number of voxel
voxel_size = infovol.voxelsize;
voxel_unit = infovol.unit;

%% INITIALIZE RESULTS (USE FOR CORRELATION)
current_phase_todo = 0;
for current_phase=1:1:number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        results_correlation(current_phase_todo).name = infovol.phasename(current_phase,1);
        % Domain length is saved (only for the volume fraction file)
        results_correlation(current_phase_todo).domain_length = voxel_number^(1/number_dimension) * voxel_size;
        results_correlation(current_phase_todo).voxel_size = voxel_size;
    end
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

% Initialization (generic)
number_phase_todo = sum(p.todo); % How many phase are we going to analyse ?
Numbervoxel_phase=zeros(number_phase_todo,1); 
timedata = zeros(number_phase_todo+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_phase_todo+1,1);
phasename_todo = cell(number_phase_todo,1);
phaselabel = zeros(number_phase_todo,1);
% Initialization (algorithm-specific)
Volumefraction_phase=zeros(number_phase_todo,1);

current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);

        % % Algorithm
        phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
        Numbervoxel_phase(current_phase_todo,1)= sum(sum(sum(Phase_microstructure==phaselabel(current_phase_todo,1) )));
        Volumefraction_phase(current_phase_todo,1) = Numbervoxel_phase(current_phase_todo,1)/voxel_number; % Defintion of volume fraction

        % % Correlation
        results_correlation(current_phase_todo).volume_fraction = Volumefraction_phase(current_phase_todo,1); % Save for correlation

        % % Time
        timedata_domain(current_phase_todo+1,1) = infovol.phasename(current_phase,1);
        timedata(current_phase_todo+1,1) = Numbervoxel_phase(current_phase_todo,1);
        timedata(current_phase_todo+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
        timedata(current_phase_todo+1,3) = toc(time_clock_start_phase); % CPU elapsed time        
    end
end
% CPU and stopwatch time - end
timedata_domain(1,1) = {'Full volume'};
timedata(1,2) = cputime-time_cpu_start_volume; % CPU elapsed time
timedata(1,3) = toc(time_clock_start_volume); % Stopwatch elapsed time
timedata_pervolume = timedata(1,:);
timedata_perphase = timedata(2:end,:);

%% TABLES
% Time
Table_time = table(timedata_domain(:,1), timedata(:,1),timedata(:,2),timedata(:,3),...
    'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Volumefraction.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_Volumefraction = table(phasename_todo,phaselabel,Volumefraction_phase,Numbervoxel_phase,...
    'VariableNames',{'Phase' 'Label' 'Volume fraction' 'Number of voxel'});%
Results_Volumefraction.Table_Volumefraction = Table_Volumefraction; % Save in main table result

%% SAVE TABLES
if opt.save.xls
    filename = 'Volume_fraction'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Volume_fraction';
    DATA_writetable.sheet(1).table=Table_Volumefraction;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp(Table_Volumefraction)
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME
%%

%% ALONG DIRECTIONS
% Initialize
for current_direction=1:1:number_dimension % Loop over all directions
    tmp = Domain_size; tmp(current_direction)=1; Results.direction(current_direction).numbervoxelslice = prod(tmp); clear tmp; % Number of voxel within each slice normal to the direction
    Results.direction(current_direction).volumefraction = zeros(Domain_size(current_direction),number_phase_todo+1); % Initialization
    Results.direction(current_direction).volumefraction(:,1) = (1:1:Domain_size(current_direction))*voxel_size; % Position along direction
end
% Calculate results
current_phase_todo = 0;
for current_phase = 1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;

        label = infovol.phaselabel(current_phase);
        % Direction 1
        for current_position=1:1:Domain_size(1) % Loop over postion
            Numbervoxel_position= sum(sum(Phase_microstructure(current_position,:,:)==label));
            Results.direction(1).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(1).numbervoxelslice;
        end
        % Direction 2
        for current_position=1:1:Domain_size(2) % Loop over postion
            Numbervoxel_position= sum(sum(Phase_microstructure(:,current_position,:)==label));
            Results.direction(2).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(2).numbervoxelslice;
        end
        if number_dimension==3 % 3D case
            % Direction 3
            for current_position=1:1:Domain_size(3) % Loop over postion
                Numbervoxel_position= sum(sum(Phase_microstructure(:,:,current_position)==label));
                Results.direction(3).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(3).numbervoxelslice;
            end
        end
    end
end

%% TABLES
% Prepare header name
clear Variable_name_table;
Variable_name_table={['Position ' voxel_unit]};
for current_phase_todo=1:1:number_phase_todo
    Variable_name_table(current_phase_todo+1)=phasename_todo(current_phase_todo);
end
% Create table
for current_direction=1:1:number_dimension % Loop over all directions
    EvolutionVolumefraction.direction(current_direction).table = array2table(Results.direction(current_direction).volumefraction,'VariableNames',Variable_name_table);
end
Results_Volumefraction.EvolutionVolumefraction = EvolutionVolumefraction; % Save in main table result

%% SAVE TABLES
if opt.save.xls
    filename = 'Volume_fraction_along_directions'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    for current_direction=1:1:number_dimension % Loop over all directions
        DATA_writetable.sheet(current_direction).name = char(infovol.directionname(current_direction));
        DATA_writetable.sheet(current_direction).table = EvolutionVolumefraction.direction(current_direction).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% FIGURES
strunit = voxel_unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
else
    axisunit = ['(' strunit ')'];
end

scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Volume fractions'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
for current_direction=1:1:number_dimension % Iterate over axe
    sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    h_title=title(['Volume fractions along ',char(infovol.directionname(current_direction))]);
    % Plot graphs
    x_ = Results.direction(current_direction).volumefraction(:,1);
    current_phase_todo = 0;
    for current_phase=1:1:number_phase % Loop over phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            y_ = Results.direction(current_direction).volumefraction(:,current_phase_todo+1);
            plot(x_,y_,'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
        end
    end
    % Axis label
    t_ = xlabel(' ');
    t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
    t_2 = axisunit;
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    ylabel('Volume fractions \epsilon');
    % Legend
    for current_phase_todo=1:1:number_phase_todo
        str_legend(current_phase_todo).name = [char(phasename_todo(current_phase_todo)) ', <\epsilon>=' num2str(Volumefraction_phase(current_phase_todo),'%1.3f')];
    end
    h_legend = legend(sub_axes,str_legend.name,'Location','best');
    % - Grid
    grid(sub_axes,opt.format.grid); % Display grid
    set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
    set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
    h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
    hold(sub_axes,'off'); % Relase figure    
end
sgtitle(Fig,'Volume fractions along directions','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    filename= 'Volume_fractions_along_direction';
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

%% FRACTAL DIMENSION (COUNTING BOX)
if p.fractal_boxcounting.todo
    fprintf('> Phase domain fractal dimension (box-counting method)\n');
    p.fractal_boxcounting.topology_dimension = number_dimension;
    p.fractal_boxcounting.plot = false;
    for current_phase_todo=1:1:number_phase_todo % Loop over phase and call box counting algorithm
        binary_array = zeros(size(Phase_microstructure));
        binary_array(Phase_microstructure==phaselabel(current_phase_todo,1))=1;
        [N(:,current_phase_todo),box_lengths,fractal_dimension(:,current_phase_todo),fractal_dimension_convergence(:,:,current_phase_todo)] = Function_fractaldimension_boxcounting(binary_array,p.fractal_boxcounting);
    end

    % Table
    Table_fractaldimension = table(phasename_todo,fractal_dimension(1,:)',fractal_dimension(2,:)',fractal_dimension(3,:)',fractal_dimension(4,:)',...
        'VariableNames',{'Phase' 'Fit from 1 to' 'Fractal dimension' 'Topology dimension' 'Fractal propensity'});
    for current_phase_todo=1:1:number_phase_todo
        Table_boxlength(current_phase_todo).t = table(fractal_dimension_convergence(:,1,current_phase_todo),fractal_dimension_convergence(:,2,current_phase_todo),fractal_dimension_convergence(:,3,current_phase_todo),...
        'VariableNames',{'Fit from 1 to' 'Fractal dimension' 'Fit norm error'});
    end
    Results_Volumefraction.Table_fractaldimension = Table_fractaldimension; % Save in main table result
    disp(Table_fractaldimension)

    % Save table
    if opt.save.xls
        filename = 'PhaseDomain_Fractaldimension_boxcounting'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name = 'Fractal dimension';
        DATA_writetable.sheet(1).table = Table_fractaldimension;
        for current_phase_todo=1:1:number_phase_todo
            DATA_writetable.sheet(1+current_phase_todo).name = char(phasename_todo(current_phase_todo));
            DATA_writetable.sheet(1+current_phase_todo).table = Table_boxlength(current_phase_todo).t;
        end           
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

    % Figure
    if p.fractal_boxcounting.topology_dimension==3
        propertyname= 'Phase volume'; % Figure name
    else
        propertyname = 'Phase surface'; % Figure name
    end
    filename = 'PhaseDomain_fractal_dimension_boxcounting';
    function_fractalfig(propertyname,filename,Current_folder, box_lengths,N,fractal_dimension_convergence,p.fractal_boxcounting.topology_dimension,number_phase,p,opt,infovol);

    % Correlation
    results_correlation(current_phase_todo).PhaseDomain_fractaldimension_boxcounting = fractal_dimension(2,current_phase_todo);
    results_correlation(current_phase_todo).PhaseDomain_fractalpropensity_boxcounting = abs(p.fractal_boxcounting.topology_dimension - fractal_dimension(2,current_phase_todo));
end

%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%

if length(p.scaling)>=2 % Check if voxel size analysis is asked
    size_choice = p.scaling;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed
    
    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1);
    property_voxelsizedependence(1,1)=voxel_size;
    property_voxelsizedependence(1,2:end)=Volumefraction_phase(:,1)';

    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fractaldimension_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1);
        fractaldimension_voxelsizedependence(1,1) = voxel_size;
        fractaldimension_voxelsizedependence(1,2:end) = fractal_dimension(2,:);
    end

    % Loop on each voxel size
    for current_iteration=1:1:number_resize

        % % Microstructure scaling
        % New voxel size
        property_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
        % Set parameters
        parameters_scaling.scaling_factor = size_choice(current_iteration);
        parameters_scaling.label_or_greylevel = 'Label';
        parameters_scaling.background = min(infovol.phaselabel);
        % Scale
        Phase_microstructure_resized = function_scaling(Phase_microstructure,parameters_scaling);

        % CPU and stopwatch time - start
        time_cpu_start_volume = cputime; % CPU start
        time_clock_start_volume = tic; % Stopwatch start
        % Number of voxel of the current resized microstructure
        voxel_number_tmp=numel(Phase_microstructure_resized);

        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over all phases
            if p.todo(current_phase)
                time_cpu_start_phase = cputime; % CPU start
                time_clock_start_phase = tic; % Stopwatch start

                current_phase_todo=current_phase_todo+1;

                % % Algorithm: SPECIFIC FOR EACH FILE
                code_tmp = infovol.phaselabel(current_phase);
                Numbervoxel_phase_tmp= sum(sum(sum(Phase_microstructure_resized==code_tmp )));
                property_voxelsizedependence(current_iteration+1,current_phase_todo+1)=Numbervoxel_phase_tmp/voxel_number_tmp;

                % % Time
                timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

            end
        end
        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];

        % % Fractal
        if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
            for current_phase_todo=1:1:number_phase_todo % Loop over phase and call box counting algorithm
                binary_array = zeros(size(Phase_microstructure_resized));
                binary_array(Phase_microstructure_resized==phaselabel(current_phase_todo,1))=1;
                [~,~,fractal_dimension_tmp, ~] = Function_fractaldimension_boxcounting(binary_array,p.fractal_boxcounting);
                fractaldimension_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
                fractaldimension_voxelsizedependence(current_iteration+1,current_phase_todo+1) = fractal_dimension_tmp(2,1);
            end      
        end

    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    property_voxelsizedependence = sortrows(property_voxelsizedependence,1);   
    
    %% EXTRAPOLATION TO 0 nm
    fprintf('> Volume fractions dependence with the voxel size\n');
    tmp = zeros(number_resize+2,number_phase_todo+1); % + 0 nm and + initial voxel size
    x=property_voxelsizedependence(:,1);
    if strcmp(p.scaling_extrapolation,'Linear')
        interpolation_voxelsize_order=1;
    elseif strcmp(p.scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order=2;
    elseif strcmp(p.scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order = min(interpolation_voxelsize_order,max_order);
    fprintf('  Extrapolation to zero voxel size: polynomial of order %i\n\n',interpolation_voxelsize_order);
    for current_phase_todo=1:1:number_phase_todo
        y=property_voxelsizedependence(:,current_phase_todo+1);
        pi = polyfit(x,y,interpolation_voxelsize_order);
        vq = polyval(pi,0);
        tmp(1,current_phase_todo+1)=vq;
        interpolation_voxelsize(current_phase_todo).pi=pi;
    end
    tmp(2:end,:) = property_voxelsizedependence;
    property_voxelsizedependence = tmp; clear tmp;
    
    if p.fractal_bertei.todo
        fractal_dimension_Bertei = zeros(number_phase_todo,3);
        % Fractal dimension according to Bertei et al., https://doi.org/10.1016/j.nanoen.2017.06.028 (Richardson/Mandelbrot formula)
        % Log(property) = m + (3-fractal dimension)*log(voxel size)
        logV = log(property_voxelsizedependence(2:end,1));
        logPs = zeros(length(logV),number_phase_todo);
        for current_phase_todo=1:1:number_phase_todo
            logP = log(property_voxelsizedependence(2:end,current_phase_todo+1));
            logPs(:,current_phase_todo)=logP;
            pf = polyfit(logV,logP,1);
            fractal_dimension_Bertei(current_phase_todo,1) = 3-pf(1);
            fractal_dimension_Bertei(current_phase_todo,2) = 3; % Topology dimension, line
            fractal_dimension_Bertei(current_phase_todo,3) = abs( fractal_dimension_Bertei(current_phase_todo,2) - fractal_dimension_Bertei(current_phase_todo,1) ); % Fractal propensity
            % Correlation
            results_correlation(current_phase_todo).PhaseDomain_fractaldimension_Mandelbrot = fractal_dimension_Bertei(current_phase_todo,1);
            results_correlation(current_phase_todo).PhaseDomain_fractalpropensity_Mandelbrot = fractal_dimension_Bertei(current_phase_todo,3);
        end
    end

    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_phase_todo=1:1:number_phase_todo
        Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
    end
    % Table
    Table_Volumefraction_voxelsizedependence = array2table(property_voxelsizedependence,...
        'VariableNames',Variable_name_table);
    if p.fractal_bertei.todo
        Variable_name_table={'Voxel size log'}; % Columns name
        for current_phase_todo=1:1:number_phase_todo
            Variable_name_table(1+current_phase_todo)={[char(phasename_todo(current_phase_todo)) ' log']};
        end
        Table_Volumefraction_voxelsizedependence_loglog = array2table([logV logPs],'VariableNames',Variable_name_table);
        Table_Fractaldimension_Bertei = table(phasename_todo(:,1),fractal_dimension_Bertei(:,1),fractal_dimension_Bertei(:,2),fractal_dimension_Bertei(:,3),...
            'VariableNames',{'Phase','Fractal dimension','Topology dimension','Fractal propensity'});   
    end    
    
    %% DISPLAY TEXT RESULTS
    disp(Table_Volumefraction_voxelsizedependence)
    if p.fractal_bertei.todo
        disp(Table_Volumefraction_voxelsizedependence_loglog)
        fprintf('Richardson/Mandelbrot formula: Log(property) = m + (3-fractal dimension)*log(voxel size))\n');
        disp(Table_Fractaldimension_Bertei)
    end    
    
    %% SAVE RESULTS
    Results_Volumefraction.voxelsizedependence = Table_Volumefraction_voxelsizedependence; % Save in main table result
    if opt.save.xls
        filename = 'Volume_fraction_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Volume_fractions';
        DATA_writetable.sheet(1).table=Table_Volumefraction_voxelsizedependence;
        if p.fractal_bertei.todo
            DATA_writetable.sheet(2).name='Log log';
            DATA_writetable.sheet(2).table=Table_Volumefraction_voxelsizedependence_loglog;
            DATA_writetable.sheet(3).name='Fractal dimension Mandelbrot';
            DATA_writetable.sheet(3).table=Table_Fractaldimension_Bertei;
        end        
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    %% FIGURES
    parameters_figure.plotlog = true; 
    parameters_figure.propertyname = 'Volume fractions';
    parameters_figure.method = [];
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence;
    parameters_figure.number_phase = number_phase_todo;
    parameters_figure.str_ylabel = 'Volume fractions';
    parameters_figure.propertynameunit = [];
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'Volume_fractions_voxel_size_dependence';
    parameters_figure.infovol = infovol;
    parameters_figure.opt = opt;
    parameters_figure.todo = p.todo;
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    

    %% FRACTAL DIMENSION
    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fprintf('> Phase domain fractal dimension (box counting) dependence with the voxel size\n');
        fprintf('  Extrapolation to zero voxel size: polynomial of order 1\n\n');

        % Sort by voxel size
        fractaldimension_voxelsizedependence = sortrows(fractaldimension_voxelsizedependence,1);

        % Extrapolation to 0 voxel size
        tmp = zeros(number_resize+2,number_phase_todo+1); % + 0 nm and + initial voxel size
        for current_phase_todo=1:1:number_phase_todo
            y=fractaldimension_voxelsizedependence(:,current_phase_todo+1);
            pi = polyfit(x,y,1);
            vq = polyval(pi,0);
            tmp(1,current_phase_todo+1)=vq;
            interpolation_voxelsize(current_phase_todo).pi=pi;           
        end
        tmp(2:end,:) = fractaldimension_voxelsizedependence;
        fractaldimension_voxelsizedependence = tmp; clear tmp;

        % Managing result
        Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
        for current_phase_todo=1:1:number_phase_todo
            Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
        end
        % Table
        Table_Fractaldimension_voxelsizedependence = array2table(fractaldimension_voxelsizedependence,...
            'VariableNames',Variable_name_table);
        fprintf('    fitted from s=1 to s=2\n');
        disp(Table_Fractaldimension_voxelsizedependence)

        % Save result
        Results_Volumefraction.Fractaldimension_voxelsizedependence = Table_Fractaldimension_voxelsizedependence; % Save in main table result
        if opt.save.xls
            filename = 'PhaseDomain_Fractaldimension_boxcounting_voxel_size_dependence'; % Filename without extension
            % Prepare the data
            clear DATA_writetable
            DATA_writetable.sheet(1).name='Fit from 1 to 2';
            DATA_writetable.sheet(1).table=Table_Fractaldimension_voxelsizedependence;
            % Save function
            Function_Writetable(Current_folder,filename,DATA_writetable)
        end

        % Correlation
        for current_phase_todo=1:1:number_phase_todo
            results_correlation(current_phase_todo).PhaseDomain_fractaldimension_boxcounting_extrapolated = fractaldimension_voxelsizedependence(1,current_phase_todo+1) ;
            results_correlation(current_phase_todo).PhaseDomain_fractalpropensity_boxcounting_extrapolated = abs(p.fractal_boxcounting.topology_dimension - fractaldimension_voxelsizedependence(1,current_phase_todo+1));
        end

        % Figure
        parameters_figure.plotlog = false; 
        parameters_figure.figname = 'Phase domain fractal dimension';
        parameters_figure.propertyname = 'Phase domain fractal dimension';
        parameters_figure.method = 'Box counting';
        parameters_figure.property_voxelsizedependence = fractaldimension_voxelsizedependence;
        parameters_figure.number_phase = number_phase_todo;
        parameters_figure.str_ylabel = 'Fractal dimension';
        parameters_figure.propertynameunit = [];
        parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
        parameters_figure.Current_folder = Current_folder;
        parameters_figure.filename = 'Phase_Fractaldimension_voxel_size_dependence';
        parameters_figure.infovol = infovol;
        parameters_figure.opt = opt;
        parameters_figure.todo = p.todo;
        Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    end

end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

if p.RVE.number_RVE>0
    for k_RVE = 1:1:p.RVE.number_RVE % Loop over all RVE asked
        RVEparameters = p.RVE.RVE(k_RVE);
        RVE(k_RVE).RVEparameters = RVEparameters; % For result structure
        if opt.save.xls || opt.save.savefig
            Sub_folder_RVE = [Current_folder RVEparameters.savename separator];
            while exist(Sub_folder_RVE,'dir')
                RVEparameters.savename = [RVEparameters.savename '_bis'];
                Sub_folder_RVE = [Current_folder RVEparameters.savename separator];
            end
            mkdir(Sub_folder_RVE);
        end
        if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
            thresholds = RVEparameters.threshold_std;
        else
            thresholds = RVEparameters.threshold_reldiff;
        end
        n_threshold = length(thresholds);

        % Nested analysis ?
        if RVEparameters.donested % Yes
            [n_nestedRVE, xcrop] = Function_NestedRVE(RVEparameters,Domain_size);
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_phase_todo,3,3); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
        else % No
            n_nestedRVE = 0;
        end

        for k_nestedRVE = 0:1:n_nestedRVE
            %% NESTED ANALYSIS
            if k_nestedRVE == 0
                Domain_size_nested = Domain_size;
            else
                x0 = xcrop(k_nestedRVE,1); x1 = xcrop(k_nestedRVE,2);
                y0 = xcrop(k_nestedRVE,3); y1 = xcrop(k_nestedRVE,4);
                z0 = xcrop(k_nestedRVE,5); z1 = xcrop(k_nestedRVE,6);
                Phase_microstructure_nested = Phase_microstructure(x0:x1,y0:y1,z0:z1);
                Domain_size_nested = size(Phase_microstructure_nested);
            end
            
            %% SUBDOMAINS
            [All_subdomain,GROUP_SUBDOMAIN, Wholevolume_size] = Function_get_Subdomains(RVEparameters, Domain_size_nested, Domain_size); % Location of all subdomains
            Wholevolume_size(1:2) = Wholevolume_size(1:2)*voxel_size;
            if RVEparameters.donested
                Result_nestedRVE(k_nestedRVE+1,1,:,:,1) = Wholevolume_size(1);
                if k_nestedRVE ==1
                    fprintf('       Nested analysis\n');
                end
                if k_nestedRVE > 0
                    fprintf('          Cropping iteration: #%i/%i (cubic root volume=%1.1f %s)\n',k_nestedRVE,n_nestedRVE,Wholevolume_size(1),infovol.unit);
                end                
                if strcmp(RVEparameters.type,'C')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'D')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'G')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'H')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2) = Wholevolume_size(2);                    
                end
            end

            % Information about subdomains
            RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
                'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Equivalent cubic length' 'Equivalent section length' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
            [number_subdomain,~] = size(All_subdomain); % The number of subdomain
            number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size

            %% ALGORITHM
            % Colunm 1 is the subdomain id
            % Colunm 2 and 3 are the sizes of the subdomain.
            Property_eachsubdomain = zeros(number_subdomain,number_phase_todo+4);
            % Property calculated for each subdomain
            for subdomain_id = 1:1:number_subdomain
                % Boundary of the subdomain
                x0 = All_subdomain(subdomain_id,7); x1 = All_subdomain(subdomain_id,8);
                y0 = All_subdomain(subdomain_id,9); y1 = All_subdomain(subdomain_id,10);
                z0 = All_subdomain(subdomain_id,11); z1 = All_subdomain(subdomain_id,12);
                clear current_subdomain;
                % Crop volume
                if k_nestedRVE == 0
                    current_subdomain = Phase_microstructure(x0:x1,y0:y1,z0:z1);
                else
                    current_subdomain = Phase_microstructure_nested(x0:x1,y0:y1,z0:z1);
                end                    

                % CPU and stopwatch time - start
                time_cpu_start_volume = cputime; % CPU start
                time_clock_start_volume = tic; % Stopwatch start
                % Number of voxel of the current resized microstructure
                voxel_number_tmp=numel(current_subdomain);

                Property_eachsubdomain(subdomain_id,1)=subdomain_id;
                % Equivalent size of the subdomain
                Property_eachsubdomain(subdomain_id,2)=All_subdomain(subdomain_id,4)*voxel_size; % Cubic root length
                Property_eachsubdomain(subdomain_id,3)=All_subdomain(subdomain_id,5)*voxel_size; % Square root length
                Property_eachsubdomain(subdomain_id,4)=All_subdomain(subdomain_id,6)*voxel_size; % Length

                current_phase_todo = 0;
                for current_phase=1:1:number_phase % Loop over all phases
                    if p.todo(current_phase)
                        time_cpu_start_phase = cputime; % CPU start
                        time_clock_start_phase = tic; % Stopwatch start
                        current_phase_todo=current_phase_todo+1;

                        % % Algorithm: SPECIFIC FOR EACH FILE
                        code_tmp = infovol.phaselabel(current_phase);
                        Numbervoxel_phase_tmp= sum(sum(sum(current_subdomain==code_tmp )));
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4)=Numbervoxel_phase_tmp/voxel_number_tmp;

                        % % Time
                        timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];
                    end
                end
                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
            end

            %% STATISTICAL ANALYSIS and RVE SIZE
            [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_phase_todo,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters);
            
            %% SAVE FOR CORRELATION
            if k_nestedRVE == 0
                for k_threshold=1:1:n_threshold
                    current_phase_todo = 0;
                    for current_phase=1:1:number_phase
                        if p.todo(current_phase)
                            current_phase_todo=current_phase_todo+1;
                            if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                if Size_RVE(1,current_phase_todo,2,1)~=0
                                    str_ = ['vf_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,1);
                                end
                                if strcmp(RVEparameters.type,'C')
                                    if Size_RVE(1,current_phase_todo,2,2)~=0
                                        str_ = ['vf_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                    end
                                elseif strcmp(RVEparameters.type,'D')
                                    if Size_RVE(1,current_phase_todo,2,2)~=0
                                        str_ = ['vf_RVE_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                    end
                                end                                
                            else
                                if Size_convergence(1,current_phase_todo,2,1)~=0
                                    str_ = ['vf_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,1);
                                end
                                if strcmp(RVEparameters.type,'G')
                                    if Size_convergence(1,current_phase_todo,2,2)~=0
                                        str_ = ['vf_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
                                    end
                                elseif strcmp(RVEparameters.type,'H')
                                    if Size_convergence(1,current_phase_todo,2,2)~=0
                                        str_ = ['vf_conv_areaFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
                                    end
                                end

                            end
                        end
                    end
                end
            end

            %% MANAGING RESULTS
            [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,RVE,k_RVE,number_phase,number_phase_todo,infovol,p);

            if RVEparameters.donested
                for k_threshold=1:1:n_threshold
                    current_phase_todo = 0;
                    for current_phase=1:1:number_phase
                        if p.todo(current_phase)
                            current_phase_todo=current_phase_todo+1;
                            if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1) = Size_RVE(k_threshold,current_phase_todo,1,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1) = Size_RVE(k_threshold,current_phase_todo,2,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1) = Size_RVE(k_threshold,current_phase_todo,3,1);
                            else
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1) = Size_convergence(k_threshold,current_phase_todo,1,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1) = Size_convergence(k_threshold,current_phase_todo,2,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1) = Size_convergence(k_threshold,current_phase_todo,3,1);
                            end
                            if strcmp(RVEparameters.type,'C')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2) = Size_RVE(k_threshold,current_phase_todo,3,2);
                            elseif strcmp(RVEparameters.type,'D')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3) = Size_RVE(k_threshold,current_phase_todo,3,2);                                
                            elseif strcmp(RVEparameters.type,'G')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3) = Size_convergence(k_threshold,current_phase_todo,3,2);
                            elseif strcmp(RVEparameters.type,'H')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2) = Size_convergence(k_threshold,current_phase_todo,3,2);                                
                            end
                        end
                    end
                end
            end

            %% TEXT DISPLAY AND SAVE RESULTS
            if k_nestedRVE == 0
                propertyname='Volume fractions';
                RVEparameters.disp_parRVE = true;
                Function_subdomains_display_and_save(RVE,k_RVE,RVEparameters,number_phase,number_phase_todo,propertyname,Sub_folder_RVE,opt,infovol,p);
            end

            %% FIGURES
            if k_nestedRVE == 0
                parameters_figure.propertyname = propertyname;
                parameters_figure.propertynameunit = [];
                parameters_figure.RVE = RVEparameters;
                parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.number_phase = number_phase;
                parameters_figure.number_phase_todo = number_phase_todo;
                parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
                parameters_figure.derivative_convergence = derivative_convergence;
                parameters_figure.relativedifference_convergence = relativedifference_convergence;
                parameters_figure.Size_RVE = Size_RVE;
                parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
                parameters_figure.Size_convergence = Size_convergence;
                parameters_figure.Wholevolume_size = Wholevolume_size;
                parameters_figure.Wholevolume_results = Volumefraction_phase;
                parameters_figure.infovol = infovol;
                parameters_figure.todo = p.todo;
                parameters_figure.opt = opt;
                Function_create_figures_RVE(parameters_figure) % Figures
            end
        end

        %% NESTED ANALYSIS RESULT
        if RVEparameters.donested
            % Table
            [RVE] = Function_nestedtable(RVE,k_RVE,RVEparameters,number_phase,propertyname,Result_nestedRVE,Sub_folder_RVE,opt,infovol,p);
            % Figure
            parameters_figure.Result_nestedRVE = Result_nestedRVE;
            Function_create_figures_nestedRVE(parameters_figure) % Figures
            % Save
            RVE(k_RVE).nestedanalysis = Result_nestedRVE;
        end        
    end    
    Results_Volumefraction.RVE.volumefractions = RVE; % Save in main table result
end

%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Volumefraction.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Volumefraction.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_Volumefraction.Table_date = Table_date;

if opt.save.xls
    % Prepare the data
    clear DATA_writetable
    % Time per volume
    DATA_writetable.sheet(1).name='Time per volume';
    DATA_writetable.sheet(1).table=Table_time_pervolume;
    % Time per phase
    DATA_writetable.sheet(2).name='Time per phase';
    DATA_writetable.sheet(2).table=Table_time_perphase;    
    % Data : Date
    DATA_writetable.sheet(3).name='Date';
    DATA_writetable.sheet(3).table=Table_date;
    % Save function
    Function_Writetable(Current_folder,'Volume_fraction_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Volume_fraction_calculation_time', 'Volume fractions', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_volume_fraction'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_volume_fraction'],'Results_Volumefraction')
end

end