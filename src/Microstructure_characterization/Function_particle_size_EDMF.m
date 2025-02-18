function [] = Function_particle_size_EDMF(Phase_microstructure,infovol,opt,p,foo1,foo2)
% Calculate Particle size fitting the cumulative function of the euclidean distance map
% Method is detailed in MATBOX documentation and in journal article https://doi.org/10.1149/1945-7111/ab913b) 
% Function_particle_size_EDMF(Phase_microstructure,infovol,opt,p,foo) - when use with the toolbox
% or
% Function_particle_size_EDMF(Phase_microstructure, labels, voxelsize, unit, removeinclusion, distancelabel) - when use as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       labels: either 'All' (all phased are characterized) or a 1D array listing the labels to characterize
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       removeinclusion, true or false: inclusion are not considered for the distance calculation
%       distancelabel, an integer: distance from labels to label=distancelabel (instead of distance from labels to ~labels). To ignore, set to -1
%       e.g.: Function_particle_size_EDMF(<your_3d_array>, 1, 0.4, 'um',false,-1); % 1 within unique(<your_3d_array>)
%       e.g.: Function_particle_size_EDMF(<your_3d_array>, 'All' , 0.4, 'um',false,-1); % Calculate for all phases
%       e.g.: Function_particle_size_EDMF(<your_3d_array>, [0, 1, 2], 0.4, 'um',true,-1); % 0,1,2 within unique(<your_3d_array>)
%       e.g.: Function_particle_size_EDMF(<your_3d_array>, 1, 0.4, 'um',true,0); % 0,1,2 within unique(<your_3d_array>)

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 6 % Case for function called as: Function_particle_size_EDMF(Phase_microstructure, labels, voxelsize, unit, removeinclusion, distancelabel, applypercluster)
        labels = infovol; clear infovol;
        voxelsize = opt; clear opt;
        unit = p; clear p;
        removeinclusion = foo1; distancelabel = foo2; clear foo1; clear foo2;
               
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['DiameterEDMF_' char(t)];
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
            if isnumeric(labels)
                if ismember(infovol.phaselabel(k),labels)
                    p.todo(k)=1;
                else
                    p.todo(k)=0;
                end
            elseif strcmp(labels,'All')
                p.todo(k)=1;
            else
                disp 'Error calling Function_particle_size_EDMF. Incorrect labels argument.'
                help Function_particle_size_EDMF
                return
            end
        end
        if ~sum(p.todo)
            disp 'Error calling Function_particle_size_EDMF. Labels argument correspond to no existing label.'
            help Function_particle_size_EDMF
            return
        end

        % Directions investigated
        for k=1:1:3
            p.direction_todo(k)=0;
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

        % Set default colormap
        p.colormap = 'jet';
        p.colormap_background = 'white';
        
        % No Voxel size dependence analysis
        p.scaling = 1;
        % No Representative Volume Element analysis
        p.RVE.number_RVE = 0;
        
    else % Incorrect number of argument
        disp 'Error calling Function_particle_size_EDMF. Wrong number of argument.'
        help Function_particle_size_EDMF
        return
    end
    
else
    % Read parameters
    removeinclusion = p.removeinclusion;
    distancelabel = p.distancelabel;
end


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Diameter_EDMF' separator];
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

%% INITIALIZE RESULTS (USE FOR CORRELATION and VISUALIZATION)
current_phase_todo = 0;
for current_phase=1:1:number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        results_correlation(current_phase_todo).name = infovol.phasename(current_phase,1);
        results_visualization(current_phase_todo).name = infovol.phasename(current_phase,1);
    end
end

%% PARAMETERS
% Note: those are not algorithm paramter (c-PSD is parameter free), but post-processing smoothing parameter to get a better looking distribution function

% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;

%%
%% ALGORITHM ON WHOLE VOLUME
%%

disp '    PARTICLE SIZE - EUCLIDEAN DISTANCE MAP FITTING (EDMF) METHOD';
disp '    ------------------------------------------------------------';
disp ' ';

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
dmap_results = zeros(number_phase_todo,5); % min, mean, max, std, std%
d50_results = zeros(number_phase_todo,1); % mean

current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
        Numbervoxel_phase(current_phase_todo,1)= sum(sum(sum(Phase_microstructure==phaselabel(current_phase_todo,1) )));

        % % Algorithm
        phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
        % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
        binary_phase=zeros(Domain_size); % Initialization
        binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase
        [distance_transform, fitted_diameter, numericalpsd(current_phase_todo).psd, analyticalpsd(current_phase_todo).psd, error_radius(current_phase_todo).val] = Function_particle_size_distancemap_Algorithm(Phase_microstructure, phaselabel(current_phase_todo,1), removeinclusion, distancelabel, voxel_size, density_fct_parameters);
        % Equivalent mean diameter
        d50_results(current_phase_todo,1) = fitted_diameter;

        % Distance map statistics
        all_distance = distance_transform(binary_phase==1);

        dmap_results(current_phase_todo,1) = min(all_distance); % Minimum
        dmap_results(current_phase_todo,2) = mean(all_distance);% Mean
        dmap_results(current_phase_todo,3) = max(all_distance); % Maximum
        dmap_results(current_phase_todo,4) = std(all_distance); % Standard deviation
        dmap_results(current_phase_todo,5) = 100*dmap_results(current_phase_todo,4)/dmap_results(current_phase_todo,2); % Relative standard deviation
    
        % Cumulative and probability density distribution functions
        PSD_results(current_phase_todo,1) = numericalpsd(current_phase_todo).psd.x50; % d50
        PSD_results(current_phase_todo,2) = numericalpsd(current_phase_todo).psd.smoothed_x50; % smoothed d50
        PSD_results(current_phase_todo,3) = numericalpsd(current_phase_todo).psd.integral_probability_density_fct;
        PSD_results(current_phase_todo,4) = numericalpsd(current_phase_todo).psd.integral_smoothed_probability_density_fct;

        % % correlation
        results_correlation(current_phase_todo).Particle_diameter_mean_EDMF = d50_results(current_phase_todo,1);
        results_correlation(current_phase_todo).Distance_to_surface_mean = dmap_results(current_phase_todo,2);
        results_correlation(current_phase_todo).Distance_to_surface_std = dmap_results(current_phase_todo,4);
        results_correlation(current_phase_todo).Distance_to_surface_relstd = dmap_results(current_phase_todo,5);

        % % Visualization
        results_visualization(current_phase_todo).distance_to_boundary = distance_transform; 

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
Results_edmf.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_dmap = table(phasename_todo,phaselabel,dmap_results(:,1),dmap_results(:,2),dmap_results(:,3),dmap_results(:,4),dmap_results(:,5),...
    'VariableNames',{'Phase' 'Label' ['Min ' voxel_unit] ['Mean ' voxel_unit] ['Max ' voxel_unit] ['Std ' voxel_unit] 'Std %'});
Table_dmap_cumulative = table(phasename_todo, phaselabel,PSD_results(:,1),PSD_results(:,2),PSD_results(:,3),PSD_results(:,4),...
    'VariableNames',{'Phase' 'Label' ['d50 (from distribution fct) ' voxel_unit] ['Smoothed d50 ' voxel_unit] 'Distribution integral' 'Smoothed distribution integral'});
Table_edmf_d50 = table(phasename_todo, phaselabel, d50_results(:,1),...
    'VariableNames',{'Name' 'Label' ['Fitted diameter ' voxel_unit]});
Results_edmf.Table_dmap = Table_dmap; % Save in main table result
Results_edmf.Table_dmap_cumulative = Table_dmap_cumulative;
Results_edmf.Table_edmf_d50 = Table_edmf_d50;

for current_phase_todo=1:1:number_phase_todo
    array_tmp(:,1) = numericalpsd(current_phase_todo).psd.cumulative_fct(:,1);
    array_tmp(:,2) = numericalpsd(current_phase_todo).psd.cumulative_fct(:,2);
    array_tmp(:,3) = numericalpsd(current_phase_todo).psd.probability_density_fct(:,2);
    if ~isempty(numericalpsd(current_phase_todo).psd.smoothed_cumulative_fct)
        array_tmp(:,4) = numericalpsd(current_phase_todo).psd.smoothed_cumulative_fct(:,1);
        array_tmp(:,5) = numericalpsd(current_phase_todo).psd.smoothed_cumulative_fct(:,2);
        array_tmp(:,6) = numericalpsd(current_phase_todo).psd.smoothed_probability_density_fct(:,2);
        Variable_name_table={['Distance ' voxel_unit] 'Cumulative function' 'Probability density distribution function' ['Distance (smoothed) ' voxel_unit] 'Cumulative function smoothed' 'Probability density function smoothed'};
    else
        Variable_name_table={['Distance ' voxel_unit] 'Cumulative function' 'Probability density distribution function'};
    end
    Table_cumulative_sizedistribution.phase(current_phase_todo).table = array2table(array_tmp,...
        'VariableNames',Variable_name_table);
    clear array_tmp
end
Results_edmf.Table_cumulative_sizedistribution = Table_cumulative_sizedistribution; % Save in main table result

%% SAVE TABLES
if opt.save.xls
    filename = 'EDMF'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Distance_from_boundary';
    DATA_writetable.sheet(1).table=Table_dmap;
    DATA_writetable.sheet(2).name='From_cumulative_fct';
    DATA_writetable.sheet(2).table=Table_dmap_cumulative;        
    DATA_writetable.sheet(3).name='Fitted_diameter';
    DATA_writetable.sheet(3).table=Table_edmf_d50;       
    % Cumulative and size distribution
    sheet_=3;
    for current_phase_todo = 1:1:number_phase_todo
        sheet_=sheet_+1;
        DATA_writetable.sheet(sheet_).name=[char(phasename_todo(current_phase_todo,1)) '_DM'];
        DATA_writetable.sheet(sheet_).table=Table_cumulative_sizedistribution.phase(current_phase_todo).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp 'Distance map, i.e distance from boundary';
disp(Table_dmap)
disp 'Distance map, i.e distance from boundary (from cumulative function)';
disp(Table_dmap_cumulative)
disp 'Fitted diameter';
disp(Table_edmf_d50)
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%
scrsz = get(0,'ScreenSize'); % Screen resolution

strunit = voxel_unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
    Dunit = '\mum';
else
    axisunit = ['(' strunit ')'];
    Dunit = voxel_unit;
end

%% CUMULATIVE FUNCTION DIAMTER FITTING
current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;

        Fig= figure;
        Fig.Name= sprintf('Distance map: fitted diameter of phase %s',char(phasename_todo(current_phase_todo,1)));
        Fig.Color='white'; % Background colour
        set(Fig, 'Position', [scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*1/2]);

        sub_axes=subplot(1,2,1,'Parent',Fig); % Create axes
        hold(sub_axes,'on'); % Active subplot
        h_title=title (' ','FontName',opt.format.fontname);
        h_title.String= {'Analytical - numerical integral difference of','the distance to surface map cumulative fct'};
        h_error=plot(error_radius(current_phase_todo).val(:,1) ,error_radius(current_phase_todo).val(:,2)); % Curves
        set(h_error,'LineStyle','-','Color','k','LineWidth',opt.format.linewidth); % Colors
        xlabel(['Sphere radius ' axisunit]); % Axis label
        ylabel('Error');
        grid(sub_axes,opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
        hold(sub_axes,'off'); % - Figure has been done

        sub_axes=subplot(1,2,2,'Parent',Fig); % Create axes
        hold(sub_axes,'on'); % Active subplot
        h_title=title (' ','FontName',opt.format.fontname);
        h_title.String= {'Distance to the surface cumulative fct','calculated for the phase and a unique sphere with a fitted diameter'};
        h_phase = plot(numericalpsd(current_phase_todo).psd.cumulative_fct(:,1), numericalpsd(current_phase_todo).psd.cumulative_fct(:,2));
        h_sphere=plot(analyticalpsd(current_phase_todo).psd(:,1),analyticalpsd(current_phase_todo).psd(:,2));
        set(h_phase,'LineStyle','-','Marker','none','Color',infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth); % Colors
        set(h_sphere,'LineStyle','-','Marker','none','Color','k','LineWidth',opt.format.linewidth);
        xlabel(['Distance to the surface ' axisunit]); % - Axis label
        ylabel('Cumulative function');
        str1_=['Numerical cumulative funtion of ' char(phasename_todo(current_phase_todo,1))]; % Legend
        str2_=['Analytical cumulative function of a sphere of diameter ' num2str(d50_results(current_phase_todo,1),'%4.2f') ' ' Dunit];
        legend(sub_axes,str1_,str2_,'Location','best');
        grid(sub_axes,opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % - Fontname and fontsize
        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
        hold(sub_axes,'off'); % - Figure has been done

        sgtitle(Fig,['Euclidean distance map cumulative function diameter fitting, ' char(phasename_todo(current_phase_todo,1))] ,'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
        if opt.save.savefig % Save figure
            filename= sprintf('Distance_to_surface_fitted_%s',char(phasename_todo(current_phase_todo,1)));
            function_savefig(Fig, Current_folder, filename, opt.save); % Call function
        end
        if opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end
end

%% CUMULATIVE AND DISTRIBUTION FUNCTIONS PLOT
parameters_distributionfigure.figureposition = [100 100 1500 800];
parameters_distributionfigure.fontname = opt.format.fontname;
parameters_distributionfigure.grid = opt.format.grid;
parameters_distributionfigure.minorgrid = opt.format.minorgrid;
parameters_distributionfigure.fullpath = Current_folder;
parameters_distributionfigure.save = opt.save.savefig;
parameters_distributionfigure.subaxe1_title = 'Cumulative function';
parameters_distributionfigure.subaxe2_title = 'Distribution function';
parameters_distributionfigure.xlabel = ['Distance to surface ' axisunit];
parameters_distributionfigure.axefontsize = opt.format.axefontsize;
parameters_distributionfigure.legendfontsize = opt.format.legendfontsize;
parameters_distributionfigure.titlefontsize = opt.format.titlefontsize;
parameters_distributionfigure.sgtitlefontsize = opt.format.sgtitlefontsize;
parameters_distributionfigure.unit = Dunit;
parameters_distributionfigure.closefig = opt.format.autoclosefig;
for current_phase_todo=1:1:number_phase_todo % Loop over all phases
    parameters_distributionfigure.figurename =  ['Distance to boundary, ' char(phasename_todo(current_phase_todo,1))];
    parameters_distributionfigure.filename = ['Distance_to_boundary_' char(phasename_todo(current_phase_todo,1))];
    parameters_distributionfigure.title = ['Distance to surface, '  char(phasename_todo(current_phase_todo,1))];
    function_probability_distribution_figure(numericalpsd(current_phase_todo).psd,parameters_distributionfigure);
end
    
%% ALONG DIRECTIONS
% Calculate the function Distance(x) defined as 1/L*int(Distance(x)*dx,0,L)=d_50
% For each direction and each phase
alongdirection_parameters.number_phase = number_phase;
alongdirection_parameters.data = results_visualization;
alongdirection_parameters.field = 'distance_to_boundary';
alongdirection_parameters.number_dimension = number_dimension;
alongdirection_parameters.Domain_size = Domain_size;
alongdirection_parameters.voxel_size = voxel_size;
alongdirection_parameters.ignore_value = 0;
alongdirection_parameters.ignore_min = true;
alongdirection_parameters.Variable_name_table = {['Position ' voxel_unit] ['mean distance ' voxel_unit] ['max distance ' voxel_unit] ['std ' voxel_unit]};
alongdirection_parameters.infovol = infovol;
alongdirection_parameters.opt = opt;
alongdirection_parameters.phasename_todo = phasename_todo;
alongdirection_parameters.todo = p.todo;
alongdirection_parameters.Table_filename = 'Distance_to_surface';
alongdirection_parameters.figure_name = 'Distance_to_surface';
alongdirection_parameters.axe_title = 'Distance to surface';
alongdirection_parameters.figure_title = 'Distance to surface along directions';
alongdirection_parameters.figure_filename = 'Distance_to_surface_along_directions';
alongdirection_parameters.ylabel = ['Distance to surface ' axisunit];
alongdirection_parameters.ylabel_unit = Dunit;
alongdirection_parameters.axisunit = axisunit;
alongdirection_parameters.legendname = 'distance';
alongdirection_parameters.mean_val = PSD_results(:,2);
alongdirection_parameters.Current_folder = Current_folder;
[Table_evolution] = Function_along_direction(alongdirection_parameters); % Call function
Results_edmf.Table_evolution = Table_evolution; % Save in main table result

%% PARTICLE DISTANCE TO SURFACE MAP
% Color map
%myColorMap = turbo(256);
str = ['myColorMap = ' p.colormap '(256);'];
eval(str);
if strcmp(p.colormap_background,'white')
    myColorMap(1,:) = 1; 
elseif strcmp(p.colormap_background,'black')
    myColorMap(1,:) = 0; 
elseif strcmp(p.colormap_background,'grey')
    myColorMap(1,:) = [0.75 0.75 0.75];      
end

for current_phase_todo=1:1:number_phase_todo % Loop over phases
    data_distance = results_visualization(current_phase_todo).distance_to_boundary;
    Fig = figure; % Create figure
    Fig.Name= ['Distance to surface, phase ' char(phasename_todo(current_phase_todo))]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
    for current_direction=1:1:number_dimension % Iterate over axe
        sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        h_title=title ({'Slice in the middle',['View normal to ' char(infovol.directionname(current_direction))]}); % Set title font
        % - Plot graphs
        if current_direction==1
            tmp=squeeze(data_distance(round(Domain_size(1)/2),:,:));
            h=image(tmp,'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',char(infovol.directionname(3)));
            t_1y = sprintf('Position along %s ',char(infovol.directionname(2)));
            set(h, 'XData', [0, Domain_size(3)*voxel_size]);
            set(h, 'YData', [0, Domain_size(2)*voxel_size]);
        elseif current_direction==2
            tmp=squeeze(data_distance(:,round(Domain_size(2)/2),:));
            h=image(tmp,'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',char(infovol.directionname(3)));
            t_1y = sprintf('Position along %s ',char(infovol.directionname(1)));
            set(h, 'XData', [0, Domain_size(3)*voxel_size]);
            set(h, 'YData', [0, Domain_size(1)*voxel_size]);
        elseif current_direction==3
            h=image(data_distance(:,:,round(Domain_size(3)/2)),'CDataMapping','scaled');
            t_1x = sprintf('Position along %s ',char(infovol.directionname(1)));
            t_1y = sprintf('Position along %s ',char(infovol.directionname(2)));
            set(h, 'XData', [0, Domain_size(1)*voxel_size]);
            set(h, 'YData', [0, Domain_size(2)*voxel_size]);
        end
        axis equal; axis tight;
%         x_value = get(sub_axes,'XTick');
%         set(sub_axes,'XtickLabel',x_value*voxel_size/1000);
%         y_value = get(sub_axes,'YTick');
%         set(sub_axes,'YtickLabel',y_value*voxel_size/1000);
        % - Axis label
        t_ = xlabel(' ');
        t_2 = axisunit;
        t_.String= [t_1x t_2]; % Sprintf does not accept greek characters
        t_ = ylabel(' ');
        t_2 = axisunit;
        t_.String= [t_1y t_2]; % Sprintf does not accept greek characters
        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
        colormap(myColorMap);
        % Create colorbar
        %colorbar('peer',sub_axes);
        h=colorbar(sub_axes);
        ylabel(h, ['Distance to surface ' axisunit]);
        set(h,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);        
        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,['Distance to surface, ' char(phasename_todo(current_phase_todo))],'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
    if opt.save.savefig % Save figure
        filename= sprintf('View_dmap_%s', char(phasename_todo(current_phase_todo)));
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end
end

%% SECTION ALONG POSITION
number_direction_todo = sum(p.direction_todo);
if sum(p.direction_todo)
    fprintf('Diameters fitted slice per slice, calculation is on-going (slow)...\n\n');

    % Calculation
    for current_phase_todo=1:1:number_phase_todo % Loop over all phases
        binary_phase=zeros(Domain_size);
        binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase
        [fitted2D(current_phase_todo).dir] = Function_particle_size_EDMF_section_Algorithm(binary_phase,p.direction_todo,removeinclusion);
    end

    % Plot
    Fig = figure; % Create figure
    Fig.Name= 'Fitted section disc diameter along directions'; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]); % Full screen figure
    current_direction_todo = 0;
    for current_direction=1:1:3
        if p.direction_todo(current_direction)
            current_direction_todo = current_direction_todo+1;
            sub_axes=subplot(1,number_direction_todo,current_direction_todo,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            h_title=title(['Section disc diameters along ',char(infovol.directionname(current_direction))]);

            current_phase_todo = 0;
            for current_phase=1:1:number_phase % Loop over phases
                if p.todo(current_phase)
                    current_phase_todo=current_phase_todo+1;
                    tmp = fitted2D(current_phase_todo).dir(current_direction);
                    streq = sprintf('D_{eq,%s}',char(infovol.directionname(current_direction)));
                    heq = plot(tmp.x*voxel_size, tmp.diameter_eqdisc(:,1)*voxel_size,'LineWidth',opt.format.linewidth,'DisplayName',[char(infovol.phasename(current_phase,1)) ', ' streq],'Color',infovol.phasecolor(current_phase,:));
                end
            end

            % Axis label
            t_ = xlabel(' ');
            t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
            t_2 = axisunit;
            t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
            ylabel(['Fitted diameters ' axisunit]);
            h_legend = legend(sub_axes,'Location','best');
            h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
            % - Grid
            grid(sub_axes,opt.format.grid); % Display grid
            set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
            h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
    end
    sgtitle(Fig,'Section disc diameters along directions','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
    if opt.save.savefig % Save figure
        filename= 'Diameters_EDMF_along_direction';
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end

    % Table
    clear Variable_name_table;
    Variable_name_table={['Position ' voxel_unit]};
    for current_phase_todo=1:1:number_phase_todo
        Variable_name_table(current_phase_todo+1)=phasename_todo(current_phase_todo);
    end
    % Create table
    for current_direction=1:1:number_dimension % Loop over all directions
        if p.direction_todo(current_direction)
            tmp = fitted2D(1).dir(current_direction);
            tmp_array = [(tmp.x)' (tmp.diameter_eqdisc)];
            for current_phase_todo=2:1:number_phase_todo
                tmp = fitted2D(current_phase_todo).dir(current_direction);
                tmp_array = [tmp_array (tmp.diameter_eqdisc)];            
            end
            EDMF_evolution.direction(current_direction).table = array2table(tmp_array*voxel_size,'VariableNames',Variable_name_table);
        end
    end
    Results_edmf.Alongdirections = EDMF_evolution; % Save in main table result

    if opt.save.xls
        filename = 'Diameters_EDMF_along_directions'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        current_direction_todo = 0;
        for current_direction=1:1:number_dimension % Loop over all directions
            if p.direction_todo(current_direction)
                current_direction_todo = current_direction_todo +1;
                DATA_writetable.sheet(current_direction_todo).name = char(infovol.directionname(current_direction));
                DATA_writetable.sheet(current_direction_todo).table = EDMF_evolution.direction(current_direction).table;
            end
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
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
    property_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1,2);
    property_voxelsizedependence(1,1,:)=voxel_size;
    property_voxelsizedependence(1,2:end,1)=d50_results(:,1)';
    property_voxelsizedependence(1,2:end,2)=dmap_results(:,2)';    

    current_phase_todo = 0;
    for current_phase=1:1:number_phase % Loop over all phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            PSDresized(current_phase_todo).iteration(1).voxelsize = voxel_size;
            PSDresized(current_phase_todo).iteration(1).psd = numericalpsd(current_phase_todo).psd;
        end
    end

    % Loop on each voxel size
    for current_iteration=1:1:number_resize

        % % Microstructure scaling
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1,:)=current_voxel_size;
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
        Current_domain_size = size(Phase_microstructure_resized);

        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over all phases
            if p.todo(current_phase)
                time_cpu_start_phase = cputime; % CPU start
                time_clock_start_phase = tic; % Stopwatch start

                current_phase_todo=current_phase_todo+1;

                % % Algorithm: SPECIFIC FOR EACH FILE
                % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
                binary_phase=zeros(Current_domain_size); % Initialization
                code_tmp = infovol.phaselabel(current_phase);
                Numbervoxel_phase_tmp= sum(sum(sum(Phase_microstructure_resized==code_tmp )));
                binary_phase(Phase_microstructure_resized == code_tmp) = 1; % Binary phase
                [distance_transform_resized, fitted_diameter_resized, numericalpsd_resized(current_phase_todo).psd, ~, ~] = Function_particle_size_distancemap_Algorithm(Phase_microstructure_resized, code_tmp, removeinclusion, distancelabel, current_voxel_size, density_fct_parameters);
                all_distance_resized = distance_transform_resized(binary_phase==1);
                property_voxelsizedependence(current_iteration+1,current_phase_todo+1,1) = fitted_diameter_resized;
                property_voxelsizedependence(current_iteration+1,current_phase_todo+1,2) = mean(all_distance_resized);
                PSDresized(current_phase_todo).iteration(current_iteration+1).voxelsize = current_voxel_size;
                PSDresized(current_phase_todo).iteration(current_iteration+1).psd = numericalpsd_resized(current_phase_todo).psd;

                % % Time
                timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

            end
        end
        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];

    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    property_voxelsizedependence(:,:,1) = sortrows(property_voxelsizedependence(:,:,1),1);  
    property_voxelsizedependence(:,:,2) = sortrows(property_voxelsizedependence(:,:,2),1);       
    
    %% EXTRAPOLATION TO 0 nm
    if strcmp(p.scaling_extrapolation,'Linear')
        interpolation_voxelsize_order=1;
    elseif strcmp(p.scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order=2;
    elseif strcmp(p.scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order = min(interpolation_voxelsize_order,max_order);

    tmp = zeros(number_resize+2,number_phase_todo+1,2); % + 0 nm and + initial voxel size
    x=property_voxelsizedependence(:,1,1);
    str_correlation(1).name = 'Particle_diameter_mean_EDMF';
    str_correlation(2).name = 'Distance_to_surface_mean';    
    fprintf('> Mean diameter dependence with the voxel size\n');
    fprintf('  Extrapolation to zero voxel size: polynomial of order %i\n\n',interpolation_voxelsize_order);
    for kparameter=1:1:2
        for current_phase_todo=1:1:number_phase_todo
            y=property_voxelsizedependence(:,current_phase_todo+1,kparameter);
            pi = polyfit(x,y,interpolation_voxelsize_order);
            vq = polyval(pi,0);
            tmp(1,current_phase_todo+1,kparameter)=vq;
            interpolation_voxelsize(current_phase_todo,kparameter).pi=pi;
            % For correlation
            results_correlation(current_phase_todo).([str_correlation(kparameter).name '_extrapolated']) = vq;            
        end
        tmp(2:end,:,kparameter) = property_voxelsizedependence(:,:,kparameter);
    end
    property_voxelsizedependence = tmp; clear tmp;

    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_phase_todo=1:1:number_phase_todo
        Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
    end
    % Table
    Table_d50_Edmf_voxelsizedependence = array2table(property_voxelsizedependence(:,:,1),...
        'VariableNames',Variable_name_table);
    Table_distance2surface_voxelsizedependence = array2table(property_voxelsizedependence(:,:,2),...
        'VariableNames',Variable_name_table);    
    
    %% DISPLAY TEXT RESULTS
    disp(Table_d50_Edmf_voxelsizedependence)
    fprintf('> Distance to surface dependence with the voxel size:\n\n');
    fprintf('  Extrapolation to zero voxel size: polynomial of order %i\n\n',interpolation_voxelsize_order);
    disp(Table_distance2surface_voxelsizedependence)

    
    %% SAVE RESULTS
    Results_edmf.voxelsizedependence_d50 = Table_d50_Edmf_voxelsizedependence; % Save in main table result
    Results_edmf.voxelsizedependence_dsurface = Table_distance2surface_voxelsizedependence; 
    if opt.save.xls
        filename = 'EMDF_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='D50_EDMF';
        DATA_writetable.sheet(1).table=Table_d50_Edmf_voxelsizedependence;
        DATA_writetable.sheet(2).name='Distance_to_surface';
        DATA_writetable.sheet(2).table=Table_distance2surface_voxelsizedependence;        
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    %% FIGURES
    parameters_figure.plotlog = false;
    parameters_figure.propertyname = 'Mean diameter';
    parameters_figure.method = 'EDMF';
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,1);
    parameters_figure.number_phase = number_phase_todo;
    parameters_figure.str_ylabel = ['D_{50} ' axisunit];
    parameters_figure.propertynameunit = Dunit;
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,1);
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'D50_Edmf_voxel_size_dependence';
    parameters_figure.infovol = infovol;
    parameters_figure.opt = opt;
    parameters_figure.todo = p.todo;
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    

    parameters_figure.propertyname = 'Mean distance to surface';
    parameters_figure.method = [];
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,2);
    parameters_figure.str_ylabel = ['Distance to surface ' axisunit];
    parameters_figure.propertynameunit = Dunit;
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,2);
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'Distance_to_surface_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures     

    %% FIGURES 2 (DISTRIBUTION)
    for current_phase_todo=1:1:number_phase_todo % Loop over all phases
        parameters_distributionfigure.figurename =  ['Distance to surface, ' char(phasename_todo(current_phase_todo,1))];
        parameters_distributionfigure.filename = ['Distance_to_surface_voxelsize_' char(phasename_todo(current_phase_todo,1))];
        parameters_distributionfigure.subaxe1_title = 'Cumulative function';
        parameters_distributionfigure.subaxe2_title = 'Distribution function';
        parameters_distributionfigure.title = ['Distance to surface voxel size dependence (C-PSD), ' char(phasename_todo(current_phase_todo,1))];
        parameters_distributionfigure.xlabel = ['Distance to surface ' axisunit];
        function_probability_distribution_size_figure(PSDresized(current_phase_todo),parameters_distributionfigure);
    end
end

%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

if p.RVE.number_RVE>0
    str_property(1).corrname = 'D50edmf';
    str_property(2).corrname = 'D2S';
    str_property(1).propertyname = 'Mean diameter';
    str_property(2).propertyname = 'Distance to surface';
    res(1).Wholevolume_results = d50_results(:,1);
    res(2).Wholevolume_results = dmap_results(:,2);

    for k_RVE = 1:1:p.RVE.number_RVE % Loop over all RVE asked
        RVEparameters = p.RVE.RVE(k_RVE);
        for property_RVE=1:1:2
            Rk(property_RVE).RVE(k_RVE).RVEparameters = RVEparameters;
        end
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
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_phase_todo,3,3, 2); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
            % 1 fitted diameter, 2 mean distance to surface
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
                Result_nestedRVE(k_nestedRVE+1,1,:,:,1,:) = Wholevolume_size(1);
                if k_nestedRVE ==1
                    fprintf('       Nested analysis\n');
                end
                if k_nestedRVE > 0
                    fprintf('          Cropping iteration: #%i/%i (cubic root volume=%1.1f %s)\n',k_nestedRVE,n_nestedRVE,Wholevolume_size(1),infovol.unit);
                end                
                if strcmp(RVEparameters.type,'C')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2,:) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'D')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3,:) = Wholevolume_size(2);                    
                elseif strcmp(RVEparameters.type,'G')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3,:) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'H')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2,:) = Wholevolume_size(2);                    
                end
            end

            % Information about subdomains
            for property_RVE=1:1:2
                Rk(property_RVE).RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
                    'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Equivalent cubic length' 'Equivalent section length' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
            end
            [number_subdomain,~] = size(All_subdomain); % The number of subdomain
            number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size

            %% ALGORITHM
            % Colunm 1 is the subdomain id
            % Colunm 2 and 3 are the sizes of the subdomain.
            Property_eachsubdomain = zeros(number_subdomain,number_phase_todo+4,2); % 1 fitted diameter, 2 mean distance to surface
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
                Current_domain_size = size(current_subdomain);                

                % CPU and stopwatch time - start
                time_cpu_start_volume = cputime; % CPU start
                time_clock_start_volume = tic; % Stopwatch start
                % Number of voxel of the current resized microstructure
                voxel_number_tmp=numel(current_subdomain);

                Property_eachsubdomain(subdomain_id,1,:)=subdomain_id;
                % Equivalent size of the subdomain
                Property_eachsubdomain(subdomain_id,2,:)=All_subdomain(subdomain_id,4)*voxel_size; % Cubic root length
                Property_eachsubdomain(subdomain_id,3,:)=All_subdomain(subdomain_id,5)*voxel_size; % Square root length
                Property_eachsubdomain(subdomain_id,4,:)=All_subdomain(subdomain_id,6)*voxel_size; % Length

                current_phase_todo = 0;
                for current_phase=1:1:number_phase % Loop over all phases
                    if p.todo(current_phase)
                        time_cpu_start_phase = cputime; % CPU start
                        time_clock_start_phase = tic; % Stopwatch start
                        current_phase_todo=current_phase_todo+1;

                        % % Algorithm: SPECIFIC FOR EACH FILE
                        code_tmp = infovol.phaselabel(current_phase);
                        Numbervoxel_phase_tmp= sum(sum(sum(current_subdomain==code_tmp )));
                        binary_phase=zeros(Current_domain_size); % Initialization
                        binary_phase(current_subdomain == code_tmp) = 1;
                        [distance_transform_subdomain, fitted_diameter_subdomain, ~, ~, ~] = Function_particle_size_distancemap_Algorithm(current_subdomain, code_tmp, removeinclusion, distancelabel, voxel_size, density_fct_parameters);
                        all_distance_subdomain = distance_transform_subdomain(binary_phase==1);
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4,1)=fitted_diameter_subdomain;
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4,2)=mean(all_distance_subdomain);

                        % % Time
                        timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

                    end
                end
                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
            end

            for property_RVE=1:1:2
                %% STATISTICAL ANALYSIS and RVE SIZE
                [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_phase_todo,GROUP_SUBDOMAIN,Property_eachsubdomain(:,:,property_RVE),voxel_size,RVEparameters);

                %% SAVE FOR CORRELATION
                if k_nestedRVE == 0
                    for k_threshold=1:1:n_threshold
                        current_phase_todo = 0;
                        for current_phase=1:1:number_phase
                            if p.todo(current_phase)
                                current_phase_todo=current_phase_todo+1;
                                if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                    if Size_RVE(1,current_phase_todo,2,1)~=0
                                        str_ = [str_property(property_RVE).corrname '_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,1);
                                    end
                                    if strcmp(RVEparameters.type,'C')
                                        if Size_RVE(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                            str_(str_=='.')='p';
                                            results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                        end
                                    elseif strcmp(RVEparameters.type,'D')
                                        if Size_RVE(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_RVE_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                            str_(str_=='.')='p';
                                            results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                        end
                                    end
                                else
                                    if Size_convergence(1,current_phase_todo,2,1)~=0
                                        str_ = [str_property(property_RVE).corrname '_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,1);
                                    end
                                    if strcmp(RVEparameters.type,'G')
                                        if Size_convergence(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                            str_(str_=='.')='p';
                                            results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
                                        end
                                    elseif strcmp(RVEparameters.type,'H')
                                        if Size_convergence(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_conv_areaFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
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
                Rk(property_RVE).RVE = Function_subdomains_manage_results(Property_eachsubdomain(:,:,property_RVE), Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,Rk(property_RVE).RVE,k_RVE,number_phase,number_phase_todo,infovol,p);

                if RVEparameters.donested
                    for k_threshold=1:1:n_threshold
                        current_phase_todo = 0;
                        for current_phase=1:1:number_phase
                            if p.todo(current_phase)
                                current_phase_todo=current_phase_todo+1;
                                if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1,property_RVE) = Size_RVE(k_threshold,current_phase_todo,1,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1,property_RVE) = Size_RVE(k_threshold,current_phase_todo,2,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1,property_RVE) = Size_RVE(k_threshold,current_phase_todo,3,1);
                                else
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1,property_RVE) = Size_convergence(k_threshold,current_phase_todo,1,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1,property_RVE) = Size_convergence(k_threshold,current_phase_todo,2,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1,property_RVE) = Size_convergence(k_threshold,current_phase_todo,3,1);
                                end
                                if strcmp(RVEparameters.type,'C')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2,property_RVE) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2,property_RVE) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2,property_RVE) = Size_RVE(k_threshold,current_phase_todo,3,2);
                                elseif strcmp(RVEparameters.type,'D')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3,property_RVE) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3,property_RVE) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3,property_RVE) = Size_RVE(k_threshold,current_phase_todo,3,2);
                                elseif strcmp(RVEparameters.type,'G')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3,property_RVE) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3,property_RVE) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3,property_RVE) = Size_convergence(k_threshold,current_phase_todo,3,2);
                                elseif strcmp(RVEparameters.type,'H')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2,property_RVE) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2,property_RVE) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2,property_RVE) = Size_convergence(k_threshold,current_phase_todo,3,2);
                                end
                            end
                        end
                    end
                end

                %% TEXT DISPLAY AND SAVE RESULTS
                if k_nestedRVE == 0
                    if property_RVE==1
                        RVEparameters.disp_parRVE = true;
                    else
                        RVEparameters.disp_parRVE = false;
                    end
                    Function_subdomains_display_and_save(Rk(property_RVE).RVE,k_RVE,RVEparameters,number_phase,number_phase_todo,str_property(property_RVE).propertyname,Sub_folder_RVE,opt,infovol,p);
                end

                %% FIGURES
                if k_nestedRVE == 0
                    parameters_figure.propertyname = str_property(property_RVE).propertyname;
                    parameters_figure.propertynameunit = Dunit;
                    parameters_figure.RVE = RVEparameters;
                    parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                    parameters_figure.savefolder = Sub_folder_RVE;
                    parameters_figure.number_phase = number_phase;
                    parameters_figure.number_phase_todo = number_phase_todo;
                    parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                    parameters_figure.Property_eachsubdomain = Property_eachsubdomain(:,:,property_RVE);
                    parameters_figure.derivative_convergence = derivative_convergence;
                    parameters_figure.relativedifference_convergence = relativedifference_convergence;
                    parameters_figure.Size_RVE = Size_RVE;
                    parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
                    parameters_figure.Size_convergence = Size_convergence;
                    parameters_figure.Wholevolume_size = Wholevolume_size;
                    parameters_figure.Wholevolume_results = res(property_RVE).Wholevolume_results;
                    parameters_figure.infovol = infovol;
                    parameters_figure.todo = p.todo;
                    parameters_figure.opt = opt;
                    Function_create_figures_RVE(parameters_figure) % Figures
                end
            end
        end

        %% NESTED ANALYSIS RESULT
        if RVEparameters.donested
            for property_RVE=1:1:2
                % Table
                [Rk(property_RVE).RVE] = Function_nestedtable(Rk(property_RVE).RVE,k_RVE,RVEparameters,number_phase,str_property(property_RVE).propertyname,Result_nestedRVE(:,:,:,:,:,property_RVE),Sub_folder_RVE,opt,infovol,p);
                % Figure
                parameters_figure.propertyname = str_property(property_RVE).propertyname;
                parameters_figure.Result_nestedRVE = Result_nestedRVE(:,:,:,:,:,property_RVE);
                Function_create_figures_nestedRVE(parameters_figure) % Figures
                % Save
                Rk(property_RVE).RVE(k_RVE).nestedanalysis = Result_nestedRVE(:,:,:,:,:,property_RVE);
            end
        end

    end
    Results_edmf.RVE.meandiameter = Rk(1).RVE; % Save in main table result
    Results_edmf.RVE.distance2surface = Rk(2).RVE;
end

%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_edmf.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_edmf.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_edmf.Table_date = Table_date;

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
    Function_Writetable(Current_folder,'D50_Edmf_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'D50_Edmf_calculation_time', 'Diameter (EMDF)', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_D50_edmf'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_D50_edmf'],'Results_edmf')
end
%% SAVE VISUALIZATION
Current_folder = [infovol.volpath 'Visualization' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Visualization_particlediameter_edmf'],'results_visualization','-v7.3');    

end