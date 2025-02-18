function [] = Charact_Volumefractions(Labels,Nanoporosity,infovol,opt,p)
% Calculate volume fractions
% Charact_Volumefractions(Labels, Nanoporosity, infovol, opt, p) - when used with the MATBOX toolbox GUI
% or
% Charact_Volumefractions(Labels, voxelsize, unit) - when used as a standalone function
% with: Labels, a 3D array: the 3D segmented volumes
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       e.g.: Charact_Volumefractions(<your_3d_array>, 0.4, 'um');

%% DEFAULT VALUES
expected_number_argument = 5;
nargin; % Number of input variable when the function is call

show_folder = false;
if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 3 % Case for function called as: Function_Volume_fractions(Labels, voxelsize, unit). Standalone use.
        voxelsize = Nanoporosity; clear Nanoporosity;
        unit = infovol; clear infovol;
        Nanoporosity = [];
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        %infovol.volumesubfolder = ['Characterization_' char(t)];
        if ispc
            infovol.mainfolder=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
        else
            infovol.mainfolder = pwd;
        end
        %infovol.volpath = fullfile(infovol.mainfolder, infovol.volumesubfolder); 
        infovol.volpath = fullfile(infovol.mainfolder, 'Characterization', char(t)); 
        
        infovol.sub = [];
        show_folder = true;

        % Set default phase information
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        infovol.phaselabel = unique(Labels);
        for k=1:1:length(infovol.phaselabel) % Loop over all unique values
            infovol.phasename(k,1) = {['Phase ' num2str(infovol.phaselabel(k))]}; % Assign phase name based on code value
            infovol.phasecolor(k,:) = allcolor(k,:);
            p.todo(k)=1;
        end
        infovol.voxelsize = voxelsize;
        infovol.unit = unit;

        % No background
        infovol.isbackground = 0;
        
        % Set default direction information
        sz = size(Labels);
        if length(sz)==2
            infovol.directionname = {'in-plane direction 1';'in-plane direction 2'};
        else
            infovol.directionname = {'in-plane direction 1';'in-plane direction 2';'through-plane direction'};
        end
        
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
        opt.format.grid = 'on'; opt.format.minorgrid = 'off';
        opt.format.autoclosefig = false;
        
        % No Voxel size dependence analysis
        p.scaling = 1;
        % No Representative Volume Element analysis
        p.RVE.number_RVE = 0;
        % Fractal analysis
        p.fractal_bertei.todo = false;     

    else % Incorrect number of argument
        disp 'Error calling Charact_Volumefractions. Wrong number of argument.'
        help Charact_Volumefractions
        return
    end
    
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
Current_folder = fullfile(infovol.volpath,'Volume_fractions',infovol.sub);
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
if show_folder
    fprintf('Results are saved in %s\n\n',Current_folder)
end

%% VOLUME INFORMATION
Domain_size = size(Labels);
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

L = length(p.todo);
if L < number_phase
    p.todo = [p.todo; ones(number_phase-L)];
elseif L > number_phase
    p.todo = p.todo(1:number_phase);
end

if infovol.isbackground
    p.todo(1) = 0;
end

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

if infovol.isbackground
    nvoxel_background = sum(sum(sum(Labels==0)));
    vf_background = nvoxel_background/voxel_number;
end

current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);

        % % Algorithm
        phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
        Numbervoxel_phase(current_phase_todo,1)= sum(sum(sum(Labels==phaselabel(current_phase_todo,1) )));
        Volumefraction_phase(current_phase_todo,1) = Numbervoxel_phase(current_phase_todo,1)/voxel_number; % Defintion of volume fraction
        if infovol.isbackground
            Volumefraction_phase(current_phase_todo,1) = Volumefraction_phase(current_phase_todo,1) / (1-vf_background);
        end

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

if infovol.isbackground
    Table_Background = table({'Background'},0,vf_background,nvoxel_background,...
        'VariableNames',{'Phase' 'Label' 'Volume fraction' 'Number of voxel'});%
end

%% SAVE TABLES
if opt.save.xls
    filename = 'Volume_fractions'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Volume_fraction';
    DATA_writetable.sheet(1).table=Table_Volumefraction;
    if infovol.isbackground
        DATA_writetable.sheet(2).name='Background';
        DATA_writetable.sheet(2).table=Table_Background;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp(Table_Volumefraction)
if infovol.isbackground
    fprintf(['Background volume fraction: ' num2str(vf_background) '\n\n'])
end
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

if infovol.isbackground
    for current_direction=1:1:number_dimension % Loop over all directions
        ResultsBackground.direction(current_direction).volumefraction = zeros(Domain_size(current_direction),1); % Initialization
    end

    % Direction 1
    for current_position=1:1:Domain_size(1) % Loop over postion
        Numbervoxel_position= sum(sum(Labels(current_position,:,:)==0));
        ResultsBackground.direction(1).volumefraction(current_position,1) = Numbervoxel_position / Results.direction(1).numbervoxelslice;
    end
    % Direction 2
    for current_position=1:1:Domain_size(2) % Loop over postion
        Numbervoxel_position= sum(sum(Labels(:,current_position,:)==0));
        ResultsBackground.direction(2).volumefraction(current_position,1) = Numbervoxel_position / Results.direction(2).numbervoxelslice;
    end
    if number_dimension==3 % 3D case
        % Direction 3
        for current_position=1:1:Domain_size(3) % Loop over postion
            Numbervoxel_position= sum(sum(Labels(:,:,current_position)==0));
            ResultsBackground.direction(3).volumefraction(current_position,1) = Numbervoxel_position / Results.direction(3).numbervoxelslice;
        end
    end

    current_phase_todo = 0;
    for current_phase = 1:1:number_phase % Loop over all phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;

            label = infovol.phaselabel(current_phase);
            % Direction 1
            for current_position=1:1:Domain_size(1) % Loop over postion
                Numbervoxel_position= sum(sum(Labels(current_position,:,:)==label));
                Results.direction(1).volumefraction(current_position,current_phase_todo+1) = (Numbervoxel_position / Results.direction(1).numbervoxelslice) / (1-ResultsBackground.direction(1).volumefraction(current_position,1));
            end
            % Direction 2
            for current_position=1:1:Domain_size(2) % Loop over postion
                Numbervoxel_position= sum(sum(Labels(:,current_position,:)==label));
                Results.direction(2).volumefraction(current_position,current_phase_todo+1) = (Numbervoxel_position / Results.direction(2).numbervoxelslice) / (1-ResultsBackground.direction(2).volumefraction(current_position,1));
            end
            if number_dimension==3 % 3D case
                % Direction 3
                for current_position=1:1:Domain_size(3) % Loop over postion
                    Numbervoxel_position= sum(sum(Labels(:,:,current_position)==label));
                    Results.direction(3).volumefraction(current_position,current_phase_todo+1) = (Numbervoxel_position / Results.direction(3).numbervoxelslice) / (1-ResultsBackground.direction(3).volumefraction(current_position,1));
                end
            end
        end
    end

else
    current_phase_todo = 0;
    for current_phase = 1:1:number_phase % Loop over all phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;

            label = infovol.phaselabel(current_phase);
            % Direction 1
            for current_position=1:1:Domain_size(1) % Loop over postion
                Numbervoxel_position= sum(sum(Labels(current_position,:,:)==label));
                Results.direction(1).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(1).numbervoxelslice;
            end
            % Direction 2
            for current_position=1:1:Domain_size(2) % Loop over postion
                Numbervoxel_position= sum(sum(Labels(:,current_position,:)==label));
                Results.direction(2).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(2).numbervoxelslice;
            end
            if number_dimension==3 % 3D case
                % Direction 3
                for current_position=1:1:Domain_size(3) % Loop over postion
                    Numbervoxel_position= sum(sum(Labels(:,:,current_position)==label));
                    Results.direction(3).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(3).numbervoxelslice;
                end
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

    if infovol.isbackground
        y_ = ResultsBackground.direction(current_direction).volumefraction;
        plot(x_,y_,'Color', infovol.phasecolor(1,:),'LineWidth',opt.format.linewidth,'LineStyle','--','DisplayName',char(infovol.phasename(current_phase,1)));
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
    if infovol.isbackground
        str_legend(current_phase_todo+1).name = ['Background, <\epsilon>=' num2str(vf_background,'%1.3f')];
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


%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%

p.scaling(p.scaling==1)=[];
if ~isempty(p.scaling)
    % Metric parameters
    pMET.fct_name = 'Volume fraction';
    pMET.fractal_bertei = p.fractal_bertei.todo;
    pMET.metric(1).name = 'Volume fraction';
    pMET.metric(1).shortname_correlation = 'vf'; 
    pMET.metric(1).unit = [];
    pMET.metric(1).result_initialres = Volumefraction_phase;
    [Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume] = ImageResolution_erroranalysis(pMET, Labels, [], voxel_size, number_phase, number_phase_todo, p.todo, p.scaling, p.scaling_extrapolation, Current_folder, opt, infovol, Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume);
end

%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

if p.RVE.number_RVE>0
    % Metric parameters
    pMET.fct_name = 'Volume fraction';
    pMET.p = p;
    pMET.metric(1).name = 'Volume fraction';    
    pMET.metric(1).shortname_correlation = 'vf'; 
    pMET.metric(1).unit = [];    
    pMET.metric(1).result_wholevolume = Volumefraction_phase;
    [Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume] = RVE_main(p.RVE.RVE, pMET, Labels, [], voxel_size, number_phase, number_phase_todo, Current_folder, opt, infovol, Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume); % Call main function
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
Current_folder = fullfile(infovol.volpath, 'Correlation');
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

fullpath = fullfile(Current_folder,'Correlation_volume_fraction');
if exist([fullpath '.mat'], 'file') == 2
    tmp = load([fullpath '.mat']);
    old = tmp.results_correlation;
    fields_old = fieldnames(old);
    current_length = length(results_correlation);
    for k=1:1:length(old)
        for k_field = 1:1:length(fields_old)
            results_correlation(current_length+k).(char(fields_old(k_field))) = old(k).(char(fields_old(k_field)));
        end
    end 
end
save(fullpath,'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = fullfile(infovol.volpath, 'Summary');
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save(fullfile(Current_folder,['Results_volume_fraction_' infovol.sub]),'Results_Volumefraction')
end

end