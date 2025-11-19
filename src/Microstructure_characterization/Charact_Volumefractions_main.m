function [] = Charact_Volumefractions_main(Labels,Nanoporosity,Wetting,infovol,opt,p)
% Calculate volume fractions
% Charact_Volumefractions_main(Labels, Nanoporosity, Wetting, infovol, opt, p) - when used with the MATBOX toolbox GUI
% or
% Charact_Volumefractions_main(Labels, voxelsize, unit, name) - when used as a standalone function
% with: Labels, a 3D array: the 3D segmented volumes
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       (optional) name, a string: folder where result is saved  
%       e.g.: Charact_Volumefractions_main(<your_3d_array>, 0.4, 'um');

%% DEFAULT VALUES IF RUN IN COMMAND LINE
expected_number_argument = 6;
nargin; % Number of input variable when the function is call

Run_in_commandline = false;
if nargin ~= expected_number_argument % Unexpected number of argument

    if nargin == 3 || nargin==4 % Case for function called as: Charact_Volumefractions_main(Labels, voxelsize, unit, name). Standalone use.
        Run_in_commandline = true;
        voxelsize = Nanoporosity; Nanoporosity = [];
        unit = Wetting; Wetting = [];
        if nargin==4
            name = infovol; infovol = [];
        end

        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        %infovol.volumesubfolder = ['Characterization_' char(t)];
        if ispc
            infovol.mainfolder=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
        else
            infovol.mainfolder = pwd;
        end
        %infovol.volpath = fullfile(infovol.mainfolder, infovol.volumesubfolder);
        if nargin==3
            infovol.volpath = fullfile(infovol.mainfolder, 'Characterization', char(t));
        else
            infovol.volpath = fullfile(infovol.mainfolder, 'Characterization', name);
        end
        infovol.sub = [];

        % Set default phase information
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        infovol.phaselabel = unique(Labels);
        for k=1:1:length(infovol.phaselabel) % Loop over all unique values
            infovol.phasename(k,1) = {['Phase ' num2str(infovol.phaselabel(k))]}; % Assign phase name based on code value
            infovol.phasecolor(k,:) = allcolor(k,:);
            p.todo(k)=1;
        end
        p.combined_todo = false;
        infovol.voxelsize = voxelsize;
        infovol.unit = unit;
        infovol.filename = 'Your array';

        % No background
        infovol.isbackground = 0;

        % Set default direction information
        sz = size(Labels);
        if length(sz)==2
            infovol.directionname = {'in-plane direction 1';'in-plane direction 2'};
            p.plotdirections = [1 1];
        else
            infovol.directionname = {'in-plane direction 1';'in-plane direction 2';'through-plane direction'};
            p.plotdirections = [1 1 1];
        end

        % Set default options
        opt.save.mat = true;
        opt.save.xls = true;
        opt.save.savefig = true;
        opt.save.fig_infig = true;
        opt.save.fig_format = {'png'};
        opt.save.png_DPI = 300;
        opt.save.Height_foroneplot = 8;
        opt.save.autoclosefig = false;
        
        % Set display options
        opt.format.fontname = 'Times New Roman';
        opt.format.axefontsize =  12;
        opt.format.legendfontsize =  10;
        opt.format.titlefontsize =  14;
        opt.format.sgtitlefontsize =  16;
        opt.format.linewidth = 2;
        opt.format.grid = 'on'; opt.format.minorgrid = 'off';
        opt.includefilenameintitle = 'off';
        opt.format.tile_spacing = 'compact';
        opt.format.layout_padding = 'compact';
        opt.format.includefilenameintitle = false;

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
if Run_in_commandline
    fprintf('Results are saved in %s\n\n',Current_folder)
end

%% VOLUME INFORMATION
sz = size(Labels);
if length(sz)==2 % Add third column with unit value
    sz = [sz 1];
    number_dimension = 2;
else
    number_dimension =3; % 3D case
end
number_phase = length(infovol.phaselabel); % Number of phase

background_label = [];
delta_back = 0;
if infovol.isbackground
    number_phase = number_phase -1;
    delta_back = 1;
    background_label = infovol.phaselabel(1);
end

voxel_number = prod(sz); % Number of voxel
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

for current_phase=1:1:number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        results_correlation(current_phase_todo).name = char(infovol.phasename(current_phase,1));
        % Domain length is saved (only for the volume fraction file)
        results_correlation(current_phase_todo).domain_length = voxel_number^(1/number_dimension) * voxel_size;
        results_correlation(current_phase_todo).voxel_size = voxel_size;
    end
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

if isempty(Nanoporosity)
    p.combined_todo = false;
end

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

% Initialization (generic)
number_phase_todo = sum(p.todo); % How many phase are we going to analyse ?
phaselabel_todo = zeros(number_phase_todo,1);
phasename_todo = cell(number_phase_todo,1);
phasecolor_todo = zeros(number_phase_todo,3);
current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        phaselabel_todo(current_phase_todo,1) = infovol.phaselabel(current_phase+delta_back);
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase+delta_back,1);
        phasecolor_todo(current_phase_todo,:) = infovol.phasecolor(current_phase+delta_back,:);
    end
end
timedata = zeros(number_phase_todo+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_phase_todo+1,1);

% Initialization (algorithm-specific)
% Phase-wise
Vf_phase_label=zeros(number_phase_todo,1);
Vf_phase_solidpart=zeros(number_phase_todo,1);
Vf_phase_porouspart_idealwetting=zeros(number_phase_todo,1);
Vf_phase_porouspart_partiallwetting=zeros(number_phase_todo,1);
% Ratio
ratio_solid=zeros(number_phase_todo,1);
ratio_poreidealwetting=zeros(number_phase_todo,1);
ratio_porepartialwetting=zeros(number_phase_todo,1);
% Statistics
stats_nanoporosity=zeros(number_phase_todo,6);
stats_wetting=zeros(number_phase_todo,6);
% Combined phase
if p.combined_todo
    nvoxel.onlypore = 0;
    nvoxel.onlysolid = 0;
    nvoxel.mixed = 0;
    nvoxel.solid = 0;
    nvoxel.pore_idealwetting = 0;
    nvoxel.pore_partialwetting = 0;
end

% Background
n_voxel_background = 0;
if infovol.isbackground
    n_voxel_background = sum(sum(sum(Labels==background_label)));
    vf_background = n_voxel_background/voxel_number;
end
n_voxel = numel(Labels) - n_voxel_background;

% Start
for current_phase_todo=1:1:number_phase_todo % Loop over all phases
    time_cpu_start_phase = cputime; % CPU start
    time_clock_start_phase = tic; % Stopwatch start

    % % Algorithm
    if Run_in_commandline
        n.n_voxel_label = sum(sum(sum( Labels==phaselabel_todo(current_phase_todo,1) )));
        Vf_phase_label(current_phase_todo,1) = n.n_voxel_label / n_voxel;
    else
        [st,vf,r,n] = Charact_Volumefractions_algorithm(Labels, phaselabel_todo(current_phase_todo,1), n_voxel, Nanoporosity, Wetting);

        Vf_phase_label(current_phase_todo,1) = vf.phase_label;
        Vf_phase_solidpart(current_phase_todo,1) = vf.phase_solid;
        Vf_phase_porouspart_idealwetting(current_phase_todo,1) = vf.phase_pore_idealwetting;
        Vf_phase_porouspart_partiallwetting(current_phase_todo,1) = vf.phase_pore_partialwetting;

        ratio_solid(current_phase_todo,1) = r.solid;
        ratio_poreidealwetting(current_phase_todo,1) = r.poreidealwetting;
        ratio_porepartialwetting(current_phase_todo,1) = r.porepartialwetting;

        stats_nanoporosity(current_phase_todo,:) = st.nanoporosity;
        stats_wetting(current_phase_todo,:) = st.wetting;

        if p.combined_todo
            nvoxel.onlypore = nvoxel.onlypore + n.onlypore;
            nvoxel.onlysolid = nvoxel.onlysolid + n.onlysolid;
            nvoxel.mixed = nvoxel.mixed + n.mixed;
            nvoxel.solid = nvoxel.solid + n.solid;
            nvoxel.pore_idealwetting = nvoxel.pore_idealwetting + n.pore_idealwetting;
            nvoxel.pore_partialwetting = nvoxel.pore_partialwetting + n.pore_partialwetting;
        end

        % % Correlation
        results_correlation(current_phase_todo).volume_fraction = Vf_phase_label(current_phase_todo,1); % Save for correlation
    end

    % Time
    timedata_domain(current_phase_todo+1,1) = phasename_todo(current_phase_todo,1);
    timedata(current_phase_todo+1,1) = n.n_voxel_label;
    timedata(current_phase_todo+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
    timedata(current_phase_todo+1,3) = toc(time_clock_start_phase); % CPU elapsed time
end

if p.combined_todo && ~isempty(Nanoporosity)
    % Categories
    vf_onlypore = nvoxel.onlypore/n_voxel;
    vf_onlysolid = nvoxel.onlysolid/n_voxel;
    vf_mixed = nvoxel.mixed/n_voxel;

    % Porosity, solid, and air
    vf_pore_idealwetting = nvoxel.pore_idealwetting/n_voxel;
    vf_solid = nvoxel.solid/n_voxel;

    current_phase_todo=current_phase_todo+1;
    results_correlation(current_phase_todo).name = 'Solid (combined)';
    results_correlation(current_phase_todo).domain_length = voxel_number^(1/number_dimension) * voxel_size;
    results_correlation(current_phase_todo).voxel_size = voxel_size;
    results_correlation(current_phase_todo).volume_fraction = vf_solid;
    current_phase_todo=current_phase_todo+1;
    results_correlation(current_phase_todo).name = 'Porosity (combined, ideal wetting)';
    results_correlation(current_phase_todo).domain_length = voxel_number^(1/number_dimension) * voxel_size;
    results_correlation(current_phase_todo).voxel_size = voxel_size;
    results_correlation(current_phase_todo).volume_fraction = vf_pore_idealwetting;

    vf_pore_partialwetting = nvoxel.pore_partialwetting/n_voxel;
    results_correlation(current_phase_todo+1).volume_fraction = vf_pore_partialwetting;
    vf_air = 1-vf_solid-vf_pore_partialwetting;

    current_phase_todo=current_phase_todo+1;
    results_correlation(current_phase_todo).name = 'Porosity (combined, partial wetting)';
    results_correlation(current_phase_todo).domain_length = voxel_number^(1/number_dimension) * voxel_size;
    results_correlation(current_phase_todo).voxel_size = voxel_size;
    results_correlation(current_phase_todo).volume_fraction = vf_pore_partialwetting;

    current_phase_todo=current_phase_todo+1;
    results_correlation(current_phase_todo).name = 'Air saturation';
    results_correlation(current_phase_todo).domain_length = voxel_number^(1/number_dimension) * voxel_size;
    results_correlation(current_phase_todo).voxel_size = voxel_size;
    results_correlation(current_phase_todo).volume_fraction = vf_air;
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

if Run_in_commandline || isempty(Nanoporosity)
    rownames = {'label';'Volume fraction (label)'};
    C = zeros(length(rownames),number_phase_todo);
    C(1,:) = phaselabel_todo';
    C(2,:) = Vf_phase_label';
    Table_Volumefraction_perlabel = array2table(C,"VariableNames",phasename_todo',"RowNames",rownames);
    Results_Volumefraction.Table_Volumefraction_perlabel = Table_Volumefraction_perlabel;
    idx_wetting_is_relevant = [];
else
    % Statistics
    rownames = {'Mean';'Median';'Mode';'Standard deviation';'Minimum';'Maximum'};
    Table_nanoporositystats_perlabel = array2table(stats_nanoporosity',"VariableNames",phasename_todo',"RowNames",rownames);
    tmp = num2cell(stats_wetting);
    wetting_is_relevant = ones(number_phase_todo,1);
    for k=1:number_phase_todo
        if stats_nanoporosity(k,6)==0 % max nanoporosity=0: i.e. pure solid, wetting is irrelevant
            tmp(k,:) = {'n/a'};
            wetting_is_relevant(k)=0;
        end
    end
    idx_wetting_is_relevant = find(wetting_is_relevant);
    phaselabel_todo_wettingisrelevant = phaselabel_todo(idx_wetting_is_relevant);
    phasename_todo_wettingisrelevant = phasename_todo(idx_wetting_is_relevant);
    phasecolor_wettingisrelevant = phasecolor_todo(idx_wetting_is_relevant,:);
    stats_wettingisrelevant = stats_wetting(idx_wetting_is_relevant,:);
    Table_wetting_perlabel = array2table(tmp',"VariableNames",phasename_todo',"RowNames",rownames);

    % Per label
    rownames = {'label';'Volume fraction (label)';'Ratio (solid)';'Volume fraction (solid)';'Ratio (pore, ideal wetting)';'Ratio (pore, partial wetting)';'Volume fraction (pore, ideal wetting)';'Volume fraction (pore, partial wetting)'};
    C = zeros(length(rownames),number_phase_todo);
    C(1,:) = phaselabel_todo';
    C(2,:) = Vf_phase_label';
    C(3,:) = ratio_solid';
    C(4,:) = Vf_phase_solidpart';
    C(5,:) = ratio_poreidealwetting';
    C(6,:) = ratio_porepartialwetting';
    C(7,:) = Vf_phase_porouspart_idealwetting';
    C(8,:) = Vf_phase_porouspart_partiallwetting';
    Table_Volumefraction_perlabel = array2table(C,"VariableNames",phasename_todo',"RowNames",rownames);
    Results_Volumefraction.Table_Volumefraction_perlabel = Table_Volumefraction_perlabel;
    
    % Combined volume fraction
    if p.combined_todo
        names = {'category: only solid';'category: only pore';'category: mixed';'Solid';'Pore (ideal wetting)';'Pore (partial wetting)';'Air'};
        vf_combined = num2cell([vf_onlysolid;vf_onlypore;vf_mixed;vf_solid;vf_pore_idealwetting;vf_pore_partialwetting;vf_air]);
        Table_Volumefraction_combined = cell2table([names,vf_combined],"VariableNames",{'Combined phase' 'Volume fraction'});
        Results_Volumefraction.Table_Volumefraction_combined = Table_Volumefraction_combined;        
    end
end

if infovol.isbackground
    Table_Background = table({'Background'},background_label,vf_background,...
        'VariableNames',{'Phase' 'Label' 'Volume fraction'});%
end

%% SAVE TABLES
if opt.save.xls
    filename = 'Volume_fractions'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    sheet=1;
    DATA_writetable.sheet(sheet).name='Volume fraction';
    DATA_writetable.sheet(sheet).table=Table_Volumefraction_perlabel;
    if ~Run_in_commandline && ~isempty(Nanoporosity)
        sheet=sheet+1;
        DATA_writetable.sheet(sheet).name='Nanoporosity';
        DATA_writetable.sheet(sheet).table=Table_nanoporositystats_perlabel;
        sheet=sheet+1;
        DATA_writetable.sheet(sheet).name='Wetting';
        DATA_writetable.sheet(sheet).table=Table_wetting_perlabel;
        if p.combined_todo
            sheet=sheet+1;
            DATA_writetable.sheet(sheet).name='Combined';
            DATA_writetable.sheet(sheet).table=Table_Volumefraction_combined;
        end
        if infovol.isbackground
            sheet=sheet+1;
            DATA_writetable.sheet(sheet).name='Background';
            DATA_writetable.sheet(sheet).table=Table_Background;
        end
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY TABLES
fprintf('> Calculated on the whole domain:\n\n');
if infovol.isbackground
    fprintf('Background: label = %i, volume fraction = %1.3f\n',background_label,vf_background);
    fprintf('Subsequent volume fractions are calculated considering all pixels except for those belonging to the background\n\n');
end
if ~Run_in_commandline && ~isempty(Nanoporosity)
    fprintf('  Nanoporosity per phase:\n');
    disp(Table_nanoporositystats_perlabel)
    fprintf('  Wetting per phase:\n');
    disp(Table_wetting_perlabel)
end
fprintf('  Volume fractions and ratios per phase:\n');
disp(Table_Volumefraction_perlabel)
if p.combined_todo
    fprintf('  Volume fractions from the combined analysis of all phases:\n');
    fprintf('     Categories: voxels are either purely solid (nanoporosity=0), purely porous (nanoporosity=1), or mixed (nanoporosity ]0,1[. Sum(category)=1\n');
    fprintf('     Solid: mean(1-Nanoporosity), Pore(ideal wetting): mean(Nanoporosity), Pore(partial wetting): mean(Nanoporosity.*wetting)\n');
    fprintf('        Solid + Pore(ideal wetting) = 1\n');
    fprintf('        Solid + Pore(partial wetting) + Air = 1\n');
    disp(Table_Volumefraction_combined)
end
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%% PLOT VOLUME FRACTIONS
pars.array_name = 'Volume fraction';
pars.array_unit = '';
pars.plots = {'bar','pie'};
pars.pie_in_percents = true;
pars.inputfilename = infovol.filename;
plot_avg_perphase(phasename_todo,phasecolor_todo,Vf_phase_label,pars,opt.format,opt.save,Current_folder,'Volume_fractions_labels');

if ~Run_in_commandline && p.combined_todo
    pars.array_name = 'Volume fraction (category)';
    pars.array_unit = '';
    pars.plots = {'bar','pie'};
    pars.pie_in_percents = true;
    pars.inputfilename = infovol.filename;
    plot_avg_perphase([{'Pure solid'} {'Mixed'} {'Pure pore'} ],[1 0 0;1 0 1;0 0 1],[vf_onlysolid vf_mixed vf_onlypore],pars,opt.format,opt.save,Current_folder,'Volume_fractions_categories');

    pars.array_name = 'Volume fraction (state)';
    pars.array_unit = '';
    pars.plots = {'bar','pie'};
    pars.pie_in_percents = true;
    pars.inputfilename = infovol.filename;
    col_solid = [0.4902 0.8549 0.3451];
    col_liquid_partial = [0.3647 0.8863 0.9059];
    col_liquid_ideal = [0.0660    0.4430    0.7450];
    col_air = [0.85 0.85 0.85];
    plot_avg_perphase([{'Solid'} {'Liquid'} {'Air'}],[col_solid; col_liquid_partial; col_air],[vf_solid vf_pore_partialwetting vf_air],pars,opt.format,opt.save,Current_folder,'Volume_fractions_state');
end

%% PLOT NANOPOROSITY AND WETTING
if ~Run_in_commandline && ~isempty(Nanoporosity)
    pars.round_value = 2;
    pars.smooth_cumulative_fct = true;    
    pars.plots = {'bar(mean and std)','pdf','cdf'}; % {'bar'}, {'bar(mean and std)'}, {'pdf'}, or {'cdf'}, or a combination (e.g., {'bar','pdf'})
    pars.inputfilename = infovol.filename;

    pars.array_name = 'Nanoporosity';
    pars.array_unit = '';
    plot_distribution_perphase(Labels,phaselabel_todo,phasename_todo,phasecolor_todo,Nanoporosity,stats_nanoporosity,pars,opt.format,opt.save,Current_folder,'Nanoporosity_labelwise');

    if ~strcmp(infovol.partial_wetting_representation,'Ideal')
        pars.array_name = 'Wetting';
        pars.array_unit = '';
        plot_distribution_perphase(Labels,phaselabel_todo_wettingisrelevant,phasename_todo_wettingisrelevant,phasecolor_wettingisrelevant,Wetting,stats_wettingisrelevant,pars,opt.format,opt.save,Current_folder,'Wetting_labelwise');
    end
end

%%
%% ALONG DIRECTIONS
%%
if number_dimension==2
    p.plotdirections(3)=0;
end
if sum(p.plotdirections)>0
    group(1).inputfilename = infovol.filename;

    group(1).yaxis_name = 'Volume fraction'; 
    group(1).yaxis_unit = '';
    group(1).yaxis_round = 3;    
    group(1).filename = 'Volume_fractions_along_axis';
    group(1).mean = Vf_phase_label;

    if ~Run_in_commandline && ~isempty(Nanoporosity)
        group(2).yaxis_name = 'Nanoporosity';
        group(2).yaxis_unit = '';
        group(2).yaxis_round = 3;
        group(2).filename = 'Nanoporosity_along_axis';
        group(2).mean = stats_nanoporosity(:,1);
        group(2).std = stats_nanoporosity(:,4);
        group(3).yaxis_name = 'Wetting';
        group(3).yaxis_unit = '';
        group(3).yaxis_round = 3;
        group(3).filename = 'Wetting_along_axis';
        group(3).mean = stats_wetting(:,1);
        group(3).std = stats_wetting(:,4);
        group(3).relevant = idx_wetting_is_relevant;
        if p.combined_todo
            group(4).yaxis_name = 'Volume fraction';
            group(4).yaxis_unit = '';
            group(4).yaxis_round = 3;
            group(4).filename = 'Volume_fractions_combined_along_axis';
            group(4).mean = [vf_solid vf_pore_idealwetting vf_pore_partialwetting vf_air];
        end
    end    

    %% CALCULATE
    dd = 0;
    for d=1:3
        if p.plotdirections(d)
            dd = dd+1;

            direction(dd).name = infovol.directionname(d);
            direction(dd).xaxis = ([1:1:sz(d)]-0.5)*infovol.voxelsize;
            direction(dd).xaxis_unit = infovol.unit;

            for x=1:sz(d)
                % Slices
                if d==1
                    sl_label = squeeze(Labels(x,:,:));
                elseif d==2
                    sl_label = squeeze(Labels(:,x,:));
                elseif d==3
                    sl_label = Labels(:,:,x);
                end
                if ~Run_in_commandline && ~isempty(Nanoporosity)
                    if d==1
                        sl_Nanoporosity = squeeze(Nanoporosity(x,:,:));
                        sl_Wetting = squeeze(Wetting(x,:,:));
                    elseif d==2
                        sl_Nanoporosity = squeeze(Nanoporosity(:,x,:));
                        sl_Wetting = squeeze(Wetting(:,x,:));
                    elseif d==3
                        sl_Nanoporosity = Nanoporosity(:,:,x);
                        sl_Wetting = Wetting(:,:,x);
                    end
                else
                    sl_Nanoporosity = [];
                    sl_Wetting = [];
                end

                voxel_number = numel(sl_label);

                % Combined phase
                if p.combined_todo
                    nvoxel_.solid = 0;
                    nvoxel_.pore_idealwetting = 0;
                    nvoxel_.pore_partialwetting = 0;
                end
                % Background
                n_voxel_background = 0;
                if infovol.isbackground
                    n_voxel_background = sum(sum(sl_label==background_label));
                    vf_background = n_voxel_background/voxel_number;
                end
                n_voxel = numel(sl_label) - n_voxel_background;

                for current_phase_todo=1:number_phase_todo
                    if Run_in_commandline
                        group(1).direction(dd).array(current_phase_todo).vals(x) = sum(sum( sl_label==phaselabel_todo(current_phase_todo,1) )) / n_voxel; % 2D array (4xn) for avg, std, min max
                        group(1).direction(dd).array(current_phase_todo).name = char(phasename_todo(current_phase_todo));
                        group(1).direction(dd).array(current_phase_todo).color = phasecolor_todo(current_phase_todo,:);
                        group(1).direction(dd).array(current_phase_todo).linestyle = '-';
                    else
                        [st_,vf_,r_,n_] = Charact_Volumefractions_algorithm(sl_label, phaselabel_todo(current_phase_todo,1), n_voxel, sl_Nanoporosity, sl_Wetting);
                        group(1).direction(dd).array(current_phase_todo).vals(x) = vf_.phase_label;
                        group(1).direction(dd).array(current_phase_todo).name = char(phasename_todo(current_phase_todo));
                        group(1).direction(dd).array(current_phase_todo).color = phasecolor_todo(current_phase_todo,:);
                        group(1).direction(dd).array(current_phase_todo).linestyle = '-';

                        if ~isempty(Nanoporosity)
                            group(2).direction(dd).array(current_phase_todo).vals(x,1) = st_.nanoporosity(1);
                            group(2).direction(dd).array(current_phase_todo).vals(x,2) = st_.nanoporosity(4);
                            group(2).direction(dd).array(current_phase_todo).vals(x,3) = st_.nanoporosity(5);
                            group(2).direction(dd).array(current_phase_todo).vals(x,4) = st_.nanoporosity(6);
                            group(2).direction(dd).array(current_phase_todo).name = char(phasename_todo(current_phase_todo));
                            group(2).direction(dd).array(current_phase_todo).color = phasecolor_todo(current_phase_todo,:);

                            group(3).direction(dd).array(current_phase_todo).vals(x,1) = st_.wetting(1);
                            group(3).direction(dd).array(current_phase_todo).vals(x,2) = st_.wetting(4);
                            group(3).direction(dd).array(current_phase_todo).vals(x,3) = st_.wetting(5);
                            group(3).direction(dd).array(current_phase_todo).vals(x,4) = st_.wetting(6);
                            group(3).direction(dd).array(current_phase_todo).name = char(phasename_todo(current_phase_todo));
                            group(3).direction(dd).array(current_phase_todo).color = phasecolor_todo(current_phase_todo,:);

                            if p.combined_todo
                                nvoxel_.solid = nvoxel_.solid + n_.solid;
                                nvoxel_.pore_idealwetting = nvoxel_.pore_idealwetting + n_.pore_idealwetting;
                                nvoxel_.pore_partialwetting = nvoxel_.pore_partialwetting + n_.pore_partialwetting;
                            end
                        end

                    end

                    if infovol.isbackground
                        group(1).direction(dd).array(current_phase_todo+1).vals(x) = vf_background;
                        group(1).direction(dd).array(current_phase_todo+1).name = 'Background';
                        group(1).direction(dd).array(current_phase_todo+1).color = [0.5 0.5 0.5];
                        group(1).direction(dd).array(current_phase_todo+1).linestyle = '--'; % not used if 2D array
                    end

                end
                if ~Run_in_commandline
                    if p.combined_todo
                        % Porosity, solid, and air
                        vf_pore_idealwetting_ = nvoxel_.pore_idealwetting/n_voxel;
                        vf_solid_ = nvoxel_.solid/n_voxel;
                        vf_pore_partialwetting_ = nvoxel_.pore_partialwetting/n_voxel;
                        vf_air_ = 1 - vf_solid_ - vf_pore_partialwetting_;
                        group(4).direction(dd).array(1).vals(x) = vf_solid_; group(4).direction(dd).array(2).vals(x) = vf_pore_idealwetting_; group(4).direction(dd).array(3).vals(x) = vf_pore_partialwetting_; group(4).direction(dd).array(4).vals(x) = vf_air_;
                        group(4).direction(dd).array(1).name = 'Solid'; group(4).direction(dd).array(2).name = 'Liquid (ideal wetting)'; group(4).direction(dd).array(3).name = 'Liquid (partial wetting)'; group(4).direction(dd).array(4).name = 'Air';
                        group(4).direction(dd).array(1).color = col_solid; group(4).direction(dd).array(2).color = col_liquid_ideal; group(4).direction(dd).array(3).color = col_liquid_partial; group(4).direction(dd).array(4).color = col_air;
                        group(4).direction(dd).array(1).linestyle = '-'; group(4).direction(dd).array(2).linestyle = ':'; group(4).direction(dd).array(3).linestyle = '-'; group(4).direction(dd).array(4).linestyle = '-';
                    end
                end
            end

        end
    end

    %% TABLE AND PLOT
    plot_along_direction(group,direction,opt.format,opt.save,Current_folder);
end

%%
%% ERROR ANALYSIS: IMAGE RESOLUTION and REPRESENTATIVITY
%%

%% SETUP METRICS
p.scaling(p.scaling==1)=[];
if ~isempty(p.scaling) || p.RVE.number_RVE>0
    % Metric parameters
    pMET.fct_name = 'Volume fraction';
    pMET.inputfilename = infovol.filename;
    pMET.fractal_bertei = p.fractal_bertei.todo;
    pMET.metric(1).filename = 'Volume_fractions';
    pMET.metric(1).name = 'Volume fraction';
    pMET.metric(1).shortname_correlation = 'vf';
    pMET.metric(1).unit = [];
    pMET.metric(1).domain_name = phasename_todo;
    pMET.metric(1).domain_label = phaselabel_todo;
    pMET.metric(1).domain_color = phasecolor_todo;
    pMET.metric(1).result_initial = Vf_phase_label;
    pMET.metric(1).scaling_extrapolation = p.scaling_extrapolation;    
    if p.combined_todo
        pMET.metric(2).filename = 'Volume_fractions_combined';
        pMET.metric(2).name = 'Volume fraction (combined)';
        pMET.metric(2).shortname_correlation = 'vf';
        pMET.metric(2).unit = [];
        pMET.metric(2).domain_name = [{'Solid'}; {'Liquid (ideal wetting)'}; {'Liquid (partial wetting)'}; {'Air'}];
        pMET.metric(2).domain_color = [col_solid; col_liquid_ideal; col_liquid_partial; col_air];
        pMET.metric(2).result_initial = [vf_solid; vf_pore_idealwetting; vf_pore_partialwetting; vf_air];
        pMET.metric(2).scaling_extrapolation = p.scaling_extrapolation;       
    end
end

%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
if ~isempty(p.scaling)
    if isempty(Nanoporosity)
        [Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume] = ImageResolution_erroranalysis(pMET, [], [], [], [], [], Labels, voxel_size, p.scaling, Current_folder, opt, infovol, Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume);
    else
        [Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume] = ImageResolution_erroranalysis(pMET, Labels, Nanoporosity, Wetting, [], [], [], voxel_size, p.scaling, Current_folder, opt, infovol, Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume);
    end
end

%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
if p.RVE.number_RVE>0
    if isempty(Nanoporosity)
        [Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume] = RVE_main(p.RVE.RVE, pMET, [], [], [], [], [], Labels, voxel_size, Current_folder, opt, infovol, Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume); % Call main function
    else
        [Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume] = RVE_main(p.RVE.RVE, pMET, Labels, Nanoporosity, Wetting, [], [], [], voxel_size, Current_folder, opt, infovol, Results_Volumefraction, results_correlation, timedata_perphase, timedata_pervolume); % Call main function
    end
end

%%
%% END
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
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Volume_fraction_calculation_time', 'Volume fractions', infovol.filename, opt.format, opt.save);

%% SAVE CORRELATION
Correlation_folder = fullfile(infovol.volpath, 'Correlation');
if exist(Correlation_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Correlation_folder);
end

if ~isempty(infovol.sub)
    fullpath = fullfile(Correlation_folder,['Correlation_volume_fraction_' infovol.sub]);
else
    fullpath = fullfile(Correlation_folder,'Correlation_volume_fraction');
end
% if exist([fullpath '.mat'], 'file') == 2
%     tmp = load([fullpath '.mat']);
%     old = tmp.results_correlation;
%     fields_old = fieldnames(old);
%     current_length = length(results_correlation);
%     for k=1:1:length(old)
%         for k_field = 1:1:length(fields_old)
%             results_correlation(current_length+k).(char(fields_old(k_field))) = old(k).(char(fields_old(k_field)));
%         end
%     end
% end
save(fullpath,'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Summary_folder = fullfile(infovol.volpath, 'Summary');
    if exist(Summary_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Summary_folder);
    end
    if ~isempty(infovol.sub)
        save(fullfile(Summary_folder,['Results_volume_fraction_' infovol.sub]),'Results_Volumefraction')
    else
        save(fullfile(Summary_folder,'Results_volume_fraction'),'Results_Volumefraction')
    end
end

end