function [] = Charact_Tortuosity_main(Labels,Nanoporosity,Wetting,Bulkdiffusivity,Bulkconductivity,infovol,opt,p)
% Calculate transport properties
% Charact_Tortuosity_main(Labels, Nanoporosity, Wetting, Bulkdiffusivity, Bulkconductivity, infovol, opt, p) - when used with the MATBOX toolbox GUI
% or
% Charact_Tortuosity_main(Array, direction, val, name) - when used as a standalone function
% with: Array, a 3D array
%       direction: 1, 2 or 3
%       val: an integer if Array is a segmented volume, total porosity ]0,1[ or [] if Array is a bulk diffusivity map
%       (optional) name, a string: folder where result is saved
%       e.g.: - Charact_Volumefractions_main(<segmented array>, 3, 1); % Effective transport properties of label 1 along direction 3
%             - Charact_Volumefractions_main(<bulkk diffuisvity map>, 1, 0.4); % Effective transport properties along direction 1, knowing total porosity is 0.4
%             - Charact_Volumefractions_main(<bulkk diffuisvity map>, 1, []);  % Effective transport properties along direction 1, not knowing total porosity
%               Call help Call_TauFactor2_multiphase for more details.
% In both cases, Deff is calculated with an homogenization calculation performed by TauFactor2
% Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
% https://doi.org/10.21105/joss.05358.

% To do:
% - Along direction
% - image res: DONE but correlation wrong
% - RVE: DONE but correlation wrong
% - summary
% - end
% - correlation
% - connected volume fraction, all volume fraction
%   bar: 3*2 bars (dir 1 all volume fraction, dir 1 connected volume fraction, dir 2,...
%   Connection definition: must be connected both, only 1st, only last,
%      in The GUI
% - direction color in the GUI
% - anisotropy
% - lower scale diffusivity stats : DONE
%   avg, std, distribution : DONE
%   lower scale diffusivity variation along direction: ON GOING - ADD PER
%   SLICE CONSIDERING ALL VOLUME and not just combined
% - Add Macmullin number: DONE
% - default values for command line
% - background (see [Deff, Tau, Bruggeman, eps] = Tortuosity_algorithm(M,p)) in the meantime
%   add crop method to start with
% - add option to select only Mc, or Deff, or p


% Background:
% -RVE: on cropped volume only
% -everything else: symmetry/inpainting: no more background
%                   iteration: iteratve on k value for background

wettedpore_combined_color = [0 0 0];
solid_combined_color = [0.5 0.5 0.5];

%% DEFAULT VALUES IF RUN IN COMMAND LINE
expected_number_argument = 8;
nargin; % Number of input variable when the function is call

Run_in_commandline = false;
if nargin ~= expected_number_argument % Unexpected number of argument

    if nargin == 3 || nargin==4 % Case for function called as: Charact_Tortuosity_main(Array, direction, val, name). Standalone use.
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
            infovol.mainfolder = winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
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
        disp 'Error calling Charact_Tortuosity_main. Wrong number of argument.'
        help Charact_Tortuosity_main
        return
    end

end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
Current_folder = fullfile(infovol.volpath,'Transport',infovol.sub);
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
    dimension = 2;
else
    dimension = 3; % 3D case
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

%%
%% BACKGROUND
%%

% ? Here or in main ?

%%
%% ALGORITHM ON WHOLE VOLUME
%%

fprintf('> Calculated on the whole domain:\n\n');

Bulkdiffusivity(isinf(Bulkdiffusivity))=0;
Bulkdiffusivity(isnan(Bulkdiffusivity))=0;
Bulkdiffusivity(Bulkdiffusivity<0)=0;
Bulkdiffusivity(Bulkdiffusivity>1)=1;

Bulkconductivity(isinf(Bulkconductivity))=0;
Bulkconductivity(isnan(Bulkconductivity))=0;
Bulkconductivity(Bulkconductivity<0)=0;
Bulkconductivity(Bulkconductivity>1)=1;

% Direction
number_direction_todo = sum(p.tau_xyz(1:dimension));
directionname_todo = cell(number_direction_todo,1);
direction_todo = zeros(number_direction_todo,1);
current_direction_todo = 0;
for current_direction=1:dimension % Loop over all phases
    if p.tau_xyz(current_direction)
        current_direction_todo = current_direction_todo + 1;
        direction_todo(current_direction_todo,1) = current_direction;
        directionname_todo(current_direction_todo,1) = infovol.directionname(current_direction,1);
    end
end

index = 0;

%% PER LABEL
L = length(p.todo);
if L < number_phase
    p.todo = [p.todo; ones(number_phase-L)];
elseif L > number_phase
    p.todo = p.todo(1:number_phase);
end

number_phase_todo = 0;
phaselabel_todo = [];
timedata_perphase_name = [];
timedata_perphase = [];
if sum(p.todo)>0
    fprintf('  Per label...\n');
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

    timedata_perphase = zeros(number_phase_todo*number_direction_todo,3);
    timedata_perphase_name = cell(number_phase_todo*number_direction_todo,1);

    % Initialization (algorithm-specific)
    % Phase-wise
    Deff = zeros(number_phase_todo,number_direction_todo);
    Mc = zeros(number_phase_todo,number_direction_todo);
    Tau = zeros(number_phase_todo,number_direction_todo);
    Bruggeman =zeros(number_phase_todo,number_direction_todo);
    eps = zeros(number_phase_todo,number_direction_todo);

    % Start
    index = 0;
    for current_phase_todo=1:number_phase_todo % Loop over all phases
        for current_direction_todo = 1:number_direction_todo % Loop over all directions
            time_cpu_start_phase = cputime; % CPU start
            time_clock_start_phase = tic; % Stopwatch start
            index = index +1;

            direction = direction_todo(current_direction_todo);
            label = phaselabel_todo(current_phase_todo);

            % % Algorithm
            if Run_in_commandline
                foo = 1

            else
                [Deff(current_phase_todo,current_direction_todo),Mc(current_phase_todo,current_direction_todo),Tau(current_phase_todo,current_direction_todo),Bruggeman(current_phase_todo,current_direction_todo),eps(current_phase_todo,current_direction_todo)] = Call_TauFactor2_binary(Labels,direction,label);
                % % Correlation
                results_correlation(index).name = {[char(phasename_todo(current_phase_todo,1)) ', ' char(directionname_todo(current_direction_todo,1))]};
                results_correlation(index).Deff = Deff(current_phase_todo,current_direction_todo);
                results_correlation(index).Mc = Mc(current_phase_todo,current_direction_todo);
                results_correlation(index).Tau = Tau(current_phase_todo,current_direction_todo);
                results_correlation(index).Bruggeman = Bruggeman(current_phase_todo,current_direction_todo);
            end

            % Time
            timedata_perphase_name(index,1) = {[char(phasename_todo(current_phase_todo,1)) ', ' char(directionname_todo(current_direction_todo,1))]};
            timedata_perphase(index,1) = sum(sum(sum( Labels==label )));
            timedata_perphase(index,2) = cputime-time_cpu_start_phase; % CPU elapsed time
            timedata_perphase(index,3) = toc(time_clock_start_phase); % CPU elapsed time
        end
    end
end

%% ON COMBINED VOLUMES
index_combined = 0;
timedata_combined_name = {};
timedata_combined = [];
if p.pore_combined_todo

    % Statistics
    BWlog = Nanoporosity.*Wetting>0;
    BW_combined_pore = single(BWlog);
    vals = Bulkdiffusivity(BWlog);
    vals = round(double(vals),4);
    vals = reshape(vals,[1 numel(vals)]);
    if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
        stats_diffusivity_combined = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
    else
        stats_diffusivity_combined = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
    end

    Label_diffusivity = double(unique(Labels(BWlog)));
    idx_diffusivity = find(ismember(double(infovol.phaselabel_semantic),Label_diffusivity));
    stats_diffusivity_label = zeros(length(idx_diffusivity),6);
    for k=1:length(idx_diffusivity)
        label = infovol.phaselabel_semantic(idx_diffusivity(k));
        BWlog = Labels==label;
        vals = Bulkdiffusivity(BWlog);
        vals = round(double(vals),4);
        vals = reshape(vals,[1 numel(vals)]);
        if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
            st = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
        else
            st = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
        end
        stats_diffusivity_label(k,:) = st;
    end

    fprintf('  Wetted pore (combined), for the domain with Nanoporosity.*Wetting>0...\n');
    if strcmp(infovol.poretransport_representation,'Bruggeman exponent p (uniform)')
        fprintf('     Local diffusivity is (Nanoporosity.*Wetting).^Bruggeman exponent\n');
    end
    pore_Deff = zeros(number_direction_todo,1);
    pore_Mc = zeros(number_direction_todo,1);
    pore_Tau = zeros(number_direction_todo,1);
    pore_Bruggeman =zeros(number_direction_todo,1);
    pore_eps =zeros(number_direction_todo,1);
    for current_direction_todo = 1:number_direction_todo % Loop over all directions
        time_cpu_start_combined = cputime; % CPU start
        time_clock_start_combined = tic; % Stopwatch start
        index = index + 1;
        index_combined = index_combined + 1;

        direction = direction_todo(current_direction_todo);

        [pore_Deff(current_direction_todo),pore_Mc(current_direction_todo),pore_Tau(current_direction_todo),pore_Bruggeman(current_direction_todo),pore_eps(current_direction_todo),n_voxel] = Call_TauFactor2_pore(Nanoporosity,Wetting,Bulkdiffusivity,direction);

        % % Correlation
        results_correlation(index).name = {['Wetted pore (combined), ' char(directionname_todo(current_direction_todo,1))]};
        results_correlation(index).Deff = pore_Deff(current_direction_todo);
        results_correlation(index).Mc = pore_Mc(current_direction_todo);
        results_correlation(index).Tau = pore_Tau(current_direction_todo);
        results_correlation(index).Bruggeman = pore_Bruggeman(current_direction_todo);

        % Time
        timedata_combined_name(index_combined,1) = {['Wetted pore (combined), ' char(directionname_todo(current_direction_todo,1))]};
        timedata_combined(index_combined,1) = n_voxel;
        timedata_combined(index_combined,2) = cputime-time_cpu_start_combined; % CPU elapsed time
        timedata_combined(index_combined,3) = toc(time_clock_start_combined); % CPU elapsed time
    end
end

if p.solid_combined_todo

    % Statistics
    BWlog = Nanoporosity<1;
    BW_combined_solid = single(BWlog);
    vals = Bulkconductivity(BWlog);
    vals = round(double(vals),4);
    vals = reshape(vals,[1 numel(vals)]);
    if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
        stats_conductivity_combined = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
    else
        stats_conductivity_combined = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
    end

    Label_conductivity = double(unique(Labels(BWlog)));
    idx_conductivity = find(ismember(double(infovol.phaselabel_semantic),Label_conductivity));
    stats_conductivity_label = zeros(length(idx_conductivity),6);
    for k=1:length(idx_conductivity)
        label = infovol.phaselabel_semantic(idx_conductivity(k));
        BWlog = Labels==label;
        vals = Bulkconductivity(BWlog);
        vals = round(double(vals),4);
        vals = reshape(vals,[1 numel(vals)]);
        if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
            st = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
        else
            st = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
        end
        stats_conductivity_label(k,:) = st;
    end

    fprintf('  Solid (combined), for the domain with (1-Nanoporosity)>0...\n');
    if strcmp(infovol.solidtransport_representation,'Bruggeman exponent p (uniform)')
        fprintf('     Local conductivity is (1-Nanoporosity).^Bruggeman exponent\n');
    end
    solid_Deff = zeros(number_direction_todo,1);
    solid_Mc = zeros(number_direction_todo,1);
    solid_Tau = zeros(number_direction_todo,1);
    solid_Bruggeman =zeros(number_direction_todo,1);
    solid_eps =zeros(number_direction_todo,1);
    for current_direction_todo = 1:number_direction_todo % Loop over all directions
        time_cpu_start_combined = cputime; % CPU start
        time_clock_start_combined = tic; % Stopwatch start
        index = index + 1;
        index_combined = index_combined + 1;

        direction = direction_todo(current_direction_todo);

        [solid_Deff(current_direction_todo),solid_Mc(current_direction_todo),solid_Tau(current_direction_todo),solid_Bruggeman(current_direction_todo),solid_eps(current_direction_todo),n_voxel] = Call_TauFactor2_solid(Nanoporosity,Bulkconductivity,direction);

        % % Correlation
        results_correlation(index).name = {['Solid (combined), ' char(directionname_todo(current_direction_todo,1))]};
        results_correlation(index).Deff = solid_Deff(current_direction_todo);
        results_correlation(index).Mc = solid_Mc(current_direction_todo);
        results_correlation(index).Tau = solid_Tau(current_direction_todo);
        results_correlation(index).Bruggeman = solid_Bruggeman(current_direction_todo);

        % Time
        timedata_combined_name(index_combined,1) = {['Solid (combined), ' char(directionname_todo(current_direction_todo,1))]};
        timedata_combined(index_combined,1) = n_voxel;
        timedata_combined(index_combined,2) = cputime-time_cpu_start_combined; % CPU elapsed time
        timedata_combined(index_combined,3) = toc(time_clock_start_combined); % CPU elapsed time
    end
end

%% TABLES
% Time
if ~isempty(timedata_perphase_name) && ~isempty(timedata_combined_name)
    timedata_name = [timedata_perphase_name; timedata_combined_name];
    timedata_vals = [timedata_perphase;timedata_combined];
elseif ~isempty(timedata_perphase_name)
    timedata_name = timedata_perphase_name;
    timedata_vals = timedata_perphase;
elseif ~isempty(timedata_combined_name)
    timedata_name = timedata_combined_name;
    timedata_vals = timedata_combined;
end
Table_time = table(timedata_name, timedata_vals(:,1),timedata_vals(:,2),timedata_vals(:,3),...
    'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Transport.Table_time = Table_time; % Save in main table result

if Run_in_commandline
    foo = 1;
else
    domainnames = [];
    domaincolors = [];
    Teff = []; Mceff = []; Taus = []; Ps = []; Epss = [];
    if number_phase_todo > 0
        domainnames = [domainnames; phasename_todo];
        domaincolors = [domaincolors; phasecolor_todo];
        Teff = [Teff;Deff];
        Mceff = [Mceff;Mc];
        Taus = [Taus;Tau];
        Ps = [Ps;Bruggeman];
        Epss = [Epss;eps];
    end
    if p.pore_combined_todo
        domainnames = [domainnames; {'Wetted pore (combined)'}];
        domaincolors = [domaincolors; wettedpore_combined_color];
        Teff = [Teff;pore_Deff'];
        Mceff = [Mceff;pore_Mc'];
        Taus = [Taus;pore_Tau'];
        Ps = [Ps;pore_Bruggeman'];
        Epss = [Epss;pore_eps'];
    end
    if p.solid_combined_todo
        domainnames = [domainnames; {'Solid (combined)'}];
        domaincolors = [domaincolors; solid_combined_color];
        Teff = [Teff;solid_Deff'];
        Mceff = [Mceff;solid_Mc'];
        Taus = [Taus;solid_Tau'];
        Ps = [Ps;solid_Bruggeman'];
        Epss = [Epss;solid_eps'];
    end

    Table_effectivetransport = array2table(Teff,"VariableNames",directionname_todo',"RowNames",domainnames);
    Table_MacMullin = array2table(Mceff,"VariableNames",directionname_todo',"RowNames",domainnames);
    Table_Tortuosity = array2table(Taus,"VariableNames",directionname_todo',"RowNames",domainnames);
    Table_Bruggeman = array2table(Ps,"VariableNames",directionname_todo',"RowNames",domainnames);
    Table_Vf = array2table(Epss,"VariableNames",directionname_todo',"RowNames",domainnames);

    Results_Transport.Table_effectivetransport = Table_effectivetransport;
    Results_Transport.Table_MacMullin = Table_MacMullin;
    Results_Transport.Table_Tortuosity = Table_Tortuosity;
    Results_Transport.Table_Bruggeman = Table_Bruggeman;
    Results_Transport.Table_Vf = Table_Vf;

    % Statistics
    if p.pore_combined_todo || p.solid_combined_todo
        rownames = {'Mean';'Median';'Mode';'Standard deviation';'Minimum';'Maximum'};
        if p.pore_combined_todo
            Table_diffusivity_perdomain = array2table([stats_diffusivity_label' stats_diffusivity_combined'],"VariableNames",[infovol.phasename(idx_diffusivity)' 'Wetted pore (combined)'],"RowNames",rownames);
        end
        if p.solid_combined_todo
            Table_conductivity_perdomain = array2table([stats_conductivity_label' stats_conductivity_combined'],"VariableNames",[infovol.phasename(idx_conductivity)' 'Solid (combined)'],"RowNames",rownames);
        end
    end
end


%% SAVE TABLES
if opt.save.xls
    filename = 'Effective_transport'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Effective transport';
    DATA_writetable.sheet(1).table=Table_effectivetransport;
    DATA_writetable.sheet(2).name='MacMullin number';
    DATA_writetable.sheet(2).table=Table_MacMullin;
    DATA_writetable.sheet(3).name='Tortuosity factor';
    DATA_writetable.sheet(3).table=Table_Tortuosity;
    DATA_writetable.sheet(4).name='Bruggeman exponent';
    DATA_writetable.sheet(4).table=Table_Bruggeman;
    DATA_writetable.sheet(5).name='Volume fractions';
    DATA_writetable.sheet(5).table=Table_Vf;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)

    if p.pore_combined_todo || p.solid_combined_todo
        clear DATA_writetable
        filename = 'Bulk_transport_statistics'; % Filename without extension
        sheet = 0;
        if p.pore_combined_todo
            sheet = sheet + 1;
            DATA_writetable.sheet(sheet).name='Diffusivity';
            DATA_writetable.sheet(sheet).table=Table_diffusivity_perdomain;
        end
        if p.solid_combined_todo
            sheet = sheet + 1;
            DATA_writetable.sheet(sheet).name='Conductivity';
            DATA_writetable.sheet(sheet).table=Table_conductivity_perdomain;
        end
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
end

%% DISPLAY TABLES
disp(' ')
if p.pore_combined_todo
    fprintf('  Bulk diffusivity per domain:\n');
    disp(Table_diffusivity_perdomain)
end
if  p.solid_combined_todo
    fprintf('  Bulk conductivity per domain:\n');
    disp(Table_conductivity_perdomain)
end
fprintf('  Effective transport:\n');
disp(Table_effectivetransport)
fprintf('  MacMulling number:\n');
disp(Table_MacMullin)
fprintf('  Tortuosity factor:\n');
disp(Table_Tortuosity)
fprintf('  Bruggeman exponent:\n');
disp(Table_Bruggeman)
fprintf('  Volume fractions:\n');
disp(Table_Vf)
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%% PLOT BULK TRANSPORT
if p.pore_combined_todo || p.solid_combined_todo
    pars.round_value = 3;
    pars.smooth_cumulative_fct = true;
    pars.plots = {'bar(mean and std)','pdf','cdf'}; % {'bar'}, {'bar(mean and std)'}, {'pdf'}, or {'cdf'}, or a combination (e.g., {'bar','pdf'})
    pars.inputfilename = infovol.filename;

    if p.pore_combined_todo
        pars.array_name = 'Bulk diffusivity';
        pars.array_unit = '';
        plot_distribution_perphase(Labels,infovol.phaselabel(idx_diffusivity),infovol.phasename(idx_diffusivity),infovol.phasecolor(idx_diffusivity,:),Bulkdiffusivity,stats_diffusivity_label,pars,opt.format,opt.save,Current_folder,'Bulk diffusivity (label)');
        if strcmp(infovol.poretransport_representation,'Diffusivity D (heterogeneous)')
            plot_distribution_perphase(BW_combined_pore,1,{'Wetted pore (combined)'},wettedpore_combined_color,Bulkdiffusivity,stats_diffusivity_combined,pars,opt.format,opt.save,Current_folder,'Bulk diffusivity (combined)');
        end
    end
    if p.solid_combined_todo
        pars.array_name = 'Bulk conductivity';
        pars.array_unit = '';
        plot_distribution_perphase(Labels,infovol.phaselabel(idx_conductivity),infovol.phasename(idx_conductivity),infovol.phasecolor(idx_conductivity,:),Bulkconductivity,stats_conductivity_label,pars,opt.format,opt.save,Current_folder,'Bulk conductivity (label)');
        if strcmp(infovol.poretransport_representation,'Diffusivity D (heterogeneous)')
            plot_distribution_perphase(BW_combined_solid,1,{'Solid (combined)'},solid_combined_color,Bulkconductivity,stats_conductivity_combined,pars,opt.format,opt.save,Current_folder,'Bulk conductivity (combined)');
        end
    end
end

%% PLOT EFFECTIVE TRANSPORT
pars.name = 'Effective transport properties';
pars.title = 'Effective transport properties D_{eff}/D_{bulk} = \epsilon / \tau = \epsilon^{p} = 1/M_{c}';
pars.array_name = {'Effective transport D_{eff}/D_{bulk}', 'MacMullin number M_{c}', 'Volume fractions \epsilon', 'Tortuosity factor \tau','Bruggeman exponent p'};
pars.vals = {Teff,Mceff,Epss,Taus,Ps};
pars.array_unit = {'','','','',''};
pars.yaxis_range = [[0,1];[0, Inf];[0 1];[0, Inf];[0, Inf]];
pars.layout = [2,3];
pars.legend = directionname_todo;
pars.inputfilename = infovol.filename;
plot_bar_multiple(domainnames,pars,opt.format,opt.save,Current_folder,'Effective_transport','Transport');

%%
%% ALONG DIRECTIONS: BULK TRANSPORT
%%

clear direction
p.plotdirections = [1 1 1];
if dimension==2
    p.plotdirections(3)=0;
end

if p.pore_combined_todo
    if sum(p.plotdirections)>0
        group(1).inputfilename = infovol.filename;
        group(1).yaxis_name = 'Bulk diffusivity';
        group(1).yaxis_unit = '';
        group(1).yaxis_round = 3;
        group(1).filename = 'Bulk_diffusivity_along_axis';
        group(1).mean = stats_diffusivity_combined(:,1);
        group(1).std = stats_diffusivity_combined(:,4);

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
                        sl_label = squeeze(BW_combined_pore(x,:,:));
                        sl_bulk = squeeze(Bulkdiffusivity(x,:,:));
                    elseif d==2
                        sl_label = squeeze(BW_combined_pore(:,x,:));
                        sl_bulk = squeeze(Bulkdiffusivity(:,x,:));
                    elseif d==3
                        sl_label = BW_combined_pore(:,:,x);
                        sl_bulk = squeeze(Bulkdiffusivity(:,:,x));
                    end

                    vals = sl_bulk(sl_label==1);
                    vals = round(double(vals),4);
                    vals = reshape(vals,[1 numel(vals)]);
                    if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
                        st = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
                    else
                        st = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
                    end
                    group(1).direction(dd).array(1).vals(x,1) = st(1);
                    group(1).direction(dd).array(1).vals(x,2) = st(4);
                    group(1).direction(dd).array(1).vals(x,3) = st(5);
                    group(1).direction(dd).array(1).vals(x,4) = st(6);
                    group(1).direction(dd).array(1).name = 'Wetted pore (combined)';
                    group(1).direction(dd).array(1).color = wettedpore_combined_color;
                  
                end

            end
        end

        %% TABLE AND PLOT
        plot_along_direction(group,direction,opt.format,opt.save,Current_folder);
    end
end

clear group
if p.solid_combined_todo
    if sum(p.plotdirections)>0
        group(1).inputfilename = infovol.filename;
        group(1).yaxis_name = 'Bulk conductivity';
        group(1).yaxis_unit = '';
        group(1).yaxis_round = 3;
        group(1).filename = 'Bulk_conductivity_along_axis';
        group(1).mean = stats_conductivity_combined(:,1);
        group(1).std = stats_conductivity_combined(:,4);

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
                        sl_label = squeeze(BW_combined_solid(x,:,:));
                        sl_bulk = squeeze(Bulkconductivity(x,:,:));
                    elseif d==2
                        sl_label = squeeze(BW_combined_solid(:,x,:));
                        sl_bulk = squeeze(Bulkconductivity(:,x,:));
                    elseif d==3
                        sl_label = BW_combined_solid(:,:,x);
                        sl_bulk = squeeze(Bulkconductivity(:,:,x));
                    end

                    vals = sl_bulk(sl_label==1);
                    vals = round(double(vals),4);
                    vals = reshape(vals,[1 numel(vals)]);
                    if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
                        st = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
                    else
                        st = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
                    end
                    group(1).direction(dd).array(1).vals(x,1) = st(1);
                    group(1).direction(dd).array(1).vals(x,2) = st(4);
                    group(1).direction(dd).array(1).vals(x,3) = st(5);
                    group(1).direction(dd).array(1).vals(x,4) = st(6);
                    group(1).direction(dd).array(1).name = 'Solid (combined)';
                    group(1).direction(dd).array(1).color = solid_combined_color;
                  
                end

            end
        end

        %% TABLE AND PLOT
        plot_along_direction(group,direction,opt.format,opt.save,Current_folder);
    end

end

%%
%% ALONG DIRECTIONS: EFFECTIVE TRANSPORT
%%
% 
% keyboard
% 
% if dimension == 2
%     p.calculate_along(3) = 0;
% end
% section_dir_todo = find(p.calculate_along>1);
% n_section_dir = length(section_dir_todo);
% 
% 
% 
% if ~isempty(section_dir_todo)
%     %% SECTIONS
%     xyz0_section = zeros(n_section_dir,3);
%     xyz1_section = zeros(n_section_dir,3);
%     for kdir = 1:n_section_dir
%         alongdir = section_dir_todo(kdir);
%         number_section = p.calculate_along(alongdir);
% 
%         % HERE
%         if alongdir == 1
%             x0 = round(linspace(1,sz(alongdir),number_section+1));
%             x1 = x0(2:end-1)-1;
%             x0(end)=[]; x1 = [x1 sz(alongdir)];
% 
%             y0 = 1; y1 = sz(2);
%             if dimension == 3
%                 z0 = 1; z1 = sz(3);
%             else
%                 z0 = 1; z1 = 1;
%             end
% 
%         elseif alongdir == 2
%             y0 = round(linspace(1,sz(alongdir),number_section+1));
%             y1 = y0(2:end-1)-1;
%             y0(end)=[]; y1 = [x1 sz(alongdir)];            
% 
%         elseif alongdir == 3
% 
% 
%         end
%     end
% 
% 
% 
%     %% ON LABELS
%     if sum(p.todo)>0
%         fprintf('  Calculated section per section along axis, per label...\n');
%         for kdir = 1:length(section_dir_todo)
%             alongdir = section_dir_todo(kdir);
%             % Initialization (algorithm-specific)
%             number_section = p.calculate_along(alongdir);
%             Deff_section = zeros(number_phase_todo,number_direction_todo,number_section);
%             Mc_section = zeros(number_phase_todo,number_direction_todo,number_section);
%             Tau_section = zeros(number_phase_todo,number_direction_todo,number_section);
%             Bruggeman_section =zeros(number_phase_todo,number_direction_todo,number_section);
%             eps_section = zeros(number_phase_todo,number_direction_todo,number_section);
% 
%             % Start
%             for current_phase_todo=1:number_phase_todo % Loop over all phases
%                 for current_direction_todo = 1:number_direction_todo % Loop over all directions
%                     for section = 1:number_section
%                         time_cpu_start_phase = cputime; % CPU start
%                         time_clock_start_phase = tic; % Stopwatch start
%                         index = index +1;
% 
%                         direction = direction_todo(current_direction_todo);
%                         label = phaselabel_todo(current_phase_todo);
% 
%                         Labels_section = 
% 
%                         % Algorithm
%                         [Deff_section(current_phase_todo,current_direction_todo,section),Mc_section(current_phase_todo,current_direction_todo,section),Tau_section(current_phase_todo,current_direction_todo,section),Bruggeman_section(current_phase_todo,current_direction_todo,section),eps_section(current_phase_todo,current_direction_todo,section)] = Call_TauFactor2_binary(Labels_section,direction,label);
% 
%                         % Time
%                         timedata_perphase_name(index,1) = {[char(phasename_todo(current_phase_todo,1)) ', ' char(directionname_todo(current_direction_todo,1)) ,', section']};
%                         timedata_perphase(index,1) = sum(sum(sum( Labels_section==label )));
%                         timedata_perphase(index,2) = cputime-time_cpu_start_phase; % CPU elapsed time
%                         timedata_perphase(index,3) = toc(time_clock_start_phase); % CPU elapsed time
%                     end
%                 end
%             end
%         end
%     end
% 
%     %% ON COMBINED DOMAINS
%
%end

%% SETUP METRICS
p.scaling(p.scaling==1)=[];
if ~isempty(p.scaling) || p.RVE.number_RVE>0
    % Metric parameters
    pMET.fct_name = 'Transport';
    pMET.inputfilename = infovol.filename;
    pMET.fractal_bertei = p.fractal_bertei.todo;
    pMET.phases_todo = p.todo;
    pMET.phaselabel_todo = phaselabel_todo;
    pMET.pore_combined_todo = p.pore_combined_todo;
    pMET.solid_combined_todo = p.solid_combined_todo;
    pMET.direction_todo = direction_todo;

    kmetric = 0;
    for dir = 1:number_direction_todo
        kmetric = kmetric+1;
        pMET.metric(kmetric).filename = ['Tortuosity_factor_ax' num2str(direction_todo(dir))];
        pMET.metric(kmetric).name = ['Tortuosity factor along ' char(directionname_todo(dir))];
        pMET.metric(kmetric).shortname_correlation = ['Tau' num2str(direction_todo(dir))];
        pMET.metric(kmetric).unit = [];
        pMET.metric(kmetric).domain_name = domainnames;
        pMET.metric(kmetric).domain_color = domaincolors;
        pMET.metric(kmetric).result_initial = Taus(:,dir);
        pMET.metric(kmetric).scaling_extrapolation = p.scaling_extrapolation.tau;
    
        kmetric = kmetric+1;
        pMET.metric(kmetric).filename = ['Bruggeman_exponent_ax' num2str(direction_todo(dir))];
        pMET.metric(kmetric).name = ['Bruggeman exponent along ' char(directionname_todo(dir))];
        pMET.metric(kmetric).shortname_correlation = ['p' num2str(direction_todo(dir))];
        pMET.metric(kmetric).unit = [];
        pMET.metric(kmetric).domain_name = domainnames;
        pMET.metric(kmetric).domain_color = domaincolors;
        pMET.metric(kmetric).result_initial = Ps(:,dir); 
        pMET.metric(kmetric).scaling_extrapolation = p.scaling_extrapolation.p;

        kmetric = kmetric+1;
        pMET.metric(kmetric).filename = ['Effective_transport_ax' num2str(direction_todo(dir))];
        pMET.metric(kmetric).name = ['Effective transport along ' char(directionname_todo(dir))];
        pMET.metric(kmetric).shortname_correlation = ['Deff' num2str(direction_todo(dir))];
        pMET.metric(kmetric).unit = [];
        pMET.metric(kmetric).domain_name = domainnames;
        pMET.metric(kmetric).domain_color = domaincolors;
        pMET.metric(kmetric).result_initial = Teff(:,dir);    
        pMET.metric(kmetric).scaling_extrapolation = p.scaling_extrapolation.Deff;
    end

end

%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
if ~isempty(p.scaling)
    if isempty(Nanoporosity)
        [Results_Transport, results_correlation, timedata_perphase, timedata_combined] = ImageResolution_erroranalysis(pMET, [], [], [], [], [], Labels, voxel_size, p.scaling, Current_folder, opt, infovol, Results_Transport, results_correlation, timedata_perphase, timedata_combined);
    else
        [Results_Transport, results_correlation, timedata_perphase, timedata_combined] = ImageResolution_erroranalysis(pMET, Labels, Nanoporosity, Wetting, Bulkdiffusivity,Bulkconductivity, [], voxel_size, p.scaling, Current_folder, opt, infovol, Results_Transport, results_correlation, timedata_perphase, timedata_combined);
    end
end


%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
if p.RVE.number_RVE>0
    if isempty(Nanoporosity)
        [Results_Transport, results_correlation, timedata_perphase, timedata_combined] = RVE_main(p.RVE.RVE, pMET, [], [], [], [], [], Labels, voxel_size, Current_folder, opt, infovol, Results_Transport, results_correlation, timedata_perphase, timedata_combined); % Call main function
    else
        [Results_Transport, results_correlation, timedata_perphase, timedata_combined] = RVE_main(p.RVE.RVE, pMET, Labels, Nanoporosity, Wetting, Bulkdiffusivity, Bulkconductivity, [], voxel_size, Current_folder, opt, infovol, Results_Transport, results_correlation, timedata_perphase, timedata_combined); % Call main function
    end
end

%% SAVE RESULTS
if opt.save.mat
    Summary_folder = fullfile(infovol.volpath, 'Summary');
    if exist(Summary_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Summary_folder);
    end
    if ~isempty(infovol.sub)
        save(fullfile(Summary_folder,['Results_Transport_' infovol.sub]),'Results_Transport')
    else
        save(fullfile(Summary_folder,'Results_Transport'),'Results_Transport')
    end
end

end