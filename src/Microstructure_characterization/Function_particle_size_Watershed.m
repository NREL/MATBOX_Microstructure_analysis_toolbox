function [] = Function_particle_size_Watershed(Phase_microstructure,infovol,opt,p)
% Calculate Particle size with watershed (immersion) method (C-PSD)
% Function_particle_size_Watershed(Phase_microstructure,infovol,opt,p) - when use with the toolbox
% or
% Function_particle_size_Watershed(Phase_microstructure,*) - when use as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 4 % Case for function called as: Function_particle_size_Watershed(Phase_microstructure,*)
        
        
    else % Incorrect number of argument
        disp 'Error calling Function_particle_size_Watershed. Wrong number of argument.'
        help Function_particle_size_Watershed
        return
    end

else
    % Read parameters
    cpsd_refining = p.correctoversegmentation_cpsd;
    sizethreshold = p.correctoversegmentation_sizethreshold;
    verbose = p.verbose;   
    visualize_2D = false;    
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Diameter_Watershed' separator];
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
% Note: those are not algorithm paramter, but post-processing smoothing parameter to get a better looking distribution function

% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;

%%
%% ALGORITHM ON WHOLE VOLUME
%%

disp '    PARTICLE SIZE - DISCRETE PARTICLE SIZE DISTRIBUTION (WATERSHED) METHOD';
disp '    ----------------------------------------------------------------------';
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
PSD_results = zeros(number_phase_todo,5); % min, mean, max, std, std%

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
        if verbose
            fprintf('Current phase: %s\n',char(phasename_todo(current_phase_todo,1)))
        end
        % Watershed algorithm
        [Label_lake] = Function_particle_size_Watershed_algorithm(binary_phase,verbose);

        % Oversegmentation correction
        [Label_lake] = Function_oversegmentation_correction(binary_phase,Label_lake, sizethreshold, cpsd_refining, verbose);

        



        % Statistics
        Particle_size = Particle_size*voxel_size;
        all_diameters = Particle_size(binary_phase==1);
        PSD_results(current_phase_todo,1) = min(all_diameters); % Minimum
        PSD_results(current_phase_todo,2) = mean(all_diameters);% Mean
        PSD_results(current_phase_todo,3) = max(all_diameters); % Maximum
        PSD_results(current_phase_todo,4) = std(all_diameters); % Standard deviation
        PSD_results(current_phase_todo,5) = 100*PSD_results(current_phase_todo,4)/PSD_results(current_phase_todo,2); % Relative standard deviation
        % Cumulative and probability density distribution functions
        [PSD(current_phase_todo).psd, ~] = Function_probability_density(all_diameters,[],density_fct_parameters);
        PSD_results(current_phase_todo,6) = PSD(current_phase_todo).psd.x50; % d50
        PSD_results(current_phase_todo,7) = PSD(current_phase_todo).psd.smoothed_x50; % smoothed d50
        PSD_results(current_phase_todo,8) = PSD(current_phase_todo).psd.integral_probability_density_fct;
        PSD_results(current_phase_todo,9) = PSD(current_phase_todo).psd.integral_smoothed_probability_density_fct;

        % % correlation
        results_correlation(current_phase_todo).Particle_diameter_mean_CPSD = PSD_results(current_phase_todo,2);
        results_correlation(current_phase_todo).Particle_diameter_std_CPSD = PSD_results(current_phase_todo,4);
        results_correlation(current_phase_todo).Particle_diameter_relstd_CPSD = PSD_results(current_phase_todo,5);
        %results_correlation(current_phase_todo).Particle_radius_mean_CPSD = PSD_results(current_phase_todo,2)/2; % Example on how to add new correlation
        results_correlation(current_phase_todo).Particle_level_details = Particle_lvl_description(current_phase_todo,1);

        % % Visualization
        results_visualization(current_phase_todo).Particle_diameter_CPSD = Particle_size;
        %results_visualization(current_phase_todo).Particle_radius_CPSD = Particle_size/2; % Example on how to add new visualization

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


end