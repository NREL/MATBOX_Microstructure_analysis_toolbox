function [microstructure3D, phase, outcome] = function_generate_ellipsoid_microstructure(domain_size,phase,Maximum_overlapping,Minimum_particle_volume_conservated,check_contiguity, stoping_conditions, do_verification, save_options)
% Generate ellipsoids-based microstructrures
% Francois Usseglio-Viretta, NREL

refresh_time_for_progression_figure = 5;% s
have_stop_conditions = ~strcmp(stoping_conditions.action,'Ignore (no stoping condition)');

%tic

%% PARAMETERS
simulate_calendering = false; % Useful to set it true to achieve high density
%check_contiguity = true; % Default, true. False useful for very elongated particles


%% DEDUCE MICROSTRUCTURE INFORMATION
%% VOLUME FRACTIONS
voxel_number = prod(domain_size); % Number of voxel
area_xy= domain_size(1)*domain_size(2); % Number of voxel per slice
number_phase = length(phase); % Number of phase
for dir_=1:1:3 % Nornalized domain's length (one value per slice)
    direction(dir_).normalized_coordinates = linspace(0,1,domain_size(dir_));
end

total_vf = 0; % Volume fraction
for current_phase=1:1:number_phase
    x = phase(current_phase).volumefraction.along_3rd_axis(1,:);
    y = phase(current_phase).volumefraction.along_3rd_axis(2,:);
    phase(current_phase).volumefraction.total=trapz(x,y); % volume fraction of current phase
    phase(current_phase).voxel_number = round( phase(current_phase).volumefraction.total*voxel_number ); % number of voxel of current phase
    total_vf = total_vf + phase(current_phase).volumefraction.total; % Sum of all volume fraction
    phase(current_phase).volumefraction.along_3rd_axis_allslices = interp1(x,y,direction(3).normalized_coordinates,'linear'); % Volume fraction per slice
end
if total_vf<=0
    error('microstructure_generation: sum of the volume fractions is <0 while it should be >0, <=1.')
elseif total_vf<1
    generate_last_phase = false; % The remaining volume will not be generated. A phase will be assigned to all the remainig voxels once all the phases have been generated
elseif total_vf==1
    generate_last_phase = true;
else
    error('microstructure_generation: sum of the volume fractions is >1 while it should be >0, <=1.')
end

%% PARTICLE SIZE
for current_phase=1:1:number_phase % Check particle size
    [tmp, ~] = size(phase(current_phase).size_histogram.along_3rd_axis);
    phase(current_phase).number_size_histogram.along_3rd_axis = tmp-1; clear tmp;  % Number of histogram
    for k=1:1:phase(current_phase).number_size_histogram.along_3rd_axis
        if sum(phase(current_phase).size_histogram.along_3rd_axis(k+1,2:end))~=100 % Sum of histogram should be 100
            error('microstructure_generation: sum of the particle size histogram is ~= 100% while it should be 100%.')
        end
    end
    is_odd = mod(phase(current_phase).size_histogram.along_3rd_axis(1,2:end),2);
    if sum(is_odd)~=numel(is_odd) % Some particle diameter is odd.
        idx = find(is_odd==0);
        phase(current_phase).size_histogram.along_3rd_axis(1,idx) = phase(current_phase).size_histogram.along_3rd_axis(1,idx)+1; % odd -> even number
        is_odd(idx) = is_odd(idx)+1; % Convert in a odd integer.
        warning('Even particle diameter (expressed in voxel length) have been converted to odd number.')
    end
end
% Histogram per slice
for current_phase=1:1:number_phase
   x = phase(current_phase).size_histogram.along_3rd_axis(2:end,1);
   y = phase(current_phase).size_histogram.along_3rd_axis(2:end,2:end);
   particlesizeprobability = interp1(x,y,direction(3).normalized_coordinates,'linear');
   if isvector(particlesizeprobability)
       particlesizeprobability=particlesizeprobability';
   end   
   phase(current_phase).size_histogram.along_3rd_axis_allslices = [phase(current_phase).size_histogram.along_3rd_axis(1,2:end); particlesizeprobability]; 
end

%% PARTICLE SIZE ELONGATION
for current_phase=1:1:number_phase % Check particle size
    [tmp1, ~] = size(phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis);
    [tmp2, ~] = size(phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis);
    phase(current_phase).number_elongation_dx_over_dy_histogram.along_3rd_axis = tmp1-1; clear tmp1;  % Number of histogram
    phase(current_phase).number_elongation_dx_over_dz_histogram.along_3rd_axis = tmp2-1; clear tmp2;  % Number of histogram
    for k=1:1:phase(current_phase).number_elongation_dx_over_dy_histogram.along_3rd_axis
        if sum(phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis(k+1,2:end))~=100 % Sum of histogram should be 100
            error('Incorrect input data: sum of the particle elongation histogram is ~= 100% while it should be 100%.')
        end
    end
    for k=1:1:phase(current_phase).number_elongation_dx_over_dz_histogram.along_3rd_axis
        if sum(phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis(k+1,2:end))~=100 % Sum of histogram should be 100
            error('Incorrect input data: sum of the particle elongation histogram is ~= 100% while it should be 100%.')
        end
    end    
end
% Histogram per slice
for current_phase=1:1:number_phase
   x = phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis(2:end,1);
   y = phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis(2:end,2:end);
   particleelongationprobability = interp1(x,y,direction(3).normalized_coordinates,'linear');
   if isvector(particleelongationprobability)
       particleelongationprobability=particleelongationprobability';
   end   
   phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis_allslices = [phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis(1,2:end); particleelongationprobability]; 
   x = phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis(2:end,1);
   y = phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis(2:end,2:end);
   particleelongationprobability = interp1(x,y,direction(3).normalized_coordinates,'linear');
   if isvector(particleelongationprobability)
       particleelongationprobability=particleelongationprobability';
   end
   phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis_allslices = [phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis(1,2:end); particleelongationprobability]; 
end
for current_phase=1:1:number_phase
    phase(current_phase).unique_dx_diameter = unique( phase(current_phase).size_histogram.along_3rd_axis(1,2:end));
    phase(current_phase).unique_dxdy_elongation = unique( phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis(1,2:end));
    phase(current_phase).unique_dxdz_elongation = unique( phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis(1,2:end));
    
    phase(current_phase).unique_dy_diameter = [];
    phase(current_phase).unique_dz_diameter = [];
    for k=1:1:length(phase(current_phase).unique_dx_diameter)
        phase(current_phase).unique_dy_diameter = [phase(current_phase).unique_dy_diameter phase(current_phase).unique_dx_diameter(k)*1./phase(current_phase).unique_dxdy_elongation];
        phase(current_phase).unique_dz_diameter = [phase(current_phase).unique_dz_diameter phase(current_phase).unique_dx_diameter(k)*1./phase(current_phase).unique_dxdz_elongation];
    end
end

%% PARTICLE SIZE ORIENTATION (ROTATION)

tolerance = 1e-3;
for current_phase=1:1:number_phase % Check rotation
    [tmp1, ~] = size(phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis);
    [tmp2, ~] = size(phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis);
    [tmp3, ~] = size(phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis);
    phase(current_phase).number_rotation_x_histogram.along_3rd_axis = tmp1-1; clear tmp1;  % Number of histogram
    phase(current_phase).number_rotation_y_histogram.along_3rd_axis = tmp2-1; clear tmp2;  % Number of histogram
    phase(current_phase).number_rotation_z_histogram.along_3rd_axis = tmp3-1; clear tmp3;  % Number of histogram
    for k=1:1:phase(current_phase).number_rotation_x_histogram.along_3rd_axis
        if abs(sum(phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis(k+1,2:end))-100)>tolerance % Sum of histogram should be 100
            error('Incorrect input data: sum of the particle rotation histogram is ~= 100% while it should be 100%.')
        end
    end
    for k=1:1:phase(current_phase).number_rotation_y_histogram.along_3rd_axis
        if abs(sum(phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis(k+1,2:end))-100)>tolerance % Sum of histogram should be 100
            error('Incorrect input data: sum of the particle rotation histogram is ~= 100% while it should be 100%.')
        end
    end   
    for k=1:1:phase(current_phase).number_rotation_z_histogram.along_3rd_axis
        if abs(sum(phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis(k+1,2:end))-100)>tolerance % Sum of histogram should be 100
            error('Incorrect input data: sum of the particle rotation histogram is ~= 100% while it should be 100%.')
        end
    end        
end

% Histogram per slice
for current_phase=1:1:number_phase
    x = phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis(2:end,1);
    y = phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis(2:end,2:end);
    particlerotationprobability = interp1(x,y,direction(3).normalized_coordinates,'linear');
    if isvector(particlerotationprobability)
        particlerotationprobability=particlerotationprobability';
    end
    phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis_allslices = [phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis(1,2:end); particlerotationprobability];
    
    x = phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis(2:end,1);
    y = phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis(2:end,2:end);
    particlerotationprobability = interp1(x,y,direction(3).normalized_coordinates,'linear');
    if isvector(particlerotationprobability)
        particlerotationprobability=particlerotationprobability';
    end
    phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis_allslices = [phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis(1,2:end); particlerotationprobability];
    
    x = phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis(2:end,1);
    y = phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis(2:end,2:end);
    particlerotationprobability = interp1(x,y,direction(3).normalized_coordinates,'linear');
    if isvector(particlerotationprobability)
        particlerotationprobability=particlerotationprobability';
    end
    phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis_allslices = [phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis(1,2:end); particlerotationprobability];
end
for current_phase=1:1:number_phase
    phase(current_phase).unique_rotation_x = unique( phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis(1,2:end));
    phase(current_phase).unique_rotation_y = unique( phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis(1,2:end));
    phase(current_phase).unique_rotation_z = unique( phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis(1,2:end));
end


%% INITIALIZE VARIABLE USE TO UPDATE VOLUME FRACTION, DIAMETER, PARTICLE ELONGATION, AND PARTICLE ORIENTATION HISTOGRAM PROBABILITY
for current_phase=1:1:number_phase
    phase(current_phase).current_volumefraction.along_3rd_axis_allslices =  phase(current_phase).volumefraction.along_3rd_axis_allslices * 0; % Initialize
    
    unique_diameter = phase(current_phase).unique_dx_diameter;
    unique_dxdy_elongation = phase(current_phase).unique_dxdy_elongation;
    unique_dxdz_elongation = phase(current_phase).unique_dxdz_elongation;    
    phase(current_phase).current_diameter_dx.along_3rd_axis_allslices =  zeros(domain_size(3), length(unique_diameter)); % Initialize
    phase(current_phase).current_elongation_dxdy.along_3rd_axis_allslices =  zeros(domain_size(3), length(unique_dxdy_elongation)); % Initialize
    phase(current_phase).current_elongation_dxdz.along_3rd_axis_allslices =  zeros(domain_size(3), length(unique_dxdz_elongation)); % Initialize
    
    phase(current_phase).current_rotation_x.along_3rd_axis_allslices =  zeros(domain_size(3), length(phase(current_phase).unique_rotation_x)); % Initialize
    phase(current_phase).current_rotation_y.along_3rd_axis_allslices =  zeros(domain_size(3), length(phase(current_phase).unique_rotation_y)); % Initialize
    phase(current_phase).current_rotation_z.along_3rd_axis_allslices =  zeros(domain_size(3), length(phase(current_phase).unique_rotation_z)); % Initialize
end

%%
%%  MICROSTRUCTURE GENERATION
%%

%% PROGRESSION RATE
if stoping_conditions.plot || have_stop_conditions
    time_progression = [0];
    volumefraction_progression = zeros(number_phase,1);  
    particle_progression = zeros(number_phase,1); 
    
    phase_being_generated = [];
    
    time_progression_generationrate = [];
    generationrate_progression = [];
    generationrate_progression_average = [];
    particle_generationrate_progression = [];
    particle_generationrate_progression_average = [];    
    
    time_progression_generationrate_avglastiterations =[];
    generationrate_progression_lastiterations = [];    
    particle_generationrate_progression_lastiterations = [];   
end


if stoping_conditions.plot
    c_ = colororder;
    
    Fig_global_progression = figure; % Create figure
    Fig_global_progression.Name= 'Algorithm generation progression and rate'; % Figure name
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_global_progression,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)/2]); % Full screen figure
    Fig_global_progression.Color='white'; % Background colour
    
    sub_axes_progression_1=subplot(1,3,1,'Parent',Fig_global_progression);
    hold(sub_axes_progression_1,'on'); % Active subplot
    h_title=title ('Volume fraction');

    for current_phase = 1:1:number_phase
        plot([0,1],[phase(current_phase).volumefraction.total phase(current_phase).volumefraction.total],'LineWidth',2,'LineStyle','--','Color',c_(current_phase,:),'DisplayName',[phase(current_phase).name ' (inputs)']);
    end
    ylim(sub_axes_progression_1,[0 +Inf])
    xlabel(sub_axes_progression_1,'Wallclock time (s)');
    ylabel(sub_axes_progression_1,'Volume fraction');
    grid(sub_axes_progression_1,'on');
    legend(sub_axes_progression_1,'Location','best');
    set(sub_axes_progression_1,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
    set(h_title,'FontSize',14)
    hold(sub_axes_progression_1,'off');
    
    sub_axes_progression_2=subplot(1,3,2,'Parent',Fig_global_progression,'XScale', 'linear', 'YScale', 'log');
    hold(sub_axes_progression_2,'on'); % Active subplo
    h_title=title ('Volume fraction generation rate');
    xlabel(sub_axes_progression_2,'Wallclock time (s)');
    ylabel(sub_axes_progression_2,'Volume fraction.s^{-1}');
    grid(sub_axes_progression_2,'on');
    legend(sub_axes_progression_2,'Location','best');
    set(sub_axes_progression_2,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
    set(h_title,'FontSize',14)    
    hold(sub_axes_progression_2,'off');
    
    sub_axes_progression_3=subplot(1,3,3,'Parent',Fig_global_progression);
    hold(sub_axes_progression_3,'on'); % Active subplo
    h_title=title ('Particle generation rate');
    xlabel(sub_axes_progression_3,'Wallclock time (s)');
    ylabel(sub_axes_progression_3,'Particle.s^{-1}');
    grid(sub_axes_progression_3,'on');
    legend(sub_axes_progression_3,'Location','best');
    set(sub_axes_progression_3,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
    set(h_title,'FontSize',14)    
    hold(sub_axes_progression_3,'off'); 
    
    sgtitle(Fig_global_progression,'Algorithm generation progression and rate','FontWeight','bold','FontSize',16,'FontName','Times new roman');
end

%%
start_generation_time = tic; % Start clockwatch
start_figure_time = tic;

% Phase are generated sequentially.
% Order your phase so that the phase i+1 contains less and smaller particles compared with the phase i
% Indeed, it is more difficult to find space for a large particles when small particles already fill the domain's space.

% Check each phase has the same number of pass
number_pass = length(phase(1).fillratio);
for current_phase=2:1:number_phase
    tmp = length(phase(current_phase).fillratio);
    if number_pass~=tmp
        error('Incorrect input data: each phase must have the same number of pass. length(phase(current_phase).fillratio) should be the same.')
    end
end

% Initialize
unassigned_id = 0; 
microstructure3D.phase = zeros(domain_size) + unassigned_id; % Phase

microstructure3D.particle_id = uint32(zeros(domain_size)); % Particle numerotation
microstructure3D.particle_volume_numbervoxel = uint32(zeros(domain_size)); % Particle volume expressed in number of voxels
microstructure3D.particle_volume_ellipsoid = zeros(domain_size); % Particle ellipsoid volume
microstructure3D.particle_equivalentspherediameter = zeros(domain_size); % Particle equivalent sphere diameter
for current_phase=1:1:number_phase
    microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy = zeros(domain_size); % Ellispoid elongation dx/dy
    microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz = zeros(domain_size); % Ellispoid elongation dx/dz
    microstructure3D.phase_(current_phase).particle_length_x = zeros(domain_size); % Ellispoid lenght along axe 1
    microstructure3D.phase_(current_phase).particle_length_y = zeros(domain_size); % Ellispoid lenght along axe 2
    microstructure3D.phase_(current_phase).particle_length_z = zeros(domain_size); % Ellispoid lenght along axe 3
    microstructure3D.phase_(current_phase).particle_rotation_x = zeros(domain_size); % Ellispoid rotation along axe 1
    microstructure3D.phase_(current_phase).particle_rotation_y = zeros(domain_size); % Ellispoid rotation along axe 2
    microstructure3D.phase_(current_phase).particle_rotation_z = zeros(domain_size); % Ellispoid rotation along axe 3
end
particle_id = 1; % Particle current numero

% Index and coordinates of all unassigned voxels
Idx_all_voxels = 1:1:voxel_number;
Idx_unassigned_voxels_pickfrom = Idx_all_voxels;
Idx_assigned_voxels=zeros(1,voxel_number);
all_Idx_assigned_voxels=[];

Idx_unassigned_voxels = Idx_all_voxels;
[I,J,K] = ind2sub(domain_size,Idx_all_voxels);
Coordinates_all_voxels = [I' J' K'];
number_unassigned_voxels = voxel_number;

% Particle info
% row i: particle id = i
preallocating_particle = 1e6;
particle_info = struct('Center_xyz',zeros(preallocating_particle,3),...
                       'Number_voxel',zeros(preallocating_particle,1),...
                       'Phase_code',zeros(preallocating_particle,1),...
                       'Phase_number',zeros(preallocating_particle,1),...
                       'Extremities',zeros(preallocating_particle,6));                   
particle_index(1:preallocating_particle) = struct('Index',0);
                
% loop pass (each new pass phase_volume_target is larger: one value per pass)
%   loop phase
% Then, distance map to fill the remaining, with competition between phase for zone attibuted to multiple phase
update_waitbar_step=0.05;
phase_volume_filled = zeros(number_phase,1); % Number of voxel assigned to each phase
phase_number_particle = zeros(number_phase,1); % Number of particle assigned to each phase
number_iteration = 0;
complete_stop = false;
outcome = 'Success !'; % Default
first_continue = true;
for current_pass=1:1:number_pass
    for current_phase=1:1:number_phase
        stop_condition_not_reached = true;
        if complete_stop
            break
        end
        % Progession bar
        unique_diameter = phase(current_phase).unique_dx_diameter;
        unique_dxdy_elongation = phase(current_phase).unique_dxdy_elongation;
        unique_dxdz_elongation = phase(current_phase).unique_dxdz_elongation;
        histogram_currentphase_z_diameter = phase(current_phase).size_histogram.along_3rd_axis_allslices(1,:); % Diameter histogram (x) of the current phase
        histogram_currentphase_z_elongation_dx_over_dy = phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis_allslices(1,:); % elongation dx/dy histogram (x) of the current phase
        histogram_currentphase_z_elongation_dx_over_dz = phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis_allslices(1,:); % elongation dx/dz histogram (x) of the current phase
        
        unique_rotation_x = phase(current_phase).unique_rotation_x;
        unique_rotation_y = phase(current_phase).unique_rotation_y;
        unique_rotation_z = phase(current_phase).unique_rotation_z;
        histogram_currentphase_z_rotation_x = phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis_allslices(1,:); % Diameter histogram (x) of the current phase
        histogram_currentphase_z_rotation_y = phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis_allslices(1,:); % elongation dx/dy histogram (x) of the current phase
        histogram_currentphase_z_rotation_z = phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis_allslices(1,:); % elongation dx/dz histogram (x) of the current phase
    
        phase_volume_target = phase(current_phase).voxel_number * phase(current_phase).fillratio(current_pass); % Number of voxel to assign for this phase before moving forward
        
        new_particle_hasbeen_generated = false;

        % Wait bar
        generation_progress = phase_volume_filled(current_phase,1)/phase_volume_target; % 0-1
        previous_generation_progress = generation_progress;
        h_waitbar = waitbar(generation_progress,['Pass:' num2str(current_pass) '/' num2str(number_pass), ' Phase ' num2str(current_phase) '/' num2str(number_phase) ' is being generated']);
        while phase_volume_filled(current_phase,1) <= phase_volume_target && stop_condition_not_reached
            number_iteration = number_iteration+1;
            
            % Progression/rate
            time_since_start = toc(start_generation_time);
            if time_since_start>stoping_conditions.maxtime
                complete_stop = true;
                outcome = ['Algorithm stalls at pass ' num2str(current_pass) ', for ' phase(current_phase).name];
                break
            end            
            
            if stoping_conditions.plot || have_stop_conditions
                % Update array for figure
                time_progression = [time_progression time_since_start];
                volumefraction_progression = [volumefraction_progression phase_volume_filled(:,1)/voxel_number];
                particle_progression = [particle_progression  phase_number_particle(:,1)];
                
                if number_iteration>1
                    phase_being_generated = [phase_being_generated current_phase];
                    time_progression_generationrate = [time_progression_generationrate time_progression(end-1)+((time_progression(end)-time_progression(end-1))/2)];
                    generationrate_progression = [generationrate_progression (volumefraction_progression(:,end)-volumefraction_progression(:,end-1)) ./ (time_progression(end)-time_progression(end-1)) ];
                    generationrate_progression_average = [generationrate_progression_average mean(generationrate_progression,2)];
                    particle_generationrate_progression = [particle_generationrate_progression (particle_progression(:,end)-particle_progression(:,end-1)) ./ (time_progression(end)-time_progression(end-1)) ];
                    particle_generationrate_progression_average = [particle_generationrate_progression_average mean(particle_generationrate_progression,2)];                    
                    
                end
                
                if number_iteration>stoping_conditions.average_on_last_n_iterations+1
                    time_progression_generationrate_avglastiterations = [time_progression_generationrate_avglastiterations time_since_start];
                    
                    % The average on a phase generation rate must be done only for iterations that correspond to this phase
                    pha = phase_being_generated(end-stoping_conditions.average_on_last_n_iterations:end);
                    val1 = generationrate_progression(:,end-stoping_conditions.average_on_last_n_iterations:end);
                    val2 = particle_generationrate_progression(:,end-stoping_conditions.average_on_last_n_iterations:end);
                    val1a = NaN(number_phase,1); val2a = NaN(number_phase,1);
                    for kk=1:1:number_phase
                        idx_p = find(pha==kk);
                        if ~isempty(idx_p) && length(idx_p)>=stoping_conditions.average_on_last_n_iterations
                            val1a(kk,1) = mean(val1(kk,idx_p),2);
                            val2a(kk,1) = mean(val2(kk,idx_p),2);
                        end
                    end
                    generationrate_progression_lastiterations = [generationrate_progression_lastiterations val1a];
                    particle_generationrate_progression_lastiterations = [particle_generationrate_progression_lastiterations val2a];
                                
                    % Check stoping condition
                    max_vf_rate = max(generationrate_progression_lastiterations(:,end));
                    max_particle_rate = max(particle_generationrate_progression_lastiterations(:,end));
                    if max_vf_rate<stoping_conditions.vfrate_threshold
                        vfrate_stoping_condition = true;
                    else
                        vfrate_stoping_condition = false;
                    end
                    if max_particle_rate<stoping_conditions.particlerate_threshold
                        particlerate_stoping_condition = true;
                    else
                        particlerate_stoping_condition = false;
                    end                    
                    if (strcmp(stoping_conditions.andor,'or') && (vfrate_stoping_condition || particlerate_stoping_condition)) || (strcmp(stoping_conditions.andor,'and') && (vfrate_stoping_condition && particlerate_stoping_condition))
                       % Reach stoping conditions
                       if strcmp(stoping_conditions.action,'Move to next phase or generation pass (and stop if last iteration)')
                           %disp 'Next'
                           if first_continue
                               outcome = ['Algorithm stalls at pass ' num2str(current_pass) ', for ' phase(current_phase).name];
                               first_continue= false;
                           else
                               outcome = [outcome '. Then, pass ' num2str(current_pass) ', for ' phase(current_phase).name];
                           end
                           stop_condition_not_reached = false; % Exit while loop and go to next phase
                       elseif strcmp(stoping_conditions.action,'Stop')
                           %disp 'Stop'
                           complete_stop = true;
                           outcome = ['Algorithm stalls at pass ' num2str(current_pass) ', for ' phase(current_phase).name];
                           break
                       end
                    end
                end
                
                % Update figure
                if stoping_conditions.plot
                    time_since_last_plot = toc(start_figure_time);
                    if time_since_last_plot>=refresh_time_for_progression_figure
                        start_figure_time = tic; % reset timer
                        
                        cla(sub_axes_progression_1); % Clear axes
                        hold(sub_axes_progression_1,'on'); % Active subplot
                        for current_phase_figure = 1:1:number_phase
                            plot(sub_axes_progression_1,[0,time_since_start],[phase(current_phase_figure).volumefraction.total phase(current_phase_figure).volumefraction.total],'LineWidth',2,'LineStyle','--','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (inputs)']);
                            plot(sub_axes_progression_1,time_progression,volumefraction_progression(current_phase_figure,:),'LineWidth',2,'LineStyle','-','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (generated)']);
                        end
                        ylim(sub_axes_progression_1,[0 +Inf])
                        xlim(sub_axes_progression_1,[0 time_since_start])
                        legend(sub_axes_progression_1,'Location','best');
                        hold(sub_axes_progression_1,'off');
                        
                        cla(sub_axes_progression_2); % Clear axes
                        hold(sub_axes_progression_2,'on'); % Active subplot
                        for current_phase_figure = 1:1:number_phase
                            if ~isempty(time_progression_generationrate)
                                plot(sub_axes_progression_2,time_progression_generationrate,generationrate_progression(current_phase_figure,:),'LineWidth',1,'LineStyle','-','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (instantaneous)']);
                                plot(sub_axes_progression_2,time_progression_generationrate,generationrate_progression_average(current_phase_figure,:),'LineWidth',2,'LineStyle','--','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (global average)']);
                            end
                            if ~isempty(time_progression_generationrate_avglastiterations)
                                plot(sub_axes_progression_2,time_progression_generationrate_avglastiterations,generationrate_progression_lastiterations(current_phase_figure,:),'LineWidth',2,'LineStyle',':','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (average of last ' num2str(stoping_conditions.average_on_last_n_iterations) ' iterations)']);
                            end
                        end
                        xlim(sub_axes_progression_2,[0 time_since_start])
                        legend(sub_axes_progression_2,'Location','best');
                        hold(sub_axes_progression_2,'off');
                        
                        cla(sub_axes_progression_3); % Clear axes
                        hold(sub_axes_progression_3,'on'); % Active subplot
                        for current_phase_figure = 1:1:number_phase
                            if ~isempty(time_progression_generationrate)
                                %plot(sub_axes_progression_3,time_progression_generationrate,particle_generationrate_progression(current_phase_figure,:),'LineWidth',1,'LineStyle','-','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (instantaneous)']);
                                plot(sub_axes_progression_3,time_progression_generationrate,particle_generationrate_progression_average(current_phase_figure,:),'LineWidth',2,'LineStyle','--','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (global average)']);
                            end
                            if ~isempty(time_progression_generationrate_avglastiterations)
                                plot(sub_axes_progression_3,time_progression_generationrate_avglastiterations,particle_generationrate_progression_lastiterations(current_phase_figure,:),'LineWidth',2,'LineStyle',':','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (average of last ' num2str(stoping_conditions.average_on_last_n_iterations) ' iterations)']);
                            end
                        end
                        xlim(sub_axes_progression_3,[0 time_since_start])
                        legend(sub_axes_progression_3,'Location','best');
                        hold(sub_axes_progression_3,'off');
                        pause(0.01); % Force refresh of the graph
                    end
                end
            end

            generation_progress = phase_volume_filled(current_phase,1)/phase_volume_target; % 0-1
            if generation_progress-previous_generation_progress > update_waitbar_step
                waitbar(generation_progress,h_waitbar);
                previous_generation_progress=generation_progress;
            end
            
            if new_particle_hasbeen_generated
                particle_id = particle_id+1; % Update
                % Update volume fractions
                phase(current_phase).current_volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max) = sum(sum(microstructure3D.phase(:,:,z_update_min:z_update_max)==phase(current_phase).code))/area_xy;
                % Update particle size histogram
                for k_diameter=1:1:length(unique_diameter)
                    current_section = zeros(z_update_max-z_update_min+1,1);
                    current_section(:) = sum(sum(microstructure3D.phase_(current_phase).particle_length_x(:,:,z_update_min:z_update_max)==unique_diameter(k_diameter)));
                    target_section = area_xy * phase(current_phase).volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max)';
                    phase(current_phase).current_diameter_dx.along_3rd_axis_allslices(z_update_min:z_update_max,k_diameter) = 100*current_section./target_section;
                end
                % Update particle elongation histogram
                for k_elongation_dxdy=1:1:length(unique_dxdy_elongation)
                    current_section = zeros(z_update_max-z_update_min+1,1);
                    current_section(:) = sum(sum(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy(:,:,z_update_min:z_update_max)==unique_dxdy_elongation(k_elongation_dxdy)));
                    target_section = area_xy * phase(current_phase).volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max)';
                    phase(current_phase).current_elongation_dxdy.along_3rd_axis_allslices(z_update_min:z_update_max,k_elongation_dxdy) = 100*current_section./target_section;
                end
                for k_elongation_dxdz=1:1:length(unique_dxdz_elongation)
                    current_section = zeros(z_update_max-z_update_min+1,1);
                    current_section(:) = sum(sum(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz(:,:,z_update_min:z_update_max)==unique_dxdz_elongation(k_elongation_dxdz)));
                    target_section = area_xy * phase(current_phase).volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max)';
                    phase(current_phase).current_elongation_dxdz.along_3rd_axis_allslices(z_update_min:z_update_max,k_elongation_dxdz) = 100*current_section./target_section;
                end
                
                % Update particle orientation histogram
                for k_otientation_x=1:1:length(unique_rotation_x)
                    current_section = zeros(z_update_max-z_update_min+1,1);
                    current_section(:) = sum(sum(microstructure3D.phase_(current_phase).particle_rotation_x(:,:,z_update_min:z_update_max)==unique_rotation_x(k_otientation_x)));
                    target_section = area_xy * phase(current_phase).volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max)';
                    phase(current_phase).current_rotation_x.along_3rd_axis_allslices(z_update_min:z_update_max,k_otientation_x) = 100*current_section./target_section;
                end
                for k_otientation_y=1:1:length(unique_rotation_y)
                    current_section = zeros(z_update_max-z_update_min+1,1);
                    current_section(:) = sum(sum(microstructure3D.phase_(current_phase).particle_rotation_y(:,:,z_update_min:z_update_max)==unique_rotation_y(k_otientation_y)));
                    target_section = area_xy * phase(current_phase).volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max)';
                    phase(current_phase).current_rotation_y.along_3rd_axis_allslices(z_update_min:z_update_max,k_otientation_y) = 100*current_section./target_section;
                end      
                for k_otientation_z=1:1:length(unique_rotation_z)
                    current_section = zeros(z_update_max-z_update_min+1,1);
                    current_section(:) = sum(sum(microstructure3D.phase_(current_phase).particle_rotation_z(:,:,z_update_min:z_update_max)==unique_rotation_z(k_otientation_z)));
                    target_section = area_xy * phase(current_phase).volumefraction.along_3rd_axis_allslices(z_update_min:z_update_max)';
                    phase(current_phase).current_rotation_z.along_3rd_axis_allslices(z_update_min:z_update_max,k_otientation_z) = 100*current_section./target_section;
                end                    
                                
                % Update unassigned voxel index
                % Note that you may use a simpler explicit method, such as I=find(microstructure3D.phase==unassigned_id), but it is much more slower.
                [all_Idx_assigned_voxels, Idx_assigned_voxels, Idx_unassigned_voxels_pickfrom] = update_unassigned_index(all_Idx_assigned_voxels, linear_indices, Idx_assigned_voxels, Idx_all_voxels, number_unassigned_voxels, Idx_unassigned_voxels_pickfrom);
                new_particle_hasbeen_generated = false; % Reset
            end
            
            % Particle center position
            rn_pos = randi([1 voxel_number],1); % Pick a random number among all the voxels
            idx_center = Idx_unassigned_voxels_pickfrom(rn_pos); % Select an indice among all the unassigned voxels
            x_center = Coordinates_all_voxels(idx_center,1);
            y_center = Coordinates_all_voxels(idx_center,2);
            z_center = Coordinates_all_voxels(idx_center,3);
            
            remaining_volume_fraction_currentphase_z = phase(current_phase).volumefraction.along_3rd_axis_allslices(z_center) - phase(current_phase).current_volumefraction.along_3rd_axis_allslices(z_center); % volume fraction of the current phase at z_center
            
            normalized_term = 0;
            for kk=1:1:number_phase
                normalized_term = normalized_term + phase(kk).volumefraction.along_3rd_axis_allslices(z_center) - phase(kk).current_volumefraction.along_3rd_axis_allslices(z_center);
            end
            
            rn = rand; % Pick a random number from 0 to 1
            if rn < remaining_volume_fraction_currentphase_z/normalized_term % Check volume fraction
                % Particle diameter
                new_probability = update_probability(phase(current_phase).size_histogram.along_3rd_axis_allslices(z_center+1,:), phase(current_phase).current_diameter_dx.along_3rd_axis_allslices(z_center,:));
                idx = choose_randomly_from_histogramprobability( new_probability ); % Diameter histogram (y) of the current phase at z_center
                current_diameter = histogram_currentphase_z_diameter(idx);
                if current_diameter==1 % One-voxel particle
                    z_update_min = z_center; z_update_max = z_center; % Specify bounds where to update probabilities
                    number_voxel_assigned = 1; % Number of new voxel for the current phase
                    number_voxel_idealshape = 1;
                    microstructure3D.phase(x_center, y_center, z_center) = phase(current_phase).code; % Insert voxel within the domain.
                    % Indices of newly assigned voxels
                    linear_indices = rn_pos;
                    % 3D arrays assignment
                    microstructure3D = ThreeD_array_assignment(microstructure3D, current_phase, linear_indices, particle_id, number_voxel_assigned, number_voxel_idealshape, 1, 1, 1, 0, 0, 0, x_center, y_center, z_center);
                    % Update
                    new_particle_hasbeen_generated = true;
                    particle_xyz_min_max = [x_center y_center z_center x_center y_center z_center];
                    
                else % n-voxel particle
                    % Particle elongation
                    new_probability = update_probability(phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis_allslices(z_center+1,:), phase(current_phase).current_elongation_dxdy.along_3rd_axis_allslices(z_center,:));
                    idx = choose_randomly_from_histogramprobability( new_probability ); % Elongation dx/dy histogram (y) of the current phase at z_center
                    current_elongation_dx_over_dy = histogram_currentphase_z_elongation_dx_over_dy(idx); % Particle elongation
                    new_probability = update_probability(phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis_allslices(z_center+1,:), phase(current_phase).current_elongation_dxdz.along_3rd_axis_allslices(z_center,:));
                    idx = choose_randomly_from_histogramprobability( new_probability ); % Elongation dx/dz histogram (y) of the current phase at z_center
                    current_elongation_dx_over_dz = histogram_currentphase_z_elongation_dx_over_dz(idx); % Particle elongation
                    dx = current_diameter; % Particle lenght along axe 1
                    dy = dx * (1/current_elongation_dx_over_dy); % Particle lenght along axe 2
                    dz = dx * (1/current_elongation_dx_over_dz); % Particle lenght along axe 3
                    binary_ellipsoid = create_ellipsoid(dx,dy,dz); % Create ellipsoid
                    
                    % Particle orientation
                    new_probability = update_probability(phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis_allslices(z_center+1,:), phase(current_phase).current_rotation_x.along_3rd_axis_allslices(z_center,:));
                    idx = choose_randomly_from_histogramprobability( new_probability ); % Rotation x histogram (y) of the current phase at z_center
                    angle_x_deg = histogram_currentphase_z_rotation_x(idx); % Particle rotation                    
                    new_probability = update_probability(phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis_allslices(z_center+1,:), phase(current_phase).current_rotation_y.along_3rd_axis_allslices(z_center,:));
                    idx = choose_randomly_from_histogramprobability( new_probability ); % Rotation y histogram (y) of the current phase at z_center
                    angle_y_deg = histogram_currentphase_z_rotation_y(idx); % Particle rotation    
                    new_probability = update_probability(phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis_allslices(z_center+1,:), phase(current_phase).current_rotation_z.along_3rd_axis_allslices(z_center,:));
                    idx = choose_randomly_from_histogramprobability( new_probability ); % Rotation z histogram (y) of the current phase at z_center
                    angle_z_deg = histogram_currentphase_z_rotation_z(idx); % Particle rotation                        
                    if dx~=dy || dx~=dz || dy~=dz % Only for ellipsoid
                        if angle_x_deg~=0 || angle_y_deg~=0 || angle_z_deg~=0 % Only if one or more rotation
                            binary_ellipsoid = rotate_domain(binary_ellipsoid,angle_x_deg, angle_y_deg, angle_z_deg); % Apply rotation matrix
                            binary_ellipsoid(binary_ellipsoid>=0.5)=1;
                            binary_ellipsoid(binary_ellipsoid~=1)=0;
                        end
                    end
                    
                    size_ellipsoid = size(binary_ellipsoid); % Dimension of the ellipsoid
                    % Subdomain and domain bounds that contain the ellipsoid
                    [x_sub_min, x_sub_max, x_min, x_max] = bounds_subdomain_domain(size_ellipsoid(1),x_center,1,domain_size(1));
                    [y_sub_min, y_sub_max, y_min, y_max] = bounds_subdomain_domain(size_ellipsoid(2),y_center,1,domain_size(2));
                    [z_sub_min, z_sub_max, z_min, z_max] = bounds_subdomain_domain(size_ellipsoid(3),z_center,1,domain_size(3));
                    subdomain = microstructure3D.phase(x_min:x_max, y_min:y_max, z_min:z_max); % Subdomain
                    subdomain_ellipsoid = binary_ellipsoid(x_sub_min:x_sub_max, y_sub_min:y_sub_max, z_sub_min:z_sub_max); % Ellipsoid subdomain
                    z_update_min = z_min; z_update_max = z_max; % Specify bounds where to update probabilities, in case particle will be generated
                    
                    % Insertion of the ellipsoid differ depending on the subdomain is empty or contains already other particles
                    unique_subdomain = unique(subdomain);
                    if length(unique_subdomain)==1 && unique_subdomain==unassigned_id % Subdomain is empty
                        if simulate_calendering==true && particle_id>1
                            new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                            continue_whileloop_and_go_to_next_randcenter = true;
                            continue
                        end
                        number_voxel_assigned = sum(sum(sum(subdomain_ellipsoid == 1))); % Number of new voxel for the current phase
                        number_voxel_idealshape = number_voxel_assigned;
                        microstructure3D.phase(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain_ellipsoid * phase(current_phase).code; % Insert ellipsoid within the domain.
                        % Indices of newly assigned voxels
                        [linear_indices] = Index_from_subdomain_to_domain(subdomain_ellipsoid, [], domain_size, [x_min y_min z_min]);
                        % 3D arrays assignment
                        microstructure3D = ThreeD_array_assignment(microstructure3D, current_phase, linear_indices, particle_id, number_voxel_assigned, number_voxel_idealshape, dx, dy, dz, angle_x_deg, angle_y_deg, angle_z_deg, x_center, y_center, z_center);
                        % Update
                        particle_xyz_min_max = [x_min y_min z_min x_max y_max z_max];
                        new_particle_hasbeen_generated = true;
                        
                    else % Subdomain contains other particles
                        subdomain_binary = subdomain;
                        subdomain_binary(subdomain_binary~=unassigned_id)=1;
                        union_subdomain_ellipsoid = subdomain_binary + subdomain_ellipsoid;
                        idx_interpenetration = find(union_subdomain_ellipsoid==2);
                        if isempty(idx_interpenetration) % There is no overlapping
                            if simulate_calendering==true && particle_id>1
                                new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                                continue_whileloop_and_go_to_next_randcenter = true;
                                continue
                            end
                            number_voxel_assigned = sum(sum(sum(subdomain_ellipsoid == 1))); % Number of new voxel for the current phase
                            number_voxel_idealshape = number_voxel_assigned;
                            microstructure3D.phase(x_min:x_max, y_min:y_max, z_min:z_max) = subdomain + (subdomain_ellipsoid * phase(current_phase).code); % Insert ellipsoid within the domain.
                            % Indices of newly assigned voxels
                            [linear_indices] = Index_from_subdomain_to_domain(subdomain_ellipsoid, [], domain_size, [x_min y_min z_min]);
                            % 3D arrays assignment
                            microstructure3D = ThreeD_array_assignment(microstructure3D, current_phase, linear_indices, particle_id, number_voxel_assigned, number_voxel_idealshape, dx, dy, dz, angle_x_deg, angle_y_deg, angle_z_deg, x_center, y_center, z_center);
                            % Update
                            particle_xyz_min_max = [x_min y_min z_min x_max y_max z_max];
                            new_particle_hasbeen_generated = true;
                            
                        else % New and existing particles are overlapping
                            continue_whileloop_and_go_to_next_randcenter = false; % Initialization
                            subdomain_id = microstructure3D.particle_id(x_min:x_max, y_min:y_max, z_min:z_max);
                            particles_id_subdomain = unique(subdomain_id);
                            particles_id_subdomain(particles_id_subdomain==0)=[];
                            n_otherparticles = length(particles_id_subdomain);
                            
                            % % First let's check these particles can overlap
                            for k_other = 1:1:n_otherparticles % Loop over all existing particles
                                other_id = particles_id_subdomain(k_other); % Get code of the other particle
                                other_phasenumber =  particle_info.Phase_number(other_id);
                                current_maximum_interpenetration = Maximum_overlapping(current_phase, other_phasenumber);
                                if current_maximum_interpenetration<=0 % No need to continue
                                    new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                                    continue_whileloop_and_go_to_next_randcenter = true;
                                    break
                                end
                            end
                            if continue_whileloop_and_go_to_next_randcenter
                                continue
                            end
                            
                            % % Then, let's check these particles do not include each other
                            [linear_indices] = Index_from_subdomain_to_domain(subdomain_ellipsoid, [], domain_size, [x_min y_min z_min]);
                            number_voxel_idealshape = length(linear_indices);
                            microstructure_tmp = zeros(domain_size);
                            microstructure_tmp(linear_indices)=1;
                            for k_other = 1:1:n_otherparticles % Loop over all existing particles
                                other_id = particles_id_subdomain(k_other); % Get code of the other particle
                                if sum(microstructure_tmp(particle_index(other_id).Index))==length(particle_index(other_id).Index)
                                    new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                                    continue_whileloop_and_go_to_next_randcenter = true;
                                    break
                                end
                            end
                            if continue_whileloop_and_go_to_next_randcenter
                                continue
                            end
                            
                            % % Then, let's check these particles overlap only to a certain extent
                            current_center = [x_center y_center z_center]; % Ellipsoid center
                            idx_all_overlapping_subdomain = [];
                            loss_voxels_fromidealshape = 0; % Initialize
                            for k_other = 1:1:n_otherparticles % Loop over all existing particles
                                other_id = particles_id_subdomain(k_other); % Get code of the other particle
                                other_phasenumber =  particle_info.Phase_number(other_id);                                
                                other_code =  particle_info.Phase_code(other_id);
                                other_center = particle_info.Center_xyz(other_id,:); % Get coordinate center of the other particle
                                tmp = find(subdomain_id(idx_interpenetration)==other_id); % Index of the overlapping region for the current particle
                                if length(tmp)>=1
                                    idx_overlapping_ellipsoid_other = idx_interpenetration(tmp); % Linear index of the overlapping between the ellispode and the current other particle, subdomain
                                    [linear_indices_overlapping_ellipsoid_other] = Index_from_subdomain_to_domain(subdomain_ellipsoid, idx_overlapping_ellipsoid_other, domain_size, [x_min y_min z_min]); % Same index, but expressed in the full domain
                                    [I_overlapping,J_overlapping,K_overlapping] = ind2sub(domain_size,linear_indices_overlapping_ellipsoid_other); % Coordinates of the overlaping region expressed in the full domain
                                    
                                    distance_overlapping_to_center = vecnorm( [I_overlapping,J_overlapping,K_overlapping]-current_center, 2, 2 ); % Euclidean distance from each voxel of the overlapping region with the particle center
                                    distance_overlapping_to_othercenter = vecnorm( [I_overlapping,J_overlapping,K_overlapping]-other_center, 2, 2 ); % Euclidean distance from each voxel of the overlapping region with the other particle center
                                    d_min_overlapping_to_centre = min(distance_overlapping_to_center); % Min and max of the previous results
                                    d_max_overlapping_to_centre = max(distance_overlapping_to_center);
                                    ratio_overlapping = (d_max_overlapping_to_centre-d_min_overlapping_to_centre)/d_max_overlapping_to_centre;
                                    d_min_overlapping_to_othercentre = min(distance_overlapping_to_othercenter); % Min and max of the previous results
                                    d_max_overlapping_to_othercentre = max(distance_overlapping_to_othercenter);
                                    ratio_overlapping_other = (d_max_overlapping_to_othercentre-d_min_overlapping_to_othercentre)/d_max_overlapping_to_centre;
                                    
                                    current_maximum_interpenetration = Maximum_overlapping(current_phase, other_phasenumber);
                                    if max(ratio_overlapping, ratio_overlapping_other)<=current_maximum_interpenetration % Check overlapping
                                        
                                        % In case the overlapping is allowed, each particle will loose this number of voxels (which them impact its original shape)
                                        % This loss is cumulative with each new overlapping.
                                        estimation_numbervoxel_foreach_particle = length(idx_overlapping_ellipsoid_other)/2; % True for sphere-sphere, estimation for ellipsoid-ellipsoid
                                        % Estimation of the other particle volume ratio if the overlapping will be allowed
                                        ratio_other = (particle_info.Number_voxel_currently_assigned(other_id) - estimation_numbervoxel_foreach_particle)/particle_info.Number_voxel_assigned_whencreated(other_id);
                                        % Estimation of the current particle volume ratio if the overlapping will be allowed
                                        loss_voxels_fromidealshape = loss_voxels_fromidealshape + estimation_numbervoxel_foreach_particle;
                                        ratio_ = (number_voxel_idealshape-loss_voxels_fromidealshape)/number_voxel_idealshape;
                                        % Check ratio of original volume is still high enough
                                        if ratio_other<=Minimum_particle_volume_conservated(other_phasenumber) || ratio_<=Minimum_particle_volume_conservated(current_phase)
                                            new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                                            continue_whileloop_and_go_to_next_randcenter = true;
                                            break
                                        else
                                            idx_all_overlapping_subdomain = [idx_all_overlapping_subdomain; idx_overlapping_ellipsoid_other];
                                        end
                                    else % Overlapping is too important
                                        new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                                        continue_whileloop_and_go_to_next_randcenter = true;
                                        break
                                    end
                                end
                            end
                            if continue_whileloop_and_go_to_next_randcenter
                                continue
                            end
                            
                            % % Then, let's distribute the overlapped region between the current particle and the others
                            subdomain_union_minus_intersection = abs(subdomain_ellipsoid - subdomain_binary); % (particle union others) - (particle intersection others)
                            [Euclidean_distance_map,Idx_distance_map] = bwdist(subdomain_union_minus_intersection);
                            index_nearest_particle = Idx_distance_map(idx_all_overlapping_subdomain);
                            id_overlapping = zeros(size(subdomain_ellipsoid));
                            id_overlapping(idx_all_overlapping_subdomain) = double(subdomain_id(index_nearest_particle)) + subdomain_ellipsoid(index_nearest_particle)*particle_id;
                            new_subdomain_id = double(subdomain_id) + subdomain_ellipsoid*particle_id; % Wrong overlapping
                            new_subdomain_id(idx_all_overlapping_subdomain) = id_overlapping(idx_all_overlapping_subdomain); % Correct overlapping
                            
                            % % Then, check particles have not been cut in pieces (discontinued particles)
                            all_particles_subdomain=[];
                            for k_particle = 1:1:n_otherparticles+1 % Loop over pre-existing particles and new particle
                                if k_particle==n_otherparticles+1
                                    id_ = particle_id; % New particle
                                else
                                    id_ = particles_id_subdomain(k_particle); % previous particle
                                end
                                binary_particle = zeros(size(new_subdomain_id)); % Initialize
                                all_particles_subdomain(k_particle).id = id_;
                                all_particles_subdomain(k_particle).idx = find(new_subdomain_id==id_);
                                binary_particle(all_particles_subdomain(k_particle).idx)=1;
                                if check_contiguity
                                    L = bwlabeln(binary_particle,6); % Check face-to-face connectivity
                                    if length(unique(L))>2 % It should be 2: background + one unique particle cluster
                                        new_particle_hasbeen_generated = false;  % Not necessary. For clarification sake.
                                        continue_whileloop_and_go_to_next_randcenter = true;
                                        break
                                    end
                                end
                                
                            end
                            if continue_whileloop_and_go_to_next_randcenter
                                continue
                            end
                            
                            % % Then, assign current ellipsoid and update other particles
                            for k_particle = 1:1:n_otherparticles+1
                                if k_particle==n_otherparticles+1 % New particle
                                    number_voxel_assigned = length(all_particles_subdomain(k_particle).idx);
                                    % Indices of newly assigned voxels
                                    [linear_indices] = Index_from_subdomain_to_domain(subdomain_ellipsoid, all_particles_subdomain(k_particle).idx, domain_size, [x_min y_min z_min]);
                                    % Phase assignment
                                    %subdomain_phase(all_particles_subdomain(k_particle).idx) = phase(current_phase).code;
                                    % 3D arrays assignment (it will overwritte overlapping region as well: thus updating the other particles)
                                    microstructure3D = ThreeD_array_assignment(microstructure3D, current_phase, linear_indices, particle_id, number_voxel_assigned, number_voxel_idealshape, dx, dy, dz, angle_x_deg, angle_y_deg, angle_z_deg, x_center, y_center, z_center);
                                    microstructure3D.phase(linear_indices) = phase(current_phase).code;
                                    % Min-max
                                    [I_particle,J_particle,K_particle] = ind2sub(domain_size,linear_indices);
                                    particle_xyz_min_max = [min(I_particle) min(J_particle) min(K_particle) max(I_particle) max(J_particle) max(K_particle)];
                                else % Previous particles
                                    % Phase assignment
                                    other_id = all_particles_subdomain(k_particle).id;
                                    other_code = particle_info.Phase_code(other_id);
                                    other_phase = particle_info.Phase_number(other_id);
                                    
                                    % Indices of the modified particle
                                    linear_indices_sub_other_before = find(subdomain_id==other_id);
                                    [linear_indices_sub2domain_other_before] = Index_from_subdomain_to_domain(subdomain_ellipsoid, linear_indices_sub_other_before, domain_size, [x_min y_min z_min]);
                                    [~,ia,~] = intersect(particle_index(other_id).Index, linear_indices_sub2domain_other_before);
                                    linear_indices_other_wo_subdomain = particle_index(other_id).Index;
                                    linear_indices_other_wo_subdomain(ia)=[];
                                    [linear_indices_sub_other_new] = Index_from_subdomain_to_domain(subdomain_ellipsoid, all_particles_subdomain(k_particle).idx, domain_size, [x_min y_min z_min]);
                                    linear_indices_other = [linear_indices_other_wo_subdomain; linear_indices_sub_other_new];
                                    
                                    % Update number of voxels of the modified particle
                                    number_voxel_other_assigned = length(linear_indices_other);
                                    % Update volume fraction
                                    phase_volume_filled(other_phase,1) = phase_volume_filled(other_phase,1) - (particle_info.Number_voxel_currently_assigned(other_id)-number_voxel_other_assigned);
                                    % Update total number of unassigned voxels
                                    number_unassigned_voxels = number_unassigned_voxels + (particle_info.Number_voxel_currently_assigned(other_id)-number_voxel_other_assigned);
                                    % Update number of voxels
                                    particle_info.Number_voxel_currently_assigned(other_id) = number_voxel_other_assigned;
                                    % Update min-max
                                    [I_particle,J_particle,K_particle] = ind2sub(domain_size,linear_indices_other);
                                    particle_info.Extremities(other_id,:) = [min(I_particle) min(J_particle) min(K_particle) max(I_particle) max(J_particle) max(K_particle)];
                                    % Update 3D arrays assignment
                                    microstructure3D = Update_ThreeD_array_assignment(microstructure3D, linear_indices_other, number_voxel_other_assigned);
                                    % Update index
                                    particle_index(other_id).Index = uint32(linear_indices_other);
                                end
                            end
                            % Update
                            new_particle_hasbeen_generated = true;
                        end
%                         if continue_whileloop_and_go_to_next_randcenter
%                             continue
%                         end
                    end
                end
            else
                new_particle_hasbeen_generated = false; % Not necessary. For clarification sake.
            end % Particle has been generated or not below this point
            
            if new_particle_hasbeen_generated % Update
                phase_volume_filled(current_phase,1) = phase_volume_filled(current_phase,1) + number_voxel_assigned;
                number_unassigned_voxels = number_unassigned_voxels - number_voxel_assigned;
                phase_number_particle(current_phase,1) = phase_number_particle(current_phase,1) + 1;
                % Particle info update
                particle_info.Center_xyz(particle_id,:) = [x_center y_center z_center];
                particle_info.Phase_code(particle_id) = phase(current_phase).code;
                particle_info.Phase_number(particle_id) = current_phase;
                particle_info.Number_voxel_idealshape(particle_id) = number_voxel_idealshape; % Do not change after particle creation
                particle_info.Number_voxel_assigned_whencreated(particle_id) = number_voxel_assigned; % Do not change after particle creation
                particle_info.Number_voxel_currently_assigned(particle_id) = number_voxel_assigned; % May change after particle creation, due to overlapping
                particle_info.Extremities(particle_id,:) = particle_xyz_min_max; % May change after particle creation, due to overlapping
                particle_index(particle_id).Index = uint32(linear_indices); % May change after particle creation, due to overlapping
            end
        end % While phase not filled enough
        delete(h_waitbar); % Close the progess bar
    end % Loop over phase
end % Loop over pass
time_generation = toc(start_generation_time);
fprintf('Generation terminated. Elapsed time: %1.1fs\n',time_generation)


% Progression figure
if stoping_conditions.plot
    % Update array for figure
    time_since_start = toc(start_generation_time);
    time_progression = [time_progression time_since_start];
    volumefraction_progression = [volumefraction_progression phase_volume_filled(:,1)/voxel_number];
    particle_progression = [particle_progression  phase_number_particle(:,1)];
    
    time_progression_generationrate = [time_progression_generationrate time_progression(end-1)+((time_progression(end)-time_progression(end-1))/2)];;
    generationrate_progression = [generationrate_progression (volumefraction_progression(:,end)-volumefraction_progression(:,end-1)) ./ (time_progression(end)-time_progression(end-1)) ];
    generationrate_progression_average = [generationrate_progression_average mean(generationrate_progression,2)];
    particle_generationrate_progression = [particle_generationrate_progression (particle_progression(:,end)-particle_progression(:,end-1)) ./ (time_progression(end)-time_progression(end-1)) ];
    particle_generationrate_progression_average = [particle_generationrate_progression_average mean(particle_generationrate_progression,2)];
    
    if number_iteration>stoping_conditions.average_on_last_n_iterations
        time_progression_generationrate_avglastiterations = [time_progression_generationrate_avglastiterations time_since_start];
        generationrate_progression_lastiterations = [generationrate_progression_lastiterations mean(generationrate_progression(:,end-stoping_conditions.average_on_last_n_iterations:end),2)];
        particle_generationrate_progression_lastiterations = [particle_generationrate_progression_lastiterations mean(particle_generationrate_progression(:,end-stoping_conditions.average_on_last_n_iterations:end),2)];
    end

    % Update figure
    cla(sub_axes_progression_1); % Clear axes
    hold(sub_axes_progression_1,'on'); % Active subplot
    for current_phase_figure = 1:1:number_phase
        plot(sub_axes_progression_1,[0,time_since_start],[phase(current_phase_figure).volumefraction.total phase(current_phase_figure).volumefraction.total],'LineWidth',2,'LineStyle','--','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (inputs)']);
        plot(sub_axes_progression_1,time_progression,volumefraction_progression(current_phase_figure,:),'LineWidth',2,'LineStyle','-','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (generated)']);
    end
    ylim(sub_axes_progression_1,[0 +Inf])
    xlim(sub_axes_progression_1,[0 time_since_start])
    legend(sub_axes_progression_1,'Location','best');
    hold(sub_axes_progression_1,'off');
    
    cla(sub_axes_progression_2); % Clear axes
    hold(sub_axes_progression_2,'on'); % Active subplot
    for current_phase_figure = 1:1:number_phase
        plot(sub_axes_progression_2,time_progression_generationrate,generationrate_progression(current_phase_figure,:),'LineWidth',1,'LineStyle','-','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (instantaneous)']);
        plot(sub_axes_progression_2,time_progression_generationrate,generationrate_progression_average(current_phase_figure,:),'LineWidth',2,'LineStyle','--','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (global average)']);
        plot(sub_axes_progression_2,time_progression_generationrate_avglastiterations,generationrate_progression_lastiterations(current_phase_figure,:),'LineWidth',2,'LineStyle',':','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (average of last ' num2str(stoping_conditions.average_on_last_n_iterations) ' iterations)']);
    end
    xlim(sub_axes_progression_2,[0 time_since_start])
    legend(sub_axes_progression_2,'Location','best');
    hold(sub_axes_progression_2,'off');
    
    cla(sub_axes_progression_3); % Clear axes
    hold(sub_axes_progression_3,'on'); % Active subplot
    for current_phase_figure = 1:1:number_phase
        %plot(sub_axes_progression_3,time_progression_generationrate,particle_generationrate_progression(current_phase_figure,:),'LineWidth',1,'LineStyle','-','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (instantaneous)']);
        plot(sub_axes_progression_3,time_progression_generationrate,particle_generationrate_progression_average(current_phase_figure,:),'LineWidth',2,'LineStyle','--','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (global average)']);
        plot(sub_axes_progression_3,time_progression_generationrate_avglastiterations,particle_generationrate_progression_lastiterations(current_phase_figure,:),'LineWidth',2,'LineStyle',':','Color',c_(current_phase_figure,:),'DisplayName',[phase(current_phase_figure).name ' (average of last ' num2str(stoping_conditions.average_on_last_n_iterations) ' iterations)']);
    end
    xlim(sub_axes_progression_3,[0 time_since_start])
    legend(sub_axes_progression_3,'Location','best');
    hold(sub_axes_progression_3,'off');
    
    if ~isempty(save_options.folder) && save_options.save_progression
        function_savefig(Fig_global_progression, save_options.folder, ['Algorithm_progression_run_' num2str(save_options.run_number)]);
    end
end

% Verification
if do_verification
    function_verification_generatedmicrostructure(microstructure3D,phase,save_options)
end



%%
%% FUNCTIONS
%%

function idx = choose_randomly_from_histogramprobability(x)
sum_x = x;
for k=2:1:length(x)
    sum_x(k) = sum_x(k)+sum_x(k-1);
end
rn = rand*100; % Pick a random number between 0 and 100
tmp = sum_x<rn; % Deduce diameter
tmp2=find(tmp==1);
if isempty(tmp2)
    idx=1;
else
    idx=max(tmp2)+1;
end
end

function [new_probability] = update_probability (initial_probability, current_sum)
diff_ = max(initial_probability-current_sum,0);
new_probability = 100*diff_/sum(diff_);
% sum(new_)==100
end

 function [sub_min, sub_max, min_, max_] = bounds_subdomain_domain(dim,center,boundmin,boundmax)
 is_odd = mod(dim,2);
 if is_odd % Odd case
     min_ = center - floor(dim/2);
     max_ = center + floor(dim/2);       
 else % Even case
     min_ = center - (dim/2 - 1);
     max_ = center + dim/2;  
 end
 if min_<boundmin
     sub_min = boundmin-min_+1;
     min_ = boundmin;
 else
     sub_min = 1;
 end
 if max_>boundmax
     sub_max = dim - (max_-boundmax);
     max_ = boundmax;
 else
     sub_max = dim;
 end
 end

function [binary_ellipsoid] = create_ellipsoid(dx,dy,dz)
% Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix

Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
[X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
binary_ellipsoid = zeros(size(X));
binary_ellipsoid(R <= 1 ) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
% Remove one voxel tip, if any?

% Visualization
%slice(Y,X,Z,binary_ellipsoid,0,0,0);
%axis equal vis3d; colorbar;
%imagesc(binary_ellipsoid(:,:,1))
end

function [binary_ellipsoid] = rotate_domain(binary_ellipsoid,angle_x_deg, angle_y_deg, angle_z_deg)
% Angle rotation in radians
angle_x_rad=deg2rad(angle_x_deg);
angle_y_rad=deg2rad(angle_y_deg);
angle_z_rad=deg2rad(angle_z_deg);

% Rotation matrix
% https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
% Transformation matrix that rotates the image around the x-axis
tx = [1 0 0 0
    0 cos(angle_x_rad)  sin(angle_x_rad) 0
    0 -sin(angle_x_rad) cos(angle_x_rad) 0
    0             0              0       1];
% Transformation matrix that rotates the image around the y-axis
ty = [cos(angle_y_rad)  0      -sin(angle_y_rad)   0
    0             1              0     0
    sin(angle_y_rad)    0       cos(angle_y_rad)   0
    0             0              0     1];
% Transformation matrix that rotates the image around the z-axis
tz = [cos(angle_z_rad) sin(angle_z_rad) 0 0
    -sin(angle_z_rad)  cos(angle_z_rad) 0 0
    0 0 1 0
    0 0 0 1];
% Then pass the matrix to the affine3d object constructor.
tx_form = affine3d(tx);
ty_form = affine3d(ty);
tz_form = affine3d(tz);

% Apply the transformation to the image
if angle_x_deg~=0
    binary_ellipsoid = imwarp(binary_ellipsoid,tx_form);
end
if angle_y_deg~=0
    binary_ellipsoid = imwarp(binary_ellipsoid,ty_form);
end
if angle_z_deg~=0
    binary_ellipsoid = imwarp(binary_ellipsoid,tz_form);
end

% Crop
idx = find(binary_ellipsoid~=0);
[I,J,K]=ind2sub(size(binary_ellipsoid),idx);
binary_ellipsoid = binary_ellipsoid(min(I):max(I),min(J):max(J),min(K):max(K));
end


function [linear_indices_domain] = Index_from_subdomain_to_domain(binary_subdomain, idx, domain_size, subdomain_min_location)
if isempty(idx)
    idx = find(binary_subdomain==1); % Index of all voxels within the subdomain
end
[Isub,Jsub, Ksub] = ind2sub(size(binary_subdomain),idx); % Coordinates of all voxels within the subdomain
Idomain = Isub + subdomain_min_location(1)-1; % Coordinates of all voxels within the full domain
Jdomain = Jsub + subdomain_min_location(2)-1;
Kdomain = Ksub + subdomain_min_location(3)-1;
linear_indices_domain = sub2ind(domain_size, Idomain, Jdomain, Kdomain);
end

function [microstructure3D] = ThreeD_array_assignment(microstructure3D, phase_number, idx, particle_id, number_voxel, number_voxel_idealshape, dx, dy, dz, anglex, angley, anglez, x_center, y_center, z_center)
microstructure3D.particle_id(idx) = particle_id; % Assign particle id
microstructure3D.particle_volume_numbervoxel(idx) = number_voxel; % Assign particle volume
microstructure3D.particle_volume_numbervoxel_idealshape(idx) = number_voxel_idealshape; % Assign particle volume
if number_voxel~=1
    volume_ideal_ellipsoid = 4/3 * pi * dx * dy * dz / 8;
    sphere_equivalentdiameter_fromidealellipsoid = 2*((3*volume_ideal_ellipsoid)/(4*pi))^(1/3);
else
    volume_ideal_ellipsoid=1;
    sphere_equivalentdiameter_fromidealellipsoid=1;
end
microstructure3D.centroid(particle_id).xyz = [x_center, y_center, z_center]; % Centroid
microstructure3D.particle_volume_idealellipsoid(idx) = volume_ideal_ellipsoid; % Assign ellipsoid volume
microstructure3D.particle_equivalentspherediameter_fromidealellipsoid(idx) = sphere_equivalentdiameter_fromidealellipsoid; % Assign particle diameter
microstructure3D.phase_(phase_number).particle_length_x(idx) = dx; % Assign ellispoid lenght along axe 1
microstructure3D.phase_(phase_number).particle_length_y(idx) = dy; % Assign ellispoid lenght along axe 2
microstructure3D.phase_(phase_number).particle_length_z(idx) = dz; % Assign ellispoid lenght along axe 3
microstructure3D.phase_(phase_number).particle_elongation_dx_over_dy(idx) = round(dx/dy,4); % Ellispoid elongation dx/dy
microstructure3D.phase_(phase_number).particle_elongation_dx_over_dz(idx) = round(dx/dz,4); % Ellispoid elongation dx/dz
microstructure3D.phase_(phase_number).particle_rotation_x(idx) = round(anglex,4); % Assign ellispoid lenght along axe 1
microstructure3D.phase_(phase_number).particle_rotation_y(idx) = round(angley,4); % Assign ellispoid lenght along axe 2
microstructure3D.phase_(phase_number).particle_rotation_z(idx) = round(anglez,4); % Assign ellispoid lenght along axe 3
end

function [microstructure3D] = Update_ThreeD_array_assignment(microstructure3D, idx, number_voxel)
microstructure3D.particle_volume_numbervoxel(idx) = number_voxel; % Assign particle volume
end


function [all_Idx_assigned_voxels, Idx_assigned_voxels, Idx_unassigned_voxels_pickfrom] = update_unassigned_index(all_Idx_assigned_voxels, linear_indices, Idx_assigned_voxels, Idx_all_voxels, number_unassigned_voxels, Idx_unassigned_voxels_pickfrom)
all_Idx_assigned_voxels=[all_Idx_assigned_voxels linear_indices']; % All assigned voxels
Idx_assigned_voxels(linear_indices)=linear_indices; % 0 + all assigned voxels
r=Idx_all_voxels-Idx_assigned_voxels; % 0 for assigned voxels, correct index for unassigned voxels
%r=unique(r); % only one 0 for the assignes voxels, correct index for unassigned voxels
%r(1)=[]; % Remove 0. Correct index for unassigned voxels
r(r==0)=[]; % Equivalent, but faster.
tmp1 = randi(number_unassigned_voxels,length(all_Idx_assigned_voxels),1)'; % Choose n number (for n assigned voxels), from 1 to the number of unassigned voxels
tmp2=r(tmp1); % For these n number, assing index of unassignes voxels
Idx_unassigned_voxels_pickfrom(all_Idx_assigned_voxels)=tmp2; % Replace the index of assigned voxel with index of unassigned voxels. So that you will only pick index of unassigned voxels
end

function [probabiliy_histogram] = Calculate_histogram_probability_perslice(array3D, direction, uniquevalue, volumefractions)

array3D = round(array3D,4);

domain_size = size(array3D);
number_unique_value = length(uniquevalue);
if direction==1
    area_inplane = domain_size(2)*domain_size(3);
elseif direction==2
     area_inplane = domain_size(1)*domain_size(3);
else 
     area_inplane = domain_size(1)*domain_size(2);
end
% Initialize
probabiliy_histogram = zeros(domain_size(direction), number_unique_value);
for current_position=1:1:domain_size(direction) % Loop over postion
    for current_value=1:1:number_unique_value % Loop over value
        probabiliy_histogram(current_position,current_value) = sum(sum(array3D(:,:,current_position)==round(uniquevalue(current_value),4)));
    end
    probabiliy_histogram(current_position,:) = 100 * probabiliy_histogram(current_position,:) /(area_inplane*volumefractions(current_position));
end
end

function [statistics_histogram] = Statistics_from_histogram(histogram, values)
% values are sorted from low to high (histogram is the weight of the values)
tmp = [values; histogram];
tmp=sortrows(tmp',1)';
values = tmp(1,:);
histogram = tmp(2:end,:);

[n_position, ~] = size(histogram);
min_ = zeros(n_position,1);
mean_ = zeros(n_position,1);
max_ = zeros(n_position,1);
std_ = zeros(n_position,1);
sum_weight = sum(histogram,2);
for k=1:1:n_position
    mean_(k) = sum(histogram(k,:) .* values)/sum_weight(k);
    idx_ = find( histogram(k,:) ~= 0);
    min_(k) = values(idx_(1));
    max_(k) = values(idx_(end));
    std_(k) = std(values, histogram(k,:));
end
statistics_histogram = [min_ mean_ max_ std_];
end


%toc
end

