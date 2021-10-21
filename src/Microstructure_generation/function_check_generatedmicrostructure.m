function [] = function_check_generatedmicrostructure(microstructure3D,phase,OPTIONS)

OPTIONS.fontname = 'Times New Roman';
OPTIONS.grid = 'on';
OPTIONS.minorgrid = 'on';
OPTIONS.save_fig = true; % always true
OPTIONS.savefig_infig = true; % always true 
OPTIONS.savefig_informat = {'png'}; % always 'png'
OPTIONS.closefigureaftercreation = false;



%% VOLUME FRACTIONS


domain_size = size(microstructure3D.phase);
number_phase = length(phase);
voxel_number = prod(domain_size);
area_xy= domain_size(1)*domain_size(2);
% Total
for current_phase = 1:1:number_phase
   phase(current_phase).volumefraction.aftergeneration = sum(sum(sum( microstructure3D.phase == phase(current_phase).code))) / voxel_number;
   phase(current_phase).volumefraction. remain = phase(current_phase).volumefraction.total - phase(current_phase).volumefraction.aftergeneration;
end

% Along direction
for current_phase = 1:1:number_phase
    vf_position = zeros(1,domain_size(3));
    for current_position=1:1:domain_size(3) % Loop over postion
        vf_position(1,current_position) = sum(sum(microstructure3D.phase(:,:,current_position)==phase(current_phase).code));
    end
    vf_position = vf_position/area_xy;
    phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration = vf_position;
    % Equivalent to phase(current_phase).current_volumefraction.along_3rd_axis_allslices
end

% Save in a table


% Figure
Fig = figure; % Create figure
Fig.Name= 'Volume fractions along 3rd axis'; % Figure name
Fig.Color='white'; % Background colour
sub_axes = axes('Parent',Fig); % Create axes
% h_title=title (' ','FontName',OPTIONS.fontname,'FontSize',OPTIONS.size); % Set title font
% h_title.String= 'Volume fractions'; % Title string
% Plot graphs
x_ = [1:1:domain_size(3)];

for k = 1:1:2
    for current_phase=1:1:number_phase % Loop over phases
        sub_axes.ColorOrderIndex = current_phase;
        
        sub_axes = subplot(2,1,k,'Parent',Fig);
        hold(sub_axes,'on')
        if k ==1
            h_obtained = plot(x_, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration,'DisplayName',[phase(current_phase).name ', achieved'],'LineWidth',OPTIONS.figure_Linewidth,'LineStyle','-');
            h_consign  = plot(x_, phase(current_phase).volumefraction.along_3rd_axis_allslices,'DisplayName',[phase(current_phase).name ', consign'],'Color', h_obtained.Color,'LineWidth',OPTIONS.figure_Linewidth,'LineStyle','--');
            xlabel('Position along third axis (voxel length)');
            ylabel('Volume fractions');
            
        else
            h_difference = plot(x_,((phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration)-(phase(current_phase).volumefraction.along_3rd_axis_allslices)),'DisplayName',[phase(current_phase).name ', difference'],'LineWidth',OPTIONS.figure_Linewidth,'Color', h_obtained.Color);
            xlabel('Position along third axis (voxel length)');
            ylabel('Volume fractions difference');
            
        end
        if strcmp(OPTIONS.grid,true)
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
        end
        legend(sub_axes,'Location','best','FontSize',OPTIONS.legend_fontsize)
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.axe_fontsize); % Fontname and fontsize
        hold(sub_axes,'off')
    end
end
Vol_fig_title = sgtitle('Comparision between consign and generated volume fraction'); 
set(Vol_fig_title,'FontSize',OPTIONS.title_fontsize)

% if OPTIONS.save_fig == true % Save figure
%     filename= 'Volume_fractions_along_3rdaxis';
%     function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
% end
% if OPTIONS.closefigureaftercreation == true
%     close(Fig); % Do not keep open figures
% end

%% PARTICLE SIZE (DX,DY,DZ), ELONGATION (DX/DY, DX/DZ), and ROTATION (RX,RY,RZ)

% Round values
for current_phase = 1:1:number_phase
    microstructure3D.phase_(current_phase).particle_length_x = round(microstructure3D.phase_(current_phase).particle_length_x,4);
    microstructure3D.phase_(current_phase).particle_length_y = round(microstructure3D.phase_(current_phase).particle_length_y,4);
    microstructure3D.phase_(current_phase).particle_length_z = round(microstructure3D.phase_(current_phase).particle_length_z,4);
    microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy = round(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy,4);
    microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz = round(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz,4);
end
% Make sure complementary volume 0 will not interfer
val_ = -10;
for current_phase = 1:1:number_phase
    microstructure3D.phase_(current_phase).particle_length_x( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_length_y( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_length_z( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_rotation_x( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_rotation_y( microstructure3D.phase~=phase(current_phase).code ) = val_;
    microstructure3D.phase_(current_phase).particle_rotation_z( microstructure3D.phase~=phase(current_phase).code ) = val_;
end

% Min-max
min_dx = 9e9; max_dx = -9e9;
min_dy = 9e9; max_dy = -9e9;
min_dz = 9e9; max_dz = -9e9;
min_rx = 9e9; max_rx = -9e9;
min_ry = 9e9; max_ry = -9e9;
min_rz = 9e9; max_rz = -9e9;
min_dxdy = 9e9; max_dxdy = -9e9;
min_dxdz = 9e9; max_dxdz = -9e9;
for current_phase = 1:1:number_phase
    idx_ = microstructure3D.phase==phase(current_phase).code;
    
    min_dx = min(min_dx, min(microstructure3D.phase_(current_phase).particle_length_x(idx_)));
    min_dy = min(min_dy, min(microstructure3D.phase_(current_phase).particle_length_y(idx_)));
    min_dz = min(min_dz, min(microstructure3D.phase_(current_phase).particle_length_z(idx_)));
    min_rx = min(min_rx, min(microstructure3D.phase_(current_phase).particle_rotation_x(idx_)));
    min_ry = min(min_ry, min(microstructure3D.phase_(current_phase).particle_rotation_y(idx_)));
    min_rz = min(min_rz, min(microstructure3D.phase_(current_phase).particle_rotation_z(idx_)));   
    min_dxdy = min(min_dxdy, min(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy(idx_)));
    min_dxdz = min(min_dxdz, min(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz(idx_))); 
    
    max_dx = max(max_dx, max(microstructure3D.phase_(current_phase).particle_length_x(idx_)));
    max_dy = max(max_dy, max(microstructure3D.phase_(current_phase).particle_length_y(idx_)));
    max_dz = max(max_dz, max(microstructure3D.phase_(current_phase).particle_length_z(idx_)));    
    max_rx = max(max_rx, max(microstructure3D.phase_(current_phase).particle_rotation_x(idx_)));
    max_ry = max(max_ry, max(microstructure3D.phase_(current_phase).particle_rotation_y(idx_)));
    max_rz = max(max_rz, max(microstructure3D.phase_(current_phase).particle_rotation_z(idx_)));      
    max_dxdy = max(max_dxdy, max(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy(idx_)));
    max_dxdz = max(max_dxdz, max(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz(idx_)));
end

% Along direction
for current_phase = 1:1:number_phase
    % dx
    phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_length_x, 3, phase(current_phase).unique_dx_diameter, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
    % dy
    phase(current_phase).diameter_dy.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_length_y, 3, phase(current_phase).unique_dy_diameter, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
    % dz
    phase(current_phase).diameter_dz.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_length_z, 3, phase(current_phase).unique_dz_diameter, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
    
    % dx_dy
    phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dy, 3, phase(current_phase).unique_dxdy_elongation, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
    % dx_dz
    phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_elongation_dx_over_dz, 3, phase(current_phase).unique_dxdz_elongation, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
     
    % rx
    phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_rotation_x, 3, phase(current_phase).unique_rotation_x, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
    % ry
    phase(current_phase).rotation_y.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_rotation_y, 3, phase(current_phase).unique_rotation_y, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
    % rz
    phase(current_phase).rotation_z.along_3rd_axis_allslice_calculatedaftergeneration = Calculate_histogram_probability_perslice(microstructure3D.phase_(current_phase).particle_rotation_z, 3, phase(current_phase).unique_rotation_z, phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration);    
end

% Statistics
for current_phase=1:1:number_phase
  % dx
  phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_dx_diameter);
  % dy
  phase(current_phase).diameter_dy.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).diameter_dy.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_dy_diameter);
  % dz
  phase(current_phase).diameter_dz.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).diameter_dz.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_dz_diameter);
  % dxdy  
  phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_dxdy_elongation);
  % dxdz
  phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_dxdz_elongation);
  % rx
  phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_rotation_x);
  % ry
  phase(current_phase).rotation_y.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).rotation_y.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_rotation_y);
  % rz
  phase(current_phase).rotation_z.along_3rd_axis_allslice_calculatedaftergeneration_stats = Statistics_from_histogram(phase(current_phase).rotation_z.along_3rd_axis_allslice_calculatedaftergeneration, phase(current_phase).unique_rotation_z);
end


%% FIGURES
% Figure (particle size consign)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle size ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    %x_ = [1:1:domain_size(3)] * voxel_size_um;
    x_ = [1:1:domain_size(3)] * 1;
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Diameter dx in-plane repartition';
            %values_ = phase(current_phase).unique_dx_diameter * voxel_size_um;
            values_ = phase(current_phase).unique_dx_diameter * 1;
            y_obtained = phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).size_histogram.along_3rd_axis_allslices(2:end,:);
            lgd_title = 'Diameter in \mum';
        elseif id_axe==2
            str_title = 'Elongation dx/dy in-plane repartition';
            values_ = phase(current_phase).unique_dxdy_elongation;
            y_obtained = phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis_allslices(2:end,:);
            lgd_title = 'Elongation ratio';
        elseif id_axe==3
            str_title = 'Elongation dx/dz in-plane repartition';
            values_ = phase(current_phase).unique_dxdz_elongation;
            y_obtained = phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis_allslices(2:end,:);
            lgd_title = 'Elongation ratio';
        end
        h_title=title (str_title,'FontName',OPTIONS.fontname,'FontSize',16); % Set title font
        %xlabel('Position along third axis (\mum)');
        xlabel('Position along third axis (voxel length)');
        ylabel('Volume repartition (%)');
        for k_value = 1:1:length(values_)
            sub_axes.ColorOrderIndex = k_value;
            h_obtained = plot(x_, y_obtained(:,k_value),'DisplayName',[num2str(values_(k_value)) ', achieved']);
            h_consign  = plot(x_, y_consign(:,k_value),'DisplayName',[num2str(values_(k_value)) ', consign']);
            set(h_obtained, 'LineWidth',2,'LineStyle','-'); % Format
            set(h_consign, 'Color', h_obtained.Color,'LineWidth',2,'LineStyle','--'); % Format
        end
        % Legend
        lgd = legend(sub_axes,'Location','best','NumColumns',2);
        title(lgd,lgd_title)
        % - Grid
        if strcmp(OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',14); % Fontname and fontsize
        h_title.FontSize = 16; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
%     if OPTIONS.save_fig == true % Save figure
%         filename= ['DiameterElongation_along_3rdaxis_consign_' phase(current_phase).name];
%         function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
%     end
%     if OPTIONS.closefigureaftercreation == true
%         close(Fig); % Do not keep open figures
%     end
end

% Figure (particle orientation consign)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle orientation ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    %x_ = [1:1:domain_size(3)] * voxel_size_um;
    x_ = [1:1:domain_size(3)] * 1;
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Rotation rx in-plane repartition';
            values_ = phase(current_phase).unique_rotation_x;
            y_obtained = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis_allslices(2:end,:);
        elseif id_axe==2
            str_title = 'Rotation ry in-plane repartition';
            values_ = phase(current_phase).unique_rotation_y;
            y_obtained = phase(current_phase).rotation_y.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis_allslices(2:end,:);
        elseif id_axe==3
            str_title = 'Rotation rz in-plane repartition';
            values_ = phase(current_phase).unique_rotation_z;
            y_obtained = phase(current_phase).rotation_z.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis_allslices(2:end,:);
        end
        h_title=title (str_title,'FontName',OPTIONS.fontname,'FontSize',16); % Set title font
        %xlabel('Position along third axis (\mum)');
        xlabel('Position along third axis (voxel length)');
        ylabel('Volume repartition (%)');
        for k_value = 1:1:length(values_)
            sub_axes.ColorOrderIndex = k_value;
            h_obtained = plot(x_, y_obtained(:,k_value),'DisplayName',[num2str(values_(k_value)) ', achieved']);
            h_consign  = plot(x_, y_consign(:,k_value),'DisplayName',[num2str(values_(k_value)) ', consign']);
            set(h_obtained, 'LineWidth',2,'LineStyle','-'); % Format
            set(h_consign, 'Color', h_obtained.Color,'LineWidth',2,'LineStyle','--'); % Format
        end
        % Legend
        lgd = legend(sub_axes,'Location','best','NumColumns',2);
        lgd_title = 'Angle in degrees';
        title(lgd,lgd_title)
        % - Grid
        if strcmp(OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',14); % Fontname and fontsize
        h_title.FontSize = 16; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
%     if OPTIONS.save_fig == true % Save figure
%         filename= ['Rotation_along_3rdaxis_consign_' phase(current_phase).name];
%         function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
%     end
%     if OPTIONS.closefigureaftercreation == true
%         close(Fig); % Do not keep open figures
%     end
end

% Figure (particle size repartition)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle size repartition ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    %x_ = [1:1:domain_size(3)] * voxel_size_um;
    x_ = [1:1:domain_size(3)] * 1;
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Diameter dx in-plane repartition';
            ylabel('Volume repartition (%)');
            %lgd_title = 'Diameter in \mum;';
            lgd_title = 'Diameter in voxel length;';
            %values_ = phase(current_phase).unique_dx_diameter * voxel_size_um;
            values_ = phase(current_phase).unique_dx_diameter * 1;
            y_obtained = phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration;
            % Legend
            clear lgd_str
            for k_value=1:1:length(values_)
                lgd_str(k_value).str = num2str(values_(k_value));
            end            
        elseif id_axe==2
            str_title = 'Elongation dx/dy in-plane repartition';
            ylabel('Volume repartition (%)');
            lgd_title = 'Elongation ratio';
            values_ = phase(current_phase).unique_dxdy_elongation;
            y_obtained = phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration;
            % Legend
            clear lgd_str
            for k_value=1:1:length(values_)
                lgd_str(k_value).str = num2str(values_(k_value));
            end            
        elseif id_axe==3
            str_title = 'Elongation dx/dz in-plane repartition';
            ylabel('Volume repartition (%)');
            lgd_title = 'Elongation ratio';
            values_ = phase(current_phase).unique_dxdz_elongation;
            y_obtained = phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration;
            % Legend
            clear lgd_str
            for k_value=1:1:length(values_)
                lgd_str(k_value).str = num2str(values_(k_value));
            end
        end
        h_title=title (str_title,'FontName',OPTIONS.fontname,'FontSize',16); % Set title font
        %xlabel('Position along third axis (\mum)');
        xlabel('Position along third axis (voxel length)');
        % Plot
        h_obtained = area(y_obtained);
        ylim([0 100])
        lgd = legend(sub_axes,lgd_str.str,'Location','best','NumColumns',2);        
        title(lgd,lgd_title)        
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',14); % Fontname and fontsize
        h_title.FontSize = 16; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
%     if OPTIONS.save_fig == true % Save figure
%         filename= ['DiameterElongation_along_3rdaxis_' phase(current_phase).name];
%         function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
%     end
%     if OPTIONS.closefigureaftercreation == true
%         close(Fig); % Do not keep open figures
%     end
end

% Figure (particle orientation repartition)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle orientation repartition ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    %x_ = [1:1:domain_size(3)] * voxel_size_um;
    x_ = [1:1:domain_size(3)] * 1;
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Rotation rx in-plane repartition';
            ylabel('Volume repartition (%)');
            values_ = phase(current_phase).unique_rotation_x;
            y_obtained = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration;
            % Legend
            clear lgd_str
            for k_value=1:1:length(values_)
                lgd_str(k_value).str = num2str(values_(k_value));
            end            
        elseif id_axe==2
            str_title = 'Rotation ry in-plane repartition';
            ylabel('Volume repartition (%)');
            values_ = phase(current_phase).unique_rotation_y;
            y_obtained = phase(current_phase).rotation_y.along_3rd_axis_allslice_calculatedaftergeneration;
            % Legend
            clear lgd_str
            for k_value=1:1:length(values_)
                lgd_str(k_value).str = num2str(values_(k_value));
            end            
        elseif id_axe==3
            str_title = 'Rotation rz in-plane repartition';
            ylabel('Volume repartition (%)');
            values_ = phase(current_phase).unique_rotation_z;
            y_obtained = phase(current_phase).rotation_z.along_3rd_axis_allslice_calculatedaftergeneration;
            % Legend
            clear lgd_str
            for k_value=1:1:length(values_)
                lgd_str(k_value).str = num2str(values_(k_value));
            end
        end
        h_title=title (str_title,'FontName',OPTIONS.fontname,'FontSize',16); % Set title font
        %xlabel('Position along third axis (\mum)');
        xlabel('Position along third axis (voxel length)');
        % Plot
        h_obtained = area(y_obtained);
        ylim([0 100])
        lgd = legend(sub_axes,lgd_str.str,'Location','best','NumColumns',2);     
        lgd_title = 'Angle in degrees';        
        title(lgd,lgd_title)        
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',14); % Fontname and fontsize
        h_title.FontSize = 16; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
%     if OPTIONS.save_fig == true % Save figure
%         filename= ['Rotation_along_3rdaxis_' phase(current_phase).name];
%         function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
%     end
%     if OPTIONS.closefigureaftercreation == true
%         close(Fig); % Do not keep open figures
%     end
end

% Figure, Particle size and elongation (stats)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle size stats ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    %x_ = [1:1:domain_size(3)] * voxel_size_um;
    x_ = [1:1:domain_size(3)] * 1;
    % Plot graphs
    for id_axe=1:1:5
        sub_axes=subplot(2,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Diameter dx';
            %ylabel('Diameter (\mum)');
            ylabel('Diameter (voxel length)');
            %y_ = phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration_stats * voxel_size_um;
            y_ = phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration_stats * 1;
        elseif id_axe==2
            str_title = 'Diameter dy';
            %ylabel('Diameter (\mum)');
            ylabel('Diameter (voxel length)');
            %y_ = phase(current_phase).diameter_dy.along_3rd_axis_allslice_calculatedaftergeneration_stats * voxel_size_um;
            y_ = phase(current_phase).diameter_dy.along_3rd_axis_allslice_calculatedaftergeneration_stats * 1;
        elseif id_axe==3
            str_title = 'Diameter dz';
            %ylabel('Diameter (\mum)');
            ylabel('Diameter (voxel length)');
            %y_ = phase(current_phase).diameter_dz.along_3rd_axis_allslice_calculatedaftergeneration_stats * voxel_size_um;
            y_ = phase(current_phase).diameter_dz.along_3rd_axis_allslice_calculatedaftergeneration_stats * 1;
        elseif id_axe==4  
            str_title = 'Elongation dx/dy';
            ylabel('Ratio');
            y_ = phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        elseif id_axe==5
            str_title = 'Elongation dx/dz';
            ylabel('Ratio');
            y_ = phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        end
        h_title=title (str_title,'FontName',OPTIONS.fontname,'FontSize',16); % Set title font
        %xlabel('Position along third axis (\mum)');
        xlabel('Position along third axis (voxel length)');
        % Plot
        % Mean
        h_mean=plot(x_,y_(:,2)); % For the legend order
        % Extremums
        h_min=plot(x_,y_(:,1));
        h_max=plot(x_,y_(:,3));
        % Colors, thickness, markers
        set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        set(h_min, 'Color', 'b','LineStyle','--','LineWidth',1);
        set(h_max, 'Color', 'r','LineStyle','--','LineWidth',1);
        % Mean with error bar (+- standard deviation)
        h_mean_witherrorbar = errorbar(x_,y_(:,2),y_(:,4));
        set(h_mean_witherrorbar, 'Color', 'k','LineWidth',1);
        h_mean=plot(x_,y_(:,2)); % Plot over the other
        set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        lgd = legend(sub_axes,'Mean diameter with +- standard deviation','Minimum diameter', 'Maximun diameter','Location','best');        
        % - Grid
        if strcmp(OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',14); % Fontname and fontsize
        h_title.FontSize = 16; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
%     if OPTIONS.save_fig == true % Save figure
%         filename= ['DiameterElongation_along_3rdaxis_stats_' phase(current_phase).name];
%         function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
%     end
%     if OPTIONS.closefigureaftercreation == true
%         close(Fig); % Do not keep open figures
%     end
end

% Figure, Particle rotation (stats)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle rotation stats ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    %x_ = [1:1:domain_size(3)] * voxel_size_um;
    x_ = [1:1:domain_size(3)] * 1;
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Rotation rx';
            y_ = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        elseif id_axe==2
            str_title = 'Rotation ry';
            y_ = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        elseif id_axe==3
            str_title = 'Rotation rz';
            y_ = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        end
        ylabel('Angle (°)');        
        h_title=title (str_title,'FontName',OPTIONS.fontname,'FontSize',16); % Set title font
        %xlabel('Position along third axis (\mum)');
        xlabel('Position along third axis (voxel length)');
        % Plot
        % Mean
        h_mean=plot(x_,y_(:,2)); % For the legend order
        % Extremums
        h_min=plot(x_,y_(:,1));
        h_max=plot(x_,y_(:,3));
        % Colors, thickness, markers
        set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        set(h_min, 'Color', 'b','LineStyle','--','LineWidth',1);
        set(h_max, 'Color', 'r','LineStyle','--','LineWidth',1);
        % Mean with error bar (+- standard deviation)
        h_mean_witherrorbar = errorbar(x_,y_(:,2),y_(:,4));
        set(h_mean_witherrorbar, 'Color', 'k','LineWidth',1);
        h_mean=plot(x_,y_(:,2)); % Plot over the other
        set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        lgd = legend(sub_axes,'Mean rotation with +- standard deviation','Minimum rotation', 'Maximun rotation','Location','best');        
        % - Grid
        if strcmp(OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',14); % Fontname and fontsize
        h_title.FontSize = 16; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
%     if OPTIONS.save_fig == true % Save figure
%         filename= ['Rotation_along_3rdaxis_stats_' phase(current_phase).name];
%         function_savefig(Fig, save_folder, filename, OPTIONS); % Call function
%     end
%     if OPTIONS.closefigureaftercreation == true
%         close(Fig); % Do not keep open figures
%     end
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

end