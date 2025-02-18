function [] = function_verification_generatedmicrostructure(microstructure3D,phase,input_porosity, current_porosity, save_options)

domain_size = size(microstructure3D.phase);
x_ = [1:1:domain_size(3)];
number_phase = length(phase);
voxel_number = prod(domain_size);
area_xy= domain_size(1)*domain_size(2);
scrsz = get(0,'ScreenSize'); % Screen resolution
c_ = [colororder; rand(1000,3)];

%% VOLUME FRACTIONS
% Total
for current_phase = 1:1:number_phase
   phase(current_phase).volumefraction.aftergeneration = sum(sum(sum( microstructure3D.phase == phase(current_phase).code))) / voxel_number;
   phase(current_phase).volumefraction.remain = phase(current_phase).volumefraction.total - phase(current_phase).volumefraction.aftergeneration;
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

% Figure
Fig = figure; % Create figure
Fig.Name= 'Volume fractions along thickness (no cropping)'; % Figure name
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)*2/3]);
Fig.Color='white'; % Background colour
for id_axe=1:1:2
    if id_axe==1
        sub_axes=subplot(2,1,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        for current_phase = 1:1:number_phase
            plot(x_/x_(end), phase(current_phase).volumefraction.along_3rd_axis_allslices,'DisplayName',[phase(current_phase).name ' (inputs)'],'Color', c_(current_phase,:),'LineWidth',2,'LineStyle','--');
            plot(x_/x_(end), phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration,'DisplayName',[phase(current_phase).name ' (generated)'],'Color', c_(current_phase,:),'LineWidth',2,'LineStyle','-');
        end
        %if plot_porosity
            plot(current_porosity(1,:), input_porosity,'DisplayName','Porosity (inputs)','Color', 'k','LineWidth',2,'LineStyle','--');
            plot(current_porosity(1,:), current_porosity(2,:),'DisplayName','Porosity (generated)','Color', 'k','LineWidth',2,'LineStyle','-');
        %end
        ylabel('Volume fractions');
    elseif id_axe==2
        sub_axes=subplot(2,1,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot        
        for current_phase = 1:1:number_phase
            plot(x_/x_(end),(phase(current_phase).volumefraction.along_3rd_axis_allslice_calculatedaftergeneration - phase(current_phase).volumefraction.along_3rd_axis_allslices),'DisplayName',[phase(current_phase).name ' (error)'],'LineWidth',2,'Color', c_(current_phase,:));
        end
        %if plot_porosity
            plot(current_porosity(1,:),(current_porosity(2,:) - input_porosity),'DisplayName','Porosity (error)','LineWidth',2,'Color', 'k');
        %end
        ylabel('Volume fraction errors');
    end
    xlabel('Normalized position along direction 3 (thickness)');
    grid(sub_axes,'on'); % Display grid
    legend('Location','best');
    set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
    hold(sub_axes,'off'); % Relase axe,'off'); % Relase axe
end
sgtitle(Fig,'Volume fractions along direction 3 (thickness)','FontWeight','bold','FontSize',16,'FontName','Times new roman');
if ~isempty(save_options.folder) && save_options.save_verification
    function_savefig(Fig, save_options.folder, ['Volumefraction_run_' num2str(save_options.run_number)]);
end


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


% % FIGURES
% Figure (particle size consign)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle size of ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*6/7 scrsz(4)*2/3]);
    % - Create axes as a subplot
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Particle diameter dx';
            values_ = phase(current_phase).unique_dx_diameter;
            y_obtained = phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).size_histogram.along_3rd_axis_allslices(2:end,:);
            lgd_title = 'Diameter in voxel length';
        elseif id_axe==2
            str_title = 'Particle elongation dx/dy';
            values_ = phase(current_phase).unique_dxdy_elongation;
            y_obtained = phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).elongation_histogram_dx_over_dy.along_3rd_axis_allslices(2:end,:);
            lgd_title = 'Elongation ratio';
        elseif id_axe==3
            str_title = 'Particle elongation dx/dz';
            values_ = phase(current_phase).unique_dxdz_elongation;
            y_obtained = phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).elongation_histogram_dx_over_dz.along_3rd_axis_allslices(2:end,:);
            lgd_title = 'Elongation ratio';
        end
        h_title=title (str_title,'FontName','Times new roman','FontSize',14); % Set title font
        xlabel('Normalized position along direction 3 (thickness)');
        ylabel('Volume repartition (%)');
        for k_value = 1:1:length(values_)
            plot(x_/x_(end), y_consign(:,k_value),'DisplayName',[num2str(values_(k_value)) ' (inputs)'], 'Color', c_(k_value,:),'LineWidth',2,'LineStyle','--');
        end   
        for k_value = 1:1:length(values_)
            plot(x_/x_(end), y_obtained(:,k_value),'DisplayName',[num2str(values_(k_value)) ' (generated)'],'Color', c_(k_value,:),'LineWidth',2,'LineStyle','-');
        end          
        % Legend
        lgd = legend(sub_axes,'Location','best','NumColumns',2);
        title(lgd,lgd_title)
        % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
        h_title.FontSize = 14; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,{['Particle diameter of ' phase(current_phase).name],'In-plane distribution'},'FontWeight','bold','FontSize',16,'FontName','Times new roman');
    if ~isempty(save_options.folder) && save_options.save_verification
        phasename = function_remove_emptyandspecialcharacter_string(phase(current_phase).name);
        function_savefig(Fig, save_options.folder, ['Particlediameter_' phasename '_run_' num2str(save_options.run_number)]);
    end
end


% Figure (particle orientation consign)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle rotation of ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*6/7 scrsz(4)*2/3]);
    % - Create axes as a subplot
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Particle rotation Rx';
            values_ = phase(current_phase).unique_rotation_x;
            y_obtained = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).orientation_histogram_angledeg_x.along_3rd_axis_allslices(2:end,:);
        elseif id_axe==2
            str_title = 'Particle rotation Ry';
            values_ = phase(current_phase).unique_rotation_y;
            y_obtained = phase(current_phase).rotation_y.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).orientation_histogram_angledeg_y.along_3rd_axis_allslices(2:end,:);
        elseif id_axe==3
            str_title = 'Particle rotation Rz';
            values_ = phase(current_phase).unique_rotation_z;
            y_obtained = phase(current_phase).rotation_z.along_3rd_axis_allslice_calculatedaftergeneration;
            y_consign = phase(current_phase).orientation_histogram_angledeg_z.along_3rd_axis_allslices(2:end,:);
        end
        h_title=title (str_title,'FontName','Times new roman','FontSize',14); % Set title font
        xlabel('Normalized position along direction 3 (thickness)');
        ylabel('Volume repartition (%)');
        for k_value = 1:1:length(values_)
            plot(x_/x_(end), y_consign(:,k_value),'DisplayName',[num2str(values_(k_value)) ' (inputs)'], 'Color', c_(k_value,:),'LineWidth',2,'LineStyle','--');
        end
        for k_value = 1:1:length(values_)
            plot(x_/x_(end), y_obtained(:,k_value),'DisplayName',[num2str(values_(k_value)) ' (generated)'],'Color', c_(k_value,:),'LineWidth',2,'LineStyle','-');
        end        
        % Legend
        lgd = legend(sub_axes,'Location','best','NumColumns',2);
        lgd_title = 'Angle in degrees';
        title(lgd,lgd_title)
        % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
        h_title.FontSize = 14; % Set title fontsize
        lgd.FontSize = 12; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,{['Particle rotation of ' phase(current_phase).name],'In-plane distribution'},'FontWeight','bold','FontSize',16,'FontName','Times new roman');
    if ~isempty(save_options.folder) && save_options.save_verification
        phasename = function_remove_emptyandspecialcharacter_string(phase(current_phase).name);
        function_savefig(Fig, save_options.folder, ['Particlerotation_' phasename '_run_' num2str(save_options.run_number)]);
    end
end

% Figure, Particle size and elongation (stats)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle size statistics of ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*6/7 scrsz(4)*2/3]);
    % - Create axes as a subplot
    % Plot graphs
    for id_axe=1:1:5
        sub_axes=subplot(2,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Diameter Dx';
            ylabel('Particle diameter (voxel length)');
            y_ = phase(current_phase).diameter_dx.along_3rd_axis_allslice_calculatedaftergeneration_stats * 1;
        elseif id_axe==2
            str_title = 'Diameter Dy';
            ylabel('Particle diameter (voxel length)');
            y_ = phase(current_phase).diameter_dy.along_3rd_axis_allslice_calculatedaftergeneration_stats * 1;
        elseif id_axe==3
            str_title = 'Diameter Dz';
            ylabel('Particle diameter (voxel length)');
            y_ = phase(current_phase).diameter_dz.along_3rd_axis_allslice_calculatedaftergeneration_stats * 1;
        elseif id_axe==4  
            str_title = 'Elongation Dx/Dy';
            ylabel('Ratio');
            y_ = phase(current_phase).elongation_dxdy.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        elseif id_axe==5
            str_title = 'Elongation Dx/Dz';
            ylabel('Ratio');
            y_ = phase(current_phase).elongation_dxdz.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        end
        h_title=title (str_title,'FontName','Times new roman','FontSize',14); % Set title font
        xlabel('Normalized position along direction 3 (thickness)');
        % Mean
        h_mean=plot(x_/x_(end),y_(:,2)); % For the legend order
        % Extremums
        h_min=plot(x_/x_(end),y_(:,1));
        h_max=plot(x_/x_(end),y_(:,3));
        % Colors, thickness, markers
        set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        set(h_min, 'Color', 'b','LineStyle','--','LineWidth',1);
        set(h_max, 'Color', 'r','LineStyle','--','LineWidth',1);

        % Mean with error bar (+- standard deviation)
        % Old approach
        % h_mean_witherrorbar = errorbar(x_/x_(end),y_(:,2),y_(:,4));
        % set(h_mean_witherrorbar, 'Color', 'k','LineWidth',1);
        % h_mean=plot(x_/x_(end),y_(:,2)); % Plot over the other
        % set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        % New approach
        std_min = (y_(:,2) - y_(:,4))';
        std_max = (y_(:,2) + y_(:,4))';
        x2 = [(x_/x_(end)), fliplr((x_/x_(end)))];
        inBetween = [std_min, fliplr(std_max)];
        h_=fill(x2, inBetween, 'k');
        set(h_,'LineStyle','none','FaceAlpha',0.25);

        lgd = legend(sub_axes,'Mean diameter with +- standard deviation','Minimum diameter', 'Maximun diameter','Location','best');        
        % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
        h_title.FontSize = 14; % Set title fontsize
        lgd.FontSize = 10; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,{['Particle diameter statistics of ' phase(current_phase).name],'In-plane distribution'},'FontWeight','bold','FontSize',16,'FontName','Times new roman');
    if ~isempty(save_options.folder) && save_options.save_verification
        phasename = function_remove_emptyandspecialcharacter_string(phase(current_phase).name);
        function_savefig(Fig, save_options.folder, ['Particlediameter_stats_' phasename '_run_' num2str(save_options.run_number)]);
    end
end

% Figure, Particle rotation (stats)
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Particle rotation stats ' phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*6/7 scrsz(4)*2/3]);
    % - Create axes as a subplot
    % Plot graphs
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            str_title = 'Particle rotation Rx';
            y_ = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        elseif id_axe==2
            str_title = 'Particle rotation Ry';
            y_ = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        elseif id_axe==3
            str_title = 'Particle rotation Rz';
            y_ = phase(current_phase).rotation_x.along_3rd_axis_allslice_calculatedaftergeneration_stats;
        end
        ylabel('Angle (°)');        
        h_title=title (str_title,'FontName','Times new roman','FontSize',14); % Set title font
        xlabel('Normalized position along direction 3 (thickness)');
        % Plot
        % Mean
        h_mean=plot(x_/x_(end),y_(:,2)); % For the legend order
        % Extremums
        h_min=plot(x_/x_(end),y_(:,1));
        h_max=plot(x_/x_(end),y_(:,3));
        % Colors, thickness, markers
        set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);
        set(h_min, 'Color', 'b','LineStyle','--','LineWidth',1);
        set(h_max, 'Color', 'r','LineStyle','--','LineWidth',1);
        
        % Mean with error bar (+- standard deviation)
        % Old approach
        % h_mean_witherrorbar = errorbar(x_/x_(end),y_(:,2),y_(:,4));
        % set(h_mean_witherrorbar, 'Color', 'k','LineWidth',1);
        % h_mean=plot(x_/x_(end),y_(:,2)); % Plot over the other
        % set(h_mean, 'Color', 'k','LineWidth',1,'LineWidth',2);

        % New approach
        std_min = (y_(:,2) - y_(:,4))';
        std_max = (y_(:,2) + y_(:,4))';
        x2 = [(x_/x_(end)), fliplr((x_/x_(end)))];
        inBetween = [std_min, fliplr(std_max)];
        h_=fill(x2, inBetween, 'k');
        set(h_,'LineStyle','none','FaceAlpha',0.25);

        lgd = legend(sub_axes,'Mean rotation with +- standard deviation','Minimum rotation', 'Maximun rotation','Location','best');        
        % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
        h_title.FontSize = 14; % Set title fontsize
        lgd.FontSize = 10; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,{['Particle rotation statistics of ' phase(current_phase).name],'In-plane distribution'},'FontWeight','bold','FontSize',16,'FontName','Times new roman');
    if ~isempty(save_options.folder) && save_options.save_verification
        phasename = function_remove_emptyandspecialcharacter_string(phase(current_phase).name);
        function_savefig(Fig, save_options.folder, ['Particlerotation_stats_' phasename '_run_' num2str(save_options.run_number)]);
    end
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
probabiliy_histogram(isnan(probabiliy_histogram)) = 0; % NaN when divided by 0 volume fraction
end

function [statistics_histogram] = Statistics_from_histogram(histogram_, values)
% values are sorted from low to high (histogram is the weight of the values)
tmp = [values; histogram_];
tmp=sortrows(tmp',1)';
values = tmp(1,:);
histogram_ = tmp(2:end,:);

[n_position, ~] = size(histogram_);
min_ = zeros(n_position,1);
mean_ = zeros(n_position,1);
max_ = zeros(n_position,1);
std_ = zeros(n_position,1);
sum_weight = sum(histogram_,2);
for k=1:1:n_position
    idx_ = find( histogram_(k,:) ~= 0);
    if ~isempty(idx_)
        mean_(k) = sum(histogram_(k,:) .* values)/sum_weight(k);
        min_(k) = values(idx_(1));
        max_(k) = values(idx_(end));
        std_(k) = std(values, histogram_(k,:));
    end
end
statistics_histogram = [min_ mean_ max_ std_];
end

end