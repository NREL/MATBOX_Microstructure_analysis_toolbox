function [D_PSD_particle_size,Label_lake] = Function_Discrete_particle_size_watershed_immersion_algoritv2(binary_phase,pars)

%% REFERENCE
% This algorithm is inspired by this publication (but the code here is original)
% Watersheds in  Digital  Spaces: An  Efficient  Algorithm  Based on  Immersion  Simulations 
% Luc  Vincent  and  Pierre  Soille
% IEEE  TRANSACTIONS  ON  PATTERN  ANALYSIS  AND  MACHINE  INTELLIGENCE, VOL.  13,  NO.  6,  JUNE  1991

cpsd_refining = pars.cpsd_refining;
downscaling_factor = pars.downscaling_factor;
details_convergence = pars.details_convergence;
visualize_2D = pars.visualize_2D;

%% DOWNSCALE
if downscaling_factor~=1
    sz = size(binary_phase);
    tmp = [sz(1)-100:1:sz(1)];
    tmp2 = mod(tmp,downscaling_factor);
    idx = max(find(tmp2==0));
    cropx = length(tmp)-idx;
    tmp = [sz(2)-100:1:sz(2)];
    tmp2 = mod(tmp,downscaling_factor);
    idx = max(find(tmp2==0));
    cropy = length(tmp)-idx;    
    binary_phase(end-cropx+1:end, :) = [];
    binary_phase(:, end-cropy+1:end) = [];

    binary_phase_initial = binary_phase;
    parameters_scaling.scaling_factor = downscaling_factor;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = 0;
    % Scale
    binary_phase = function_scaling(binary_phase,parameters_scaling);
end

%% DOMAIN SIZE
% Domain size of the microstructure
Domain_size=size(binary_phase);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
else
    visualize_2D = false;
end
scrsz = get(0,'ScreenSize'); % Screen resolution

if visualize_2D
    lake_id_RGB_color = 0.1 + (0.9-0.1).*rand(1e6,3);
    lake_id_grey_color = sum(lake_id_RGB_color,2)/3;
    desktop=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
    idx_source=[]; existingvalues=[];
    % Edges
    [index_border_phase,~,~,~] = Function_identify_edges(binary_phase);    
    
    % Prepare Video
    video_handle = VideoWriter([[desktop '\'] 'Illustration_2D_watershed'],'mpeg-4');
    set(video_handle,'Quality',100); % Set video quality
    set(video_handle,'FrameRate',5); % Set video framerate
    open(video_handle)
    Fig_video = figure; % Create figure
    Fig_video.Name= 'Illustration of two-dimensional watershed method'; % Figure name
    Fig_video.Color='white'; % Background colour
    set(Fig_video,'position',scrsz); % Full screen figure
    sgtitle(Fig_video,'Animation of two-dimensional watershed method','FontName','Times New Roman','FontSize',18);
    axe_video(1).s = subplot(1,2,1,'Parent',Fig_video); % Create axes 
    axe_video(2).s = subplot(1,2,2,'Parent',Fig_video);
    y_video_l = []; y_video_2 = []; num_voxel = sum(sum(sum(binary_phase==1)));
    
    % Prepare illustration 1
    Fig_illustration = figure; % Create figure
    Fig_illustration.Name= 'Illustration of two-dimensional watershed method'; % Figure name
    Fig_illustration.Color='white'; % Background colour
    set(Fig_illustration,'position',scrsz); % Full screen figure
    sgtitle(Fig_illustration,'Illustration of two-dimensional watershed method','FontName','Times New Roman','FontSize',18);
    sub_axes_illustration(1).s=subplot(2,2,1,'Parent',Fig_illustration); % Create axes    
    hold(sub_axes_illustration(1).s,'on'); % Active subplot
    imagesc(sub_axes_illustration(1).s,binary_phase);
    colormap(sub_axes_illustration(1).s,'gray');
    axis(sub_axes_illustration(1).s,'tight'); % Fit the axes box
    axis(sub_axes_illustration(1).s,'equal'); % Aspect ratio is 1:1
    % Label
    xlabel('Voxels');
    ylabel('Voxels');
    set(sub_axes_illustration(1).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
    title('Binary image','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(1).s);
    hold(sub_axes_illustration(1).s,'off'); % Active subplot
end


%% CALCULATE EUCLIDEAN DISTANCE MAP
% We use the matlab built-in function bwdist to calculate the distance transform (or distance map) with an euclidean distance
% bwdist(BW) : For each pixel in BW, the distance transform assigns a number that is the distance between that pixel and the nearest nonzero pixel of BW
% then, bwdist(~binary_phase,'euclidean') :
Distance_map=bwdist(~binary_phase,'euclidean');

% Save time
%length(unique(Distance_map))
%Distance_map = round(Distance_map);
length(unique(Distance_map))

% Topography analogy:
Distance_map = -Distance_map; % Inverse the distance sign: 0 is the reference altitude, minimum is the lower altitude point
Distance_map(~binary_phase) = +Inf; % A negative infinite value is attributed to the complementary

if visualize_2D
    sub_axes_illustration(2).s=subplot(2,2,2,'Parent',Fig_illustration); % Create axes  
    hold(sub_axes_illustration(2).s,'on'); % Active subplot
    imagesc(sub_axes_illustration(2).s,Distance_map);
    custom_map = jet; custom_map(end,:)=[1 1 1];
    colormap(sub_axes_illustration(2).s,custom_map);
    axis(sub_axes_illustration(2).s,'tight'); % Fit the axes box
    axis(sub_axes_illustration(2).s,'equal'); % Aspect ratio is 1:1
    colorbar
    % Label
    xlabel('Voxels');
    ylabel('Voxels');
    set(sub_axes_illustration(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
    title('Altitude map','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(2).s);
    box(sub_axes_illustration(2).s,'on')
    hold(sub_axes_illustration(2).s,'off'); % Active subplot
end

%% SORTING THE DISTANCE MAP VALUES
all_altitude = unique(Distance_map); % Get all altitudes
all_altitude = sortrows(all_altitude,1); % Matlab has already sorted this array from -inf to +inf. Let's be sure of this
all_altitude(end)=[]; % Remove last value value (= +Inf)
number_altitude = length(all_altitude); % Number of altitude level

%% LABELING LAKE ITERATIVELY
Label_lake = zeros(Domain_size(1), Domain_size(2), Domain_size(3)); % Initialisation
if details_convergence
    disp 'Flooding simulation in progress...(%)';
end
next_step_display=1;
t_altitudeiteration = tic;
for current_altitude_iteration = 1:1:number_altitude % We start from the lower position and go up
    % Progress of the calculation in percents
    progress_force_calculation = 100*current_altitude_iteration/number_altitude;
    if progress_force_calculation>=next_step_display && details_convergence
        fprintf (' %4.1f',progress_force_calculation);
        next_step_display=next_step_display+1;
    end
    current_altitude = all_altitude(current_altitude_iteration); % Current altitude value
    
    % Binary image of this level of altitude
    binary_altitude_map = zeros(Domain_size(1), Domain_size(2), Domain_size(3)); % Initialisation
    binary_altitude_map(Distance_map<=current_altitude)=1; % Mark 1 when voxel are located at the current altitude OR LOWER: LAKES ARE BEING FILLING FROM THEIR DEEPEST POINT
    
    % Label the current level image
    % L = bwlabeln(BW,6)
    % L is containing labels for the connected components in BW
    % L is the same size as BW
    % The elements of L are integer values greater than or equal to 0.
    % The pixels labeled 0 are the background. The pixels labeled 1 make up one connected cluster. The pixels labeled 2 make up a second connected cluster, and so on
    Label_current_altitude = bwlabeln(binary_altitude_map,6);
    
    unique_cluster = unique(Label_current_altitude); % Number of different cluster identified
    unique_cluster(1)=[]; % Remove 0 (complementary phase)
    number_cluster = length(unique_cluster); % Number of different cluster

    % We label the different lake
    if current_altitude_iteration>1
        for current_iteration = 1:1:number_cluster % Loop over all the cluster identified
            Current_label = unique_cluster(current_iteration); % Current label id
            binary_cluster = zeros(Domain_size(1), Domain_size(2), Domain_size(3)); % Create binary cluster
            binary_cluster(Label_current_altitude==Current_label)=1;
            index_cluster = find(binary_cluster==1); % Index voxels of the cluster
            intersection = binary_cluster.*Label_lake; % Calculate intersection with the already identified lake
            index_intersection = find(intersection~=0); % Index of the intersected voxels
            number_intersection = length(index_intersection); % Number of voxels intersected
            
            % % % Four cases are possible
            
            % % Case 1 and 2: intersection is empty AND...
            if number_intersection==0
                % Binary image of the existing lake
                binary_existing_lake = zeros(Domain_size(1), Domain_size(2), Domain_size(3));
                binary_existing_lake(Label_lake~=0)=1;
                % Distance bewteen all non attributes voxels and the existing lakes
                [Dist_existing_lake,Index_nearest_lake] = bwdist(binary_existing_lake,'chessboard');
                % Restrict result to the current new cluster
                Dist_existing_lake=Dist_existing_lake.*binary_cluster;
                % All distance
                unique_distance = unique(Dist_existing_lake);
                % Remove 0
                unique_distance(1)=[];
                % Minimum distance
                min_distance = unique_distance(1);
                if min_distance>2
                    % Case 1: AND... far from existing lake: this cluster is a new lake
                    % Set new id for this new lake
                    max_lake_id = max(unique(Label_lake))+1;
                    % Attribute this new id to the cluster
                    Label_lake(index_cluster)=max_lake_id;
                else
                    % Case 2: AND... 'adjacent' with one or more existing lake
                    % Index of voxels of the cluster close to the existing lake
                    index_near = find(Dist_existing_lake==min_distance);
                    % Id of these lake
                    index_adjacent_lake = Index_nearest_lake(index_near);
                    id_adjacent_lake = Label_lake(index_adjacent_lake);
                    % Number of adjacent lake
                    number_adjacent_lake = length(id_adjacent_lake);
                    if number_adjacent_lake==1
                        % Simple case: Attribute all voxels of the label to this unique lake
                        Label_lake(index_cluster)=id_adjacent_lake(1);
                    else
                        % Voxel will be attributed to the closest lake
                        for current_voxel=1:1:length(index_cluster)
                            choose_lake = zeros(number_adjacent_lake,8);
                            [x_tmp,y_tmp,z_tmp] = ind2sub(Domain_size,index_cluster(current_voxel));
                            choose_lake(:,1)=x_tmp;
                            choose_lake(:,2)=y_tmp;
                            choose_lake(:,3)=z_tmp;
                            for current_lake=1:1:number_adjacent_lake
                                [x_tmp,y_tmp,z_tmp] = ind2sub(Domain_size,index_adjacent_lake(current_lake));
                                choose_lake(current_lake,4)=x_tmp;
                                choose_lake(current_lake,5)=y_tmp;
                                choose_lake(current_lake,6)=z_tmp;
                                choose_lake(current_lake,7)=Label_lake(index_adjacent_lake(current_lake));
                            end
                            % Calculate distance
                            choose_lake(:,8)= (choose_lake(:,1)-choose_lake(:,4)).^2 + (choose_lake(:,2)-choose_lake(:,5)).^2 + (choose_lake(:,3)-choose_lake(:,5)).^2;
                            % Sort by increasing distance
                            choose_lake = sortrows(choose_lake,8);
                            % Attribute the current voxel
                            Label_lake(index_cluster(current_voxel))=choose_lake(1,7);
                        end
                    end
                end
            end
            
            % % Case 3 and 4: intersection is not empty AND...
            if number_intersection>0
                % Does intersected voxels belong to different lake?
                different_lake = unique(Label_lake(index_intersection));
                % Number different lake
                number_different_lake = length(different_lake);
                % % Case 3: AND... intersection with only one existing lake
                if number_different_lake==1
                    % Attribute all voxels of the label to this unique lake
                    Label_lake(index_cluster)=different_lake;
                end
                % % Case 4: AND... intersection with more than 1 existing lake
                if number_different_lake>1
                    % Each voxel of the cluster will be attributed to the closest lake
                    % % Step 1: Get subdomain that contains the cluster
                    % Get all coordinates
                    [C1,C2,C3] = ind2sub(Domain_size,index_cluster);
                    % Min, max coordinates
                    x_min = min(C1); y_min = min(C2); z_min = min(C3);
                    x_max = max(C1); y_max = max(C2); z_max = max(C3);
                    % Extract subdomain for the binary cluster
                    subdomain_cluster = zeros(x_max-x_min+1,y_max-y_min+1,z_max-z_min+1);
                    subdomain_cluster = binary_cluster(x_min:x_max,y_min:y_max,z_min:z_max);
                    % Extract subdomain for the lake
                    subdomain_lake = zeros(x_max-x_min+1,y_max-y_min+1,z_max-z_min+1);
                    subdomain_lake = Label_lake(x_min:x_max,y_min:y_max,z_min:z_max);
                    % Extract subdomain for the intersection
                    subdomain_intersection = zeros(x_max-x_min+1,y_max-y_min+1,z_max-z_min+1);
                    subdomain_intersection = intersection(x_min:x_max,y_min:y_max,z_min:z_max);
                    % % Step 2: Prepare the matrix that will calculate the shortest distance
                    subdomain_distance = subdomain_intersection;
                    subdomain_distance (subdomain_distance~=0)=1;
                    % Step 3: Calculate shortest distance from voxel to lake
                    [~,Distance_map_nearest_lake] = bwdist(subdomain_distance,'euclidean');
                    % Step 4: Mark the voxels that does not belong to the current cluster
                    Distance_map_nearest_lake=double(Distance_map_nearest_lake);
                    Distance_map_nearest_lake(subdomain_cluster==0)=-1;
                    % Step 5: Attribute to each voxel of the cluster the nearest lake
                    attraction_point = unique(Distance_map_nearest_lake);
                    number_attraction_point = length(attraction_point);
                    for current_=1:1:number_attraction_point
                        current_attraction_point = attraction_point(current_);
                        if current_attraction_point~= -1
                            % Find all voxels attracted by this point
                            index_attracted = find(Distance_map_nearest_lake==current_attraction_point);
                            % Attribute these voxel to the corresponding lake
                            subdomain_lake(index_attracted)=subdomain_lake(current_attraction_point);
                        end
                    end
                    % Step 6: Re-introduced the subomain to the whole domain
                    Label_lake(x_min:x_max,y_min:y_max,z_min:z_max)=subdomain_lake;
                end
            end
        end
    else
        % Firt iteration?
        % Then all the different cluster at the minimum voxel identifies
        % different lake (their lower altitude)
        Label_lake = Label_current_altitude;
    end
    
    % Update illustration
    if visualize_2D
        % Initializaion
        slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
        slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
        slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
        slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
        slice_grey = zeros(Domain_size(1),Domain_size(2)); % Blue color map
        currentvalues = unique(Label_lake);
        for k=1:1:length(currentvalues)
            idx=find(Label_lake==currentvalues(k));
            if isempty( intersect(existingvalues,currentvalues(k)) ) && currentvalues(k)~=0
                idx_source=[idx_source; idx];
                existingvalues=unique([existingvalues currentvalues(k)]);
            end
            if currentvalues(k)==0
                slice_r(idx) = 1; slice_g(idx) = 1; slice_b(idx) = 1; % Background
                slice_grey(idx) = 1;
            else
                slice_r(idx) = lake_id_RGB_color(k,1); slice_g(idx) = lake_id_RGB_color(k,2); slice_b(idx) = lake_id_RGB_color(k,3); % Background                
                slice_grey(idx) = lake_id_grey_color(k);
            end
        end
        % Border
        slice_r(index_border_phase) = 0; slice_g(index_border_phase) = 0; slice_b(index_border_phase) = 0;
        slice_grey(index_border_phase) = 0;        
        % Source
        slice_r(idx_source) = 1; slice_g(idx_source) = 0; slice_b(idx_source) = 0;
        slice_color(:,:,1)=slice_r; % Attribute RGB color
        slice_color(:,:,2)=slice_g;
        slice_color(:,:,3)=slice_b;
        hold(axe_video(1).s,'on');
        slice_image = image(slice_color,'parent',axe_video(1).s); % Display the slice
        imagesc(axe_video(1).s,Distance_map, 'AlphaData', .25);
        colormap(axe_video(1).s,custom_map);
        axis(axe_video(1).s,'tight'); % Fit the axes box
        axis(axe_video(1).s,'equal'); % Aspect ratio is 1:1
        % Label
        xlabel(axe_video(1).s,'Voxels');
        ylabel(axe_video(1).s,'Voxels');
        set(axe_video(1).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title({'Filling the volume from lake sources',['Number of bassins identified: ' num2str(length(currentvalues)-1,'%i')]},'FontName','Times New Roman','FontSize',18,'Parent',axe_video(1).s);        
        set(axe_video(1).s,'YDir','normal')
        hold(axe_video(1).s,'off');
        
        % Progression
        volumeprogression = 100*sum(sum(sum(Label_lake~=0)))/num_voxel;
        cla(axe_video(2).s); % Clear axe
        grid(axe_video(2).s,'on');
        yyaxis(axe_video(2).s,'left');
        x_video_l = all_altitude(1:current_altitude_iteration);
        y_video_l = [y_video_l length(currentvalues)-1];
        plot(axe_video(2).s, x_video_l, y_video_l,'LineStyle','-','Marker','o','MarkerSize',12,'LineWidth',2);
        xlabel(axe_video(2).s,'Altitude (voxel length)');
        xlim(axe_video(2).s,[min(x_video_l) 0]);
        ylabel(axe_video(2).s,'Number of bassins');
        yyaxis(axe_video(2).s,'right');
        y_video_2 = [y_video_2 volumeprogression];
        plot(axe_video(2).s, x_video_l, y_video_2,'LineStyle','-','Marker','x','MarkerSize',12,'LineWidth',2);
        ylabel(axe_video(2).s,'Phase volume (%)');        
        set(axe_video(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title('Progression','FontName','Times New Roman','FontSize',18,'Parent',axe_video(2).s);        
        % Video
        stored_frame(current_altitude_iteration) = getframe(Fig_video);
        writeVideo(video_handle,stored_frame(current_altitude_iteration))        
    end    

end
toc(t_altitudeiteration)


if visualize_2D
    close(video_handle) % Close video   
    
    background = 0;
    edgewithbackground = false;
    
    sub_axes_illustration(3).s=subplot(2,2,3,'Parent',Fig_illustration); % Create axes  
    hold(sub_axes_illustration(3).s,'on'); % Active subplot
    slice_image = image(slice_color,'parent',sub_axes_illustration(3).s); % Display the slice
    axis(sub_axes_illustration(3).s,'tight'); % Fit the axes box
    axis(sub_axes_illustration(3).s,'equal'); % Aspect ratio is 1:1
    % Label
    xlabel(sub_axes_illustration(3).s,'Voxels');
    ylabel(sub_axes_illustration(3).s,'Voxels');
    set(sub_axes_illustration(3).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
    title(['Number of bassins identified: ' num2str(length(currentvalues)-1,'%i')],'FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(3).s);
    box(sub_axes_illustration(3).s,'on')
    hold(sub_axes_illustration(3).s,'off'); % Active subplot    
    
    % Find watershedlines
    [index_border_label,~,~,~] = Function_identify_labelsedges(Label_lake, background, edgewithbackground);
    sub_axes_illustration(4).s=subplot(2,2,4,'Parent',Fig_illustration); % Create axes  
    hold(sub_axes_illustration(4).s,'on'); % Active subplot
    tmp = slice_grey; tmp(index_border_label) = 1;    
    slice_color(:,:,1)=tmp; % Attribute RGB color
    tmp = slice_grey; tmp(index_border_label) = 0;    
    slice_color(:,:,2)=tmp;
    tmp = slice_grey; tmp(index_border_label) = 0;    
    slice_color(:,:,3)=tmp;
    slice_image = image(slice_color,'parent',sub_axes_illustration(4).s); % Display the slice
    axis(sub_axes_illustration(4).s,'tight'); % Fit the axes box
    axis(sub_axes_illustration(4).s,'equal'); % Aspect ratio is 1:1
    % Label
    xlabel(sub_axes_illustration(4).s,'Voxels');
    ylabel(sub_axes_illustration(4).s,'Voxels');
    set(sub_axes_illustration(4).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
    title('Watershed lines','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(4).s);
    box(sub_axes_illustration(4).s,'on')
    hold(sub_axes_illustration(4).s,'off'); % Active subplot       
    function_savefig(Fig_illustration, [desktop '\'], 'Illustration of two-dimensional watershed method')
end

% Get all discrete particle id
unique_particle = unique(Label_lake);
unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
% Get number of particle
number_particle = length(unique_particle);
if details_convergence
    fprintf ('Initial number of discrete particle identified: %i\n',number_particle);
end

%% CORRECT IDENTIFICATION

tic_cpsdcorrection = tic;
chess_pattern_removal = false; % Not required for Watershed method
P1 = []; P2 = []; P3 = []; % Not required for Watershed method
[Label_lake] = Function_correct_DPSD_identification(binary_phase, Label_lake, cpsd_refining, details_convergence, chess_pattern_removal, P1, P2, P3, visualize_2D);
toc(tic_cpsdcorrection);

%% UPSCALING -- VERY SLOW
if downscaling_factor~=1
    tic_upscalingcorrection = tic;
    parameters_scaling.scaling_factor = 1/downscaling_factor;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = 0;
    % Scale
    Label_lake_upscaled = function_scaling(uint16(Label_lake),parameters_scaling); % -- VERY SLOW
    Domain_size = size(Label_lake_upscaled);
    Label_lake = Label_lake_upscaled;

    % Remove values outside binary_phase_initial
    Label_lake(binary_phase_initial==0)=0;

    % Is there some missing values ?
    tmp = zeros(size(Label_lake));
    tmp(Label_lake==0)=1;
    missingpoints = double(binary_phase_initial).*double(tmp);
    id_missingpoints = find(missingpoints);
    if ~isempty(id_missingpoints)>0 % Yes, assign to nearest label lake
        [~,idx] = bwdist(~tmp);
        Label_lake(id_missingpoints) = Label_lake(idx(id_missingpoints));
    end
    if length(Domain_size)==2
        Domain_size = [Domain_size 1];
    end
    toc(tic_upscalingcorrection);
end


%% FINAL DISCRETE PARTICLE SIZE
% Initialisation
tic_size = tic;
D_PSD_particle_size = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% Unique id
unique_particle = unique(Label_lake);
unique_particle(1)=[];
% Number of different particle
number_particle = length(unique_particle);
% loop over all particles
for current_=1:1:number_particle
    % Get id_
    particle_id = unique_particle(current_);
    % Find all voxels of the particle
    index_particle=find(Label_lake==particle_id);
    % Number of voxel
    number_voxel_particle = length(index_particle);
    volume_ = number_voxel_particle;
    area_ = number_voxel_particle;
    % Equivalent diameter size
    if Domain_size(3)>1
        % 3D case
        equivalent_diameter_size= 2 * ((3*volume_  /(4*pi))^(1/3));
    else
        % 2D case
        equivalent_diameter_size= 2 * ((area_/pi)^(1/2));
    end
    D_PSD_particle_size(index_particle)= equivalent_diameter_size;
end
toc(tic_size);

%% VISUALIZATION 
if visualize_2D % Only if 2D
    parameters_figure=[];
    Function_visualize_particle_identification(binary_phase, Label_lake, D_PSD_particle_size, parameters_figure);
end

end