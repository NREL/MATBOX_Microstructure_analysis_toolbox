function [D_PSD_particle_size, Label_lake, Repulsive_matrix, Repulsive_matrix_sign] = Function_Discrete_particle_size_PCRF_algorithm_v4(binary_phase,cpsd_refining,downscaling_factor, erosion_distance, percluster, details_convergence,visualize_2D,display_erosion,exponent_law,max_dist)

%% REFERENCE
% Quantitative Relationships Between Pore Tortuosity, Pore Topology, and Solid Particle Morphology Using a Novel Discrete Particle Size Algorithm
% Francois L. E. Usseglio-Viretta,Donal P. Finegan, Andrew Colclasure, Thomas M. M. Heenan, Daniel Abraham, Paul Shearing, and Kandler Smith,
% Journal of The Electrochemical Society, 2020 167 100513

sz = size(binary_phase);
dimension = length(sz);
if dimension==2
    sz = [sz 1];
end
 

%% DOWNSCALE
if downscaling_factor~=1
    tmp = [sz(1)-100:1:sz(1)];
    tmp2 = mod(tmp,downscaling_factor);
    idx = max(find(tmp2==0));
    cropx = length(tmp)-idx;
    tmp = [sz(2)-100:1:sz(2)];
    tmp2 = mod(tmp,downscaling_factor);
    idx = max(find(tmp2==0));
    cropy = length(tmp)-idx;
    if dimension==2
        binary_phase(end-cropx+1:end, :) = [];
        binary_phase(:, end-cropy+1:end) = [];
    else
        tmp = [sz(3)-100:1:sz(3)];
        tmp2 = mod(tmp,downscaling_factor);
        idx = max(find(tmp2==0));
        cropz = length(tmp)-idx;
        binary_phase(end-cropx+1:end,:,:) = [];
        binary_phase(:,end-cropy+1:end,:) = [];
        binary_phase(:,:,end-cropz+1:end) = [];
    end

    binary_phase_initial = binary_phase;
    parameters_scaling.scaling_factor = downscaling_factor;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = 0;
    % Scale
    binary_phase = function_scaling(binary_phase,parameters_scaling);
    sz = size(binary_phase);
    if dimension==2
        sz = [sz 1];
    end
end
if dimension==2
    zdisplay=1;
else
    zdisplay = round(sz(3)/2);
end


if erosion_distance > 0
    dilatation_distance = 0;
    distance_method = 'Euclidean';
    tolerance = 0.1;
    [binary_phase] = Function_morphologyopening_erosiondilatation_algorithm(binary_phase,erosion_distance,tolerance,dilatation_distance,distance_method); % Fill crack, fill gap btw adjacent particles
    if display_erosion
        figure; imagesc(binary_phase(:,:,zdisplay)); axis equal; colormap copper; pause(0.1);
    end
    % dilatation_distance = erosion_distance;
    % distance_method = 'Euclidean';
    % tolerance = 0.1;
    % [binary_phase] = Function_morphologyopening_erosiondilatation_algorithm(binary_phase,0,tolerance,dilatation_distance,distance_method); % Fill crack, fill gap btw adjacent particles
    % figure
    % imagesc(binary_phase); axis equal; colormap copper; pause(0.1);    
    % [binary_phase] = Function_morphologyopening_erosiondilatation_algorithm(binary_phase,2*erosion_distance,tolerance,0,distance_method); % Accentue bottleneck
    % figure
    % imagesc(binary_phase); axis equal; colormap copper; pause(0.1);
end

%% ANALYSE PER CLUSTER
if percluster
    if dimension==2
        Clusters = bwlabel(binary_phase,4);
    else
        Clusters = bwlabeln(binary_phase,6);
    end
    uni_clusters = unique(Clusters);
    uni_clusters(uni_clusters==0)=[];
    ncluster = length(uni_clusters);
    Label_lake_full = zeros(sz);
    Repulsive_matrix_full = zeros(sz(1),sz(2),sz(3),3);
    Repulsive_matrix_sign_full = zeros(sz(1),sz(2),sz(3),3);
else
    ncluster=1;
end

for kcluster = 1:1:ncluster
    if percluster
        idx=find(Clusters==uni_clusters(kcluster));
        [IX,IY,IZ]=ind2sub(sz,idx);
        x_mincluster = min(IX); x_maxcluster = max(IX);
        y_mincluster = min(IY); y_maxcluster = max(IY);
        if dimension==2
            z_mincluster = 1; z_maxcluster = 1;
            binary_phase = Clusters(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster);
        else
            z_mincluster = min(IZ); z_maxcluster = max(IZ);
            binary_phase = Clusters(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster);
        end
        binary_phase(binary_phase~=uni_clusters(kcluster))=0;
        binary_phase(binary_phase~=0)=1;
    end

    % Add a 0 layer
    if length(binary_phase)>1
        sz_tmp = size(binary_phase);
        if length(sz_tmp)==dimension
            tmp = zeros(sz_tmp+2);
            if dimension==2
                tmp(2:end-1,2:end-1) = binary_phase;
            else
                tmp(2:end-1,2:end-1,2:end-1) = binary_phase;
            end
        else
            tmp = zeros(sz_tmp(1)+2,sz_tmp(2)+2,3);
            tmp(2:end-1,2:end-1,2) = binary_phase(:,:,1);
        end
    else
        if dimension==2
            tmp = zeros(3,3);
            tmp(2,2) = binary_phase;
        else
            tmp = zeros(3,3,3);
            tmp(2,2,2) = binary_phase;
        end
    end
    binary_phase = tmp;

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
        video_handle = VideoWriter([[desktop '\'] 'Illustration_2D_PCRF'],'mpeg-4');
        set(video_handle,'Quality',100); % Set video quality
        set(video_handle,'FrameRate',5); % Set video framerate
        open(video_handle)
        Fig_video = figure; % Create figure
        Fig_video.Name= 'Illustration of two-dimensional PCRF method'; % Figure name
        Fig_video.Color='white'; % Background colour
        set(Fig_video,'position',scrsz); % Full screen figure
        sgtitle(Fig_video,'Animation of two-dimensional PCRF method','FontName','Times New Roman','FontSize',18);
        axe_video(1).s = subplot(1,2,1,'Parent',Fig_video); % Create axes
        axe_video(2).s = subplot(1,2,2,'Parent',Fig_video);
        y_video_l = []; y_video_2 = []; num_voxel = sum(sum(sum(binary_phase==1)));
        current_frame = 0;

        % Prepare illustration 1
        Fig_illustration = figure; % Create figure
        Fig_illustration.Name= 'Illustration of two-dimensional PCRF method'; % Figure name
        Fig_illustration.Color='white'; % Background colour
        set(Fig_illustration,'position',scrsz); % Full screen figure
        sgtitle(Fig_illustration,'Illustration of two-dimensional PCRF method','FontName','Times New Roman','FontSize',18);
        sub_axes_illustration(1).s=subplot(2,3,1,'Parent',Fig_illustration); % Create axes
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

    % % For the investigated phase
    % Number of voxel
    number_voxel_phase = sum(sum(sum(binary_phase==1)));
    % Find index for the phase
    index_phase = find(binary_phase==1);
    % Get all coordinates
    [P1,P2,P3] = ind2sub(Domain_size,index_phase);

    % % For the complementary phase
    % Number of voxel
    number_voxel_complementary = sum(sum(sum(binary_phase==0)));
    % Find index for the phase
    index_complementary_phase = find(binary_phase==0);
    % Get all coordinates
    [C1,C2,C3] = ind2sub(Domain_size,index_complementary_phase);

    % Matrix that contains the indice
    index_array=zeros(Domain_size(1),Domain_size(2),Domain_size(3),2);
    for current_voxel_phase=1:1:number_voxel_phase
        x_=P1(current_voxel_phase);
        y_=P2(current_voxel_phase);
        z_=P3(current_voxel_phase);
        index_array(x_,y_,z_,1)=index_phase(current_voxel_phase);
    end

    %% STEP 1: IDENTIFY THE EDGE
    % % Identify voxel of the complemtary phase that belong at the border
    % They will be marked 2 in the binary phase

    for current_voxel=1:1:number_voxel_complementary
        % Get back voxel coordinate
        x_=C1(current_voxel); y_=C2(current_voxel); z_=C3(current_voxel);

        % Check x-
        if x_>1
            if binary_phase(x_-1,y_,z_)==1
                binary_phase(x_,y_,z_)=2;
            end
        end
        % Check y-
        if y_>1
            if binary_phase(x_,y_-1,z_)==1
                binary_phase(x_,y_,z_)=2;
            end
        end
        % Check z-
        if z_>1
            if binary_phase(x_,y_,z_-1)==1
                binary_phase(x_,y_,z_)=2;
            end
        end
        % Check x+
        if x_<Domain_size(1)
            if binary_phase(x_+1,y_,z_)==1
                binary_phase(x_,y_,z_)=2;
            end
        end
        % Check y+
        if y_<Domain_size(2)
            if binary_phase(x_,y_+1,z_)==1
                binary_phase(x_,y_,z_)=2;
            end
        end
        % Check z+
        if z_<Domain_size(3)
            if binary_phase(x_,y_,z_+1)==1
                binary_phase(x_,y_,z_)=2;
            end
        end

    end
    % Number of voxel at the border
    number_voxel_border = sum(sum(sum(binary_phase==2)));
    % Find index for the border voxel
    index_border_phase = find(binary_phase==2);
    % Get all coordinates
    [B1,B2,B3] = ind2sub(Domain_size,index_border_phase);

    %% STEP 2: EVALUATE REPULSIVE FORCE INDUCED BY EACH VOXEL AT THE BORDER

    % Initialise repulsive force F matrix
    % There are three components F=Fx.ex + Fy.ey + Fz.ez
    Repulsive_matrix = zeros(Domain_size(1),Domain_size(2),Domain_size(3),3);
    % Initialise
    history_test=zeros(2*max_dist,3);

    % Voxel are sorting according to their distance with the current borde voxel
    % data_=zeros(number_voxel_close,11);
    % data_(:,1)=r12;
    % data_(:,2)=P_1; data_(:,3)=P_2; data_(:,4)=P_3;
    % data_(:,5)=n_step;
    % data_(:,6)=dx_; data_(:,7)=dy_; data_(:,8)=dz_;
    % data_(:,9)=dx_12; data_(:,10)=dy_12; data_(:,11)=dz_12;

    if details_convergence
        disp 'Calculation of the repulsive force in progress...(%)';
    end
    next_step_display=1;
    % Loop over all voxel at the edge
    for current_voxel_border=1:1:number_voxel_border
        % Progress of the calculation in percents
        progress_force_calculation = 100*current_voxel_border/number_voxel_border;
        if progress_force_calculation>=next_step_display && details_convergence
            fprintf (' %4.1f',progress_force_calculation);
            next_step_display=next_step_display+1;
        end
        % Get back voxel coordinate
        x_b=B1(current_voxel_border); y_b=B2(current_voxel_border); z_b=B3(current_voxel_border);

        % Reduce the index only to the voxel close to this border
        % Range of the analysis
        x_min = max(1,x_b-max_dist); y_min = max(1,y_b-max_dist); z_min = max(1,z_b-max_dist);
        x_max = min(Domain_size(1),x_b+max_dist); y_max = min(Domain_size(2),y_b+max_dist); z_max = min(Domain_size(3),z_b+max_dist);
        % Crop the index array
        index_close = index_array(x_min:x_max,y_min:y_max,z_min:z_max,1);
        % Unique index
        unique_index = unique(index_close);
        % Remove 0
        unique_index(1)=[];
        % Number of voxel
        number_voxel_close = length(unique_index);

        % Data will be stored here
        data_=zeros( number_voxel_close, 11);

        % Get the voxel coordinates
        [data_(:,2),data_(:,3),data_(:,4)] = ind2sub(Domain_size,unique_index);

        % Convert in our metrics
        x_b=x_b-0.5;y_b=y_b-0.5;z_b=z_b-0.5;

        % Convert in our metrics and evaluate distance along each direction
        data_(:,9) = data_(:,2)-0.5-x_b;
        data_(:,10) = data_(:,3)-0.5-y_b;
        data_(:,11) = data_(:,4)-0.5-z_b;
        % Evaluate euclidean distance
        data_(:,1) = (data_(:,9).^2+data_(:,10).^2+data_(:,11).^2).^(0.5);

        % Number of step required to test voxels between them
        data_(:,5) = floor(data_(:,1));

        % Delta position
        data_(:,6) = data_(:,9)./data_(:,5);
        data_(:,7) = data_(:,10)./data_(:,5);
        data_(:,8) = data_(:,11)./data_(:,5);

        % Sort with r12 deacreasing
        data_ = sortrows(data_,-1);

        % For each voxel of the phase
        for current_voxel_phase=1:1:number_voxel_close

            % Get back voxel coordinate
            x_pp=data_(current_voxel_phase,2); y_pp=data_(current_voxel_phase,3); z_pp=data_(current_voxel_phase,4);

            % Voxel whom values have been already calculated (as being between
            % a previous tested voxel and the voxel border) are not recalculated
            if index_array(x_pp,y_pp,z_pp,2)~=current_voxel_border

                % Test loop
                test_=1; iteration_test=0;

                while (test_==1 && (iteration_test+1)<=data_(current_voxel_phase,5))
                    iteration_test=iteration_test+1;
                    % Tested position
                    x_tested = x_b+iteration_test*data_(current_voxel_phase,6);
                    y_tested = y_b+iteration_test*data_(current_voxel_phase,7);
                    z_tested = z_b+iteration_test*data_(current_voxel_phase,8);
                    % Converted in voxels coordinates
                    x_tested=floor(x_tested)+1;
                    y_tested=floor(y_tested)+1;
                    z_tested=floor(z_tested)+1;
                    % Check voxel
                    test_ = binary_phase(x_tested,y_tested,z_tested);
                    % Conserve history of tested voxels
                    history_test(iteration_test,1)=x_tested;
                    history_test(iteration_test,2)=y_tested;
                    history_test(iteration_test,3)=z_tested;
                end

                if test_==1
                    % % Clear line of view between the voxel exists
                    force_=1/(data_(current_voxel_phase,1)^exponent_law); % total force >0
                    force_1 = force_ *data_(current_voxel_phase,9)/data_(current_voxel_phase,1); % Component along direction 1, with the sign
                    force_2 = force_ *data_(current_voxel_phase,10)/data_(current_voxel_phase,1); % Component along direction 2, with the sign
                    force_3 = force_ *data_(current_voxel_phase,11)/data_(current_voxel_phase,1); % Component along direction 3, with the sign
                    Repulsive_matrix(x_pp,y_pp,z_pp,1) = Repulsive_matrix(x_pp,y_pp,z_pp,1) + force_1;
                    Repulsive_matrix(x_pp,y_pp,z_pp,2) = Repulsive_matrix(x_pp,y_pp,z_pp,2) + force_2;
                    Repulsive_matrix(x_pp,y_pp,z_pp,3) = Repulsive_matrix(x_pp,y_pp,z_pp,3) + force_3;
                    index_array(x_pp, y_pp, z_pp,2)=current_voxel_border;
                    % % intermediate voxels
                    % Get intermediate coordinates
                    x_intermediate = history_test(1:data_(current_voxel_phase,5),1);
                    y_intermediate = history_test(1:data_(current_voxel_phase,5),2);
                    z_intermediate = history_test(1:data_(current_voxel_phase,5),3);

                    % 1 9 10 11
                    v9 = x_intermediate-0.5-x_b; v10 = y_intermediate-0.5-y_b; v11 = z_intermediate-0.5-z_b;
                    v1 = (v9.^2+v10.^2+v11.^2).^(0.5);
                    for marker_=1:1:data_(current_voxel_phase,5)
                        % Get intermediate coordinates
                        x_int = x_intermediate(marker_,1);
                        y_int = y_intermediate(marker_,1);
                        z_int = z_intermediate(marker_,1);
                        if index_array(x_int, y_int, z_int,2)~=current_voxel_border
                            % Force
                            force_=1/(v1(marker_,1)^exponent_law); % total force >0
                            force_1 = force_ *v9(marker_,1)/v1(marker_,1); % Component along direction 1, with the sign
                            force_2 = force_ *v10(marker_,1)/v1(marker_,1); % Component along direction 2, with the sign
                            force_3 = force_ *v11(marker_,1)/v1(marker_,1); % Component along direction 3, with the sign
                            % Update force
                            Repulsive_matrix(x_int,y_int,z_int,1) = Repulsive_matrix(x_int,y_int,z_int,1) + force_1;
                            Repulsive_matrix(x_int,y_int,z_int,2) = Repulsive_matrix(x_int,y_int,z_int,2) + force_2;
                            Repulsive_matrix(x_int,y_int,z_int,3) = Repulsive_matrix(x_int,y_int,z_int,3) + force_3;
                            % Mark these voxels as calculated
                            index_array(x_int, y_int, z_int,2)=current_voxel_border;
                        end
                    end

                end
            end
        end

    end
    if details_convergence
        fprintf (' End of calculation\n');
    end

    Repulsive_matrix_sign = Repulsive_matrix;
    Repulsive_matrix_sign(:,:,:,1) = sign(Repulsive_matrix(:,:,:,1));
    Repulsive_matrix_sign(:,:,:,2) = sign(Repulsive_matrix(:,:,:,2));
    Repulsive_matrix_sign(:,:,:,3) = sign(Repulsive_matrix(:,:,:,3));

    if visualize_2D
        sub_axes_illustration(2).s=subplot(2,3,2,'Parent',Fig_illustration); % Create axes
        hold(sub_axes_illustration(2).s,'on'); % Active subplot
        tmp = Repulsive_matrix(:,:,1,1); tmp(binary_phase==0)=NaN;
        imagesc(sub_axes_illustration(2).s,tmp);
        custom_map = jet; custom_map(end,:)=[1 1 1];
        colormap(sub_axes_illustration(2).s,custom_map);
        axis(sub_axes_illustration(2).s,'tight'); % Fit the axes box
        axis(sub_axes_illustration(2).s,'equal'); % Aspect ratio is 1:1
        colorbar
        % Label
        xlabel('Voxels');
        ylabel('Voxels');
        set(sub_axes_illustration(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title('F(x)','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(2).s);
        box(sub_axes_illustration(2).s,'on')
        hold(sub_axes_illustration(2).s,'off'); % Active subplot

        sub_axes_illustration(3).s=subplot(2,3,3,'Parent',Fig_illustration); % Create axes
        hold(sub_axes_illustration(3).s,'on'); % Active subplot
        tmp = Repulsive_matrix(:,:,1,2); tmp(binary_phase==0)=NaN;
        imagesc(sub_axes_illustration(3).s,tmp);
        custom_map = jet; custom_map(end,:)=[1 1 1];
        colormap(sub_axes_illustration(3).s,custom_map);
        axis(sub_axes_illustration(3).s,'tight'); % Fit the axes box
        axis(sub_axes_illustration(3).s,'equal'); % Aspect ratio is 1:1
        colorbar
        % Label
        xlabel('Voxels');
        ylabel('Voxels');
        set(sub_axes_illustration(3).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title('F(y)','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(3).s);
        box(sub_axes_illustration(3).s,'on')
        hold(sub_axes_illustration(3).s,'off'); % Active subplot

        sub_axes_illustration(3).s=subplot(2,3,5,'Parent',Fig_illustration); % Create axes
        hold(sub_axes_illustration(3).s,'on'); % Active subplot
        tmp = sign(Repulsive_matrix(:,:,1,1)); tmp(binary_phase==0)=0;
        imagesc(sub_axes_illustration(3).s,tmp);
        custom_map = [0 0 1; 0 0 0; 1 0 0];
        colormap(sub_axes_illustration(3).s,custom_map);
        axis(sub_axes_illustration(3).s,'tight'); % Fit the axes box
        axis(sub_axes_illustration(3).s,'equal'); % Aspect ratio is 1:1
        colorbar
        % Label
        xlabel('Voxels');
        ylabel('Voxels');
        set(sub_axes_illustration(3).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title('sign(F(x))','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(3).s);
        box(sub_axes_illustration(3).s,'on')
        hold(sub_axes_illustration(3).s,'off'); % Active subplot

        sub_axes_illustration(4).s=subplot(2,3,6,'Parent',Fig_illustration); % Create axes
        hold(sub_axes_illustration(4).s,'on'); % Active subplot
        tmp = sign(Repulsive_matrix(:,:,1,2)); tmp(binary_phase==0)=0;
        imagesc(sub_axes_illustration(4).s,tmp);
        custom_map = [0 0 1; 0 0 0; 1 0 0];
        colormap(sub_axes_illustration(4).s,custom_map);
        axis(sub_axes_illustration(4).s,'tight'); % Fit the axes box
        axis(sub_axes_illustration(4).s,'equal'); % Aspect ratio is 1:1
        colorbar
        % Label
        xlabel('Voxels');
        ylabel('Voxels');
        set(sub_axes_illustration(4).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title('sign(F(y))','FontName','Times New Roman','FontSize',14,'Parent',sub_axes_illustration(4).s);
        box(sub_axes_illustration(4).s,'on')
        hold(sub_axes_illustration(4).s,'off'); % Active subplot
        function_savefig(Fig_illustration,desktop,'\Field')
    end


    %% STEP3: DROP CHARGED PARTICLE WITHIN THE PHASE AND CALCULATE THEIR MOVEMENT

    % What is called 'seed' in this program:
    % - it is simply the point where the particle will stop its movement
    % - it is not used later as a starting point for the particle growing:
    % there is no particle growing in this algorithm
    % - actual role of seed in this algorithm: distinct partcicle whom seed is adjecent are merged later
    % - it is used as an intermediary result, that role will be later replaced
    % with the center of mass of the identified particle

    % Discrete particle give the id number of the discrete particle.

    % The path used by a particle define a particle
    % If the path of a particle meet the path of a previous particle, then both
    % path define the same particle.

    % Initialisation
    Label_lake = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    seed_particle = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    % Discret particle id numerotation
    particle_id=0;
    % Loop over all voxels
    for current_voxel=1:1:number_voxel_phase

        % Get back voxel coordinate
        x_=P1(current_voxel); y_=P2(current_voxel); z_=P3(current_voxel);
        % Check if the voxel has been already attributed to a particle
        if Label_lake(x_,y_,z_)==0
            % Not yet attributed!
            % Temporary particle numerotation
            particle_id_tmp = particle_id+1;
            Label_lake(x_,y_,z_)=particle_id_tmp;
            % Store path
            path_particle = [x_,y_,z_];
            % Determine next move
            movement=1;
            while movement==1
                % Evaluate local force repulsive
                % Column1: direction
                % Column2: force value
                % Column3: sign
                local_force = [1 abs(Repulsive_matrix(x_,y_,z_,1)) sign(Repulsive_matrix(x_,y_,z_,1));2 abs(Repulsive_matrix(x_,y_,z_,2)) sign(Repulsive_matrix(x_,y_,z_,2));3 abs(Repulsive_matrix(x_,y_,z_,3)) sign(Repulsive_matrix(x_,y_,z_,3));];
                % Sort according to force value
                local_force = sortrows(local_force,-2);
                for direction=1:1:3
                    % Direction
                    if local_force(direction,1)==1
                        diff_=[1 0 0];
                    elseif local_force(direction,1)==2
                        diff_=[0 1 0];
                    elseif local_force(direction,1)==3
                        diff_=[0 0 1];
                    end
                    % Sens
                    diff_=diff_*local_force(direction,3);
                    % New position
                    new_ = [x_,y_,z_]+diff_;
                    % Check if new position is within the domain
                    if new_(1)>=1 && new_(2)>=1 && new_(3)>=1 && new_(1)<=Domain_size(1) && new_(2)<=Domain_size(2) && new_(3)<=Domain_size(3)
                        % Check if new position belongs to the phase
                        if binary_phase(new_(1),new_(2),new_(3))==1
                            % Actualise position
                            x_=new_(1); y_=new_(2); z_=new_(3);
                        end
                    end
                end
                % Check if new position has been already attributed
                if Label_lake(x_,y_,z_)==0
                    % Not attributed!
                    Label_lake(x_,y_,z_)=particle_id_tmp;
                    % Store path
                    path_particle = [path_particle; [x_,y_,z_]];
                    %disp 'New';
                elseif Label_lake(x_,y_,z_)==particle_id_tmp
                    % We back to one of out our previous position, or we didnt move
                    % We have identified a new particle!
                    movement=0;
                    new_particle=1;
                    particle_id=particle_id+1;
                    seed_particle(x_,y_,z_)=particle_id;
                    % Store path
                    path_particle = [path_particle; [x_,y_,z_]];
                    % The seed particle is defined as the two last position
                    % (because the drop moves vice versa from one to another)
                    n_path=size(path_particle,1);
                    if n_path>=2
                        x_previous = path_particle(n_path-1,1);
                        y_previous = path_particle(n_path-1,2);
                        z_previous = path_particle(n_path-1,3);
                        if seed_particle(x_previous,y_previous,z_previous)==0
                            seed_particle(x_previous,y_previous,z_previous)=particle_id;
                        end
                    end
                    %disp 'Backtrack or no movement';
                else
                    % We have back to an particle already identified
                    movement=0;
                    % We did not have identified a new particle!
                    new_particle=0;
                    %disp 'Already identified';
                end



                if visualize_2D % Uncomment to see step by step. Although, visualization is very slow.
                    %                 [length_trajectory,~] = size(path_particle);
                    %
                    %                 % Initializaion
                    %                 slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
                    %                 slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
                    %                 slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
                    %                 slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
                    %                 slice_grey = zeros(Domain_size(1),Domain_size(2)); % Blue color map
                    %                 currentvalues = unique(Label_lake);
                    %                 for k=1:1:length(currentvalues)
                    %                     idx=find(Label_lake==currentvalues(k));
                    %                     if isempty( intersect(existingvalues,currentvalues(k)) ) && currentvalues(k)~=0
                    %                         idx_source=[idx_source; idx];
                    %                         existingvalues=unique([existingvalues currentvalues(k)]);
                    %                     end
                    %                     if currentvalues(k)==0
                    %                         slice_r(idx) = 1; slice_g(idx) = 1; slice_b(idx) = 1; % Background
                    %                         slice_grey(idx) = 1;
                    %                     else
                    %                         slice_r(idx) = lake_id_RGB_color(k,1); slice_g(idx) = lake_id_RGB_color(k,2); slice_b(idx) = lake_id_RGB_color(k,3); % Background
                    %                         slice_grey(idx) = lake_id_grey_color(k);
                    %                     end
                    %                 end
                    %                 % Border
                    %                 slice_r(index_border_phase) = 0; slice_g(index_border_phase) = 0; slice_b(index_border_phase) = 0;
                    %                 slice_grey(index_border_phase) = 0;
                    %                 % Source
                    %                 %slice_r(idx_source) = 1; slice_g(idx_source) = 0; slice_b(idx_source) = 0;
                    %                 slice_color(:,:,1)=slice_r; % Attribute RGB color
                    %                 slice_color(:,:,2)=slice_g;
                    %                 slice_color(:,:,3)=slice_b;
                    %
                    %                 hold(axe_video(1).s,'on');
                    %                 slice_image = image(slice_color,'parent',axe_video(1).s); % Display the slice
                    %                 tmp = sign(Repulsive_matrix(:,:,1,1)); tmp(binary_phase==0)=0; tmp(Label_lake~=0)=0;
                    %                 imagesc(axe_video(1).s,tmp, 'AlphaData', .25);
                    %                 colormap(axe_video(1).s,custom_map);
                    %                 axis(axe_video(1).s,'tight'); % Fit the axes box
                    %                 axis(axe_video(1).s,'equal'); % Aspect ratio is 1:1
                    %                 % Label
                    %                 xlabel(axe_video(1).s,'Voxels');
                    %                 ylabel(axe_video(1).s,'Voxels');
                    %                 set(axe_video(1).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
                    %                 title({['Current trajectory length: ' num2str(length_trajectory)],['Number of particles identified: ' num2str(length(currentvalues)-1,'%i')]},'FontName','Times New Roman','FontSize',18,'Parent',axe_video(1).s);
                    %                 set(axe_video(1).s,'YDir','normal')
                    %                 hold(axe_video(1).s,'off');
                    %
                    %                 hold(axe_video(2).s,'on');
                    %                 slice_image = image(slice_color,'parent',axe_video(2).s); % Display the slice
                    %                 tmp = sign(Repulsive_matrix(:,:,1,2)); tmp(binary_phase==0)=0; tmp(Label_lake~=0)=0;
                    %                 imagesc(axe_video(2).s,tmp, 'AlphaData', .25);
                    %                 colormap(axe_video(2).s,custom_map);
                    %                 axis(axe_video(2).s,'tight'); % Fit the axes box
                    %                 axis(axe_video(2).s,'equal'); % Aspect ratio is 1:1
                    %                 % Label
                    %                 xlabel(axe_video(2).s,'Voxels');
                    %                 ylabel(axe_video(2).s,'Voxels');
                    %                 set(axe_video(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
                    %                 title({['Current trajectory length: ' num2str(length_trajectory)],['Number of particles identified: ' num2str(length(currentvalues)-1,'%i')]},'FontName','Times New Roman','FontSize',18,'Parent',axe_video(2).s);
                    %                 set(axe_video(2).s,'YDir','normal')
                    %                 hold(axe_video(2).s,'off');
                    %
                    %                 % Video
                    %                 current_frame = current_frame+1;
                    %                 stored_frame(current_frame) = getframe(Fig_video);
                    %                 writeVideo(video_handle,stored_frame(current_frame))
                end

            end
            if new_particle==0
                % Remap the path with the founded particle
                % Particle id
                particle_found = Label_lake(x_,y_,z_);
                % Number of voxel
                n_voxel=size(path_particle,1);
                for k=1:1:n_voxel
                    Label_lake(path_particle(k,1),path_particle(k,2),path_particle(k,3))=particle_found;
                end
            end
        end
    end
    % Get all discrete particle id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
    % Get number of particle
    number_particle = length(unique_particle);
    if details_convergence==1
        fprintf ('Initial number of discrete particle identified: %i\n',number_particle);
    end

    %% STEP4: ADJACENTS SEEDS ARE MERGED

    % The merging range is limited to one cube of 3*3*3 voxel lenght centered on the investigated seed
    % If the user want it, it can change the merging range, but then, it is not
    % more a non-parametric method.
    merge_range=1;
    % Position of all seeds
    % Find index of all seeds
    index_seed = find(seed_particle~=0);
    % Get their coordinates
    [S1,S2,S3] = ind2sub(Domain_size,index_seed);
    % Number of all seeds
    n_seeds = length(index_seed);
    % Loop over all seeds
    for current_seed=1:1:n_seeds
        % Get back voxel coordinate
        x_=S1(current_seed); y_=S2(current_seed); z_=S3(current_seed);
        % Seed id
        seed_id = seed_particle(x_,y_,z_);
        % Initialize position
        pos_x=x_; pos_y=y_; pos_z=z_;
        for current_range=1:1:merge_range
            % Determine  all possible positions within the domains
            if x_>=current_range+1
                pos_x=[x_-current_range:x_];
            end
            if y_>=current_range+1
                pos_y=[y_-1:y_];
            end
            if z_>=current_range+1
                pos_z=[z_-1:z_];
            end
            if x_<=Domain_size(1)-current_range
                pos_x=[pos_x x_+1];
            end
            if y_<=Domain_size(2)-current_range
                pos_y=[pos_y y_+1];
            end
            if z_<=Domain_size(3)-current_range
                pos_z=[pos_z z_+1];
            end
        end
        % Check each of these positions
        for ix = 1:1:length(pos_x)
            x2 = pos_x(ix);
            for iy = 1:1:length(pos_y)
                y2 = pos_y(iy);
                for iz = 1:1:length(pos_z)
                    z2 = pos_z(iz);
                    if (seed_particle(x2,y2,z2)~=0 && seed_particle(x2,y2,z2)~=seed_id) % Check within investigated phase and different particle
                        % Then replace particle id
                        adjacent_id = seed_particle(x2,y2,z2);
                        index_adjacent_particle = find(seed_particle==adjacent_id);
                        seed_particle(index_adjacent_particle)=seed_id;
                        index_adjacent_particle2 = find(Label_lake==adjacent_id);
                        Label_lake(index_adjacent_particle2)=seed_id;
                    end
                end
            end
        end
    end
    % Re-order seed id
    unique_seed = unique(seed_particle);
    number_seed = length(unique_seed);
    for seed_id=2:1:number_seed % 0 is the complementary phase
        id_to_be_replaced = unique_seed(seed_id);
        index_ = find(seed_particle==id_to_be_replaced);
        seed_particle(index_)=seed_id-1;
        index_ = find(Label_lake==id_to_be_replaced);
        Label_lake(index_)=seed_id-1;
    end
    % Get all discrete particle id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
    % Get number of particle
    number_particle = length(unique_particle);
    if details_convergence==1
        fprintf ('Adjacent seeds have been merged: number of discrete particle: %i\n',number_particle);
    end

    %% TEST
    unique_particle_1 = unique(Label_lake);
    unique_particle_2 = unique(seed_particle);
    if unique_particle_1~=unique_particle_2
        disp 'Error on particle number'
    end
    index_seed = find(seed_particle~=0);
    for current_index=1:1:length(index_seed)
        id_seed = seed_particle(index_seed(current_index));
        id_discrete = Label_lake(index_seed(current_index));
        if id_seed~=id_discrete
            disp 'Error on particle id'
        end
    end
   

    %% CORRECT IDENTIFICATION
    chess_pattern_removal = true; % Required for PCRF method
    tmp_binary_phase = zeros(size(binary_phase));
    tmp_binary_phase(binary_phase==1)=1;
    [Label_lake] = Function_correct_DPSD_identification(tmp_binary_phase, Label_lake, cpsd_refining, details_convergence, chess_pattern_removal, P1, P2, P3, visualize_2D);
    
    %% REMOVE THIN LAYER
    if dimension == 2
        Label_lake = Label_lake(2:end-1,2:end-1);
        Repulsive_matrix = Repulsive_matrix(2:end-1,2:end-1,1,:);
        Repulsive_matrix_sign = Repulsive_matrix_sign(2:end-1,2:end-1,1,:);
    else
        Label_lake = Label_lake(2:end-1,2:end-1,2:end-1);
        Repulsive_matrix = Repulsive_matrix(2:end-1,2:end-1,1,2:end-1,1,:);
        Repulsive_matrix_sign = Repulsive_matrix_sign(2:end-1,2:end-1,1,2:end-1,1,:);
    end

    %% PER CLUSTER
    if percluster
        m = max(max(max(Label_lake_full)));
        Label_lake = Label_lake + m;
        Label_lake(Label_lake==m)=0;

        Label_lake_full(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster) = Label_lake_full(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster) + Label_lake;
        Repulsive_matrix_full(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster,:) = Repulsive_matrix_full(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster,:) + Repulsive_matrix;
        Repulsive_matrix_sign_full(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster,:) = Repulsive_matrix_sign_full(x_mincluster:x_maxcluster,y_mincluster:y_maxcluster,z_mincluster:z_maxcluster,:) + Repulsive_matrix_sign;

        % % Randomize cluster
        % [Label_lake_full,n] = randomize_labels(Label_lake_full);
        % % Figure
        % tmp=randi(255,n+1,3)/255;
        % tmp(1,:)=[1.0 1.0 1.0];
        % figure; imagesc(Label_lake_full); axis equal; axis tight; colormap(tmp);
        % pause(0.1);
    end
    
end

if percluster % Remove temporary arrays
    Label_lake = Label_lake_full;
    clear Label_lake_full
    Repulsive_matrix = Repulsive_matrix_full;
    clear Repulsive_matrix_full
    Repulsive_matrix_sign = Repulsive_matrix_sign_full;
    clear Repulsive_matrix_sign_full    
end

%% UPSCALING
if downscaling_factor~=1
    % Randomize cluster
    [Label_lake_r,n] = randomize_labels(Label_lake);
    % Figure
    tmp=randi(255,n+1,3)/255;
    tmp(1,:)=[1.0 1.0 1.0];
    if display_erosion
        figure; imagesc(Label_lake_r(:,:,zdisplay)); axis equal; axis tight; colormap(tmp); pause(0.1);
    end

    parameters_scaling.scaling_factor = 1/downscaling_factor;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = 0;
    % Scale
    Label_lake_upscaled = function_scaling(uint16(Label_lake),parameters_scaling);
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

    chess_pattern_removal = false;
    % Find index for the phase
    index_phase = find(binary_phase_initial==1);
    % Get all coordinates
    [P1,P2,P3] = ind2sub(Domain_size,index_phase);
    [Label_lake] = Function_correct_DPSD_identification(binary_phase_initial, Label_lake, cpsd_refining, details_convergence, chess_pattern_removal, P1, P2, P3, visualize_2D);
end

%% FINAL DISCRETE PARTICLE SIZE
if length(Domain_size)==2
    Domain_size = [Domain_size 1];
end
% Initialisation
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

%% VISUALIZATION
if visualize_2D % Only if 2D
    parameters_figure=[];
    Function_visualize_particle_identification(binary_phase, Label_lake, D_PSD_particle_size, parameters_figure);
end

%% FUNCTIONS
function [Crnd,n] = randomize_labels(C)
uni = unique(C);
uni(uni==0)=[]; % remove background
n = length(uni);
rnd = randperm(n); % random permutation of the integers from 1 to n without repeating elements.
Crnd = zeros(size(C));
for k=1:1:n
    Crnd(C==uni(k)) = rnd(k);
end
end

end

