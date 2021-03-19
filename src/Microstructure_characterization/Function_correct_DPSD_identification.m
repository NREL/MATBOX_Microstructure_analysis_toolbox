function [Label_lake] = Function_correct_DPSD_identification(binary_phase, Label_lake, cpsd_refining, details_convergence, chess_pattern_removal, P1, P2, P3, visualize_2D)

% ITERATIVE PROCESS - NON SPECIFIC TO THE WATERSHED-IMMERSION METHOD

% If chess_pattern_removal true
% STEP 1a: remove chess pattern usually induced by the method
% STEP 1b: remove the remaining chess pattern not replaced in step 1

% Commom for both chess_pattern_removal true or false
% STEP 1c: search for non-continguous discrete particle and correct them
% STEP 2: treat 1 voxel size particle
% STEP 3: calculate D-PSD size
% STEP 4: calculate C-PSD size from another algorithm (will be done only for the first iteration)
% STEP 5: treat voxel whom D-PSD size < C-PSD size

% Initialize iteration
iteration_=0;
number_change=1e9;

% Domain size of the microstructure
Domain_size=size(binary_phase);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end

% Initialise Discrete_particle checksum and list
% Checksum_state0_Label_lake = DataHash(Label_lake);
% That may induce an OUT OF MEMORY error if Label_lake is too large
% Instead, we are calculating the hash slice by slice
Checksum_state0_Label_lake=[];
for current_=1:1:Domain_size(3)
    current_checksum = DataHash(Label_lake(:,:,current_));
    Checksum_state0_Label_lake=[Checksum_state0_Label_lake current_checksum];
end
Checksum_list.state(1).hash=Checksum_state0_Label_lake;

% Checksum (hash) on Label_lake is used to verify that the iterative process does not bring us back to a previous state and run indefinitely
% The function used to check the Hask is DataHash, from Jan Simon, downloaded from Matlab File Exchange (Copyright (c) 2016, Jan Simon, All rights reserved)
% Copyright notice, list of conditions, and disclaimer are in the third party licences folder of the repo

if visualize_2D && cpsd_refining
    lake_id_RGB_color = 0.1 + (0.9-0.1).*rand(1e6,3);
    lake_id_grey_color = sum(lake_id_RGB_color,2)/3;
    % Prepare Video
    scrsz = get(0,'ScreenSize'); % Screen resolution
    desktop=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
    video_handle = VideoWriter([[desktop '\'] 'Illustration_2D_correction_DPSD'],'mpeg-4');
    set(video_handle,'Quality',100); % Set video quality
    set(video_handle,'FrameRate',1); % Set video framerate
    open(video_handle)
    Fig_video = figure; % Create figure
    Fig_video.Name= 'Illustration of two-dimensional particle identification correction'; % Figure name
    Fig_video.Color='white'; % Background colour
    set(Fig_video,'position',scrsz); % Full screen figure
    hsg = sgtitle(Fig_video,'Animation of two-dimensional particle identification correction','FontName','Times New Roman','FontSize',18);
    axe_video(1).s = subplot(2,3,1,'Parent',Fig_video); % Create axes
    axe_video(2).s = subplot(2,3,2,'Parent',Fig_video);
    axe_video(3).s = subplot(2,3,3,'Parent',Fig_video);
    axe_video(4).s = subplot(2,3,4,'Parent',Fig_video);
    axe_video(5).s = subplot(2,3,5,'Parent',Fig_video);
    custom_jet_colormap = jet; custom_jet_colormap(1,:)=[1 1 1];
    colormap(axe_video(1).s,custom_jet_colormap);
    colormap(axe_video(2).s,custom_jet_colormap);
    
    % Prepare illustration
    Fig_illustration = figure; % Create figure
    Fig_illustration.Name= 'Illustration of two-dimensional particle identification correction'; % Figure name
    Fig_illustration.Color='white'; % Background colour
    set(Fig_illustration,'position',scrsz); % Full screen figure
    sgtitle(Fig_illustration,'Illustration of two-dimensional particle identification correction','FontName','Times New Roman','FontSize',18);
    axe_illustration(1).s = subplot(1,3,1,'Parent',Fig_illustration); % Create axes
    axe_illustration(2).s = subplot(1,3,2,'Parent',Fig_illustration);
    axe_illustration(3).s = subplot(1,3,3,'Parent',Fig_illustration);
    all_idx_smaller_than_cpsd=[];
    
%     if chess_pattern_removal
%         % Prepare illustration
%         Fig_chess = figure; % Create figure
%         Fig_chess.Name= 'Illustration of two-dimensional chess pattern removal'; % Figure name
%         Fig_chess.Color='white'; % Background colour
%         set(Fig_chess,'position',scrsz); % Full screen figure
%         sgtitle(Fig_chess,'Illustration of two-dimensional chess pattern removal','FontName','Times New Roman','FontSize',18);
%         axe_chess(1).s = subplot(1,2,1,'Parent',Fig_chess); % Create axes
%         axe_chess(2).s = subplot(1,2,2,'Parent',Fig_chess);
%         all_idx_chess=[];
%     end
    
end

while number_change~=0
    iteration_=iteration_+1;
    if details_convergence==1
        fprintf ('Current iteration: %i\n',iteration_);
    end
    if visualize_2D && cpsd_refining
        set(hsg,'String',{'Animation of two-dimensional particle identification correction',['Iteration ' num2str(iteration_,'%i')]});
        
        % Edge detection
        background = 0;
        edgewithbackground = false;
        
        slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
        slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
        slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
        slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
        slice_grey = zeros(Domain_size(1),Domain_size(2)); % Blue color map
        currentvalues = unique(Label_lake);
        for k=1:1:length(currentvalues)
            idx=find(Label_lake==currentvalues(k));
            if currentvalues(k)==0
                slice_r(idx) = 1; slice_g(idx) = 1; slice_b(idx) = 1; % Background
                slice_grey(idx) = 1;
            else
                slice_r(idx) = lake_id_RGB_color(k,1); slice_g(idx) = lake_id_RGB_color(k,2); slice_b(idx) = lake_id_RGB_color(k,3); % Background
                slice_grey(idx) = lake_id_grey_color(k);
            end
        end
        % Edges
        [index_border_phase,~,~,~] = Function_identify_edges(binary_phase);
        slice_r(index_border_phase) = 0; slice_g(index_border_phase) = 0; slice_b(index_border_phase) = 0;
        slice_grey(index_border_phase) = 0;
        % Watershed lines
        [index_border_label,~,~,~] = Function_identify_labelsedges(Label_lake, background, edgewithbackground);
        tmp = slice_grey; tmp(index_border_label) = 1;
        slice_color(:,:,1)=tmp; % Attribute RGB color
        tmp = slice_grey; tmp(index_border_label) = 0;
        slice_color(:,:,2)=tmp;
        tmp = slice_grey; tmp(index_border_label) = 0;
        slice_color(:,:,3)=tmp;
        slice_image = image(slice_color,'parent',axe_video(4).s); % Display the slice
        axis(axe_video(4).s,'tight'); % Fit the axes box
        axis(axe_video(4).s,'equal'); % Aspect ratio is 1:1
        set(axe_video(4).s,'YDir','normal')
        % Label
        xlabel(axe_video(4).s,'Voxels');
        ylabel(axe_video(4).s,'Voxels');
        set(axe_video(4).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title('Particle identification with watershed lines (before)','FontName','Times New Roman','FontSize',14,'Parent',axe_video(4).s);
        box(axe_video(4).s,'on')
        hold(axe_video(4).s,'off'); % Active subplot
        if iteration_==1
            slice_image = image(slice_color,'parent',axe_illustration(1).s); % Display the slice
            axis(axe_illustration(1).s,'tight'); % Fit the axes box
            axis(axe_illustration(1).s,'equal'); % Aspect ratio is 1:1
            set(axe_illustration(1).s,'YDir','normal')
            % Label
            xlabel(axe_illustration(1).s,'Voxels');
            ylabel(axe_illustration(1).s,'Voxels');
            set(axe_illustration(1).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
            title('Particle identification with watershed lines (before)','FontName','Times New Roman','FontSize',14,'Parent',axe_illustration(1).s);
            box(axe_illustration(1).s,'on')
            hold(axe_illustration(1).s,'off'); % Active subplot
        end
    end
    
    if chess_pattern_removal
        number_voxel_phase = sum(sum(sum(binary_phase==1)));
        %% STEP1: DISCRETE PARTICLES ARE ANALYSED. FOR WHICH A CHESS PATTERN IS IDENTIFIED, THEY ARE MEGRED WITH THE NEAREST PARTICLE THAT DOES NOT EXIBIT THIS PATTERN
        
        % Chess pattern definition (x+,x-,y+,y-,z+,z- must be different with the value in x,y,z)
        %      ~A
        %   ~A  A ~A -> voxel A belongs to the chess
        %      ~A
        
        % This is an iterative process
        % It will go on until there are no more chess pattern or no more change are induced
        sav_index_chess=[];
        number_voxel_chess = 1;
        iteration_chess=0;
        while number_voxel_chess~=0
            iteration_chess=iteration_chess+1;
            % Initialise chess pattern
            chess_pattern = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
            if Domain_size(3)==1
                for current_voxel_phase=1:1:number_voxel_phase
                    % Get back voxel coordinate
                    x_p=P1(current_voxel_phase); y_p=P2(current_voxel_phase); z_p=P3(current_voxel_phase);
                    % Chess pattern analysis only for voxels on the 'subdomains' Domain_size-1
                    if x_p>1 && y_p>1 && x_p<Domain_size(1) && y_p<Domain_size(2)
                        % Check chess pattern
                        id_= Label_lake(x_p,y_p,z_p);
                        % Along plane 12
                        if (Label_lake(x_p+1,y_p,z_p)~=id_ && Label_lake(x_p-1,y_p,z_p)~=id_ && Label_lake(x_p,y_p+1,z_p)~=id_ && Label_lake(x_p,y_p-1,z_p)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                    end
                end
            else
                for current_voxel_phase=1:1:number_voxel_phase
                    % Get back voxel coordinate
                    x_p=P1(current_voxel_phase); y_p=P2(current_voxel_phase); z_p=P3(current_voxel_phase);
                    % Check chess pattern
                    id_= Label_lake(x_p,y_p,z_p);
                    % Chess pattern analysis only for voxels on the 'subdomains' Domain_size-1
                    if x_p>1 && y_p>1 && z_p>1 && x_p<Domain_size(1) && y_p<Domain_size(2) && z_p<Domain_size(3)
                        % Along plane 12
                        if (Label_lake(x_p+1,y_p,z_p)~=id_ && Label_lake(x_p-1,y_p,z_p)~=id_ && Label_lake(x_p,y_p+1,z_p)~=id_ && Label_lake(x_p,y_p-1,z_p)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                        % Along plane 13
                        if (Label_lake(x_p+1,y_p,z_p)~=id_ && Label_lake(x_p-1,y_p,z_p)~=id_ && Label_lake(x_p,y_p,z_p+1)~=id_ && Label_lake(x_p,y_p,z_p-1)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                        % Along plane 23
                        if (Label_lake(x_p,y_p+1,z_p)~=id_ && Label_lake(x_p,y_p-1,z_p)~=id_ && Label_lake(x_p,y_p,z_p+1)~=id_ && Label_lake(x_p,y_p,z_p-1)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                        % At the border
                    elseif (x_p==1 || x_p==Domain_size(1)) && (y_p>1 && z_p>1 && y_p<Domain_size(2) && z_p<Domain_size(3))
                        % Along plane 23
                        if (Label_lake(x_p,y_p+1,z_p)~=id_ && Label_lake(x_p,y_p-1,z_p)~=id_ && Label_lake(x_p,y_p,z_p+1)~=id_ && Label_lake(x_p,y_p,z_p-1)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                    elseif (y_p==1 || y_p==Domain_size(2)) && (x_p>1 && z_p>1 && x_p<Domain_size(1) && z_p<Domain_size(3))
                        % Along plane 13
                        if (Label_lake(x_p+1,y_p,z_p)~=id_ && Label_lake(x_p-1,y_p,z_p)~=id_ && Label_lake(x_p,y_p,z_p+1)~=id_ && Label_lake(x_p,y_p,z_p-1)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                    elseif (z_p==1 || z_p==Domain_size(3)) && (x_p>1 && y_p>1 && x_p<Domain_size(1) && y_p<Domain_size(2))
                        % Along plane 12
                        if (Label_lake(x_p+1,y_p,z_p)~=id_ && Label_lake(x_p-1,y_p,z_p)~=id_ && Label_lake(x_p,y_p+1,z_p)~=id_ && Label_lake(x_p,y_p-1,z_p)~=id_)
                            chess_pattern(x_p,y_p,z_p)=1;
                        end
                    end
                end
            end
            
            % For each voxel on the chess pattern, they will be attributed to the nearest particle id
            % Position of all chess voxels
            index_chess = find(chess_pattern==1);
            if length(index_chess)==length(sav_index_chess) && isequal(index_chess,sav_index_chess)
                % This mean iterative process has come to an end: this iteration
                % start has the same statut with the previous iteration start
                break
            end
            if iteration_chess>1
                if length(index_chess)>=length(sav_index_chess)
                    % No more amelioration
                    break
                end
            end
            sav_index_chess = index_chess; % Will be compared in the next iteration
            
                        
%             if iteration_==1 && visualize_2D && chess_pattern_removal
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
%                 slice_r(idx_source) = 1; slice_g(idx_source) = 0; slice_b(idx_source) = 0;
%                 slice_color(:,:,1)=slice_r; % Attribute RGB color
%                 slice_color(:,:,2)=slice_g;
%                 slice_color(:,:,3)=slice_b;
%                 
%                 hold(axe_chess(1).s,'on');
%                 slice_image = image(slice_color,'parent',axe_chess(1).s); % Display the slice
%                 axis(axe_chess(1).s,'tight'); % Fit the axes box
%                 axis(axe_chess(1).s,'equal'); % Aspect ratio is 1:1
%                 % Label
%                 xlabel(axe_chess(1).s,'Voxels');
%                 ylabel(axe_chess(1).s,'Voxels');
%                 set(axe_chess(1).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
%                 title({'Before chess pattern removal',['Number of particles identified: ' num2str(length(currentvalues)-1,'%i')]},'FontName','Times New Roman','FontSize',18,'Parent',axe_chess(1).s);
%                 set(axe_chess(1).s,'YDir','normal')
%                 hold(axe_chess(1).s,'off');    
%                                 
%                 tmp = slice_grey; tmp(index_chess) = 0;
%                 slice_color(:,:,1)=tmp; % Attribute RGB color
%                 tmp = slice_grey; tmp(index_chess) = 0;
%                 slice_color(:,:,2)=tmp;
%                 tmp = slice_grey; tmp(index_chess) = 1;
%                 slice_color(:,:,3)=tmp;
%                 slice_image = image(slice_color,'parent',axe_chess(2).s); % Display the slice
%                 axis(axe_chess(2).s,'tight'); % Fit the axes box
%                 axis(axe_chess(2).s,'equal'); % Aspect ratio is 1:1
%                 % Label
%                 xlabel(axe_chess(2).s,'Voxels');
%                 ylabel(axe_chess(2).s,'Voxels');
%                 set(axe_chess(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
%                 title('Chess pattern detection','FontName','Times New Roman','FontSize',18,'Parent',axe_chess(2).s);
%                 set(axe_chess(2).s,'YDir','normal')
%                 hold(axe_chess(2).s,'off');                  
%             end

            % Coordinnated of the voxels that belong to the chess pattern
            [Ch1,Ch2,Ch3] = ind2sub(Domain_size,index_chess);
            % Number of voxel that belong to the chess pattern
            number_voxel_chess = length(Ch1);
            if details_convergence==1
                fprintf ('   Number of voxel with a chess pattern: %i\n',number_voxel_chess);
            end
            % Calculate the center of mass of all particles, without taking into
            % consideration their voxels that belong to the chess pattern
            % Initialise matrix
            mass_center = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
            % Get all discrete particle id
            unique_particle = unique(Label_lake);
            unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
            number_particle = length(unique_particle);
            % Loop over all different particle id
            for k=1:1:number_particle
                current_id=unique_particle(k);
                % Find all voxels of this particle
                index_mass = find(Label_lake==current_id);
                % Remove from the analyis all these voxel that belong to the chess pattern
                actual_index_mass=[];
                n_m=length(index_mass);
                for kk=1:1:n_m
                    current_voxel_index=index_mass(kk);
                    index_test = find(index_chess==current_voxel_index);
                    if isempty(index_test)
                        % This voxel does not belong to the chess pattern and will
                        % be considered for the mass center calculation
                        actual_index_mass = [actual_index_mass current_voxel_index];
                    end
                end
                % Get the mass center
                if ~isempty(actual_index_mass) % To forecast the case this particle is made only of a chess pattern
                    % Get all voxel coordinates
                    [M1,M2,M3] = ind2sub(Domain_size,actual_index_mass);
                    % Deduce the mass center
                    x_mass=round(mean(M1)); y_mass=round(mean(M2)); z_mass=round(mean(M3));
                    mass_center(x_mass,y_mass,z_mass)=current_id;
                else
                    %disp empty
                end
            end
            
            % Mass center position
            % Get all voxels of the mass center
            index_mass_position = find(mass_center~=0);
            % Number of mass center
            n_mass_center_=length(index_mass_position);
            % Get their coordinates
            [MC1,MC2,MC3] = ind2sub(Domain_size,index_mass_position);
            mass_center_position = zeros(n_mass_center_,5);
            for k_=1:1:n_mass_center_
                % Position
                xx_=MC1(k_); yy_=MC2(k_); zz_=MC3(k_);
                mass_center_position(k_,1)=mass_center(xx_,yy_,zz_); % Id
                mass_center_position(k_,2)=xx_; % Position x
                mass_center_position(k_,3)=yy_; % Position y
                mass_center_position(k_,4)=zz_; % Position z
            end
            
            % Loop over all voxel that belong to the chess pattern
            for current_voxel_chess=1:1:number_voxel_chess
                % Get back voxel coordinate
                x_ch=Ch1(current_voxel_chess); y_ch=Ch2(current_voxel_chess); z_ch=Ch3(current_voxel_chess);
                % Calculate distance between each mass center and the current voxel
                voxel_pos=zeros(n_mass_center_,3);
                voxel_pos(:,1)=x_ch;
                voxel_pos(:,2)=y_ch;
                voxel_pos(:,3)=z_ch;
                mass_center_position(:,5) = ((mass_center_position(:,2)-voxel_pos(:,1)).^2 + (mass_center_position(:,3)-voxel_pos(:,2)).^2 + (mass_center_position(:,4)-voxel_pos(:,3)).^2).^(0.5);
                % Sort by increasing distance
                mass_center_position = sortrows(mass_center_position,5);
                % This iteration of mass center will be potentially modified for the analysis of each different voxel
                % Each new voxel start with an unchanged mass center matrix
                mass_center_position_modified = mass_center_position;
                
                keep_trying=1;
                while keep_trying==1
                    if isempty(mass_center_position_modified)
                        % That mean we have tried all the mass center. This voxel is likely an isolated voxel we will deal with in the next step
                        break
                    end
                    % Identify the nearest mass center for the voxel investigated
                    nearest_particle_id = mass_center_position_modified(1,1);
                    % We will check there is a clear line of view between this mass center and the chess voxel. If not, the next nearest seed will be
                    % considered and so on (through a modification of the mass center matrix, stored in mass_center_modified)
                    % Get back mass center coordinate
                    x_m = mass_center_position_modified(1,2); y_m = mass_center_position_modified(1,3); z_m = mass_center_position_modified(1,4);
                    % Convert in our metrics
                    x_ch_=x_ch-0.5;y_ch_=y_ch-0.5;z_ch_=z_ch-0.5;
                    x_m_=x_m-0.5;y_m_=y_m-0.5;z_m_=z_m-0.5;
                    % Evaluate distance
                    dx_12 = x_m_-x_ch_; dy_12 = y_m_-y_ch_; dz_12 = z_m_-z_ch_;
                    r12 = sqrt(dx_12^2+dy_12^2+dz_12^2);
                    % Number of step required to test voxels between them
                    n_step = floor(r12);
                    % Delta position
                    dx_ = dx_12/n_step;
                    dy_ = dy_12/n_step;
                    dz_ = dz_12/n_step;
                    % Test loop
                    test_=1; iteration_test=0;
                    while (test_==1 && (iteration_test+1)<=n_step)
                        iteration_test=iteration_test+1;
                        % Tested position
                        x_tested = x_ch_+iteration_test*dx_;
                        y_tested = y_ch_+iteration_test*dy_;
                        z_tested = z_ch_+iteration_test*dz_;
                        % Converted in voxels coordinates
                        x_tested=floor(x_tested)+1;
                        y_tested=floor(y_tested)+1;
                        z_tested=floor(z_tested)+1;
                        % Check voxel
                        test_ = binary_phase(x_tested,y_tested,z_tested);
                    end
                    if test_==1
                        % Clear line of view between the voxel exists
                        % Attribute the chess voxel with this id
                        Label_lake(x_ch,y_ch,z_ch)=nearest_particle_id;
                        % Sucees for this voxel!
                        keep_trying=0;
                    else
                        % We need to try the next nearest particle
                        mass_center_position_modified(1,:)=[];
                    end
                end
            end
            % Re-order seed id
            unique_id = unique(Label_lake);
            number_id = length(unique_id);
            for current_id=2:1:number_id % 0 is the complementary phase
                id_to_be_replaced = unique_id(current_id);
                index_ = find(Label_lake==id_to_be_replaced);
                Label_lake(index_)=current_id-1;
            end
        end
        % Get all discrete particle id
        unique_particle = unique(Label_lake);
        unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
        % Get number of particle
        number_particle = length(unique_particle);
        if details_convergence==1
            fprintf ('   - Chess pattern removed (step1):     number of discrete particle: %i\n',number_particle);
        end
        
        
        
        %% STEP 2: REMOVE THE LAST CHESS PATTERN, IF ANY
        
        % Some voxels identified as belonging to a chess pattern may still be there
        % These particular voxels are often sandwiched between two particle of the
        % complementary phases, that block the line of view
        % We attribute them according to their immediate surrunding
        
        % While loop initialisation
        sav_index_chess=[];
        number_voxel_chess=1;
        while number_voxel_chess~=0
            
            % Find index of all remaining voxel that still belong to the chess pattern
            index_chess = find(chess_pattern==1);
            % Check while loop stop
            if length(index_chess)==length(sav_index_chess) && isequal(index_chess,sav_index_chess)
                % This mean iterative process has come to an end: this iteration
                % start has the same statut with the previous iteration start
                break
            end
            sav_index_chess = index_chess; % Will be compared in the next iteration
            
            % Get all their coordinate
            [Ch1,Ch2,Ch3] = ind2sub(Domain_size,index_chess);
            % Number of these voxels
            number_voxel_chess = length(Ch1);
            % Loop over all voxels that belong to the chess pattern
            for current_voxel_chess=1:1:number_voxel_chess % Note: if n_voxel=0, this loop is ignored by Matlalb
                % Get back coordinate
                x_ch=Ch1(current_voxel_chess); y_ch=Ch2(current_voxel_chess); z_ch=Ch3(current_voxel_chess);
                % Get surrounding voxel value
                surrounding_values=[];
                % x-
                if x_ch>1
                    x_test=x_ch-1; y_test=y_ch; z_test=z_ch;
                    if binary_phase(x_test,y_test,z_test)==1 && chess_pattern(x_test,y_test,z_test)==0
                        % This adjecent voxel belong to the phase and does not belong to a chess pattern
                        surrounding_values=[surrounding_values Label_lake(x_test,y_test,z_test)];
                    end
                end
                % y-
                if y_ch>1
                    x_test=x_ch; y_test=y_ch-1; z_test=z_ch;
                    if binary_phase(x_test,y_test,z_test)==1 && chess_pattern(x_test,y_test,z_test)==0
                        % This adjecent voxel belong to the phase and does not belong to a chess pattern
                        surrounding_values=[surrounding_values Label_lake(x_test,y_test,z_test)];
                    end
                end
                % z-
                if z_ch>1
                    x_test=x_ch; y_test=y_ch; z_test=z_ch-1;
                    if binary_phase(x_test,y_test,z_test)==1 && chess_pattern(x_test,y_test,z_test)==0
                        % This adjecent voxel belong to the phase and does not belong to a chess pattern
                        surrounding_values=[surrounding_values Label_lake(x_test,y_test,z_test)];
                    end
                end
                % x+
                if x_ch<Domain_size(1)
                    x_test=x_ch+1; y_test=y_ch; z_test=z_ch;
                    if binary_phase(x_test,y_test,z_test)==1 && chess_pattern(x_test,y_test,z_test)==0
                        % This adjecent voxel belong to the phase and does not belong to a chess pattern
                        surrounding_values=[surrounding_values Label_lake(x_test,y_test,z_test)];
                    end
                end
                % y+
                if y_ch<Domain_size(2)
                    x_test=x_ch; y_test=y_ch+1; z_test=z_ch;
                    if binary_phase(x_test,y_test,z_test)==1 && chess_pattern(x_test,y_test,z_test)==0
                        % This adjecent voxel belong to the phase and does not belong to a chess pattern
                        surrounding_values=[surrounding_values Label_lake(x_test,y_test,z_test)];
                    end
                end
                % z+
                if z_ch<Domain_size(3)
                    x_test=x_ch; y_test=y_ch; z_test=z_ch+1;
                    if binary_phase(x_test,y_test,z_test)==1 && chess_pattern(x_test,y_test,z_test)==0
                        % This adjecent voxel belong to the phase and does not belong to a chess pattern
                        surrounding_values=[surrounding_values Label_lake(x_test,y_test,z_test)];
                    end
                end
                % Replace with the id that is preominent
                if ~isempty(surrounding_values)
                    % Unique value at the surrounding
                    unique_surrounding = unique(surrounding_values);
                    % Get the value that is preominemt
                    sum_surrounding=zeros(length(unique_surrounding),1);
                    for n_=1:1:length(unique_surrounding)
                        sum_surrounding(n_)=sum(surrounding_values==unique_surrounding(n_));
                    end
                    max_index=find(sum_surrounding==max(sum_surrounding));
                    max_index=max_index(1); % In case there is equality
                    preominent_value = unique_surrounding(max_index);
                    % Replace chess voxel
                    Label_lake(x_ch,y_ch,z_ch)=preominent_value;
                    % Update chess pattern matrix
                    chess_pattern(x_ch,y_ch,z_ch)=0;
                end
            end
        end
        % Get all discrete particle id
        unique_particle = unique(Label_lake);
        unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
        % Get number of particle
        number_particle = length(unique_particle);
        if details_convergence==1
            fprintf ('   - Chess pattern removed (step2):     number of discrete particle: %i\n',number_particle);
        end
        
        
        
    end
    
    
    
    %% STEP 1c: SEARCH FOR NON-CONTINGUOUS PARTICLE
    
    % We use the Matlab built-in function for detetecting and labeling non continguous particle/cluster
    % L = bwlabeln(BW,6)
    % L is containing labels for the connected components in BW
    % L is the same size as BW
    % The elements of L are integer values greater than or equal to 0.
    % The pixels labeled 0 are the background. The pixels labeled 1 make up one connected cluster. The pixels labeled 2 make up a second connected cluster, and so on
    
    % Get all unique id particle
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove id 0 (complementary phase)
    % Get number of different id
    number_id = length(unique_particle);
    % Loop over id
    for current_=1:1:number_id
        % Get id
        current_id = unique_particle(current_);
        % Create binary matrix for this particle id
        binary_particle = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
        binary_particle(Label_lake==current_id)=1;
        % Label particle
        L_particle = bwlabeln(binary_particle,6);
        % Unique label
        unique_cluster = unique(L_particle);
        unique_cluster(1)=[]; % remobe background cluster id
        % Number of distinct cluster
        cluster_number=length(unique_cluster);
        % Treat the case where there are more than 1 cluster
        if cluster_number>1
            % The largest cluster will be unchanged
            % All the other cluster will be treated as follow:
            % Case a) If they are surrounded by ONLY voxels of the complementary phase, then this cluster is considered as a new particle
            % Case b) otherwise, voxels of this cluster will be attributed to the particle id of the nearest voxel that belong to the phase and
            % that does not belong to this particular cluster
            cluster_size=zeros(cluster_number,2);
            % Column: cluster id, size
            for current_cluster=1:1:cluster_number
                cluster_id = unique_cluster(current_cluster);
                cluster_size(current_cluster,1)=cluster_id;
                cluster_size(current_cluster,2)=sum(sum(sum(L_particle==cluster_id)));
            end
            % Sort per cluster size (decreasing order)
            cluster_size = sortrows(cluster_size,-2);
            % We will deal only with the second, third etc. cluster
            for analysed_cluster=2:1:cluster_number
                % Get cluster id
                cluster_id = cluster_size(analysed_cluster,1);
                % Create cluster binary matrix
                binary_cluster = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
                binary_cluster(L_particle==cluster_id)=1;
                % Find all voxels of this cluster
                index_cluster=find(binary_cluster==1);
                % Get all their coordinate
                [Cl1,Cl2,Cl3] = ind2sub(Domain_size,index_cluster);
                % Number of voxel
                number_voxel_cluster = length(Cl1);
                % % Identify voxel of the complemtary that belong at the border
                % They will be marked 2 in the binary phase
                surrounding_values=[];
                pos_surrounding=zeros(1,5);
                for current_voxel=1:1:number_voxel_cluster
                    % Get back voxel coordinate
                    x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                    % Check x-
                    if x_>1
                        if Label_lake(x_-1,y_,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_-1,y_,z_)];
                            if Label_lake(x_-1,y_,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_-1,y_,z_,0,Label_lake(x_-1,y_,z_)]];
                            end
                        end
                    end
                    % Check y-
                    if y_>1
                        if Label_lake(x_,y_-1,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_-1,z_)];
                            if Label_lake(x_,y_-1,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_-1,z_,0,Label_lake(x_,y_-1,z_)]];
                            end
                        end
                    end
                    % Check z-
                    if z_>1
                        if Label_lake(x_,y_,z_-1)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_,z_-1)];
                            if Label_lake(x_,y_,z_-1)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_,z_-1,0,Label_lake(x_,y_,z_-1)]];
                            end
                        end
                    end
                    % Check x+
                    if x_<Domain_size(1)
                        if Label_lake(x_+1,y_,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_+1,y_,z_)];
                            if Label_lake(x_+1,y_,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_+1,y_,z_,0,Label_lake(x_+1,y_,z_)]];
                            end
                        end
                    end
                    % Check y+
                    if y_<Domain_size(2)
                        if Label_lake(x_,y_+1,z_)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_+1,z_)];
                            if Label_lake(x_,y_+1,z_)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_+1,z_,0,Label_lake(x_,y_+1,z_)]];
                            end
                        end
                    end
                    % Check z+
                    if z_<Domain_size(3)
                        if Label_lake(x_,y_,z_+1)~=current_id
                            surrounding_values=[surrounding_values Label_lake(x_,y_,z_+1)];
                            if Label_lake(x_,y_,z_+1)~=0
                                pos_surrounding=[pos_surrounding; [x_,y_,z_+1,0,Label_lake(x_,y_,z_+1)]];
                            end
                        end
                    end
                end
                % Remove fist 0 line
                pos_surrounding(1,:)=[];
                [n_surrounding,~]=size(pos_surrounding);
                % Replace accoring to case a or case b
                if length(unique(surrounding_values))==1 && surrounding_values(1)==0
                    % Case a:
                    unique_particle = unique(Label_lake);
                    max_id = max(unique_particle);
                    new_id = max_id+1;
                    Label_lake(index_cluster)=new_id;
                else
                    % Case b:
                    % Calculate minimum distance to nearest voxel
                    for current_voxel=1:1:number_voxel_cluster
                        % Get back voxel coordinate
                        x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                        % Calculate all distance betwenen the voxel and the surrounding voxels
                        for voxel_surround=1:1:n_surrounding
                            x_s=pos_surrounding(voxel_surround,1);
                            y_s=pos_surrounding(voxel_surround,2);
                            z_s=pos_surrounding(voxel_surround,3);
                            dist_surround = sqrt((x_s-x_)^2 + (y_s-y_)^2 + (z_s-z_)^2);
                            pos_surrounding(voxel_surround,4)=dist_surround;
                        end
                        % Sort by increasing order
                        pos_surrounding = sortrows(pos_surrounding,4);
                        % Get minimum distance
                        min_dist = pos_surrounding(1,4);
                        % Get all index with this distance
                        index_min = find(pos_surrounding(:,4)==min_dist);
                        % Get all id
                        list_id=[];
                        for kkk=1:1:length(index_min)
                            list_id = [list_id pos_surrounding(kkk,5)];
                        end
                        % Get unique id
                        unique_surrounding = unique(list_id);
                        % Get the value that is preominemt
                        sum_surrounding=zeros(length(unique_surrounding),1);
                        for n_=1:1:length(unique_surrounding)
                            sum_surrounding(n_)=sum(list_id==unique_surrounding(n_));
                        end
                        max_index=find(sum_surrounding==max(sum_surrounding));
                        max_index=max_index(1); % In case there is equality
                        preominent_value = unique_surrounding(max_index);
                        % Replace value
                        Label_lake(x_,y_,z_)=preominent_value;
                    end
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
        fprintf ('   - Search for non connected particle: number of discrete particle: %i\n',number_particle);
    end
    
    %% STEP 2: CASE OF THE X VOXEL SIZE PARTICLE
    
    % X voxel-size particle will be removed
    % Since I do NOT want parameter, i will consider only the case X=1
    critical_minimal_size = 1;
    
    % If you really want to use it, you should consider iterate from smaller particle to
    % larger particle (modify the algorithm in consequence)
    
    % Get all unique id particle
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove id 0 (complementary phase)
    % Get number of different id
    number_id = length(unique_particle);
    % Loop over id
    for current_=1:1:number_id
        % Get id
        current_id = unique_particle(current_);
        % Size
        D_PSD_particle_size = sum(sum(sum(Label_lake==current_id)));
        if D_PSD_particle_size<=critical_minimal_size
            index_particle=find(Label_lake==current_id);
            % Get all their coordinate
            [Cl1,Cl2,Cl3] = ind2sub(Domain_size,index_particle);
            % Number of voxel
            number_voxel_particle = length(Cl1);
            % % Identify voxel of the complemtary that belong at the border
            surrounding_values=[];
            pos_surrounding=zeros(1,5);
            % Loop over all voxel of the particle
            for current_voxel=1:1:number_voxel_particle
                % Get back voxel coordinate
                x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                % Check x-
                if x_>1
                    if Label_lake(x_-1,y_,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_-1,y_,z_)];
                        if Label_lake(x_-1,y_,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_-1,y_,z_,0,Label_lake(x_-1,y_,z_)]];
                        end
                    end
                end
                % Check y-
                if y_>1
                    if Label_lake(x_,y_-1,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_-1,z_)];
                        if Label_lake(x_,y_-1,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_-1,z_,0,Label_lake(x_,y_-1,z_)]];
                        end
                    end
                end
                % Check z-
                if z_>1
                    if Label_lake(x_,y_,z_-1)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_,z_-1)];
                        if Label_lake(x_,y_,z_-1)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_,z_-1,0,Label_lake(x_,y_,z_-1)]];
                        end
                    end
                end
                % Check x+
                if x_<Domain_size(1)
                    if Label_lake(x_+1,y_,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_+1,y_,z_)];
                        if Label_lake(x_+1,y_,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_+1,y_,z_,0,Label_lake(x_+1,y_,z_)]];
                        end
                    end
                end
                % Check y+
                if y_<Domain_size(2)
                    if Label_lake(x_,y_+1,z_)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_+1,z_)];
                        if Label_lake(x_,y_+1,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_+1,z_,0,Label_lake(x_,y_+1,z_)]];
                        end
                    end
                end
                % Check z+
                if z_<Domain_size(3)
                    if Label_lake(x_,y_,z_+1)~=current_id
                        surrounding_values=[surrounding_values Label_lake(x_,y_,z_+1)];
                        if Label_lake(x_,y_,z_+1)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_,z_+1,0,Label_lake(x_,y_,z_+1)]];
                        end
                    end
                end
            end
            
            % Remove fist 0 line
            pos_surrounding(1,:)=[];
            [n_surrounding,~]=size(pos_surrounding);
            
            % Replace accoring to case a or case b
            if length(unique(surrounding_values))==1 && surrounding_values(1)==0
                % Case a: is removed
                Label_lake(index_particle)=0;
            else
                % Case b:
                % Calculate minimum distance to nearest voxel
                for current_voxel=1:1:number_voxel_particle
                    % Get back voxel coordinate
                    x_=Cl1(current_voxel); y_=Cl2(current_voxel); z_=Cl3(current_voxel);
                    % Calculate all distance betwenen the voxel and the surrounding voxels
                    for voxel_surround=1:1:n_surrounding
                        x_s=pos_surrounding(voxel_surround,1);
                        y_s=pos_surrounding(voxel_surround,2);
                        z_s=pos_surrounding(voxel_surround,3);
                        dist_surround = sqrt((x_s-x_)^2 + (y_s-y_)^2 + (z_s-z_)^2);
                        pos_surrounding(voxel_surround,4)=dist_surround;
                    end
                    % Sort by increasing order
                    pos_surrounding = sortrows(pos_surrounding,4);
                    % Get minimum distance
                    min_dist = pos_surrounding(1,4);
                    % Get all index with this distance
                    index_min = find(pos_surrounding(:,4)==min_dist);
                    % Get all id
                    list_id=[];
                    for kkk=1:1:length(index_min)
                        list_id = [list_id pos_surrounding(kkk,5)];
                    end
                    % Get unique id
                    unique_surrounding = unique(list_id);
                    % Get the value that is preominemt
                    sum_surrounding=zeros(length(unique_surrounding),1);
                    for n_=1:1:length(unique_surrounding)
                        sum_surrounding(n_)=sum(list_id==unique_surrounding(n_));
                    end
                    max_index=find(sum_surrounding==max(sum_surrounding));
                    max_index=max_index(1); % In case there is equality
                    preominent_value = unique_surrounding(max_index);
                    % Replace value
                    Label_lake(x_,y_,z_)=preominent_value;
                end
            end
            
        end
    end
    % Re-order seed id
    unique_id = unique(Label_lake);
    number_id = length(unique_id);
    for current_id=2:1:number_id % 0 is the complementary phase
        id_to_be_replaced = unique_id(current_id);
        index_ = find(Label_lake==id_to_be_replaced);
        Label_lake(index_)=current_id-1;
    end
    % Get all discrete particle id
    unique_particle = unique(Label_lake);
    unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
    % Get number of particle
    number_particle = length(unique_particle);
    if details_convergence==1
        fprintf ('   - One-voxel size particle treated:   number of discrete particle: %i\n',number_particle);
    end
    
    
    if cpsd_refining==1
        %% STEP 3: Discrete-PSD
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
        if details_convergence==1
            fprintf ('   - Discrete particle size calculated\n');
        end
        if visualize_2D && cpsd_refining
            hold(axe_video(2).s,'on'); % Active subplot
            imagesc(axe_video(2).s,D_PSD_particle_size);
            axis(axe_video(2).s,'tight'); % Fit the axes box
            axis(axe_video(2).s,'equal'); % Aspect ratio is 1:1
            colorbar(axe_video(2).s);
            % Label
            xlabel(axe_video(2).s,'Voxels');
            ylabel(axe_video(2).s,'Voxels');
            set(axe_video(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
            title('Discrete particle size','FontName','Times New Roman','FontSize',14,'Parent',axe_video(2).s);
            box(axe_video(2).s,'on')
            hold(axe_video(2).s,'off'); % Active subplot
        end
        
        %% STEP 4: Continuum-PSD
        % It will be calculated only one time
        if iteration_==1
            % Run the C-PSD algorithm
            [C_PSD_particle_size] = Function_particle_size_CPSD_Algorithm(binary_phase);
            if visualize_2D && cpsd_refining
                hold(axe_video(1).s,'on'); % Active subplot
                imagesc(axe_video(1).s,C_PSD_particle_size);
                axis(axe_video(1).s,'tight'); % Fit the axes box
                axis(axe_video(1).s,'equal'); % Aspect ratio is 1:1
                colorbar(axe_video(1).s);
                % Label
                xlabel(axe_video(1).s,'Voxels');
                ylabel(axe_video(1).s,'Voxels');
                set(axe_video(1).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
                title('Continuum spherical particle size','FontName','Times New Roman','FontSize',14,'Parent',axe_video(1).s);
                box(axe_video(1).s,'on')
                hold(axe_video(1).s,'off'); % Active subplot
            end
        end
        
        
        
        %% STEP 5: COMPARE SIZE BETWEEN C-PSD AND D-PSD
        % Initialise
        smaller_than_cpsd = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
        % Set =1 when D-PSD is smaller than C-PSD
        smaller_than_cpsd(D_PSD_particle_size<C_PSD_particle_size)=1;
        % % Two cases can explaines D-PSD < C-PSD
        % Case a) The particle is at the domain's edge: C-PSD algorithm does not exibit
        % edge effect compare to the D-PSD that will decrease the particle size at
        % the border. Then, there is nothing to do for these particles
        % Case b) these voxels are likely stuck between two large particles. We will replace them
        
        % L = bwlabeln(BW,6)
        % L is containing labels for the connected components in BW
        % L is the same size as BW
        % The elements of L are integer values greater than or equal to 0.
        % The pixels labeled 0 are the background. The pixels labeled 1 make up one connected cluster. The pixels labeled 2 make up a second connected cluster, and so on
        % Label zone
        L_zone = bwlabeln(smaller_than_cpsd,6);
        
        % % Step1: remove from the analysis zone located at the domain's edge
        % Unique zone
        unique_zone = unique(Label_lake);
        unique_zone(1)=[]; % remove background
        % Number of distinct zone
        zone_number=length(unique_zone);
        % Loop over all zones
        for current_=1:1:zone_number
            % Get id
            current_zone = unique_zone(current_);
            % Get all voxels
            index_zone = find(Label_lake==current_zone);
            % Get all coordinates
            [Z1,Z2,Z3] = ind2sub(Domain_size,index_zone);
            % Get min max
            x_min=min(Z1); y_min=min(Z2); z_min=min(Z3);
            x_max=max(Z1); y_max=max(Z2); z_max=max(Z3);
            % Check if the zone is located at the border
            if Domain_size(3)>1
                % 3D case
                at_the_border = (x_min==1 || y_min==1 || z_min==1 || x_max==Domain_size(1) || y_max==Domain_size(2) || z_max==Domain_size(3));
            else
                % 2D case
                at_the_border = (x_min==1 || y_min==1 || x_max==Domain_size(1) || y_max==Domain_size(2));
            end
            % Remove the zone from the analysis if true
            if at_the_border==1
                L_zone(index_zone)=0;
                smaller_than_cpsd(index_zone)=0;
            end
        end
        
        if visualize_2D && cpsd_refining
            idx_smaller_than_cpsd = find(smaller_than_cpsd==1);
            hold(axe_video(3).s,'on'); % Active subplot
            tmp = slice_grey; tmp(idx_smaller_than_cpsd) = 0;
            slice_color(:,:,1)=tmp; % Attribute RGB color
            tmp = slice_grey; tmp(idx_smaller_than_cpsd) = 0;
            slice_color(:,:,2)=tmp;
            tmp = slice_grey; tmp(idx_smaller_than_cpsd) = 1;
            slice_color(:,:,3)=tmp;
            slice_image = image(slice_color,'parent',axe_video(3).s); % Display the slice
            axis(axe_video(3).s,'tight'); % Fit the axes box
            axis(axe_video(3).s,'equal'); % Aspect ratio is 1:1
            set(axe_video(3).s,'YDir','normal')
            % Label
            xlabel(axe_video(3).s,'Voxels');
            ylabel(axe_video(3).s,'Voxels');
            set(axe_video(3).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
            title('Region with d-PSD(x) < c-PSD(x)','FontName','Times New Roman','FontSize',14,'Parent',axe_video(3).s);
            box(axe_video(3).s,'on')
            hold(axe_video(3).s,'off'); % Active subplot
            
            all_idx_smaller_than_cpsd = [all_idx_smaller_than_cpsd; idx_smaller_than_cpsd];
            tmp = slice_grey; tmp(all_idx_smaller_than_cpsd) = 0;
            slice_color(:,:,1)=tmp; % Attribute RGB color
            tmp = slice_grey; tmp(all_idx_smaller_than_cpsd) = 0;
            slice_color(:,:,2)=tmp;
            tmp = slice_grey; tmp(all_idx_smaller_than_cpsd) = 1;
            slice_color(:,:,3)=tmp;
            slice_image = image(slice_color,'parent',axe_illustration(2).s); % Display the slice
            axis(axe_illustration(2).s,'tight'); % Fit the axes box
            axis(axe_illustration(2).s,'equal'); % Aspect ratio is 1:1
            set(axe_illustration(2).s,'YDir','normal')
            % Label
            xlabel(axe_illustration(2).s,'Voxels');
            ylabel(axe_illustration(2).s,'Voxels');
            set(axe_illustration(2).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
            title('Region with d-PSD(x) < c-PSD(x)','FontName','Times New Roman','FontSize',14,'Parent',axe_illustration(2).s);
            box(axe_illustration(2).s,'on')
            hold(axe_illustration(2).s,'off'); % Active subplot
            
        end
        
        % % Step 2: replace the remaining voxel with nearest voxel
        % Unique zone
        unique_zone = unique(L_zone);
        unique_zone(1)=[]; % remove background
        % Number of distinct zone
        zone_number=length(unique_zone);
        % Loop over all zones
        for current_=1:1:zone_number
            % Get id
            current_zone = unique_zone(current_);
            % Get all voxels
            index_zone = find(L_zone==current_zone);
            % Get all coordinates
            [Z1,Z2,Z3] = ind2sub(Domain_size,index_zone);
            % Number of voxel
            number_voxel_zone  = length(Z1);
            % % Identify voxel of the complemtary that belong at the border
            pos_surrounding=zeros(1,5);
            for current_voxel=1:1:number_voxel_zone
                % Get back voxel coordinate
                x_=Z1(current_voxel); y_=Z2(current_voxel); z_=Z3(current_voxel);
                % Check x-
                if x_>1
                    if L_zone(x_-1,y_,z_)~=current_zone
                        if Label_lake(x_-1,y_,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_-1,y_,z_,0,Label_lake(x_-1,y_,z_)]];
                        end
                    end
                end
                % Check y-
                if y_>1
                    if L_zone(x_,y_-1,z_)~=current_zone
                        if Label_lake(x_,y_-1,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_-1,z_,0,Label_lake(x_,y_-1,z_)]];
                        end
                    end
                end
                % Check z-
                if z_>1
                    if L_zone(x_,y_,z_-1)~=current_zone
                        if Label_lake(x_,y_,z_-1)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_,z_-1,0,Label_lake(x_,y_,z_-1)]];
                        end
                    end
                end
                % Check x+
                if x_<Domain_size(1)
                    if L_zone(x_+1,y_,z_)~=current_zone
                        if Label_lake(x_+1,y_,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_+1,y_,z_,0,Label_lake(x_+1,y_,z_)]];
                        end
                    end
                end
                % Check y+
                if y_<Domain_size(2)
                    if L_zone(x_,y_+1,z_)~=current_zone
                        if Label_lake(x_,y_+1,z_)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_+1,z_,0,Label_lake(x_,y_+1,z_)]];
                        end
                    end
                end
                % Check z+
                if z_<Domain_size(3)
                    if L_zone(x_,y_,z_+1)~=current_zone
                        if Label_lake(x_,y_,z_+1)~=0
                            pos_surrounding=[pos_surrounding; [x_,y_,z_+1,0,Label_lake(x_,y_,z_+1)]];
                        end
                    end
                end
            end
            % Remove fist 0 line
            pos_surrounding(1,:)=[];
            [n_surrounding,~]=size(pos_surrounding);
            
            if n_surrounding>0
                % Calculate minimum distance to nearest voxel
                for current_voxel=1:1:number_voxel_zone
                    % Get back voxel coordinate
                    x_=Z1(current_voxel); y_=Z2(current_voxel); z_=Z3(current_voxel);
                    % Calculate all distance betwenen the voxel and the surrounding voxels
                    for voxel_surround=1:1:n_surrounding
                        x_s=pos_surrounding(voxel_surround,1);
                        y_s=pos_surrounding(voxel_surround,2);
                        z_s=pos_surrounding(voxel_surround,3);
                        dist_surround = sqrt((x_s-x_)^2 + (y_s-y_)^2 + (z_s-z_)^2);
                        pos_surrounding(voxel_surround,4)=dist_surround;
                    end
                    % Sort by increasing order
                    pos_surrounding = sortrows(pos_surrounding,4);
                    % Get minimum distance
                    min_dist = pos_surrounding(1,4);
                    % Get all index with this distance
                    index_min = find(pos_surrounding(:,4)==min_dist);
                    % Get all id
                    list_id=[];
                    for kkk=1:1:length(index_min)
                        list_id = [list_id pos_surrounding(kkk,5)];
                    end
                    % Get unique id
                    unique_surrounding = unique(list_id);
                    % Get the value that is preominemt
                    sum_surrounding=zeros(length(unique_surrounding),1);
                    for n_=1:1:length(unique_surrounding)
                        sum_surrounding(n_)=sum(list_id==unique_surrounding(n_));
                    end
                    max_index=find(sum_surrounding==max(sum_surrounding));
                    max_index=max_index(1); % In case there is equality
                    preominent_value = unique_surrounding(max_index);
                    % Replace value
                    Label_lake(x_,y_,z_)=preominent_value;
                end
            end
        end
        % Re-order seed id
        unique_id = unique(Label_lake);
        number_id = length(unique_id);
        for current_id=2:1:number_id % 0 is the complementary phase
            id_to_be_replaced = unique_id(current_id);
            index_ = find(Label_lake==id_to_be_replaced);
            Label_lake(index_)=current_id-1;
        end
        % Get all discrete particle id
        unique_particle = unique(Label_lake);
        unique_particle(1)=[]; % Remove the 0, allocated to the complementary phase
        % Get number of particle
        number_particle = length(unique_particle);
        if details_convergence==1
            fprintf ('   - C-PSD & D-PSD analysis:            number of discrete particle: %i\n',number_particle);
        end
        
        
        if visualize_2D
            slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
            slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
            slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
            slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
            slice_grey = zeros(Domain_size(1),Domain_size(2)); % Blue color map
            currentvalues = unique(Label_lake);
            for k=1:1:length(currentvalues)
                idx=find(Label_lake==currentvalues(k));
                if currentvalues(k)==0
                    slice_r(idx) = 1; slice_g(idx) = 1; slice_b(idx) = 1; % Background
                    slice_grey(idx) = 1;
                else
                    slice_r(idx) = lake_id_RGB_color(k,1); slice_g(idx) = lake_id_RGB_color(k,2); slice_b(idx) = lake_id_RGB_color(k,3); % Background
                    slice_grey(idx) = lake_id_grey_color(k);
                end
            end
            % Edges
            [index_border_phase,~,~,~] = Function_identify_edges(binary_phase);
            slice_r(index_border_phase) = 0; slice_g(index_border_phase) = 0; slice_b(index_border_phase) = 0;
            slice_grey(index_border_phase) = 0;
            % Watershed lines
            % Edge detection
            background = 0;
            edgewithbackground = false;
            [index_border_label,~,~,~] = Function_identify_labelsedges(Label_lake, background, edgewithbackground);
            tmp = slice_grey; tmp(index_border_label) = 1;
            slice_color(:,:,1)=tmp; % Attribute RGB color
            tmp = slice_grey; tmp(index_border_label) = 0;
            slice_color(:,:,2)=tmp;
            tmp = slice_grey; tmp(index_border_label) = 0;
            slice_color(:,:,3)=tmp;
            slice_image = image(slice_color,'parent',axe_video(5).s); % Display the slice
            axis(axe_video(5).s,'tight'); % Fit the axes box
            axis(axe_video(5).s,'equal'); % Aspect ratio is 1:1
            set(axe_video(5).s,'YDir','normal')
            % Label
            xlabel(axe_video(5).s,'Voxels');
            ylabel(axe_video(5).s,'Voxels');
            set(axe_video(5).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
            title('Particle identification with watershed lines (after)','FontName','Times New Roman','FontSize',14,'Parent',axe_video(5).s);
            box(axe_video(5).s,'on')
            hold(axe_video(5).s,'off'); % Active subplot
            % Video
            stored_frame(iteration_) = getframe(Fig_video);
            writeVideo(video_handle,stored_frame(iteration_))
            
            slice_image = image(slice_color,'parent',axe_illustration(3).s); % Display the slice
            axis(axe_illustration(3).s,'tight'); % Fit the axes box
            axis(axe_illustration(3).s,'equal'); % Aspect ratio is 1:1
            set(axe_illustration(3).s,'YDir','normal')
            % Label
            xlabel(axe_illustration(3).s,'Voxels');
            ylabel(axe_illustration(3).s,'Voxels');
            set(axe_illustration(3).s,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
            title('Particle identification with watershed lines (after)','FontName','Times New Roman','FontSize',14,'Parent',axe_illustration(3).s);
            box(axe_illustration(3).s,'on')
            hold(axe_illustration(3).s,'off'); % Active subplot
            
        end
        
    end
    
    %% CHANGE FROM PREVIOUS ITERATION: EXIT WHILE LOOP TEST
    % Compare with previous iteration resutt
    if iteration_>1
        number_change=sum(sum(sum(previous_discrete_particle~=Label_lake)));
        if details_convergence==1
            fprintf ('   CONVERGENCE CHECK: Number of voxels changed compared with previous iteration: %i\n',number_change);
        end
    end
    % For the next iteration
    previous_discrete_particle = Label_lake;
    
    %% CHECKSUM
    
    % Current Discrete_particle checksum
    % Checksum_currentstate_Label_lake = DataHash(Label_lake);
    % That may induce an OUT OF MEMORY error if Label_lake is too large
    % Instead, we are calculating the hash slice by slice
    Checksum_currentstate_Label_lake=[];
    for current_=1:1:Domain_size(3)
        current_checksum = DataHash(Label_lake(:,:,current_));
        Checksum_currentstate_Label_lake=[Checksum_currentstate_Label_lake current_checksum];
    end
    
    % Compare with previous checksum
    exit_loop_checksum=0;
    for previous_state=1:1:iteration_
        Checksum_previousstate_Label_lake = Checksum_list.state(previous_state).hash;
        check_checksum = isequal(Checksum_previousstate_Label_lake,Checksum_currentstate_Label_lake);
        if check_checksum==1
            exit_loop_checksum=1;
        end
    end
    % Save current checksum in the list
    Checksum_list.state(iteration_+1).hash=Checksum_currentstate_Label_lake;
    % Exit loop due to checksum equality
    if (exit_loop_checksum==1 && number_change>0)
        if details_convergence==1
            disp 'Iterative process has been founded to back to a previous state. Iterations are stopped.'
        end
        break
    end
    
end
if visualize_2D && cpsd_refining
    close(video_handle) % Close video
    function_savefig(Fig_illustration, [desktop '\'], 'Illustration of two-dimensional particle identification correction')
end

end

