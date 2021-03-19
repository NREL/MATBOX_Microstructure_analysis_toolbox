function [] = Function_connectivity(Phase_microstructure, PROPERTY, OPTIONS, INFO)
%Calculate connectivity
% Function_connectivity(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_connectivity(array, voxelsize) - when use as a standalone function

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 2 % Case for function called as: Function_particle_size_CPSD(Phase_microstructure, voxel_size)
        voxel_size = PROPERTY;
        clear PROPERTY OPTIONS
        
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        INFO.resultfoldername = ['Volume_characterization_' char(t)];
        desktop=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
        OPTIONS.mainsavefolder = [desktop '\' INFO.resultfoldername];           
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        
        % Set default phase information
        unique_code = unique(Phase_microstructure);
        INFO.number_phase = length(unique_code);
        INFO.voxel_number = numel(Phase_microstructure);
        for k=1:1:INFO.number_phase % Loop over all unique values
            INFO.phase(k).code = unique_code(k); % Assing code
            INFO.phasename(k,1) = {['Phase ' num2str(INFO.phase(k).code)]}; % Assign phase name based on code value
            INFO.phase(k).filename = ['Phase_' num2str(INFO.phase(k).code)]; % Assign phase name based on code value
            INFO.phase(k).name = INFO.phase(k).filename;
            INFO.phase(k).color = allcolor(k,:);
        end
        INFO.voxel_size = voxel_size; % nanometers
        
        % Set default direction information
        INFO.direction(1).name = 'in-plane direction 1';
        INFO.direction(2).name = 'in-plane direction 2';
        INFO.direction(3).name = 'through-plane direction';
        INFO.direction(1).filename = 'inplane_direction_1';
        INFO.direction(2).filename = 'inplane_direction_2';
        INFO.direction(3).filename = 'throughplane_direction';
        
        % Set default options
        OPTIONS.save_resultsmat = true;
        OPTIONS.save_xls = true;
        OPTIONS.save_fig = true;
        OPTIONS.savefig_infig = true;
        OPTIONS.savefig_informat = {'png'};
        
        % Set display options
        OPTIONS.fontname = 'Times New Roman';
        OPTIONS.displaytext = true;
        OPTIONS.closefigureaftercreation = false;
        OPTIONS.grid = 'on'; OPTIONS.minorgrid = 'on';
        OPTIONS.Linewidth = 2;
        OPTIONS.Fontsize_axe =  12;
        OPTIONS.Fontsize_legend =  12;
        OPTIONS.Fontsize_title =  14;
        
        % No Voxel size dependence analysis and RVE analysis
        PROPERTY.connectivity.voxel_size_dependence.todo = false;
        PROPERTY.connectivity.number_RVE = 0;
        PROPERTY.connectivity.todo_slice = [ {false}; {false}; {false}];
        PROPERTY.connectivity.sliceparameter = [{5};{5};{5}];
        PROPERTY.connectivity.slicechoice = 1;        
        
    else % Incorrect number of argument
        disp 'Error calling Function_connectivity. Wrong number of argument.'
        help Function_connectivity
    end
    
end


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Connectivity\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Connectivity/'; % Result are saved in this subfolder
end
Current_folder = [main_folder Sub_folder]; 
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% INFO
Domain_size = size(Phase_microstructure);
if length(Domain_size)==2 % Add third column with unit value
    Domain_size = [Domain_size 1];
end
if min(Domain_size) == 1 % 2D case
    number_dimension = 2;
else
    number_dimension =3; % 3D case
end
number_phase = INFO.number_phase; % Number of phase
voxel_number = INFO.voxel_number; % Number of voxel
voxel_size = INFO.voxel_size; % Voxel size nm

%% INITIALIZE RESULTS (USE FOR CORRELATION and VISUALIZATION)
for current_phase=1:1:number_phase
    results_correlation(current_phase).name = INFO.phase(current_phase).name;
    results_visualization(current_phase).name = INFO.phase(current_phase).name;
end

%% DISPLAY
if OPTIONS.displaytext==true
    disp '    CONNECTIVITY';
    disp '    ------------';
    disp ' ';
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

% Voxels and clusters are marked with:
connected_id = 1; 
isolated_id = 2;
unknown_id = 3;

%% CALCULATION
% Initialization
connectivity_illdefined = zeros(number_phase,2); % number voxels, volume fraction
connectivity_LargestIsolatedUnknown = zeros(number_phase,3); % largest, isolated, unknown
connectivity_TransportFace2Face = zeros(number_phase,9); % Connected, isolated, unknown for each direction
connectivity_TransportFromFace1 = zeros(number_phase,9); % Connected, isolated, unknown for each direction
connectivity_TransportFromFace2 = zeros(number_phase,9); % Connected, isolated, unknown for each direction


Time_measure=[];
for current_phase=1:1:number_phase % Loop over all phases
    time_cpu_start = cputime; % CPU start
    tic; % Stopwatch start
    
    code_tmp = INFO.phase(current_phase).code;
    binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3)); % Initialization
    binary_phase(Phase_microstructure == code_tmp) = 1;
    
    [Connectivity_structure] = Function_Connectivity_Algorithm(binary_phase, voxel_size, connected_id, unknown_id, isolated_id); % Call algorithm
    
    % Number of voxel of the current resized microstructure
    voxel_number_tmp=sum(sum(sum( binary_phase==1)));
    % CPU and stopwatch time - end
    time_cpu_elapsed = cputime-time_cpu_start;
    time_stopwatch_elapsed = toc;
    Time_tmp = [voxel_number_tmp time_cpu_elapsed time_stopwatch_elapsed];
    Time_measure = [Time_measure;Time_tmp];
    
    connectivity_illdefined(current_phase,1) = Connectivity_structure.illdetailed_cluster_number;
    connectivity_illdefined(current_phase,2) = Connectivity_structure.illdetailed_cluster_phasevolumefraction;
    
    connectivity_LargestIsolatedUnknown(current_phase,1) = Connectivity_structure.Clusters_LargestIsolatedUnknown.main_cluster_phasefraction;
    connectivity_LargestIsolatedUnknown(current_phase,2) = Connectivity_structure.Clusters_LargestIsolatedUnknown.isolated_cluster_phasefraction;
    connectivity_LargestIsolatedUnknown(current_phase,3) = Connectivity_structure.Clusters_LargestIsolatedUnknown.unknown_cluster_phasefraction;
    
    for direction=1:1:3
        connectivity_TransportFace2Face(current_phase,(direction-1)*3+1) = Connectivity_structure.Clusters_TransportFace2Face.direction(direction).connected_cluster_phasefraction;
        connectivity_TransportFace2Face(current_phase,(direction-1)*3+2) = Connectivity_structure.Clusters_TransportFace2Face.direction(direction).isolated_cluster_phasefraction;
        connectivity_TransportFace2Face(current_phase,(direction-1)*3+3) = Connectivity_structure.Clusters_TransportFace2Face.direction(direction).unknown_cluster_phasefraction;

        connectivity_TransportFromFace1(current_phase,(direction-1)*3+1) = Connectivity_structure.Clusters_TransportFromFace1.direction(direction).connected_cluster_phasefraction;
        connectivity_TransportFromFace1(current_phase,(direction-1)*3+2) = Connectivity_structure.Clusters_TransportFromFace1.direction(direction).isolated_cluster_phasefraction;
        connectivity_TransportFromFace1(current_phase,(direction-1)*3+3) = Connectivity_structure.Clusters_TransportFromFace1.direction(direction).unknown_cluster_phasefraction;

        connectivity_TransportFromFace2(current_phase,(direction-1)*3+1) = Connectivity_structure.Clusters_TransportFromFace2.direction(direction).connected_cluster_phasefraction;
        connectivity_TransportFromFace2(current_phase,(direction-1)*3+2) = Connectivity_structure.Clusters_TransportFromFace2.direction(direction).isolated_cluster_phasefraction;
        connectivity_TransportFromFace2(current_phase,(direction-1)*3+3) = Connectivity_structure.Clusters_TransportFromFace2.direction(direction).unknown_cluster_phasefraction;
    end
    
    % Save (correlation, visualization)
    results_correlation(current_phase).main_cluster_phasefraction = connectivity_LargestIsolatedUnknown(current_phase,1);
    results_correlation(current_phase).isolated_cluster_phasefraction = connectivity_LargestIsolatedUnknown(current_phase,2);
    results_correlation(current_phase).unknown_cluster_phasefraction = connectivity_LargestIsolatedUnknown(current_phase,3);
    
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_1 =  connectivity_TransportFace2Face(current_phase,1);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_1 = connectivity_TransportFace2Face(current_phase,2);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_1 = connectivity_TransportFace2Face(current_phase,3); 
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_2 =  connectivity_TransportFace2Face(current_phase,4);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_2 = connectivity_TransportFace2Face(current_phase,5);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_2 = connectivity_TransportFace2Face(current_phase,6);  
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_3 =  connectivity_TransportFace2Face(current_phase,7);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_3 = connectivity_TransportFace2Face(current_phase,8);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_3 = connectivity_TransportFace2Face(current_phase,9);      

    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_1_fromface1 =  connectivity_TransportFromFace1(current_phase,1);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_1_face1 = connectivity_TransportFromFace1(current_phase,2);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_1_face1 = connectivity_TransportFromFace1(current_phase,3); 
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_2_fromface1 =  connectivity_TransportFromFace1(current_phase,4);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_2_face1 = connectivity_TransportFromFace1(current_phase,5);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_2_face1 = connectivity_TransportFromFace1(current_phase,6);  
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_3_fromface1 =  connectivity_TransportFromFace1(current_phase,7);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_3_face1 = connectivity_TransportFromFace1(current_phase,8);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_3_face1 = connectivity_TransportFromFace1(current_phase,9);  

    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_1_fromface2 =  connectivity_TransportFromFace2(current_phase,1);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_1_face2 = connectivity_TransportFromFace2(current_phase,2);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_1_face2 = connectivity_TransportFromFace2(current_phase,3); 
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_2_fromface2 =  connectivity_TransportFromFace2(current_phase,4);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_2_face2 = connectivity_TransportFromFace2(current_phase,5);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_2_face2 = connectivity_TransportFromFace2(current_phase,6);  
    results_correlation(current_phase).Cluster_phasefraction_connected_along_direction_3_fromface2 =  connectivity_TransportFromFace2(current_phase,7);
    %results_correlation(current_phase).Cluster_phasefraction_isolated_along_direction_3_face2 = connectivity_TransportFromFace2(current_phase,8);
    %results_correlation(current_phase).Cluster_phasefraction_unknown_along_direction_3_face2 = connectivity_TransportFromFace2(current_phase,9);  
    
    results_visualization(current_phase).Clusters_sortedpersize = Connectivity_structure.Clusters_sortedpersize_3Darray;
    results_visualization(current_phase).Clusters_LargestIsolatedUnknown = Connectivity_structure.Clusters_LargestIsolatedUnknown.array;    
    results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_1 = Connectivity_structure.Clusters_TransportFace2Face.direction(1).array;        
    results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_2 = Connectivity_structure.Clusters_TransportFace2Face.direction(2).array;        
    results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_3 = Connectivity_structure.Clusters_TransportFace2Face.direction(3).array;      
    results_visualization(current_phase).Clusters_ConnectedFromFace1_along_direction_1 = Connectivity_structure.Clusters_TransportFromFace1.direction(1).array;        
    results_visualization(current_phase).Clusters_ConnectedFromFace1_along_direction_2 = Connectivity_structure.Clusters_TransportFromFace1.direction(2).array;        
    results_visualization(current_phase).Clusters_ConnectedFromFace1_along_direction_3 = Connectivity_structure.Clusters_TransportFromFace1.direction(3).array;     
    results_visualization(current_phase).Clusters_ConnectedFromFace2_along_direction_1 = Connectivity_structure.Clusters_TransportFromFace2.direction(1).array;        
    results_visualization(current_phase).Clusters_ConnectedFromFace2_along_direction_2 = Connectivity_structure.Clusters_TransportFromFace2.direction(2).array;        
    results_visualization(current_phase).Clusters_ConnectedFromFace2_along_direction_3 = Connectivity_structure.Clusters_TransportFromFace2.direction(3).array;   
end


%% TABLES
% Time
Table_time = table(Time_measure(:,1)*1e-6,Time_measure(:,2),Time_measure(:,3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_connectivity.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_connectivity = table(INFO.phasename,connectivity_LargestIsolatedUnknown(:,1),connectivity_LargestIsolatedUnknown(:,2),connectivity_LargestIsolatedUnknown(:,3),...
    'VariableNames',{'Name' 'Main_cluster' 'Isolated_clusters' 'Unknown_clusters'});
Table_face2face = table(INFO.phasename,connectivity_TransportFace2Face(:,1),connectivity_TransportFace2Face(:,2),connectivity_TransportFace2Face(:,3),connectivity_TransportFace2Face(:,4),connectivity_TransportFace2Face(:,5),connectivity_TransportFace2Face(:,6),connectivity_TransportFace2Face(:,7),connectivity_TransportFace2Face(:,8),connectivity_TransportFace2Face(:,9),...
    'VariableNames',{'Name' 'Dir1_connected' 'Dir1_isolated' 'Dir1_unknown' 'Dir2_connected' 'Dir2_isolated' 'Dir2_unknown' 'Dir3_connected' 'Dir3_isolated' 'Dir3_unknown'});
Table_fromface1 = table(INFO.phasename,connectivity_TransportFromFace1(:,1),connectivity_TransportFromFace1(:,2),connectivity_TransportFromFace1(:,3),connectivity_TransportFromFace1(:,4),connectivity_TransportFromFace1(:,5),connectivity_TransportFromFace1(:,6),connectivity_TransportFromFace1(:,7),connectivity_TransportFromFace1(:,8),connectivity_TransportFromFace1(:,9),...
    'VariableNames',{'Name' 'Dir1_connected' 'Dir1_isolated' 'Dir1_unknown' 'Dir2_connected' 'Dir2_isolated' 'Dir2_unknown' 'Dir3_connected' 'Dir3_isolated' 'Dir3_unknown'});
Table_fromface2 = table(INFO.phasename,connectivity_TransportFromFace2(:,1),connectivity_TransportFromFace2(:,2),connectivity_TransportFromFace2(:,3),connectivity_TransportFromFace2(:,4),connectivity_TransportFromFace2(:,5),connectivity_TransportFromFace2(:,6),connectivity_TransportFromFace2(:,7),connectivity_TransportFromFace2(:,8),connectivity_TransportFromFace2(:,9),...
    'VariableNames',{'Name' 'Dir1_connected' 'Dir1_isolated' 'Dir1_unknown' 'Dir2_connected' 'Dir2_isolated' 'Dir2_unknown' 'Dir3_connected' 'Dir3_isolated' 'Dir3_unknown'});
Table_illdefinedclusters = table(INFO.phasename,connectivity_illdefined(:,1),connectivity_illdefined(:,2),...
    'VariableNames',{'Name' 'Number_1_voxel_size_cluster' 'Phase_volume_fraction'});
Results_connectivity.Table_connectivity = Table_connectivity; % Save in main table result
Results_connectivity.Table_face2face = Table_face2face;
Results_connectivity.Table_fromface1 = Table_fromface1;
Results_connectivity.Table_fromface2 = Table_fromface2;
Results_connectivity.Table_illdefinedclusters = Table_illdefinedclusters;

%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Connectivity'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Results';
    DATA_writetable.sheet(1).table=Table_connectivity;
    DATA_writetable.sheet(2).name='Face2Face';
    DATA_writetable.sheet(2).table=Table_face2face;  
    DATA_writetable.sheet(3).name='FromFaceStart';
    DATA_writetable.sheet(3).table=Table_fromface1;  
    DATA_writetable.sheet(4).name='FromFaceEnd';
    DATA_writetable.sheet(4).table=Table_fromface2;        
    DATA_writetable.sheet(5).name='Illdefinedclusters';
    DATA_writetable.sheet(5).table=Table_illdefinedclusters;       
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
if OPTIONS.displaytext==true
    fprintf('> Calculated on the whole domain:\n\n');
    disp 'Connectivity, % of the phase';
    disp(Table_connectivity)
    disp 'Connectivity from extremes faces, % of phase';
    disp(Table_face2face)
    disp 'Connectivity from first face, % of phase';
    disp(Table_fromface1)
    disp 'Connectivity from last face, % of phase';
    disp(Table_fromface2)
    disp 'Ill defined cluster (one-voxel size)';
    disp(Table_illdefinedclusters)    
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_time)
end


%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%
scrsz = get(0,'ScreenSize'); % Screen resolution

%% CONNECTIVITY MAP
for current_phase=1:1:number_phase % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Connectivity, phase ' INFO.phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]); % Full screen figure
    id_axe=0;
    for view_direction=1:1:number_dimension % Iterate over axe
        for k=1:1:4
            if k==1
                str_='Volume connectivity';
                data_tmp = results_visualization(current_phase).Clusters_LargestIsolatedUnknown;
            else
                str_=sprintf('Connectivity along direction %i',k-1);
                if k==2
                    data_tmp = results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_1;
                elseif k==3
                    data_tmp = results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_2;
                elseif k==4
                    data_tmp = results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_3;
                end
            end
            id_axe=id_axe+1;
            sub_axes=subplot(number_dimension,4,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            h_title=title ({str_,'Slice in the middle',['View normal to ' INFO.direction(view_direction).name]}); % Set title font
            % - Plot graphs
            if view_direction==1
                tmp=squeeze(data_tmp(round(Domain_size(1)/2),:,:));
                h=image(tmp,'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',INFO.direction(3).name);
                t_1y = sprintf('Position along %s ',INFO.direction(2).name);
                set(h, 'XData', [0, Domain_size(3)*voxel_size/1000]);
                set(h, 'YData', [0, Domain_size(2)*voxel_size/1000]);
            elseif view_direction==2
                tmp=squeeze(data_tmp(:,round(Domain_size(2)/2),:));
                h=image(tmp,'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',INFO.direction(3).name);
                t_1y = sprintf('Position along %s ',INFO.direction(1).name);
                set(h, 'XData', [0, Domain_size(3)*voxel_size/1000]);
                set(h, 'YData', [0, Domain_size(1)*voxel_size/1000]);                
            elseif view_direction==3
                h=image(data_tmp(:,:,round(Domain_size(3)/2)),'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',INFO.direction(1).name);
                t_1y = sprintf('Position along %s ',INFO.direction(2).name);
                set(h, 'XData', [0, Domain_size(1)*voxel_size/1000]);
                set(h, 'YData', [0, Domain_size(2)*voxel_size/1000]);                
            end
            axis equal; axis tight;
            % - Axis label
            t_ = xlabel(' ');
            t_2 = '(\mum)';
            t_.String= [t_1x t_2]; % Sprintf does not accept greek characters
            t_ = ylabel(' ');
            t_2 = '(\mum)';
            t_.String= [t_1y t_2]; % Sprintf does not accept greek characters
            set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
            % Color map
            myColorMap = jet(256);
            myColorMap(1,:) = 1; % background color is white
            colormap(myColorMap);
            caxis([0 3])
            % Create colorbar
            %colorbar('peer',sub_axes);
            h=colorbar(sub_axes);
            set(h,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe);
            h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
    end
    sgtitle(Fig,['Cluster connectivity, connected(1), isolated(2), and unknown(3), ' INFO.phase(current_phase).name],'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
    if OPTIONS.save_fig == true % Save figure
        filename= sprintf('View_connectivity_%s', INFO.phase(current_phase).filename);
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end
clear tmp data_tmp;

%% ALONG DIRECTIONS
% Initialize
for current_direction=1:1:number_dimension % Loop over all directions
    tmp = Domain_size; tmp(current_direction)=1; Results.direction(current_direction).numbervoxelslice = prod(tmp); clear tmp; % Number of voxel within each slice normal to the direction
    Results.direction(current_direction).connectivity = zeros(Domain_size(current_direction),number_phase+1,3); % Initialization
    Results.direction(current_direction).connectivity_F2F1 = zeros(Domain_size(current_direction),number_phase+1,3); % Initialization
    Results.direction(current_direction).connectivity_F2F2 = zeros(Domain_size(current_direction),number_phase+1,3); % Initialization
    Results.direction(current_direction).connectivity_F2F3 = zeros(Domain_size(current_direction),number_phase+1,3); % Initialization
    for k=1:1:3
        Results.direction(current_direction).connectivity(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size/1000)'; % Axis (position along direction)
        Results.direction(current_direction).connectivity_F2F1(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size/1000)'; % Axis (position along direction)
        Results.direction(current_direction).connectivity_F2F2(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size/1000)'; % Axis (position along direction)
        Results.direction(current_direction).connectivity_F2F3(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size/1000)'; % Axis (position along direction)
    end
end
% Calculate results
for current_phase = 1:1:number_phase % Loop over all phases
    Vtmp = results_visualization(current_phase).Clusters_LargestIsolatedUnknown;
    F1tmp = results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_1;
    F2tmp =results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_2;
    F3tmp = results_visualization(current_phase).Clusters_ConnectedFace2Face_along_direction_3;
    for c=1:1:3
        % Direction 1
        for current_position=1:1:Domain_size(1) % Loop over postion
            n= sum(sum(Vtmp(current_position,:,:)==c)); phasevoxel = sum(sum(Vtmp(current_position,:,:)>0));
            Results.direction(1).connectivity(current_position,current_phase+1,c) = n / phasevoxel;
            n= sum(sum(F1tmp(current_position,:,:)==c)); 
            Results.direction(1).connectivity_F2F1(current_position,current_phase+1,c) = n / phasevoxel;       
            n= sum(sum(F2tmp(current_position,:,:)==c));
            Results.direction(1).connectivity_F2F2(current_position,current_phase+1,c) = n / phasevoxel;     
            n= sum(sum(F3tmp(current_position,:,:)==c));
            Results.direction(1).connectivity_F2F3(current_position,current_phase+1,c) = n / phasevoxel;                 
        end
        % Direction 2
        for current_position=1:1:Domain_size(2) % Loop over postion
            n= sum(sum(Vtmp(:,current_position,:)==c));  phasevoxel = sum(sum(Vtmp(:,current_position,:)>0));
            Results.direction(2).connectivity(current_position,current_phase+1,c) = n / phasevoxel;
            n= sum(sum(F1tmp(:,current_position,:)==c));
            Results.direction(2).connectivity_F2F1(current_position,current_phase+1,c) = n / phasevoxel;    
            n= sum(sum(F2tmp(:,current_position,:)==c));
            Results.direction(2).connectivity_F2F2(current_position,current_phase+1,c) = n / phasevoxel;      
            n= sum(sum(F3tmp(:,current_position,:)==c));
            Results.direction(2).connectivity_F2F3(current_position,current_phase+1,c) = n / phasevoxel;                  
        end
        % Direction 3
        for current_position=1:1:Domain_size(3) % Loop over postion
            n= sum(sum(Vtmp(:,:,current_position)==c));  phasevoxel = sum(sum(Vtmp(:,:,current_position)>0));
            Results.direction(3).connectivity(current_position,current_phase+1,c) = n / phasevoxel;
            n= sum(sum(F1tmp(:,:,current_position)==c));
            Results.direction(3).connectivity_F2F1(current_position,current_phase+1,c) = n / phasevoxel;      
            n= sum(sum(F2tmp(:,:,current_position)==c));
            Results.direction(3).connectivity_F2F2(current_position,current_phase+1,c) = n / phasevoxel;      
            n= sum(sum(F3tmp(:,:,current_position)==c));
            Results.direction(3).connectivity_F2F3(current_position,current_phase+1,c) = n / phasevoxel;                  
        end
    end
end
clear Vtmp F1tmp F2tmp F3tmp

%% TABLES

% Prepare header name
clear Variable_name_table;
Variable_name_table(1)={'Position_um'};
for current_phase=1:1:number_phase
    Variable_name_table(current_phase+1)={INFO.phase(current_phase).filename};
end
% Create table
for c=1:1:3
    for current_direction=1:1:number_dimension % Loop over all directions
        EvolutionConnectivity.type(1).direction(current_direction).connectivity(c).table = array2table(Results.direction(current_direction).connectivity(:,:,c),'VariableNames',Variable_name_table);
        EvolutionConnectivity.type(2).direction(current_direction).connectivity(c).table = array2table(Results.direction(current_direction).connectivity_F2F1(:,:,c),'VariableNames',Variable_name_table);
        EvolutionConnectivity.type(3).direction(current_direction).connectivity(c).table = array2table(Results.direction(current_direction).connectivity_F2F2(:,:,c),'VariableNames',Variable_name_table);
        EvolutionConnectivity.type(4).direction(current_direction).connectivity(c).table = array2table(Results.direction(current_direction).connectivity_F2F3(:,:,c),'VariableNames',Variable_name_table);
    end
end
Results_connectivity.EvolutionConnectivity = EvolutionConnectivity.type(1); % Save in main table result
Results_connectivity.EvolutionConnectivity_F2F1 = EvolutionConnectivity.type(2); % Save in main table result
Results_connectivity.EvolutionConnectivity_F2F2 = EvolutionConnectivity.type(3); % Save in main table result
Results_connectivity.EvolutionConnectivity_F2F3 = EvolutionConnectivity.type(4); % Save in main table result

%% SAVE TABLES
if OPTIONS.save_xls==true
    for current_direction=1:1:number_dimension % Loop over all directions
        filename = ['Connectivity_along_d' num2str(current_direction) '_' INFO.direction(current_direction).filename]; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheer_number = 0;
        for connectivity_type=1:1:4 % Loop over all connectivity type
            if connectivity_type==1
                str_a='Volume';
            elseif connectivity_type==2
                str_a='Face2Face_normal_1';
            elseif connectivity_type==3
                str_a='Face2Face_normal_2';
            elseif connectivity_type==4
                str_a='Face2Face_normal_3';
            end
            for c=1:1:3 % Loop over connectivity id
                if c==1
                    str_b='Connected';
                elseif c==2
                    str_b='Isolated';
                elseif c==3
                    str_b='Unknown';
                end
                sheer_number=sheer_number+1;
                DATA_writetable.sheet(sheer_number).name = [str_a '_' str_b];
                DATA_writetable.sheet(sheer_number).table = EvolutionConnectivity.type(connectivity_type).direction(current_direction).connectivity(c).table;
            end
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
end

%% FIGURES
for current_phase=1:1:number_phase
    Fig = figure; % Create figure
    Fig.Name= ['Connectivity along directions, ' INFO.phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]); % Full screen figure
    id_axe=0;
    for current_direction=1:1:number_dimension % Iterate over axe
        for c=1:1:3 % Loop over connectivity id
            if c==1
                str_id='Connected cluster(s)';
            elseif c==2
                str_id='Isolated clusters';
            elseif c==3
                str_id='Unknown clusters';
            end
            id_axe=id_axe+1;
            sub_axes=subplot(number_dimension,3,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            h_title=title (' '); % Set title font
            h_title.String = {['Connectivity along ' INFO.direction(current_direction).name],str_id}; % Set title font
            % Plot graphs
            for connectivity_type=1:1:4 % Loop over phases
                if connectivity_type==1
                    x_ = Results.direction(current_direction).connectivity(:,1,c);
                    y_ = Results.direction(current_direction).connectivity(:,current_phase+1,c);
                elseif connectivity_type==2
                    x_ = Results.direction(current_direction).connectivity_F2F1(:,1,c);
                    y_ = Results.direction(current_direction).connectivity_F2F1(:,current_phase+1,c);
                elseif connectivity_type==3
                    x_ = Results.direction(current_direction).connectivity_F2F2(:,1,c);
                    y_ = Results.direction(current_direction).connectivity_F2F2(:,current_phase+1,c);
                elseif connectivity_type==4
                    x_ = Results.direction(current_direction).connectivity_F2F3(:,1,c);                    
                    y_ = Results.direction(current_direction).connectivity_F2F3(:,current_phase+1,c);
                end
                h=plot(x_,100*y_); % Percents
                set(h, 'LineWidth',OPTIONS.Linewidth); % Colors
            end
            % Axis label
            t_ = xlabel(' ');
            t_1 = sprintf('Position along %s ',INFO.direction(current_direction).name);
            t_2 = '(\mum)';
            t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
            ylabel('Connnectivity (%)');
            % Legend
            h_legend = legend(sub_axes,{'Volume ',['Face2Face normal ' INFO.direction(1).name],['Face2Face normal ' INFO.direction(2).name],['Face2Face normal ' INFO.direction(3).name]},'Location','best');
            % - Grid
            if strcmp(OPTIONS.grid,'on')
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
            end
            set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
            h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
            h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
    end
    sgtitle(Fig,Fig.Name,'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
    if OPTIONS.save_fig == true % Save figure
        filename= Fig.Name;
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end

%% ALONG DIRECTIONS (THICK SLICE PER THICK SLICE)
for current_direction=1:1:number_dimension
    if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
        if PROPERTY.connectivity.slicechoice==1
            number_thickslice = round(cell2mat(PROPERTY.connectivity.sliceparameter(current_direction)));
        else
            thickness_slice = cell2mat(PROPERTY.connectivity.sliceparameter(current_direction));
            number_thickslice = round(Domain_size(current_direction)*voxel_size/1000 / thickness_slice);
        end        
        number_thickslice = max([1 number_thickslice]);
        position_slice = round(linspace(1,Domain_size(current_direction),number_thickslice+1));
        for current_phase=1:1:number_phase
            code_tmp = INFO.phase(current_phase).code; % The code of the phase
            binary_phase = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
            binary_phase( Phase_microstructure==code_tmp )=1;
            
            val_ = zeros(number_thickslice,15); % Start position, end position, porosity,  Connected, isolated, unknow (x4)
            for current_slice=1:1:number_thickslice
                if current_direction==1
                    x0 = position_slice(current_slice); x1 = position_slice(current_slice+1);
                    y0 = 1; y1 = Domain_size(2);
                    z0 = 1; z1 = Domain_size(3);
                elseif current_direction==2
                    x0 = 1; x1 = Domain_size(1);
                    y0 = position_slice(current_slice); y1 = position_slice(current_slice+1);
                    z0 = 1; z1 = Domain_size(3);
                elseif current_direction==3
                    x0 = 1; x1 = Domain_size(1);
                    y0 = 1; y1 = Domain_size(2);
                    z0 = position_slice(current_slice); z1 = position_slice(current_slice+1);
                end
                val_(current_slice,1) = position_slice(current_slice);
                val_(current_slice,2) = position_slice(current_slice+1);
                thick_slice = binary_phase(x0:x1,y0:y1,z0:z1);
                
                voxel_number = sum(sum(sum(thick_slice==1)));
                val_(current_slice,3) = voxel_number/numel(thick_slice);
                
                time_cpu_start = cputime; % CPU time start
                tic; % Stopwatch start
                
                [Connectivity_structure_tmp] = Function_Connectivity_Algorithm(thick_slice, voxel_size, connected_id, unknown_id, isolated_id); % Call algorithm
                val_(current_slice,4) = Connectivity_structure_tmp.Clusters_LargestIsolatedUnknown.main_cluster_phasefraction;
                val_(current_slice,5) = Connectivity_structure_tmp.Clusters_LargestIsolatedUnknown.isolated_cluster_phasefraction;
                val_(current_slice,6) = Connectivity_structure_tmp.Clusters_LargestIsolatedUnknown.unknown_cluster_phasefraction;
                
                val_(current_slice,7) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(1).connected_cluster_phasefraction;
                val_(current_slice,8) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(1).isolated_cluster_phasefraction;
                val_(current_slice,9) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(1).unknown_cluster_phasefraction;     
                
                val_(current_slice,10) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(2).connected_cluster_phasefraction;
                val_(current_slice,11) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(2).isolated_cluster_phasefraction;
                val_(current_slice,12) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(2).unknown_cluster_phasefraction;     
                
                val_(current_slice,13) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(3).connected_cluster_phasefraction;
                val_(current_slice,14) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(3).isolated_cluster_phasefraction;
                val_(current_slice,15) = Connectivity_structure_tmp.Clusters_TransportFace2Face.direction(3).unknown_cluster_phasefraction;                
                
                % CPU and stopwatch time - end
                time_cpu_elapsed = cputime-time_cpu_start;
                time_stopwatch_elapsed = toc;
                Time_tmp = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
                Time_measure = [Time_measure;Time_tmp];
            end
            % Convert to micrometers
            val_(:,1) = val_(:,1)*voxel_size/1000;
            val_(:,2) = val_(:,2)*voxel_size/1000;
            val_(1,1)=0;
            % Table
            tmp_ = table(val_(:,1),val_(:,2),val_(:,3),val_(:,4),val_(:,5),val_(:,6),val_(:,7),val_(:,8),val_(:,9),val_(:,10),val_(:,11),val_(:,12),val_(:,13),val_(:,14),val_(:,15),...
                'VariableNames',{'Start' 'End' 'Porosity' 'Connected_volume' 'Isolated_volume' 'Unknown_volume' 'Connected_Face2Face_normal_1' 'Isolated_Face2Face_normal_1' 'Unknown_Face2Face_normal_1'  'Connected_Face2Face_normal_2' 'Isolated_Face2Face_normal_2' 'Unknown_Face2Face_normal_2' 'Connected_Face2Face_normal_3' 'Isolated_Face2Face_normal_3' 'Unknown_Face2Face_normal_3'});
            Connectivityslice.direction(current_direction).phase(current_phase).table = tmp_;
        end
    end
end

if logical(cell2mat(PROPERTY.connectivity.todo_slice(1))) || logical(cell2mat(PROPERTY.connectivity.todo_slice(2))) || logical(cell2mat(PROPERTY.connectivity.todo_slice(3)))
    Results_connectivity.Table_Connectivity_slice = Connectivityslice; % Save in main table result
    
    % SAVE TABLES
    if OPTIONS.save_xls==true
        filename = 'Connectivity_slice'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheet_=0;
        for current_phase=1:1:number_phase
            for current_direction=1:1:number_dimension
                if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
                    sheet_=sheet_+1;
                    DATA_writetable.sheet(sheet_).name=[INFO.phase(current_phase).filename ' ' INFO.direction(current_direction).filename];
                    DATA_writetable.sheet(sheet_).table=Results_connectivity.Table_Connectivity_slice.direction(current_direction).phase(current_phase).table;
                end
            end
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    % FIGURES
    Color_order = get(0, 'DefaultAxesColorOrder');
    for current_phase=1:1:number_phase
        code_tmp = INFO.phase(current_phase).code; % Phase code
        volumefraction_phase =  sum(sum(sum(Phase_microstructure==code_tmp)))/numel(Phase_microstructure);
        clear str_legend
        Fig = figure; % Create figure
        Fig.Name= ['Connectivity per thick slice, ' INFO.phase(current_phase).name];
        Fig.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig,'position',scrsz); % Full screen figure
        % - Create axes as a subplot
        for id_axe=1:1:5
            k_legend=0;
            sub_axes=subplot(2,3,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            if id_axe==1
                h_title=title ('Volume fraction \epsilon'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_connectivity.Table_Connectivity_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,3);
                            y_(2*k,1)   = t_(k,3);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [volumefraction_phase volumefraction_phase]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(volumefraction_phase,'%1.3f')];
                    end
                end
                ylabel('Volume fraction \epsilon');
            elseif id_axe==2
                h_title=title ('Connectivity: volume'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_connectivity.Table_Connectivity_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,4);
                            y_(2*k,1)   = t_(k,4);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [connectivity_LargestIsolatedUnknown(current_phase,1) connectivity_LargestIsolatedUnknown(current_phase,1)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(connectivity_LargestIsolatedUnknown(current_phase,1),'%1.3f')];
                    end
                    ylabel('Connectivity (%)');
                end
            elseif id_axe==3
                h_title=title (['Connectivity: Face2Face normal ' INFO.direction(1).name]); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_connectivity.Table_Connectivity_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,7);
                            y_(2*k,1)   = t_(k,7);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [connectivity_TransportFace2Face(current_phase,(1-1)*3+1) connectivity_TransportFace2Face(current_phase,(1-1)*3+1)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(connectivity_TransportFace2Face(current_phase,(1-1)*3+1),'%1.3f')];
                    end
                    ylabel('Connectivity (%)');
                end
            elseif id_axe==4
                h_title=title (['Connectivity: Face2Face normal ' INFO.direction(2).name]); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_connectivity.Table_Connectivity_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,10);
                            y_(2*k,1)   = t_(k,10);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [connectivity_TransportFace2Face(current_phase,(2-1)*3+1) connectivity_TransportFace2Face(current_phase,(2-1)*3+1)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(connectivity_TransportFace2Face(current_phase,(2-1)*3+1),'%1.3f')];
                    end
                    ylabel('Connectivity (%)');
                end
            elseif id_axe==5
                h_title=title (['Connectivity: Face2Face normal ' INFO.direction(3).name]); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.connectivity.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_connectivity.Table_Connectivity_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,13);
                            y_(2*k,1)   = t_(k,13);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [connectivity_TransportFace2Face(current_phase,(3-1)*3+1) connectivity_TransportFace2Face(current_phase,(3-1)*3+1)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(connectivity_TransportFace2Face(current_phase,(3-1)*3+1),'%1.3f')];
                    end
                    ylabel('Connectivity (%)');
                end
            end
            xlabel('Position along domain''s direction (\mum)');
            h_legend = legend(sub_axes,str_legend.name,'Location','best');
            % - Grid
            if strcmp(OPTIONS.grid,'on')
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
            end
            set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
            h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
            h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
        sgtitle(Fig,['Connectivity along directions, ' INFO.phase(current_phase).name] ,'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
        if OPTIONS.save_fig == true % Save figure
            function_savefig(Fig, Current_folder, ['Connectivity_slice_' INFO.phase(current_phase).filename], OPTIONS); % Call function
        end
        if OPTIONS.closefigureaftercreation == true
            close(Fig); % Do not keep open figures
        end
    end    
end


%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%
interpolation_voxelsize_order=1;

if PROPERTY.connectivity.voxel_size_dependence.todo % Check if voxel size analysis is asked
    size_choice = PROPERTY.connectivity.voxel_size_dependence.voxel;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed

    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_phase+1,5);
    property_voxelsizedependence(1,1,:)=voxel_size;
    property_voxelsizedependence(1,2:end,1)=connectivity_LargestIsolatedUnknown(:,1)';
    property_voxelsizedependence(1,2:end,2)=connectivity_TransportFace2Face(:,(1-1)*3+1)';
    property_voxelsizedependence(1,2:end,3)=connectivity_TransportFace2Face(:,(2-1)*3+1)';
    property_voxelsizedependence(1,2:end,4)=connectivity_TransportFace2Face(:,(3-1)*3+1)';
    property_voxelsizedependence(1,2:end,5)=connectivity_illdefined(:,2)';
    for current_iteration=1:1:number_resize
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1,:)=current_voxel_size;
        % Microstructure resized
        [Phase_microstructure_resized] = function_scale_array(Phase_microstructure, voxel_size, current_voxel_size, INFO.phaseinfo);        
        Current_domain_size=size(Phase_microstructure_resized);

        for current_phase=1:1:number_phase
            % CPU and stopwatch time - start
            time_cpu_start = cputime;
            tic;
            
            code_tmp = INFO.phase(current_phase).code; % Code of the phase
            binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Binary phase
            binary_phase(Phase_microstructure_resized == code_tmp) = 1;
            
            [Connectivity_structure_resized] = Function_Connectivity_Algorithm(binary_phase, current_voxel_size, connected_id, unknown_id, isolated_id); % Call algorithm

            property_voxelsizedependence(current_iteration+1,current_phase+1,1)=Connectivity_structure_resized.Clusters_LargestIsolatedUnknown.main_cluster_phasefraction;
            property_voxelsizedependence(current_iteration+1,current_phase+1,2)=Connectivity_structure_resized.Clusters_TransportFace2Face.direction(1).connected_cluster_phasefraction;
            property_voxelsizedependence(current_iteration+1,current_phase+1,3)=Connectivity_structure_resized.Clusters_TransportFace2Face.direction(2).connected_cluster_phasefraction;
            property_voxelsizedependence(current_iteration+1,current_phase+1,4)=Connectivity_structure_resized.Clusters_TransportFace2Face.direction(3).connected_cluster_phasefraction;
            property_voxelsizedependence(current_iteration+1,current_phase+1,5)=Connectivity_structure_resized.illdetailed_cluster_phasevolumefraction;
            
            % Number of voxel of the current resized microstructure
            voxel_number_tmp=sum(sum(sum( binary_phase==1)));
            % CPU and stopwatch time - end
            time_cpu_elapsed = cputime-time_cpu_start;
            time_stopwatch_elapsed = toc;
            Time_tmp = [voxel_number_tmp time_cpu_elapsed time_stopwatch_elapsed];
            Time_measure = [Time_measure;Time_tmp];            
        end

    end
    clear Phase_microstructure_resized Connectivity_structure_resized;
    
    %% EXTRAPOLATION TO 0 nm
    str_correlation(1).name = 'main_cluster_phasefraction';
    str_correlation(2).name = 'Cluster_phasefraction_connected_along_direction_1';
    str_correlation(3).name = 'Cluster_phasefraction_connected_along_direction_2';
    str_correlation(4).name = 'Cluster_phasefraction_connected_along_direction_3';
    tmp = zeros(number_resize+2,number_phase+1,5);
    for property_resized = 1:1:5
        x=property_voxelsizedependence(:,1,property_resized);
        for current_phase=1:1:number_phase
            y=property_voxelsizedependence(:,current_phase+1,property_resized);
            p = polyfit(x,y,interpolation_voxelsize_order);
            vq = polyval(p,0);
            tmp(1,current_phase+1,property_resized)=vq;
            interpolation_voxelsize(current_phase,property_resized).p=p;
            % For correlation
            if property_resized<5
                results_correlation(current_phase).([str_correlation(property_resized).name '_extrapolated']) = vq;
            end
        end
        tmp(2:end,:,property_resized) = property_voxelsizedependence(:,:,property_resized);
    end
    property_voxelsizedependence = tmp; clear tmp;
        
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={'Voxel_size_nm'}; % Columns name
    for current_phase=1:1:number_phase
        Variable_name_table(1+current_phase)={INFO.phase(current_phase).filename};
    end
    % Table
    Table_maincluster_voxelsizedependence = array2table(property_voxelsizedependence(:,:,1),...
        'VariableNames',Variable_name_table);
    Table_illdefinedcluster_voxelsizedependence = array2table(property_voxelsizedependence(:,:,5),...
        'VariableNames',Variable_name_table);
    
    k_=0;
    for current_direction=1:1:number_dimension
        for current_phase=1:1:number_phase
            k_=k_+1;
            Variable_name_table(1+k_)={[INFO.phase(current_phase).filename '_' INFO.direction(current_direction).filename]};
        end
    end    
    % table
    Table_connectivity_Face2Face_voxelsizedependence = array2table([property_voxelsizedependence(:,:,2) property_voxelsizedependence(:,2:end,3) property_voxelsizedependence(:,2:end,4)],...
        'VariableNames',Variable_name_table);
    
    %% DISPLAY TEXT RESULTS
    if (OPTIONS.displaytext==1)
        fprintf('> Connectivity dependence with the voxel size:\n\n');
        disp(Table_maincluster_voxelsizedependence)
        fprintf('> Connectivity from extremes faces dependence with the voxel size:\n\n');
        disp(Table_connectivity_Face2Face_voxelsizedependence)        
        fprintf('> Ill defined cluster dependence with the voxel size:\n\n');
        disp(Table_illdefinedcluster_voxelsizedependence)        
    end
        
    %% SAVE RESULTS
    Results_connectivity.voxelsizedependence_maincluster = Table_maincluster_voxelsizedependence; % Save in main table result
    Results_connectivity.voxelsizedependence_Face2Faceconnectivity = Table_connectivity_Face2Face_voxelsizedependence; % Save in main table result    
    Results_connectivity.voxelsizedependence_illdefined = Table_illdefinedcluster_voxelsizedependence; % Save in main table result    
    if OPTIONS.save_xls==true
        filename = 'Connectivity_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Main_cluster';
        DATA_writetable.sheet(1).table=Table_maincluster_voxelsizedependence;
        DATA_writetable.sheet(2).name='Face2Face_connectivity';
        DATA_writetable.sheet(2).table=Table_connectivity_Face2Face_voxelsizedependence;     
        DATA_writetable.sheet(3).name='Illdefined_clusters';
        DATA_writetable.sheet(3).table=Table_illdefinedcluster_voxelsizedependence;          
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end        
        
    %% FIGURES
    parameters_figure.number_phase = number_phase;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.INFO = INFO;
    parameters_figure.OPTIONS = OPTIONS;
    
    parameters_figure.propertyname = 'Main cluster connectivity';
    parameters_figure.method = '';
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,1);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,1);
    parameters_figure.str_ylabel = 'Connectivity (%)';
    parameters_figure.propertynameunit = '(%)';
    parameters_figure.filename = 'Main_cluster_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures

    parameters_figure.propertyname = 'Ill defined cluster';
    parameters_figure.method = '';
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,5);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,5);
    parameters_figure.str_ylabel = 'Volume fraction relative to phase (%)';
    parameters_figure.propertynameunit = '(%)';
    parameters_figure.filename = 'Illdefinedclusters_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    
    
    parameters_figure.propertyname = 'Face2face connectivity';
    parameters_figure.method = '';
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,2);
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,2);
    parameters_figure.interpolation_voxelsize2 = interpolation_voxelsize(:,3);
    parameters_figure.property_voxelsizedependence2 = property_voxelsizedependence(:,:,3);
    parameters_figure.interpolation_voxelsize3 = interpolation_voxelsize(:,4);
    parameters_figure.property_voxelsizedependence3 = property_voxelsizedependence(:,:,4);
    parameters_figure.str_ylabel = 'Connectivity (%)';
    parameters_figure.propertynameunit = '(%)';
    parameters_figure.filename = 'Face2Face_connectivity_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    
end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) ANALYSIS
%%

% Note: direction 4 is used to store volume connectivty. Directions 1,2,3
% correspond to face2face connectivity along their respective direction.

if PROPERTY.connectivity.number_RVE>0
    for n_RVE=1:1:PROPERTY.connectivity.number_RVE % Loop over all RVE asked
        RVEparameters.name = PROPERTY.connectivity.RVE(n_RVE).name;
        RVEparameters.savename = PROPERTY.connectivity.RVE(n_RVE).savename;
        RVEparameters.type = PROPERTY.connectivity.RVE(n_RVE).type;
        RVEparameters.divisions = PROPERTY.connectivity.RVE(n_RVE).divisions;
        RVEparameters.subs2 = PROPERTY.connectivity.RVE(n_RVE).subs2;
        RVEparameters.subs4 = PROPERTY.connectivity.RVE(n_RVE).subs4;
        RVEparameters.Aspectratio = PROPERTY.connectivity.RVE(n_RVE).Aspectratio;
        if  strcmp(PROPERTY.connectivity.RVE(n_RVE).type,'A')
            RVEparameters.Aspectratio_name = [num2str(Domain_size(1)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(2)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(3)/Domain_size(3),'%1.3f\t')];
        elseif strcmp(PROPERTY.connectivity.RVE(n_RVE).type,'B') || strcmp(PROPERTY.connectivity.RVE(n_RVE).type,'D')
            RVEparameters.Aspectratio_name = [num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
        end
        RVEparameters.Constantdirection = PROPERTY.connectivity.RVE(n_RVE).Constantdirection;
        RVEparameters.Growthdirection = PROPERTY.connectivity.RVE(n_RVE).Growthdirection;
        RVEparameters.Growthperstep = PROPERTY.connectivity.RVE(n_RVE).Growthperstep;
        RVEparameters.Growthrelativeto = PROPERTY.connectivity.RVE(n_RVE).Growthrelativeto;
        RVEparameters.threshold_std = PROPERTY.connectivity.RVE(n_RVE).threshold_std;
        RVEparameters.threshold_numbersubvolumes = PROPERTY.connectivity.RVE(n_RVE).threshold_numbersubvolumes;
        RVEparameters.firstuniquevolume_size = PROPERTY.connectivity.RVE(n_RVE).firstuniquevolume_size;
        RVEparameters.firstuniquevolume_unit = PROPERTY.connectivity.RVE(n_RVE).firstuniquevolume_unit;
        
        if OPTIONS.save_xls==true || OPTIONS.save_fig == true
            if ispc
                Sub_folder_RVE = [Current_folder RVEparameters.savename '\'];
            else
                Sub_folder_RVE = [Current_folder RVEparameters.savename '/'];
            end
            if exist(Sub_folder_RVE,'dir')==0 % Folder existence is checked, and created if necessary
                mkdir(Sub_folder_RVE);
            end
        end
        
        [All_subdomain,GROUP_SUBDOMAIN, Wholevolume_size] = Function_get_Subdomains(RVEparameters, Domain_size, voxel_size); % Location of all subdomains
        Wholevolume_size(1:2) = Wholevolume_size(1:2)*voxel_size/1000; % Size are set in micrometer
        
        % Information about subdomains
        tmp = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),...
            'VariableNames',{'Subdomain_Id' 'Group_Id' 'Number_subdomain' 'Equivalent_cubic_length' 'Section_length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
         for property_RVE=1:1:3
            for current_direction=1:1:4
                res(property_RVE).direction(current_direction).RVE(n_RVE).RVEparameters = RVEparameters;
                res(property_RVE).direction(current_direction).RVE(n_RVE).info = tmp;
            end
        end
        
        [number_subdomain,~] = size(All_subdomain); % The number of subdomain
        number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size
        
        %% ALGORITHM
        % Initialisation
        % Colunm 1 is the subdomain id
        % Colunm 2 and 3 are the sizes of the subdomain.
        Property_eachsubdomain = zeros(number_subdomain,number_phase+3,3,4); %1 Connected, 2 Isolated, 3 Unknown
        % Property calculated for each subdomain
        for subdomain_id = 1:1:number_subdomain
            % Boundary of the subdomain
            x0 = All_subdomain(subdomain_id,6); x1 = All_subdomain(subdomain_id,7);
            y0 = All_subdomain(subdomain_id,8); y1 = All_subdomain(subdomain_id,9);
            z0 = All_subdomain(subdomain_id,10); z1 = All_subdomain(subdomain_id,11);
            clear current_subdomain;
            current_subdomain = Phase_microstructure(x0:x1,y0:y1,z0:z1);
            Current_domain_size = size(current_subdomain);
            current_voxel_number = numel(current_subdomain);
            Property_eachsubdomain(subdomain_id,1,:,:)=subdomain_id;
            % Equivalent size of the subdomain
            Property_eachsubdomain(subdomain_id,2,:,:)=All_subdomain(subdomain_id,4)*voxel_size/1000; % Size are set in micrometer
            Property_eachsubdomain(subdomain_id,3,:,:)=All_subdomain(subdomain_id,5)*voxel_size/1000;            
            for current_phase=1:1:number_phase % Loop over all phases
                time_cpu_start = cputime; % CPU time start
                tic; % Stopwatch start                
                
                code_tmp = INFO.phase(current_phase).code; % Code of the phase
                binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Initialization
                binary_phase(current_subdomain == code_tmp) = 1;
                
                voxel_number = sum(sum(sum( binary_phase==1 )));
                
                [Connectivity_structure_rve] = Function_Connectivity_Algorithm(binary_phase, voxel_size, connected_id, unknown_id, isolated_id); % Call algorithm
                Property_eachsubdomain(subdomain_id,current_phase+3,1,4) = Connectivity_structure_rve.Clusters_LargestIsolatedUnknown.main_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,1,1) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(1).connected_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,1,2) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(2).connected_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,1,3) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(3).connected_cluster_phasefraction;
                
                Property_eachsubdomain(subdomain_id,current_phase+3,2,4) = Connectivity_structure_rve.Clusters_LargestIsolatedUnknown.isolated_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,2,1) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(1).isolated_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,2,2) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(2).isolated_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,2,3) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(3).isolated_cluster_phasefraction;
                
                Property_eachsubdomain(subdomain_id,current_phase+3,3,4) = Connectivity_structure_rve.Clusters_LargestIsolatedUnknown.unknown_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,3,1) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(1).unknown_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,3,2) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(2).unknown_cluster_phasefraction;
                Property_eachsubdomain(subdomain_id,current_phase+3,3,3) = Connectivity_structure_rve.Clusters_TransportFace2Face.direction(3).unknown_cluster_phasefraction;                
                
                % CPU and stopwatch time - end
                time_cpu_elapsed = cputime-time_cpu_start;
                time_stopwatch_elapsed = toc;
                Time_tmp = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
                Time_measure = [Time_measure;Time_tmp];
                
            end
        end
        
        %% STATISTICAL ANALYSIS and RVE SIZE
        str_RVE(1).name = 'ConnectedConnectivity_RVE_';
        str_RVE(2).name = 'IsolatedConnectivity_RVE_';
        str_RVE(3).name = 'UnknownConnectivity_RVE_';
        str_RVE(1).propertyname = 'Connected cluster';
        str_RVE(2).propertyname = 'Isolated cluster';  
        str_RVE(3).propertyname = 'Unknown cluster';      
        str_RVE(1).propertynameunit = '';
        str_RVE(2).propertynameunit = '';
        str_RVE(3).propertynameunit = '';
        
        res(1).Wholevolume_results = [connectivity_TransportFace2Face(:,(1-1)*3+1) connectivity_TransportFace2Face(:,(2-1)*3+1) connectivity_TransportFace2Face(:,(3-1)*3+1) connectivity_LargestIsolatedUnknown(:,1)];
        res(2).Wholevolume_results = [connectivity_TransportFace2Face(:,(1-1)*3+2) connectivity_TransportFace2Face(:,(2-1)*3+2) connectivity_TransportFace2Face(:,(3-1)*3+2) connectivity_LargestIsolatedUnknown(:,2)];
        res(3).Wholevolume_results = [connectivity_TransportFace2Face(:,(1-1)*3+3) connectivity_TransportFace2Face(:,(2-1)*3+3) connectivity_TransportFace2Face(:,(3-1)*3+3) connectivity_LargestIsolatedUnknown(:,3)];

        INFO.direction(4).filename ='Volume';
        INFO.direction(4).name ='Volume';
        
        for property_RVE=1:1:3
            for current_direction=1:1:4
                [res(property_RVE).direction(current_direction).Property_subdomains_statistics, res(property_RVE).direction(current_direction).Size_RVE] = Function_subdomains_statistical_analysis(number_group_size,number_phase,GROUP_SUBDOMAIN,Property_eachsubdomain(:,:,property_RVE,current_direction),voxel_size,RVEparameters);
                % For correlation
                for current_phase=1:1:number_phase
                    if max(res(property_RVE).direction(current_direction).Size_RVE(2,:,1))~=0
                        results_correlation(current_phase).([str_RVE(property_RVE).name INFO.direction(current_direction).filename '_'  RVEparameters.savename]) = res(property_RVE).direction(current_direction).Size_RVE(2,current_phase,1);
                    end
                    if strcmp(RVEparameters.type,'C')
                        if max(res(property_RVE).direction(current_direction).Size_RVE(2,:,2))~=0
                            results_correlation(current_phase).([str_RVE(property_RVE).name INFO.direction(current_direction).filename '_length_' RVEparameters.savename]) = res(property_RVE).direction(current_direction).Size_RVE(2,current_phase,2);
                        end
                    end
                end
                
                %% MANAGING RESULTS
                [res(property_RVE).direction(current_direction).RVE] = Function_subdomains_manage_results(Property_eachsubdomain(:,:,property_RVE,current_direction),...
                                                                                                          res(property_RVE).direction(current_direction).Property_subdomains_statistics,...
                                                                                                          RVEparameters,...
                                                                                                          res(property_RVE).direction(current_direction).RVE,...
                                                                                                          res(property_RVE).direction(current_direction).Size_RVE,...
                                                                                                          n_RVE,number_phase,INFO);
                
                %% TEXT DISPLAY AND SAVE RESULTS
                propertyname=[str_RVE(property_RVE).propertyname ' ' INFO.direction(current_direction).name];
                if property_RVE>1 || current_direction>1
                    INFO.showrveparameters=false;
                else
                    INFO.showrveparameters=true;
                end
                Function_subdomains_display_and_save(OPTIONS,INFO,res(property_RVE).direction(current_direction).RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)
                
                %% FIGURES
                parameters_figure.propertyname = propertyname;
                parameters_figure.propertynameunit = str_RVE(property_RVE).propertynameunit;
                parameters_figure.RVE = RVEparameters;
                parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.OPTIONS = OPTIONS;
                parameters_figure.INFO = INFO;
                parameters_figure.number_phase = number_phase;
                parameters_figure.Property_subdomains_statistics = res(property_RVE).direction(current_direction).Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain(:,:,property_RVE,current_direction);
                parameters_figure.Size_RVE = res(property_RVE).direction(current_direction).Size_RVE;
                parameters_figure.Wholevolume_results = res(property_RVE).Wholevolume_results(:,current_direction);
                parameters_figure.Wholevolume_size =  Wholevolume_size;
                Function_create_figures_RVE(parameters_figure) % Figures
            end
        end
        
    end
    Results_connectivity.RVE.connected = res(1).direction; % Save in main table result
    Results_connectivity.RVE.isolated = res(2).direction;
    Results_connectivity.RVE.unknown = res(3).direction;
end

%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Connectivity_calculation_time_per_phase';
[Results_connectivity] = function_time_figure(Time_measure,date_start, date_end, Results_connectivity, Current_folder, filename, 'Connectivity', OPTIONS);
 
%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_connectivity'],'Results_connectivity')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Correlation_connectivity'],'results_correlation')
    Sub_folder = 'Visualization\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Visualization_connectivity'],'results_visualization','-v7.3');    
end


end

