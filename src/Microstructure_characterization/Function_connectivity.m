function [] = Function_connectivity(Phase_microstructure,infovol,opt,p,foo1,foo2)
% Calculate connectivity (percolation)
% Function_connectivity(Phase_microstructure, infovol, opt, p) - when used with the MATBOX toolbox
% or
% Function_connectivity(Phase_microstructure, voxelsize, unit, ind, union, coupled) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       e.g.: Function_Volume_fractions(<your_3d_array>, 0.4, 'um');

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 6 % Case for function called as: Function_connectivity(Phase_microstructure, voxelsize, unit, ind, union, coupled). Standalone use.
        voxelsize = infovol; clear infovol;
        unit = opt; clear opt;
        phasetodo = p; clear p;
        unionphasetodo = foo1; clear foo1;
        coupledphasetodo = foo2; clear foo2;
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['Connectivity_' char(t)];
        if ispc
            infovol.mainfolder=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
            separator = '\';
        else
            infovol.mainfolder = pwd;
            separator = '/';
        end
        infovol.volpath = [infovol.mainfolder separator infovol.volumesubfolder separator];    

        % Set default phase information
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        infovol.phaselabel = unique(Phase_microstructure);
        for k=1:1:length(infovol.phaselabel) % Loop over all unique values
            infovol.phasename(k,1) = {['Phase ' num2str(infovol.phaselabel(k))]}; % Assign phase name based on code value
            infovol.phasecolor(k,:) = allcolor(k,:);
            p.todo(k)=1;
        end
        infovol.voxelsize = voxelsize;
        infovol.unit = unit;
        
        % Set default direction information
        infovol.directionname = {'in-plane direction 1';'in-plane direction 2';'through-plane direction'};
        
        % Set default options
        opt.save.mat = true;
        opt.save.xls = true;
        opt.save.savefig = true;
        opt.save.fig_infig = true;
        opt.save.fig_format = {'png'};
        
        % Set display options
        opt.format.fontname = 'Times New Roman';
        opt.format.axefontsize =  12;
        opt.format.legendfontsize =  10;
        opt.format.titlefontsize =  14;
        opt.format.sgtitlefontsize =  16;
        opt.format.linewidth = 2;
        opt.format.grid = 'on'; opt.format.minorgrid = 'on';
        opt.format.autoclosefig = false;
        
        % No Voxel size dependence analysis
        p.scaling = 1;
        % No Representative Volume Element analysis
        p.RVE.number_RVE = 0;
  
    else % Incorrect number of argument
        disp 'Error calling Function_connectivity. Wrong number of argument.'
        help Function_connectivity
        return
    end
    
else
    % Read parameters
    phasetodo = p.ind; % Individual phase
    unionphasetodo = p.union; % Union phase
    coupledphasetodo = p.coupled; % Couple phase
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Connectivity' separator];
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% TO BE REMOVED LATER
%parameters_scaling.scaling_factor = 2;
%parameters_scaling.label_or_greylevel = 'Label';
%parameters_scaling.background = 0;
% Scale
%Phase_microstructure = function_scaling(Phase_microstructure,parameters_scaling);
%infovol.voxelsize = infovol.voxelsize/2;

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

%%
%% ALGORITHM ON WHOLE VOLUME
%%

disp '    CONNECTIVITY';
disp '    ------------';
disp ' ';

% Voxels and clusters are marked with:
connected_id = 1; 
isolated_id = 2;
unknown_id = 3;

%% CALCULATION
% % Indiviual phase and union of phases
if strcmp(phasetodo,'None') || isempty(phasetodo)
    number_phase_todo = 0;
else
    number_phase_todo = length(phasetodo);
end
% Union of phases
if strcmp(unionphasetodo,'None') || strcmp(unionphasetodo,'n/a') || isempty(unionphasetodo)
    numberunion_todo=0;
else
    id=find(unionphasetodo==';');
    if isempty(id)
        uniontodo(1).label = str2num(unionphasetodo);
        numberunion_todo = 1;
    else
        numberunion_todo = length(id)+1;
        for kunion=1:1:numberunion_todo
            if kunion==1
                uniontodo(kunion).label = str2num(  unionphasetodo(1:id(1)-1) );
            elseif kunion==length(id)+1
                uniontodo(kunion).label = str2num(  unionphasetodo(id(kunion-1)+1:end) );
            else
                uniontodo(kunion).label = str2num(  unionphasetodo(id(kunion-1)+1:id(kunion)-1) );
            end
        end
    end
end
number_domain_todo = number_phase_todo + numberunion_todo;

if number_domain_todo>0
    time_cpu_start_volume = cputime; % CPU start
    time_clock_start_volume = tic; % Stopwatch start

    % Initialization (generic)
    Numbervoxel_phase=zeros(number_domain_todo,1);
    timedata = zeros(number_domain_todo+1,3); timedata(1,1) = voxel_number;
    timedata_domain = cell(number_domain_todo+1,1);
    phasename_todo = cell(number_domain_todo,1);

    % Initialization (algorithm-specific)
    connectivity_illdefined = zeros(number_domain_todo,2); % number voxels, volume fraction
    connectivity_LargestIsolatedUnknown = zeros(number_domain_todo,3); % largest, isolated, unknown
    connectivity_TransportFace2Face = zeros(number_domain_todo,9); % Connected, isolated, unknown for each direction
    connectivity_TransportFromFace1 = zeros(number_domain_todo,9); % Connected, isolated, unknown for each direction
    connectivity_TransportFromFace2 = zeros(number_domain_todo,9); % Connected, isolated, unknown for each direction

    for current_domain_todo=1:1:number_domain_todo % Loop over all domain to do
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        % Select label and provide names
        if current_domain_todo<=number_phase_todo
            labels = phasetodo(current_domain_todo);
            id = find(infovol.phaselabel == labels);
            phasename_todo(current_domain_todo,1) = infovol.phasename(id,1);
        else
            labels = uniontodo(current_domain_todo-number_phase_todo).label;
            tmp = [];
            for k=1:1:length(labels)
                id = find(infovol.phaselabel == labels(k));
                if k==1
                    tmp = [char(infovol.phasename(id,1))];
                else
                    tmp = [tmp ' and ' char(infovol.phasename(id,1))];
                end
            end
            phasename_todo(current_domain_todo,1) = {tmp};
        end

        % Create binary phase
        binary_phase = zeros(Domain_size);
        for k=1:1:length(labels)
            binary_phase( Phase_microstructure==labels(k) ) = 1;
        end
        Numbervoxel_phase(current_domain_todo,1) = sum(sum(sum(binary_phase)));

        [Connectivity_structure] = Function_Connectivity_Algorithm(binary_phase, voxel_size, connected_id, unknown_id, isolated_id); % Call algorithm

        % Ill connected: 1-voxel size cluster
        connectivity_illdefined(current_domain_todo,1) = Connectivity_structure.illdetailed_cluster_number;
        connectivity_illdefined(current_domain_todo,2) = Connectivity_structure.illdetailed_cluster_phasevolumefraction;
        % Largest cluster, sum of isolated clusters, sum of unknwon clusters
        connectivity_LargestIsolatedUnknown(current_domain_todo,1) = Connectivity_structure.Clusters_LargestIsolatedUnknown.main_cluster_phasefraction;
        connectivity_LargestIsolatedUnknown(current_domain_todo,2) = Connectivity_structure.Clusters_LargestIsolatedUnknown.isolated_cluster_phasefraction;
        connectivity_LargestIsolatedUnknown(current_domain_todo,3) = Connectivity_structure.Clusters_LargestIsolatedUnknown.unknown_cluster_phasefraction;
        % Oriented and directional clusters
        for direction=1:1:3
            connectivity_TransportFace2Face(current_domain_todo,(direction-1)*3+1) = Connectivity_structure.Clusters_TransportFace2Face.direction(direction).connected_cluster_phasefraction;
            connectivity_TransportFace2Face(current_domain_todo,(direction-1)*3+2) = Connectivity_structure.Clusters_TransportFace2Face.direction(direction).isolated_cluster_phasefraction;
            connectivity_TransportFace2Face(current_domain_todo,(direction-1)*3+3) = Connectivity_structure.Clusters_TransportFace2Face.direction(direction).unknown_cluster_phasefraction;
            connectivity_TransportFromFace1(current_domain_todo,(direction-1)*3+1) = Connectivity_structure.Clusters_TransportFromFace1.direction(direction).connected_cluster_phasefraction;
            connectivity_TransportFromFace1(current_domain_todo,(direction-1)*3+2) = Connectivity_structure.Clusters_TransportFromFace1.direction(direction).isolated_cluster_phasefraction;
            connectivity_TransportFromFace1(current_domain_todo,(direction-1)*3+3) = Connectivity_structure.Clusters_TransportFromFace1.direction(direction).unknown_cluster_phasefraction;
            connectivity_TransportFromFace2(current_domain_todo,(direction-1)*3+1) = Connectivity_structure.Clusters_TransportFromFace2.direction(direction).connected_cluster_phasefraction;
            connectivity_TransportFromFace2(current_domain_todo,(direction-1)*3+2) = Connectivity_structure.Clusters_TransportFromFace2.direction(direction).isolated_cluster_phasefraction;
            connectivity_TransportFromFace2(current_domain_todo,(direction-1)*3+3) = Connectivity_structure.Clusters_TransportFromFace2.direction(direction).unknown_cluster_phasefraction;
        end

        % Correlation
        results_correlation(current_domain_todo).name = phasename_todo(current_domain_todo,1);
        results_correlation(current_domain_todo).main_cluster_phasefraction = connectivity_LargestIsolatedUnknown(current_domain_todo,1);
        results_correlation(current_domain_todo).isolated_cluster_phasefraction = connectivity_LargestIsolatedUnknown(current_domain_todo,2);
        results_correlation(current_domain_todo).unknown_cluster_phasefraction = connectivity_LargestIsolatedUnknown(current_domain_todo,3);

        % Uncomment if you want to also correlate isolated and unknown clusters

        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_1 =  connectivity_TransportFace2Face(current_domain_todo,1);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_1 = connectivity_TransportFace2Face(current_domain_todo,2);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_1 = connectivity_TransportFace2Face(current_domain_todo,3);
        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_2 =  connectivity_TransportFace2Face(current_domain_todo,4);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_2 = connectivity_TransportFace2Face(current_domain_todo,5);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_2 = connectivity_TransportFace2Face(current_domain_todo,6);
        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_3 =  connectivity_TransportFace2Face(current_domain_todo,7);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_3 = connectivity_TransportFace2Face(current_domain_todo,8);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_3 = connectivity_TransportFace2Face(current_domain_todo,9);

        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_1_fromface1 =  connectivity_TransportFromFace1(current_domain_todo,1);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_1_face1 = connectivity_TransportFromFace1(current_domain_todo,2);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_1_face1 = connectivity_TransportFromFace1(current_domain_todo,3);
        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_2_fromface1 =  connectivity_TransportFromFace1(current_domain_todo,4);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_2_face1 = connectivity_TransportFromFace1(current_domain_todo,5);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_2_face1 = connectivity_TransportFromFace1(current_domain_todo,6);
        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_3_fromface1 =  connectivity_TransportFromFace1(current_domain_todo,7);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_3_face1 = connectivity_TransportFromFace1(current_domain_todo,8);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_3_face1 = connectivity_TransportFromFace1(current_domain_todo,9);

        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_1_fromface2 =  connectivity_TransportFromFace2(current_domain_todo,1);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_1_face2 = connectivity_TransportFromFace2(current_domain_todo,2);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_1_face2 = connectivity_TransportFromFace2(current_domain_todo,3);
        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_2_fromface2 =  connectivity_TransportFromFace2(current_domain_todo,4);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_2_face2 = connectivity_TransportFromFace2(current_domain_todo,5);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_2_face2 = connectivity_TransportFromFace2(current_domain_todo,6);
        results_correlation(current_domain_todo).Cluster_phasefraction_connected_along_direction_3_fromface2 =  connectivity_TransportFromFace2(current_domain_todo,7);
        %results_correlation(current_domain_todo).Cluster_phasefraction_isolated_along_direction_3_face2 = connectivity_TransportFromFace2(current_domain_todo,8);
        %results_correlation(current_domain_todo).Cluster_phasefraction_unknown_along_direction_3_face2 = connectivity_TransportFromFace2(current_domain_todo,9);

        results_visualization(current_domain_todo).name = phasename_todo(current_domain_todo,1);
        results_visualization(current_domain_todo).Clusters_sortedpersize = Connectivity_structure.Clusters_sortedpersize_3Darray;
        results_visualization(current_domain_todo).Clusters_LargestIsolatedUnknown = Connectivity_structure.Clusters_LargestIsolatedUnknown.array;
        results_visualization(current_domain_todo).Clusters_ConnectedFace2Face_along_direction_1 = Connectivity_structure.Clusters_TransportFace2Face.direction(1).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFace2Face_along_direction_2 = Connectivity_structure.Clusters_TransportFace2Face.direction(2).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFace2Face_along_direction_3 = Connectivity_structure.Clusters_TransportFace2Face.direction(3).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFromFace1_along_direction_1 = Connectivity_structure.Clusters_TransportFromFace1.direction(1).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFromFace1_along_direction_2 = Connectivity_structure.Clusters_TransportFromFace1.direction(2).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFromFace1_along_direction_3 = Connectivity_structure.Clusters_TransportFromFace1.direction(3).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFromFace2_along_direction_1 = Connectivity_structure.Clusters_TransportFromFace2.direction(1).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFromFace2_along_direction_2 = Connectivity_structure.Clusters_TransportFromFace2.direction(2).array;
        results_visualization(current_domain_todo).Clusters_ConnectedFromFace2_along_direction_3 = Connectivity_structure.Clusters_TransportFromFace2.direction(3).array;

        % % Time
        timedata_domain(current_domain_todo+1,1) = phasename_todo(current_domain_todo,1);
        timedata(current_domain_todo+1,1) = Numbervoxel_phase(current_domain_todo,1);
        timedata(current_domain_todo+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
        timedata(current_domain_todo+1,3) = toc(time_clock_start_phase); % CPU elapsed time
    end
    % CPU and stopwatch time - end
    timedata_domain(1,1) = {'Full volume'};
    timedata(1,2) = cputime-time_cpu_start_volume; % CPU elapsed time
    timedata(1,3) = toc(time_clock_start_volume); % Stopwatch elapsed time
    timedata_pervolume = timedata(1,:);
    timedata_perphase = timedata(2:end,:);
end

% % Coupled phases
[number_coupledphase_todo,~] = size(coupledphasetodo);
if number_coupledphase_todo>0
    time_cpu_start_volume = cputime; % CPU start
    time_clock_start_volume = tic; % Stopwatch start

    % Initialization (generic)
    Numbervoxel_coupledphase=zeros(number_coupledphase_todo,1);
    timedata_coupledphase = zeros(number_coupledphase_todo+1,3); timedata_coupledphase(1,1) = voxel_number;
    timedata_coupledphase_domain = cell(number_coupledphase_todo+1,1);
    coupledphasename_todo = cell(number_coupledphase_todo,1);

    % Initialization (algorithm-specific)
    support_phase_connectivity = zeros(number_coupledphase_todo,1);
    connectivity_throughsupport_catIII = zeros(number_coupledphase_todo,1);
    connectivity_direct_catIV = zeros(number_coupledphase_todo,1);
    connectivity_catIII_IV = zeros(number_coupledphase_todo,1);

    for current_coupledphase_todo = 1:1:number_coupledphase_todo
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_coupledphase = coupledphasetodo(current_coupledphase_todo,:);

        label_main = current_coupledphase(1);
        label_support = current_coupledphase(2);
        id_main = find(infovol.phaselabel == label_main);
        id_support = find(infovol.phaselabel == label_support);
        namestr = [char(infovol.phasename(id_main,1)) ' supported by ' char(infovol.phasename(id_support,1))];
        coupledphasename_todo(current_coupledphase_todo) = {namestr};

        % Support phase connectivity
        binary_phase = zeros(Domain_size);
        binary_phase( Phase_microstructure==label_support ) = 1;
        Numbervoxel_coupledphase(current_coupledphase_todo) = sum(sum(sum(binary_phase)));
        [Connectivity_structure] = Function_Connectivity_Algorithm(binary_phase, voxel_size, connected_id, unknown_id, isolated_id); % Call algorithm
        if strcmp(p.coupled_direction,'Direction 1') && strcmp(p.coupled_orientation,'Connected to first face')
            tmp = Connectivity_structure.Clusters_TransportFromFace1.direction(1).array;
            str_coupled = [char(infovol.directionname(1)) ', connected to first face'];
        elseif strcmp(p.coupled_direction,'Direction 1') && strcmp(p.coupled_orientation,'Connected to last face')
            tmp = Connectivity_structure.Clusters_TransportFromFace2.direction(1).array;
            str_coupled = [char(infovol.directionname(1)) ', connected to last face'];
        elseif strcmp(p.coupled_direction,'Direction 1') && strcmp(p.coupled_orientation,'Connected to both faces')
            tmp = Connectivity_structure.Clusters_TransportFace2Face.direction(1).array;
            str_coupled = [char(infovol.directionname(1)) ', connected to both face'];
            
        elseif strcmp(p.coupled_direction,'Direction 2') && strcmp(p.coupled_orientation,'Connected to first face')
            tmp = Connectivity_structure.Clusters_TransportFromFace1.direction(2).array;
            str_coupled = [char(infovol.directionname(2)) ', connected to first face'];            
        elseif strcmp(p.coupled_direction,'Direction 2') && strcmp(p.coupled_orientation,'Connected to last face')
            tmp = Connectivity_structure.Clusters_TransportFromFace2.direction(2).array;
            str_coupled = [char(infovol.directionname(2)) ', connected to last face'];            
        elseif strcmp(p.coupled_direction,'Direction 2') && strcmp(p.coupled_orientation,'Connected to both faces')
            tmp = Connectivity_structure.Clusters_TransportFace2Face.direction(2).array;
            str_coupled = [char(infovol.directionname(2)) ', connected to both face'];           

        elseif strcmp(p.coupled_direction,'Direction 3') && strcmp(p.coupled_orientation,'Connected to first face')
            tmp = Connectivity_structure.Clusters_TransportFromFace1.direction(3).array;
            str_coupled = [char(infovol.directionname(3)) ', connected to first face'];            
        elseif strcmp(p.coupled_direction,'Direction 3') && strcmp(p.coupled_orientation,'Connected to last face')
            tmp = Connectivity_structure.Clusters_TransportFromFace2.direction(3).array;
            str_coupled = [char(infovol.directionname(3)) ', connected to last face'];                       
        elseif strcmp(p.coupled_direction,'Direction 3') && strcmp(p.coupled_orientation,'Connected to both faces')
            tmp = Connectivity_structure.Clusters_TransportFace2Face.direction(3).array;
            str_coupled = [char(infovol.directionname(3)) ', connected to both face'];                        
        end        
        ConnectedSupport = zeros(Domain_size); ConnectedSupport(tmp==connected_id)=1;
        support_phase_connectivity(current_coupledphase_todo) = 100 * sum(sum(sum(ConnectedSupport))) / sum(sum(sum(binary_phase)));

        % Main phase connectivity
        binary_phase = zeros(Domain_size);
        binary_phase( Phase_microstructure==label_main ) = 1;
        Numbervoxel_coupledphase(current_coupledphase_todo) = Numbervoxel_coupledphase(current_coupledphase_todo) + sum(sum(sum(binary_phase)));
        MainPhaseClusters = bwlabeln(binary_phase,6);
        Id_clusters = unique(MainPhaseClusters);
        MainPhaseClusters_and_connectedsupport = MainPhaseClusters;
        MainPhaseClusters_and_connectedsupport(ConnectedSupport==1) = max(Id_clusters)+1;

        % Connectivity matrix
        background=0;
        interface_complemetaryvolume=1;
        interface_anotherlabel=2;
        [Connectivity_matrix, ~] = Function_connectivitymatrix(MainPhaseClusters_and_connectedsupport, background, interface_complemetaryvolume, interface_anotherlabel);

        % Clusters of main phase in contact with edge through connected cluster of support phase
        row = Connectivity_matrix(end,:);
        idx = find(row>0);
        idcluster_contact_supportphase = Connectivity_matrix(1,idx);
        idcluster_contact_supportphase(1)=[]; idcluster_contact_supportphase(end)=[]; % Renove support
        s=0;
        MainsPhaseConnected_throughSupport_catIII = zeros(Domain_size);
        for kcluster=1:1:length(idcluster_contact_supportphase)
            MainsPhaseConnected_throughSupport_catIII( MainPhaseClusters==idcluster_contact_supportphase(kcluster) ) = 1;
        end
        connectivity_throughsupport_catIII(current_coupledphase_todo) = 100 * sum(sum(sum(MainsPhaseConnected_throughSupport_catIII))) / sum(sum(sum(binary_phase)));

        % Clusters of main phase in contact with edge directly
        MainsPhaseConnected_directly_catIV = zeros(Domain_size);
        for kcluster=2:1:length(Id_clusters) % Skip background
            idp = find(MainPhaseClusters==Id_clusters(kcluster));
            [IX,IY,IZ]=ind2sub(Domain_size,idp);
            if strcmp(p.coupled_direction,'Direction 1') && strcmp(p.coupled_orientation,'Connected to first face') && min(IX)==1
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 1') && strcmp(p.coupled_orientation,'Connected to last face') && max(IX)==Domain_size(1)
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 1') && strcmp(p.coupled_orientation,'Connected to both faces') && min(IX)==1 && max(IX)==Domain_size(1)
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 2') && strcmp(p.coupled_orientation,'Connected to first face') && min(IY)==1
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 2') && strcmp(p.coupled_orientation,'Connected to last face') && max(IY)==Domain_size(2) 
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 2') && strcmp(p.coupled_orientation,'Connected to both faces') && min(IY)==1 && max(IY)==Domain_size(2) 
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 3') && strcmp(p.coupled_orientation,'Connected to first face') && min(IZ)==1
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 3') && strcmp(p.coupled_orientation,'Connected to last face') && max(IZ)==Domain_size(3)
                MainsPhaseConnected_directly_catIV(idp)=1;
            elseif strcmp(p.coupled_direction,'Direction 3') && strcmp(p.coupled_orientation,'Connected to both faces') && min(IZ)==1 && max(IZ)==Domain_size(3)
                MainsPhaseConnected_directly_catIV(idp)=1;
            end 

        end
        connectivity_direct_catIV(current_coupledphase_todo) = 100 * sum(sum(sum(MainsPhaseConnected_directly_catIV))) / sum(sum(sum(binary_phase)));

        % Clusters of main phase in contact with edge both ways
        MainsPhaseConnected_catIII_IV = MainsPhaseConnected_throughSupport_catIII + MainsPhaseConnected_directly_catIV;
        MainsPhaseConnected_catIII_IV(MainsPhaseConnected_catIII_IV~=0)=1;
        connectivity_catIII_IV(current_coupledphase_todo) = 100 * sum(sum(sum(MainsPhaseConnected_catIII_IV))) / sum(sum(sum(binary_phase)));

        % Array: 0 background, 1 support phase (connected), 2 support phase (not connected), 3 main phase (connected through support), 4 main phase (connected directly), 5 main phase (not connected)
        Supported_connectivity = zeros(Domain_size);
        Supported_connectivity( Phase_microstructure==label_support ) = 2;
        Supported_connectivity( ConnectedSupport==1 ) = 1;
        Supported_connectivity( binary_phase==1 ) = 5;      
        Supported_connectivity( MainsPhaseConnected_throughSupport_catIII==1 ) = 3;          
        Supported_connectivity( MainsPhaseConnected_directly_catIV==1 ) = 4;          

        % Correlation
        results_correlation(number_domain_todo+current_coupledphase_todo).name = {[namestr ' ' str_coupled]};
        results_correlation(number_domain_todo+current_coupledphase_todo).support_phase_connectivity =  support_phase_connectivity(current_coupledphase_todo);
        results_correlation(number_domain_todo+current_coupledphase_todo).connectivity_throughsupport_catIII =  connectivity_throughsupport_catIII(current_coupledphase_todo);
        results_correlation(number_domain_todo+current_coupledphase_todo).connectivity_direct_catIV =  connectivity_direct_catIV(current_coupledphase_todo);
        results_correlation(number_domain_todo+current_coupledphase_todo).connectivity_catIII_IV =  connectivity_catIII_IV(current_coupledphase_todo);

        % Visualization
        results_visualization(number_domain_todo+current_coupledphase_todo).name = {[namestr ' ' str_coupled]};
        results_visualization(number_domain_todo+current_coupledphase_todo).Supported_connectivity = Supported_connectivity;

        % % Time
        timedata_coupledphase_domain(current_coupledphase_todo+1,1) = {[namestr ' ' str_coupled]};
        timedata_coupledphase(current_coupledphase_todo+1,1) = Numbervoxel_coupledphase(current_coupledphase_todo);
        timedata_coupledphase(current_coupledphase_todo+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
        timedata_coupledphase(current_coupledphase_todo+1,3) = toc(time_clock_start_phase); % CPU elapsed time        

    end
    % CPU and stopwatch time - end
    timedata_coupledphase_domain(1,1) = {'Full volume'};
    timedata_coupledphase(1,2) = cputime-time_cpu_start_volume; % CPU elapsed time
    timedata_coupledphase(1,3) = toc(time_clock_start_volume); % Stopwatch elapsed time
    timedata_pervolume = timedata_coupledphase(1,:);
    timedata_perphase = timedata_coupledphase(2:end,:);
end

%% TABLES

if number_domain_todo>0
    % Time
    Table_time = table(timedata_domain(:,1), timedata(:,1),timedata(:,2),timedata(:,3),...
        'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
    Results_connectivity.Table_time = Table_time; % Save in main table result

    % Result calculated on whole volume
    Table_connectivity = table(phasename_todo,connectivity_LargestIsolatedUnknown(:,1),connectivity_LargestIsolatedUnknown(:,2),connectivity_LargestIsolatedUnknown(:,3),...
        'VariableNames',{'Phase(s)' 'Main cluster' 'Isolated clusters' 'Unknown clusters'});%
    Table_face2face = table(phasename_todo,connectivity_TransportFace2Face(:,1),connectivity_TransportFace2Face(:,2),connectivity_TransportFace2Face(:,3),connectivity_TransportFace2Face(:,4),connectivity_TransportFace2Face(:,5),connectivity_TransportFace2Face(:,6),connectivity_TransportFace2Face(:,7),connectivity_TransportFace2Face(:,8),connectivity_TransportFace2Face(:,9),...
        'VariableNames',{'Phase(s)' 'Dir1 connected' 'Dir1 isolated' 'Dir1 unknown' 'Dir2 connected' 'Dir2 isolated' 'Dir2 unknown' 'Dir3 connected' 'Dir3 isolated' 'Dir3 unknown'});
    Table_fromface1 = table(phasename_todo,connectivity_TransportFromFace1(:,1),connectivity_TransportFromFace1(:,2),connectivity_TransportFromFace1(:,3),connectivity_TransportFromFace1(:,4),connectivity_TransportFromFace1(:,5),connectivity_TransportFromFace1(:,6),connectivity_TransportFromFace1(:,7),connectivity_TransportFromFace1(:,8),connectivity_TransportFromFace1(:,9),...
        'VariableNames',{'Phase(s)' 'Dir1 connected' 'Dir1 isolated' 'Dir1 unknown' 'Dir2 connected' 'Dir2 isolated' 'Dir2 unknown' 'Dir3 connected' 'Dir3 isolated' 'Dir3 unknown'});
    Table_fromface2 = table(phasename_todo,connectivity_TransportFromFace2(:,1),connectivity_TransportFromFace2(:,2),connectivity_TransportFromFace2(:,3),connectivity_TransportFromFace2(:,4),connectivity_TransportFromFace2(:,5),connectivity_TransportFromFace2(:,6),connectivity_TransportFromFace2(:,7),connectivity_TransportFromFace2(:,8),connectivity_TransportFromFace2(:,9),...
        'VariableNames',{'Phase(s)' 'Dir1 connected' 'Dir1 isolated' 'Dir1 unknown' 'Dir2 connected' 'Dir2 isolated' 'Dir2 unknown' 'Dir3 connected' 'Dir3 isolated' 'Dir3 unknown'});
    Table_illdefinedclusters = table(phasename_todo,connectivity_illdefined(:,1),connectivity_illdefined(:,2),...
        'VariableNames',{'Phase(s)' 'Number one-voxel size cluster' 'Phase percentage'});
    Results_connectivity.Table_connectivity = Table_connectivity; % Save in main table result
    Results_connectivity.Table_face2face = Table_face2face;
    Results_connectivity.Table_fromface1 = Table_fromface1;
    Results_connectivity.Table_fromface2 = Table_fromface2;
    Results_connectivity.Table_illdefinedclusters = Table_illdefinedclusters;
end

if number_coupledphase_todo>0
    % Time
    Table_coupledphase_time = table(timedata_coupledphase_domain(:,1), timedata_coupledphase(:,1),timedata_coupledphase(:,2),timedata_coupledphase(:,3),...
        'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
    Results_connectivity.Table_coupledphase_time = Table_coupledphase_time; % Save in main table result

    coupled_direction = cell(number_coupledphase_todo,1);
    coupled_direction(:,1)={str_coupled};
    Table_coupledconnectivity = table(coupled_direction, coupledphasename_todo,support_phase_connectivity,connectivity_throughsupport_catIII,connectivity_direct_catIV,connectivity_catIII_IV,...
        'VariableNames',{'Direction' 'Phases' 'Support phase connectivity' 'Main phase connected through support phase' 'Main phase directly connected' 'Main phase connected both ways'});
    Results_connectivity.Table_coupledconnectivity = Table_coupledconnectivity;
end

%% SAVE TABLES
if opt.save.xls
    filename = 'Connectivity'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    sheetnumber = 1;
    if number_domain_todo>0
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
        sheetnumber = 6;
    end
    if number_coupledphase_todo>0
        DATA_writetable.sheet(sheetnumber).name='Coupled connectivity';
        DATA_writetable.sheet(sheetnumber).table=Table_coupledconnectivity;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
if number_domain_todo>0
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
if number_coupledphase_todo>0
    disp 'Coupled Connectivity';
    disp(Table_coupledconnectivity)
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_coupledphase_time)    
end

% keyboard
% 
% %%
% %% ADDITIONAL RESULTS ON THE WHOLE VOLUME
% %%
% 
% %% ALONG DIRECTIONS
% if number_domain_todo>0
%     % Initialize
%     for current_direction=1:1:number_dimension % Loop over all directions
%         tmp = Domain_size; tmp(current_direction)=1; Results.direction(current_direction).numbervoxelslice = prod(tmp); clear tmp; % Number of voxel within each slice normal to the direction
%         Results.direction(current_direction).connectivity = zeros(Domain_size(current_direction),number_domain_todo+1,3); % Initialization
%         Results.direction(current_direction).connectivity_F2F1 = zeros(Domain_size(current_direction),number_domain_todo+1,3); % Initialization
%         Results.direction(current_direction).connectivity_F2F2 = zeros(Domain_size(current_direction),number_domain_todo+1,3); % Initialization
%         Results.direction(current_direction).connectivity_F2F3 = zeros(Domain_size(current_direction),number_domain_todo+1,3); % Initialization
%         for k=1:1:3
%             Results.direction(current_direction).connectivity(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size)'; % Axis (position along direction)
%             Results.direction(current_direction).connectivity_F2F1(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size)'; % Axis (position along direction)
%             Results.direction(current_direction).connectivity_F2F2(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size)'; % Axis (position along direction)
%             Results.direction(current_direction).connectivity_F2F3(:,1,k) = ((1:1:Domain_size(current_direction))*voxel_size)'; % Axis (position along direction)
%         end
%     end
%     % Calculate results
%     for current_domain_todo = 1:1:number_domain_todo % Loop over all phases
%         Vtmp = results_visualization(current_domain_todo).Clusters_LargestIsolatedUnknown;
%         F1tmp = results_visualization(current_domain_todo).Clusters_ConnectedFace2Face_along_direction_1;
%         F2tmp = results_visualization(current_domain_todo).Clusters_ConnectedFace2Face_along_direction_2;
%         F3tmp = results_visualization(current_domain_todo).Clusters_ConnectedFace2Face_along_direction_3;
%         for c=1:1:3
%             % Direction 1
%             for current_position=1:1:Domain_size(1) % Loop over postion
%                 n= sum(sum(Vtmp(current_position,:,:)==c)); phasevoxel = sum(sum(Vtmp(current_position,:,:)>0));
%                 Results.direction(1).connectivity(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F1tmp(current_position,:,:)==c));
%                 Results.direction(1).connectivity_F2F1(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F2tmp(current_position,:,:)==c));
%                 Results.direction(1).connectivity_F2F2(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F3tmp(current_position,:,:)==c));
%                 Results.direction(1).connectivity_F2F3(current_position,current_domain_todo+1,c) = n / phasevoxel;
%             end
%             % Direction 2
%             for current_position=1:1:Domain_size(2) % Loop over postion
%                 n= sum(sum(Vtmp(:,current_position,:)==c));  phasevoxel = sum(sum(Vtmp(:,current_position,:)>0));
%                 Results.direction(2).connectivity(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F1tmp(:,current_position,:)==c));
%                 Results.direction(2).connectivity_F2F1(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F2tmp(:,current_position,:)==c));
%                 Results.direction(2).connectivity_F2F2(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F3tmp(:,current_position,:)==c));
%                 Results.direction(2).connectivity_F2F3(current_position,current_domain_todo+1,c) = n / phasevoxel;
%             end
%             % Direction 3
%             for current_position=1:1:Domain_size(3) % Loop over postion
%                 n= sum(sum(Vtmp(:,:,current_position)==c));  phasevoxel = sum(sum(Vtmp(:,:,current_position)>0));
%                 Results.direction(3).connectivity(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F1tmp(:,:,current_position)==c));
%                 Results.direction(3).connectivity_F2F1(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F2tmp(:,:,current_position)==c));
%                 Results.direction(3).connectivity_F2F2(current_position,current_domain_todo+1,c) = n / phasevoxel;
%                 n= sum(sum(F3tmp(:,:,current_position)==c));
%                 Results.direction(3).connectivity_F2F3(current_position,current_domain_todo+1,c) = n / phasevoxel;
%             end
%         end
%     end
%     clear Vtmp F1tmp F2tmp F3tmp
% end
% 
% 
% % Initialize
% for current_direction=1:1:number_dimension % Loop over all directions
%     tmp = Domain_size; tmp(current_direction)=1; Results.direction(current_direction).numbervoxelslice = prod(tmp); clear tmp; % Number of voxel within each slice normal to the direction
%     Results.direction(current_direction).volumefraction = zeros(Domain_size(current_direction),number_phase_todo+1); % Initialization
%     Results.direction(current_direction).volumefraction(:,1) = (1:1:Domain_size(current_direction))*voxel_size; % Position along direction
% end
% % Calculate results
% current_phase_todo = 0;
% for current_phase = 1:1:number_phase % Loop over all phases
%     if p.todo(current_phase)
%         current_phase_todo=current_phase_todo+1;
% 
%         label = infovol.phaselabel(current_phase);
%         % Direction 1
%         for current_position=1:1:Domain_size(1) % Loop over postion
%             Numbervoxel_position= sum(sum(Phase_microstructure(current_position,:,:)==label));
%             Results.direction(1).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(1).numbervoxelslice;
%         end
%         % Direction 2
%         for current_position=1:1:Domain_size(2) % Loop over postion
%             Numbervoxel_position= sum(sum(Phase_microstructure(:,current_position,:)==label));
%             Results.direction(2).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(2).numbervoxelslice;
%         end
%         if number_dimension==3 % 3D case
%             % Direction 3
%             for current_position=1:1:Domain_size(3) % Loop over postion
%                 Numbervoxel_position= sum(sum(Phase_microstructure(:,:,current_position)==label));
%                 Results.direction(3).volumefraction(current_position,current_phase_todo+1) = Numbervoxel_position / Results.direction(3).numbervoxelslice;
%             end
%         end
%     end
% end
% 
% %% TABLES
% % Prepare header name
% clear Variable_name_table;
% Variable_name_table={['Position ' voxel_unit]};
% for current_phase_todo=1:1:number_phase_todo
%     Variable_name_table(current_phase_todo+1)=phasename_todo(current_phase_todo);
% end
% % Create table
% for current_direction=1:1:number_dimension % Loop over all directions
%     EvolutionVolumefraction.direction(current_direction).table = array2table(Results.direction(current_direction).volumefraction,'VariableNames',Variable_name_table);
% end
% Results_Volumefraction.EvolutionVolumefraction = EvolutionVolumefraction; % Save in main table result
% 
% %% SAVE TABLES
% if opt.save.xls
%     filename = 'Volume_fraction_along_directions'; % Filename without extension
%     % Prepare the data
%     clear DATA_writetable
%     for current_direction=1:1:number_dimension % Loop over all directions
%         DATA_writetable.sheet(current_direction).name = char(infovol.directionname(current_direction));
%         DATA_writetable.sheet(current_direction).table = EvolutionVolumefraction.direction(current_direction).table;
%     end
%     % Save function
%     Function_Writetable(Current_folder,filename,DATA_writetable)
% end
% 
% %% FIGURES
% strunit = voxel_unit;
% if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
%     axisunit = '(\mum)';
% else
%     axisunit = ['(' strunit ')'];
% end
% 
% scrsz = get(0,'ScreenSize'); % Screen resolution
% Fig = figure; % Create figure
% Fig.Name= 'Volume fractions'; % Figure name
% Fig.Color='white'; % Background colour
% set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
% for current_direction=1:1:number_dimension % Iterate over axe
%     sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
%     hold(sub_axes,'on'); % Active subplot
%     h_title=title(['Volume fractions along ',char(infovol.directionname(current_direction))]);
%     % Plot graphs
%     x_ = Results.direction(current_direction).volumefraction(:,1);
%     current_phase_todo = 0;
%     for current_phase=1:1:number_phase % Loop over phases
%         if p.todo(current_phase)
%             current_phase_todo=current_phase_todo+1;
%             y_ = Results.direction(current_direction).volumefraction(:,current_phase_todo+1);
%             plot(x_,y_,'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
%         end
%     end
%     % Axis label
%     t_ = xlabel(' ');
%     t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
%     t_2 = axisunit;
%     t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
%     ylabel('Volume fractions');
%     % Legend
%     for current_phase_todo=1:1:number_phase_todo
%         str_legend(current_phase_todo).name = [char(phasename_todo(current_phase_todo)) ', <\epsilon>=' num2str(Volumefraction_phase(current_phase_todo),'%1.3f')];
%     end
%     h_legend = legend(sub_axes,str_legend.name,'Location','best');
%     % - Grid
%     grid(sub_axes,opt.format.grid); % Display grid
%     set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
%     set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
%     h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
%     h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
%     hold(sub_axes,'off'); % Relase figure    
% end
% sgtitle(Fig,'Volume fractions along directions','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
% if opt.save.savefig % Save figure
%     filename= 'Volume_fractions_along_direction';
%     function_savefig(Fig, Current_folder, filename, opt.save); % Call function
% end
% if opt.format.autoclosefig
%     close(Fig); % Do not keep open figures
% end
% 
% 
% %%
% %% IMAGE RESOLUTION SENSITIVITY ANALYSIS
% %%
% 
% if length(p.scaling)>=2 % Check if voxel size analysis is asked
%     size_choice = p.scaling;
%     size_choice = sort(size_choice);
%     size_choice(size_choice==1)=[];
%     number_resize=length(size_choice); % Number of different voxel size that will be analyzed
%     
%     %% CALCULATION
%     % Initialization
%     property_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1);
%     property_voxelsizedependence(1,1)=voxel_size;
%     property_voxelsizedependence(1,2:end)=Volumefraction_phase(:,1)';
% 
%     if p.fractal.todo && p.fractal.voxelsize
%         fractaldimension_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1);
%         fractaldimension_voxelsizedependence(1,1) = voxel_size;
%         fractaldimension_voxelsizedependence(1,2:end) = fractal_dimension(2,:);
%     end
% 
%     % Loop on each voxel size
%     for current_iteration=1:1:number_resize
% 
%         % % Microstructure scaling
%         % New voxel size
%         property_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
%         % Set parameters
%         parameters_scaling.scaling_factor = size_choice(current_iteration);
%         parameters_scaling.label_or_greylevel = 'Label';
%         parameters_scaling.background = min(infovol.phaselabel);
%         % Scale
%         Phase_microstructure_resized = function_scaling(Phase_microstructure,parameters_scaling);
% 
%         % CPU and stopwatch time - start
%         time_cpu_start_volume = cputime; % CPU start
%         time_clock_start_volume = tic; % Stopwatch start
%         % Number of voxel of the current resized microstructure
%         voxel_number_tmp=numel(Phase_microstructure_resized);
% 
%         current_phase_todo = 0;
%         for current_phase=1:1:number_phase % Loop over all phases
%             if p.todo(current_phase)
%                 time_cpu_start_phase = cputime; % CPU start
%                 time_clock_start_phase = tic; % Stopwatch start
% 
%                 current_phase_todo=current_phase_todo+1;
% 
%                 % % Algorithm: SPECIFIC FOR EACH FILE
%                 code_tmp = infovol.phaselabel(current_phase);
%                 Numbervoxel_phase_tmp= sum(sum(sum(Phase_microstructure_resized==code_tmp )));
%                 property_voxelsizedependence(current_iteration+1,current_phase_todo+1)=Numbervoxel_phase_tmp/voxel_number_tmp;
% 
%                 % % Time
%                 timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];
% 
%             end
%         end
%         % CPU and stopwatch time - end
%         timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
% 
%         % % Fractal
%         if p.fractal.todo && p.fractal.voxelsize
%             for current_phase_todo=1:1:number_phase_todo % Loop over phase and call box counting algorithm
%                 binary_array = zeros(size(Phase_microstructure_resized));
%                 binary_array(Phase_microstructure_resized==phaselabel(current_phase_todo,1))=1;
%                 [~,~,fractal_dimension_tmp, ~] = Function_fractaldimension_boxcounting(binary_array,p.fractal);
%                 fractaldimension_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
%                 fractaldimension_voxelsizedependence(current_iteration+1,current_phase_todo+1) = fractal_dimension_tmp(2,1);
%             end      
%         end
% 
%     end
%     clear Phase_microstructure_resized;
% 
%     % Sort per voxel size
%     property_voxelsizedependence = sortrows(property_voxelsizedependence,1);   
%     
%     %% EXTRAPOLATION TO 0 nm
%     fprintf('> Volume fractions dependence with the voxel size\n');
%     tmp = zeros(number_resize+2,number_phase_todo+1); % + 0 nm and + initial voxel size
%     x=property_voxelsizedependence(:,1);
%     if strcmp(p.scaling_extrapolation,'Linear')
%         interpolation_voxelsize_order=1;
%     elseif strcmp(p.scaling_extrapolation,'Quadratic')
%         interpolation_voxelsize_order=2;
%     elseif strcmp(p.scaling_extrapolation,'Cubic')
%         interpolation_voxelsize_order=3;
%     end
%     max_order = length(p.scaling)-1;
%     interpolation_voxelsize_order = min(interpolation_voxelsize_order,max_order);
%     fprintf('  Extrapolation to zero voxel size: polynomial of order %i\n\n',interpolation_voxelsize_order);
%     for current_phase_todo=1:1:number_phase_todo
%         y=property_voxelsizedependence(:,current_phase_todo+1);
%         pi = polyfit(x,y,interpolation_voxelsize_order);
%         vq = polyval(pi,0);
%         tmp(1,current_phase_todo+1)=vq;
%         interpolation_voxelsize(current_phase_todo).pi=pi;
%     end
%     tmp(2:end,:) = property_voxelsizedependence;
%     property_voxelsizedependence = tmp; clear tmp;
%     
%     %% MANAGING RESULTS
%     % Results are saved in a table
%     Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
%     for current_phase_todo=1:1:number_phase_todo
%         Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
%     end
%     % Table
%     Table_Volumefraction_voxelsizedependence = array2table(property_voxelsizedependence,...
%         'VariableNames',Variable_name_table);
%     
%     %% DISPLAY TEXT RESULTS
%     disp(Table_Volumefraction_voxelsizedependence)
%     
%     %% SAVE RESULTS
%     Results_Volumefraction.voxelsizedependence = Table_Volumefraction_voxelsizedependence; % Save in main table result
%     if opt.save.xls
%         filename = 'Volume_fraction_voxel_size_dependence'; % Filename without extension
%         % Prepare the data
%         clear DATA_writetable
%         DATA_writetable.sheet(1).name='Volume_fractions';
%         DATA_writetable.sheet(1).table=Table_Volumefraction_voxelsizedependence;
%         % Save function
%         Function_Writetable(Current_folder,filename,DATA_writetable)
%     end
%     
%     %% FIGURES
%     parameters_figure.plotlog = false; 
%     parameters_figure.propertyname = 'Volume fractions';
%     parameters_figure.method = [];
%     parameters_figure.property_voxelsizedependence = property_voxelsizedependence;
%     parameters_figure.number_phase = number_phase_todo;
%     parameters_figure.str_ylabel = 'Volume fractions';
%     parameters_figure.propertynameunit = [];
%     parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
%     parameters_figure.Current_folder = Current_folder;
%     parameters_figure.filename = 'Volume_fractions_voxel_size_dependence';
%     parameters_figure.infovol = infovol;
%     parameters_figure.opt = opt;
%     parameters_figure.todo = p.todo;
%     Function_create_figure_voxelsizedependence(parameters_figure) % Figures    
% 
%     %% FRACTAL DIMENSION
%     if p.fractal.todo && p.fractal.voxelsize
%         fprintf('> Phase domain fractal dimension dependence with the voxel size\n');
%         fprintf('  Extrapolation to zero voxel size: polynomial of order 1\n\n');
% 
%         % Sort by voxel size
%         fractaldimension_voxelsizedependence = sortrows(fractaldimension_voxelsizedependence,1);
% 
%         % Extrapolation to 0 voxel size
%         tmp = zeros(number_resize+2,number_phase_todo+1); % + 0 nm and + initial voxel size
%         for current_phase_todo=1:1:number_phase_todo
%             y=fractaldimension_voxelsizedependence(:,current_phase_todo+1);
%             pi = polyfit(x,y,1);
%             vq = polyval(pi,0);
%             tmp(1,current_phase_todo+1)=vq;
%             interpolation_voxelsize(current_phase_todo).pi=pi;           
%         end
%         tmp(2:end,:) = fractaldimension_voxelsizedependence;
%         fractaldimension_voxelsizedependence = tmp; clear tmp;
% 
%         % Managing result
%         Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
%         for current_phase_todo=1:1:number_phase_todo
%             Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
%         end
%         % Table
%         Table_Fractaldimension_voxelsizedependence = array2table(fractaldimension_voxelsizedependence,...
%             'VariableNames',Variable_name_table);
%         fprintf('    fitted from s=1 to s=2\n');
%         disp(Table_Fractaldimension_voxelsizedependence)
% 
%         % Save result
%         Results_Volumefraction.Fractaldimension_voxelsizedependence = Table_Fractaldimension_voxelsizedependence; % Save in main table result
%         if opt.save.xls
%             filename = 'PhaseDomain_Fractaldimension_voxel_size_dependence'; % Filename without extension
%             % Prepare the data
%             clear DATA_writetable
%             DATA_writetable.sheet(1).name='Fit from 1 to 2';
%             DATA_writetable.sheet(1).table=Table_Fractaldimension_voxelsizedependence;
%             % Save function
%             Function_Writetable(Current_folder,filename,DATA_writetable)
%         end
% 
%         % Correlation
%         for current_phase_todo=1:1:number_phase_todo
%             results_correlation(current_phase_todo).PhaseDomain_fractaldimension_extrapolated = fractaldimension_voxelsizedependence(1,current_phase_todo+1) ;
%             results_correlation(current_phase_todo).PhaseDomain_fractalpropensity_extrapolated = abs(p.fractal.topology_dimension - fractaldimension_voxelsizedependence(1,current_phase_todo+1));
%         end
% 
%         % Figure
%         parameters_figure.plotlog = false; 
%         parameters_figure.figname = 'Phase domain fractal dimension';
%         parameters_figure.propertyname = 'Phase domain fractal dimension';
%         parameters_figure.method = 'Box counting';
%         parameters_figure.property_voxelsizedependence = fractaldimension_voxelsizedependence;
%         parameters_figure.number_phase = number_phase_todo;
%         parameters_figure.str_ylabel = 'Fractal dimension';
%         parameters_figure.propertynameunit = [];
%         parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
%         parameters_figure.Current_folder = Current_folder;
%         parameters_figure.filename = 'Phase_Fractaldimension_voxel_size_dependence';
%         parameters_figure.infovol = infovol;
%         parameters_figure.opt = opt;
%         parameters_figure.todo = p.todo;
%         Function_create_figure_voxelsizedependence(parameters_figure) % Figures
%     end
% 
% end
% 
% 
% %%
% %% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
% %%
% 
% if p.RVE.number_RVE>0
%     for k_RVE = 1:1:p.RVE.number_RVE % Loop over all RVE asked
%         RVEparameters = p.RVE.RVE(k_RVE);
%         RVE(k_RVE).RVEparameters = RVEparameters; % For result structure
%         if opt.save.xls || opt.save.savefig
%             Sub_folder_RVE = [Current_folder RVEparameters.savename separator];
%             while exist(Sub_folder_RVE,'dir')
%                 RVEparameters.savename = [RVEparameters.savename '_bis'];
%                 Sub_folder_RVE = [Current_folder RVEparameters.savename separator];
%             end
%             mkdir(Sub_folder_RVE);
%         end
%         if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C')
%             thresholds = RVEparameters.threshold_std;
%         else
%             thresholds = RVEparameters.threshold_reldiff;
%         end
%         n_threshold = length(thresholds);
% 
%         % Nested analysis ?
%         if RVEparameters.donested % Yes
%             [n_nestedRVE, xcrop] = Function_NestedRVE(RVEparameters,Domain_size);
%             Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_phase_todo,3,3); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
%         else % No
%             n_nestedRVE = 0;
%         end
% 
%         for k_nestedRVE = 0:1:n_nestedRVE
% 
%             %% NESTED ANALYSIS
%             if k_nestedRVE == 0
%                 Domain_size_nested = Domain_size;
%             else
%                 x0 = xcrop(k_nestedRVE,1); x1 = xcrop(k_nestedRVE,2);
%                 y0 = xcrop(k_nestedRVE,3); y1 = xcrop(k_nestedRVE,4);
%                 z0 = xcrop(k_nestedRVE,5); z1 = xcrop(k_nestedRVE,6);
%                 Phase_microstructure_nested = Phase_microstructure(x0:x1,y0:y1,z0:z1);
%                 Domain_size_nested = size(Phase_microstructure_nested);
%             end
%             
%             %% SUBDOMAINS
%             [All_subdomain,GROUP_SUBDOMAIN, Wholevolume_size] = Function_get_Subdomains(RVEparameters, Domain_size_nested, Domain_size); % Location of all subdomains
%             Wholevolume_size(1:2) = Wholevolume_size(1:2)*voxel_size;
%             if RVEparameters.donested
%                 Result_nestedRVE(k_nestedRVE+1,1,:,:,1) = Wholevolume_size(1);
%                 if k_nestedRVE ==1
%                     fprintf('       Nested analysis\n');
%                 end
%                 if k_nestedRVE > 0
%                     fprintf('          Cropping iteration: #%i/%i (cubic root volume=%1.1f %s)\n',k_nestedRVE,n_nestedRVE,Wholevolume_size(1),infovol.unit);
%                 end                
%                 if strcmp(RVEparameters.type,'C')
%                     Result_nestedRVE(k_nestedRVE+1,1,:,:,2) = Wholevolume_size(2);
%                 elseif strcmp(RVEparameters.type,'F')
%                     Result_nestedRVE(k_nestedRVE+1,1,:,:,3) = Wholevolume_size(2);
%                 end
%             end
% 
%             % Information about subdomains
%             RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
%                 'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Equivalent cubic length' 'Equivalent section length' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
%             [number_subdomain,~] = size(All_subdomain); % The number of subdomain
%             number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size
% 
%             %% ALGORITHM
%             % Colunm 1 is the subdomain id
%             % Colunm 2 and 3 are the sizes of the subdomain.
%             Property_eachsubdomain = zeros(number_subdomain,number_phase_todo+4);
%             % Property calculated for each subdomain
%             for subdomain_id = 1:1:number_subdomain
%                 % Boundary of the subdomain
%                 x0 = All_subdomain(subdomain_id,7); x1 = All_subdomain(subdomain_id,8);
%                 y0 = All_subdomain(subdomain_id,9); y1 = All_subdomain(subdomain_id,10);
%                 z0 = All_subdomain(subdomain_id,11); z1 = All_subdomain(subdomain_id,12);
%                 clear current_subdomain;
%                 % Crop volume
%                 if k_nestedRVE == 0
%                     current_subdomain = Phase_microstructure(x0:x1,y0:y1,z0:z1);
%                 else
%                     current_subdomain = Phase_microstructure_nested(x0:x1,y0:y1,z0:z1);
%                 end                    
% 
%                 % CPU and stopwatch time - start
%                 time_cpu_start_volume = cputime; % CPU start
%                 time_clock_start_volume = tic; % Stopwatch start
%                 % Number of voxel of the current resized microstructure
%                 voxel_number_tmp=numel(current_subdomain);
% 
%                 Property_eachsubdomain(subdomain_id,1)=subdomain_id;
%                 % Equivalent size of the subdomain
%                 Property_eachsubdomain(subdomain_id,2)=All_subdomain(subdomain_id,4)*voxel_size; % Cubic root length
%                 Property_eachsubdomain(subdomain_id,3)=All_subdomain(subdomain_id,5)*voxel_size; % Square root length
%                 Property_eachsubdomain(subdomain_id,4)=All_subdomain(subdomain_id,6)*voxel_size; % Length
% 
%                 current_phase_todo = 0;
%                 for current_phase=1:1:number_phase % Loop over all phases
%                     if p.todo(current_phase)
%                         time_cpu_start_phase = cputime; % CPU start
%                         time_clock_start_phase = tic; % Stopwatch start
%                         current_phase_todo=current_phase_todo+1;
% 
%                         % % Algorithm: SPECIFIC FOR EACH FILE
%                         code_tmp = infovol.phaselabel(current_phase);
%                         Numbervoxel_phase_tmp= sum(sum(sum(current_subdomain==code_tmp )));
%                         Property_eachsubdomain(subdomain_id,current_phase_todo+4)=Numbervoxel_phase_tmp/voxel_number_tmp;
% 
%                         % % Time
%                         timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];
%                     end
%                 end
%                 % CPU and stopwatch time - end
%                 timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
%             end
% 
%             %% STATISTICAL ANALYSIS and RVE SIZE
%             [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_phase_todo,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters);
%             
%             %% SAVE FOR CORRELATION
%             if k_nestedRVE == 0
%                 for k_threshold=1:1:n_threshold
%                     current_phase_todo = 0;
%                     for current_phase=1:1:number_phase
%                         if p.todo(current_phase)
%                             current_phase_todo=current_phase_todo+1;
%                             if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C')
%                                 if Size_RVE(1,current_phase_todo,2,1)~=0
%                                     str_ = ['vf_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
%                                     str_(str_=='.')='p';
%                                     results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,1);
%                                 end
%                                 if strcmp(RVEparameters.type,'C')
%                                     if Size_RVE(1,current_phase_todo,2,2)~=0
%                                         str_ = ['vf_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
%                                         str_(str_=='.')='p';
%                                         results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
%                                     end
%                                 end
%                             else
%                                 if Size_convergence(1,current_phase_todo,2,1)~=0
%                                     str_ = ['vf_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
%                                     str_(str_=='.')='p';
%                                     results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,1);
%                                 end
%                                 if strcmp(RVEparameters.type,'F')
%                                     if Size_convergence(1,current_phase_todo,2,2)~=0
%                                         str_ = ['vf_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
%                                         str_(str_=='.')='p';
%                                         results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
% 
%             %% MANAGING RESULTS
%             [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,RVE,k_RVE,number_phase,number_phase_todo,infovol,p);
% 
%             if RVEparameters.donested
%                 for k_threshold=1:1:n_threshold
%                     current_phase_todo = 0;
%                     for current_phase=1:1:number_phase
%                         if p.todo(current_phase)
%                             current_phase_todo=current_phase_todo+1;
%                             if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C')
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1) = Size_RVE(k_threshold,current_phase_todo,1,1);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1) = Size_RVE(k_threshold,current_phase_todo,2,1);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1) = Size_RVE(k_threshold,current_phase_todo,3,1);
%                             else
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1) = Size_convergence(k_threshold,current_phase_todo,1,1);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1) = Size_convergence(k_threshold,current_phase_todo,2,1);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1) = Size_convergence(k_threshold,current_phase_todo,3,1);
%                             end
%                             if strcmp(RVEparameters.type,'C')
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2) = Size_RVE(k_threshold,current_phase_todo,1,2);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2) = Size_RVE(k_threshold,current_phase_todo,2,2);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2) = Size_RVE(k_threshold,current_phase_todo,3,2);
%                             elseif strcmp(RVEparameters.type,'F')
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3) = Size_convergence(k_threshold,current_phase_todo,1,2);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3) = Size_convergence(k_threshold,current_phase_todo,2,2);
%                                 Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3) = Size_convergence(k_threshold,current_phase_todo,3,2);
%                             end
%                         end
%                     end
%                 end
%             end
% 
%             %% TEXT DISPLAY AND SAVE RESULTS
%             if k_nestedRVE == 0
%                 propertyname='Volume fractions';
%                 RVEparameters.disp_parRVE = true;
%                 Function_subdomains_display_and_save(RVE,k_RVE,RVEparameters,number_phase,number_phase_todo,propertyname,Sub_folder_RVE,opt,infovol,p);
%             end
% 
%             %% FIGURES
%             if k_nestedRVE == 0
%                 parameters_figure.propertyname = propertyname;
%                 parameters_figure.propertynameunit = [];
%                 parameters_figure.RVE = RVEparameters;
%                 parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
%                 parameters_figure.savefolder = Sub_folder_RVE;
%                 parameters_figure.number_phase = number_phase;
%                 parameters_figure.number_phase_todo = number_phase_todo;
%                 parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
%                 parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
%                 parameters_figure.derivative_convergence = derivative_convergence;
%                 parameters_figure.relativedifference_convergence = relativedifference_convergence;
%                 parameters_figure.Size_RVE = Size_RVE;
%                 parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
%                 parameters_figure.Size_convergence = Size_convergence;
%                 parameters_figure.Wholevolume_size = Wholevolume_size;
%                 parameters_figure.Wholevolume_results = Volumefraction_phase;
%                 parameters_figure.infovol = infovol;
%                 parameters_figure.todo = p.todo;
%                 parameters_figure.opt = opt;
%                 Function_create_figures_RVE(parameters_figure) % Figures
%             end
%         end
% 
%         %% NESTED ANALYSIS RESULT
%         if RVEparameters.donested
%             % Table
%             [RVE] = Function_nestedtable(RVE,k_RVE,RVEparameters,number_phase,propertyname,Result_nestedRVE,Sub_folder_RVE,opt,infovol,p);
%             % Figure
%             parameters_figure.Result_nestedRVE = Result_nestedRVE;
%             Function_create_figures_nestedRVE(parameters_figure) % Figures
%             % Save
%             RVE(k_RVE).nestedanalysis = Result_nestedRVE;
%         end        
%     end    
%     Results_Volumefraction.RVE.volumefractions = RVE; % Save in main table result
% end

%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_connectivity.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_connectivity.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_connectivity.Table_date = Table_date;

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
    Function_Writetable(Current_folder,'Connectivity_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Connectivity_calculation_time', 'Connectivity', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_connectivity'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_connectivity'],'Results_connectivity')
end

%% SAVE VISUALIZATION
Current_folder = [infovol.volpath 'Visualization' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Visualization_connectivity'],'results_visualization','-v7.3');    


end