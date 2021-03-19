classdef Microstructure_generation_additives_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        Microstructure_additive_phase_generation  matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        LoadvolumeandvolumefractionTab  matlab.ui.container.Tab
        Load_instructions               matlab.ui.control.Label
        Load_title                      matlab.ui.control.Label
        Load_LoadvolumeButton           matlab.ui.control.Button
        Load_VolumeloaedLabel           matlab.ui.control.Label
        Load_ExistingphasesUITable      matlab.ui.control.Table
        Load_ExistingphasesLabel        matlab.ui.control.Label
        Load_NewphaseparameterLabel     matlab.ui.control.Label
        VolumefractiondenseLabel_3      matlab.ui.control.Label
        Load_volumefraction_dense       matlab.ui.control.NumericEditField
        NanoporosityLabel_3             matlab.ui.control.Label
        Load_nanoporosity               matlab.ui.control.NumericEditField
        TargetvolumefractionLabel_3     matlab.ui.control.Label
        Load_targetvolumefraction       matlab.ui.control.NumericEditField
        AdditiveidLabel                 matlab.ui.control.Label
        Load_Additive_id                matlab.ui.control.NumericEditField
        VisualizeloadedmicrostructureButton  matlab.ui.control.Button
        BetweenneighboursparticlesTab   matlab.ui.container.Tab
        Bridge_instructions             matlab.ui.control.Label
        Bridge_title                    matlab.ui.control.Label
        Bridge_NewphaseparameterLabel   matlab.ui.control.Label
        TargetvolumefractionLabel       matlab.ui.control.Label
        Bridge_targetvolumefraction     matlab.ui.control.NumericEditField
        Bridge_parameters_UITable       matlab.ui.control.Table
        Bridge_generation               matlab.ui.control.Button
        ProgressionGaugeLabel           matlab.ui.control.Label
        ProgressionGauge                matlab.ui.control.LinearGauge
        TextArea                        matlab.ui.control.TextArea
        AchievedEditFieldLabel          matlab.ui.control.Label
        Bridged_AchievedEditField       matlab.ui.control.NumericEditField
        EnergycriterionTab              matlab.ui.container.Tab
        PD_title                        matlab.ui.control.Label
        PD_instructions                 matlab.ui.control.Label
        PD_NewphaseparameterLabel       matlab.ui.control.Label
        PD_parameters_UITable           matlab.ui.control.Table
        AchievedEditFieldLabel_2        matlab.ui.control.Label
        PD_AchievedEditField            matlab.ui.control.NumericEditField
        TargetvolumefractionLabel_4     matlab.ui.control.Label
        PD_targetvolumefraction         matlab.ui.control.NumericEditField
        PD_generation                   matlab.ui.control.Button
        PD_TextArea                     matlab.ui.control.TextArea
        ProgressionGauge_2Label         matlab.ui.control.Label
        PD_ProgressionGauge             matlab.ui.control.LinearGauge
        OutcomeTab                      matlab.ui.container.Tab
        O_title                         matlab.ui.control.Label
        O_savebutton                    matlab.ui.control.Button
        O_vf_UITable                    matlab.ui.control.Table
        O_instructions                  matlab.ui.control.Label
        O_vf_Label                      matlab.ui.control.Label
        O_visualize                     matlab.ui.control.Button
        O_Label                         matlab.ui.control.Label
        O_FilenameLabel                 matlab.ui.control.Label
        O_filename_save                 matlab.ui.control.EditField
        O_savefolderButton              matlab.ui.control.Button
        O_savefolderLabel               matlab.ui.control.Label
        AboutTab                        matlab.ui.container.Tab
        O_title_2                       matlab.ui.control.Label
        About_TextArea                  matlab.ui.control.TextArea
        QuotationinstructionsTextAreaLabel  matlab.ui.control.Label
        About_Quotationinstructions     matlab.ui.control.TextArea
        About_Logo_NREL                 matlab.ui.control.Image
        OpendocumentationButton         matlab.ui.control.Button
        Github_ushbutton                matlab.ui.control.Button
        LinksLabel                      matlab.ui.control.Label
        Mistry_article                  matlab.ui.control.Button
    end

    
    properties (Access = public)
        Microstructure=[]; % Microstructure initialized empty
        Microstructure_initial=[]; % Microstructure initialized empty
        History_operations=[]; History_parameters=[]; History_elapsedtime=[]; % Log
        Savefolder=[]; % save folder
        filesave_ini=[]; % file savename
    end
    
    methods (Access = private)
        
        function [] = Turn_on_off(app,status_GUI,tab)
            child_handles = allchild(tab); % Get all children
            [n_children,~]=size(child_handles);
            for k=1:1:n_children
                child_handles(k).Visible=status_GUI;
                if isprop(child_handles(k),'Enable')
                    child_handles(k).Enable=status_GUI;
                end
            end
            
            if strcmp(tab.Title,'OutcomeTab')
                app.O_savebutton.Enable='off';
                app.O_visualize.Enable='off';
            end            
        end
        
        
        function [] = check_turnonoff_tab(app)
            idx = find( cell2mat(app.Load_ExistingphasesUITable.Data(:,3))==1 );
            target_volume_fraction = app.Load_targetvolumefraction.Value;
            if ~isempty(idx) && target_volume_fraction<=cell2mat(app.Load_ExistingphasesUITable.Data(idx,2))
                statut_tab = 'on';
            else
                statut_tab = 'off';
            end
            app.Turn_on_off(statut_tab,app.BetweenneighboursparticlesTab)
            app.Turn_on_off(statut_tab,app.EnergycriterionTab)    
            app.Turn_on_off(statut_tab,app.OutcomeTab) 
        end
        
        function [] = finished_generation(app)
            % Update volume fractions
            phases_before = unique(app.Microstructure_initial);
            phases_after = unique(app.Microstructure);
            phases = unique([phases_before; phases_after]);
            number_phase = length(phases);
            tabledata=zeros(number_phase,3);
            for current_phase=1:1:number_phase
                tabledata(current_phase,1)=phases(current_phase);
                tabledata(current_phase,2)=sum(sum(sum( app.Microstructure_initial==phases(current_phase)  ))) / numel(app.Microstructure_initial);
                tabledata(current_phase,3)=sum(sum(sum( app.Microstructure==phases(current_phase)  ))) / numel(app.Microstructure);
            end
            app.O_vf_UITable.Data = tabledata;
            app.O_savebutton.Enable='on';
            app.O_visualize.Enable='on';
            app.Turn_on_off('on',app.OutcomeTab)
            app.TabGroup.SelectedTab = app.OutcomeTab; % Go to outcome tab
        end
        
        function currentDir = getcurrentdir(app)
            if isdeployed % Stand-alone mode.
                [status, result] = system('path');
                currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
            else % MATLAB mode.
                currentDir = pwd;
            end
        end            
                    
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.Turn_on_off('off',app.BetweenneighboursparticlesTab)
            app.Turn_on_off('off',app.EnergycriterionTab)
            app.Turn_on_off('off',app.OutcomeTab)
        end

        % Button pushed function: Load_LoadvolumeButton
        function Load_LoadvolumeButtonPushed(app, event)
            % Open dialog box to choose file path
            [FileName,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)'},'Select a .tif file');
            if FileName==0
                % User clicked cancel button or closed the dialog box
                app.Load_VolumeloaedLabel.Text ='No volume loaded. Please load a volume with pre-existing phase(s).';
                app.O_savefolderLabel.Text = 'No save folder selected. Please select a folder.';
                app.Microstructure=[]; % Reset Microstructure
                app.Microstructure_initial=[]; % Reset Microstructure
                app.History_operations=[]; app.History_parameters=[]; app.History_elapsedtime=[]; % Reset log
                app.Savefolder=[];
                app.filesave_ini=[];
                app.Turn_on_off('off',app.BetweenneighboursparticlesTab)
                app.Turn_on_off('off',app.EnergycriterionTab)
                app.Turn_on_off('off',app.OutcomeTab)
                app.Load_ExistingphasesUITable.Data = [];
                
                app.VisualizeloadedmicrostructureButton.Enable = 'off';
                app.Load_volumefraction_dense.Enable = 'off'; app.VolumefractiondenseLabel_3.Enable = 'off';
                app.Load_nanoporosity.Enable = 'off'; app.NanoporosityLabel_3.Enable = 'off';
                app.Load_Additive_id.Enable = 'off'; app.AdditiveidLabel.Enable = 'off';
                
            else
                [Microstructure_local,outcome] = function_load_tif([PathName FileName]);
                if outcome % Success to import
                    
                    % Update global variables
                    app.Microstructure=Microstructure_local; % Update Microstructure
                    app.Microstructure_initial = app.Microstructure;
                    
                    %% Load tab
                    % Update save folder
                    app.Savefolder = PathName;
                    app.Load_VolumeloaedLabel.Text = [PathName FileName];
                    app.O_savefolderLabel.Text = app.Savefolder;
                    [~,tmp,~] = fileparts(FileName); % Remove .tif
                    app.filesave_ini = tmp;
                    app.Turn_on_off('on',app.BetweenneighboursparticlesTab)
                    app.Turn_on_off('on',app.EnergycriterionTab)
                    % Volume fraction and background
                    phases = unique(Microstructure_local);
                    number_phase = length(phases);
                    tabledata=cell(number_phase,3);
                    for current_phase=1:1:number_phase
                        tabledata(current_phase,1)={num2str(phases(current_phase))};
                        tabledata(current_phase,2)={sum(sum(sum( Microstructure_local==phases(current_phase)  ))) / numel(Microstructure_local)};
                        if phases(current_phase)==0
                            tabledata(current_phase,3)={true}; % By default, background is phase 0
                        else
                            tabledata(current_phase,3)={false};
                        end
                        
                    end
                    app.Load_ExistingphasesUITable.Data = tabledata;
                    
                    %% Between neighbours particles
                    tabledata=cell(3,3);
                    tabledata(1,:)={'Randomize cluster selection','[0,1]',0};
                    tabledata(2,:)={'Minimum distance from boundary in voxel length','[0,+Inf]',1};
                    tabledata(3,:)={'Erosion depth in voxel length','[0,+Inf] ',1};
                    app.Bridge_parameters_UITable.Data = tabledata;
                    
                    %% Preferential deposition
                    tabledata=cell(2,3);
                    tabledata(1,:)={'Morphology parameter w','[0.001,0.999]',0.5};
                    tabledata(2,:)={'Fraction of additive phase, relative to target volume fraction, added at each step ','[0,1]',0.01};
                    app.PD_parameters_UITable.Data = tabledata;                    
                    
                    app.VisualizeloadedmicrostructureButton.Enable = 'on';
                    app.Load_volumefraction_dense.Enable = 'on';  app.VolumefractiondenseLabel_3.Enable = 'on';
                    app.Load_nanoporosity.Enable = 'on'; app.NanoporosityLabel_3.Enable = 'on';
                    app.Load_Additive_id.Enable = 'on'; app.AdditiveidLabel.Enable = 'on';
                    app.Load_Additive_id.Value = double(max(phases)+1);
                    
                    app.Load_volumefraction_denseValueChanged;
                    
                else
                    app.Load_VolumeloaedLabel.Text ='Error: volume cannot be loaded! File is not a .tif file';
                    app.Load_savefolderLabel.Text = 'No save folder selected. Please select a folder.';
                    app.Microstructure=[]; % Reset Microstructure
                    app.Microstructure_initial=[]; % Reset Microstructure
                    app.History_operations=[]; app.History_parameters=[]; app.History_elapsedtime=[]; % Reset log
                    app.Savefolder=[];
                    app.filesave_ini=[];
                    app.Turn_on_off('off',app.BetweenneighboursparticlesTab)
                    app.Turn_on_off('off',app.EnergycriterionTab)
                    app.Turn_on_off('off',app.OutcomeTab)
                    app.Load_ExistingphasesUITable.Data = [];
                    app.VisualizeloadedmicrostructureButton.Enable = 'off';
                    app.Load_volumefraction_dense.Enable = 'off';  app.VolumefractiondenseLabel_3.Enable = 'off';
                    app.Load_nanoporosity.Enable = 'off'; app.NanoporosityLabel_3.Enable = 'off';
                    app.Load_Additive_id.Enable = 'off'; app.AdditiveidLabel.Enable = 'off';
                end
            end
            
            
        end

        % Cell edit callback: Load_ExistingphasesUITable
        function Load_ExistingphasesUITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(2)==3
                tabledata = app.Load_ExistingphasesUITable.Data;
                tabledata(:,3) = {false};
                if newData==1
                    tabledata(indices(1),3)={true};
                end
                app.Load_ExistingphasesUITable.Data = tabledata;
            end
            app.Load_volumefraction_denseValueChanged;
        end

        % Button pushed function: Bridge_generation
        function Bridge_generationButtonPushed(app, event)
            app.Turn_on_off('off',app.OutcomeTab) 
                        
            % Prepare parameters
            idx = find( cell2mat(app.Load_ExistingphasesUITable.Data(:,3))==1 );
            app.Microstructure = app.Microstructure_initial; % Reset
            background_id = str2num(cell2mat(app.Load_ExistingphasesUITable.Data(idx,1)));
            target_volume_fraction = app.Load_targetvolumefraction.Value;
            additive_id = app.Load_Additive_id.Value;
            randomize_level = cell2mat(app.Bridge_parameters_UITable.Data(1,3));
            minimum_dist = cell2mat(app.Bridge_parameters_UITable.Data(2,3));
            erosiondilatation_depth = cell2mat(app.Bridge_parameters_UITable.Data(3,3));
            
            % You could use it as a function as
            %function [app.Microstructure] = generate_additive(app.Microstructure,background_id,additive_id,target_volume_fraction,minimum_dist,randomize_level,erosiondilatation_depth)
            
            time_cpu_start = cputime; % CPU start
            tic; % Stopwatch start
            
            app.ProgressionGauge.Value = 0; pause(0.001); % Update progression bar
            app.TextArea.Value = {'Calculating','distance to boundary...'}; pause(0.01);
            
            binary_phase=zeros(size(app.Microstructure)); % Initialization
            binary_phase(app.Microstructure == background_id) = 1;
            [distance_bnd2bnd] = Function_particle_size_CPSD_Algorithm(binary_phase); % Calculate distance from boundary to boundary
            
%            % Visualization step-by-step
%             figure
%             axis equal tight
%             imagesc(distance_bnd2bnd(:,:,2))  
%             set(gca,'YDir','normal')
%             custommap = jet; custommap(1,:)=[0 0 0];
%             colormap(custommap)
            
            distance_bnd2bnd(distance_bnd2bnd==0)=9e9; % Remove 0 (background) from the potential additive location
            distance_bnd2bnd(distance_bnd2bnd<=minimum_dist*sqrt(3))=9e9; % Remove small features from the potential additive location
            
            app.TextArea.Value = {'Generating','additive phase...'}; pause(0.01);
                        
            ntot = numel(app.Microstructure);
            current_additive_volumefraction = 0; % Initialize
            increment = 1;
            iter=0;
            d=0; % Initialize
            while current_additive_volumefraction < target_volume_fraction
                d=d+increment; % Increment distance
                iter=iter+1;
                if d>minimum_dist*sqrt(3)
                    potential_additives = distance_bnd2bnd<=d; % Select all voxels that could be assigned to additive
                    n = sum(sum(sum(potential_additives))); % How many of them ?
                    if n>0 % At least one
                        if randomize_level==0
                            if n/ntot <= target_volume_fraction-current_additive_volumefraction % All potential voxel are not enough to match target volume fraction
                                app.Microstructure(potential_additives)=additive_id;
                            else % Randomize among possible location is mandatory not to overshoot the target volume fraction
                                [L,ncluster] = bwlabeln(potential_additives,18); % Identify cluster among these potential locations
                                cluster_possibility=1:1:ncluster;
                                while current_additive_volumefraction <= target_volume_fraction
                                    cluster_choice=randi(length(cluster_possibility)); % Random selection of cluster
                                    app.Microstructure(L==cluster_possibility(cluster_choice))=additive_id;
                                    current_additive_volumefraction = sum(sum(sum(app.Microstructure==additive_id)))/ntot; % Update current volume fraction
                                    cluster_possibility(cluster_choice)=[]; % Avoid re-assignment
                                end
                            end
                        else % Use a portion of them, per cluster randomly chosen
                            [L,ncluster] = bwlabeln(potential_additives,18); % Identify cluster among these potential locations
                            cluster_possibility=1:1:ncluster;
                            volume_increment=0; % initalize
                            while volume_increment <= (1-randomize_level)*n
                                cluster_choice=randi(length(cluster_possibility)); % Random selection of cluster
                                loc = L==cluster_possibility(cluster_choice);
                                app.Microstructure(loc)=additive_id;
                                volume_increment = volume_increment + sum(sum(sum( loc )));
                                cluster_possibility(cluster_choice)=[]; % Avoid re-assignment
                            end
                        end
                    end
                end
                distance_bnd2bnd(app.Microstructure==additive_id)=9e9; % Remove already assigned additive from the potential additive new location
                current_additive_volumefraction = sum(sum(sum(app.Microstructure==additive_id)))/ntot; % Update current volume fraction
                
                app.ProgressionGauge.Value = 100*current_additive_volumefraction/target_volume_fraction; pause(0.001); % Update progression bar
                
                progress(iter,1) = current_additive_volumefraction;
                if iter>10 && length(unique(progress(end-5:end,1)))==1
                    progress
                    break
                end
                
%                %Visualization step-by-step
%                 figure
%                 axis equal tight
%                 imagesc(app.Microstructure(:,:,2))
%                 set(gca,'YDir','normal')
%                 custommap = copper; custommap(1,:)=[0 0 0];
%                 colormap(custommap)
                
                
            end
            app.ProgressionGauge.Value = 100; pause(0.001); % Update progression bar
                        
            if erosiondilatation_depth>0 % Optional step
                app.TextArea.Value = {'Erosion of additive phase...'}; pause(0.01);
                binary_phase=zeros(size(app.Microstructure)); % Initialization
                binary_phase(app.Microstructure == additive_id) = 1;
                [distance_bnd2bnd] = Function_particle_size_CPSD_Algorithm(binary_phase); % Calculate distance from boundary to boundary
                distance_bnd2bnd(distance_bnd2bnd==0)=9e9; % Remove 0 (background) from the potential additive location
                app.Microstructure(distance_bnd2bnd<(erosiondilatation_depth*sqrt(3)+0.0001))=background_id;                
            end
             
% %             Visualization step-by-step
%             figure
%             axis equal tight
%             imagesc(app.Microstructure(:,:,2))
%             set(gca,'YDir','normal')
%             custommap = copper; custommap(1,:)=[0 0 0];
%             colormap(custommap)
            
            if minimum_dist>0
                app.TextArea.Value = {'Correcting minimum','distance...'}; pause(0.01);
                binary_phase=zeros(size(app.Microstructure)); % Initialization
                binary_phase(app.Microstructure == background_id) = 1;
                [distance_bnd2bnd] = Function_particle_size_CPSD_Algorithm(binary_phase); % Calculate distance from boundary to boundary
                distance_bnd2bnd(distance_bnd2bnd==0)=9e9; % Remove 0 (background) from the potential additive location
                app.Microstructure(distance_bnd2bnd<minimum_dist*sqrt(3)+0.0001)=additive_id;
            end
            
% %             Visualization step-by-step
%             figure
%             axis equal tight
%             imagesc(app.Microstructure(:,:,2))
%             set(gca,'YDir','normal')
%             custommap = copper; custommap(1,:)=[0 0 0];
%             colormap(custommap)            
    

            delta = round(sum(sum(sum( app.Microstructure==additive_id  ))) - numel(app.Microstructure)*target_volume_fraction)
            erosion_cluster = true;
            while delta>0
                %erosion_cluster
                %delta
                app.TextArea.Value = {'Resorbing volume fraction','difference...'}; pause(0.01);
                binary_phase=zeros(size(app.Microstructure)); % Create binary matrix
                if erosion_cluster==false
                    % Random voxel erosion
                    binary_phase(app.Microstructure==background_id) = 1;
                    distance_map = bwdist(binary_phase,'chessboard'); % Distance map
                    distance_map(app.Microstructure ~= additive_id)=2; % Consider only additive
                    idx = find(distance_map<2);
                    if length(idx)<=delta
                        app.Microstructure(idx)=background_id;
                    else
                        tmp = randi(length(idx),delta,1);
                        app.Microstructure(idx(tmp))=background_id;
                    end
                else
                    % Random cluster erosion
                    binary_phase(app.Microstructure==additive_id) = 1;
                    [L,~] = bwlabeln(binary_phase,18);
                    % Determine cluster size
                    [C,~,ic] = unique(L);
                    a_counts = accumarray(ic,1);
                    size_cluster = [C, a_counts];
                    size_cluster = sortrows(size_cluster,2); % Sort in ascending order
                    while size_cluster(1,2)<delta && delta>0 
                        app.Microstructure(L==size_cluster(1,1))=background_id; % Smallest cluster is removed
                        delta = round(sum(sum(sum( app.Microstructure==additive_id  ))) - numel(app.Microstructure)*target_volume_fraction);
                        size_cluster(1,:)=[]; % Smallest cluster is removed from the list 
                    end
                    erosion_cluster = false; % Switch to voxel erosion
                end
                delta = round(sum(sum(sum( app.Microstructure==additive_id  ))) - numel(app.Microstructure)*target_volume_fraction);
            end
           
              
% %             Visualization step-by-step
%             figure
%             axis equal tight
%             imagesc(app.Microstructure(:,:,2))
%             set(gca,'YDir','normal')
%             custommap = copper; custommap(1,:)=[0 0 0];
%             colormap(custommap)               
            
            % CPU and stopwatch time - end
            time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
            time_stopwatch_elapsed = toc; % Stopwatch elapsed time
            
            app.TextArea.Value = {[num2str(ntot,'%i'),' voxels (domain)'],['Stopwatch time:',num2str(time_stopwatch_elapsed,'%1.1f'), 's'],['CPU time:',num2str(time_cpu_elapsed,'%1.1f'), 's']}; pause(0.01);
            app.Bridged_AchievedEditField.Value = sum(sum(sum( app.Microstructure==additive_id  ))) / numel(app.Microstructure);
            
            % Done
            app.O_filename_save.Value = [app.filesave_ini '_bridge.tif'];
            app.finished_generation
        end

        % Cell edit callback: Bridge_parameters_UITable
        function Bridge_parameters_UITableCellEdit(app, event)
            % Check parameters
            tabledata = app.Bridge_parameters_UITable.Data;
            v=cell2mat(tabledata(1,3));
            if isnan(v) || v<0 || v>1
                tabledata(1,3) = {0};
            end            
            v=cell2mat(tabledata(2,3));
            if isnan(v) || v<0
                tabledata(2,3) = {1};
            end                  
            v=cell2mat(tabledata(3,3));
            if isnan(v) || v<0 || v~=round(v)
                tabledata(3,3) = {1};
            end           
            app.Bridge_parameters_UITable.Data = tabledata;
        end

        % Value changed function: Load_nanoporosity, 
        % Load_volumefraction_dense
        function Load_volumefraction_denseValueChanged(app, event)
            volumefraction_dense = app.Load_volumefraction_dense.Value;
            nanoporosity = app.Load_nanoporosity.Value;
            target_volume_fraction = volumefraction_dense * 1/(1-nanoporosity);
            if ~isempty(app.Load_ExistingphasesUITable.Data)
                idx = find( cell2mat(app.Load_ExistingphasesUITable.Data(:,3))==1 );
                if ~isempty(idx)
                    max_vf = cell2mat(app.Load_ExistingphasesUITable.Data(idx,2));
                else
                    max_vf = 9e9;
                end
                app.Load_targetvolumefraction.Value = min([max_vf, target_volume_fraction]);
            else
                app.Load_targetvolumefraction.Value = target_volume_fraction;
            end
            app.Bridge_targetvolumefraction.Value = app.Load_targetvolumefraction.Value;
            app.PD_targetvolumefraction.Value = app.Load_targetvolumefraction.Value;
        end

        % Value changed function: Load_Additive_id
        function Load_Additive_idValueChanged(app, event)
            additive_id = app.Load_Additive_id.Value;
            existing_id = str2num(cell2mat(app.Load_ExistingphasesUITable.Data(:,1)));
            if sum(ismember(existing_id,additive_id))>0
                app.Load_Additive_id.Value = max(existing_id)+1;
            end
        end

        % Button pushed function: O_savebutton
        function O_savebuttonButtonPushed(app, event)
            function_save_tif(app.Microstructure, [app.Savefolder app.O_filename_save.Value]);
        end

        % Button pushed function: PD_generation
        function PD_generationButtonPushed(app, event)
            app.Turn_on_off('off',app.OutcomeTab) 
                        
            % Prepare parameters
            idx = find( cell2mat(app.Load_ExistingphasesUITable.Data(:,3))==1 );
            app.Microstructure = app.Microstructure_initial; % Reset
            background_id = str2num(cell2mat(app.Load_ExistingphasesUITable.Data(idx,1)));
            target_volume_fraction = app.Load_targetvolumefraction.Value;
            additive_id = app.Load_Additive_id.Value;
            morphology_parameter = cell2mat(app.PD_parameters_UITable.Data(1,3));
            fractionstep = cell2mat(app.PD_parameters_UITable.Data(2,3));
             
            time_cpu_start = cputime; % CPU start
            tic; % Stopwatch start
            
            app.PD_ProgressionGauge.Value = 0; pause(0.001); % Update progression bar
            app.PD_TextArea.Value = {'Calculating...'}; pause(0.001);
            
            % Prepare data set
            I  = zeros(size(app.Microstructure));
            I(app.Microstructure~=background_id) = 1; % 0, background, 1 other phases, 2 additive phase
            
            %__________________________________________________________________________
            % Programmed by     : Aashutosh Mistry (START)
            % Created on        : Jul 05, 2020
            % Article ref.      : https://doi.org/10.1021/acsami.7b17771
            %   Mistry, Smith & Mukherjee (2018) Secondary Phase Stochastics in
            %   Lithium-Ion Battery Electrodes, ACS Appl. Mater. Interfaces 10(7) pp.
            %   6317-6326
            
            [Nx, Ny, Nz] = size(I);                     % voxel ranges
            N = Nx*Ny*Nz;                               % total voxels
            Nsimultaneous = round(N*target_volume_fraction*fractionstep); % how many new CBD voxels are in a single run
            % 1D index for each of the 3D indices
            Pid = zeros(Nx,Ny,Nz);
            for k=1:Nz
                for j=1:Ny
                    for i=1:Nx
                        Pid(i,j,k) = Index1D(i,j,k, Nx,Ny,Nz);
                    end
                end
            end
            
            % generating a neighbor list - ensuring periodicity of boundaries _________
            [Eid, Wid, Nid, Sid, Uid, Lid] = Neighbors(Nx,Ny,Nz);
            % multiple rounds of deposition
            e2 = 0; % Initialization
            while e2<target_volume_fraction
                I = AddCBDwMorphology_step (I, morphology_parameter, Pid, Eid,Wid,Nid,Sid,Uid,Lid, Nsimultaneous);
                e2 = sum(sum(sum( I==2 )))/N;  % new CBD volume fraction
                app.PD_AchievedEditField.Value = e2;
                app.PD_ProgressionGauge.Value = 100*e2/target_volume_fraction; pause(0.001); % Update progression bar
            end
            app.Microstructure(I==2)= additive_id;
            % Programmed by     : Aashutosh Mistry (END)
            %__________________________________________________________________________
            
            % CPU and stopwatch time - end
            time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
            time_stopwatch_elapsed = toc; % Stopwatch elapsed time
            
            app.PD_TextArea.Value = {[num2str(N,'%i'),' voxels (domain)'],['Stopwatch time:',num2str(time_stopwatch_elapsed,'%1.1f'), 's'],['CPU time:',num2str(time_cpu_elapsed,'%1.1f'), 's']}; pause(0.01);
            app.PD_AchievedEditField.Value = sum(sum(sum( app.Microstructure==additive_id  ))) / numel(app.Microstructure);
            
            % Done
            app.O_filename_save.Value = [app.filesave_ini '_w' num2str(morphology_parameter,'%1.3f') '.tif'];
            app.finished_generation            
        end

        % Button pushed function: O_visualize
        function O_visualizeButtonPushed(app, event)
            direction_name = {'Normal to Axis 1';'Normal to Axis 2';'Normal to Axis 3'};
            Microstructure_basic_visualization_interface(app.Microstructure, direction_name);            
        end

        % Button pushed function: 
        % VisualizeloadedmicrostructureButton
        function VisualizeloadedmicrostructureButtonPushed(app, event)
            direction_name = {'Normal to Axis 1';'Normal to Axis 2';'Normal to Axis 3'};
            Microstructure_basic_visualization_interface(app.Microstructure_initial, direction_name);
        end

        % Cell edit callback: PD_parameters_UITable
        function PD_parameters_UITableCellEdit(app, event)
            % Check parameters
            tabledata = app.PD_parameters_UITable.Data;
            v=cell2mat(tabledata(1,3));
            if isnan(v) || v<0.001 || v>0.999
                tabledata(1,3) = {0.5};
            end            
            v=cell2mat(tabledata(2,3));
            if isnan(v) || v<0 || v>1
                tabledata(2,3) = {0.1};
            end                  
            app.PD_parameters_UITable.Data = tabledata;            
        end

        % Button pushed function: O_savefolderButton
        function O_savefolderButtonPushed(app, event)
            if isempty(app.Savefolder)
                selpath = uigetdir(pwd,'Select save folder');
            else
                selpath = uigetdir(app.Savefolder,'Select save folder');
            end
            if selpath~=0
                if ispc
                    selpath=[selpath '\'];
                else
                    selpath=[selpath '/'];
                end
                app.Savefolder=selpath;
                app.O_savefolderLabel.Text = app.Savefolder;
            end
        end

        % Button pushed function: OpendocumentationButton
        function OpendocumentationButtonPushed(app, event)
            % path = matlab.desktop.editor.getActiveFilename; % Path of active file (but does not work for app file)
            path = app.getcurrentdir;
            if ispc
                separation_folder = '\';
            else
                separation_folder = '/';
            end
            higherlevelfolder = extractBetween(path,path(1:5),['MATBOX_Microstructure_analysis_toolbox' separation_folder],'Boundaries','inclusive');
            documentation_path = [char(higherlevelfolder) 'Documentation' separation_folder 'NREL_MATBOX_Microstructure_analysis_toolbox_documentation.pdf'];
            if exist(documentation_path,'file')
                open(documentation_path);
            else
                disp 'MATLAB did not find the file NREL_MATBOX_Microstructure_analysis_toolbox_documentation.pdf'.
                disp 'Default location is \MATBOX_Microstructure_analysis_toolbox\Documentation\';
            end 
        end

        % Image clicked function: About_Logo_NREL
        function About_Logo_NRELImageClicked(app, event)
            url = 'https://www.nrel.gov/transportation/';
            web(url)
        end

        % Button pushed function: Github_ushbutton
        function Github_ushbuttonPushed(app, event)
            url = 'https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox';
            web(url)
        end

        % Button pushed function: Mistry_article
        function Mistry_articlePushed(app, event)
            url = 'https://doi.org/10.1021/acsami.7b17771';
            web(url)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create Microstructure_additive_phase_generation and hide until all components are created
            app.Microstructure_additive_phase_generation = uifigure('Visible', 'off');
            app.Microstructure_additive_phase_generation.Position = [100 100 958 590];
            app.Microstructure_additive_phase_generation.Name = 'Additive phase generation';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.Microstructure_additive_phase_generation);
            app.TabGroup.TabLocation = 'left';
            app.TabGroup.Position = [2 1 957 590];

            % Create LoadvolumeandvolumefractionTab
            app.LoadvolumeandvolumefractionTab = uitab(app.TabGroup);
            app.LoadvolumeandvolumefractionTab.Title = 'Load volume and volume fraction';

            % Create Load_instructions
            app.Load_instructions = uilabel(app.LoadvolumeandvolumefractionTab);
            app.Load_instructions.FontAngle = 'italic';
            app.Load_instructions.Position = [13 505 718 42];
            app.Load_instructions.Text = {'Instructions: load a volume, select background and choose your target volume fraction for the additive phase.'; 'Then move to the tab that corresponds to the type of additive/method you want to generate.'};

            % Create Load_title
            app.Load_title = uilabel(app.LoadvolumeandvolumefractionTab);
            app.Load_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Load_title.HorizontalAlignment = 'center';
            app.Load_title.FontWeight = 'bold';
            app.Load_title.Position = [13 557 718 22];
            app.Load_title.Text = 'Volume with pre-existing phases';

            % Create Load_LoadvolumeButton
            app.Load_LoadvolumeButton = uibutton(app.LoadvolumeandvolumefractionTab, 'push');
            app.Load_LoadvolumeButton.ButtonPushedFcn = createCallbackFcn(app, @Load_LoadvolumeButtonPushed, true);
            app.Load_LoadvolumeButton.Position = [13 473 100 22];
            app.Load_LoadvolumeButton.Text = 'Load volume';

            % Create Load_VolumeloaedLabel
            app.Load_VolumeloaedLabel = uilabel(app.LoadvolumeandvolumefractionTab);
            app.Load_VolumeloaedLabel.Position = [129 473 602 22];
            app.Load_VolumeloaedLabel.Text = 'No volume loaded. Please load a volume with pre-existing phase(s).';

            % Create Load_ExistingphasesUITable
            app.Load_ExistingphasesUITable = uitable(app.LoadvolumeandvolumefractionTab);
            app.Load_ExistingphasesUITable.ColumnName = {'Id'; 'Volume fraction'; 'Background'};
            app.Load_ExistingphasesUITable.RowName = {};
            app.Load_ExistingphasesUITable.ColumnEditable = [false false true];
            app.Load_ExistingphasesUITable.CellEditCallback = createCallbackFcn(app, @Load_ExistingphasesUITableCellEdit, true);
            app.Load_ExistingphasesUITable.Position = [13 283 381 142];

            % Create Load_ExistingphasesLabel
            app.Load_ExistingphasesLabel = uilabel(app.LoadvolumeandvolumefractionTab);
            app.Load_ExistingphasesLabel.FontWeight = 'bold';
            app.Load_ExistingphasesLabel.Position = [13 432 405 22];
            app.Load_ExistingphasesLabel.Text = 'Existing phases. Select background where additive will be generated';

            % Create Load_NewphaseparameterLabel
            app.Load_NewphaseparameterLabel = uilabel(app.LoadvolumeandvolumefractionTab);
            app.Load_NewphaseparameterLabel.FontWeight = 'bold';
            app.Load_NewphaseparameterLabel.Position = [13 240 137 22];
            app.Load_NewphaseparameterLabel.Text = 'Additive phase';

            % Create VolumefractiondenseLabel_3
            app.VolumefractiondenseLabel_3 = uilabel(app.LoadvolumeandvolumefractionTab);
            app.VolumefractiondenseLabel_3.HorizontalAlignment = 'right';
            app.VolumefractiondenseLabel_3.Enable = 'off';
            app.VolumefractiondenseLabel_3.Position = [13 211 132 22];
            app.VolumefractiondenseLabel_3.Text = 'Volume fraction (dense)';

            % Create Load_volumefraction_dense
            app.Load_volumefraction_dense = uieditfield(app.LoadvolumeandvolumefractionTab, 'numeric');
            app.Load_volumefraction_dense.Limits = [0 1];
            app.Load_volumefraction_dense.ValueChangedFcn = createCallbackFcn(app, @Load_volumefraction_denseValueChanged, true);
            app.Load_volumefraction_dense.Enable = 'off';
            app.Load_volumefraction_dense.Position = [158 211 102 22];
            app.Load_volumefraction_dense.Value = 0.1;

            % Create NanoporosityLabel_3
            app.NanoporosityLabel_3 = uilabel(app.LoadvolumeandvolumefractionTab);
            app.NanoporosityLabel_3.HorizontalAlignment = 'right';
            app.NanoporosityLabel_3.Enable = 'off';
            app.NanoporosityLabel_3.Position = [13 182 76 22];
            app.NanoporosityLabel_3.Text = 'Nanoporosity';

            % Create Load_nanoporosity
            app.Load_nanoporosity = uieditfield(app.LoadvolumeandvolumefractionTab, 'numeric');
            app.Load_nanoporosity.Limits = [0 1];
            app.Load_nanoporosity.ValueChangedFcn = createCallbackFcn(app, @Load_volumefraction_denseValueChanged, true);
            app.Load_nanoporosity.Enable = 'off';
            app.Load_nanoporosity.Position = [158 182 102 22];

            % Create TargetvolumefractionLabel_3
            app.TargetvolumefractionLabel_3 = uilabel(app.LoadvolumeandvolumefractionTab);
            app.TargetvolumefractionLabel_3.HorizontalAlignment = 'right';
            app.TargetvolumefractionLabel_3.Position = [13 153 124 22];
            app.TargetvolumefractionLabel_3.Text = 'Target volume fraction';

            % Create Load_targetvolumefraction
            app.Load_targetvolumefraction = uieditfield(app.LoadvolumeandvolumefractionTab, 'numeric');
            app.Load_targetvolumefraction.Limits = [0 1];
            app.Load_targetvolumefraction.Editable = 'off';
            app.Load_targetvolumefraction.Position = [158 153 102 22];
            app.Load_targetvolumefraction.Value = 0.1;

            % Create AdditiveidLabel
            app.AdditiveidLabel = uilabel(app.LoadvolumeandvolumefractionTab);
            app.AdditiveidLabel.HorizontalAlignment = 'right';
            app.AdditiveidLabel.Enable = 'off';
            app.AdditiveidLabel.Position = [76 121 61 22];
            app.AdditiveidLabel.Text = 'Additive id';

            % Create Load_Additive_id
            app.Load_Additive_id = uieditfield(app.LoadvolumeandvolumefractionTab, 'numeric');
            app.Load_Additive_id.Limits = [0 Inf];
            app.Load_Additive_id.RoundFractionalValues = 'on';
            app.Load_Additive_id.ValueChangedFcn = createCallbackFcn(app, @Load_Additive_idValueChanged, true);
            app.Load_Additive_id.Enable = 'off';
            app.Load_Additive_id.Position = [158 121 102 22];
            app.Load_Additive_id.Value = 2;

            % Create VisualizeloadedmicrostructureButton
            app.VisualizeloadedmicrostructureButton = uibutton(app.LoadvolumeandvolumefractionTab, 'push');
            app.VisualizeloadedmicrostructureButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizeloadedmicrostructureButtonPushed, true);
            app.VisualizeloadedmicrostructureButton.Enable = 'off';
            app.VisualizeloadedmicrostructureButton.Position = [408 403 182 22];
            app.VisualizeloadedmicrostructureButton.Text = 'Visualize loaded microstructure';

            % Create BetweenneighboursparticlesTab
            app.BetweenneighboursparticlesTab = uitab(app.TabGroup);
            app.BetweenneighboursparticlesTab.Title = 'Between neighbours particles';
            app.BetweenneighboursparticlesTab.ForegroundColor = [0 0 1];

            % Create Bridge_instructions
            app.Bridge_instructions = uilabel(app.BetweenneighboursparticlesTab);
            app.Bridge_instructions.FontAngle = 'italic';
            app.Bridge_instructions.Position = [13 525 718 22];
            app.Bridge_instructions.Text = 'Instructions: Choose algorithms parameters and click generate. Once done, you can go to the outcome tab for viusalization and saving.';

            % Create Bridge_title
            app.Bridge_title = uilabel(app.BetweenneighboursparticlesTab);
            app.Bridge_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Bridge_title.HorizontalAlignment = 'center';
            app.Bridge_title.FontWeight = 'bold';
            app.Bridge_title.Position = [13 557 718 22];
            app.Bridge_title.Text = 'Additive phase will be generated so that they will connect neighbours particles, effectively forming bridges';

            % Create Bridge_NewphaseparameterLabel
            app.Bridge_NewphaseparameterLabel = uilabel(app.BetweenneighboursparticlesTab);
            app.Bridge_NewphaseparameterLabel.FontWeight = 'bold';
            app.Bridge_NewphaseparameterLabel.Position = [13 493 137 22];
            app.Bridge_NewphaseparameterLabel.Text = 'New phase parameters';

            % Create TargetvolumefractionLabel
            app.TargetvolumefractionLabel = uilabel(app.BetweenneighboursparticlesTab);
            app.TargetvolumefractionLabel.HorizontalAlignment = 'right';
            app.TargetvolumefractionLabel.Position = [515 493 124 22];
            app.TargetvolumefractionLabel.Text = 'Target volume fraction';

            % Create Bridge_targetvolumefraction
            app.Bridge_targetvolumefraction = uieditfield(app.BetweenneighboursparticlesTab, 'numeric');
            app.Bridge_targetvolumefraction.Limits = [0 1];
            app.Bridge_targetvolumefraction.Editable = 'off';
            app.Bridge_targetvolumefraction.Position = [660 493 71 22];
            app.Bridge_targetvolumefraction.Value = 0.1;

            % Create Bridge_parameters_UITable
            app.Bridge_parameters_UITable = uitable(app.BetweenneighboursparticlesTab);
            app.Bridge_parameters_UITable.ColumnName = {'Parameter'; 'Range'; 'Value'};
            app.Bridge_parameters_UITable.RowName = {};
            app.Bridge_parameters_UITable.ColumnEditable = [false false true];
            app.Bridge_parameters_UITable.CellEditCallback = createCallbackFcn(app, @Bridge_parameters_UITableCellEdit, true);
            app.Bridge_parameters_UITable.Position = [13 335 437 151];

            % Create Bridge_generation
            app.Bridge_generation = uibutton(app.BetweenneighboursparticlesTab, 'push');
            app.Bridge_generation.ButtonPushedFcn = createCallbackFcn(app, @Bridge_generationButtonPushed, true);
            app.Bridge_generation.BackgroundColor = [0.0745 0.6235 1];
            app.Bridge_generation.FontSize = 16;
            app.Bridge_generation.FontWeight = 'bold';
            app.Bridge_generation.FontColor = [1 1 1];
            app.Bridge_generation.Enable = 'off';
            app.Bridge_generation.Position = [13 13 96 87];
            app.Bridge_generation.Text = 'Generate';

            % Create ProgressionGaugeLabel
            app.ProgressionGaugeLabel = uilabel(app.BetweenneighboursparticlesTab);
            app.ProgressionGaugeLabel.HorizontalAlignment = 'center';
            app.ProgressionGaugeLabel.Position = [122 46 70 22];
            app.ProgressionGaugeLabel.Text = 'Progression';

            % Create ProgressionGauge
            app.ProgressionGauge = uigauge(app.BetweenneighboursparticlesTab, 'linear');
            app.ProgressionGauge.MajorTicks = [0 10 20 30 40 50 60 70 80 90 100];
            app.ProgressionGauge.Position = [203 13 351 87];

            % Create TextArea
            app.TextArea = uitextarea(app.BetweenneighboursparticlesTab);
            app.TextArea.Position = [573 13 158 87];
            app.TextArea.Value = {'Not running'};

            % Create AchievedEditFieldLabel
            app.AchievedEditFieldLabel = uilabel(app.BetweenneighboursparticlesTab);
            app.AchievedEditFieldLabel.HorizontalAlignment = 'right';
            app.AchievedEditFieldLabel.Position = [515 466 55 22];
            app.AchievedEditFieldLabel.Text = 'Achieved';

            % Create Bridged_AchievedEditField
            app.Bridged_AchievedEditField = uieditfield(app.BetweenneighboursparticlesTab, 'numeric');
            app.Bridged_AchievedEditField.Editable = 'off';
            app.Bridged_AchievedEditField.Position = [660 466 71 22];

            % Create EnergycriterionTab
            app.EnergycriterionTab = uitab(app.TabGroup);
            app.EnergycriterionTab.Title = 'Energy criterion';
            app.EnergycriterionTab.ForegroundColor = [0 0 1];

            % Create PD_title
            app.PD_title = uilabel(app.EnergycriterionTab);
            app.PD_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.PD_title.HorizontalAlignment = 'center';
            app.PD_title.FontWeight = 'bold';
            app.PD_title.Position = [13 557 718 22];
            app.PD_title.Text = 'Additive phase preferenital deposition is controlled through a energy criterion';

            % Create PD_instructions
            app.PD_instructions = uilabel(app.EnergycriterionTab);
            app.PD_instructions.FontAngle = 'italic';
            app.PD_instructions.Position = [13 525 718 22];
            app.PD_instructions.Text = 'Instructions: Choose algorithms parameters and click generate. Once done, you can go to the outcome tab for viusalization and saving.';

            % Create PD_NewphaseparameterLabel
            app.PD_NewphaseparameterLabel = uilabel(app.EnergycriterionTab);
            app.PD_NewphaseparameterLabel.FontWeight = 'bold';
            app.PD_NewphaseparameterLabel.Position = [13 493 137 22];
            app.PD_NewphaseparameterLabel.Text = 'New phase parameters';

            % Create PD_parameters_UITable
            app.PD_parameters_UITable = uitable(app.EnergycriterionTab);
            app.PD_parameters_UITable.ColumnName = {'Parameter'; 'Range'; 'Value'};
            app.PD_parameters_UITable.RowName = {};
            app.PD_parameters_UITable.ColumnEditable = [false false true];
            app.PD_parameters_UITable.CellEditCallback = createCallbackFcn(app, @PD_parameters_UITableCellEdit, true);
            app.PD_parameters_UITable.Position = [13 335 718 151];

            % Create AchievedEditFieldLabel_2
            app.AchievedEditFieldLabel_2 = uilabel(app.EnergycriterionTab);
            app.AchievedEditFieldLabel_2.HorizontalAlignment = 'right';
            app.AchievedEditFieldLabel_2.Position = [13 272 55 22];
            app.AchievedEditFieldLabel_2.Text = 'Achieved';

            % Create PD_AchievedEditField
            app.PD_AchievedEditField = uieditfield(app.EnergycriterionTab, 'numeric');
            app.PD_AchievedEditField.Editable = 'off';
            app.PD_AchievedEditField.Position = [158 272 71 22];

            % Create TargetvolumefractionLabel_4
            app.TargetvolumefractionLabel_4 = uilabel(app.EnergycriterionTab);
            app.TargetvolumefractionLabel_4.HorizontalAlignment = 'right';
            app.TargetvolumefractionLabel_4.Position = [13 299 124 22];
            app.TargetvolumefractionLabel_4.Text = 'Target volume fraction';

            % Create PD_targetvolumefraction
            app.PD_targetvolumefraction = uieditfield(app.EnergycriterionTab, 'numeric');
            app.PD_targetvolumefraction.Limits = [0 1];
            app.PD_targetvolumefraction.Editable = 'off';
            app.PD_targetvolumefraction.Position = [158 299 71 22];
            app.PD_targetvolumefraction.Value = 0.1;

            % Create PD_generation
            app.PD_generation = uibutton(app.EnergycriterionTab, 'push');
            app.PD_generation.ButtonPushedFcn = createCallbackFcn(app, @PD_generationButtonPushed, true);
            app.PD_generation.BackgroundColor = [0.0745 0.6235 1];
            app.PD_generation.FontSize = 16;
            app.PD_generation.FontWeight = 'bold';
            app.PD_generation.FontColor = [1 1 1];
            app.PD_generation.Enable = 'off';
            app.PD_generation.Position = [13 13 96 87];
            app.PD_generation.Text = 'Generate';

            % Create PD_TextArea
            app.PD_TextArea = uitextarea(app.EnergycriterionTab);
            app.PD_TextArea.Position = [573 13 158 87];
            app.PD_TextArea.Value = {'Not running'};

            % Create ProgressionGauge_2Label
            app.ProgressionGauge_2Label = uilabel(app.EnergycriterionTab);
            app.ProgressionGauge_2Label.HorizontalAlignment = 'center';
            app.ProgressionGauge_2Label.Position = [122 46 70 22];
            app.ProgressionGauge_2Label.Text = 'Progression';

            % Create PD_ProgressionGauge
            app.PD_ProgressionGauge = uigauge(app.EnergycriterionTab, 'linear');
            app.PD_ProgressionGauge.MajorTicks = [0 10 20 30 40 50 60 70 80 90 100];
            app.PD_ProgressionGauge.Position = [203 13 351 87];

            % Create OutcomeTab
            app.OutcomeTab = uitab(app.TabGroup);
            app.OutcomeTab.Title = 'Outcome';

            % Create O_title
            app.O_title = uilabel(app.OutcomeTab);
            app.O_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.O_title.HorizontalAlignment = 'center';
            app.O_title.FontWeight = 'bold';
            app.O_title.Position = [13 557 718 22];
            app.O_title.Text = 'Check generation result and save';

            % Create O_savebutton
            app.O_savebutton = uibutton(app.OutcomeTab, 'push');
            app.O_savebutton.ButtonPushedFcn = createCallbackFcn(app, @O_savebuttonButtonPushed, true);
            app.O_savebutton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.O_savebutton.FontSize = 16;
            app.O_savebutton.FontWeight = 'bold';
            app.O_savebutton.FontColor = [1 1 1];
            app.O_savebutton.Enable = 'off';
            app.O_savebutton.Position = [17 15 96 38];
            app.O_savebutton.Text = 'Save';

            % Create O_vf_UITable
            app.O_vf_UITable = uitable(app.OutcomeTab);
            app.O_vf_UITable.ColumnName = {'Id'; 'Inital'; 'After generation'};
            app.O_vf_UITable.RowName = {};
            app.O_vf_UITable.Position = [13 335 381 151];

            % Create O_instructions
            app.O_instructions = uilabel(app.OutcomeTab);
            app.O_instructions.FontAngle = 'italic';
            app.O_instructions.Position = [13 525 718 22];
            app.O_instructions.Text = 'Instructions: Verify volume fractions are those expected after generation. Select save folder and write filename, then save.';

            % Create O_vf_Label
            app.O_vf_Label = uilabel(app.OutcomeTab);
            app.O_vf_Label.FontWeight = 'bold';
            app.O_vf_Label.Position = [13 493 137 22];
            app.O_vf_Label.Text = 'Volume fractions';

            % Create O_visualize
            app.O_visualize = uibutton(app.OutcomeTab, 'push');
            app.O_visualize.ButtonPushedFcn = createCallbackFcn(app, @O_visualizeButtonPushed, true);
            app.O_visualize.BackgroundColor = [0.502 0.502 0.502];
            app.O_visualize.FontSize = 16;
            app.O_visualize.FontWeight = 'bold';
            app.O_visualize.FontColor = [1 1 1];
            app.O_visualize.Enable = 'off';
            app.O_visualize.Position = [408 448 96 38];
            app.O_visualize.Text = 'Visualize';

            % Create O_Label
            app.O_Label = uilabel(app.OutcomeTab);
            app.O_Label.Position = [13 286 718 28];
            app.O_Label.Text = {'If you are not satisifed with the generated additive phase, you can try with other parameters or methods.'; 'The microstructure is reset each time you click on a ''Generate'' button.'};

            % Create O_FilenameLabel
            app.O_FilenameLabel = uilabel(app.OutcomeTab);
            app.O_FilenameLabel.HorizontalAlignment = 'right';
            app.O_FilenameLabel.Position = [13 63 86 22];
            app.O_FilenameLabel.Text = 'Save Filename';

            % Create O_filename_save
            app.O_filename_save = uieditfield(app.OutcomeTab, 'text');
            app.O_filename_save.Position = [129 63 602 22];
            app.O_filename_save.Value = 'none';

            % Create O_savefolderButton
            app.O_savefolderButton = uibutton(app.OutcomeTab, 'push');
            app.O_savefolderButton.ButtonPushedFcn = createCallbackFcn(app, @O_savefolderButtonPushed, true);
            app.O_savefolderButton.Position = [13 92 100 22];
            app.O_savefolderButton.Text = 'Save folder';

            % Create O_savefolderLabel
            app.O_savefolderLabel = uilabel(app.OutcomeTab);
            app.O_savefolderLabel.Position = [129 92 602 22];
            app.O_savefolderLabel.Text = 'No save folder selected. Please select a folder.';

            % Create AboutTab
            app.AboutTab = uitab(app.TabGroup);
            app.AboutTab.Title = 'About';

            % Create O_title_2
            app.O_title_2 = uilabel(app.AboutTab);
            app.O_title_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.O_title_2.HorizontalAlignment = 'center';
            app.O_title_2.FontWeight = 'bold';
            app.O_title_2.Position = [13 557 718 22];
            app.O_title_2.Text = 'About the additive generation module';

            % Create About_TextArea
            app.About_TextArea = uitextarea(app.AboutTab);
            app.About_TextArea.Position = [13 403 718 141];
            app.About_TextArea.Value = {'This module enables you to generate an additive phase within the backgroud phase of a microstructure, without modifying the existing solid phases.'; ''; 'The generation can be performed through two different approaches:'; ''; '1) Additive is preferentially located between neighboring particles.'; ''; '2) Additive is preferentially located either at the surface of existing particle surface or at the surface of already deposited additive, based on a user-defined energy criteria.'};

            % Create QuotationinstructionsTextAreaLabel
            app.QuotationinstructionsTextAreaLabel = uilabel(app.AboutTab);
            app.QuotationinstructionsTextAreaLabel.HorizontalAlignment = 'right';
            app.QuotationinstructionsTextAreaLabel.FontWeight = 'bold';
            app.QuotationinstructionsTextAreaLabel.Position = [13 363 134 22];
            app.QuotationinstructionsTextAreaLabel.Text = 'Quotation instructions';

            % Create About_Quotationinstructions
            app.About_Quotationinstructions = uitextarea(app.AboutTab);
            app.About_Quotationinstructions.Position = [162 240 569 148];
            app.About_Quotationinstructions.Value = {'- For any additive generated with this module:'; 'F. L. E. Usseglio-Viretta et al., MATBOX: An Open-source Microstructure Analysis Toolbox for microstructure generation, segmentation, characterization, visualization, correlation, and meshing, SoftwareX, in preparation'; ''; '- If you generate the phase using the energy criterion, please ALSO quote:'; 'A. N. Mistry, K. Smith, and P. P. Mukherjee, Secondary Phase Stochastics in Lithium-Ion Battery Electrodes, ACS Appl. Mater. Interfaces 10(7) pp. 6317-6326 (2018), https://doi.org/10.1021/acsami.7b17771'; ''};

            % Create About_Logo_NREL
            app.About_Logo_NREL = uiimage(app.AboutTab);
            app.About_Logo_NREL.ImageClickedFcn = createCallbackFcn(app, @About_Logo_NRELImageClicked, true);
            app.About_Logo_NREL.Position = [13 32 264 100];
            app.About_Logo_NREL.ImageSource = 'logo_NREL.png';

            % Create OpendocumentationButton
            app.OpendocumentationButton = uibutton(app.AboutTab, 'push');
            app.OpendocumentationButton.ButtonPushedFcn = createCallbackFcn(app, @OpendocumentationButtonPushed, true);
            app.OpendocumentationButton.BackgroundColor = [0.8 0.8 0.8];
            app.OpendocumentationButton.Position = [584 187 147 38];
            app.OpendocumentationButton.Text = 'Open documentation';

            % Create Github_ushbutton
            app.Github_ushbutton = uibutton(app.AboutTab, 'push');
            app.Github_ushbutton.ButtonPushedFcn = createCallbackFcn(app, @Github_ushbuttonPushed, true);
            app.Github_ushbutton.HorizontalAlignment = 'left';
            app.Github_ushbutton.BackgroundColor = [0.8 0.8 0.8];
            app.Github_ushbutton.Position = [295 84 441 40];
            app.Github_ushbutton.Text = 'Github repository: MATBOX';

            % Create LinksLabel
            app.LinksLabel = uilabel(app.AboutTab);
            app.LinksLabel.FontWeight = 'bold';
            app.LinksLabel.Position = [13 134 37 22];
            app.LinksLabel.Text = 'Links';

            % Create Mistry_article
            app.Mistry_article = uibutton(app.AboutTab, 'push');
            app.Mistry_article.ButtonPushedFcn = createCallbackFcn(app, @Mistry_articlePushed, true);
            app.Mistry_article.HorizontalAlignment = 'left';
            app.Mistry_article.BackgroundColor = [0.8 0.8 0.8];
            app.Mistry_article.Position = [295 39 441 40];
            app.Mistry_article.Text = 'Journal article: Secondary Phase Stochastics in Lithium-Ion Battery Electrodes';

            % Show the figure after all components are created
            app.Microstructure_additive_phase_generation.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Microstructure_generation_additives_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.Microstructure_additive_phase_generation)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.Microstructure_additive_phase_generation)
        end
    end
end