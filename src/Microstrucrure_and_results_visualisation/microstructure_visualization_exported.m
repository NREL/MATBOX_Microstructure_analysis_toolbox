classdef microstructure_visualization_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        Visualizationmodule_UIFigure    matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        InstructionsTab                 matlab.ui.container.Tab
        Main_instructions_title         matlab.ui.control.Label
        instructions_1                  matlab.ui.control.Label
        instructions_2                  matlab.ui.control.Label
        instructions_3                  matlab.ui.control.Label
        ViewvolumesTab                  matlab.ui.container.Tab
        Checkseg_title_2                matlab.ui.control.Label
        Viewvolume_Selecttif_Button     matlab.ui.control.Button
        Checkseg_instructions_2         matlab.ui.control.Label
        View_UITable                    matlab.ui.control.Table
        View_Labelloaded                matlab.ui.control.Label
        View_createfigure_button        matlab.ui.control.Button
        View_3Dfigure_button            matlab.ui.control.Button
        IdleLabel                       matlab.ui.control.Label
        VoxelsizeEditFieldLabel_2       matlab.ui.control.Label
        View_VoxelsizeEditField         matlab.ui.control.NumericEditField
        UnitEditFieldLabel_2            matlab.ui.control.Label
        View_UnitEditField              matlab.ui.control.EditField
        ChecksegmentationTab            matlab.ui.container.Tab
        Checkseg_title                  matlab.ui.control.Label
        Checkseg_Importtif_grey_Button  matlab.ui.control.Button
        Checkseg_Importtif_segmented_Button  matlab.ui.control.Button
        Checkseg_instructions           matlab.ui.control.Label
        Checkseg_createfigure_button    matlab.ui.control.Button
        Checkseg_Greyvolumename         matlab.ui.control.Label
        Checkseg_Segmentedvolumename    matlab.ui.control.Label
        VoxelsizeEditFieldLabel         matlab.ui.control.Label
        Checkseg_VoxelsizeEditField     matlab.ui.control.NumericEditField
        UnitEditFieldLabel              matlab.ui.control.Label
        Checkseg_UnitEditField          matlab.ui.control.EditField
        ViewcharacterizationresultsTab  matlab.ui.container.Tab
        Res_title                       matlab.ui.control.Label
        Res_Selectfolder_Button         matlab.ui.control.Button
        Res_instructions                matlab.ui.control.Label
        Res_Labelloaded                 matlab.ui.control.Label
        Res_createfigure_button         matlab.ui.control.Button
        Res_3Dfigure_button             matlab.ui.control.Button
        VoxelsizeEditFieldLabel_3       matlab.ui.control.Label
        Res_VoxelsizeEditField          matlab.ui.control.NumericEditField
        UnitEditFieldLabel_3            matlab.ui.control.Label
        Res_UnitEditField               matlab.ui.control.EditField
        PropertytoplotDropDownLabel     matlab.ui.control.Label
        PropertytoplotDropDown          matlab.ui.control.DropDown
        ForthephaseDropDownLabel        matlab.ui.control.Label
        ForthephaseDropDown             matlab.ui.control.DropDown
        Res_ErrormessageLabel           matlab.ui.control.Label
        PropertyunitLabel               matlab.ui.control.Label
        Res_ResUnitEditField            matlab.ui.control.EditField
    end

    
    properties (Access = private)
        % View volumes
        allpath = {};
        uniquevolume = [];
        volumes = [];
        
        % Check segmentation
        is_loaded_segmentedvolume = false;
        is_loaded_greyvolume = false;
        is_domainsize_coherent = false;
        segmentedvolume = [];
        greyvolume = [];
        
        % Characterization result
        results_filenames = [];
        current_result_folder = [];
        arrayresult = [];
        resultname = [];
        result_id_mat =[];
    end
    
    methods (Access = private)
        
        

    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.View_UITable.Visible = 'off'; app.View_UITable.Enable = 'off';
            app.Res_ErrormessageLabel.Visible = 'off';
        end

        % Button pushed function: Checkseg_Importtif_grey_Button
        function Checkseg_Importtif_grey_ButtonPushed(app, event)
            % Open dialog box to choose file path
            str_dialogbox = 'Select grey-level volume';
            [FileName,PathName,~] = uigetfile({'*.tif;*.tiff','Tif image (*.tif, *.tiff)'},str_dialogbox);
            if FileName==0
                % User clicked cancel button or closed the dialog box
                set(app.Checkseg_Greyvolumename,'FontColor','k','Text','No grey-level volume loaded');
                app.is_loaded_greyvolume = false;
            else
                [app.greyvolume, outcome] = function_load_tif([PathName FileName], 'uint16' ); % Load new volume
                if outcome % Success to import
                    app.is_loaded_greyvolume = true;
                    % Check dimension
                    if app.is_loaded_greyvolume && app.is_loaded_segmentedvolume
                        size_segmented = size(app.segmentedvolume);
                        dim_segmented = length(size_segmented);
                        size_grey = size(app.greyvolume);
                        dim_grey = length(size_grey);
                        if dim_segmented==dim_grey && sum(size_segmented==size_grey)==dim_segmented
                            app.is_domainsize_coherent=true;
                            set(app.Checkseg_Greyvolumename,'FontColor','k','Text',FileName);
                        else
                            app.is_domainsize_coherent=false;
                            set(app.Checkseg_Greyvolumename,'FontColor','r','Text','Volumes have different fields of view!');
                        end
                    else
                        set(app.Checkseg_Greyvolumename,'FontColor','k','Text',FileName);
                    end
                else
                    set(app.Checkseg_Greyvolumename,'FontColor','r','Text','Failed to load the volume!');
                    app.is_loaded_greyvolume = false;
                end
            end
            if app.is_loaded_greyvolume && app.is_loaded_segmentedvolume && app.is_domainsize_coherent
                set(app.Checkseg_createfigure_button,'enable','on');
            else
                set(app.Checkseg_createfigure_button,'enable','off');
            end
        end

        % Button pushed function: 
        % Checkseg_Importtif_segmented_Button
        function Checkseg_Importtif_segmented_ButtonPushed(app, event)
            % Open dialog box to choose file path
            str_dialogbox = 'Select segmented volume';
            [FileName,PathName,~] = uigetfile({'*.tif;*.tiff','Tif image (*.tif, *.tiff)'},str_dialogbox);
            if FileName==0
                % User clicked cancel button or closed the dialog box
                set(app.Checkseg_Segmentedvolumename,'FontColor','k','Text','No segmented volume loaded');
                app.is_loaded_segmentedvolume = false;
            else
                [app.segmentedvolume, outcome] = function_load_tif([PathName FileName], 'uint8' ); % Load new volume
                if outcome % Success to import
                    app.is_loaded_segmentedvolume = true;
                    % Check dimension
                    if app.is_loaded_greyvolume && app.is_loaded_segmentedvolume
                        size_segmented = size(app.segmentedvolume);
                        dim_segmented = length(size_segmented);
                        size_grey = size(app.greyvolume);
                        dim_grey = length(size_grey);    
                        if dim_segmented==dim_grey && sum(size_segmented==size_grey)==dim_segmented
                            app.is_domainsize_coherent=true;
                            set(app.Checkseg_Segmentedvolumename,'FontColor','k','Text',FileName);
                        else
                            app.is_domainsize_coherent=false;
                            set(app.Checkseg_Segmentedvolumename,'FontColor','r','Text','Volumes have different fields of view!');
                        end
                    else
                        set(app.Checkseg_Segmentedvolumename,'FontColor','k','Text',FileName);
                    end
                else
                    set(app.Checkseg_Segmentedvolumename,'FontColor','r','Text','Failed to load the volume!');
                    app.is_loaded_segmentedvolume = false;
                end
            end
            if app.is_loaded_greyvolume && app.is_loaded_segmentedvolume && app.is_domainsize_coherent
                set(app.Checkseg_createfigure_button,'enable','on');
            else
                set(app.Checkseg_createfigure_button,'enable','off');
            end
        end

        % Button pushed function: Checkseg_createfigure_button
        function Checkseg_createfigure_buttonButtonPushed(app, event)
            microstructure_visualization_segmentation(app.greyvolume, app.segmentedvolume, app.Checkseg_Greyvolumename.Text, app.Checkseg_Segmentedvolumename.Text, app.Checkseg_VoxelsizeEditField.Value, app.Checkseg_UnitEditField.Value);
        end

        % Button pushed function: Viewvolume_Selecttif_Button
        function Viewvolume_Selecttif_ButtonPushed(app, event)
            flg=true;   % set the logic variable to start the loop
            nfile=0; files={};
            app.allpath = {}; % reset
            while flg
                str_dialogbox = 'Select One or More Files. Click cancel to stop importing more files';
                [FileName,PathName,~] = uigetfile({'*.tif;*.tiff','Tif image (*.tif, *.tiff)'},str_dialogbox,'MultiSelect', 'on');
                if isequal(FileName,0), flg=0; break; end   %  no file selected; quit
                if iscell(FileName) % More than 1 file selected
                    [~,n]=size(FileName);
                    for k=1:1:n
                        nfile=nfile+1;
                        files{nfile,1} = [char(FileName(k))];
                        app.allpath{nfile,1} = [char(PathName) char(FileName(k))];
                    end
                else
                    nfile=nfile+1;
                    files{nfile,1} = [char(FileName)];
                    app.allpath{nfile,1} = [char(PathName) char(FileName)];
                end
            end
            if ~isempty(app.allpath)
                if nfile>1
                    % Import
                    app.IdleLabel.Text = 'Importing...'; pause(0.01);                    
                    allimport = true;
                    for k=1:1:nfile
                        [app.volumes(k).array, outcome] = function_load_tif(app.allpath{k,1}, 'uint16' ); % Load new volume
                        if outcome % Success to import
                            set(app.View_Labelloaded,'FontColor','k','Text',[files{k,:} ' imported.']);  pause(0.01);
                        else
                            set(app.View_Labelloaded,'FontColor','r','Text',[files{k,:} ' import failed!']);
                            allimport = false;
                            break                            
                        end
                    end
                    app.IdleLabel.Text = 'Idle'; pause(0.01);
                    if allimport
                        set(app.View_Labelloaded,'FontColor','k','Text',[num2str(nfile) ' volumes imported.']);
                        app.View_UITable.Data = [files num2cell(ones(nfile,1))];
                        app.View_UITable.ColumnFormat = {'char' 'logical'};
                        app.View_UITable.Visible = 'on'; app.View_UITable.Enable = 'on';
                        app.View_createfigure_button.Enable = 'on';
                    else
                        app.View_UITable.Visible = 'off'; app.View_UITable.Enable = 'off';
                        app.View_createfigure_button.Enable = 'off';
                    end
                    app.View_3Dfigure_button.Enable = 'off';
                else
                    app.View_UITable.Visible = 'off'; app.View_UITable.Enable = 'off';
                    % Import
                    app.IdleLabel.Text = 'Importing...'; pause(0.01);
                    [app.uniquevolume, outcome] = function_load_tif(app.allpath{nfile,1}, 'uint16' ); % Load new volume
                    if outcome % Success to import
                        set(app.View_Labelloaded,'FontColor','k','Text',files{1,:});
                        app.View_createfigure_button.Enable = 'on';
                        app.View_3Dfigure_button.Enable = 'on';
                    else
                        set(app.View_Labelloaded,'FontColor','r','Text','Failed to load the volume!');
                        app.View_createfigure_button.Enable = 'off';
                        app.View_3Dfigure_button.Enable = 'off';
                    end
                    app.IdleLabel.Text = 'Idle'; pause(0.01);
                end
            else
                app.View_Labelloaded.Text = 'No volume selected';
                app.View_createfigure_button.Enable = 'off';
                app.View_3Dfigure_button.Enable = 'off';
                app.View_UITable.Visible = 'off'; app.View_UITable.Enable = 'off';
            end
            
            
            
        end

        % Button pushed function: View_createfigure_button
        function View_createfigure_buttonButtonPushed(app, event)
            [nfile,~]=size(app.allpath);
            if nfile==1 % 1 volume visualization
                voxelsize = app.View_VoxelsizeEditField.Value;
                voxelsize_unitname = app.View_UnitEditField.Value;
                arrayunit = '';
                microstructure_visualization_slices(app.uniquevolume, app.View_Labelloaded.Text, voxelsize, voxelsize_unitname, arrayunit);
            else 
                choice = cell2mat(app.View_UITable.Data(:,2));
                if sum(choice)>1 % n volume visualization
                    % Build the string require to call the visualization function
                    str = 'Microstructure_comparison_visualization_interface(';
                    first_volume = true;
                    for k=1:1:length(choice)
                        if choice(k)
                            if first_volume
                                str = [str 'app.volumes(' num2str(k) ').array'];
                                first_volume = false;
                            else
                                str = [str ',app.volumes(' num2str(k) ').array'];
                            end
                        end
                    end
                    str = [str ');'];
                    eval(str);
                else % 1 volume visualization
                    voxelsize = app.View_VoxelsizeEditField.Value;
                    voxelsize_unitname = app.View_UnitEditField.Value;
                    arrayunit = '';
                    idx = find(choice==1);
                    microstructure_visualization_slices(app.volumes(idx).array, app.View_UITable.Data(idx,1), voxelsize, voxelsize_unitname, arrayunit);
                end
            end            
        end

        % Button pushed function: View_3Dfigure_button
        function View_3Dfigure_buttonButtonPushed(app, event)
            Fig_3D=figure;
            col_ = bone;
            volshow(app.uniquevolume,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');
        end

        % Button pushed function: Res_Selectfolder_Button
        function Res_Selectfolder_ButtonPushed(app, event)
            str_dialogbox = 'Select the main folder of the volume you want to investigate'; % Set string of the dialog box
            volume_main_folder = uigetdir(pwd,str_dialogbox); % Open dialog box to choose file path
            if volume_main_folder==0
                % User clicked cancel button or closed the dialog box
                app.Res_Labelloaded.Text = 'No folder selected';
                app.PropertytoplotDropDown.Enable = 'off'; app.PropertytoplotDropDownLabel.Enable = 'off';
                app.ForthephaseDropDown.Enable = 'off'; app.ForthephaseDropDownLabel.Enable = 'off';
                app.Res_ErrormessageLabel.Visible = 'off';
                app.Res_createfigure_button.Enable = 'off';
                app.Res_3Dfigure_button.Enable = 'off';
            else
                if ispc
                    app.current_result_folder = [volume_main_folder '\' 'Visualization' '\'];
                    tmp_ = find(volume_main_folder=='\');
                else
                    app.current_result_folder = [volume_main_folder '/' 'Visualization' '/'];
                    tmp_ = find(volume_main_folder=='/');
                end                
                volume_name = volume_main_folder(tmp_(end)+1:length(volume_main_folder));
                % Find all mat files that start with 'Visualization'
                MyFolderInfo=dir(app.current_result_folder);
                n_files=length(MyFolderInfo);
                detected_files=0;
                app.results_filenames={};
                for k_file=1:1:n_files
                    [~,filename,ext] = fileparts(MyFolderInfo(k_file).name);
                    if strcmp(ext,'.mat')
                        sub_name = filename(1:13);
                        if strcmp(sub_name,'Visualization')
                            detected_files=detected_files+1;
                            app.results_filenames(detected_files)={[filename '.mat']};
                        end
                    end
                end
                if detected_files
                    app.Res_ErrormessageLabel.Visible = 'off'; app.Res_ErrormessageLabel.Text = 'Error message';
                    app.Res_Labelloaded.Text = volume_name; % Set name
                    number_results = length(app.results_filenames); % Number of results to summarize
                    fields=[];
                    for current_matfile=1:1:number_results % Loop over all saved results
                        current_result=char(app.results_filenames(current_matfile)); % Load result
                        pathname = [app.current_result_folder current_result];
                        datamat = load(pathname);
                        if isfield(datamat,'results_visualization')
                            data = datamat.results_visualization;
                            fields = [fields; fieldnames(data)];
                        end
                    end
                    fields = unique(fields);
                    field_results={};
                    kk=0;
                    for k=1:1:length(fields)
                        if ~strcmp(char(fields(k)),'name')
                            kk=kk+1;
                            field_results(kk,1)=fields(k);
                        end
                    end
                    app.PropertytoplotDropDownLabel.Enable = 'on';
                    set(app.PropertytoplotDropDown,'Items', field_results,'enable','on');
                    app.PropertytoplotDropDownValueChanged;
                else
                    app.Res_ErrormessageLabel.Visible = 'on'; app.Res_ErrormessageLabel.Text = 'Wrong folder! ...\Main_folder_to_select\Visualization\Visualization*.mat';
                    app.Res_Labelloaded.Text = 'No folder selected';
                    app.PropertytoplotDropDown.Enable = 'off'; app.PropertytoplotDropDownLabel.Enable = 'off';
                    app.ForthephaseDropDown.Enable = 'off'; app.ForthephaseDropDownLabel.Enable = 'off';
                    app.Res_ErrormessageLabel.Visible = 'off';
                    app.Res_createfigure_button.Enable = 'off';
                    app.Res_3Dfigure_button.Enable = 'off';
                end
            end
        end

        % Button pushed function: Res_createfigure_button
        function Res_createfigure_buttonButtonPushed(app, event)
            voxelsize = app.Res_VoxelsizeEditField.Value;
            voxelsize_unitname = app.Res_UnitEditField.Value;
            arrayunit = app.Res_ResUnitEditField.Value;
            microstructure_visualization_slices(app.arrayresult, app.resultname, voxelsize, voxelsize_unitname, arrayunit);
        end

        % Button pushed function: Res_3Dfigure_button
        function Res_3Dfigure_buttonButtonPushed(app, event)
            Fig_3D=figure;
            col_ = turbo;
            %volshow(app.arrayresult,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');
            volshow(app.arrayresult,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','MaximumIntensityProjection');
        end

        % Value changed function: PropertytoplotDropDown
        function PropertytoplotDropDownValueChanged(app, event)
            property_string = app.PropertytoplotDropDown.Value;
            number_results = length(app.results_filenames); % Number of results to summarize
            app.result_id_mat=[];
            for current_matfile=1:1:number_results % Loop over all saved results
                current_result=char(app.results_filenames(current_matfile)); % Load result
                pathname = [app.current_result_folder current_result];
                datamat = load(pathname);
                if isfield(datamat,'results_visualization')
                    data = datamat.results_visualization;
                    fields = fieldnames(data);
                    for k=1:1:length(fields)
                        if strcmp( char(fields(k)),property_string)
                            app.result_id_mat=current_matfile;
                            [~, number_phase] = size(data);
                            phase_results={};
                            for kk=1:1:number_phase
                                phase_results(kk,1)={data(kk).name};
                            end
                            set(app.ForthephaseDropDown,'Items', phase_results,'enable','on');
                            app.ForthephaseDropDownLabel.Enable = 'on';
                            app.ForthephaseDropDownValueChanged;
                        end
                    end
                end
            end
            app.ForthephaseDropDown;
        end

        % Value changed function: ForthephaseDropDown
        function ForthephaseDropDownValueChanged(app, event)
            phase_string = app.ForthephaseDropDown.Value;
            volume_name = app.Res_Labelloaded.Text;
            property_string = app.PropertytoplotDropDown.Value;
            app.resultname = [volume_name ', ' property_string ', ' phase_string];     
            current_result=char(app.results_filenames(app.result_id_mat)); % Load result
            pathname = [app.current_result_folder current_result];
            datamat = load(pathname);
            data = datamat.results_visualization;
            % Select phase
            [~, number_phase] = size(data);
            for k=1:1:number_phase
                if strcmp(data(k).name, phase_string)
                    k_phase=k;
                end
            end
            % Select array
            app.arrayresult = data(k_phase).(property_string);
            app.Res_createfigure_button.Enable = 'on';
            if length(size(app.arrayresult))==3
                app.Res_3Dfigure_button.Enable = 'on';
            else
                app.Res_3Dfigure_button.Enable = 'off';
            end
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create Visualizationmodule_UIFigure and hide until all components are created
            app.Visualizationmodule_UIFigure = uifigure('Visible', 'off');
            app.Visualizationmodule_UIFigure.Position = [100 100 858 590];
            app.Visualizationmodule_UIFigure.Name = 'Visualization module';
            app.Visualizationmodule_UIFigure.Icon = 'Icon_visualization.png';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.Visualizationmodule_UIFigure);
            app.TabGroup.TabLocation = 'left';
            app.TabGroup.Position = [1 1 858 590];

            % Create InstructionsTab
            app.InstructionsTab = uitab(app.TabGroup);
            app.InstructionsTab.Title = 'Instructions';

            % Create Main_instructions_title
            app.Main_instructions_title = uilabel(app.InstructionsTab);
            app.Main_instructions_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Main_instructions_title.HorizontalAlignment = 'center';
            app.Main_instructions_title.FontWeight = 'bold';
            app.Main_instructions_title.Position = [9 559 657 22];
            app.Main_instructions_title.Text = 'Main instructions';

            % Create instructions_1
            app.instructions_1 = uilabel(app.InstructionsTab);
            app.instructions_1.Position = [9 527 657 22];
            app.instructions_1.Text = '- If you want to visualize one or several tif files, go to the "View volume(s)" tab';

            % Create instructions_2
            app.instructions_2 = uilabel(app.InstructionsTab);
            app.instructions_2.Position = [9 490 657 22];
            app.instructions_2.Text = '- If you want to evaluate visusally the relevance of a segmentation, go to the "Check segmentation" tab';

            % Create instructions_3
            app.instructions_3 = uilabel(app.InstructionsTab);
            app.instructions_3.Position = [9 419 657 56];
            app.instructions_3.Text = {'- If you want to visualize a voxel-wise property (such as particle diameter, particle lable, cluster connectivity etc.)'; 'go to the the "View characterization results" tab. Note that you need to run the characterization module first to generate'; 'the results to visualize in this tab. If you do not find the property you want to visualize you can easily add it (see'; 'the documentation for details).'};

            % Create ViewvolumesTab
            app.ViewvolumesTab = uitab(app.TabGroup);
            app.ViewvolumesTab.Title = 'View volume(s)';

            % Create Checkseg_title_2
            app.Checkseg_title_2 = uilabel(app.ViewvolumesTab);
            app.Checkseg_title_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Checkseg_title_2.HorizontalAlignment = 'center';
            app.Checkseg_title_2.FontWeight = 'bold';
            app.Checkseg_title_2.Position = [9 559 657 22];
            app.Checkseg_title_2.Text = 'Select one or several volume files';

            % Create Viewvolume_Selecttif_Button
            app.Viewvolume_Selecttif_Button = uibutton(app.ViewvolumesTab, 'push');
            app.Viewvolume_Selecttif_Button.ButtonPushedFcn = createCallbackFcn(app, @Viewvolume_Selecttif_ButtonPushed, true);
            app.Viewvolume_Selecttif_Button.BackgroundColor = [0.8 0.8 0.8];
            app.Viewvolume_Selecttif_Button.Position = [9 479 200 33];
            app.Viewvolume_Selecttif_Button.Text = 'Click to select one or several files';

            % Create Checkseg_instructions_2
            app.Checkseg_instructions_2 = uilabel(app.ViewvolumesTab);
            app.Checkseg_instructions_2.Position = [9 527 657 22];
            app.Checkseg_instructions_2.Text = 'Select file(s), enter a voxel size and then create figure. Volumes can have different dimension.';

            % Create View_UITable
            app.View_UITable = uitable(app.ViewvolumesTab);
            app.View_UITable.ColumnName = {'Files'; 'Select volumes to plot'};
            app.View_UITable.ColumnWidth = {'1x', 'auto'};
            app.View_UITable.RowName = {};
            app.View_UITable.ColumnEditable = [false true];
            app.View_UITable.Position = [9 284 657 185];

            % Create View_Labelloaded
            app.View_Labelloaded = uilabel(app.ViewvolumesTab);
            app.View_Labelloaded.Position = [219 485 447 22];
            app.View_Labelloaded.Text = 'No volume selected';

            % Create View_createfigure_button
            app.View_createfigure_button = uibutton(app.ViewvolumesTab, 'push');
            app.View_createfigure_button.ButtonPushedFcn = createCallbackFcn(app, @View_createfigure_buttonButtonPushed, true);
            app.View_createfigure_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.View_createfigure_button.FontSize = 14;
            app.View_createfigure_button.FontWeight = 'bold';
            app.View_createfigure_button.FontColor = [1 1 1];
            app.View_createfigure_button.Enable = 'off';
            app.View_createfigure_button.Position = [9 24 164 33];
            app.View_createfigure_button.Text = '2D slices visualization';

            % Create View_3Dfigure_button
            app.View_3Dfigure_button = uibutton(app.ViewvolumesTab, 'push');
            app.View_3Dfigure_button.ButtonPushedFcn = createCallbackFcn(app, @View_3Dfigure_buttonButtonPushed, true);
            app.View_3Dfigure_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.View_3Dfigure_button.FontSize = 14;
            app.View_3Dfigure_button.FontWeight = 'bold';
            app.View_3Dfigure_button.FontColor = [1 1 1];
            app.View_3Dfigure_button.Enable = 'off';
            app.View_3Dfigure_button.Position = [188 24 133 33];
            app.View_3Dfigure_button.Text = '3D visualization';

            % Create IdleLabel
            app.IdleLabel = uilabel(app.ViewvolumesTab);
            app.IdleLabel.Position = [336 29 330 22];
            app.IdleLabel.Text = 'Idle';

            % Create VoxelsizeEditFieldLabel_2
            app.VoxelsizeEditFieldLabel_2 = uilabel(app.ViewvolumesTab);
            app.VoxelsizeEditFieldLabel_2.HorizontalAlignment = 'right';
            app.VoxelsizeEditFieldLabel_2.Position = [9 69 60 22];
            app.VoxelsizeEditFieldLabel_2.Text = 'Voxel size';

            % Create View_VoxelsizeEditField
            app.View_VoxelsizeEditField = uieditfield(app.ViewvolumesTab, 'numeric');
            app.View_VoxelsizeEditField.Limits = [1e-05 Inf];
            app.View_VoxelsizeEditField.Position = [84 69 63 22];
            app.View_VoxelsizeEditField.Value = 1;

            % Create UnitEditFieldLabel_2
            app.UnitEditFieldLabel_2 = uilabel(app.ViewvolumesTab);
            app.UnitEditFieldLabel_2.HorizontalAlignment = 'right';
            app.UnitEditFieldLabel_2.Position = [162 69 27 22];
            app.UnitEditFieldLabel_2.Text = 'Unit';

            % Create View_UnitEditField
            app.View_UnitEditField = uieditfield(app.ViewvolumesTab, 'text');
            app.View_UnitEditField.Position = [209 69 63 22];
            app.View_UnitEditField.Value = 'um';

            % Create ChecksegmentationTab
            app.ChecksegmentationTab = uitab(app.TabGroup);
            app.ChecksegmentationTab.Title = 'Check segmentation';

            % Create Checkseg_title
            app.Checkseg_title = uilabel(app.ChecksegmentationTab);
            app.Checkseg_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Checkseg_title.HorizontalAlignment = 'center';
            app.Checkseg_title.FontWeight = 'bold';
            app.Checkseg_title.Position = [9 559 657 22];
            app.Checkseg_title.Text = 'Select grey level and segemented volume files to compare them';

            % Create Checkseg_Importtif_grey_Button
            app.Checkseg_Importtif_grey_Button = uibutton(app.ChecksegmentationTab, 'push');
            app.Checkseg_Importtif_grey_Button.ButtonPushedFcn = createCallbackFcn(app, @Checkseg_Importtif_grey_ButtonPushed, true);
            app.Checkseg_Importtif_grey_Button.BackgroundColor = [0.8 0.8 0.8];
            app.Checkseg_Importtif_grey_Button.Position = [9 479 200 33];
            app.Checkseg_Importtif_grey_Button.Text = 'Click to import grey-level image';

            % Create Checkseg_Importtif_segmented_Button
            app.Checkseg_Importtif_segmented_Button = uibutton(app.ChecksegmentationTab, 'push');
            app.Checkseg_Importtif_segmented_Button.ButtonPushedFcn = createCallbackFcn(app, @Checkseg_Importtif_segmented_ButtonPushed, true);
            app.Checkseg_Importtif_segmented_Button.BackgroundColor = [0.8 0.8 0.8];
            app.Checkseg_Importtif_segmented_Button.Position = [9 434 200 33];
            app.Checkseg_Importtif_segmented_Button.Text = 'Click to import segmented image';

            % Create Checkseg_instructions
            app.Checkseg_instructions = uilabel(app.ChecksegmentationTab);
            app.Checkseg_instructions.Position = [9 527 657 22];
            app.Checkseg_instructions.Text = 'Select both files, enter a voxel size and then create figure. Volumes must share same dimension.';

            % Create Checkseg_createfigure_button
            app.Checkseg_createfigure_button = uibutton(app.ChecksegmentationTab, 'push');
            app.Checkseg_createfigure_button.ButtonPushedFcn = createCallbackFcn(app, @Checkseg_createfigure_buttonButtonPushed, true);
            app.Checkseg_createfigure_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Checkseg_createfigure_button.FontSize = 14;
            app.Checkseg_createfigure_button.FontWeight = 'bold';
            app.Checkseg_createfigure_button.FontColor = [1 1 1];
            app.Checkseg_createfigure_button.Enable = 'off';
            app.Checkseg_createfigure_button.Position = [9 352 200 33];
            app.Checkseg_createfigure_button.Text = 'Create figure';

            % Create Checkseg_Greyvolumename
            app.Checkseg_Greyvolumename = uilabel(app.ChecksegmentationTab);
            app.Checkseg_Greyvolumename.Position = [219 485 447 22];
            app.Checkseg_Greyvolumename.Text = 'No volume loaded';

            % Create Checkseg_Segmentedvolumename
            app.Checkseg_Segmentedvolumename = uilabel(app.ChecksegmentationTab);
            app.Checkseg_Segmentedvolumename.Position = [219 440 447 22];
            app.Checkseg_Segmentedvolumename.Text = 'No volume loaded';

            % Create VoxelsizeEditFieldLabel
            app.VoxelsizeEditFieldLabel = uilabel(app.ChecksegmentationTab);
            app.VoxelsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelsizeEditFieldLabel.Position = [9 400 60 22];
            app.VoxelsizeEditFieldLabel.Text = 'Voxel size';

            % Create Checkseg_VoxelsizeEditField
            app.Checkseg_VoxelsizeEditField = uieditfield(app.ChecksegmentationTab, 'numeric');
            app.Checkseg_VoxelsizeEditField.Limits = [1e-05 Inf];
            app.Checkseg_VoxelsizeEditField.Position = [84 400 63 22];
            app.Checkseg_VoxelsizeEditField.Value = 1;

            % Create UnitEditFieldLabel
            app.UnitEditFieldLabel = uilabel(app.ChecksegmentationTab);
            app.UnitEditFieldLabel.HorizontalAlignment = 'right';
            app.UnitEditFieldLabel.Position = [161 400 27 22];
            app.UnitEditFieldLabel.Text = 'Unit';

            % Create Checkseg_UnitEditField
            app.Checkseg_UnitEditField = uieditfield(app.ChecksegmentationTab, 'text');
            app.Checkseg_UnitEditField.Position = [208 400 63 22];
            app.Checkseg_UnitEditField.Value = 'um';

            % Create ViewcharacterizationresultsTab
            app.ViewcharacterizationresultsTab = uitab(app.TabGroup);
            app.ViewcharacterizationresultsTab.Title = 'View characterization results';

            % Create Res_title
            app.Res_title = uilabel(app.ViewcharacterizationresultsTab);
            app.Res_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Res_title.HorizontalAlignment = 'center';
            app.Res_title.FontWeight = 'bold';
            app.Res_title.Position = [9 559 657 22];
            app.Res_title.Text = 'Select one folder. You need to run the microstructure characterization module first before';

            % Create Res_Selectfolder_Button
            app.Res_Selectfolder_Button = uibutton(app.ViewcharacterizationresultsTab, 'push');
            app.Res_Selectfolder_Button.ButtonPushedFcn = createCallbackFcn(app, @Res_Selectfolder_ButtonPushed, true);
            app.Res_Selectfolder_Button.BackgroundColor = [0.8 0.8 0.8];
            app.Res_Selectfolder_Button.Position = [9 479 200 33];
            app.Res_Selectfolder_Button.Text = 'Click to select a  folder';

            % Create Res_instructions
            app.Res_instructions = uilabel(app.ViewcharacterizationresultsTab);
            app.Res_instructions.Position = [9 521 657 28];
            app.Res_instructions.Text = {'Select the main folder generated by the microstructure characteriation module, then choose the property to plot and for'; 'which phase. Enter a voxel size and then create the figure.'};

            % Create Res_Labelloaded
            app.Res_Labelloaded = uilabel(app.ViewcharacterizationresultsTab);
            app.Res_Labelloaded.Position = [219 485 447 22];
            app.Res_Labelloaded.Text = 'No folder selected';

            % Create Res_createfigure_button
            app.Res_createfigure_button = uibutton(app.ViewcharacterizationresultsTab, 'push');
            app.Res_createfigure_button.ButtonPushedFcn = createCallbackFcn(app, @Res_createfigure_buttonButtonPushed, true);
            app.Res_createfigure_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Res_createfigure_button.FontSize = 14;
            app.Res_createfigure_button.FontWeight = 'bold';
            app.Res_createfigure_button.FontColor = [1 1 1];
            app.Res_createfigure_button.Enable = 'off';
            app.Res_createfigure_button.Position = [9 24 164 33];
            app.Res_createfigure_button.Text = '2D slices visualization';

            % Create Res_3Dfigure_button
            app.Res_3Dfigure_button = uibutton(app.ViewcharacterizationresultsTab, 'push');
            app.Res_3Dfigure_button.ButtonPushedFcn = createCallbackFcn(app, @Res_3Dfigure_buttonButtonPushed, true);
            app.Res_3Dfigure_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Res_3Dfigure_button.FontSize = 14;
            app.Res_3Dfigure_button.FontWeight = 'bold';
            app.Res_3Dfigure_button.FontColor = [1 1 1];
            app.Res_3Dfigure_button.Enable = 'off';
            app.Res_3Dfigure_button.Position = [188 24 133 33];
            app.Res_3Dfigure_button.Text = '3D visualization';

            % Create VoxelsizeEditFieldLabel_3
            app.VoxelsizeEditFieldLabel_3 = uilabel(app.ViewcharacterizationresultsTab);
            app.VoxelsizeEditFieldLabel_3.HorizontalAlignment = 'right';
            app.VoxelsizeEditFieldLabel_3.Position = [9 69 60 22];
            app.VoxelsizeEditFieldLabel_3.Text = 'Voxel size';

            % Create Res_VoxelsizeEditField
            app.Res_VoxelsizeEditField = uieditfield(app.ViewcharacterizationresultsTab, 'numeric');
            app.Res_VoxelsizeEditField.Limits = [1e-05 Inf];
            app.Res_VoxelsizeEditField.Position = [84 69 63 22];
            app.Res_VoxelsizeEditField.Value = 1;

            % Create UnitEditFieldLabel_3
            app.UnitEditFieldLabel_3 = uilabel(app.ViewcharacterizationresultsTab);
            app.UnitEditFieldLabel_3.HorizontalAlignment = 'right';
            app.UnitEditFieldLabel_3.Position = [162 69 27 22];
            app.UnitEditFieldLabel_3.Text = 'Unit';

            % Create Res_UnitEditField
            app.Res_UnitEditField = uieditfield(app.ViewcharacterizationresultsTab, 'text');
            app.Res_UnitEditField.Position = [209 69 63 22];
            app.Res_UnitEditField.Value = 'um';

            % Create PropertytoplotDropDownLabel
            app.PropertytoplotDropDownLabel = uilabel(app.ViewcharacterizationresultsTab);
            app.PropertytoplotDropDownLabel.HorizontalAlignment = 'right';
            app.PropertytoplotDropDownLabel.Enable = 'off';
            app.PropertytoplotDropDownLabel.Position = [9 440 87 22];
            app.PropertytoplotDropDownLabel.Text = 'Property to plot';

            % Create PropertytoplotDropDown
            app.PropertytoplotDropDown = uidropdown(app.ViewcharacterizationresultsTab);
            app.PropertytoplotDropDown.Items = {};
            app.PropertytoplotDropDown.ValueChangedFcn = createCallbackFcn(app, @PropertytoplotDropDownValueChanged, true);
            app.PropertytoplotDropDown.Enable = 'off';
            app.PropertytoplotDropDown.Position = [111 440 236 22];
            app.PropertytoplotDropDown.Value = {};

            % Create ForthephaseDropDownLabel
            app.ForthephaseDropDownLabel = uilabel(app.ViewcharacterizationresultsTab);
            app.ForthephaseDropDownLabel.HorizontalAlignment = 'right';
            app.ForthephaseDropDownLabel.Enable = 'off';
            app.ForthephaseDropDownLabel.Position = [9 413 80 22];
            app.ForthephaseDropDownLabel.Text = 'For the phase';

            % Create ForthephaseDropDown
            app.ForthephaseDropDown = uidropdown(app.ViewcharacterizationresultsTab);
            app.ForthephaseDropDown.Items = {};
            app.ForthephaseDropDown.ValueChangedFcn = createCallbackFcn(app, @ForthephaseDropDownValueChanged, true);
            app.ForthephaseDropDown.Enable = 'off';
            app.ForthephaseDropDown.Position = [111 413 236 22];
            app.ForthephaseDropDown.Value = {};

            % Create Res_ErrormessageLabel
            app.Res_ErrormessageLabel = uilabel(app.ViewcharacterizationresultsTab);
            app.Res_ErrormessageLabel.Position = [9 379 657 22];
            app.Res_ErrormessageLabel.Text = 'Error message';

            % Create PropertyunitLabel
            app.PropertyunitLabel = uilabel(app.ViewcharacterizationresultsTab);
            app.PropertyunitLabel.HorizontalAlignment = 'right';
            app.PropertyunitLabel.Position = [287 69 74 22];
            app.PropertyunitLabel.Text = 'Property unit';

            % Create Res_ResUnitEditField
            app.Res_ResUnitEditField = uieditfield(app.ViewcharacterizationresultsTab, 'text');
            app.Res_ResUnitEditField.Position = [381 69 63 22];

            % Show the figure after all components are created
            app.Visualizationmodule_UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = microstructure_visualization_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.Visualizationmodule_UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.Visualizationmodule_UIFigure)
        end
    end
end