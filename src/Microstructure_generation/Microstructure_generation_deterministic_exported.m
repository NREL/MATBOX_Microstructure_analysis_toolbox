classdef Microstructure_generation_deterministic_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        DeterministicgenerationUIFigure  matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        InstructionsTab                matlab.ui.container.Tab
        Instructions_title             matlab.ui.control.Label
        Save_Label                     matlab.ui.control.Label
        ClicktoselectsavefolderButton  matlab.ui.control.Button
        Save_folder_text               matlab.ui.control.Label
        Label                          matlab.ui.control.Label
        PeriodicspheresTab             matlab.ui.container.Tab
        Spheres_title                  matlab.ui.control.Label
        Sphere_GenerateButton          matlab.ui.control.Button
        Sphere_VisualizeButton         matlab.ui.control.Button
        Sphere_SaveButton              matlab.ui.control.Button
        StatutTextArea_2Label_2        matlab.ui.control.Label
        Spheres_StatutTextArea         matlab.ui.control.TextArea
        Coil_Label_2                   matlab.ui.control.Label
        Spheres_parameter_UITable      matlab.ui.control.Table
        Spheres_TextArea               matlab.ui.control.TextArea
        Image2                         matlab.ui.control.Image
        ExampleLabel_2                 matlab.ui.control.Label
        LayersTab                      matlab.ui.container.Tab
        Layer_title                    matlab.ui.control.Label
        CoilserpentineTab              matlab.ui.container.Tab
        Coil_title                     matlab.ui.control.Label
        Coil_Label                     matlab.ui.control.Label
        Coil_GenerateButton            matlab.ui.control.Button
        Coil_VisualizeButton           matlab.ui.control.Button
        Coil_SaveButton                matlab.ui.control.Button
        StatutTextArea_2Label          matlab.ui.control.Label
        Coil_StatutTextArea            matlab.ui.control.TextArea
        Coil_parameter_UITable         matlab.ui.control.Table
        TextArea                       matlab.ui.control.TextArea
        Image                          matlab.ui.control.Image
        ExampleLabel                   matlab.ui.control.Label
        Coil_CheckBox                  matlab.ui.control.CheckBox
        FractalTab                     matlab.ui.container.Tab
        Fractal_title                  matlab.ui.control.Label
        WhyfractalAreyounotcuriousabouttheirtortuosityfactorLabel  matlab.ui.control.Label
        SelectfractalDropDownLabel     matlab.ui.control.Label
        SelectfractalDropDown          matlab.ui.control.DropDown
        Fractal_parameter_UITable      matlab.ui.control.Table
        Fractal_GenerateButton         matlab.ui.control.Button
        Fractal_VisualizeButton        matlab.ui.control.Button
        Fractal_SaveButton             matlab.ui.control.Button
        StatutTextAreaLabel            matlab.ui.control.Label
        StatutTextArea                 matlab.ui.control.TextArea
        Fractal_HelpTextArea           matlab.ui.control.TextArea
        Rectangular_ChannelsTab        matlab.ui.container.Tab
        Rchannel_title                 matlab.ui.control.Label
        Channel_instructionsLabel      matlab.ui.control.Label
        Channel_GenerateButton         matlab.ui.control.Button
        Channel_VisualizeButton        matlab.ui.control.Button
        Channel_SaveButton             matlab.ui.control.Button
        StatutTextArea_3Label          matlab.ui.control.Label
        Channel_StatutTextArea         matlab.ui.control.TextArea
        PatternDropDownLabel           matlab.ui.control.Label
        Channel_PatternDropDown        matlab.ui.control.DropDown
        Channel_UITable                matlab.ui.control.Table
        Channel_rectangular1D_Image    matlab.ui.control.Image
        Channel_apply_text             matlab.ui.control.Label
        ApplypatternonLabel            matlab.ui.control.Label
        Channel_applyPatternDropDown   matlab.ui.control.DropDown
        ClicktoloadmicrostructureButton  matlab.ui.control.Button
        ApplyChannel_UITable           matlab.ui.control.Table
    end

    
    properties (Access = private)
        % Saving options
        mainsavefoldersave = []; savefolder_is_defined = false;
        % Geometry
        loaded_microstructure = [];
        microstructure = [];
        name_microstructure = []; name_microstructure_loaded=[];
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Coil/serpentine
            parameters = {'Coil total diameter';'Coil thickness';'Coil length';'Rotation speed';'No-coil region normalized thickness';'No-coil region normalized length'};
            values = {'40';'20';'100';'1';'0.5';'0.2'};
            range = {'>1, integer';'>1, integer';'>1, integer';'>=0, float';'>=0, float';'>=0, float'};
            app.Coil_parameter_UITable.Data = [parameters values range];
            % Spheres
            parameters = {'Number of particles along through plane';'Number of particles along in plane';'In plane periodicity';'Particle diameter';'Normalized distance between center';'No-particle region normalized through-plane length';'No-particle region normalized in-plane length'};
            values = {'6';'3';'1';'50';'0.9';'0.2';'0.2';};
            range = {'>1, integer';'>1, integer';'1 (true) or 0 (false)';'>1, integer';'[0,1], float';'>=0, float';'>=0, float'};
            app.Spheres_parameter_UITable.Data = [parameters values range];   
            % Channels
            app.Channel_PatternDropDownValueChanged
            app.Channel_applyPatternDropDownValueChanged
        end

        % Button pushed function: ClicktoselectsavefolderButton
        function ClicktoselectsavefolderButtonPushed(app, event)
            str_dialogbox = 'Select location where the save folder will be created';
            resultfolder_location = uigetdir(matlabroot,str_dialogbox); % Open dialog box
            if resultfolder_location==0
                % User clicked cancel button or closed the dialog box
                set(app.Save_folder_text,'Text','Save folder location: NOT DEFINED','FontColor','r');
                app.savefolder_is_defined = false;
            else
                % Update savefolder
                if ispc
                    app.mainsavefoldersave = [resultfolder_location '\'];
                else
                    app.mainsavefoldersave = [resultfolder_location '/'];
                end
                app.savefolder_is_defined = true;
                set(app.Save_folder_text,'Text',['Save folder location: ' app.mainsavefoldersave],'FontColor','k');
            end
        end

        % Button pushed function: Fractal_GenerateButton
        function Fractal_GenerateButtonPushed(app, event)
            fractal_choice = app.SelectfractalDropDown.Value;
            app.Fractal_VisualizeButton.Enable = 'off';app.Fractal_SaveButton.Enable = 'off';
            if strcmp(fractal_choice,'Menger sponge')
                tmp = app.Fractal_parameter_UITable.Data;
                number_of_iteration = str2num(cell2mat(tmp(1,2)));
                Cube_subdivision = str2num(cell2mat(tmp(2,2)));
                app.StatutTextArea.Value = 'Generating Menger sponge...'; pause(0.1);
                [app.microstructure] = generate_menger_sponge(number_of_iteration,Cube_subdivision);
                app.StatutTextArea.Value = ['Done! Volume dimension is ' num2str(size(app.microstructure))];
                app.name_microstructure = ['MengerSponge_Iter' num2str(number_of_iteration) '_Subdivision' num2str(Cube_subdivision)];
            end
            app.Fractal_VisualizeButton.Enable = 'on';
            if app.savefolder_is_defined
                app.Fractal_SaveButton.Enable = 'on';
            end            
        end

        % Button pushed function: Coil_VisualizeButton, 
        % Fractal_VisualizeButton, Sphere_VisualizeButton
        function Fractal_VisualizeButtonPushed(app, event)
            Microstructure_basic_visualization_interface(app.microstructure);                    
        end

        % Button pushed function: Fractal_SaveButton
        function Fractal_SaveButtonPushed(app, event)
            if app.savefolder_is_defined
                app.StatutTextArea.Value = 'Saving...'; pause(0.1);
                function_save_tif(app.microstructure,[app.mainsavefoldersave app.name_microstructure '.tif']);
                app.StatutTextArea.Value = 'Done!';
            end
        end

        % Value changed function: SelectfractalDropDown
        function SelectfractalDropDownValueChanged(app, event)
            fractal_choice = app.SelectfractalDropDown.Value;
            if strcmp(fractal_choice,'Menger sponge')
                parameters = {'Number of iteration';'Cube subdivision'};
                values = {'3';'2'};
                range = {'>1, integer';'>1, integer'};
                app.Fractal_parameter_UITable.Data = [parameters values range];
                app.Fractal_HelpTextArea.Value = 'Number of voxel along each axe is: (3^Number of iteration)*Cube subdivision';
                app.Fractal_HelpTextArea.Visible = 'on';
                app.Fractal_GenerateButton.Enable = 'on'; app.Fractal_VisualizeButton.Enable = 'off';app.Fractal_SaveButton.Enable = 'off';
            end
            
        end

        % Button pushed function: Coil_SaveButton
        function Coil_SaveButtonPushed(app, event)
            if app.savefolder_is_defined
                app.Coil_StatutTextArea.Value = 'Saving...'; pause(0.1);
                function_save_tif(app.microstructure,[app.mainsavefoldersave app.name_microstructure '.tif']);
                app.Coil_StatutTextArea.Value = 'Done!';
            end
        end

        % Button pushed function: Coil_GenerateButton
        function Coil_GenerateButtonPushed(app, event)
            app.Coil_VisualizeButton.Enable = 'off';app.Coil_SaveButton.Enable = 'off';
                tmp = app.Coil_parameter_UITable.Data;
                coil_diameter = str2num(cell2mat(tmp(1,2)));
                coil_thickness = str2num(cell2mat(tmp(2,2)));
                coil_length = str2num(cell2mat(tmp(3,2)));
                rotation_speed = str2num(cell2mat(tmp(4,2)));
                additional_distance_inplane = str2num(cell2mat(tmp(5,2)));
                additional_distance_throughplane = str2num(cell2mat(tmp(6,2)));
                coil_is_background = app.Coil_CheckBox.Value;
                app.Coil_StatutTextArea.Value = 'Generating serptentine...'; pause(0.1);
                [app.microstructure] = generate_coil(coil_diameter,coil_thickness,coil_length,rotation_speed,additional_distance_throughplane,additional_distance_inplane, coil_is_background);
                app.Coil_StatutTextArea.Value = ['Done! Volume dimension is ' num2str(size(app.microstructure))];
                app.name_microstructure = ['Serptentine_Diameter' num2str(coil_diameter) '_Thickness' num2str(coil_thickness) '_Length' num2str(coil_length) '_Rotation' num2str(rotation_speed)];
            app.Coil_VisualizeButton.Enable = 'on';
            if app.savefolder_is_defined
                app.Coil_SaveButton.Enable = 'on';
            end
        end

        % Button pushed function: Sphere_SaveButton
        function Sphere_SaveButtonPushed(app, event)
            app.Spheres_StatutTextArea.Value = 'Saving...'; pause(0.1);
            function_save_tif(app.microstructure,[app.mainsavefoldersave app.name_microstructure '.tif']);
            app.Spheres_StatutTextArea.Value = 'Done!';
        end

        % Button pushed function: Sphere_GenerateButton
        function Sphere_GenerateButtonPushed(app, event)
            app.Sphere_VisualizeButton.Enable = 'off';app.Sphere_SaveButton.Enable = 'off';
            tmp = app.Spheres_parameter_UITable.Data;
            p.number_particles_through_plane = str2num(cell2mat(tmp(1,2)));
            p.number_particles_in_plane = str2num(cell2mat(tmp(2,2)));
            p.cut_particle_inplane_extremities = str2num(cell2mat(tmp(3,2)));           
            p.particle_diameter = str2num(cell2mat(tmp(4,2)));   
            p.distance_between_center = str2num(cell2mat(tmp(5,2)));   
            p.additional_distance_throughplane = str2num(cell2mat(tmp(6,2)));   
            p.additional_distance_inplane = str2num(cell2mat(tmp(7,2)));   
            app.Spheres_StatutTextArea.Value = 'Generating spheres...'; pause(0.1);
            [app.microstructure] = generate_spheres(p);
            app.Spheres_StatutTextArea.Value = ['Done! Volume dimension is ' num2str(size(app.microstructure))];
            app.name_microstructure = ['Spheres_nTP' num2str(p.number_particles_through_plane) '_nIP' num2str(p.number_particles_in_plane) '_Periodic' num2str(p.particle_diameter) '_Diameter' num2str(p.particle_diameter)];
            app.Sphere_VisualizeButton.Enable = 'on';
            if app.savefolder_is_defined
                app.Sphere_SaveButton.Enable = 'on';
            end
        end

        % Value changed function: Channel_PatternDropDown
        function Channel_PatternDropDownValueChanged(app, event)
            pattern_choice = app.Channel_PatternDropDown.Value;
            if strcmp(pattern_choice,'Rectangular channels along one direction')
                parameters = {'Porous matrix half-width w1/2';'Channel top half width w2/2';'Channel bottom normalized width w';'Channel normalized thikness t';'Domain''s edge crops porous matrix (1) or channel (0)';'Crop domain''s edge for periodicity';'Channel depth along axis';'Channel in-plane axis'};
                values = {'20';'5';'0.5';'0.5';'1';'1';'3';'1'};
                range = {'>1, integer';'>1, integer';'[0,1], float';'[0,1], float';'1 or 0';'1 (true) or 0 (false)';'1, 2 or 3';'1, 2 or 3, but different from previous parameter'};
                app.Channel_UITable.Data = [parameters values range];
                app.Channel_rectangular1D_Image.Visible = 'on';
            end
           
            
                
            
        end

        % Button pushed function: ClicktoloadmicrostructureButton
        function ClicktoloadmicrostructureButtonPushed(app, event)
            
            [FileName,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)'},'Select volume');
            if FileName==0
                % User clicked cancel button or closed the dialog box
                app.loaded_microstructure = []; app.microstructure = [];
            else
                app.Channel_StatutTextArea.Value = 'Loading... please wait'; pause(0.01);
                [app.loaded_microstructure,outcome] = function_load_tif([PathName FileName]);
                if outcome % Success
                    app.Channel_StatutTextArea.Value = 'Done!';
                    app.Channel_GenerateButton.Enable = 'on'; app.Channel_VisualizeButton.Enable = 'on';
                else
                    app.Channel_StatutTextArea.Value = 'File could not be loaded!';
                    app.loaded_microstructure = []; app.microstructure = [];
                    app.Channel_GenerateButton.Enable = 'off'; app.Channel_VisualizeButton.Enable = 'off';
                end
            end            
            
            
           
            
            
        end

        % Value changed function: Channel_applyPatternDropDown
        function Channel_applyPatternDropDownValueChanged(app, event)
            apply_choice = app.Channel_applyPatternDropDown.Value;
            if strcmp(apply_choice,'Heterogenous microstructure (cut through particles)')
                app.ClicktoloadmicrostructureButton.Visible = 'on'; app.ClicktoloadmicrostructureButton.Enable = 'on';
                app.ApplyChannel_UITable.Visible = 'off'; app.ApplyChannel_UITable.Enable = 'off';
                app.Channel_GenerateButton.Enable = 'off'; app.Channel_VisualizeButton.Enable = 'off';
            elseif strcmp(apply_choice,'Heterogenous microstructure (remove particles)')
                app.ClicktoloadmicrostructureButton.Visible = 'on'; app.ClicktoloadmicrostructureButton.Enable = 'on';
                app.ApplyChannel_UITable.Visible = 'on'; app.ApplyChannel_UITable.Enable = 'on';
                parameters = {'Percentage of particles that intersect with macropore channels interface to remove'};
                values = {'50'};
                range = {'[0,100], float'};
                app.ApplyChannel_UITable.Data = [parameters values range];                
                app.Channel_GenerateButton.Enable = 'off'; app.Channel_VisualizeButton.Enable = 'off';
            elseif strcmp(apply_choice,'Homogenous medium') 
                app.ClicktoloadmicrostructureButton.Visible = 'off'; app.ClicktoloadmicrostructureButton.Enable = 'off';
                parameters = {'Number of voxel along direction 1';'along direction 2';'along direction 3'};
                values = {'100';'100';'50'};
                range = {'>1, integer';'>1, integer';'>1, integer'};
                app.ApplyChannel_UITable.Data = [parameters values range];                   
                app.ApplyChannel_UITable.Visible = 'on'; app.ApplyChannel_UITable.Enable = 'on';
                app.Channel_GenerateButton.Enable = 'on'; app.Channel_VisualizeButton.Enable = 'off';
            else
                app.ApplyChannel_UITable.Visible = 'off'; app.ApplyChannel_UITable.Enable = 'off';
                app.ClicktoloadmicrostructureButton.Visible = 'off'; app.ClicktoloadmicrostructureButton.Enable = 'off';
                app.Channel_GenerateButton.Enable = 'off'; app.Channel_VisualizeButton.Enable = 'off';
            end
            
        end

        % Button pushed function: Channel_VisualizeButton
        function Channel_VisualizeButtonPushed(app, event)
            if isempty(app.microstructure)
                Microstructure_basic_visualization_interface(app.loaded_microstructure); 
            else
                Microstructure_comparison_visualization_interface(app.loaded_microstructure, app.microstructure);
            end            
        end

        % Button pushed function: Channel_GenerateButton
        function Channel_GenerateButtonPushed(app, event)
            app.Channel_StatutTextArea.Value = 'Generating channels...'; pause(0.1);
            p.apply = char(app.Channel_applyPatternDropDown.Value);
            if strcmp(p.apply,'Heterogenous microstructure (remove particles)')
                p.percentremove = str2num(cell2mat(app.ApplyChannel_UITable.Data(1,2)));
            elseif strcmp(p.apply,'Homogenous medium')
                n1 = str2num(cell2mat(app.ApplyChannel_UITable.Data(1,2)));
                n2 = str2num(cell2mat(app.ApplyChannel_UITable.Data(2,2)));
                n3 = str2num(cell2mat(app.ApplyChannel_UITable.Data(3,2)));
                app.loaded_microstructure = ones(n1,n2,n3);
            end
            
            pattern_choice = app.Channel_PatternDropDown.Value;
            tmp = app.Channel_UITable.Data(:,2);
            if strcmp(pattern_choice,'Rectangular channels along one direction') 
                p.half_w1 = str2num(cell2mat(tmp(1)));
                p.half_w2 = str2num(cell2mat(tmp(2)));
                p.w = str2num(cell2mat(tmp(3)));
                p.t = str2num(cell2mat(tmp(4)));
                p.crop_porous1_channel0 = str2num(cell2mat(tmp(5)));
                p.crop_periodicity = str2num(cell2mat(tmp(6)));
                p.thickness_axis = str2num(cell2mat(tmp(7)));
                p.inplane_axis = str2num(cell2mat(tmp(8)));
                [app.loaded_microstructure, app.microstructure] = generate_rectangular_1D_channels(app.loaded_microstructure, p);
                app.name_microstructure = ['Channel_R1D_w1_' num2str(2*p.half_w1) '_w2_' num2str(2*p.half_w2) '_w_' num2str(p.w) '_t_' num2str(p.t) '_Crop_' num2str(p.crop_porous1_channel0) '_Periodic_' num2str(p.crop_periodicity)];
                app.name_microstructure_loaded = 'Channel_R1D_baseline';
            end
            app.Channel_VisualizeButton.Enable = 'on'; app.Channel_SaveButton.Enable = 'on';
            app.Channel_StatutTextArea.Value = ['Done! Volume dimension is ' num2str(size(app.microstructure))]; 
        end

        % Button pushed function: Channel_SaveButton
        function Channel_SaveButtonPushed(app, event)
            app.Channel_StatutTextArea.Value = 'Saving...'; pause(0.1);
            function_save_tif(app.loaded_microstructure,[app.mainsavefoldersave app.name_microstructure_loaded '.tif']);
            function_save_tif(app.microstructure,[app.mainsavefoldersave app.name_microstructure '.tif']);
            app.Channel_StatutTextArea.Value = 'Done!';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create DeterministicgenerationUIFigure and hide until all components are created
            app.DeterministicgenerationUIFigure = uifigure('Visible', 'off');
            app.DeterministicgenerationUIFigure.Position = [100 100 850 591];
            app.DeterministicgenerationUIFigure.Name = 'Deterministic generation';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.DeterministicgenerationUIFigure);
            app.TabGroup.TabLocation = 'left';
            app.TabGroup.Position = [1 1 850 591];

            % Create InstructionsTab
            app.InstructionsTab = uitab(app.TabGroup);
            app.InstructionsTab.Title = 'Instructions';

            % Create Instructions_title
            app.Instructions_title = uilabel(app.InstructionsTab);
            app.Instructions_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_title.HorizontalAlignment = 'center';
            app.Instructions_title.FontWeight = 'bold';
            app.Instructions_title.Position = [8 559 709 22];
            app.Instructions_title.Text = 'Instructions';

            % Create Save_Label
            app.Save_Label = uilabel(app.InstructionsTab);
            app.Save_Label.Position = [8 505 709 42];
            app.Save_Label.Text = {'Choose a save folder and then go to the tab that correspond to the type of microstructure you want to generate.'; 'Blue tabs create a bi-phase (background and generated solid) geometry from scratch.'; 'Red tabs add a phase on top of an existing geometry.'};

            % Create ClicktoselectsavefolderButton
            app.ClicktoselectsavefolderButton = uibutton(app.InstructionsTab, 'push');
            app.ClicktoselectsavefolderButton.ButtonPushedFcn = createCallbackFcn(app, @ClicktoselectsavefolderButtonPushed, true);
            app.ClicktoselectsavefolderButton.BackgroundColor = [0.8 0.8 0.8];
            app.ClicktoselectsavefolderButton.Position = [8 472 166 26];
            app.ClicktoselectsavefolderButton.Text = 'Click to select save folder';

            % Create Save_folder_text
            app.Save_folder_text = uilabel(app.InstructionsTab);
            app.Save_folder_text.FontColor = [1 0 0];
            app.Save_folder_text.Position = [184 474 533 22];
            app.Save_folder_text.Text = 'Save folder location: NOT DEFINED';

            % Create Label
            app.Label = uilabel(app.InstructionsTab);
            app.Label.FontColor = [1 0 0];
            app.Label.Position = [8 45 417 22];
            app.Label.Text = 'This submodule is not yet completed (all functionnality are not implemented)';

            % Create PeriodicspheresTab
            app.PeriodicspheresTab = uitab(app.TabGroup);
            app.PeriodicspheresTab.Title = 'Periodic spheres';
            app.PeriodicspheresTab.ForegroundColor = [0 0 1];

            % Create Spheres_title
            app.Spheres_title = uilabel(app.PeriodicspheresTab);
            app.Spheres_title.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Spheres_title.HorizontalAlignment = 'center';
            app.Spheres_title.FontWeight = 'bold';
            app.Spheres_title.Position = [8 559 709 22];
            app.Spheres_title.Text = 'Spheres-based geometries';

            % Create Sphere_GenerateButton
            app.Sphere_GenerateButton = uibutton(app.PeriodicspheresTab, 'push');
            app.Sphere_GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @Sphere_GenerateButtonPushed, true);
            app.Sphere_GenerateButton.BackgroundColor = [0.8 0.8 0.8];
            app.Sphere_GenerateButton.Position = [8 11 100 40];
            app.Sphere_GenerateButton.Text = 'Generate';

            % Create Sphere_VisualizeButton
            app.Sphere_VisualizeButton = uibutton(app.PeriodicspheresTab, 'push');
            app.Sphere_VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @Fractal_VisualizeButtonPushed, true);
            app.Sphere_VisualizeButton.BackgroundColor = [0.8 0.8 0.8];
            app.Sphere_VisualizeButton.Enable = 'off';
            app.Sphere_VisualizeButton.Position = [128 11 100 40];
            app.Sphere_VisualizeButton.Text = 'Visualize';

            % Create Sphere_SaveButton
            app.Sphere_SaveButton = uibutton(app.PeriodicspheresTab, 'push');
            app.Sphere_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Sphere_SaveButtonPushed, true);
            app.Sphere_SaveButton.BackgroundColor = [0.8 0.8 0.8];
            app.Sphere_SaveButton.Enable = 'off';
            app.Sphere_SaveButton.Position = [248 11 100 40];
            app.Sphere_SaveButton.Text = 'Save';

            % Create StatutTextArea_2Label_2
            app.StatutTextArea_2Label_2 = uilabel(app.PeriodicspheresTab);
            app.StatutTextArea_2Label_2.HorizontalAlignment = 'right';
            app.StatutTextArea_2Label_2.Position = [366 27 37 22];
            app.StatutTextArea_2Label_2.Text = 'Statut';

            % Create Spheres_StatutTextArea
            app.Spheres_StatutTextArea = uitextarea(app.PeriodicspheresTab);
            app.Spheres_StatutTextArea.Position = [418 11 299 40];
            app.Spheres_StatutTextArea.Value = {'Idle'};

            % Create Coil_Label_2
            app.Coil_Label_2 = uilabel(app.PeriodicspheresTab);
            app.Coil_Label_2.FontAngle = 'italic';
            app.Coil_Label_2.Position = [8 524 709 28];
            app.Coil_Label_2.Text = {'Sphere-based geometries are relevant for comparing electrochemical lithium ion battery micro-scale models with their macro-scale'; 'counterpart as the latters are typically relying on a spherical assumptions.'};

            % Create Spheres_parameter_UITable
            app.Spheres_parameter_UITable = uitable(app.PeriodicspheresTab);
            app.Spheres_parameter_UITable.ColumnName = {'Parameter'; 'Value'; 'Range'};
            app.Spheres_parameter_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Spheres_parameter_UITable.RowName = {};
            app.Spheres_parameter_UITable.ColumnEditable = [false true false];
            app.Spheres_parameter_UITable.Position = [8 321 411 188];

            % Create Spheres_TextArea
            app.Spheres_TextArea = uitextarea(app.PeriodicspheresTab);
            app.Spheres_TextArea.Position = [434 419 247 90];
            app.Spheres_TextArea.Value = {'Particle diameter in number of voxels. Distance between particle center and no particle region distance are normalized with particle diameter.'; ''};

            % Create Image2
            app.Image2 = uiimage(app.PeriodicspheresTab);
            app.Image2.Position = [434 234 218 170];
            app.Image2.ImageSource = 'Periodic_spheres.png';

            % Create ExampleLabel_2
            app.ExampleLabel_2 = uilabel(app.PeriodicspheresTab);
            app.ExampleLabel_2.FontAngle = 'italic';
            app.ExampleLabel_2.Position = [592 368 52 22];
            app.ExampleLabel_2.Text = 'Example';

            % Create LayersTab
            app.LayersTab = uitab(app.TabGroup);
            app.LayersTab.Title = 'Layers';
            app.LayersTab.ForegroundColor = [0 0 1];

            % Create Layer_title
            app.Layer_title = uilabel(app.LayersTab);
            app.Layer_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Layer_title.HorizontalAlignment = 'center';
            app.Layer_title.FontWeight = 'bold';
            app.Layer_title.Position = [8 559 709 22];
            app.Layer_title.Text = 'Layers-based geometries';

            % Create CoilserpentineTab
            app.CoilserpentineTab = uitab(app.TabGroup);
            app.CoilserpentineTab.Title = 'Coil/serpentine';
            app.CoilserpentineTab.ForegroundColor = [0 0 1];

            % Create Coil_title
            app.Coil_title = uilabel(app.CoilserpentineTab);
            app.Coil_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Coil_title.HorizontalAlignment = 'center';
            app.Coil_title.FontWeight = 'bold';
            app.Coil_title.Position = [8 559 709 22];
            app.Coil_title.Text = 'Serpentine geometry';

            % Create Coil_Label
            app.Coil_Label = uilabel(app.CoilserpentineTab);
            app.Coil_Label.FontAngle = 'italic';
            app.Coil_Label.Position = [8 524 709 28];
            app.Coil_Label.Text = {'Serpentine geometies are very sinous - thus with high tortuosity. They are interesting to test depletion scenario in electrochemical model'; 'for electrolytem while having a small geometry then being relevant for debugging.'};

            % Create Coil_GenerateButton
            app.Coil_GenerateButton = uibutton(app.CoilserpentineTab, 'push');
            app.Coil_GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @Coil_GenerateButtonPushed, true);
            app.Coil_GenerateButton.BackgroundColor = [0.8 0.8 0.8];
            app.Coil_GenerateButton.Position = [8 11 100 40];
            app.Coil_GenerateButton.Text = 'Generate';

            % Create Coil_VisualizeButton
            app.Coil_VisualizeButton = uibutton(app.CoilserpentineTab, 'push');
            app.Coil_VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @Fractal_VisualizeButtonPushed, true);
            app.Coil_VisualizeButton.BackgroundColor = [0.8 0.8 0.8];
            app.Coil_VisualizeButton.Enable = 'off';
            app.Coil_VisualizeButton.Position = [128 11 100 40];
            app.Coil_VisualizeButton.Text = 'Visualize';

            % Create Coil_SaveButton
            app.Coil_SaveButton = uibutton(app.CoilserpentineTab, 'push');
            app.Coil_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Coil_SaveButtonPushed, true);
            app.Coil_SaveButton.BackgroundColor = [0.8 0.8 0.8];
            app.Coil_SaveButton.Enable = 'off';
            app.Coil_SaveButton.Position = [248 11 100 40];
            app.Coil_SaveButton.Text = 'Save';

            % Create StatutTextArea_2Label
            app.StatutTextArea_2Label = uilabel(app.CoilserpentineTab);
            app.StatutTextArea_2Label.HorizontalAlignment = 'right';
            app.StatutTextArea_2Label.Position = [366 27 37 22];
            app.StatutTextArea_2Label.Text = 'Statut';

            % Create Coil_StatutTextArea
            app.Coil_StatutTextArea = uitextarea(app.CoilserpentineTab);
            app.Coil_StatutTextArea.Position = [418 11 299 40];
            app.Coil_StatutTextArea.Value = {'Idle'};

            % Create Coil_parameter_UITable
            app.Coil_parameter_UITable = uitable(app.CoilserpentineTab);
            app.Coil_parameter_UITable.ColumnName = {'Parameter'; 'Value'; 'Range'};
            app.Coil_parameter_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Coil_parameter_UITable.RowName = {};
            app.Coil_parameter_UITable.ColumnEditable = [false true false];
            app.Coil_parameter_UITable.Position = [8 342 368 167];

            % Create TextArea
            app.TextArea = uitextarea(app.CoilserpentineTab);
            app.TextArea.Position = [391 409 290 100];
            app.TextArea.Value = {'Total diameter should be > 2 times the coil thickness (dimension in number of voxel). ''No coil region'' parameters are normalized with coil thickness and length. Coil angle is increased per voxel increment along coil length by rotation speed * 2 * pi/coil_diameter.'};

            % Create Image
            app.Image = uiimage(app.CoilserpentineTab);
            app.Image.Position = [392 220 218 170];
            app.Image.ImageSource = 'Coil.png';

            % Create ExampleLabel
            app.ExampleLabel = uilabel(app.CoilserpentineTab);
            app.ExampleLabel.FontAngle = 'italic';
            app.ExampleLabel.Position = [550 351 52 22];
            app.ExampleLabel.Text = 'Example';

            % Create Coil_CheckBox
            app.Coil_CheckBox = uicheckbox(app.CoilserpentineTab);
            app.Coil_CheckBox.Text = {'Assign coil to 0 and background to 1.'; 'Coil will be then the pore domain.'};
            app.Coil_CheckBox.Position = [8 302 223 28];

            % Create FractalTab
            app.FractalTab = uitab(app.TabGroup);
            app.FractalTab.Title = 'Fractal';
            app.FractalTab.ForegroundColor = [0 0 1];

            % Create Fractal_title
            app.Fractal_title = uilabel(app.FractalTab);
            app.Fractal_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Fractal_title.HorizontalAlignment = 'center';
            app.Fractal_title.FontWeight = 'bold';
            app.Fractal_title.Position = [8 559 709 22];
            app.Fractal_title.Text = 'Fractal geometries';

            % Create WhyfractalAreyounotcuriousabouttheirtortuosityfactorLabel
            app.WhyfractalAreyounotcuriousabouttheirtortuosityfactorLabel = uilabel(app.FractalTab);
            app.WhyfractalAreyounotcuriousabouttheirtortuosityfactorLabel.FontAngle = 'italic';
            app.WhyfractalAreyounotcuriousabouttheirtortuosityfactorLabel.Position = [8 530 506 22];
            app.WhyfractalAreyounotcuriousabouttheirtortuosityfactorLabel.Text = 'Why fractal ? Are you not curious about their tortuosity factor ? (And fractal are cool anyway)';

            % Create SelectfractalDropDownLabel
            app.SelectfractalDropDownLabel = uilabel(app.FractalTab);
            app.SelectfractalDropDownLabel.HorizontalAlignment = 'right';
            app.SelectfractalDropDownLabel.Position = [8 493 75 22];
            app.SelectfractalDropDownLabel.Text = 'Select fractal';

            % Create SelectfractalDropDown
            app.SelectfractalDropDown = uidropdown(app.FractalTab);
            app.SelectfractalDropDown.Items = {'', 'Menger sponge'};
            app.SelectfractalDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectfractalDropDownValueChanged, true);
            app.SelectfractalDropDown.Position = [98 493 133 22];
            app.SelectfractalDropDown.Value = '';

            % Create Fractal_parameter_UITable
            app.Fractal_parameter_UITable = uitable(app.FractalTab);
            app.Fractal_parameter_UITable.ColumnName = {'Parameter'; 'Value'; 'Range'};
            app.Fractal_parameter_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Fractal_parameter_UITable.RowName = {};
            app.Fractal_parameter_UITable.ColumnEditable = [false true false];
            app.Fractal_parameter_UITable.Position = [8 293 302 185];

            % Create Fractal_GenerateButton
            app.Fractal_GenerateButton = uibutton(app.FractalTab, 'push');
            app.Fractal_GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @Fractal_GenerateButtonPushed, true);
            app.Fractal_GenerateButton.BackgroundColor = [0.8 0.8 0.8];
            app.Fractal_GenerateButton.Enable = 'off';
            app.Fractal_GenerateButton.Position = [8 11 100 40];
            app.Fractal_GenerateButton.Text = 'Generate';

            % Create Fractal_VisualizeButton
            app.Fractal_VisualizeButton = uibutton(app.FractalTab, 'push');
            app.Fractal_VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @Fractal_VisualizeButtonPushed, true);
            app.Fractal_VisualizeButton.BackgroundColor = [0.8 0.8 0.8];
            app.Fractal_VisualizeButton.Enable = 'off';
            app.Fractal_VisualizeButton.Position = [128 11 100 40];
            app.Fractal_VisualizeButton.Text = 'Visualize';

            % Create Fractal_SaveButton
            app.Fractal_SaveButton = uibutton(app.FractalTab, 'push');
            app.Fractal_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Fractal_SaveButtonPushed, true);
            app.Fractal_SaveButton.BackgroundColor = [0.8 0.8 0.8];
            app.Fractal_SaveButton.Enable = 'off';
            app.Fractal_SaveButton.Position = [248 11 100 40];
            app.Fractal_SaveButton.Text = 'Save';

            % Create StatutTextAreaLabel
            app.StatutTextAreaLabel = uilabel(app.FractalTab);
            app.StatutTextAreaLabel.HorizontalAlignment = 'right';
            app.StatutTextAreaLabel.Position = [366 27 37 22];
            app.StatutTextAreaLabel.Text = 'Statut';

            % Create StatutTextArea
            app.StatutTextArea = uitextarea(app.FractalTab);
            app.StatutTextArea.Position = [418 11 299 40];
            app.StatutTextArea.Value = {'Idle'};

            % Create Fractal_HelpTextArea
            app.Fractal_HelpTextArea = uitextarea(app.FractalTab);
            app.Fractal_HelpTextArea.Editable = 'off';
            app.Fractal_HelpTextArea.Visible = 'off';
            app.Fractal_HelpTextArea.Position = [325 372 392 103];

            % Create Rectangular_ChannelsTab
            app.Rectangular_ChannelsTab = uitab(app.TabGroup);
            app.Rectangular_ChannelsTab.Title = 'Channels';
            app.Rectangular_ChannelsTab.ForegroundColor = [1 0 0];

            % Create Rchannel_title
            app.Rchannel_title = uilabel(app.Rectangular_ChannelsTab);
            app.Rchannel_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Rchannel_title.HorizontalAlignment = 'center';
            app.Rchannel_title.FontWeight = 'bold';
            app.Rchannel_title.Position = [8 559 709 22];
            app.Rchannel_title.Text = 'Macro-pore channels ';

            % Create Channel_instructionsLabel
            app.Channel_instructionsLabel = uilabel(app.Rectangular_ChannelsTab);
            app.Channel_instructionsLabel.FontAngle = 'italic';
            app.Channel_instructionsLabel.Position = [8 510 709 42];
            app.Channel_instructionsLabel.Text = {'Structured electrodes with macropore channels (refered in the litterature as dual pore or secondary pore network) enhance'; 'the effective through-plane ionic diffusion, enabling fast-charging and/or ultra thick electrodes. Select first the channel pattern and'; 'enter the associated parameters in the dedicated table.'};

            % Create Channel_GenerateButton
            app.Channel_GenerateButton = uibutton(app.Rectangular_ChannelsTab, 'push');
            app.Channel_GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @Channel_GenerateButtonPushed, true);
            app.Channel_GenerateButton.BackgroundColor = [0.8 0.8 0.8];
            app.Channel_GenerateButton.Enable = 'off';
            app.Channel_GenerateButton.Position = [8 11 100 40];
            app.Channel_GenerateButton.Text = 'Generate';

            % Create Channel_VisualizeButton
            app.Channel_VisualizeButton = uibutton(app.Rectangular_ChannelsTab, 'push');
            app.Channel_VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @Channel_VisualizeButtonPushed, true);
            app.Channel_VisualizeButton.BackgroundColor = [0.8 0.8 0.8];
            app.Channel_VisualizeButton.Enable = 'off';
            app.Channel_VisualizeButton.Position = [128 11 100 40];
            app.Channel_VisualizeButton.Text = 'Visualize';

            % Create Channel_SaveButton
            app.Channel_SaveButton = uibutton(app.Rectangular_ChannelsTab, 'push');
            app.Channel_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Channel_SaveButtonPushed, true);
            app.Channel_SaveButton.BackgroundColor = [0.8 0.8 0.8];
            app.Channel_SaveButton.Enable = 'off';
            app.Channel_SaveButton.Position = [248 11 100 40];
            app.Channel_SaveButton.Text = 'Save';

            % Create StatutTextArea_3Label
            app.StatutTextArea_3Label = uilabel(app.Rectangular_ChannelsTab);
            app.StatutTextArea_3Label.HorizontalAlignment = 'right';
            app.StatutTextArea_3Label.Position = [366 27 37 22];
            app.StatutTextArea_3Label.Text = 'Statut';

            % Create Channel_StatutTextArea
            app.Channel_StatutTextArea = uitextarea(app.Rectangular_ChannelsTab);
            app.Channel_StatutTextArea.Position = [418 11 299 40];
            app.Channel_StatutTextArea.Value = {'Idle'};

            % Create PatternDropDownLabel
            app.PatternDropDownLabel = uilabel(app.Rectangular_ChannelsTab);
            app.PatternDropDownLabel.HorizontalAlignment = 'right';
            app.PatternDropDownLabel.Position = [8 479 44 22];
            app.PatternDropDownLabel.Text = 'Pattern';

            % Create Channel_PatternDropDown
            app.Channel_PatternDropDown = uidropdown(app.Rectangular_ChannelsTab);
            app.Channel_PatternDropDown.Items = {'Rectangular channels along one direction'};
            app.Channel_PatternDropDown.ValueChangedFcn = createCallbackFcn(app, @Channel_PatternDropDownValueChanged, true);
            app.Channel_PatternDropDown.Position = [67 479 281 22];
            app.Channel_PatternDropDown.Value = 'Rectangular channels along one direction';

            % Create Channel_UITable
            app.Channel_UITable = uitable(app.Rectangular_ChannelsTab);
            app.Channel_UITable.ColumnName = {'Parameter'; 'Value'; 'Range'};
            app.Channel_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Channel_UITable.RowName = {};
            app.Channel_UITable.ColumnEditable = [false true false];
            app.Channel_UITable.Position = [8 254 506 220];

            % Create Channel_rectangular1D_Image
            app.Channel_rectangular1D_Image = uiimage(app.Rectangular_ChannelsTab);
            app.Channel_rectangular1D_Image.Position = [516 304 204 169];
            app.Channel_rectangular1D_Image.ImageSource = 'Channels.png';

            % Create Channel_apply_text
            app.Channel_apply_text = uilabel(app.Rectangular_ChannelsTab);
            app.Channel_apply_text.FontAngle = 'italic';
            app.Channel_apply_text.Position = [8 210 709 28];
            app.Channel_apply_text.Text = {'Choose if you want to apply channels on an homogenous material, or on top of an existing microstructure. In the latter case you'; 'must also indicate if you want the channels to cut through particles (load phase Id) or remove them entirely (load particles id).'};

            % Create ApplypatternonLabel
            app.ApplypatternonLabel = uilabel(app.Rectangular_ChannelsTab);
            app.ApplypatternonLabel.HorizontalAlignment = 'right';
            app.ApplypatternonLabel.Position = [8 183 93 22];
            app.ApplypatternonLabel.Text = 'Apply pattern on';

            % Create Channel_applyPatternDropDown
            app.Channel_applyPatternDropDown = uidropdown(app.Rectangular_ChannelsTab);
            app.Channel_applyPatternDropDown.Items = {'', 'Heterogenous microstructure (cut through particles)', 'Heterogenous microstructure (remove particles)', 'Homogenous medium'};
            app.Channel_applyPatternDropDown.ValueChangedFcn = createCallbackFcn(app, @Channel_applyPatternDropDownValueChanged, true);
            app.Channel_applyPatternDropDown.Position = [116 183 349 22];
            app.Channel_applyPatternDropDown.Value = '';

            % Create ClicktoloadmicrostructureButton
            app.ClicktoloadmicrostructureButton = uibutton(app.Rectangular_ChannelsTab, 'push');
            app.ClicktoloadmicrostructureButton.ButtonPushedFcn = createCallbackFcn(app, @ClicktoloadmicrostructureButtonPushed, true);
            app.ClicktoloadmicrostructureButton.BackgroundColor = [0.8 0.8 0.8];
            app.ClicktoloadmicrostructureButton.Enable = 'off';
            app.ClicktoloadmicrostructureButton.Position = [475 181 166 26];
            app.ClicktoloadmicrostructureButton.Text = 'Click to load microstructure';

            % Create ApplyChannel_UITable
            app.ApplyChannel_UITable = uitable(app.Rectangular_ChannelsTab);
            app.ApplyChannel_UITable.ColumnName = {'Parameter'; 'Value'; 'Range'};
            app.ApplyChannel_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.ApplyChannel_UITable.RowName = {};
            app.ApplyChannel_UITable.ColumnEditable = [false true false];
            app.ApplyChannel_UITable.Position = [8 66 506 110];

            % Show the figure after all components are created
            app.DeterministicgenerationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Microstructure_generation_deterministic_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.DeterministicgenerationUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.DeterministicgenerationUIFigure)
        end
    end
end