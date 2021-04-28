classdef Main_menu_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        MATBOXMicrostructureAnalysisToolboxUIFigure  matlab.ui.Figure
        MATBOXNRELMicrostructureAnalysisToolboxLabel  matlab.ui.control.Label
        Version091March182021Label      matlab.ui.control.Label
        SelectthemoduleassociatedwithyouractivityLabel  matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        MicrostructuregenerationButton  matlab.ui.control.Button
        MicrostructurecharacterizationButton_2  matlab.ui.control.Button
        ROIFilteringandsegmentationButton_2  matlab.ui.control.Button
        VisualizationButton             matlab.ui.control.Button
        CorrelationButton               matlab.ui.control.Button
        MeshingButton                   matlab.ui.control.Button
        Label_2                         matlab.ui.control.Label
        AboutButton                     matlab.ui.control.Button
        RepositoryButton                matlab.ui.control.Button
        DocumentationButton             matlab.ui.control.Button
        About_Logo_NREL                 matlab.ui.control.Image
        Generation_DropDown             matlab.ui.control.DropDown
    end

    
    methods (Access = private)
        
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

        % Button pushed function: MicrostructuregenerationButton
        function MicrostructuregenerationButtonPushed(app, event)
            choice = app.Generation_DropDown.Value;
            if strcmp(choice,'Ellipsoid-based stochastic')
                microstructure_generation_ellipsoids_GUI
            elseif strcmp(choice,'Additives (stochastic or energetic method)')
                Microstructure_generation_additives
            elseif strcmp(choice,'Deterministic geometries')
                Microstructure_generation_deterministic
            end
        end

        % Button pushed function: 
        % ROIFilteringandsegmentationButton_2
        function ROIFilteringandsegmentationButton_2Pushed(app, event)
            Segmentation
        end

        % Button pushed function: 
        % MicrostructurecharacterizationButton_2
        function MicrostructurecharacterizationButton_2Pushed(app, event)
            microstructure_characterization_GUI
        end

        % Button pushed function: VisualizationButton
        function VisualizationButtonPushed(app, event)
            %microstructure_visualization_GUI
            microstructure_visualization
        end

        % Button pushed function: CorrelationButton
        function CorrelationButtonPushed(app, event)
            microstructure_correlation_GUI
        end

        % Button pushed function: MeshingButton
        function MeshingButtonPushed(app, event)
            microstructure_meshing
        end

        % Image clicked function: About_Logo_NREL
        function About_Logo_NRELImageClicked(app, event)
            web('https://www.nrel.gov/transportation/energy-storage.html');
        end

        % Button pushed function: AboutButton
        function AboutButtonPushed(app, event)
            web('https://www.nrel.gov/transportation/data-tools.html');
        end

        % Button pushed function: RepositoryButton
        function RepositoryButtonPushed(app, event)
            web('https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox');
        end

        % Button pushed function: DocumentationButton
        function DocumentationButtonPushed(app, event)
            Find_file('NREL_MATBOX_Microstructure_analysis_toolbox_documentation.pdf','MATBOX_Microstructure_analysis_toolbox','Default location is \MATBOX_Microstructure_analysis_toolbox\Documentation\');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create MATBOXMicrostructureAnalysisToolboxUIFigure and hide until all components are created
            app.MATBOXMicrostructureAnalysisToolboxUIFigure = uifigure('Visible', 'off');
            app.MATBOXMicrostructureAnalysisToolboxUIFigure.AutoResizeChildren = 'off';
            app.MATBOXMicrostructureAnalysisToolboxUIFigure.Position = [100 100 439 862];
            app.MATBOXMicrostructureAnalysisToolboxUIFigure.Name = 'MATBOX Microstructure Analysis Toolbox';
            app.MATBOXMicrostructureAnalysisToolboxUIFigure.Icon = 'Icon.png';
            app.MATBOXMicrostructureAnalysisToolboxUIFigure.Resize = 'off';

            % Create MATBOXNRELMicrostructureAnalysisToolboxLabel
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel = uilabel(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel.BackgroundColor = [0.9333 0.7882 0];
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel.HorizontalAlignment = 'center';
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel.FontSize = 22;
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel.FontWeight = 'bold';
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel.Position = [1 788 439 75];
            app.MATBOXNRELMicrostructureAnalysisToolboxLabel.Text = {'MATBOX'; 'NREL Microstructure Analysis Toolbox '};

            % Create Version091March182021Label
            app.Version091March182021Label = uilabel(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.Version091March182021Label.HorizontalAlignment = 'center';
            app.Version091March182021Label.FontAngle = 'italic';
            app.Version091March182021Label.Position = [1 761 439 22];
            app.Version091March182021Label.Text = 'Version: 0.91 - March 18, 2021';

            % Create SelectthemoduleassociatedwithyouractivityLabel
            app.SelectthemoduleassociatedwithyouractivityLabel = uilabel(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.SelectthemoduleassociatedwithyouractivityLabel.HorizontalAlignment = 'center';
            app.SelectthemoduleassociatedwithyouractivityLabel.FontSize = 18;
            app.SelectthemoduleassociatedwithyouractivityLabel.FontWeight = 'bold';
            app.SelectthemoduleassociatedwithyouractivityLabel.Position = [12 727 417 22];
            app.SelectthemoduleassociatedwithyouractivityLabel.Text = 'Select the module associated with your activity';

            % Create Label
            app.Label = uilabel(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.Label.FontSize = 8;
            app.Label.Position = [17 711 408 24];
            app.Label.Text = '___________________________________________________________________________________________';

            % Create MicrostructuregenerationButton
            app.MicrostructuregenerationButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.MicrostructuregenerationButton.ButtonPushedFcn = createCallbackFcn(app, @MicrostructuregenerationButtonPushed, true);
            app.MicrostructuregenerationButton.Icon = 'Icon_generation.png';
            app.MicrostructuregenerationButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.MicrostructuregenerationButton.FontSize = 18;
            app.MicrostructuregenerationButton.FontColor = [1 1 1];
            app.MicrostructuregenerationButton.Tooltip = {'Numerical '; 'generation of microstructures, both active materials and additive phases.'};
            app.MicrostructuregenerationButton.Position = [12 628 417 79];
            app.MicrostructuregenerationButton.Text = 'Microstructure generation';

            % Create MicrostructurecharacterizationButton_2
            app.MicrostructurecharacterizationButton_2 = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.MicrostructurecharacterizationButton_2.ButtonPushedFcn = createCallbackFcn(app, @MicrostructurecharacterizationButton_2Pushed, true);
            app.MicrostructurecharacterizationButton_2.Icon = 'Icon_Characterization.png';
            app.MicrostructurecharacterizationButton_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.MicrostructurecharacterizationButton_2.FontSize = 18;
            app.MicrostructurecharacterizationButton_2.FontColor = [1 1 1];
            app.MicrostructurecharacterizationButton_2.Tooltip = {'Calculate classic properties for battery and fuel cells macromodels: volume fraction, connectivity, specific surface area, particle diameter and tortuosity factor. Choose among a set of different numerical methods to finely evaluate properties value. Perform Representative Volume Element (RVE) analysis and image resolution sensitivity analysis to ckeck the calculation and estimate the error. Calculate morphology parameters: particle shpericity and particle elongation. Calculate geometric tortuosity with the microstructure graph representation.'};
            app.MicrostructurecharacterizationButton_2.Position = [12 440 417 79];
            app.MicrostructurecharacterizationButton_2.Text = 'Microstructure characterization';

            % Create ROIFilteringandsegmentationButton_2
            app.ROIFilteringandsegmentationButton_2 = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.ROIFilteringandsegmentationButton_2.ButtonPushedFcn = createCallbackFcn(app, @ROIFilteringandsegmentationButton_2Pushed, true);
            app.ROIFilteringandsegmentationButton_2.Icon = 'Icon_ROI.png';
            app.ROIFilteringandsegmentationButton_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.ROIFilteringandsegmentationButton_2.FontSize = 18;
            app.ROIFilteringandsegmentationButton_2.FontColor = [1 1 1];
            app.ROIFilteringandsegmentationButton_2.Tooltip = {'Region of interest selection (crop, background detection, rotation) and up/down scaling. Evaluate image quality, enhance contrast, and apply noise-reduction filter to help further segmentation. Apply basic global or local segmentation (manual and Otsu) and calculate microstructure parameters sensitivity with the segmentation threshold.'};
            app.ROIFilteringandsegmentationButton_2.Position = [12 534 417 79];
            app.ROIFilteringandsegmentationButton_2.Text = 'ROI, Filtering, and segmentation';

            % Create VisualizationButton
            app.VisualizationButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.VisualizationButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizationButtonPushed, true);
            app.VisualizationButton.Icon = 'Icon_visualization.png';
            app.VisualizationButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.VisualizationButton.FontSize = 18;
            app.VisualizationButton.FontColor = [1 1 1];
            app.VisualizationButton.Tooltip = {'Visualize microstructures, with 3D orthogonal slice representation. Compare grey-level and segmented volumes using overlay. Load pixel-wise defined results from the microstructure characterisation module and visualize them.'};
            app.VisualizationButton.Position = [12 346 417 79];
            app.VisualizationButton.Text = 'Visualization';

            % Create CorrelationButton
            app.CorrelationButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.CorrelationButton.ButtonPushedFcn = createCallbackFcn(app, @CorrelationButtonPushed, true);
            app.CorrelationButton.Icon = 'Correlation.png';
            app.CorrelationButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.CorrelationButton.FontSize = 18;
            app.CorrelationButton.FontColor = [1 1 1];
            app.CorrelationButton.Tooltip = {'Load microstructure parameters calculated with the microstructure characterization module on different volumes and plot them as a function of each other to estimate their correlation.'};
            app.CorrelationButton.Position = [12 252 417 79];
            app.CorrelationButton.Text = 'Correlation';

            % Create MeshingButton
            app.MeshingButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.MeshingButton.ButtonPushedFcn = createCallbackFcn(app, @MeshingButtonPushed, true);
            app.MeshingButton.Icon = 'Icon_Meshing.png';
            app.MeshingButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.MeshingButton.FontSize = 18;
            app.MeshingButton.FontColor = [1 1 1];
            app.MeshingButton.Tooltip = {'Generate tetrahedron-based mesh(es) of multi-phases (pore, and solids) and multi-domains (current collectors, anode, separator, cathode). Choose between structured and unstructured mesh methods.'};
            app.MeshingButton.Position = [12 158 417 79];
            app.MeshingButton.Text = 'Meshing';

            % Create Label_2
            app.Label_2 = uilabel(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.Label_2.FontSize = 8;
            app.Label_2.Position = [17 119 408 24];
            app.Label_2.Text = '___________________________________________________________________________________________';

            % Create AboutButton
            app.AboutButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.AboutButton.ButtonPushedFcn = createCallbackFcn(app, @AboutButtonPushed, true);
            app.AboutButton.BackgroundColor = [0.8 0.8 0.8];
            app.AboutButton.FontSize = 18;
            app.AboutButton.Position = [288 83 137 29];
            app.AboutButton.Text = 'About';

            % Create RepositoryButton
            app.RepositoryButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.RepositoryButton.ButtonPushedFcn = createCallbackFcn(app, @RepositoryButtonPushed, true);
            app.RepositoryButton.BackgroundColor = [0.8 0.8 0.8];
            app.RepositoryButton.FontSize = 18;
            app.RepositoryButton.Position = [288 48 137 29];
            app.RepositoryButton.Text = 'Repository';

            % Create DocumentationButton
            app.DocumentationButton = uibutton(app.MATBOXMicrostructureAnalysisToolboxUIFigure, 'push');
            app.DocumentationButton.ButtonPushedFcn = createCallbackFcn(app, @DocumentationButtonPushed, true);
            app.DocumentationButton.BackgroundColor = [0.8 0.8 0.8];
            app.DocumentationButton.FontSize = 18;
            app.DocumentationButton.Position = [288 12 137 29];
            app.DocumentationButton.Text = 'Documentation';

            % Create About_Logo_NREL
            app.About_Logo_NREL = uiimage(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.About_Logo_NREL.ImageClickedFcn = createCallbackFcn(app, @About_Logo_NRELImageClicked, true);
            app.About_Logo_NREL.Position = [12 12 264 100];
            app.About_Logo_NREL.ImageSource = 'logo_NREL.png';

            % Create Generation_DropDown
            app.Generation_DropDown = uidropdown(app.MATBOXMicrostructureAnalysisToolboxUIFigure);
            app.Generation_DropDown.Items = {'Ellipsoid-based stochastic', 'Additives (stochastic or energetic method)', 'Deterministic geometries'};
            app.Generation_DropDown.Position = [243 632 183 22];
            app.Generation_DropDown.Value = 'Ellipsoid-based stochastic';

            % Show the figure after all components are created
            app.MATBOXMicrostructureAnalysisToolboxUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Main_menu_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.MATBOXMicrostructureAnalysisToolboxUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MATBOXMicrostructureAnalysisToolboxUIFigure)
        end
    end
end