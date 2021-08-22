classdef Microstructure_generation_stochastic_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        EllipsoidbasedstochasticgenerationmoduleUIFigure  matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        Instructions                    matlab.ui.container.Tab
        NREL_logo                       matlab.ui.control.Image
        Github_button                   matlab.ui.control.Button
        Doc_button                      matlab.ui.control.Button
        SelectasolidphaseandsetthenumberofLabel_19  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_18  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_17  matlab.ui.control.Label
        Hyperlink                       matlab.ui.control.Hyperlink
        SelectasolidphaseandsetthenumberofLabel_16  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_15  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_14  matlab.ui.control.Label
        RunthegenerationalgorithmLabel  matlab.ui.control.Label
        StepsperformedafterthegenerationprogressLabel  matlab.ui.control.Label
        SetalgorithmadditionalparametersLabel  matlab.ui.control.Label
        SetpropertiesdefinedwithanhistogramdistributionLabel  matlab.ui.control.Label
        SetpropertiesdefinedwithascalarLabel  matlab.ui.control.Label
        Instructions_title              matlab.ui.control.Label
        SaveTab                         matlab.ui.container.Tab
        SaveinputsoutpuscomparisonifselectedCheckBox  matlab.ui.control.CheckBox
        ThenselectwhatinformationyouwanttosaveLabel  matlab.ui.control.Label
        SelectwhereyouwanttosaveinputsandresultsLabel  matlab.ui.control.Label
        NosavefolderselectedLabel       matlab.ui.control.Label
        Save_folder_button              matlab.ui.control.Button
        ForphaseandparticlelablelsaveforDropDown  matlab.ui.control.DropDown
        ForphaseandparticlelablelsaveforDropDownLabel  matlab.ui.control.Label
        SaveparticleadditionalinformationCheckBox  matlab.ui.control.CheckBox
        SavecharacterizationresultsifselectedCheckBox  matlab.ui.control.CheckBox
        SavealgorithmprogressionplotifselectedCheckBox  matlab.ui.control.CheckBox
        AdditionalLabel                 matlab.ui.control.Label
        OutputsLabel                    matlab.ui.control.Label
        InputsLabel                     matlab.ui.control.Label
        SaveinmatinputsofthegenerationalgorithmfunctionCheckBox  matlab.ui.control.CheckBox
        Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox  matlab.ui.control.CheckBox
        Savephaselabelinformation3Dtifstackfileuint8formatCheckBox  matlab.ui.control.CheckBox
        Instructions_title_3            matlab.ui.control.Label
        VolumefractionsTab              matlab.ui.container.Tab
        Phase_vf_save                   matlab.ui.control.Button
        Phase_vf_plot                   matlab.ui.control.Button
        PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel  matlab.ui.control.Label
        Volumefraction_UITable          matlab.ui.control.Table
        ThensetthesolidvolumefractionsLabel  matlab.ui.control.Label
        NumberSlice_volumefraction      matlab.ui.control.NumericEditField
        NumberofsolidphaseEditField_2Label  matlab.ui.control.Label
        NumberofsolidphaseEditFieldLabel  matlab.ui.control.Label
        NumberofsolidphaseEditField     matlab.ui.control.NumericEditField
        VoxelsizenmEditFieldLabel       matlab.ui.control.Label
        VoxelsizenmEditField            matlab.ui.control.NumericEditField
        Instructions_title_2            matlab.ui.control.Label
        SetdimensionsLabel              matlab.ui.control.Label
        AlgothmisfasterforLabel         matlab.ui.control.Label
        Phase_domainUITable             matlab.ui.control.Table
        AlgothmisfasterforLabel_2       matlab.ui.control.Label
        SetNumberofsolidphasesLabel     matlab.ui.control.Label
        Phase_LabelnameUItable          matlab.ui.control.Table
        AlgothmisfasterforLabel_3       matlab.ui.control.Label
        SetvolumefractionsalongthemicrostructurethicknessLabel  matlab.ui.control.Label
        AlgothmisfasterforLabel_4       matlab.ui.control.Label
        AlgothmisfasterforLabel_5       matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        DiametersTab                    matlab.ui.container.Tab
        PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_2  matlab.ui.control.Label
        Diameter_save                   matlab.ui.control.Button
        Diameter_plot                   matlab.ui.control.Button
        DxDz_UITable                    matlab.ui.control.Table
        DxDy_UITable                    matlab.ui.control.Table
        AlgothmisfasterforLabel_7       matlab.ui.control.Label
        Rowii1isthediameterhistogramLabel  matlab.ui.control.Label
        Diameters_UITable               matlab.ui.control.Table
        SelectasolidphaseandsetthenumberofLabel_2  matlab.ui.control.Label
        InstructionsTextArea            matlab.ui.control.TextArea
        InstructionsTextAreaLabel       matlab.ui.control.Label
        Diameters_main_UITable          matlab.ui.control.Table
        AlgothmisfasterforLabel_6       matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel  matlab.ui.control.Label
        Selectphase_diameters_DropDown  matlab.ui.control.DropDown
        SelectphaseDropDownLabel        matlab.ui.control.Label
        Instructions_diameter           matlab.ui.control.Label
        RotationsTab                    matlab.ui.container.Tab
        PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_3  matlab.ui.control.Label
        Rotation_save                   matlab.ui.control.Button
        Rotation_plot                   matlab.ui.control.Button
        Rotation_z_UITable              matlab.ui.control.Table
        Rotation_y_UITable              matlab.ui.control.Table
        AlgothmisfasterforLabel_9       matlab.ui.control.Label
        Rowii1isthediameterhistogramLabel_2  matlab.ui.control.Label
        Rotation_x_UITable              matlab.ui.control.Table
        SelectasolidphaseandsetthenumberofLabel_4  matlab.ui.control.Label
        InstructionsTextArea_2          matlab.ui.control.TextArea
        InstructionsTextArea_2Label     matlab.ui.control.Label
        Orientation_main_UITable        matlab.ui.control.Table
        AlgothmisfasterforLabel_8       matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_3  matlab.ui.control.Label
        Selectphase_orientation_DropDown  matlab.ui.control.DropDown
        SelectphaseDropDownLabel_2      matlab.ui.control.Label
        Instructions_orientation        matlab.ui.control.Label
        ParticlevisualizationTab        matlab.ui.container.Tab
        Plot_example_particle_2         matlab.ui.control.Button
        Plot_example_particle           matlab.ui.control.Button
        Parameter_example_UITable       matlab.ui.control.Table
        AlgothmisfasterforLabel_10      matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_5  matlab.ui.control.Label
        Instructions_Visualization      matlab.ui.control.Label
        OverlappingTab                  matlab.ui.container.Tab
        Overlapping_save                matlab.ui.control.Button
        Label_5                         matlab.ui.control.Label
        DisableparticlecontiguitycheckCheckBox  matlab.ui.control.CheckBox
        Minimumvolume_overlapping_UITable  matlab.ui.control.Table
        Label_4                         matlab.ui.control.Label
        overlapping_UITable             matlab.ui.control.Table
        Overlapping_Image               matlab.ui.control.Image
        Label_3                         matlab.ui.control.Label
        Label_2                         matlab.ui.control.Label
        Instructions_overlapping        matlab.ui.control.Label
        OrderTab                        matlab.ui.container.Tab
        Order_save                      matlab.ui.control.Button
        Order_UITable                   matlab.ui.control.Table
        Label_6                         matlab.ui.control.Label
        NumberofgenerationpassEditField  matlab.ui.control.NumericEditField
        NumberofgenerationpassEditFieldLabel  matlab.ui.control.Label
        PseudocodeTextArea              matlab.ui.control.TextArea
        PseudocodeTextAreaLabel         matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_7  matlab.ui.control.Label
        Image                           matlab.ui.control.Image
        SelectasolidphaseandsetthenumberofLabel_6  matlab.ui.control.Label
        Instructions_order              matlab.ui.control.Label
        StopconditionsTab               matlab.ui.container.Tab
        Label_8                         matlab.ui.control.Label
        MaximumwallclocktimesEditField  matlab.ui.control.NumericEditField
        MaximumwallclocktimesEditFieldLabel  matlab.ui.control.Label
        AlgothmisfasterforLabel_11      matlab.ui.control.Label
        stopcondition_action_DropDown   matlab.ui.control.DropDown
        thenLabel                       matlab.ui.control.Label
        iterationsLabel                 matlab.ui.control.Label
        whenaveragedonthelastEditField  matlab.ui.control.NumericEditField
        whenaveragedonthelastEditFieldLabel  matlab.ui.control.Label
        stopcondition_andor_DropDown    matlab.ui.control.DropDown
        particles1goesbelowEditField    matlab.ui.control.NumericEditField
        particles1goesbelowEditFieldLabel  matlab.ui.control.Label
        Ifvolumefractions1goesbelowEditField  matlab.ui.control.NumericEditField
        Ifvolumefractions1goesbelowEditFieldLabel  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_9  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_8  matlab.ui.control.Label
        PlotalgorithmprogressionCheckBox  matlab.ui.control.CheckBox
        Instructions_order_2            matlab.ui.control.Label
        UpscaleoptionalTab              matlab.ui.container.Tab
        ApplyscalingalsotoparticlelabelslowCheckBox  matlab.ui.control.CheckBox
        ApplyscalingtophaselabelfastCheckBox  matlab.ui.control.CheckBox
        Label_7                         matlab.ui.control.Label
        Image2                          matlab.ui.control.Image
        Scalingfactor1upscaling1downscaling1noscalingEditField  matlab.ui.control.NumericEditField
        Scalingfactor1upscaling1downscaling1noscalingEditFieldLabel  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_13  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_12  matlab.ui.control.Label
        Instructions_upscaling          matlab.ui.control.Label
        VerificationandCharacterizationTab  matlab.ui.container.Tab
        CompareoutputsversusinputsCheckBox  matlab.ui.control.CheckBox
        SelectasolidphaseandsetthenumberofLabel_22  matlab.ui.control.Label
        CharacterizeDropDown            matlab.ui.control.DropDown
        CharacterizeDropDownLabel       matlab.ui.control.Label
        PoretortuosityfactorslowforlargevolumesCheckBox  matlab.ui.control.CheckBox
        ConnectivityslowforlargevolumesCheckBox  matlab.ui.control.CheckBox
        VolumefractionsCheckBox         matlab.ui.control.CheckBox
        SelectasolidphaseandsetthenumberofLabel_21  matlab.ui.control.Label
        Instructions_characterization   matlab.ui.control.Label
        GenerationTab                   matlab.ui.container.Tab
        outcometime_UITable             matlab.ui.control.Table
        AlgothmisfasterforLabel_16      matlab.ui.control.Label
        AlgothmisfasterforLabel_15      matlab.ui.control.Label
        Tortuosity_UITable              matlab.ui.control.Table
        AlgothmisfasterforLabel_14      matlab.ui.control.Label
        Volume_fraction_UITable         matlab.ui.control.Table
        AlgothmisfasterforLabel_13      matlab.ui.control.Label
        AlgothmisfasterforLabel_12      matlab.ui.control.Label
        Connectivity_UITable            matlab.ui.control.Table
        SelectasolidphaseandsetthenumberofLabel_11  matlab.ui.control.Label
        NumberofrunsEditField           matlab.ui.control.NumericEditField
        NumberofrunsEditFieldLabel      matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_10  matlab.ui.control.Label
        Instructions_generation         matlab.ui.control.Label
        Generate_Button                 matlab.ui.control.Button
        VisualizationTab                matlab.ui.container.Tab
        DviewcolormapDropDown           matlab.ui.control.DropDown
        DviewcolormapDropDownLabel      matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_26  matlab.ui.control.Label
        Plot_particle_label_3D          matlab.ui.control.Button
        Plot_phase_label_3D             matlab.ui.control.Button
        DviewLabel                      matlab.ui.control.Label
        DslicesLabel                    matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_25  matlab.ui.control.Label
        SelectasolidphaseandsetthenumberofLabel_24  matlab.ui.control.Label
        Plot_particle_label             matlab.ui.control.Button
        Plot_phase_label                matlab.ui.control.Button
        SelectasolidphaseandsetthenumberofLabel_23  matlab.ui.control.Label
        Visualization_UITable           matlab.ui.control.Table
        Instructions_generation_3       matlab.ui.control.Label
        HowtoquoteTab                   matlab.ui.container.Tab
        TextArea                        matlab.ui.control.TextArea
        SelectasolidphaseandsetthenumberofLabel_20  matlab.ui.control.Label
        Instructions_generation_2       matlab.ui.control.Label
    end

    
    properties (Access = private)
        voxel_size = []; % Voxel size (before upscaling)
        domain_size = [];
        number_solidphase = [];
        currentphase_diameter_setup = 1;
        currentphase_elongation_setup = 1;
        sav_diameter = [];
        sav_elongation = [];
        
        % Input of the generation algorithm
        phase = []; 
        Maximum_overlapping = [];
        Minimum_particle_volume_conservated = [];
        disable_contiguituy_check = [];        
        
        % Output
        generation_result = [];

        % Save options
        Savefolder = [];
        
        

        
    end
    
    methods (Access = private)
        
        function [] = update_Phase_Labelname(app)
            str_name = {};
            for k=1:1:app.number_solidphase
                str_name = [str_name; ['Phase ' num2str(k)]];
            end
            app.Phase_LabelnameUItable.Data = [ num2cell([1:1:app.number_solidphase]') num2cell([1:1:app.number_solidphase]') str_name ];
        end
                
        function [] = update_volumefraction_table(app)
            c_ = {'Normalized position along direction 3'};
            for k=1:1:app.number_solidphase
                c_(k+1,1) = app.Phase_LabelnameUItable.Data(k,3);
            end
            c_(k+2,1) = {'Porosity'};
            app.Volumefraction_UITable.ColumnName = c_;
            pos_ = round( linspace(0,1,app.NumberSlice_volumefraction.Value)' ,3);
            vf = round( ones(app.NumberSlice_volumefraction.Value, app.number_solidphase+1) .* 1/(app.number_solidphase+1) , 3);
            app.Volumefraction_UITable.Data = [num2cell(pos_) num2cell(vf)];
            
            edit_ = ones(1,app.number_solidphase+2); edit_(1,end) = 0; % Can edit all except porosity
            app.Volumefraction_UITable.ColumnEditable = logical(edit_);            
        end
        
        function [] = initialize_diameter_table(app)
            c_ = {'Normalized position along direction 3'};
            number_diameter = cell2mat(app.Diameters_main_UITable.Data(1,3));
            for k=1:1:number_diameter
                c_(k+1) = {['Dx ' num2str(k)]};
            end
            c_(k+2) = {'Total'};     
            number_slice = cell2mat(app.Diameters_main_UITable.Data(1,2));
            d = zeros(number_slice+1,number_diameter+2)+NaN;
            d(1,2:end-1) = 17:2:17+2*(number_diameter-1);
            d(2:end,1) = round(linspace(0,1,number_slice), 3);
            d(2:end,2:end-1) = ones(number_slice,number_diameter) *round(100/number_diameter,3);
            if number_diameter>1
                d(2:end,end-1) = 100 - (number_diameter-1) * round(100/number_diameter,3); % Make sure total is exactly 100
            end
            d(2:end,end) = sum(d(2:end,2:end-1),2);
            app.Diameters_UITable.ColumnName = c_;
            app.Diameters_UITable.Data = d;
            edit_ = ones(1,number_diameter+2); edit_(1,end) = 0;
            app.Diameters_UITable.ColumnEditable = logical(edit_);   
        end
        
        function [] = initialize_dxdy_table(app)
            c_ = {'Normalized position along direction 3'};
            number_ratio = cell2mat(app.Diameters_main_UITable.Data(2,3));
            for k=1:1:number_ratio
                c_(k+1) = {['Dx/Dy ' num2str(k)]};
            end
            c_(k+2) = {'Total'};     
            number_slice = cell2mat(app.Diameters_main_UITable.Data(2,2));
            d = zeros(number_slice+1,number_ratio+2)+NaN;
            d(1,2:end-1) = 1:0.5:1+0.5*(number_ratio-1);
            d(2:end,1) = round(linspace(0,1,number_slice), 3);
            d(2:end,2:end-1) = ones(number_slice,number_ratio) *round(100/number_ratio,3);
            if number_ratio>1
                d(2:end,end-1) = 100 - (number_ratio-1) * round(100/number_ratio,3); % Make sure total is exactly 100
            end
            d(2:end,end) = sum(d(2:end,2:end-1),2);
            app.DxDy_UITable.ColumnName = c_;
            app.DxDy_UITable.Data = d;       
            edit_ = ones(1,number_ratio+2); edit_(1,end) = 0;
            app.DxDy_UITable.ColumnEditable = logical(edit_);               
        end
        
        function [] = initialize_dxdz_table(app)
            c_ = {'Normalized position along direction 3'};
            number_ratio = cell2mat(app.Diameters_main_UITable.Data(3,3));
            for k=1:1:number_ratio
                c_(k+1) = {['Dx/Dz ' num2str(k)]};
            end
            c_(k+2) = {'Total'};     
            number_slice = cell2mat(app.Diameters_main_UITable.Data(3,2));
            d = zeros(number_slice+1,number_ratio+2)+NaN;
            d(1,2:end-1) = 1:0.5:1+0.5*(number_ratio-1);
            d(2:end,1) = round(linspace(0,1,number_slice), 3);
            d(2:end,2:end-1) = ones(number_slice,number_ratio) *round(100/number_ratio,3);
            if number_ratio>1
                d(2:end,end-1) = 100 - (number_ratio-1) * round(100/number_ratio,3); % Make sure total is exactly 100
            end
            d(2:end,end) = sum(d(2:end,2:end-1),2);
            app.DxDz_UITable.ColumnName = c_;
            app.DxDz_UITable.Data = d;              
            edit_ = ones(1,number_ratio+2); edit_(1,end) = 0;
            app.DxDz_UITable.ColumnEditable = logical(edit_); 
        end
        
        function [] = initialize_Ox_table(app)
            c_ = {'Normalized position along direction 3'};
            number_orientation = cell2mat(app.Orientation_main_UITable.Data(1,3));
            for k=1:1:number_orientation
                c_(k+1) = {['Rx ' num2str(k)]};
            end
            c_(k+2) = {'Total'};     
            number_slice = cell2mat(app.Orientation_main_UITable.Data(1,2));
            d = zeros(number_slice+1,number_orientation+2)+NaN;
            initial_degree = linspace(0,180,number_orientation+1);
            d(1,2:end-1) = initial_degree(2:end); % do not take 0
            d(2:end,1) = round(linspace(0,1,number_slice), 3);
            d(2:end,2:end-1) = ones(number_slice,number_orientation) *round(100/number_orientation,3);
            if number_orientation>1
                d(2:end,end-1) = 100 - (number_orientation-1) * round(100/number_orientation,3); % Make sure total is exactly 100
            end
            d(2:end,end) = sum(d(2:end,2:end-1),2);
            app.Rotation_x_UITable.ColumnName = c_;
            app.Rotation_x_UITable.Data = d;
            edit_ = ones(1,number_orientation+2); edit_(1,end) = 0;
            app.Rotation_x_UITable.ColumnEditable = logical(edit_);              
        end
        
        function [] = initialize_Oy_table(app)
            c_ = {'Normalized position along direction 3'};
            number_orientation = cell2mat(app.Orientation_main_UITable.Data(2,3));
            for k=1:1:number_orientation
                c_(k+1) = {['Ry ' num2str(k)]};
            end
            c_(k+2) = {'Total'};     
            number_slice = cell2mat(app.Orientation_main_UITable.Data(2,2));
            d = zeros(number_slice+1,number_orientation+2)+NaN;
            initial_degree = linspace(0,180,number_orientation+1);
            d(1,2:end-1) = initial_degree(2:end); % do not take 0
            d(2:end,1) = round(linspace(0,1,number_slice), 3);
            d(2:end,2:end-1) = ones(number_slice,number_orientation) *round(100/number_orientation,3);
            if number_orientation>1
                d(2:end,end-1) = 100 - (number_orientation-1) * round(100/number_orientation,3); % Make sure total is exactly 100
            end
            d(2:end,end) = sum(d(2:end,2:end-1),2);
            app.Rotation_y_UITable.ColumnName = c_;
            app.Rotation_y_UITable.Data = d;
            edit_ = ones(1,number_orientation+2); edit_(1,end) = 0;
            app.Rotation_y_UITable.ColumnEditable = logical(edit_);               
        end
        
        function [] = initialize_Oz_table(app)
            c_ = {'Normalized position along direction 3'};
            number_orientation = cell2mat(app.Orientation_main_UITable.Data(3,3));
            for k=1:1:number_orientation
                c_(k+1) = {['Rz ' num2str(k)]};
            end
            c_(k+2) = {'Total'};     
            number_slice = cell2mat(app.Orientation_main_UITable.Data(3,2));
            d = zeros(number_slice+1,number_orientation+2)+NaN;
            initial_degree = linspace(0,180,number_orientation+1);
            d(1,2:end-1) = initial_degree(2:end); % do not take 0
            d(2:end,1) = round(linspace(0,1,number_slice), 3);
            d(2:end,2:end-1) = ones(number_slice,number_orientation) *round(100/number_orientation,3);
            if number_orientation>1
                d(2:end,end-1) = 100 - (number_orientation-1) * round(100/number_orientation,3); % Make sure total is exactly 100
            end
            d(2:end,end) = sum(d(2:end,2:end-1),2);
            app.Rotation_z_UITable.ColumnName = c_;
            app.Rotation_z_UITable.Data = d;
            edit_ = ones(1,number_orientation+2); edit_(1,end) = 0;
            app.Rotation_z_UITable.ColumnEditable = logical(edit_);              
        end        
        
        
        function [] = Initialize_order_table(app)
            number_pass = app.NumberofgenerationpassEditField.Value;
            n_phase = length(app.Phase_LabelnameUItable.Data(:,3));
            for k=1:1:number_pass
                r_(k,1) = {['Pass ' num2str(k)]};
            end
            app.Order_UITable.ColumnName = app.Phase_LabelnameUItable.Data(:,3);
            app.Order_UITable.RowName = r_;
            order_ = zeros(number_pass,n_phase);
            for k=1:1:n_phase
                order_(:,k) = linspace(0.25,1,number_pass)';
            end
            app.Order_UITable.Data = order_;
            edit_ = ones(1,n_phase);
            app.Order_UITable.ColumnEditable = logical(edit_);
        end
        
        function Res = Connectivity_mainclusterpercent(app,M,label)
            % Simplified version of the function_connectivity of the characterization module
            binary_phase=zeros(size(M)); % Initialization
            binary_phase(M == label) = 1;
            connected_id = 1;
            isolated_id = 2;
            unknown_id = 3;
            [Connectivity_structure] = Function_Connectivity_Algorithm(binary_phase, 1, connected_id, unknown_id, isolated_id); % Call algorithm
            Res = Connectivity_structure.Clusters_LargestIsolatedUnknown.main_cluster_phasefraction;
        end

        
        function t = Build_table_fromUItab(app,uit)
            d = uit.Data;
            [~,n_column] = size(d);
            str = 'table(d(:,1)';
            for k=2:1:n_column
                str = [str ',d(:,' num2str(k) ')'];
            end
            str = [str ',''VariableNames'',(uit.ColumnName)'');'];
            t = eval(str);
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Phase tab
            app.voxel_size = app.VoxelsizenmEditField.Value;
            app.domain_size = [100 100 100];
            app.Phase_domainUITable.Data = [ num2cell([1;2;3]) num2cell(app.domain_size') num2cell([app.domain_size.*app.voxel_size./1000]') ];
            % Phase name
            app.number_solidphase = app.NumberofsolidphaseEditField.Value;
            app.update_Phase_Labelname;
            % Volume fraction table
            app.update_volumefraction_table;
            % Particle visualization table
            app.Parameter_example_UITable.Data = [{'Diameter Dx';'Diameter ratio Dx/Dy';'Diameter ratio Dx/Dz';'Rotation normal with x-axis Rx';'Rotation normal with y-axis Ry';'Rotation normal with z-axis Rz'} num2cell([17;1;1;180;180;180]) {'Number of voxel';'[]';'[]';'Degree';'Degree';'Degree'} {'>1, integer';']0,+inf[';']0,+inf[';']0,180]';']0,180]';']0,180]'} ];
            % Visualization table
            app.Visualization_UITable.ColumnFormat = {'numeric','logical','logical'};
            
            
        end

        % Cell edit callback: Phase_domainUITable
        function Phase_domainUITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(2)==2 % Modified number of voxel
                if isnan(newData) || newData<1 || round(newData)~=newData % Incorrect input
                    app.Phase_domainUITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                    f = msgbox('Invalid Value: must be integer >=1', 'Error','warn');
                else
                    app.Phase_domainUITable.Data(:,3) = num2cell([ cell2mat(app.Phase_domainUITable.Data(:,2)).*app.voxel_size./1000]');
                end


            end



        end

        % Value changed function: VoxelsizenmEditField
        function VoxelsizenmEditFieldValueChanged(app, event)
            value = app.VoxelsizenmEditField.Value;
            if isnan(value) || value<0 % Incorrect input
                app.VoxelsizenmEditField.Value = app.voxel_size;
                f = msgbox('Invalid Value: must be >0', 'Error','warn');
            else
                app.voxel_size = value;
                app.Phase_domainUITable.Data(:,3) = num2cell([ cell2mat(app.Phase_domainUITable.Data(:,2)).*app.voxel_size./1000]');
            end
            
            
        end

        % Value changed function: NumberofsolidphaseEditField
        function NumberofsolidphaseEditFieldValueChanged(app, event)
            value = app.NumberofsolidphaseEditField.Value;
            if isnan(value) || value<1 || round(value)~=value % Incorrect input
                app.NumberofsolidphaseEditField.Value = app.number_solidphase;
                f = msgbox('Invalid Value: must be integer >=1', 'Error','warn');
            else
                app.number_solidphase = value;
                app.update_Phase_Labelname;
                app.update_volumefraction_table;
            end
            
        end

        % Cell edit callback: Phase_LabelnameUItable
        function Phase_LabelnameUItableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(2)==2 % Modified label
                if isnan(newData) || newData<1 || round(newData)~=newData % Incorrect input
                    app.Phase_LabelnameUItable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                    f = msgbox({'Invalid Value: must be integer >=1','Label 0 is reserved for background'}, 'Error','warn');
                else
                    tentative = cell2mat(app.Phase_LabelnameUItable.Data(:, 2));
                    if length(tentative)~=length(unique(tentative)) % Not unique
                        app.Phase_LabelnameUItable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                        f = msgbox('Invalid Value: all label must be different', 'Error','warn');
                    else
                        app.update_volumefraction_table;
                    end
                end
            elseif indices(2)==3 % Modified name
                app.update_volumefraction_table;
            end
        end

        % Value changed function: NumberSlice_volumefraction
        function NumberSlice_volumefractionValueChanged(app, event)
            value = app.NumberSlice_volumefraction.Value;
            if isnan(value) || value<2 || round(value)~=value % Incorrect input
                f = msgbox('Invalid Value: must be integer >=2', 'Error','warn');
            else
                app.update_volumefraction_table;
            end  
        end

        % Cell edit callback: Volumefraction_UITable
        function Volumefraction_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if isnan(newData)
                f = msgbox('Invalid Value: must be float between 0 and 1', 'Error','warn');
                app.Volumefraction_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
            else
                pos_ = round(cell2mat(app.Volumefraction_UITable.Data(:,1)), 3);
                vf_solid = round(cell2mat(app.Volumefraction_UITable.Data(:,2:end-1)), 3);
                % Check position
                if min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.Volumefraction_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end
                % Check solid volume fraction
                if min(min(vf_solid))<0 || max(max(vf_solid))>1 || sum(sum(vf_solid,2)>=1)>=1
                    f = msgbox({'Invalid Value: must be float between 0 and 1','and sum of solid volume fraction below 1'}, 'Error','warn');
                    app.Volumefraction_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end
                % Round
                app.Volumefraction_UITable.Data(:,1) = num2cell( round( cell2mat(app.Volumefraction_UITable.Data(:,1)), 3));
                app.Volumefraction_UITable.Data(:,2:end-1) = num2cell( round( cell2mat(app.Volumefraction_UITable.Data(:,2:end-1)), 3));
                % Deduce porosity
                vf_solid = cell2mat(app.Volumefraction_UITable.Data(:,2:end-1));
                vf_pore = 1-sum(vf_solid,2);
                app.Volumefraction_UITable.Data(:,end) = num2cell(vf_pore);
            end
            

        end

        % Button pushed function: Phase_vf_plot
        function Phase_vf_plotButtonPushed(app, event)
            % Data
            pos_ = cell2mat(app.Volumefraction_UITable.Data(:,1));
            vf_solid = cell2mat(app.Volumefraction_UITable.Data(:,2:end-1));
            vf_pore = cell2mat(app.Volumefraction_UITable.Data(:,end));
            phase_name = app.Phase_LabelnameUItable.Data(:,3);
            [~,n_phase] = size(vf_solid);
            
            Fig = figure; % Create figure
            Fig.Name= 'Volume fractions input'; % Figure name
            Fig.Color='white'; % Background colour
            ax_=axes('parent',Fig);
            h_title = title('Volume fractions input');
            hold(ax_,'on'); % Active subplot
            for k=1:1:n_phase
                plot(pos_,vf_solid(:,k),'LineWidth',2,'LineStyle','-','MarkerSize',12,'Marker','+','DisplayName',char(phase_name(k)));
            end
            plot(pos_,vf_pore,'LineWidth',2,'LineStyle','--','MarkerSize',12,'Marker','+','Color','k','DisplayName','Porosity');
            xlabel('Normalized position along direction 3 (thickness)');
            ylabel('Volume fraction');
            legend('Location','best');
            grid(ax_,'on'); % Display grid
            set(ax_,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
            h_title.FontSize = 14; % Set title fontsize
            hold(ax_,'off'); % Relase axe
        end

        % Button pushed function: Phase_vf_save
        function Phase_vf_saveButtonPushed(app, event)
            % Save in phase
            pos_ = cell2mat(app.Volumefraction_UITable.Data(:,1));
            vf_solid = cell2mat(app.Volumefraction_UITable.Data(:,2:end-1));
            phase_name = app.Phase_LabelnameUItable.Data(:,3);
            phase_label = app.Phase_LabelnameUItable.Data(:,2);
            [~,n_phase] = size(vf_solid);
            for k = 1:1:n_phase
                app.phase(k).name = char(phase_name(k));
                app.phase(k).code = cell2mat(phase_label(k));
                app.phase(k).volumefraction.along_3rd_axis = [pos_ vf_solid(:,k)]';
            end                    
            
            app.domain_size = cell2mat(app.Phase_domainUITable.Data(:,2)');
            % % Diameter tab
            % Prepare diameter main tab
            app.Selectphase_diameters_DropDown.Items = phase_name;
            app.Selectphase_diameters_DropDown.Value = phase_name(1);
            app.Diameters_main_UITable.Data = [{'Diameter Dx';'Diameter ratio Dx/Dy';'Diameter ratio Dx/Dz'}  num2cell([2;2;2]) num2cell([1;1;1])];
            % Diameter and elongation tab
            app.initialize_diameter_table;
            app.initialize_dxdy_table;
            app.initialize_dxdz_table;
            % Save default values
            for k=1:1:n_phase
                app.currentphase_diameter_setup = k;
                app.Diameter_saveButtonPushed;
            end
            app.currentphase_diameter_setup = 1;
            % Enable GUI
            app.Selectphase_diameters_DropDown.Enable = 'on'; app.SelectphaseDropDownLabel.Enable = 'on';
            app.Diameters_main_UITable.Enable = 'on';
            app.Diameters_UITable.Enable = 'on'; app.DxDz_UITable.Enable = 'on'; app.DxDy_UITable.Enable = 'on';
            app.Diameter_plot.Enable = 'on'; app.Diameter_save.Enable = 'on';
            
            % % Orientations tab
            % Prepare orientation main tab
            app.Selectphase_orientation_DropDown.Items = phase_name;
            app.Selectphase_orientation_DropDown.Value = phase_name(1);
            app.Orientation_main_UITable.Data = [{'Rotation normal to x-axis (Rx)';'Rotation normal to y-axis (Ry)';'Rotation normal to z-axis (Rz)'}  num2cell([2;2;2]) num2cell([1;1;1])];
            % Diameter and elongation tab
            app.initialize_Ox_table;
            app.initialize_Oy_table;
            app.initialize_Oz_table;
            % Save default values
            for k=1:1:n_phase
                app.currentphase_elongation_setup = k;
                app.Rotation_saveButtonPushed;
            end
            app.currentphase_elongation_setup = 1;
            % Enable GUI
            app.Selectphase_orientation_DropDown.Enable = 'on'; app.SelectphaseDropDownLabel_2.Enable = 'on';
            app.Orientation_main_UITable.Enable = 'on';
            app.Rotation_x_UITable.Enable = 'on'; app.Rotation_y_UITable.Enable = 'on'; app.Rotation_z_UITable.Enable = 'on';
            app.Rotation_plot.Enable = 'on'; app.Rotation_save.Enable = 'on';            
            
            % % Overlapping tab
            % Prepare tab
            app.overlapping_UITable.ColumnName = app.Phase_LabelnameUItable.Data(:,3);
            app.overlapping_UITable.RowName = app.Phase_LabelnameUItable.Data(:,3);
            app.overlapping_UITable.Data = zeros(n_phase,n_phase)+0.5;
            edit_ = ones(1,n_phase);
            app.overlapping_UITable.ColumnEditable = logical(edit_);                       
            app.overlapping_UITable.Enable = 'on'; app.overlapping_UITable.Visible = 'on';
            app.Minimumvolume_overlapping_UITable.ColumnName = app.Phase_LabelnameUItable.Data(:,3);
            app.Minimumvolume_overlapping_UITable.Data = zeros(1,n_phase)+0.5;
            app.Minimumvolume_overlapping_UITable.ColumnEditable = logical(edit_);   
            app.Minimumvolume_overlapping_UITable.Enable = 'on'; app.Minimumvolume_overlapping_UITable.Visible = 'on';
            app.Overlapping_save.Enable = 'on';
            % Save default value
            app.Overlapping_saveButtonPushed;
            
            % % Order tab
            app.Initialize_order_table;
            app.Order_UITable.Enable = 'on'; app.Order_UITable.Visible = 'on';
            app.Order_save.Enable = 'on';
            % Save default value
            app.Order_saveButtonPushed;            
            
            % % Generation tab
            app.Volume_fraction_UITable.ColumnName = [{'Run#'}, {'Pores'} app.Phase_LabelnameUItable.Data(:,3)'];
            app.Connectivity_UITable.ColumnName = [{'Run#'}, {'Pores'} app.Phase_LabelnameUItable.Data(:,3)', {'Union of solid phases'}];
                        
            app.Generate_Button.Enable = 'on';
            
            % % Move to diameter tab  
            app.TabGroup.SelectedTab = app.DiametersTab;           
        end

        % Value changed function: Selectphase_diameters_DropDown
        function Selectphase_diameters_DropDownValueChanged(app, event)
            value = app.Selectphase_diameters_DropDown.Value;
            app.currentphase_diameter_setup = find(ismember(app.Selectphase_diameters_DropDown.Items,value));
            % Re-load previous state
            app.Diameters_main_UITable.Data = app.sav_diameter(app.currentphase_diameter_setup).D_data;
            app.Diameters_UITable.Data = app.sav_diameter(app.currentphase_diameter_setup).Dx_data;
            app.Diameters_UITable.ColumnName = app.sav_diameter(app.currentphase_diameter_setup).Dx_ColumnName;
            app.DxDy_UITable.Data = app.sav_diameter(app.currentphase_diameter_setup).DxDy_data;
            app.DxDy_UITable.ColumnName = app.sav_diameter(app.currentphase_diameter_setup).DxDy_ColumnName;
            app.DxDz_UITable.Data = app.sav_diameter(app.currentphase_diameter_setup).DxDz_data;
            app.DxDz_UITable.ColumnName = app.sav_diameter(app.currentphase_diameter_setup).DxDz_ColumnName;              
        end

        % Cell edit callback: Diameters_main_UITable
        function Diameters_main_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(2)==2 % Number of slices
                if isnan(newData) || newData<2 || round(newData)~=newData
                    f = msgbox('Invalid Value: must be integer >=2', 'Error','warn');
                    app.Diameters_main_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end
            elseif indices(2)==3 % Number of values
                if isnan(newData) || newData<1 || round(newData)~=newData
                    f = msgbox('Invalid Value: must be integer >=1', 'Error','warn');
                    app.Diameters_main_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end                
            end
            if indices(1)==1
                app.initialize_diameter_table;
            elseif indices(1)==2
                app.initialize_dxdy_table;
            elseif indices(1)==3
                app.initialize_dxdz_table;
            end            
        end

        % Cell edit callback: Diameters_UITable
        function Diameters_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 && indices(2)>1 % Diameter
                diam = app.Diameters_UITable.Data(1,2:end-1);
                if isnan(newData) || newData<1 || round(newData)~=newData || sum(diam==sort(diam))~=length(diam) || length(diam)~=length(unique(diam) )
                    f = msgbox({'Invalid Value: must be integer >1','Diameters must be increasing and with no repeated value'}, 'Error','warn');
                    app.Diameters_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)==1 % Position
                pos_ = app.Diameters_UITable.Data(2:end,1);
                if isnan(newData) || min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.Diameters_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)>1 % Percentage
                rowp = round(app.Diameters_UITable.Data(indices(1),2:end-1), 3); % Round it
                if isnan(newData) || newData<0 || newData>100
                    f = msgbox({'Invalid Value: must be float between 0 and 100','with sum equal to 100'}, 'Error','warn');
                    app.Diameters_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                else
                    number_diameter = cell2mat(app.Diameters_main_UITable.Data(1,3));
                    if number_diameter>1 && abs(sum(rowp)-100)<0.1
                        rowp(end) = 100-sum(rowp(1:end-1)); % Make sure 100%
                    end
                    app.Diameters_UITable.Data(indices(1),2:end-1) = rowp;
                    app.Diameters_UITable.Data(indices(1),end) = sum(rowp);
                    if sum(rowp)~=100
                        app.Diameter_plot.Enable = 'off';
                        app.Diameter_save.Enable = 'off';
                    else
                        app.Diameter_plot.Enable = 'on';
                        app.Diameter_save.Enable = 'on'; 
                    end
                end
            end
            app.Diameters_UITable.Data(1,1) = NaN;
            app.Diameters_UITable.Data(1,end) = NaN;
        end

        % Cell edit callback: DxDy_UITable
        function DxDy_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 && indices(2)>1 % Ratio
                ratio = app.DxDy_UITable.Data(1,2:end-1);
                if isnan(newData) || newData<0 || sum(ratio==sort(ratio))~=length(ratio) || length(ratio)~=length(unique(ratio) )
                    f = msgbox({'Invalid Value: must be >0','Ratio must be increasing and with no repeated value'}, 'Error','warn');
                    app.DxDy_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)==1 % Position
                pos_ = app.DxDy_UITable.Data(2:end,1);
                if isnan(newData) || min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.DxDy_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)>1 % Percentage
                rowp = round(app.DxDy_UITable.Data(indices(1),2:end-1), 3); % Round it
                if isnan(newData) || newData<0 || newData>100
                    f = msgbox({'Invalid Value: must be float between 0 and 100','with sum equal to 100'}, 'Error','warn');
                    app.DxDy_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                else
                    number_ratio = cell2mat(app.Diameters_main_UITable.Data(2,3));
                    if number_ratio>1 && abs(sum(rowp)-100)<0.1
                        rowp(end) = 100-sum(rowp(1:end-1)); % Make sure 100%
                    end
                    app.DxDy_UITable.Data(indices(1),2:end-1) = rowp;
                    app.DxDy_UITable.Data(indices(1),end) = sum(rowp);
                    if sum(rowp)~=100
                        app.Diameter_plot.Enable = 'off';
                        app.Diameter_save.Enable = 'off';
                    else
                        app.Diameter_plot.Enable = 'on';
                        app.Diameter_save.Enable = 'on'; 
                    end
                end
            end
            app.DxDy_UITable.Data(1,1) = NaN;
            app.DxDy_UITable.Data(1,end) = NaN;            
        end

        % Cell edit callback: DxDz_UITable
        function DxDz_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 && indices(2)>1 % Ratio
                ratio = app.DxDz_UITable.Data(1,2:end-1);
                if isnan(newData) || newData<0 || sum(ratio==sort(ratio))~=length(ratio) || length(ratio)~=length(unique(ratio) )
                    f = msgbox({'Invalid Value: must be >0','Ratio must be increasing and with no repeated value'}, 'Error','warn');
                    app.DxDz_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)==1 % Position
                pos_ = app.DxDz_UITable.Data(2:end,1);
                if isnan(newData) || min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.DxDz_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)>1 % Percentage
                rowp = round(app.DxDz_UITable.Data(indices(1),2:end-1), 3); % Round it
                if isnan(newData) || newData<0 || newData>100
                    f = msgbox({'Invalid Value: must be float between 0 and 100','with sum equal to 100'}, 'Error','warn');
                    app.DxDz_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                else
                    number_ratio = cell2mat(app.Diameters_main_UITable.Data(3,3));
                    if number_ratio>1 && abs(sum(rowp)-100)<0.1
                        rowp(end) = 100-sum(rowp(1:end-1)); % Make sure 100%
                    end
                    app.DxDz_UITable.Data(indices(1),2:end-1) = rowp;
                    app.DxDz_UITable.Data(indices(1),end) = sum(rowp);
                    if sum(rowp)~=100
                        app.Diameter_plot.Enable = 'off';
                        app.Diameter_save.Enable = 'off';
                    else
                        app.Diameter_plot.Enable = 'on';
                        app.Diameter_save.Enable = 'on'; 
                    end
                end
            end
            app.DxDz_UITable.Data(1,1) = NaN;
            app.DxDz_UITable.Data(1,end) = NaN;             
        end

        % Button pushed function: Diameter_plot
        function Diameter_plotButtonPushed(app, event)
            Fig = figure; % Create figure
            Fig.Name= 'Particle diameter and elongation inputs'; % Figure name
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)/2]); % Full screen figure
            Fig.Color='white'; % Background colour
            for id_axe=1:1:3
                if id_axe==1
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig);
                    hold(sub_axes,'on'); % Active subplot
                    h_title=title ('Diameter Dx (voxel length)');
                    x = app.Diameters_UITable.Data(2:end,1);
                    y = app.Diameters_UITable.Data(2:end,2:end-1);
                    A1 = area(x,y);
                    [~,number_diameter] = size(y);
                    for k=1:1:number_diameter
                        set(A1(k),'DisplayName',['Dx ' num2str(k) ' = ' num2str(app.Diameters_UITable.Data(1,k+1))]);
                    end
                elseif id_axe==2
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig);
                    hold(sub_axes,'on'); % Active subplo
                    h_title=title ('Diameter ratio Dx/Dy');
                    x = app.DxDy_UITable.Data(2:end,1);
                    y = app.DxDy_UITable.Data(2:end,2:end-1);  
                    A2 = area(x,y);
                    [~,number_ratio] = size(y);
                    for k=1:1:number_ratio
                        set(A2(k),'DisplayName',['Dx/Dy ' num2str(k) ' = ' num2str(app.DxDy_UITable.Data(1,k+1))]);
                    end                    
                elseif id_axe==3
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig);
                    hold(sub_axes,'on'); % Active subplo
                    h_title=title ('Diameter ratio Dx/Dz');
                    x = app.DxDz_UITable.Data(2:end,1);
                    y = app.DxDz_UITable.Data(2:end,2:end-1);    
                    A3 = area(x,y);
                    [~,number_ratio] = size(y);
                    for k=1:1:number_ratio
                        set(A3(k),'DisplayName',['Dx/Dz ' num2str(k) ' = ' num2str(app.DxDz_UITable.Data(1,k+1))]);
                    end                          
                end
                xlabel('Normalized position along direction 3 (thickness)');
                ylabel('%');
                legend('Location','best');
                h_maintitle=sgtitle('Particle diameter and elongation inputs');
                set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
                set(h_title,'FontSize',14)
                set(h_maintitle,'FontSize',16,'FontName','Times new roman')
                hold(sub_axes,'off'); % Relase axe,'off'); % Relase axe
            end     
        end

        % Button pushed function: Diameter_save
        function Diameter_saveButtonPushed(app, event)
            % Save in phase
            Dz = app.Diameters_UITable.Data;
            DxDy = app.DxDy_UITable.Data;
            DxDz = app.DxDz_UITable.Data;
            app.phase(app.currentphase_diameter_setup).size_histogram.along_3rd_axis = (Dz(:,1:end-1));
            app.phase(app.currentphase_diameter_setup).elongation_histogram_dx_over_dy.along_3rd_axis = DxDy(:,1:end-1);
            app.phase(app.currentphase_diameter_setup).elongation_histogram_dx_over_dz.along_3rd_axis = DxDz(:,1:end-1);
            % Save for re-load
            app.sav_diameter(app.currentphase_diameter_setup).D_data = app.Diameters_main_UITable.Data;
            app.sav_diameter(app.currentphase_diameter_setup).Dx_data = app.Diameters_UITable.Data;
            app.sav_diameter(app.currentphase_diameter_setup).Dx_ColumnName = app.Diameters_UITable.ColumnName;
            app.sav_diameter(app.currentphase_diameter_setup).DxDy_data = app.DxDy_UITable.Data;
            app.sav_diameter(app.currentphase_diameter_setup).DxDy_ColumnName = app.DxDy_UITable.ColumnName;
            app.sav_diameter(app.currentphase_diameter_setup).DxDz_data = app.DxDz_UITable.Data;
            app.sav_diameter(app.currentphase_diameter_setup).DxDz_ColumnName = app.DxDz_UITable.ColumnName;            

        end

        % Value changed function: Selectphase_orientation_DropDown
        function Selectphase_orientation_DropDownValueChanged(app, event)
            value = app.Selectphase_orientation_DropDown.Value;
            app.currentphase_elongation_setup = find(ismember(app.Selectphase_orientation_DropDown.Items,value));
            % Re-load previous state
            app.Orientation_main_UITable.Data = app.sav_elongation(app.currentphase_elongation_setup).R_data;
            app.Rotation_x_UITable.Data = app.sav_elongation(app.currentphase_elongation_setup).Rx_data;
            app.Rotation_x_UITable.ColumnName = app.sav_elongation(app.currentphase_elongation_setup).Rx_ColumnName;
            app.Rotation_y_UITable.Data = app.sav_elongation(app.currentphase_elongation_setup).Ry_data;
            app.Rotation_y_UITable.ColumnName = app.sav_elongation(app.currentphase_elongation_setup).Ry_ColumnName;
            app.Rotation_z_UITable.Data = app.sav_elongation(app.currentphase_elongation_setup).Rz_data;
            app.Rotation_z_UITable.ColumnName = app.sav_elongation(app.currentphase_elongation_setup).Rz_ColumnName;               
        end

        % Cell edit callback: Orientation_main_UITable
        function Orientation_main_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(2)==2 % Number of slices
                if isnan(newData) || newData<2 || round(newData)~=newData
                    f = msgbox('Invalid Value: must be integer >=2', 'Error','warn');
                    app.Orientation_main_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end
            elseif indices(2)==3 % Number of values
                if isnan(newData) || newData<1 || round(newData)~=newData
                    f = msgbox('Invalid Value: must be integer >=1', 'Error','warn');
                    app.Diameters_maiOrientation_main_UITablen_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end                
            end
            if indices(1)==1
                app.initialize_Ox_table;
            elseif indices(1)==2
                app.initialize_Oy_table;
            elseif indices(1)==3
                app.initialize_Oz_table;
            end     
            
        end

        % Cell edit callback: Rotation_x_UITable
        function Rotation_x_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 && indices(2)>1 % Rotation
                rot_ = app.Rotation_x_UITable.Data(1,2:end-1);
                if isnan(newData) || newData<=0 || newData>180 || sum(rot_==sort(rot_))~=length(rot_) || length(rot_)~=length(unique(rot_) )
                    f = msgbox({'Invalid Value: must be >0, <=180','Rotations must be increasing and with no repeated value'}, 'Error','warn');
                    app.Rotation_x_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)==1 % Position
                pos_ = app.Rotation_x_UITable.Data(2:end,1);
                if isnan(newData) || min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.Rotation_x_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)>1 % Percentage
                rowp = round(app.Rotation_x_UITable.Data(indices(1),2:end-1), 3); % Round it
                if isnan(newData) || newData<0 || newData>100
                    f = msgbox({'Invalid Value: must be float between 0 and 100','with sum equal to 100'}, 'Error','warn');
                    app.Rotation_x_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                else
                    number_rotation = cell2mat(app.Orientation_main_UITable.Data(1,3));
                    if number_rotation>1 && abs(sum(rowp)-100)<0.1
                        rowp(end) = 100-sum(rowp(1:end-1)); % Make sure 100%
                    end
                    app.Rotation_x_UITable.Data(indices(1),2:end-1) = rowp;
                    app.Rotation_x_UITable.Data(indices(1),end) = sum(rowp);
                    if sum(rowp)~=100
                        app.Rotation_plot.Enable = 'off';
                        app.Rotation_save.Enable = 'off';
                    else
                        app.Rotation_plot.Enable = 'on';
                        app.Rotation_save.Enable = 'on'; 
                    end
                end
            end
            app.Rotation_x_UITable.Data(1,1) = NaN;
            app.Rotation_x_UITable.Data(1,end) = NaN;            
        end

        % Cell edit callback: Rotation_y_UITable
        function Rotation_y_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 && indices(2)>1 % Rotation
                rot_ = app.Rotation_y_UITable.Data(1,2:end-1);
                if isnan(newData) || newData<=0 || newData>180 || sum(rot_==sort(rot_))~=length(rot_) || length(rot_)~=length(unique(rot_) )
                    f = msgbox({'Invalid Value: must be >0, <=180','Rotations must be increasing and with no repeated value'}, 'Error','warn');
                    app.Rotation_y_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)==1 % Position
                pos_ = app.Rotation_y_UITable.Data(2:end,1);
                if isnan(newData) || min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.Rotation_y_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)>1 % Percentage
                rowp = round(app.Rotation_y_UITable.Data(indices(1),2:end-1), 3); % Round it
                if isnan(newData) || newData<0 || newData>100
                    f = msgbox({'Invalid Value: must be float between 0 and 100','with sum equal to 100'}, 'Error','warn');
                    app.Rotation_y_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                else
                    number_rotation = cell2mat(app.Orientation_main_UITable.Data(2,3));
                    if number_rotation>1 && abs(sum(rowp)-100)<0.1
                        rowp(end) = 100-sum(rowp(1:end-1)); % Make sure 100%
                    end
                    app.Rotation_y_UITable.Data(indices(1),2:end-1) = rowp;
                    app.Rotation_y_UITable.Data(indices(1),end) = sum(rowp);
                    if sum(rowp)~=100
                        app.Rotation_plot.Enable = 'off';
                        app.Rotation_save.Enable = 'off';
                    else
                        app.Rotation_plot.Enable = 'on';
                        app.Rotation_save.Enable = 'on'; 
                    end
                end
            end
            app.Rotation_y_UITable.Data(1,1) = NaN;
            app.Rotation_y_UITable.Data(1,end) = NaN;              
        end

        % Cell edit callback: Rotation_z_UITable
        function Rotation_z_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 && indices(2)>1 % Rotation
                rot_ = app.Rotation_z_UITable.Data(1,2:end-1);
                if isnan(newData) || newData<=0 || newData>180 || sum(rot_==sort(rot_))~=length(rot_) || length(rot_)~=length(unique(rot_) )
                    f = msgbox({'Invalid Value: must be >0, <=180','Rotations must be increasing and with no repeated value'}, 'Error','warn');
                    app.Rotation_z_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)==1 % Position
                pos_ = app.Rotation_z_UITable.Data(2:end,1);
                if isnan(newData) || min(pos_)<0 || max(pos_)>1 || pos_(1)~=0 || pos_(end)~=1 || sum(pos_==sort(pos_))~=length(pos_) || length(pos_)~=length(unique(pos_))
                    f = msgbox({'Invalid Value: must be float between 0 and 1','increasing, with first value=0, last value=1','and with no repeated value'}, 'Error','warn');
                    app.Rotation_z_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                end
            elseif indices(1)>1 && indices(2)>1 % Percentage
                rowp = round(app.Rotation_z_UITable.Data(indices(1),2:end-1), 3); % Round it
                if isnan(newData) || newData<0 || newData>100
                    f = msgbox({'Invalid Value: must be float between 0 and 100','with sum equal to 100'}, 'Error','warn');
                    app.Rotation_z_UITable.Data(indices(1), indices(2)) = event.PreviousData;
                else
                    number_rotation = cell2mat(app.Orientation_main_UITable.Data(3,3));
                    if number_rotation>1 && abs(sum(rowp)-100)<0.1
                        rowp(end) = 100-sum(rowp(1:end-1)); % Make sure 100%
                    end
                    app.Rotation_z_UITable.Data(indices(1),2:end-1) = rowp;
                    app.Rotation_z_UITable.Data(indices(1),end) = sum(rowp);
                    if sum(rowp)~=100
                        app.Rotation_plot.Enable = 'off';
                        app.Rotation_save.Enable = 'off';
                    else
                        app.Rotation_plot.Enable = 'on';
                        app.Rotation_save.Enable = 'on'; 
                    end
                end
            end
            app.Rotation_z_UITable.Data(1,1) = NaN;
            app.Rotation_z_UITable.Data(1,end) = NaN;             
        end

        % Button pushed function: Rotation_plot
        function Rotation_plotButtonPushed(app, event)
            Fig = figure; % Create figure
            Fig.Name= 'Particle rotation inputs'; % Figure name
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)/2]); % Full screen figure
            Fig.Color='white'; % Background colour
            for id_axe=1:1:3
                if id_axe==1
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig);
                    hold(sub_axes,'on'); % Active subplot
                    h_title=title ('Rotation Rx (degrees)');
                    x = app.Rotation_x_UITable.Data(2:end,1);
                    y = app.Rotation_x_UITable.Data(2:end,2:end-1);
                    A1 = area(x,y);
                    [~,number_diameter] = size(y);
                    for k=1:1:number_diameter
                        set(A1(k),'DisplayName',['Rx ' num2str(k) ' = ' num2str(app.Rotation_x_UITable.Data(1,k+1))]);
                    end
                elseif id_axe==2
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig);
                    hold(sub_axes,'on'); % Active subplo
                    h_title=title ('Rotation Ry (degrees)');
                    x = app.Rotation_y_UITable.Data(2:end,1);
                    y = app.Rotation_y_UITable.Data(2:end,2:end-1);
                    A2 = area(x,y);
                    [~,number_ratio] = size(y);
                    for k=1:1:number_ratio
                        set(A2(k),'DisplayName',['Ry ' num2str(k) ' = ' num2str(app.Rotation_y_UITable.Data(1,k+1))]);
                    end
                elseif id_axe==3
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig);
                    hold(sub_axes,'on'); % Active subplo
                    h_title=title ('Rotation Rz (degrees)');
                    x = app.Rotation_z_UITable.Data(2:end,1);
                    y = app.Rotation_z_UITable.Data(2:end,2:end-1);
                    A3 = area(x,y);
                    [~,number_ratio] = size(y);
                    for k=1:1:number_ratio
                        set(A3(k),'DisplayName',['Rz ' num2str(k) ' = ' num2str(app.Rotation_z_UITable.Data(1,k+1))]);
                    end
                end
                xlabel('Normalized position along direction 3 (thickness)');
                ylabel('%');
                legend('Location','best');
                h_maintitle=sgtitle('Particle rotation inputs');
                set(sub_axes,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
                set(h_title,'FontSize',14)
                set(h_maintitle,'FontSize',16,'FontName','Times new roman')
                hold(sub_axes,'off'); % Relase axe,'off'); % Relase axe
            end
        end

        % Button pushed function: Rotation_save
        function Rotation_saveButtonPushed(app, event)
            % Save in phase
            Rx = app.Rotation_x_UITable.Data;
            Ry = app.Rotation_y_UITable.Data;
            Rz = app.Rotation_z_UITable.Data;
            app.phase(app.currentphase_elongation_setup).orientation_histogram_angledeg_x.along_3rd_axis = Rx(:,1:end-1);
            app.phase(app.currentphase_elongation_setup).orientation_histogram_angledeg_y.along_3rd_axis = Ry(:,1:end-1);
            app.phase(app.currentphase_elongation_setup).orientation_histogram_angledeg_z.along_3rd_axis = Rz (:,1:end-1);
            % Save for re-load
            app.sav_elongation(app.currentphase_elongation_setup).R_data = app.Orientation_main_UITable.Data;
            app.sav_elongation(app.currentphase_elongation_setup).Rx_data = app.Rotation_x_UITable.Data;
            app.sav_elongation(app.currentphase_elongation_setup).Rx_ColumnName = app.Rotation_x_UITable.ColumnName;
            app.sav_elongation(app.currentphase_elongation_setup).Ry_data = app.Rotation_y_UITable.Data;
            app.sav_elongation(app.currentphase_elongation_setup).Ry_ColumnName = app.Rotation_y_UITable.ColumnName;
            app.sav_elongation(app.currentphase_elongation_setup).Rz_data = app.Rotation_z_UITable.Data;
            app.sav_elongation(app.currentphase_elongation_setup).Rz_ColumnName = app.Rotation_z_UITable.ColumnName;             
        end

        % Cell edit callback: Parameter_example_UITable
        function Parameter_example_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1 % Diameter
                if isnan(newData) || newData<1 || round(newData)~=newData % Incorrect input
                    f = msgbox('Invalid Value: must be integer >=1', 'Error','warn');    
                    app.Parameter_example_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end
            elseif indices(1)==2 || indices(1)==3 % Ratio
                if isnan(newData) || newData<=0 % Incorrect input
                    f = msgbox('Invalid Value: must be >0', 'Error','warn');
                    app.Parameter_example_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end
            else % Rotation
                if isnan(newData) || newData<=0 || newData>180 % Incorrect input
                    f = msgbox('Invalid Value: must be >0 <=180', 'Error','warn');
                    app.Parameter_example_UITable.Data(indices(1), indices(2)) = num2cell(event.PreviousData);
                end                
            end           
        end

        % Button pushed function: Plot_example_particle
        function Plot_example_particleButtonPushed(app, event)
            % Get inputs
            Dx = cell2mat(app.Parameter_example_UITable.Data(1,2));
            DxDy = cell2mat(app.Parameter_example_UITable.Data(2,2));
            DxDz = cell2mat(app.Parameter_example_UITable.Data(3,2));
            Rx = cell2mat(app.Parameter_example_UITable.Data(4,2));
            Ry = cell2mat(app.Parameter_example_UITable.Data(5,2));
            Rz = cell2mat(app.Parameter_example_UITable.Data(6,2));
            % Generate ellipsoids
            Dy = Dx/DxDy;
            Dz = Dx/DxDz;
            [binary_ellipsoid] = create_ellipsoid(Dx,Dy,Dz);
            % Rotation
            if Rx~=180 || Ry~=180 || Rz~=180
                [binary_ellipsoid] = rotate_domain(binary_ellipsoid,Rx, Ry, Rz);
            end
            % Patch
            [f,v] = function_patch_3Darray(binary_ellipsoid);
            
            Fig = figure; % Create figure
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)*2/3]); % Full screen figure            
            Fig.Name= 'Particle example'; % Figure name
            Fig.Color='white'; % Background colour
            ax_=axes('parent',Fig);
            h_title = title('Particle example');
            hold(ax_,'on'); % Active subplot
            %linelength = (1/3)*min([Dx,Dy,Dz]);
            %plot3([0 linelength],[0 0],[0 0],'DisplayName','Axe 1','Parent',ax_);
            %plot3([0 0],[0 linelength],[0 0],'DisplayName','Axe 2','Parent',ax_);
            %plot3([0 0],[0 0],[0 linelength],'DisplayName','Axe 3','Parent',ax_);
            col=v(:,3); % Z-axis is color altitude
            patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','interp','DisplayName','Particle','Parent',ax_);
            h=colorbar(ax_);
            ylabel(h, 'z position');
            xlabel('x - 1st axis')
            ylabel('y - 2nd axis')
            zlabel('z - 3rd axis (thickness)')
            axis(ax_,'equal'); axis(ax_,'tight');
            %legend(ax_,{'x - 1st axis';'y - 2nd axis';'z - 3rd axis (thickness)';'Particle'},'Location','best');
            view(ax_,3)
            set(ax_,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
            set(h,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
            h_title.FontSize = 14; % Set title fontsize   
            hold(ax_,'off'); % Relase axe
        end

        % Button pushed function: Plot_example_particle_2
        function Plot_example_particle_2ButtonPushed(app, event)
            % Get inputs
            Dx = cell2mat(app.Parameter_example_UITable.Data(1,2));
            DxDy = cell2mat(app.Parameter_example_UITable.Data(2,2));
            DxDz = cell2mat(app.Parameter_example_UITable.Data(3,2));
            Rx = cell2mat(app.Parameter_example_UITable.Data(4,2));
            Ry = cell2mat(app.Parameter_example_UITable.Data(5,2));
            Rz = cell2mat(app.Parameter_example_UITable.Data(6,2));
            % Generate ellipsoids
            Dy = Dx/DxDy;
            Dz = Dx/DxDz;
            [binary_ellipsoid] = create_ellipsoid(Dx,Dy,Dz);
            % Rotation
            if Rx~=180 || Ry~=180 || Rz~=180
                [binary_ellipsoid] = rotate_domain(binary_ellipsoid,Rx, Ry, Rz);
            end
            % Plot
            Fig_3D=figure;
            volshow(binary_ellipsoid,'Parent', Fig_3D,'BackgroundColor','w','Colormap',copper);            
        end

        % Cell edit callback: overlapping_UITable
        function overlapping_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if isnan(newData) || newData<0 % Incorrect input
                f = msgbox('Invalid Value: must be >0', 'Error','warn');
                app.overlapping_UITable.Data(indices(1), indices(2)) = event.PreviousData;
            else
                if indices(1)~=indices(2) % Enfore symmetry
                    app.overlapping_UITable.Data(indices(2), indices(1)) = app.overlapping_UITable.Data(indices(1), indices(2));
                end
            end
        end

        % Cell edit callback: Minimumvolume_overlapping_UITable
        function Minimumvolume_overlapping_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if isnan(newData) || newData<0 || newData>1  % Incorrect input
                f = msgbox('Invalid Value: must be [0,1]', 'Error','warn');
                app.Minimumvolume_overlapping_UITable.Data(indices(1), indices(2)) = event.PreviousData;
            end            
            
        end

        % Button pushed function: Overlapping_save
        function Overlapping_saveButtonPushed(app, event)
            app.Maximum_overlapping = app.overlapping_UITable.Data;
            app.Minimum_particle_volume_conservated = app.Minimumvolume_overlapping_UITable.Data;    
            app.disable_contiguituy_check = app.DisableparticlecontiguitycheckCheckBox.Value;
            % Move to next tab
            app.TabGroup.SelectedTab = app.OrderTab;   
        end

        % Button pushed function: Order_save
        function Order_saveButtonPushed(app, event)
            % Get data
            fillratio = app.Order_UITable.Data;
            [~,number_phase] = size(fillratio);
            % Save
            for k = 1:1:number_phase
                app.phase(k).fillratio = fillratio(:,k);
            end
            % Move to next tab
            app.TabGroup.SelectedTab = app.StopconditionsTab;  
            
        end

        % Value changed function: NumberofgenerationpassEditField
        function NumberofgenerationpassEditFieldValueChanged(app, event)
              app.Initialize_order_table;   
        end

        % Cell edit callback: Order_UITable
        function Order_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            column_ = app.Order_UITable.Data(:,indices(2));
            if isnan(newData) || newData<0 || newData>1 || column_(end)<1 || sum(column_==sort(column_))~=length(column_)
                f = msgbox({'Invalid Value: must be [0,1]','and increasing (or stable) from pass i to pass i+1','Fill ratio for last pass must be 1'}, 'Error','warn');
                app.Order_UITable.Data(indices(1), indices(2)) = event.PreviousData;
            end
        end

        % Button pushed function: Generate_Button
        function Generate_ButtonPushed(app, event)
            % Clear table
            app.Volume_fraction_UITable.Data = [];
            app.outcometime_UITable.Data = [];
            app.Connectivity_UITable.Data = [];
            app.Tortuosity_UITable.Data = [];
            app.Visualization_UITable.Data = [];
            
            % Stop condition options
            stoping_conditions.plot = app.PlotalgorithmprogressionCheckBox.Value;
            stoping_conditions.vfrate_threshold = app.Ifvolumefractions1goesbelowEditField.Value;
            stoping_conditions.andor = app.stopcondition_andor_DropDown.Value;
            stoping_conditions.particlerate_threshold = app.particles1goesbelowEditField.Value;
            stoping_conditions.average_on_last_n_iterations = app.whenaveragedonthelastEditField.Value;
            stoping_conditions.action = app.stopcondition_action_DropDown.Value;
            stoping_conditions.maxtime = app.MaximumwallclocktimesEditField.Value;
                        
            % Verification and characterization options
            do_verification = app.CompareoutputsversusinputsCheckBox.Value;
            postprocessing.volumefractions = app.VolumefractionsCheckBox.Value;
            postprocessing.connectivity = app.ConnectivityslowforlargevolumesCheckBox.Value;
            postprocessing.tortuosity = app.PoretortuosityfactorslowforlargevolumesCheckBox.Value;
            postprocessing.perform_on = app.CharacterizeDropDown.Value;
            postprocessing.save = app.SavecharacterizationresultsifselectedCheckBox.Value;
            
            % Scaling options
            scaling_factor = 1/app.Scalingfactor1upscaling1downscaling1noscalingEditField.Value;
            scale_phaselabel = app.ApplyscalingtophaselabelfastCheckBox.Value;
            scale_particlelabel = app.ApplyscalingalsotoparticlelabelslowCheckBox.Value;
            
            % Visualization tab
            app.Plot_phase_label.Enable = 'off';
            app.Plot_particle_label.Enable = 'off';      
            
            % Save options (within algorithm)
            save_options.folder = app.Savefolder;
            save_options.save_progression = app.SavealgorithmprogressionplotifselectedCheckBox.Value;
            save_options.save_verification = app.SaveinputsoutpuscomparisonifselectedCheckBox.Value;
            
            % Save inputs
            save_input_mat = app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox.Value;
 
            % Save outputs
            save_phaselabel = app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox.Value;
            save_particlelabel = app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox.Value;
            save_original_andor_scaled = app.ForphaseandparticlelablelsaveforDropDown.Value;
            save_additionalinfo = app.SaveparticleadditionalinformationCheckBox.Value;
            
            number_run = app.NumberofrunsEditField.Value;
            wallclocktime = zeros(number_run,3);
            for k_run=1:1:number_run
                save_options.run_number = k_run;

                tic
                % Call generation algorithm
                [microstructure3D, current_phaseinfo, outcome] = function_generate_ellipsoid_microstructure(app.domain_size,app.phase,app.Maximum_overlapping,app.Minimum_particle_volume_conservated,~app.disable_contiguituy_check, stoping_conditions, do_verification, save_options);
                wallclocktime(k_run,1) = toc;
                % Save for later visualization
                app.generation_result(k_run).microstructure3D_phaselabel = microstructure3D.phase;
                app.generation_result(k_run).microstructure3D_particlelabel = microstructure3D.particle_id;

                % Upscaling
                if scaling_factor~=1
                    tic
                    % Set parameters
                    parameters_scaling.scaling_factor = scaling_factor;
                    parameters_scaling.label_or_greylevel = 'Label';
                    parameters_scaling.background = 0;
                    % Scale
                    if scale_phaselabel
                        app.generation_result(k_run).microstructure3D_phaselabel_scaled = function_scaling(app.generation_result(k_run).microstructure3D_phaselabel,parameters_scaling);
                    end
                    if scale_particlelabel
                        app.generation_result(k_run).microstructure3D_particlelabel_scaled = function_scaling(app.generation_result(k_run).microstructure3D_particlelabel,parameters_scaling);
                    end
                    wallclocktime(k_run,2) = toc;
                end

                % % Characterization
                tic
                
                % Volume fractions
                if postprocessing.volumefractions
                    vf = zeros(1,app.number_solidphase+1);
                    if scaling_factor~=1 && strcmp(postprocessing.perform_on,'volume after upscaling')
                        vf(1) = sum(sum(sum( app.generation_result(k_run).microstructure3D_phaselabel_scaled==0 )));
                        for k_phase = 1:1:app.number_solidphase
                            vf(k_phase+1) = sum(sum(sum( app.generation_result(k_run).microstructure3D_phaselabel_scaled == app.phase(k_phase).code )));
                        end
                        vf = vf/numel(app.generation_result(k_run).microstructure3D_phaselabel_scaled);
                    else
                        vf(1) = sum(sum(sum( app.generation_result(k_run).microstructure3D_phaselabel==0 )));
                        for k_phase = 1:1:app.number_solidphase
                            vf(k_phase+1) = sum(sum(sum( app.generation_result(k_run).microstructure3D_phaselabel == app.phase(k_phase).code )));
                        end
                        vf = vf/numel(app.generation_result(k_run).microstructure3D_phaselabel);
                    end
                    if k_run==1
                        app.Volume_fraction_UITable.Data = [k_run, vf];
                    else
                        app.Volume_fraction_UITable.Data(k_run,:) = [k_run, vf];
                    end
                end
                
                % Connectivity
                if postprocessing.connectivity
                    conn = zeros(1,app.number_solidphase+2);
                    if scaling_factor~=1 && strcmp(postprocessing.perform_on,'volume after upscaling')
                        conn(1) = app.Connectivity_mainclusterpercent(app.generation_result(k_run).microstructure3D_phaselabel_scaled,0);
                        for k_phase = 1:1:app.number_solidphase
                            conn(k_phase+1) = app.Connectivity_mainclusterpercent(app.generation_result(k_run).microstructure3D_phaselabel_scaled,app.phase(k_phase).code);
                        end
                        tmp = zeros(size(app.generation_result(k_run).microstructure3D_phaselabel_scaled));
                        tmp(app.generation_result(k_run).microstructure3D_phaselabel_scaled~=0)=1;
                        conn(end) = app.Connectivity_mainclusterpercent(tmp,1);
                    else
                        conn(1) = app.Connectivity_mainclusterpercent(app.generation_result(k_run).microstructure3D_phaselabel,0);
                        for k_phase = 1:1:app.number_solidphase
                            conn(k_phase+1) = app.Connectivity_mainclusterpercent(app.generation_result(k_run).microstructure3D_phaselabel,app.phase(k_phase).code);
                        end       
                        tmp = zeros(size(app.generation_result(k_run).microstructure3D_phaselabel));
                        tmp(app.generation_result(k_run).microstructure3D_phaselabel~=0)=1;
                        conn(end) = app.Connectivity_mainclusterpercent(tmp,1);                        
                    end
                    if k_run==1
                        app.Connectivity_UITable.Data = [k_run, conn];
                    else
                        app.Connectivity_UITable.Data(k_run,:) = [k_run, conn];
                    end
                end
                
                % Pore tortuosity factor
                if postprocessing.tortuosity
                    tau = zeros(1,3);
                    if scaling_factor~=1 && strcmp(postprocessing.perform_on,'volume after upscaling')
                        binary_phase = zeros(size(app.generation_result(k_run).microstructure3D_phaselabel_scaled));
                        binary_phase(app.generation_result(k_run).microstructure3D_phaselabel_scaled==0)=1;
                    else
                        binary_phase = zeros(size(app.generation_result(k_run).microstructure3D_phaselabel));
                        binary_phase(app.generation_result(k_run).microstructure3D_phaselabel==0)=1;
                    end
                    tau(1) = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]).Tau_W1.Tau;
                    tau(2) = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]).Tau_W2.Tau;
                    tau(3) = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]).Tau_W3.Tau;
                    if k_run==1
                        app.Tortuosity_UITable.Data = [k_run, tau];
                    else
                        app.Tortuosity_UITable.Data(k_run,:) = [k_run, tau];
                    end
                end
                
                % End of characterization
                wallclocktime(k_run,3) = toc;
                
                % Outcome and wallclock time
                if k_run==1
                    app.outcometime_UITable.Data = [outcome, {num2str(wallclocktime(k_run,1),'%1.1f')}, {num2str(wallclocktime(k_run,2),'%1.1f')}, {num2str(wallclocktime(k_run,3),'%1.1f')}];
                else
                    app.outcometime_UITable.Data(k_run,:) = [outcome, {num2str(wallclocktime(k_run,1),'%1.1f')}, {num2str(wallclocktime(k_run,2),'%1.1f')}, {num2str(wallclocktime(k_run,3),'%1.1f')}];
                end                  
                
                % Save phase and particle label
                if ~isempty(app.Savefolder)
                    if save_phaselabel
                        if strcmp(save_original_andor_scaled,'volume before upscaling') || strcmp(save_original_andor_scaled,'both')
                            function_save_tif( uint8(app.generation_result(k_run).microstructure3D_phaselabel), [app.Savefolder 'Phaselabel_run_' num2str(k_run) '.tif']);
                        end
                        if scaling_factor~=1 && (strcmp(save_original_andor_scaled,'volume after upscaling') || strcmp(save_original_andor_scaled,'both'))
                            function_save_tif( uint8(app.generation_result(k_run).microstructure3D_phaselabel_scaled), [app.Savefolder 'Phaselabel_scaled_run_' num2str(k_run) '.tif']);
                        end
                    end
                    if save_particlelabel
                        if strcmp(save_original_andor_scaled,'volume before upscaling') || strcmp(save_original_andor_scaled,'both')
                            function_save_tif( uint16(app.generation_result(k_run).microstructure3D_particlelabel), [app.Savefolder 'Particlelabel_run_' num2str(k_run) '.tif']);
                        end
                        if scaling_factor~=1 && scale_particlelabel && (strcmp(save_original_andor_scaled,'volume after upscaling') || strcmp(save_original_andor_scaled,'both'))
                            function_save_tif( uint16(app.generation_result(k_run).microstructure3D_particlelabel_scaled), [app.Savefolder 'Particlelabel_scaled_run_' num2str(k_run) '.tif']);
                        end
                    end                    
                end
                
                % Save additional info
                if ~isempty(app.Savefolder)
                    if save_additionalinfo
                        save([app.Savefolder 'Additionalinfo_run_' num2str(k_run) '.mat'],'microstructure3D','-mat');
                    end
                end
                
                % % Visualization tab
                if scaling_factor==1
                    if k_run==1
                        app.Visualization_UITable.Data = [k_run, false, true];
                    else
                        app.Visualization_UITable.Data(k_run,:) = [k_run, false, true];
                    end 
                else
                    if k_run==1
                        app.Visualization_UITable.Data = [k_run, false, true];
                        app.Visualization_UITable.Data(2,:) = [k_run, true, true];
                    else
                        app.Visualization_UITable.Data(2*k_run-1,:) = [k_run, false, true];
                        app.Visualization_UITable.Data(2*k_run,:) = [k_run, true, true];
                    end                     
                end
                         
            end
            
            % Save characterization
            if ~isempty(app.Savefolder) && postprocessing.save
                filename = 'Characterization'; % Filename without extension
                % Prepare the data
                clear DATA_writetable

                if postprocessing.volumefractions
                    t = app.Build_table_fromUItab(app.Volume_fraction_UITable);
                    DATA_writetable.sheet(1).name='Volume_fractions';
                    DATA_writetable.sheet(1).table=t;
                end

                if postprocessing.connectivity
                    t = app.Build_table_fromUItab(app.Connectivity_UITable);
                    DATA_writetable.sheet(2).name='Connectivity';
                    DATA_writetable.sheet(2).table=t;
                end

                if postprocessing.tortuosity
                    t = app.Build_table_fromUItab(app.Tortuosity_UITable);
                    DATA_writetable.sheet(3).name='Tortuosity_factor';
                    DATA_writetable.sheet(3).table=t;
                end
                
                t = app.Build_table_fromUItab(app.outcometime_UITable);
                DATA_writetable.sheet(4).name='Outcome_time';
                DATA_writetable.sheet(4).table=t;                

                % Save function
                Function_Writetable(app.Savefolder,filename,DATA_writetable)
            end

            % Save inputs
            if ~isempty(app.Savefolder)
                if save_input_mat
                    domain_size = app.domain_size;
                    phase = app.phase;
                    Maximum_overlapping = app.Maximum_overlapping;
                    Minimum_particle_volume_conservated = app.Minimum_particle_volume_conservated;
                    check_contiguity = ~app.disable_contiguituy_check;
                    stopingconditions = stoping_conditions;
                    stopingconditions.plot = true; % Overwritte
                    doverification = true; % Overwritte
                    saveoptions.folder = []; % Overwritte
                    save([app.Savefolder 'Inputs.mat'],'domain_size','phase','Maximum_overlapping','Minimum_particle_volume_conservated','check_contiguity','stopingconditions', 'doverification', 'saveoptions','-mat');
                end
            end
            
            % Visualization tab
            app.Plot_phase_label.Enable = 'on';
            app.Plot_particle_label.Enable = 'on';
            if number_run==1
                app.Plot_phase_label_3D.Enable = 'on';
                app.Plot_particle_label_3D.Enable = 'on';
            else
                app.Plot_phase_label_3D.Enable = 'off';
                app.Plot_particle_label_3D.Enable = 'off';    
            end
                        
        end

        % Button pushed function: Doc_button
        function Doc_buttonButtonPushed(app, event)
            Find_file('NREL_MATBOX_Microstructure_analysis_toolbox_documentation.pdf','MATBOX_Microstructure_analysis_toolbox','Default location is \MATBOX_Microstructure_analysis_toolbox\Documentation\');                        
        end

        % Button pushed function: Github_button
        function Github_buttonButtonPushed(app, event)
            url = 'https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/';
            web(url)
        end

        % Image clicked function: NREL_logo
        function NREL_logoClicked(app, event)
            url = 'https://www.nrel.gov/transportation/energy-storage.html/';
            web(url)            
        end

        % Button pushed function: Plot_phase_label
        function Plot_phase_labelButtonPushed(app, event)
            runs = app.Visualization_UITable.Data(:,1);
            scaled = app.Visualization_UITable.Data(:,2);
            choices = app.Visualization_UITable.Data(:,3);
            idx = find(choices==1);
            if ~isempty(idx)
                if length(idx)==1 % One volume to visualize
                    if scaled(idx)
                        microstructure_visualization_slices(app.generation_result(runs(idx)).microstructure3D_phaselabel_scaled, ['Generated volume #' num2str(runs(idx))], 1, 'voxel', [])
                    else
                        microstructure_visualization_slices(app.generation_result(runs(idx)).microstructure3D_phaselabel, ['Generated volume #' num2str(runs(idx)) ', scaled'], 1, 'voxel', [])
                    end
                else % Several volumes to visualize
                    str = 'Microstructure_comparison_visualization_interface(';
                    for k=1:1:length(idx)
                        if scaled(idx(k))
                            str = [str 'app.generation_result(runs(idx(' num2str(k) '))).microstructure3D_phaselabel_scaled,'];
                        else
                            str = [str 'app.generation_result(runs(idx(' num2str(k) '))).microstructure3D_phaselabel,'];
                        end
                    end
                    str(end)=[]; % Remove last ,
                    str = [str ');'];
                    eval(str);
                end
            end            
        end

        % Button pushed function: Plot_particle_label
        function Plot_particle_labelButtonPushed(app, event)
            runs = app.Visualization_UITable.Data(:,1);
            scaled = app.Visualization_UITable.Data(:,2);
            choices = app.Visualization_UITable.Data(:,3);
            idx = find(choices==1);
            if ~isempty(idx)
                if length(idx)==1 % One volume to visualize
                    if scaled(idx)
                        microstructure_visualization_slices(app.generation_result(runs(idx)).microstructure3D_particlelabel_scaled, ['Generated volume #' num2str(runs(idx))], 1, 'voxel', [])
                    else
                        microstructure_visualization_slices(app.generation_result(runs(idx)).microstructure3D_particlelabel, ['Generated volume #' num2str(runs(idx)) ', scaled'], 1, 'voxel', [])
                    end
                else % Several volumes to visualize
                    str = 'Microstructure_comparison_visualization_interface(';
                    for k=1:1:length(idx)
                        if scaled(idx(k))
                            str = [str 'app.generation_result(runs(idx(' num2str(k) '))).microstructure3D_particlelabel_scaled,'];
                        else
                            str = [str 'app.generation_result(runs(idx(' num2str(k) '))).microstructure3D_particlelabel,'];
                        end
                    end
                    str(end)=[]; % Remove last ,
                    str = [str ');'];
                    eval(str);
                end
            end
        end

        % Button pushed function: Save_folder_button
        function Save_folder_buttonButtonPushed(app, event)
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
                app.NosavefolderselectedLabel.Text = app.Savefolder;
                statut = 'on';
            else
                app.Savefolder = [];
                app.NosavefolderselectedLabel.Text = 'No save folder selected';
                statut = 'off';
            end
                app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox.Enable = statut;
                app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox.Enable = statut;
                app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox.Enable = statut;
                app.ForphaseandparticlelablelsaveforDropDown.Enable = statut;
                app.ForphaseandparticlelablelsaveforDropDownLabel.Enable = statut;
                app.SaveparticleadditionalinformationCheckBox.Enable = statut;
                app.SavealgorithmprogressionplotifselectedCheckBox.Enable = statut;
                app.SavecharacterizationresultsifselectedCheckBox.Enable = statut;           
                app.SaveinputsoutpuscomparisonifselectedCheckBox.Enable = statut;
            
        end

        % Cell edit callback: Visualization_UITable
        function Visualization_UITableCellEdit(app, event)
            choices = app.Visualization_UITable.Data(:,3);
            idx = find(choices==1);
            if length(idx)==1
                app.Plot_phase_label_3D.Enable = 'on';
                app.Plot_particle_label_3D.Enable = 'on';
            else
                app.Plot_phase_label_3D.Enable = 'off';
                app.Plot_particle_label_3D.Enable = 'off';    
            end
            if ~isempty(idx)
                app.Plot_phase_label.Enable = 'on';
                app.Plot_particle_label.Enable = 'on';                
            else
                app.Plot_phase_label.Enable = 'off';
                app.Plot_particle_label.Enable = 'off';
            end
        end

        % Button pushed function: Plot_phase_label_3D
        function Plot_phase_label_3DButtonPushed(app, event)
            runs = app.Visualization_UITable.Data(:,1);
            scaled = app.Visualization_UITable.Data(:,2);
            choices = app.Visualization_UITable.Data(:,3);
            idx = find(choices==1);
            Fig_3D=figure;
            col_ = eval(app.DviewcolormapDropDown.Value);
            if scaled(idx)
                volshow(app.generation_result(runs(idx)).microstructure3D_phaselabel_scaled,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');
            else
                volshow(app.generation_result(runs(idx)).microstructure3D_phaselabel,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');
            end
        end

        % Button pushed function: Plot_particle_label_3D
        function Plot_particle_label_3DButtonPushed(app, event)
            runs = app.Visualization_UITable.Data(:,1);
            scaled = app.Visualization_UITable.Data(:,2);
            choices = app.Visualization_UITable.Data(:,3);
            idx = find(choices==1);
            Fig_3D=figure;
            col_ = eval(app.DviewcolormapDropDown.Value);
            if scaled(idx)
                volshow(app.generation_result(runs(idx)).microstructure3D_particlelabel_scaled,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');
            else
                volshow(app.generation_result(runs(idx)).microstructure3D_particlelabel,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');
            end            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create EllipsoidbasedstochasticgenerationmoduleUIFigure and hide until all components are created
            app.EllipsoidbasedstochasticgenerationmoduleUIFigure = uifigure('Visible', 'off');
            app.EllipsoidbasedstochasticgenerationmoduleUIFigure.Position = [100 100 1110 775];
            app.EllipsoidbasedstochasticgenerationmoduleUIFigure.Name = 'Ellipsoid-based stochastic generation module';
            app.EllipsoidbasedstochasticgenerationmoduleUIFigure.Icon = 'Icon_generation.png';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.EllipsoidbasedstochasticgenerationmoduleUIFigure);
            app.TabGroup.Tooltip = {'Run the generation algorithm'};
            app.TabGroup.TabLocation = 'left';
            app.TabGroup.Position = [1 -3 1110 779];

            % Create Instructions
            app.Instructions = uitab(app.TabGroup);
            app.Instructions.Tooltip = {'Instructions'};
            app.Instructions.Title = 'Instructions';

            % Create Instructions_title
            app.Instructions_title = uilabel(app.Instructions);
            app.Instructions_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_title.HorizontalAlignment = 'center';
            app.Instructions_title.FontWeight = 'bold';
            app.Instructions_title.Position = [11 746 977 22];
            app.Instructions_title.Text = 'Main instructions';

            % Create SetpropertiesdefinedwithascalarLabel
            app.SetpropertiesdefinedwithascalarLabel = uilabel(app.Instructions);
            app.SetpropertiesdefinedwithascalarLabel.FontWeight = 'bold';
            app.SetpropertiesdefinedwithascalarLabel.FontColor = [0.2902 0.451 0.0863];
            app.SetpropertiesdefinedwithascalarLabel.Position = [11 349 219 22];
            app.SetpropertiesdefinedwithascalarLabel.Text = '   Set properties defined with a scalar';

            % Create SetpropertiesdefinedwithanhistogramdistributionLabel
            app.SetpropertiesdefinedwithanhistogramdistributionLabel = uilabel(app.Instructions);
            app.SetpropertiesdefinedwithanhistogramdistributionLabel.FontWeight = 'bold';
            app.SetpropertiesdefinedwithanhistogramdistributionLabel.FontColor = [0 0 1];
            app.SetpropertiesdefinedwithanhistogramdistributionLabel.Position = [11 327 331 22];
            app.SetpropertiesdefinedwithanhistogramdistributionLabel.Text = '   Set properties defined with an histogram (distribution)';

            % Create SetalgorithmadditionalparametersLabel
            app.SetalgorithmadditionalparametersLabel = uilabel(app.Instructions);
            app.SetalgorithmadditionalparametersLabel.FontWeight = 'bold';
            app.SetalgorithmadditionalparametersLabel.FontColor = [1 0.302 0];
            app.SetalgorithmadditionalparametersLabel.Position = [11 305 221 22];
            app.SetalgorithmadditionalparametersLabel.Text = '   Set algorithm additional parameters';

            % Create StepsperformedafterthegenerationprogressLabel
            app.StepsperformedafterthegenerationprogressLabel = uilabel(app.Instructions);
            app.StepsperformedafterthegenerationprogressLabel.FontWeight = 'bold';
            app.StepsperformedafterthegenerationprogressLabel.FontColor = [0.4941 0.1843 0.5569];
            app.StepsperformedafterthegenerationprogressLabel.Position = [11 283 662 22];
            app.StepsperformedafterthegenerationprogressLabel.Text = '   Set parameters for tasks performed after the generation algorithm is finished, but that we configure in advance';

            % Create RunthegenerationalgorithmLabel
            app.RunthegenerationalgorithmLabel = uilabel(app.Instructions);
            app.RunthegenerationalgorithmLabel.FontWeight = 'bold';
            app.RunthegenerationalgorithmLabel.FontColor = [1 0 0];
            app.RunthegenerationalgorithmLabel.Position = [11 261 183 22];
            app.RunthegenerationalgorithmLabel.Text = '   Run the generation algorithm';

            % Create SelectasolidphaseandsetthenumberofLabel_14
            app.SelectasolidphaseandsetthenumberofLabel_14 = uilabel(app.Instructions);
            app.SelectasolidphaseandsetthenumberofLabel_14.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_14.Position = [11 709 976 22];
            app.SelectasolidphaseandsetthenumberofLabel_14.Text = {'Stochastic (i.e., based on randomness) microstructure generation algorithm. Main features are listed below:'; ''};

            % Create SelectasolidphaseandsetthenumberofLabel_15
            app.SelectasolidphaseandsetthenumberofLabel_15 = uilabel(app.Instructions);
            app.SelectasolidphaseandsetthenumberofLabel_15.FontAngle = 'italic';
            app.SelectasolidphaseandsetthenumberofLabel_15.Position = [11 615 976 98];
            app.SelectasolidphaseandsetthenumberofLabel_15.Text = {'- n phase microstructure'; '- Ellipsoid-based, with particle overlapping control (so that you can overcome particle packing density limit)'; '- Control of volume fractions, particle diameters, particle elongation, and particle orientation/rotation.'; '  Properties are defined along the volume thickness, so that you can generate graded microstructures.'; '  Particle size and morphology aredefined with histograms, so that you can generate in-plane heterogeneities.'; '- Particle generation order (see dedicated tab for more details).'};

            % Create SelectasolidphaseandsetthenumberofLabel_16
            app.SelectasolidphaseandsetthenumberofLabel_16 = uilabel(app.Instructions);
            app.SelectasolidphaseandsetthenumberofLabel_16.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_16.Position = [11 578 976 22];
            app.SelectasolidphaseandsetthenumberofLabel_16.Text = 'You can see examples in the Git hub repository using the link below:';

            % Create Hyperlink
            app.Hyperlink = uihyperlink(app.Instructions);
            app.Hyperlink.URL = 'https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/tree/master/Data_example';
            app.Hyperlink.Position = [11 556 750 22];
            app.Hyperlink.Text = 'Example of numerically generated microstructures (scroll to the bottom to see pictures, tif files are also availabe in a subfolder)';

            % Create SelectasolidphaseandsetthenumberofLabel_17
            app.SelectasolidphaseandsetthenumberofLabel_17 = uilabel(app.Instructions);
            app.SelectasolidphaseandsetthenumberofLabel_17.FontAngle = 'italic';
            app.SelectasolidphaseandsetthenumberofLabel_17.Position = [11 491 994 65];
            app.SelectasolidphaseandsetthenumberofLabel_17.Text = {'As-generated microstructures look obviously numerically generated. The aim of the algorithm is to generate microstructures with very well defined properties, not to generate volumes'; 'that look like obtained from imaging, i.e., prioritize control over fidelity.'; 'Combined with the MATBOX characterization and meshing module, it enables establishing correlations between microstructure parameters (e.g., effective diffusion coefficient as'; 'function of particle elongation) or between microstructure parameters and battery performances (using the mesh in a microscale electrochemical model).'};

            % Create SelectasolidphaseandsetthenumberofLabel_18
            app.SelectasolidphaseandsetthenumberofLabel_18 = uilabel(app.Instructions);
            app.SelectasolidphaseandsetthenumberofLabel_18.HorizontalAlignment = 'center';
            app.SelectasolidphaseandsetthenumberofLabel_18.FontSize = 14;
            app.SelectasolidphaseandsetthenumberofLabel_18.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_18.Position = [85 444 830 32];
            app.SelectasolidphaseandsetthenumberofLabel_18.Text = {'Therefore, the intented original use of this algorithm is to generate a variety of microstructures for design space analysis'; 'to identify promising set of parameters for a given application.'};

            % Create SelectasolidphaseandsetthenumberofLabel_19
            app.SelectasolidphaseandsetthenumberofLabel_19 = uilabel(app.Instructions);
            app.SelectasolidphaseandsetthenumberofLabel_19.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_19.Position = [11 376 976 28];
            app.SelectasolidphaseandsetthenumberofLabel_19.Text = {'How to use it: simply follow instructions for each tab, from top to bottom.'; 'Each tab is color-coded as indicated below:'};

            % Create Doc_button
            app.Doc_button = uibutton(app.Instructions, 'push');
            app.Doc_button.ButtonPushedFcn = createCallbackFcn(app, @Doc_buttonButtonPushed, true);
            app.Doc_button.BackgroundColor = [0.0745 0.6235 1];
            app.Doc_button.FontSize = 16;
            app.Doc_button.FontWeight = 'bold';
            app.Doc_button.Position = [11 26 211 40];
            app.Doc_button.Text = 'Open documentation';

            % Create Github_button
            app.Github_button = uibutton(app.Instructions, 'push');
            app.Github_button.ButtonPushedFcn = createCallbackFcn(app, @Github_buttonButtonPushed, true);
            app.Github_button.BackgroundColor = [0.0745 0.6235 1];
            app.Github_button.FontSize = 16;
            app.Github_button.FontWeight = 'bold';
            app.Github_button.Position = [242 26 211 40];
            app.Github_button.Text = 'Check Git hub repository';

            % Create NREL_logo
            app.NREL_logo = uiimage(app.Instructions);
            app.NREL_logo.ImageClickedFcn = createCallbackFcn(app, @NREL_logoClicked, true);
            app.NREL_logo.Position = [688 26 299 121];
            app.NREL_logo.ImageSource = 'logo_NREL.png';

            % Create SaveTab
            app.SaveTab = uitab(app.TabGroup);
            app.SaveTab.Title = 'Save';

            % Create Instructions_title_3
            app.Instructions_title_3 = uilabel(app.SaveTab);
            app.Instructions_title_3.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_title_3.HorizontalAlignment = 'center';
            app.Instructions_title_3.FontWeight = 'bold';
            app.Instructions_title_3.Position = [11 746 977 22];
            app.Instructions_title_3.Text = 'Saving options';

            % Create Savephaselabelinformation3Dtifstackfileuint8formatCheckBox
            app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox = uicheckbox(app.SaveTab);
            app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox.Enable = 'off';
            app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox.Text = 'Save phase label information (3D tif stack file, uint8 format)';
            app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox.Position = [11 489 341 22];
            app.Savephaselabelinformation3Dtifstackfileuint8formatCheckBox.Value = true;

            % Create Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox
            app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox = uicheckbox(app.SaveTab);
            app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox.Enable = 'off';
            app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox.Text = 'Save particle label information (3D tif stack file, uint16 format)';
            app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox.Position = [11 462 354 22];
            app.Saveparticlelabelinformation3Dtifstackfileuint16formatCheckBox.Value = true;

            % Create SaveinmatinputsofthegenerationalgorithmfunctionCheckBox
            app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox = uicheckbox(app.SaveTab);
            app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox.Enable = 'off';
            app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox.Text = 'Save in .mat (inputs of the generation algorithm function, so that you can re-run the algorithm in command line if needed)';
            app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox.Position = [11 553 678 22];
            app.SaveinmatinputsofthegenerationalgorithmfunctionCheckBox.Value = true;

            % Create InputsLabel
            app.InputsLabel = uilabel(app.SaveTab);
            app.InputsLabel.Position = [11 580 42 22];
            app.InputsLabel.Text = 'Inputs:';

            % Create OutputsLabel
            app.OutputsLabel = uilabel(app.SaveTab);
            app.OutputsLabel.Position = [11 516 51 22];
            app.OutputsLabel.Text = 'Outputs:';

            % Create AdditionalLabel
            app.AdditionalLabel = uilabel(app.SaveTab);
            app.AdditionalLabel.Position = [11 371 62 22];
            app.AdditionalLabel.Text = 'Additional:';

            % Create SavealgorithmprogressionplotifselectedCheckBox
            app.SavealgorithmprogressionplotifselectedCheckBox = uicheckbox(app.SaveTab);
            app.SavealgorithmprogressionplotifselectedCheckBox.Enable = 'off';
            app.SavealgorithmprogressionplotifselectedCheckBox.Text = 'Save algorithm progression plot (if selected in the stop conditions tab)';
            app.SavealgorithmprogressionplotifselectedCheckBox.Position = [11 344 399 22];
            app.SavealgorithmprogressionplotifselectedCheckBox.Value = true;

            % Create SavecharacterizationresultsifselectedCheckBox
            app.SavecharacterizationresultsifselectedCheckBox = uicheckbox(app.SaveTab);
            app.SavecharacterizationresultsifselectedCheckBox.Enable = 'off';
            app.SavecharacterizationresultsifselectedCheckBox.Text = 'Save characterization results (if selected in the verification and characterization tab)';
            app.SavecharacterizationresultsifselectedCheckBox.Position = [11 317 475 22];
            app.SavecharacterizationresultsifselectedCheckBox.Value = true;

            % Create SaveparticleadditionalinformationCheckBox
            app.SaveparticleadditionalinformationCheckBox = uicheckbox(app.SaveTab);
            app.SaveparticleadditionalinformationCheckBox.Enable = 'off';
            app.SaveparticleadditionalinformationCheckBox.Text = 'Save particle additional information: particle centroid, diameter, elongation (.mat). Only for volume before upscaling';
            app.SaveparticleadditionalinformationCheckBox.Position = [11 408 646 22];
            app.SaveparticleadditionalinformationCheckBox.Value = true;

            % Create ForphaseandparticlelablelsaveforDropDownLabel
            app.ForphaseandparticlelablelsaveforDropDownLabel = uilabel(app.SaveTab);
            app.ForphaseandparticlelablelsaveforDropDownLabel.HorizontalAlignment = 'right';
            app.ForphaseandparticlelablelsaveforDropDownLabel.Enable = 'off';
            app.ForphaseandparticlelablelsaveforDropDownLabel.Position = [11 435 224 22];
            app.ForphaseandparticlelablelsaveforDropDownLabel.Text = '     For phase and particle lablel, save for';

            % Create ForphaseandparticlelablelsaveforDropDown
            app.ForphaseandparticlelablelsaveforDropDown = uidropdown(app.SaveTab);
            app.ForphaseandparticlelablelsaveforDropDown.Items = {'volume before upscaling', 'volume after upscaling', 'both'};
            app.ForphaseandparticlelablelsaveforDropDown.Enable = 'off';
            app.ForphaseandparticlelablelsaveforDropDown.Position = [250 435 172 22];
            app.ForphaseandparticlelablelsaveforDropDown.Value = 'volume before upscaling';

            % Create Save_folder_button
            app.Save_folder_button = uibutton(app.SaveTab, 'push');
            app.Save_folder_button.ButtonPushedFcn = createCallbackFcn(app, @Save_folder_buttonButtonPushed, true);
            app.Save_folder_button.BackgroundColor = [0.0745 0.6235 1];
            app.Save_folder_button.FontSize = 14;
            app.Save_folder_button.FontWeight = 'bold';
            app.Save_folder_button.FontColor = [1 1 1];
            app.Save_folder_button.Position = [11 671 187 33];
            app.Save_folder_button.Text = 'Click to select save folder';

            % Create NosavefolderselectedLabel
            app.NosavefolderselectedLabel = uilabel(app.SaveTab);
            app.NosavefolderselectedLabel.FontAngle = 'italic';
            app.NosavefolderselectedLabel.Position = [11 644 977 22];
            app.NosavefolderselectedLabel.Text = 'No save folder selected';

            % Create SelectwhereyouwanttosaveinputsandresultsLabel
            app.SelectwhereyouwanttosaveinputsandresultsLabel = uilabel(app.SaveTab);
            app.SelectwhereyouwanttosaveinputsandresultsLabel.FontWeight = 'bold';
            app.SelectwhereyouwanttosaveinputsandresultsLabel.Position = [11 709 304 22];
            app.SelectwhereyouwanttosaveinputsandresultsLabel.Text = '1) Select where you want to save inputs and results';

            % Create ThenselectwhatinformationyouwanttosaveLabel
            app.ThenselectwhatinformationyouwanttosaveLabel = uilabel(app.SaveTab);
            app.ThenselectwhatinformationyouwanttosaveLabel.FontWeight = 'bold';
            app.ThenselectwhatinformationyouwanttosaveLabel.Position = [11 607 289 22];
            app.ThenselectwhatinformationyouwanttosaveLabel.Text = '2) Then select what information you want to save';

            % Create SaveinputsoutpuscomparisonifselectedCheckBox
            app.SaveinputsoutpuscomparisonifselectedCheckBox = uicheckbox(app.SaveTab);
            app.SaveinputsoutpuscomparisonifselectedCheckBox.Enable = 'off';
            app.SaveinputsoutpuscomparisonifselectedCheckBox.Text = 'Save inputs/outpus comparison (if selected in the verification and characterization tab)';
            app.SaveinputsoutpuscomparisonifselectedCheckBox.Position = [11 290 489 22];
            app.SaveinputsoutpuscomparisonifselectedCheckBox.Value = true;

            % Create VolumefractionsTab
            app.VolumefractionsTab = uitab(app.TabGroup);
            app.VolumefractionsTab.Tooltip = {'Phase and volume fractions'};
            app.VolumefractionsTab.Title = 'Volume fractions';
            app.VolumefractionsTab.ForegroundColor = [0.2431 0.4118 0.0275];

            % Create Label
            app.Label = uilabel(app.VolumefractionsTab);
            app.Label.FontColor = [0.851 0.3255 0.098];
            app.Label.Position = [406 578 581 42];
            app.Label.Text = {'If you ask for a high density microstructure, you will have to allow particle overlapping (specified in the'; 'tab of same name) otherwise packing density will block the algorithm progression.'};

            % Create AlgothmisfasterforLabel_5
            app.AlgothmisfasterforLabel_5 = uilabel(app.VolumefractionsTab);
            app.AlgothmisfasterforLabel_5.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_5.Position = [406 625 600 79];
            app.AlgothmisfasterforLabel_5.Text = {'- Microstructure parameter gradations (volume fraction and particle size and morphology) are defined along'; 'the microstructure thickness, which is by convention in this code the direction 3 (see left table).'; '- Properties are defined on some slices, and interpolated in between so that user does not have to specify'; 'values for each slices.'; '- In the table below you can specify volume fractions of the solid phases (porosity is deduced) on some slices.'};

            % Create AlgothmisfasterforLabel_4
            app.AlgothmisfasterforLabel_4 = uilabel(app.VolumefractionsTab);
            app.AlgothmisfasterforLabel_4.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_4.Position = [11 588 355 28];
            app.AlgothmisfasterforLabel_4.Text = {'Algorithm is adimentional. Voxel size is an optional value to'; 'simply convert dimension from voxel length to a physical length.'};

            % Create SetvolumefractionsalongthemicrostructurethicknessLabel
            app.SetvolumefractionsalongthemicrostructurethicknessLabel = uilabel(app.VolumefractionsTab);
            app.SetvolumefractionsalongthemicrostructurethicknessLabel.FontWeight = 'bold';
            app.SetvolumefractionsalongthemicrostructurethicknessLabel.Position = [406 709 433 22];
            app.SetvolumefractionsalongthemicrostructurethicknessLabel.Text = '3) Set volume fractions along the microstructure thickness.';

            % Create AlgothmisfasterforLabel_3
            app.AlgothmisfasterforLabel_3 = uilabel(app.VolumefractionsTab);
            app.AlgothmisfasterforLabel_3.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_3.Position = [11 394 355 22];
            app.AlgothmisfasterforLabel_3.Text = 'You can then choose the phase label and name.';

            % Create Phase_LabelnameUItable
            app.Phase_LabelnameUItable = uitable(app.VolumefractionsTab);
            app.Phase_LabelnameUItable.ColumnName = {'Phase#'; 'Label'; 'Name'};
            app.Phase_LabelnameUItable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Phase_LabelnameUItable.RowName = {};
            app.Phase_LabelnameUItable.ColumnEditable = [false true true];
            app.Phase_LabelnameUItable.CellEditCallback = createCallbackFcn(app, @Phase_LabelnameUItableCellEdit, true);
            app.Phase_LabelnameUItable.Position = [11 19 355 370];

            % Create SetNumberofsolidphasesLabel
            app.SetNumberofsolidphasesLabel = uilabel(app.VolumefractionsTab);
            app.SetNumberofsolidphasesLabel.FontWeight = 'bold';
            app.SetNumberofsolidphasesLabel.Position = [11 448 186 22];
            app.SetNumberofsolidphasesLabel.Text = '2) Set Number of solid phase(s)';

            % Create AlgothmisfasterforLabel_2
            app.AlgothmisfasterforLabel_2 = uilabel(app.VolumefractionsTab);
            app.AlgothmisfasterforLabel_2.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_2.Position = [167 626 25 22];
            app.AlgothmisfasterforLabel_2.Text = '';

            % Create Phase_domainUITable
            app.Phase_domainUITable = uitable(app.VolumefractionsTab);
            app.Phase_domainUITable.ColumnName = {'Direction'; 'Number of voxel'; 'Length (um)'};
            app.Phase_domainUITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Phase_domainUITable.RowName = {};
            app.Phase_domainUITable.ColumnEditable = [false true false];
            app.Phase_domainUITable.CellEditCallback = createCallbackFcn(app, @Phase_domainUITableCellEdit, true);
            app.Phase_domainUITable.Position = [11 480 355 103];

            % Create AlgothmisfasterforLabel
            app.AlgothmisfasterforLabel = uilabel(app.VolumefractionsTab);
            app.AlgothmisfasterforLabel.FontAngle = 'italic';
            app.AlgothmisfasterforLabel.Position = [11 648 355 56];
            app.AlgothmisfasterforLabel.Text = {'Algorithm is faster for larger (particle volume / domain volume)'; 'ratio. It is then recommended to start with a high ratio to setup'; 'properly your parameters, and only then moving to a larger'; 'domain when you are satisfied with your early generation results.'};

            % Create SetdimensionsLabel
            app.SetdimensionsLabel = uilabel(app.VolumefractionsTab);
            app.SetdimensionsLabel.FontWeight = 'bold';
            app.SetdimensionsLabel.Position = [11 709 108 22];
            app.SetdimensionsLabel.Text = '1) Set dimensions';

            % Create Instructions_title_2
            app.Instructions_title_2 = uilabel(app.VolumefractionsTab);
            app.Instructions_title_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_title_2.HorizontalAlignment = 'center';
            app.Instructions_title_2.FontWeight = 'bold';
            app.Instructions_title_2.Position = [11 746 977 22];
            app.Instructions_title_2.Text = 'Specify domain size, number of phases, and volume fractions';

            % Create VoxelsizenmEditField
            app.VoxelsizenmEditField = uieditfield(app.VolumefractionsTab, 'numeric');
            app.VoxelsizenmEditField.Limits = [1e-09 Inf];
            app.VoxelsizenmEditField.ValueChangedFcn = createCallbackFcn(app, @VoxelsizenmEditFieldValueChanged, true);
            app.VoxelsizenmEditField.Position = [117 621 40 22];
            app.VoxelsizenmEditField.Value = 500;

            % Create VoxelsizenmEditFieldLabel
            app.VoxelsizenmEditFieldLabel = uilabel(app.VolumefractionsTab);
            app.VoxelsizenmEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelsizenmEditFieldLabel.Position = [11 621 91 22];
            app.VoxelsizenmEditFieldLabel.Text = 'Voxel size  (nm)';

            % Create NumberofsolidphaseEditField
            app.NumberofsolidphaseEditField = uieditfield(app.VolumefractionsTab, 'numeric');
            app.NumberofsolidphaseEditField.Limits = [1 Inf];
            app.NumberofsolidphaseEditField.RoundFractionalValues = 'on';
            app.NumberofsolidphaseEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofsolidphaseEditFieldValueChanged, true);
            app.NumberofsolidphaseEditField.Position = [152 421 45 22];
            app.NumberofsolidphaseEditField.Value = 1;

            % Create NumberofsolidphaseEditFieldLabel
            app.NumberofsolidphaseEditFieldLabel = uilabel(app.VolumefractionsTab);
            app.NumberofsolidphaseEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofsolidphaseEditFieldLabel.Position = [11 421 126 22];
            app.NumberofsolidphaseEditFieldLabel.Text = 'Number of solid phase';

            % Create NumberofsolidphaseEditField_2Label
            app.NumberofsolidphaseEditField_2Label = uilabel(app.VolumefractionsTab);
            app.NumberofsolidphaseEditField_2Label.HorizontalAlignment = 'right';
            app.NumberofsolidphaseEditField_2Label.Position = [406 551 422 22];
            app.NumberofsolidphaseEditField_2Label.Text = 'First specify the number of slices for which you want to define volume fractions';

            % Create NumberSlice_volumefraction
            app.NumberSlice_volumefraction = uieditfield(app.VolumefractionsTab, 'numeric');
            app.NumberSlice_volumefraction.Limits = [2 Inf];
            app.NumberSlice_volumefraction.RoundFractionalValues = 'on';
            app.NumberSlice_volumefraction.ValueChangedFcn = createCallbackFcn(app, @NumberSlice_volumefractionValueChanged, true);
            app.NumberSlice_volumefraction.HorizontalAlignment = 'center';
            app.NumberSlice_volumefraction.Position = [839 551 45 22];
            app.NumberSlice_volumefraction.Value = 2;

            % Create ThensetthesolidvolumefractionsLabel
            app.ThensetthesolidvolumefractionsLabel = uilabel(app.VolumefractionsTab);
            app.ThensetthesolidvolumefractionsLabel.Position = [407 518 580 28];
            app.ThensetthesolidvolumefractionsLabel.Text = {'Then set the solid volume fractions input and plot them. Notice in the plot the interpolation for the volume'; 'fractions between the slices you have manually specified.'};

            % Create Volumefraction_UITable
            app.Volumefraction_UITable = uitable(app.VolumefractionsTab);
            app.Volumefraction_UITable.ColumnName = {'Normalized position along direction 3'; 'Phase 1'; 'Porosity'};
            app.Volumefraction_UITable.RowName = {};
            app.Volumefraction_UITable.ColumnEditable = [false true false];
            app.Volumefraction_UITable.CellEditCallback = createCallbackFcn(app, @Volumefraction_UITableCellEdit, true);
            app.Volumefraction_UITable.Position = [407 71 581 442];

            % Create PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel = uilabel(app.VolumefractionsTab);
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel.FontWeight = 'bold';
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel.Position = [407 19 241 28];
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel.Text = {'4) Plot (optional) to check your input and'; '(mandatory) save them to move on.'};

            % Create Phase_vf_plot
            app.Phase_vf_plot = uibutton(app.VolumefractionsTab, 'push');
            app.Phase_vf_plot.ButtonPushedFcn = createCallbackFcn(app, @Phase_vf_plotButtonPushed, true);
            app.Phase_vf_plot.BackgroundColor = [0.0745 0.6235 1];
            app.Phase_vf_plot.FontSize = 14;
            app.Phase_vf_plot.FontWeight = 'bold';
            app.Phase_vf_plot.FontColor = [1 1 1];
            app.Phase_vf_plot.Position = [715 19 158 33];
            app.Phase_vf_plot.Text = 'Plot volume fractions';

            % Create Phase_vf_save
            app.Phase_vf_save = uibutton(app.VolumefractionsTab, 'push');
            app.Phase_vf_save.ButtonPushedFcn = createCallbackFcn(app, @Phase_vf_saveButtonPushed, true);
            app.Phase_vf_save.BackgroundColor = [0.0745 0.6235 1];
            app.Phase_vf_save.FontSize = 14;
            app.Phase_vf_save.FontWeight = 'bold';
            app.Phase_vf_save.FontColor = [1 1 1];
            app.Phase_vf_save.Position = [891 19 97 33];
            app.Phase_vf_save.Text = 'Save';

            % Create DiametersTab
            app.DiametersTab = uitab(app.TabGroup);
            app.DiametersTab.Tooltip = {'Particle diameter'};
            app.DiametersTab.Title = 'Diameters';
            app.DiametersTab.BackgroundColor = 'none';
            app.DiametersTab.ForegroundColor = [0 0 1];

            % Create Instructions_diameter
            app.Instructions_diameter = uilabel(app.DiametersTab);
            app.Instructions_diameter.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_diameter.HorizontalAlignment = 'center';
            app.Instructions_diameter.FontWeight = 'bold';
            app.Instructions_diameter.Position = [11 746 977 22];
            app.Instructions_diameter.Text = 'For each solid phase, set the particle diameter distributions along the microstructure thickness';

            % Create SelectphaseDropDownLabel
            app.SelectphaseDropDownLabel = uilabel(app.DiametersTab);
            app.SelectphaseDropDownLabel.HorizontalAlignment = 'right';
            app.SelectphaseDropDownLabel.Enable = 'off';
            app.SelectphaseDropDownLabel.Position = [11 659 75 22];
            app.SelectphaseDropDownLabel.Text = 'Select phase';

            % Create Selectphase_diameters_DropDown
            app.Selectphase_diameters_DropDown = uidropdown(app.DiametersTab);
            app.Selectphase_diameters_DropDown.Items = {};
            app.Selectphase_diameters_DropDown.ValueChangedFcn = createCallbackFcn(app, @Selectphase_diameters_DropDownValueChanged, true);
            app.Selectphase_diameters_DropDown.Enable = 'off';
            app.Selectphase_diameters_DropDown.Position = [101 659 188 22];
            app.Selectphase_diameters_DropDown.Value = {};

            % Create SelectasolidphaseandsetthenumberofLabel
            app.SelectasolidphaseandsetthenumberofLabel = uilabel(app.DiametersTab);
            app.SelectasolidphaseandsetthenumberofLabel.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel.Position = [11 709 974 22];
            app.SelectasolidphaseandsetthenumberofLabel.Text = '1) Select a solid phase and set the number of slices for which you want to define particle diameter/elongation, and also how many values for each of these properties.';

            % Create AlgothmisfasterforLabel_6
            app.AlgothmisfasterforLabel_6 = uilabel(app.DiametersTab);
            app.AlgothmisfasterforLabel_6.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_6.Position = [11 686 417 18];
            app.AlgothmisfasterforLabel_6.Text = 'Your choice will format the table in the section 2).';

            % Create Diameters_main_UITable
            app.Diameters_main_UITable = uitable(app.DiametersTab);
            app.Diameters_main_UITable.ColumnName = {'Parameters'; 'Number of slices'; 'Number of values'};
            app.Diameters_main_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Diameters_main_UITable.RowName = {};
            app.Diameters_main_UITable.ColumnEditable = [false true true];
            app.Diameters_main_UITable.CellEditCallback = createCallbackFcn(app, @Diameters_main_UITableCellEdit, true);
            app.Diameters_main_UITable.Enable = 'off';
            app.Diameters_main_UITable.Position = [11 550 417 104];

            % Create InstructionsTextAreaLabel
            app.InstructionsTextAreaLabel = uilabel(app.DiametersTab);
            app.InstructionsTextAreaLabel.HorizontalAlignment = 'right';
            app.InstructionsTextAreaLabel.Position = [462 680 67 22];
            app.InstructionsTextAreaLabel.Text = 'Instructions';

            % Create InstructionsTextArea
            app.InstructionsTextArea = uitextarea(app.DiametersTab);
            app.InstructionsTextArea.Position = [544 550 442 154];
            app.InstructionsTextArea.Value = {'- Each phase is setup by default with a unique particle diameter and with no particle elongation (i.e., spheres).'; '- You can set up your custom diameter and elongation distributions for each phase. You will have to repeat all steps for each phase.'; '- Before moving to a next phase you have to click on the save button to keep your choice.'; '- Diameters and elongations are defined for each normalized position along the thickness not by a single scalar (as done for volume fractions in the previous tab) but with an histogram.'; '- You can visualize particle morphology in the particle viusalization tab for testing.'};

            % Create SelectasolidphaseandsetthenumberofLabel_2
            app.SelectasolidphaseandsetthenumberofLabel_2 = uilabel(app.DiametersTab);
            app.SelectasolidphaseandsetthenumberofLabel_2.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_2.Position = [11 518 1005 22];
            app.SelectasolidphaseandsetthenumberofLabel_2.Text = '2) Set diameter and elongation. If you set n values in the table above, you will have  *1, *2,... *n. Set their value in the first row. Set their distribution in the subsequent rows.';

            % Create Diameters_UITable
            app.Diameters_UITable = uitable(app.DiametersTab);
            app.Diameters_UITable.ColumnName = {'Normalized position along direction 3'; 'Dx 1'; 'Total'};
            app.Diameters_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Diameters_UITable.RowName = {};
            app.Diameters_UITable.ColumnEditable = [false true false];
            app.Diameters_UITable.CellEditCallback = createCallbackFcn(app, @Diameters_UITableCellEdit, true);
            app.Diameters_UITable.Enable = 'off';
            app.Diameters_UITable.Position = [11 334 817 152];

            % Create Rowii1isthediameterhistogramLabel
            app.Rowii1isthediameterhistogramLabel = uilabel(app.DiametersTab);
            app.Rowii1isthediameterhistogramLabel.FontAngle = 'italic';
            app.Rowii1isthediameterhistogramLabel.Position = [835 398 150 84];
            app.Rowii1isthediameterhistogramLabel.Text = {'Row i (i>1) is the diameter'; 'histogram expressed in'; 'percentages. Total must be'; '100%.'; ''; 'Ignore NaN values'};

            % Create AlgothmisfasterforLabel_7
            app.AlgothmisfasterforLabel_7 = uilabel(app.DiametersTab);
            app.AlgothmisfasterforLabel_7.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_7.Position = [11 491 977 22];
            app.AlgothmisfasterforLabel_7.Text = 'First table is for diameter Dx (set values in voxel length). Second and third tables are diameter ratios (i.e., equal to 1 for a sphere).';

            % Create DxDy_UITable
            app.DxDy_UITable = uitable(app.DiametersTab);
            app.DxDy_UITable.ColumnName = {'Normalized position along direction 3'; 'Dx/Dy 1'; 'Total'};
            app.DxDy_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.DxDy_UITable.RowName = {};
            app.DxDy_UITable.ColumnEditable = [false true false];
            app.DxDy_UITable.CellEditCallback = createCallbackFcn(app, @DxDy_UITableCellEdit, true);
            app.DxDy_UITable.Enable = 'off';
            app.DxDy_UITable.Position = [11 177 817 152];

            % Create DxDz_UITable
            app.DxDz_UITable = uitable(app.DiametersTab);
            app.DxDz_UITable.ColumnName = {'Normalized position along direction 3'; 'Dx/Dz 1'; 'Total'};
            app.DxDz_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.DxDz_UITable.RowName = {};
            app.DxDz_UITable.ColumnEditable = [false true false];
            app.DxDz_UITable.CellEditCallback = createCallbackFcn(app, @DxDz_UITableCellEdit, true);
            app.DxDz_UITable.Enable = 'off';
            app.DxDz_UITable.Position = [11 20 817 152];

            % Create Diameter_plot
            app.Diameter_plot = uibutton(app.DiametersTab, 'push');
            app.Diameter_plot.ButtonPushedFcn = createCallbackFcn(app, @Diameter_plotButtonPushed, true);
            app.Diameter_plot.BackgroundColor = [0.0745 0.6235 1];
            app.Diameter_plot.FontSize = 14;
            app.Diameter_plot.FontWeight = 'bold';
            app.Diameter_plot.FontColor = [1 1 1];
            app.Diameter_plot.Enable = 'off';
            app.Diameter_plot.Position = [851 76 116 40];
            app.Diameter_plot.Text = {'Plot'; 'distributions'};

            % Create Diameter_save
            app.Diameter_save = uibutton(app.DiametersTab, 'push');
            app.Diameter_save.ButtonPushedFcn = createCallbackFcn(app, @Diameter_saveButtonPushed, true);
            app.Diameter_save.BackgroundColor = [0.0745 0.6235 1];
            app.Diameter_save.FontSize = 14;
            app.Diameter_save.FontWeight = 'bold';
            app.Diameter_save.FontColor = [1 1 1];
            app.Diameter_save.Enable = 'off';
            app.Diameter_save.Position = [851 33 116 33];
            app.Diameter_save.Text = 'Save';

            % Create PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_2
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_2 = uilabel(app.DiametersTab);
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_2.FontWeight = 'bold';
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_2.Position = [835 129 168 98];
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_2.Text = {'3) Plot (optional) to check'; 'your inputs and (mandatory)'; 'save them before moving'; 'to the next phase.'; '4) Repeat steps 1,2,3 for'; 'each phase.'; '5) Move to next tab.'};

            % Create RotationsTab
            app.RotationsTab = uitab(app.TabGroup);
            app.RotationsTab.Tooltip = {'Particle rotation (orientation)'};
            app.RotationsTab.Title = 'Rotations';
            app.RotationsTab.BackgroundColor = 'none';
            app.RotationsTab.ForegroundColor = [0 0 1];

            % Create Instructions_orientation
            app.Instructions_orientation = uilabel(app.RotationsTab);
            app.Instructions_orientation.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_orientation.HorizontalAlignment = 'center';
            app.Instructions_orientation.FontWeight = 'bold';
            app.Instructions_orientation.Position = [11 746 977 22];
            app.Instructions_orientation.Text = 'For each solid phase, set the particle rotation distributions along the microstructure thickness';

            % Create SelectphaseDropDownLabel_2
            app.SelectphaseDropDownLabel_2 = uilabel(app.RotationsTab);
            app.SelectphaseDropDownLabel_2.HorizontalAlignment = 'right';
            app.SelectphaseDropDownLabel_2.Enable = 'off';
            app.SelectphaseDropDownLabel_2.Position = [11 659 75 22];
            app.SelectphaseDropDownLabel_2.Text = 'Select phase';

            % Create Selectphase_orientation_DropDown
            app.Selectphase_orientation_DropDown = uidropdown(app.RotationsTab);
            app.Selectphase_orientation_DropDown.Items = {};
            app.Selectphase_orientation_DropDown.ValueChangedFcn = createCallbackFcn(app, @Selectphase_orientation_DropDownValueChanged, true);
            app.Selectphase_orientation_DropDown.Enable = 'off';
            app.Selectphase_orientation_DropDown.Position = [101 659 100 22];
            app.Selectphase_orientation_DropDown.Value = {};

            % Create SelectasolidphaseandsetthenumberofLabel_3
            app.SelectasolidphaseandsetthenumberofLabel_3 = uilabel(app.RotationsTab);
            app.SelectasolidphaseandsetthenumberofLabel_3.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_3.Position = [11 709 909 22];
            app.SelectasolidphaseandsetthenumberofLabel_3.Text = '1) Select a solid phase and set the number of slices for which you want to define particle rotations, and also how many values for each of these properties.';

            % Create AlgothmisfasterforLabel_8
            app.AlgothmisfasterforLabel_8 = uilabel(app.RotationsTab);
            app.AlgothmisfasterforLabel_8.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_8.Position = [11 686 417 18];
            app.AlgothmisfasterforLabel_8.Text = 'Your choice will format the table in the section 2).';

            % Create Orientation_main_UITable
            app.Orientation_main_UITable = uitable(app.RotationsTab);
            app.Orientation_main_UITable.ColumnName = {'Parameters'; 'Number of slices'; 'Number of values'};
            app.Orientation_main_UITable.RowName = {};
            app.Orientation_main_UITable.ColumnEditable = [false true true];
            app.Orientation_main_UITable.CellEditCallback = createCallbackFcn(app, @Orientation_main_UITableCellEdit, true);
            app.Orientation_main_UITable.Enable = 'off';
            app.Orientation_main_UITable.Position = [11 550 453 104];

            % Create InstructionsTextArea_2Label
            app.InstructionsTextArea_2Label = uilabel(app.RotationsTab);
            app.InstructionsTextArea_2Label.HorizontalAlignment = 'right';
            app.InstructionsTextArea_2Label.Position = [462 680 67 22];
            app.InstructionsTextArea_2Label.Text = 'Instructions';

            % Create InstructionsTextArea_2
            app.InstructionsTextArea_2 = uitextarea(app.RotationsTab);
            app.InstructionsTextArea_2.Position = [544 550 442 154];
            app.InstructionsTextArea_2.Value = {'- Each phase is setup by default with a unique particle rotation (or orientation)'; 'for each cartesian axis x,y,z.'; '- You can set up your custom rotation distributions for each phase. You will have to repeat all steps for each phase.'; '- Before moving to a next phase you have to click on the save button to keep your choice.'; '- Rotations are defined for each normalized position along the thickness not by a single scalar but with an histogram (as done for particle diameters in the previous tab).'; '- You can visualize particle morphology in the particle viusalization tab for testing.'};

            % Create SelectasolidphaseandsetthenumberofLabel_4
            app.SelectasolidphaseandsetthenumberofLabel_4 = uilabel(app.RotationsTab);
            app.SelectasolidphaseandsetthenumberofLabel_4.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_4.Position = [11 518 935 22];
            app.SelectasolidphaseandsetthenumberofLabel_4.Text = '2) Set Orientations. If you set n values in the table above, you will have  *1, *2,... *n. Set their value in the first row. Set their distribution in the subsequent rows.';

            % Create Rotation_x_UITable
            app.Rotation_x_UITable = uitable(app.RotationsTab);
            app.Rotation_x_UITable.ColumnName = {'Normalized position along direction 3'; 'Dx 1'; 'Total'};
            app.Rotation_x_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Rotation_x_UITable.RowName = {};
            app.Rotation_x_UITable.ColumnEditable = [false true false];
            app.Rotation_x_UITable.CellEditCallback = createCallbackFcn(app, @Rotation_x_UITableCellEdit, true);
            app.Rotation_x_UITable.Enable = 'off';
            app.Rotation_x_UITable.Position = [11 334 817 152];

            % Create Rowii1isthediameterhistogramLabel_2
            app.Rowii1isthediameterhistogramLabel_2 = uilabel(app.RotationsTab);
            app.Rowii1isthediameterhistogramLabel_2.FontAngle = 'italic';
            app.Rowii1isthediameterhistogramLabel_2.Position = [835 398 150 84];
            app.Rowii1isthediameterhistogramLabel_2.Text = {'Row i (i>1) is the rotation'; 'histogram expressed in'; 'percentages. Total must be'; '100%.'; ''; 'Ignore NaN values'};

            % Create AlgothmisfasterforLabel_9
            app.AlgothmisfasterforLabel_9 = uilabel(app.RotationsTab);
            app.AlgothmisfasterforLabel_9.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_9.Position = [11 491 977 22];
            app.AlgothmisfasterforLabel_9.Text = 'First, second, and third table are for rotation normal to x-axis,y-axis, and z-axis, respectively. Values are in degrees. Do not use 0, as 180 degree is used instead.';

            % Create Rotation_y_UITable
            app.Rotation_y_UITable = uitable(app.RotationsTab);
            app.Rotation_y_UITable.ColumnName = {'Normalized position along direction 3'; 'Dx/Dy 1'; 'Total'};
            app.Rotation_y_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Rotation_y_UITable.RowName = {};
            app.Rotation_y_UITable.ColumnEditable = [false true false];
            app.Rotation_y_UITable.CellEditCallback = createCallbackFcn(app, @Rotation_y_UITableCellEdit, true);
            app.Rotation_y_UITable.Enable = 'off';
            app.Rotation_y_UITable.Position = [11 177 817 152];

            % Create Rotation_z_UITable
            app.Rotation_z_UITable = uitable(app.RotationsTab);
            app.Rotation_z_UITable.ColumnName = {'Normalized position along direction 3'; 'Dx/Dz 1'; 'Total'};
            app.Rotation_z_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Rotation_z_UITable.RowName = {};
            app.Rotation_z_UITable.ColumnEditable = [false true false];
            app.Rotation_z_UITable.CellEditCallback = createCallbackFcn(app, @Rotation_z_UITableCellEdit, true);
            app.Rotation_z_UITable.Enable = 'off';
            app.Rotation_z_UITable.Position = [11 20 817 152];

            % Create Rotation_plot
            app.Rotation_plot = uibutton(app.RotationsTab, 'push');
            app.Rotation_plot.ButtonPushedFcn = createCallbackFcn(app, @Rotation_plotButtonPushed, true);
            app.Rotation_plot.BackgroundColor = [0.0745 0.6235 1];
            app.Rotation_plot.FontSize = 14;
            app.Rotation_plot.FontWeight = 'bold';
            app.Rotation_plot.FontColor = [1 1 1];
            app.Rotation_plot.Enable = 'off';
            app.Rotation_plot.Position = [851 76 116 40];
            app.Rotation_plot.Text = {'Plot'; 'distributions'};

            % Create Rotation_save
            app.Rotation_save = uibutton(app.RotationsTab, 'push');
            app.Rotation_save.ButtonPushedFcn = createCallbackFcn(app, @Rotation_saveButtonPushed, true);
            app.Rotation_save.BackgroundColor = [0.0745 0.6235 1];
            app.Rotation_save.FontSize = 14;
            app.Rotation_save.FontWeight = 'bold';
            app.Rotation_save.FontColor = [1 1 1];
            app.Rotation_save.Enable = 'off';
            app.Rotation_save.Position = [851 33 116 33];
            app.Rotation_save.Text = 'Save';

            % Create PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_3
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_3 = uilabel(app.RotationsTab);
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_3.FontWeight = 'bold';
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_3.Position = [835 129 168 98];
            app.PlotoptionaltocheckyourinputandmandatorysavethemtomoveonLabel_3.Text = {'3) Plot (optional) to check'; 'your inputs and (mandatory)'; 'save them before moving'; 'to the next phase.'; '4) Repeat steps 1,2,3 for'; 'each phase.'; '5) Move to next tab.'};

            % Create ParticlevisualizationTab
            app.ParticlevisualizationTab = uitab(app.TabGroup);
            app.ParticlevisualizationTab.Tooltip = {'Particle visualization'};
            app.ParticlevisualizationTab.Title = 'Particle visualization';
            app.ParticlevisualizationTab.ForegroundColor = [0 0 1];

            % Create Instructions_Visualization
            app.Instructions_Visualization = uilabel(app.ParticlevisualizationTab);
            app.Instructions_Visualization.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_Visualization.HorizontalAlignment = 'center';
            app.Instructions_Visualization.FontWeight = 'bold';
            app.Instructions_Visualization.Position = [11 746 977 22];
            app.Instructions_Visualization.Text = 'Here you can test the diameter and rotation parameters for a single particle';

            % Create SelectasolidphaseandsetthenumberofLabel_5
            app.SelectasolidphaseandsetthenumberofLabel_5 = uilabel(app.ParticlevisualizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_5.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_5.Position = [11 709 348 22];
            app.SelectasolidphaseandsetthenumberofLabel_5.Text = '1) Enter parameters in the table below and plot the particle.';

            % Create AlgothmisfasterforLabel_10
            app.AlgothmisfasterforLabel_10 = uilabel(app.ParticlevisualizationTab);
            app.AlgothmisfasterforLabel_10.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_10.Position = [11 682 919 22];
            app.AlgothmisfasterforLabel_10.Text = 'This is only for visualization. Parameters entered here are not used elsewhere. The aim of this tab is to familiarize the user with the parameters of the two previous tabs.';

            % Create Parameter_example_UITable
            app.Parameter_example_UITable = uitable(app.ParticlevisualizationTab);
            app.Parameter_example_UITable.ColumnName = {'Parameter'; 'Value'; 'Unit'; 'Range'};
            app.Parameter_example_UITable.RowName = {};
            app.Parameter_example_UITable.ColumnEditable = [false true false false];
            app.Parameter_example_UITable.CellEditCallback = createCallbackFcn(app, @Parameter_example_UITableCellEdit, true);
            app.Parameter_example_UITable.Position = [11 512 608 165];

            % Create Plot_example_particle
            app.Plot_example_particle = uibutton(app.ParticlevisualizationTab, 'push');
            app.Plot_example_particle.ButtonPushedFcn = createCallbackFcn(app, @Plot_example_particleButtonPushed, true);
            app.Plot_example_particle.BackgroundColor = [0.0745 0.6235 1];
            app.Plot_example_particle.FontSize = 14;
            app.Plot_example_particle.FontWeight = 'bold';
            app.Plot_example_particle.FontColor = [1 1 1];
            app.Plot_example_particle.Position = [11 462 317 40];
            app.Plot_example_particle.Text = {'Plot example particle'; 'with patch function (slow but more accurate)'};

            % Create Plot_example_particle_2
            app.Plot_example_particle_2 = uibutton(app.ParticlevisualizationTab, 'push');
            app.Plot_example_particle_2.ButtonPushedFcn = createCallbackFcn(app, @Plot_example_particle_2ButtonPushed, true);
            app.Plot_example_particle_2.BackgroundColor = [0.0745 0.6235 1];
            app.Plot_example_particle_2.FontSize = 14;
            app.Plot_example_particle_2.FontWeight = 'bold';
            app.Plot_example_particle_2.FontColor = [1 1 1];
            app.Plot_example_particle_2.Position = [11 412 317 40];
            app.Plot_example_particle_2.Text = {'Plot example particle'; 'with volview (fast but less accurate)'};

            % Create OverlappingTab
            app.OverlappingTab = uitab(app.TabGroup);
            app.OverlappingTab.Tooltip = {'Particle overlapping'};
            app.OverlappingTab.Title = 'Overlapping';
            app.OverlappingTab.ForegroundColor = [1 0.302 0];

            % Create Instructions_overlapping
            app.Instructions_overlapping = uilabel(app.OverlappingTab);
            app.Instructions_overlapping.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_overlapping.HorizontalAlignment = 'center';
            app.Instructions_overlapping.FontWeight = 'bold';
            app.Instructions_overlapping.Position = [11 746 977 22];
            app.Instructions_overlapping.Text = 'Set particle overlapping paramters';

            % Create Label_2
            app.Label_2 = uilabel(app.OverlappingTab);
            app.Label_2.Position = [11 675 977 56];
            app.Label_2.Text = {'Theoretical maximum density for unisize sphere packing is 74% [1], and for random spatial distribution of unisize spheres 63.4% [2].'; 'Overlapping enables to overcome this packing density limit.'; '[1] A. Bezdek and W. Kuperberg, Arxiv (2010).'; '[2] C. Song, P. Wang, and H. A. Makse, Nature, 453, 629–632 (2008).'};

            % Create Label_3
            app.Label_3 = uilabel(app.OverlappingTab);
            app.Label_3.FontWeight = 'bold';
            app.Label_3.Position = [11 647 977 22];
            app.Label_3.Text = 'You can control particle overlapping between particles of phase i and particle of phase j through the overlapping symmetric matrix.';

            % Create Overlapping_Image
            app.Overlapping_Image = uiimage(app.OverlappingTab);
            app.Overlapping_Image.Position = [11 258 550 374];
            app.Overlapping_Image.ImageSource = 'Particle_overlapping.png';

            % Create overlapping_UITable
            app.overlapping_UITable = uitable(app.OverlappingTab);
            app.overlapping_UITable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.overlapping_UITable.RowName = {};
            app.overlapping_UITable.CellEditCallback = createCallbackFcn(app, @overlapping_UITableCellEdit, true);
            app.overlapping_UITable.Enable = 'off';
            app.overlapping_UITable.Visible = 'off';
            app.overlapping_UITable.Position = [573 258 415 374];

            % Create Label_4
            app.Label_4 = uilabel(app.OverlappingTab);
            app.Label_4.FontWeight = 'bold';
            app.Label_4.Position = [8 187 995 56];
            app.Label_4.Text = {'The above condition is checked for each new particle. However, cumulative overlapping may degrade too much particle morphology, even if each indidivual overlapping'; 'passes the above condition. You can specify below the (normalized) minimum volume each particle must preserve from its initial volume when generated for each phase.'; 'If an overlapping would remove particle volume so that its remaining volume would go below this value, the particle is not generated and the algorithm goes to the next'; 'iteration.'; ''};

            % Create Minimumvolume_overlapping_UITable
            app.Minimumvolume_overlapping_UITable = uitable(app.OverlappingTab);
            app.Minimumvolume_overlapping_UITable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.Minimumvolume_overlapping_UITable.RowName = {};
            app.Minimumvolume_overlapping_UITable.CellEditCallback = createCallbackFcn(app, @Minimumvolume_overlapping_UITableCellEdit, true);
            app.Minimumvolume_overlapping_UITable.Enable = 'off';
            app.Minimumvolume_overlapping_UITable.Visible = 'off';
            app.Minimumvolume_overlapping_UITable.Position = [8 115 977 62];

            % Create DisableparticlecontiguitycheckCheckBox
            app.DisableparticlecontiguitycheckCheckBox = uicheckbox(app.OverlappingTab);
            app.DisableparticlecontiguitycheckCheckBox.Text = 'Disable particle contiguity check';
            app.DisableparticlecontiguitycheckCheckBox.Position = [7 45 194 22];

            % Create Label_5
            app.Label_5 = uilabel(app.OverlappingTab);
            app.Label_5.FontWeight = 'bold';
            app.Label_5.Position = [8 72 986 28];
            app.Label_5.Text = {'For very elongated particles, overlapping may cut particles. By default contiguity is verified and if a particle is cut into several parts due to the newly generated particle,'; 'the latter is canceled and the algorithm goes to the next iteration. You can disable it (not recommended) if you have difficutly generating elongated particles.'};

            % Create Overlapping_save
            app.Overlapping_save = uibutton(app.OverlappingTab, 'push');
            app.Overlapping_save.ButtonPushedFcn = createCallbackFcn(app, @Overlapping_saveButtonPushed, true);
            app.Overlapping_save.BackgroundColor = [0.0745 0.6235 1];
            app.Overlapping_save.FontSize = 14;
            app.Overlapping_save.FontWeight = 'bold';
            app.Overlapping_save.FontColor = [1 1 1];
            app.Overlapping_save.Enable = 'off';
            app.Overlapping_save.Position = [869 20 116 33];
            app.Overlapping_save.Text = 'Save';

            % Create OrderTab
            app.OrderTab = uitab(app.TabGroup);
            app.OrderTab.Tooltip = {'Particle generation order'};
            app.OrderTab.Title = 'Order';
            app.OrderTab.ForegroundColor = [1 0.302 0];

            % Create Instructions_order
            app.Instructions_order = uilabel(app.OrderTab);
            app.Instructions_order.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_order.HorizontalAlignment = 'center';
            app.Instructions_order.FontWeight = 'bold';
            app.Instructions_order.Position = [11 746 977 22];
            app.Instructions_order.Text = 'Set particle generation order';

            % Create SelectasolidphaseandsetthenumberofLabel_6
            app.SelectasolidphaseandsetthenumberofLabel_6 = uilabel(app.OrderTab);
            app.SelectasolidphaseandsetthenumberofLabel_6.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_6.Position = [11 689 1007 42];
            app.SelectasolidphaseandsetthenumberofLabel_6.Text = {'Algorithm generation may stall if inputs are ill-defined. One example would be try generating a phase with two very different particle diameters: large and tiny particles.'; 'While finding place for the tiny particles is easy, it may become impossible to find place for the large ones, especially as packing density is increasing.'; 'A workaround would be to first generate all the large particles, and only after the small ones, i.e., controlling the particle generation order, as illustrated in the figure below. '};

            % Create Image
            app.Image = uiimage(app.OrderTab);
            app.Image.Position = [165 591 700 88];
            app.Image.ImageSource = 'Particle_generation_order.png';

            % Create SelectasolidphaseandsetthenumberofLabel_7
            app.SelectasolidphaseandsetthenumberofLabel_7 = uilabel(app.OrderTab);
            app.SelectasolidphaseandsetthenumberofLabel_7.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_7.Position = [11 543 975 28];
            app.SelectasolidphaseandsetthenumberofLabel_7.Text = {'You can specify below the number of generation pass and the target ratio matrix r for each phase and pass iteration. It reads as follow, for the value in row i, column j:'; 'For the pass iteration i, generate for the phase j, r(i,j) times the target volume fraction of the phase j.'};

            % Create PseudocodeTextAreaLabel
            app.PseudocodeTextAreaLabel = uilabel(app.OrderTab);
            app.PseudocodeTextAreaLabel.HorizontalAlignment = 'right';
            app.PseudocodeTextAreaLabel.Position = [11 514 76 22];
            app.PseudocodeTextAreaLabel.Text = 'Pseudo-code';

            % Create PseudocodeTextArea
            app.PseudocodeTextArea = uitextarea(app.OrderTab);
            app.PseudocodeTextArea.Position = [102 469 522 69];
            app.PseudocodeTextArea.Value = {'Loop over pass'; '   Loop over phase'; '      While volume fraction(phase) ≤ target ratio (pass) x target volume fraction(phase) '; '         Generate new particle'};

            % Create NumberofgenerationpassEditFieldLabel
            app.NumberofgenerationpassEditFieldLabel = uilabel(app.OrderTab);
            app.NumberofgenerationpassEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofgenerationpassEditFieldLabel.Position = [11 421 150 22];
            app.NumberofgenerationpassEditFieldLabel.Text = 'Number of generation pass';

            % Create NumberofgenerationpassEditField
            app.NumberofgenerationpassEditField = uieditfield(app.OrderTab, 'numeric');
            app.NumberofgenerationpassEditField.Limits = [1 Inf];
            app.NumberofgenerationpassEditField.RoundFractionalValues = 'on';
            app.NumberofgenerationpassEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofgenerationpassEditFieldValueChanged, true);
            app.NumberofgenerationpassEditField.Position = [176 421 42 22];
            app.NumberofgenerationpassEditField.Value = 1;

            % Create Label_6
            app.Label_6 = uilabel(app.OrderTab);
            app.Label_6.FontAngle = 'italic';
            app.Label_6.Position = [11 20 817 140];
            app.Label_6.Text = {'For the example mentioned above (the tiny and large particles), you would need to assign the tiny and large particles to different phases.'; 'The parameter would be: 2 phases, 1 or 2 generation pass, with the following fill ratio matrix:'; ''; '              Large particles  | Small particles'; 'Pass 1 |             1                       1         '; ''; '              Small particles  | Large particles'; 'Pass 1 |             0                       1         '; 'Pass 2 |             1                       1'; 'Note that for the second pass, the generation for the large particle is skipped as the target has been already reached at the end of the first pass.'};

            % Create Order_UITable
            app.Order_UITable = uitable(app.OrderTab);
            app.Order_UITable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.Order_UITable.RowName = {};
            app.Order_UITable.CellEditCallback = createCallbackFcn(app, @Order_UITableCellEdit, true);
            app.Order_UITable.Enable = 'off';
            app.Order_UITable.Visible = 'off';
            app.Order_UITable.Position = [11 176 613 235];

            % Create Order_save
            app.Order_save = uibutton(app.OrderTab, 'push');
            app.Order_save.ButtonPushedFcn = createCallbackFcn(app, @Order_saveButtonPushed, true);
            app.Order_save.BackgroundColor = [0.0745 0.6235 1];
            app.Order_save.FontSize = 14;
            app.Order_save.FontWeight = 'bold';
            app.Order_save.FontColor = [1 1 1];
            app.Order_save.Enable = 'off';
            app.Order_save.Position = [869 20 116 33];
            app.Order_save.Text = 'Save';

            % Create StopconditionsTab
            app.StopconditionsTab = uitab(app.TabGroup);
            app.StopconditionsTab.Tooltip = {'Algorithm stop conditions (stalling)'};
            app.StopconditionsTab.Title = 'Stop conditions';
            app.StopconditionsTab.ForegroundColor = [1 0.302 0];

            % Create Instructions_order_2
            app.Instructions_order_2 = uilabel(app.StopconditionsTab);
            app.Instructions_order_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_order_2.HorizontalAlignment = 'center';
            app.Instructions_order_2.FontWeight = 'bold';
            app.Instructions_order_2.Position = [11 746 977 22];
            app.Instructions_order_2.Text = 'Algorithm progression and stop conditions';

            % Create PlotalgorithmprogressionCheckBox
            app.PlotalgorithmprogressionCheckBox = uicheckbox(app.StopconditionsTab);
            app.PlotalgorithmprogressionCheckBox.Text = 'Plot algorithm progression';
            app.PlotalgorithmprogressionCheckBox.Position = [11 662 161 22];
            app.PlotalgorithmprogressionCheckBox.Value = true;

            % Create SelectasolidphaseandsetthenumberofLabel_8
            app.SelectasolidphaseandsetthenumberofLabel_8 = uilabel(app.StopconditionsTab);
            app.SelectasolidphaseandsetthenumberofLabel_8.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_8.Position = [11 689 976 42];
            app.SelectasolidphaseandsetthenumberofLabel_8.Text = {'Algorithm may stall if inputs are ill-defined or unrealistic (e.g., asking for a density higher than the packing density limit wihout particle overlapping).'; 'You can plot the volume fraction, volume fraction rate, and particle generation rate as function of wallclock time to evaluate if a long generation time is due to a too small'; '(particle volume / domain volume) ratio or due to unrealistic inputs.'};

            % Create SelectasolidphaseandsetthenumberofLabel_9
            app.SelectasolidphaseandsetthenumberofLabel_9 = uilabel(app.StopconditionsTab);
            app.SelectasolidphaseandsetthenumberofLabel_9.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_9.Position = [11 619 976 23];
            app.SelectasolidphaseandsetthenumberofLabel_9.Text = 'You can specify stop conditions in case the algorithm is stalling, to avoid losing time.';

            % Create Ifvolumefractions1goesbelowEditFieldLabel
            app.Ifvolumefractions1goesbelowEditFieldLabel = uilabel(app.StopconditionsTab);
            app.Ifvolumefractions1goesbelowEditFieldLabel.HorizontalAlignment = 'right';
            app.Ifvolumefractions1goesbelowEditFieldLabel.Position = [11 592 181 22];
            app.Ifvolumefractions1goesbelowEditFieldLabel.Text = 'If volume fraction.s-1 goes below';

            % Create Ifvolumefractions1goesbelowEditField
            app.Ifvolumefractions1goesbelowEditField = uieditfield(app.StopconditionsTab, 'numeric');
            app.Ifvolumefractions1goesbelowEditField.Limits = [0 Inf];
            app.Ifvolumefractions1goesbelowEditField.Position = [207 592 105 22];
            app.Ifvolumefractions1goesbelowEditField.Value = 0.0001;

            % Create particles1goesbelowEditFieldLabel
            app.particles1goesbelowEditFieldLabel = uilabel(app.StopconditionsTab);
            app.particles1goesbelowEditFieldLabel.HorizontalAlignment = 'right';
            app.particles1goesbelowEditFieldLabel.Position = [392 592 128 22];
            app.particles1goesbelowEditFieldLabel.Text = 'particle.s-1 goes below';

            % Create particles1goesbelowEditField
            app.particles1goesbelowEditField = uieditfield(app.StopconditionsTab, 'numeric');
            app.particles1goesbelowEditField.Limits = [0 Inf];
            app.particles1goesbelowEditField.Position = [535 592 105 22];
            app.particles1goesbelowEditField.Value = 0.1;

            % Create stopcondition_andor_DropDown
            app.stopcondition_andor_DropDown = uidropdown(app.StopconditionsTab);
            app.stopcondition_andor_DropDown.Items = {'and', 'or'};
            app.stopcondition_andor_DropDown.Position = [322 592 60 22];
            app.stopcondition_andor_DropDown.Value = 'or';

            % Create whenaveragedonthelastEditFieldLabel
            app.whenaveragedonthelastEditFieldLabel = uilabel(app.StopconditionsTab);
            app.whenaveragedonthelastEditFieldLabel.HorizontalAlignment = 'right';
            app.whenaveragedonthelastEditFieldLabel.Position = [650 592 146 22];
            app.whenaveragedonthelastEditFieldLabel.Text = 'when averaged on the last';

            % Create whenaveragedonthelastEditField
            app.whenaveragedonthelastEditField = uieditfield(app.StopconditionsTab, 'numeric');
            app.whenaveragedonthelastEditField.Limits = [0 Inf];
            app.whenaveragedonthelastEditField.RoundFractionalValues = 'on';
            app.whenaveragedonthelastEditField.Position = [811 592 53 22];
            app.whenaveragedonthelastEditField.Value = 50;

            % Create iterationsLabel
            app.iterationsLabel = uilabel(app.StopconditionsTab);
            app.iterationsLabel.Position = [874 592 54 22];
            app.iterationsLabel.Text = 'iterations';

            % Create thenLabel
            app.thenLabel = uilabel(app.StopconditionsTab);
            app.thenLabel.Position = [11 565 29 22];
            app.thenLabel.Text = 'then';

            % Create stopcondition_action_DropDown
            app.stopcondition_action_DropDown = uidropdown(app.StopconditionsTab);
            app.stopcondition_action_DropDown.Items = {'Ignore (no stoping condition)', 'Move to next phase or generation pass (and stop if last iteration)', 'Stop'};
            app.stopcondition_action_DropDown.Position = [54 565 465 22];
            app.stopcondition_action_DropDown.Value = 'Ignore (no stoping condition)';

            % Create AlgothmisfasterforLabel_11
            app.AlgothmisfasterforLabel_11 = uilabel(app.StopconditionsTab);
            app.AlgothmisfasterforLabel_11.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_11.Position = [11 538 755 22];
            app.AlgothmisfasterforLabel_11.Text = 'The threshold on the volume fraction rate and particle generation rate are applied on the maximum calculated among the different phases.';

            % Create MaximumwallclocktimesEditFieldLabel
            app.MaximumwallclocktimesEditFieldLabel = uilabel(app.StopconditionsTab);
            app.MaximumwallclocktimesEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumwallclocktimesEditFieldLabel.FontWeight = 'bold';
            app.MaximumwallclocktimesEditFieldLabel.Position = [11 496 164 22];
            app.MaximumwallclocktimesEditFieldLabel.Text = 'Maximum wallclock time (s)';

            % Create MaximumwallclocktimesEditField
            app.MaximumwallclocktimesEditField = uieditfield(app.StopconditionsTab, 'numeric');
            app.MaximumwallclocktimesEditField.Limits = [0 Inf];
            app.MaximumwallclocktimesEditField.FontWeight = 'bold';
            app.MaximumwallclocktimesEditField.Position = [190 496 64 22];
            app.MaximumwallclocktimesEditField.Value = 3600;

            % Create Label_8
            app.Label_8 = uilabel(app.StopconditionsTab);
            app.Label_8.FontAngle = 'italic';
            app.Label_8.Position = [264 496 680 22];
            app.Label_8.Text = 'Algorithm will always stop after reaching this generation time. Use high value if you want to ignore this last stoping condition.';

            % Create UpscaleoptionalTab
            app.UpscaleoptionalTab = uitab(app.TabGroup);
            app.UpscaleoptionalTab.Tooltip = {'Upscale (optional)'};
            app.UpscaleoptionalTab.Title = 'Upscale (optional)';
            app.UpscaleoptionalTab.BackgroundColor = [0.9412 0.9412 0.9412];
            app.UpscaleoptionalTab.ForegroundColor = [0.4941 0.1843 0.5569];

            % Create Instructions_upscaling
            app.Instructions_upscaling = uilabel(app.UpscaleoptionalTab);
            app.Instructions_upscaling.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_upscaling.HorizontalAlignment = 'center';
            app.Instructions_upscaling.FontWeight = 'bold';
            app.Instructions_upscaling.Position = [11 746 977 22];
            app.Instructions_upscaling.Text = 'Up(down) microstructure scaling';

            % Create SelectasolidphaseandsetthenumberofLabel_12
            app.SelectasolidphaseandsetthenumberofLabel_12 = uilabel(app.UpscaleoptionalTab);
            app.SelectasolidphaseandsetthenumberofLabel_12.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_12.Position = [11 703 826 28];
            app.SelectasolidphaseandsetthenumberofLabel_12.Text = {'Algorithm generation may take some time, especially if you are targeting a high density and/or a low (particle volume / domain volume) ratio.'; 'One way to speed up the process is to generate with a coarse representation and then to upscale.'};

            % Create SelectasolidphaseandsetthenumberofLabel_13
            app.SelectasolidphaseandsetthenumberofLabel_13 = uilabel(app.UpscaleoptionalTab);
            app.SelectasolidphaseandsetthenumberofLabel_13.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_13.Position = [11 580 938 22];
            app.SelectasolidphaseandsetthenumberofLabel_13.Text = 'Note that in-house upscale algorithm limits aliasing, so interface will be better represented compared with a simple MATLAB "imresize3" with "neareast" option';

            % Create Scalingfactor1upscaling1downscaling1noscalingEditFieldLabel
            app.Scalingfactor1upscaling1downscaling1noscalingEditFieldLabel = uilabel(app.UpscaleoptionalTab);
            app.Scalingfactor1upscaling1downscaling1noscalingEditFieldLabel.HorizontalAlignment = 'right';
            app.Scalingfactor1upscaling1downscaling1noscalingEditFieldLabel.Position = [11 676 336 22];
            app.Scalingfactor1upscaling1downscaling1noscalingEditFieldLabel.Text = 'Scaling factor (>1: upscaling, <1: downscaling, =1 no scaling)';

            % Create Scalingfactor1upscaling1downscaling1noscalingEditField
            app.Scalingfactor1upscaling1downscaling1noscalingEditField = uieditfield(app.UpscaleoptionalTab, 'numeric');
            app.Scalingfactor1upscaling1downscaling1noscalingEditField.Position = [362 676 46 22];
            app.Scalingfactor1upscaling1downscaling1noscalingEditField.Value = 1;

            % Create Image2
            app.Image2 = uiimage(app.UpscaleoptionalTab);
            app.Image2.Position = [11 25 974 550];
            app.Image2.ImageSource = 'Scaling.png';

            % Create Label_7
            app.Label_7 = uilabel(app.UpscaleoptionalTab);
            app.Label_7.FontAngle = 'italic';
            app.Label_7.Position = [423 676 423 22];
            app.Label_7.Text = 'Number of voxel (see "volume fraction tab") is multiplied by the scaling factor.';

            % Create ApplyscalingtophaselabelfastCheckBox
            app.ApplyscalingtophaselabelfastCheckBox = uicheckbox(app.UpscaleoptionalTab);
            app.ApplyscalingtophaselabelfastCheckBox.Enable = 'off';
            app.ApplyscalingtophaselabelfastCheckBox.Text = 'Apply scaling to phase label (fast)';
            app.ApplyscalingtophaselabelfastCheckBox.Position = [11 649 202 22];
            app.ApplyscalingtophaselabelfastCheckBox.Value = true;

            % Create ApplyscalingalsotoparticlelabelslowCheckBox
            app.ApplyscalingalsotoparticlelabelslowCheckBox = uicheckbox(app.UpscaleoptionalTab);
            app.ApplyscalingalsotoparticlelabelslowCheckBox.Text = 'Apply scaling also to particle label (slow)';
            app.ApplyscalingalsotoparticlelabelslowCheckBox.Position = [11 622 239 22];

            % Create VerificationandCharacterizationTab
            app.VerificationandCharacterizationTab = uitab(app.TabGroup);
            app.VerificationandCharacterizationTab.Tooltip = {'Verification and Characterization'};
            app.VerificationandCharacterizationTab.Title = 'Verification and Characterization';
            app.VerificationandCharacterizationTab.ForegroundColor = [0.4941 0.1843 0.5569];

            % Create Instructions_characterization
            app.Instructions_characterization = uilabel(app.VerificationandCharacterizationTab);
            app.Instructions_characterization.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_characterization.HorizontalAlignment = 'center';
            app.Instructions_characterization.FontWeight = 'bold';
            app.Instructions_characterization.Position = [11 746 977 22];
            app.Instructions_characterization.Text = 'Verification and basic microstructure characterization';

            % Create SelectasolidphaseandsetthenumberofLabel_21
            app.SelectasolidphaseandsetthenumberofLabel_21 = uilabel(app.VerificationandCharacterizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_21.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_21.Position = [11 628 614 28];
            app.SelectasolidphaseandsetthenumberofLabel_21.Text = {'You can perform simple characterization to quickly check the generated volume.'; 'Note that a more in-depth characterization can be perfomed using the MATBOX characterization module.'};

            % Create VolumefractionsCheckBox
            app.VolumefractionsCheckBox = uicheckbox(app.VerificationandCharacterizationTab);
            app.VolumefractionsCheckBox.Enable = 'off';
            app.VolumefractionsCheckBox.Text = 'Volume fractions';
            app.VolumefractionsCheckBox.Position = [11 601 111 22];
            app.VolumefractionsCheckBox.Value = true;

            % Create ConnectivityslowforlargevolumesCheckBox
            app.ConnectivityslowforlargevolumesCheckBox = uicheckbox(app.VerificationandCharacterizationTab);
            app.ConnectivityslowforlargevolumesCheckBox.Text = 'Connectivity (slow for large volumes)';
            app.ConnectivityslowforlargevolumesCheckBox.Position = [11 574 219 22];

            % Create PoretortuosityfactorslowforlargevolumesCheckBox
            app.PoretortuosityfactorslowforlargevolumesCheckBox = uicheckbox(app.VerificationandCharacterizationTab);
            app.PoretortuosityfactorslowforlargevolumesCheckBox.Text = 'Pore tortuosity factor (normalized effective diffusion coefficient = porosity / tortuosity factor)';
            app.PoretortuosityfactorslowforlargevolumesCheckBox.Position = [11 547 513 22];

            % Create CharacterizeDropDownLabel
            app.CharacterizeDropDownLabel = uilabel(app.VerificationandCharacterizationTab);
            app.CharacterizeDropDownLabel.HorizontalAlignment = 'right';
            app.CharacterizeDropDownLabel.Position = [11 520 77 22];
            app.CharacterizeDropDownLabel.Text = 'Characterize ';

            % Create CharacterizeDropDown
            app.CharacterizeDropDown = uidropdown(app.VerificationandCharacterizationTab);
            app.CharacterizeDropDown.Items = {'volume before upscaling', 'volume after upscaling'};
            app.CharacterizeDropDown.Position = [103 520 168 22];
            app.CharacterizeDropDown.Value = 'volume before upscaling';

            % Create SelectasolidphaseandsetthenumberofLabel_22
            app.SelectasolidphaseandsetthenumberofLabel_22 = uilabel(app.VerificationandCharacterizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_22.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_22.Position = [11 703 984 28];
            app.SelectasolidphaseandsetthenumberofLabel_22.Text = {'Inputs (volume fractions, particle diameter and rotation) can be compared with the obtained geometries. Note that no post-processing or characterization is performed:'; 'we simply use the properties assigned to each phase during the generation process.'};

            % Create CompareoutputsversusinputsCheckBox
            app.CompareoutputsversusinputsCheckBox = uicheckbox(app.VerificationandCharacterizationTab);
            app.CompareoutputsversusinputsCheckBox.Text = 'Compare outputs versus inputs';
            app.CompareoutputsversusinputsCheckBox.Position = [11 676 188 22];
            app.CompareoutputsversusinputsCheckBox.Value = true;

            % Create GenerationTab
            app.GenerationTab = uitab(app.TabGroup);
            app.GenerationTab.Title = 'Generation';
            app.GenerationTab.ForegroundColor = [1 0 0];

            % Create Generate_Button
            app.Generate_Button = uibutton(app.GenerationTab, 'push');
            app.Generate_Button.ButtonPushedFcn = createCallbackFcn(app, @Generate_ButtonPushed, true);
            app.Generate_Button.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Generate_Button.FontSize = 16;
            app.Generate_Button.FontWeight = 'bold';
            app.Generate_Button.Enable = 'off';
            app.Generate_Button.Position = [554 641 211 40];
            app.Generate_Button.Text = 'Run generation algorithm';

            % Create Instructions_generation
            app.Instructions_generation = uilabel(app.GenerationTab);
            app.Instructions_generation.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_generation.HorizontalAlignment = 'center';
            app.Instructions_generation.FontWeight = 'bold';
            app.Instructions_generation.Position = [11 746 977 22];
            app.Instructions_generation.Text = 'Once all parameters have been set, run the generation algorithm';

            % Create SelectasolidphaseandsetthenumberofLabel_10
            app.SelectasolidphaseandsetthenumberofLabel_10 = uilabel(app.GenerationTab);
            app.SelectasolidphaseandsetthenumberofLabel_10.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_10.Position = [11 709 531 22];
            app.SelectasolidphaseandsetthenumberofLabel_10.Text = '1) Set number of runs. Algorithm will be run sequentially. Valuable to check reproducibility';

            % Create NumberofrunsEditFieldLabel
            app.NumberofrunsEditFieldLabel = uilabel(app.GenerationTab);
            app.NumberofrunsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofrunsEditFieldLabel.Position = [11 682 88 22];
            app.NumberofrunsEditFieldLabel.Text = 'Number of runs';

            % Create NumberofrunsEditField
            app.NumberofrunsEditField = uieditfield(app.GenerationTab, 'numeric');
            app.NumberofrunsEditField.Limits = [1 Inf];
            app.NumberofrunsEditField.RoundFractionalValues = 'on';
            app.NumberofrunsEditField.Position = [114 682 37 22];
            app.NumberofrunsEditField.Value = 1;

            % Create SelectasolidphaseandsetthenumberofLabel_11
            app.SelectasolidphaseandsetthenumberofLabel_11 = uilabel(app.GenerationTab);
            app.SelectasolidphaseandsetthenumberofLabel_11.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_11.Position = [11 650 538 22];
            app.SelectasolidphaseandsetthenumberofLabel_11.Text = '2) Click to run the algorithms. Tables below will be updated each time a run is finished.';

            % Create Connectivity_UITable
            app.Connectivity_UITable = uitable(app.GenerationTab);
            app.Connectivity_UITable.ColumnName = {'Run#'; 'Pores'; 'Union of solid phases'};
            app.Connectivity_UITable.RowName = {};
            app.Connectivity_UITable.Position = [11 44 470 260];

            % Create AlgothmisfasterforLabel_12
            app.AlgothmisfasterforLabel_12 = uilabel(app.GenerationTab);
            app.AlgothmisfasterforLabel_12.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_12.Position = [11 617 660 22];
            app.AlgothmisfasterforLabel_12.Text = 'Volume fractions';

            % Create AlgothmisfasterforLabel_13
            app.AlgothmisfasterforLabel_13 = uilabel(app.GenerationTab);
            app.AlgothmisfasterforLabel_13.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_13.Position = [8 309 470 28];
            app.AlgothmisfasterforLabel_13.Text = {'Connectivity (largest cluster vol%. For more in-depth connectivity analysis,'; 'please use the characterization module.'};

            % Create Volume_fraction_UITable
            app.Volume_fraction_UITable = uitable(app.GenerationTab);
            app.Volume_fraction_UITable.ColumnName = {'Run#'; 'Pores'; 'Phase 1'};
            app.Volume_fraction_UITable.RowName = {};
            app.Volume_fraction_UITable.Position = [11 352 563 260];

            % Create AlgothmisfasterforLabel_14
            app.AlgothmisfasterforLabel_14 = uilabel(app.GenerationTab);
            app.AlgothmisfasterforLabel_14.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_14.Position = [518 312 470 22];
            app.AlgothmisfasterforLabel_14.Text = 'Pore tortuosity factor';

            % Create Tortuosity_UITable
            app.Tortuosity_UITable = uitable(app.GenerationTab);
            app.Tortuosity_UITable.ColumnName = {'Run#'; 'Direction 1'; 'Direction 2'; 'Direction 3 (thickness)'};
            app.Tortuosity_UITable.RowName = {};
            app.Tortuosity_UITable.Position = [518 44 470 260];

            % Create AlgothmisfasterforLabel_15
            app.AlgothmisfasterforLabel_15 = uilabel(app.GenerationTab);
            app.AlgothmisfasterforLabel_15.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_15.Position = [17 12 977 22];
            app.AlgothmisfasterforLabel_15.Text = 'For more in-depth characterization, please use the characterization module.';

            % Create AlgothmisfasterforLabel_16
            app.AlgothmisfasterforLabel_16 = uilabel(app.GenerationTab);
            app.AlgothmisfasterforLabel_16.FontAngle = 'italic';
            app.AlgothmisfasterforLabel_16.Position = [817 618 172 22];
            app.AlgothmisfasterforLabel_16.Text = 'Outcome and wallclock time (s)';

            % Create outcometime_UITable
            app.outcometime_UITable = uitable(app.GenerationTab);
            app.outcometime_UITable.ColumnName = {'Outcome'; 'Generation'; 'Scaling'; 'Characterization'};
            app.outcometime_UITable.RowName = {};
            app.outcometime_UITable.Position = [611 352 377 260];

            % Create VisualizationTab
            app.VisualizationTab = uitab(app.TabGroup);
            app.VisualizationTab.Title = 'Visualization';

            % Create Instructions_generation_3
            app.Instructions_generation_3 = uilabel(app.VisualizationTab);
            app.Instructions_generation_3.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_generation_3.HorizontalAlignment = 'center';
            app.Instructions_generation_3.FontWeight = 'bold';
            app.Instructions_generation_3.Position = [11 746 977 22];
            app.Instructions_generation_3.Text = 'Visualize as-generated microstructures';

            % Create Visualization_UITable
            app.Visualization_UITable = uitable(app.VisualizationTab);
            app.Visualization_UITable.ColumnName = {'Run#'; 'Scaled'; 'Visualize'};
            app.Visualization_UITable.RowName = {};
            app.Visualization_UITable.ColumnEditable = [false false true];
            app.Visualization_UITable.CellEditCallback = createCallbackFcn(app, @Visualization_UITableCellEdit, true);
            app.Visualization_UITable.Position = [11 438 355 266];

            % Create SelectasolidphaseandsetthenumberofLabel_23
            app.SelectasolidphaseandsetthenumberofLabel_23 = uilabel(app.VisualizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_23.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_23.Position = [11 709 441 22];
            app.SelectasolidphaseandsetthenumberofLabel_23.Text = '1) Select what volumes you would like to visualize (one or several at once).';

            % Create Plot_phase_label
            app.Plot_phase_label = uibutton(app.VisualizationTab, 'push');
            app.Plot_phase_label.ButtonPushedFcn = createCallbackFcn(app, @Plot_phase_labelButtonPushed, true);
            app.Plot_phase_label.BackgroundColor = [0.0745 0.6235 1];
            app.Plot_phase_label.FontSize = 14;
            app.Plot_phase_label.FontWeight = 'bold';
            app.Plot_phase_label.FontColor = [1 1 1];
            app.Plot_phase_label.Enable = 'off';
            app.Plot_phase_label.Position = [82 368 122 24];
            app.Plot_phase_label.Text = 'Plot phase label';

            % Create Plot_particle_label
            app.Plot_particle_label = uibutton(app.VisualizationTab, 'push');
            app.Plot_particle_label.ButtonPushedFcn = createCallbackFcn(app, @Plot_particle_labelButtonPushed, true);
            app.Plot_particle_label.BackgroundColor = [0.0745 0.6235 1];
            app.Plot_particle_label.FontSize = 14;
            app.Plot_particle_label.FontWeight = 'bold';
            app.Plot_particle_label.FontColor = [1 1 1];
            app.Plot_particle_label.Enable = 'off';
            app.Plot_particle_label.Position = [209 368 131 24];
            app.Plot_particle_label.Text = 'Plot particle label';

            % Create SelectasolidphaseandsetthenumberofLabel_24
            app.SelectasolidphaseandsetthenumberofLabel_24 = uilabel(app.VisualizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_24.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_24.Position = [11 396 355 22];
            app.SelectasolidphaseandsetthenumberofLabel_24.Text = '2) Then choose to plot either phase label or particle label';

            % Create SelectasolidphaseandsetthenumberofLabel_25
            app.SelectasolidphaseandsetthenumberofLabel_25 = uilabel(app.VisualizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_25.FontAngle = 'italic';
            app.SelectasolidphaseandsetthenumberofLabel_25.Position = [345 369 355 22];
            app.SelectasolidphaseandsetthenumberofLabel_25.Text = 'Note: functions from the visualization module are called.';

            % Create DslicesLabel
            app.DslicesLabel = uilabel(app.VisualizationTab);
            app.DslicesLabel.Position = [21 369 58 22];
            app.DslicesLabel.Text = '2D slices:';

            % Create DviewLabel
            app.DviewLabel = uilabel(app.VisualizationTab);
            app.DviewLabel.Position = [27 337 52 22];
            app.DviewLabel.Text = '3D view:';

            % Create Plot_phase_label_3D
            app.Plot_phase_label_3D = uibutton(app.VisualizationTab, 'push');
            app.Plot_phase_label_3D.ButtonPushedFcn = createCallbackFcn(app, @Plot_phase_label_3DButtonPushed, true);
            app.Plot_phase_label_3D.BackgroundColor = [0.0745 0.6235 1];
            app.Plot_phase_label_3D.FontSize = 14;
            app.Plot_phase_label_3D.FontWeight = 'bold';
            app.Plot_phase_label_3D.FontColor = [1 1 1];
            app.Plot_phase_label_3D.Enable = 'off';
            app.Plot_phase_label_3D.Position = [82 336 122 24];
            app.Plot_phase_label_3D.Text = 'Plot phase label';

            % Create Plot_particle_label_3D
            app.Plot_particle_label_3D = uibutton(app.VisualizationTab, 'push');
            app.Plot_particle_label_3D.ButtonPushedFcn = createCallbackFcn(app, @Plot_particle_label_3DButtonPushed, true);
            app.Plot_particle_label_3D.BackgroundColor = [0.0745 0.6235 1];
            app.Plot_particle_label_3D.FontSize = 14;
            app.Plot_particle_label_3D.FontWeight = 'bold';
            app.Plot_particle_label_3D.FontColor = [1 1 1];
            app.Plot_particle_label_3D.Enable = 'off';
            app.Plot_particle_label_3D.Position = [209 336 131 24];
            app.Plot_particle_label_3D.Text = 'Plot particle label';

            % Create SelectasolidphaseandsetthenumberofLabel_26
            app.SelectasolidphaseandsetthenumberofLabel_26 = uilabel(app.VisualizationTab);
            app.SelectasolidphaseandsetthenumberofLabel_26.FontAngle = 'italic';
            app.SelectasolidphaseandsetthenumberofLabel_26.Position = [345 337 355 22];
            app.SelectasolidphaseandsetthenumberofLabel_26.Text = 'Note: 3D view is available only when one volume is selected';

            % Create DviewcolormapDropDownLabel
            app.DviewcolormapDropDownLabel = uilabel(app.VisualizationTab);
            app.DviewcolormapDropDownLabel.HorizontalAlignment = 'right';
            app.DviewcolormapDropDownLabel.Position = [82 309 101 22];
            app.DviewcolormapDropDownLabel.Text = '3D view colormap';

            % Create DviewcolormapDropDown
            app.DviewcolormapDropDown = uidropdown(app.VisualizationTab);
            app.DviewcolormapDropDown.Items = {'parula', 'turbo', 'jet', 'gray', 'bone', 'copper'};
            app.DviewcolormapDropDown.Position = [198 309 100 22];
            app.DviewcolormapDropDown.Value = 'parula';

            % Create HowtoquoteTab
            app.HowtoquoteTab = uitab(app.TabGroup);
            app.HowtoquoteTab.Title = 'How to quote';

            % Create Instructions_generation_2
            app.Instructions_generation_2 = uilabel(app.HowtoquoteTab);
            app.Instructions_generation_2.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_generation_2.HorizontalAlignment = 'center';
            app.Instructions_generation_2.FontWeight = 'bold';
            app.Instructions_generation_2.Position = [11 746 977 22];
            app.Instructions_generation_2.Text = 'How to quote';

            % Create SelectasolidphaseandsetthenumberofLabel_20
            app.SelectasolidphaseandsetthenumberofLabel_20 = uilabel(app.HowtoquoteTab);
            app.SelectasolidphaseandsetthenumberofLabel_20.FontWeight = 'bold';
            app.SelectasolidphaseandsetthenumberofLabel_20.Position = [11 709 976 22];
            app.SelectasolidphaseandsetthenumberofLabel_20.Text = 'If you use this algorithm or report results produced thanks to this algorithm, please quote it as follow:';

            % Create TextArea
            app.TextArea = uitextarea(app.HowtoquoteTab);
            app.TextArea.Editable = 'off';
            app.TextArea.Position = [11 644 978 60];
            app.TextArea.Value = {'F. L. E. Usseglio-Viretta, P. Patel, E. Bernhardt, A Mistry, P. P. Mukherjee, J. Allen, S. J. Cooper, and K. Smith'; 'MATBOX: An Open-source Microstructure Analysis Toolbox for microstructure generation, segmentation, characterization, visualization, correlation, and meshing'; 'SoftwareX, in preparation'};

            % Show the figure after all components are created
            app.EllipsoidbasedstochasticgenerationmoduleUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Microstructure_generation_stochastic_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.EllipsoidbasedstochasticgenerationmoduleUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.EllipsoidbasedstochasticgenerationmoduleUIFigure)
        end
    end
end