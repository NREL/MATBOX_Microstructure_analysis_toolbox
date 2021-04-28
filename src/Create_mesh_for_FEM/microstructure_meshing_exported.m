classdef microstructure_meshing_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        MeshingmoduleUIFigure           matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        InstructionsTab                 matlab.ui.container.Tab
        Instructions_title              matlab.ui.control.Label
        Instructions_Label_3            matlab.ui.control.TextArea
        Instructions_Label_1            matlab.ui.control.Label
        Instructions_Label_2            matlab.ui.control.Label
        Instructions_Label_4            matlab.ui.control.Label
        Instructions_Label_5            matlab.ui.control.TextArea
        Instructions_Label_6            matlab.ui.control.Label
        Instructions_Label_7            matlab.ui.control.TextArea
        FolderandsaveoptionsTab         matlab.ui.container.Tab
        Saveoptions_title               matlab.ui.control.Label
        Save_Label                      matlab.ui.control.Label
        ClicktoselectsavefolderButton   matlab.ui.control.Button
        Save_folder_text                matlab.ui.control.Label
        TabisNOTsetupcorrectlyLabel     matlab.ui.control.Label
        Folder_Lamp                     matlab.ui.control.Lamp
        Setup_UITable                   matlab.ui.control.Table
        SelectmicrostructuresTab        matlab.ui.container.Tab
        SelectMicrostructure_Title      matlab.ui.control.Label
        OneuniquemicrostructureLabel    matlab.ui.control.Label
        Halfcellorfullcellselectatleast2volumesLabel  matlab.ui.control.Label
        PolycrystallinearchitectureBetaLabel  matlab.ui.control.Label
        ImporttifButton_Unique          matlab.ui.control.Button
        ImporttifButton_Cell            matlab.ui.control.Button
        ImporttifButton_Poly            matlab.ui.control.Button
        LeftelectrodeLabel              matlab.ui.control.Label
        SeparatorLabel                  matlab.ui.control.Label
        RightelectrodeLabel             matlab.ui.control.Label
        ChoosedomainDropDownLabel       matlab.ui.control.Label
        ChoosedomainDropDown            matlab.ui.control.DropDown
        HomogenousButton_Cell           matlab.ui.control.Button
        orLabel                         matlab.ui.control.Label
        LengthinvoxelEditFieldLabel     matlab.ui.control.Label
        LengthinvoxelEditField          matlab.ui.control.NumericEditField
        SavemicrostructureButton        matlab.ui.control.Button
        Label                           matlab.ui.control.Label
        TiffileloadedLabel              matlab.ui.control.Label
        Microstructure_vf_UITable       matlab.ui.control.Table
        VolumefractionsLabel            matlab.ui.control.Label
        Select_VisualizemicrostructureButton  matlab.ui.control.Button
        Select_Lamp                     matlab.ui.control.Lamp
        TabisNOTsetupcorrectlyLabel_2   matlab.ui.control.Label
        Select_ROI_UITable              matlab.ui.control.Table
        SelecttheregionofinterestROILabel  matlab.ui.control.Label
        SelecttheregionofinterestROIoptionalLabel_2  matlab.ui.control.Label
        Select_Orientation_UITable      matlab.ui.control.Table
        Select_Rescale_UITable          matlab.ui.control.Table
        SelecttheregionofinterestROIoptionalLabel_3  matlab.ui.control.Label
        Rows12unswapedButton            matlab.ui.control.StateButton
        Rows13unswapedButton            matlab.ui.control.StateButton
        Rows23unswapedButton            matlab.ui.control.StateButton
        Row1unflippedButton             matlab.ui.control.StateButton
        Row2unflippedButton             matlab.ui.control.StateButton
        Row3unflippedButton             matlab.ui.control.StateButton
        ScalingfactorEditFieldLabel     matlab.ui.control.Label
        ScalingfactorEditField          matlab.ui.control.NumericEditField
        VoxelsizebeforeEditFieldLabel   matlab.ui.control.Label
        VoxelsizebeforeEditField        matlab.ui.control.NumericEditField
        VoxelsizeafterEditFieldLabel    matlab.ui.control.Label
        VoxelsizeafterEditField         matlab.ui.control.NumericEditField
        ApplyscalingButton              matlab.ui.control.Button
        BackgroundidEditFieldLabel      matlab.ui.control.Label
        BackgroundidEditField           matlab.ui.control.NumericEditField
        CC_left                         matlab.ui.control.Label
        CC_right                        matlab.ui.control.Label
        DimensioncompatibilityTab       matlab.ui.container.Tab
        Compatibility_Title             matlab.ui.control.Label
        Comp_voxelsize_text             matlab.ui.control.Label
        Comp_Lamp                       matlab.ui.control.Lamp
        TabisNOTsetupcorrectlyLabel_3   matlab.ui.control.Label
        Label_2                         matlab.ui.control.Label
        Comp_voxelsize_text_2           matlab.ui.control.Label
        Dim_TP_UITable                  matlab.ui.control.Table
        Comp_voxelsize_text_3           matlab.ui.control.Label
        VoxelsizeEditFieldLabel         matlab.ui.control.Label
        VoxelsizeEditField              matlab.ui.control.NumericEditField
        ActualizeresetdimensionsButton  matlab.ui.control.Button
        Actualizedim_text               matlab.ui.control.Label
        Comp_voxelsize_text_4           matlab.ui.control.Label
        Dim_IP1_UITable                 matlab.ui.control.Table
        Comp_voxelsize_text_5           matlab.ui.control.Label
        Dim_IP2_UITable                 matlab.ui.control.Table
        Dim_Autocrop_1_Button           matlab.ui.control.Button
        Dim_Autocrop_2_Button           matlab.ui.control.Button
        Label_3                         matlab.ui.control.Label
        SavecellButton                  matlab.ui.control.Button
        Select_VisualizecellButton      matlab.ui.control.Button
        Inplanedimension1arematchingLabel  matlab.ui.control.Label
        Inplanedimension2arematchingLabel  matlab.ui.control.Label
        Dim_YoucanskipthistabLabel      matlab.ui.control.Label
        micrometersLabel                matlab.ui.control.Label
        MorphologyopeningTab            matlab.ui.container.Tab
        PerformalloperationssequentiallyLabel  matlab.ui.control.Label
        MO_Title                        matlab.ui.control.Label
        MO_text_1                       matlab.ui.control.Label
        MO_SavecellButton               matlab.ui.control.Button
        MO_text_2                       matlab.ui.control.Label
        MO_VisualizecellButton          matlab.ui.control.Button
        MO_vf_UITable                   matlab.ui.control.Table
        MO_VolumefractionsLabel         matlab.ui.control.Label
        MO_operations_UITable           matlab.ui.control.Table
        MO_VolumefractionsLabel_2       matlab.ui.control.Label
        MO_text_4                       matlab.ui.control.Label
        MO_text_3                       matlab.ui.control.Label
        ErosiondistanceinvoxellengthEditFieldLabel  matlab.ui.control.Label
        ErosiondistanceinvoxellengthEditField  matlab.ui.control.NumericEditField
        DilatationdistanceinvoxellengthEditFieldLabel  matlab.ui.control.Label
        DilatationdistanceinvoxellengthEditField  matlab.ui.control.NumericEditField
        NumberofiterationEditFieldLabel  matlab.ui.control.Label
        ErosiondilatationNumberofiterationEditField  matlab.ui.control.NumericEditField
        PerformDropDownLabel            matlab.ui.control.Label
        PerformErosionDropDown          matlab.ui.control.DropDown
        ErositiondilatationisalwaysperformedfirstLabel  matlab.ui.control.Label
        MO_text_5                       matlab.ui.control.Label
        DetectandremoveverticeverticeandedgeedgeconnectionsLabel  matlab.ui.control.Label
        PerformDropDown_2Label          matlab.ui.control.Label
        PerformCleanConnectionDropDown  matlab.ui.control.DropDown
        MO_text_6                       matlab.ui.control.Label
        GrainvolumefilterinnumberofvoxelEditFieldLabel  matlab.ui.control.Label
        GrainvolumefilterinnumberofvoxelEditField  matlab.ui.control.NumericEditField
        MO_DoButton                     matlab.ui.control.Button
        MO_UndoButton                   matlab.ui.control.Button
        MO_SaveButton                   matlab.ui.control.Button
        MO_statut_Label                 matlab.ui.control.Label
        MO_text_7                       matlab.ui.control.Label
        MO_PreservebackgroundCheckBox   matlab.ui.control.CheckBox
        MO_UnionsoliduniqueclusterCheckBox  matlab.ui.control.CheckBox
        MO_PoreuniqueclusterCheckBox    matlab.ui.control.CheckBox
        MO_SkipButton                   matlab.ui.control.Button
        MO_UndoAllButton                matlab.ui.control.Button
        CheckgraincontiguityandfixitifdiscontinousslowCheckBox  matlab.ui.control.CheckBox
        AssemblecellTab                 matlab.ui.container.Tab
        Assemble_Title                  matlab.ui.control.Label
        Assemble_Instructions           matlab.ui.control.Label
        Assemble_Lamp                   matlab.ui.control.Lamp
        TabisNOTsetupcorrectlyLabel_4   matlab.ui.control.Label
        Assemble_UITable                matlab.ui.control.Table
        Assemble_Label_4                matlab.ui.control.Label
        Assemble_SavecellButton         matlab.ui.control.Button
        Assemble_VisualizecellButton    matlab.ui.control.Button
        Label_4                         matlab.ui.control.Label
        UpdateButton                    matlab.ui.control.Button
        UpdateButton_label              matlab.ui.control.Label
        Group_UITable                   matlab.ui.control.Table
        UpdateButton_label_2            matlab.ui.control.Label
        Meshgenerationchoice            matlab.ui.container.Tab
        Meshgen_Title                   matlab.ui.control.Label
        Meshgen_Iso2mesh                matlab.ui.control.Button
        Meshgen_Cuboidmesh              matlab.ui.control.Button
        Label_5                         matlab.ui.control.Label
        Label_6                         matlab.ui.control.Label
        Label_7                         matlab.ui.control.Label
        Iso2meshoptionsTab              matlab.ui.container.Tab
        Iso2mesh_Title                  matlab.ui.control.Label
        CGAL_UITable                    matlab.ui.control.Table
        Iso2mesh_label                  matlab.ui.control.Label
        Iso2mesh_label_2                matlab.ui.control.Label
        MethodDropDownLabel             matlab.ui.control.Label
        Surfacemesh_MethodDropDown      matlab.ui.control.DropDown
        Smoothing_UITable               matlab.ui.control.Table
        Iso2mesh_label_3                matlab.ui.control.Label
        Meshdensity_UITable             matlab.ui.control.Table
        Label_12                        matlab.ui.control.Label
        Iso2meshversionnumberLabel      matlab.ui.control.Label
        StructuredmeshoptionsTab        matlab.ui.container.Tab
        Regularmesh_Title               matlab.ui.control.Label
        EachvoxelwillberepresentedwithDropDownLabel  matlab.ui.control.Label
        EachvoxelwillberepresentedwithDropDown  matlab.ui.control.DropDown
        Regularmesh_Image5              matlab.ui.control.Image
        Regularmesh_Image24             matlab.ui.control.Image
        Regularmesh_Image6              matlab.ui.control.Image
        CreatemeshTab                   matlab.ui.container.Tab
        Createmesh_Title                matlab.ui.control.Label
        Createmesh_button               matlab.ui.control.Button
        Meshchoice_UITable              matlab.ui.control.Table
        Label_8                         matlab.ui.control.Label
        Mesh_VisualizeButton            matlab.ui.control.Button
        Cellindexstartsat0CheckBox      matlab.ui.control.CheckBox
        Createmesh_choicemesh           matlab.ui.control.Label
        Mesh_SaveButton                 matlab.ui.control.Button
        StatutTextAreaLabel             matlab.ui.control.Label
        Createmesh_StatutTextArea       matlab.ui.control.TextArea
        CreateMesh_text_1               matlab.ui.control.Label
        LightingDropDownLabel           matlab.ui.control.Label
        LightingDropDown                matlab.ui.control.DropDown
        MeshedgesDropDownLabel          matlab.ui.control.Label
        MeshedgesDropDown               matlab.ui.control.DropDown
        ColormapDropDownLabel           matlab.ui.control.Label
        ColormapDropDown                matlab.ui.control.DropDown
        Transparency01EditFieldLabel    matlab.ui.control.Label
        Transparency01EditField         matlab.ui.control.NumericEditField
        ShowDropDownLabel               matlab.ui.control.Label
        ShowDropDown                    matlab.ui.control.DropDown
        SavepngCheckBox                 matlab.ui.control.CheckBox
        SavefigCheckBox                 matlab.ui.control.CheckBox
        Label_9                         matlab.ui.control.Label
        GridCheckBox                    matlab.ui.control.CheckBox
        AxislabelCheckBox               matlab.ui.control.CheckBox
        Meshformat_UITable              matlab.ui.control.Table
        CreateMesh_text_2               matlab.ui.control.Label
        csvandmatLabel                  matlab.ui.control.Label
        Face_CheckBox                   matlab.ui.control.CheckBox
        ColorisDropDownLabel            matlab.ui.control.Label
        ColorisDropDown                 matlab.ui.control.DropDown
        CalculatemeshqualityTab         matlab.ui.container.Tab
        Meshquality_Title               matlab.ui.control.Label
        Label_11                        matlab.ui.control.Label
        Meshquality_JoeLiu_button       matlab.ui.control.Button
        AboutTab                        matlab.ui.container.Tab
        About_title                     matlab.ui.control.Label
        About_TextArea                  matlab.ui.control.TextArea
        QuotationinstructionsTextAreaLabel  matlab.ui.control.Label
        About_Quotationinstructions     matlab.ui.control.TextArea
        About_Logo_NREL                 matlab.ui.control.Image
        OpendocumentationButton         matlab.ui.control.Button
        Iso2MeshwebsiteButton           matlab.ui.control.Button
        Label_10                        matlab.ui.control.Label
    end

    
    properties (Access = private)
        % Default saving options
        mainsavefoldersave = []; savefolder_is_defined = false;
        save_mat = true;
        save_csv = false;
        save_msh = false;
        save_inp = false;
        save_stl = false;
        save_binarystl = false;
        save_parameters = true;
        save_loaded_tif = true;
        save_modified_tif = true;
                
        % Microstructure type
        micro_unique = false;
        micro_cell = false;
        micro_poly = false;
        
        % Format
        datatype_array = 'uint8';
        
        % Domain
        imported_domain = [0 0 0 0 0];
        Choosedomain = [];
        Domain =[];
        loaded_array = [];
        beforescaling_array = [];
        current_array = [];
        roi_sav = [];
        % Domain in cell
        Domain_in_cell = [];
        cellarray_firstassemble = []; % Cell after dimension (for half or full cell), microstructure for unique or polychristalline choice
        cellarray_secondassemble = []; % Cell afer morphology opening
        MO_Microstructure = []; MO_Microstructure_saved = []; % Cell during morphology opening
        cell_before_reassignemnt = []; % Cell after morphology opening
        cell_after_reassignement = []; % Cell after assemble
        cell_group_after_reassignement = []; % Cell group id after assemble
        Microstructure_dimension = []; % Voxel size; number of voxe along axe 1; number of voxe along axe 2; number of voxe along axe 3;
        Separator_bounds = []; % X min and max coordinates of separator if any
        
        domainbound=[]; % Domain location
        % Morphology opening
        n_mo_step = 3;
        morphology_operation=[];
        
        % Mesh options
        meshmethod = [];
        do_one_mesh_per_phase = false; do_one_mesh_per_group = true; do_one_mesh_wholegeometry = true;
        indexstart_zero = [];
        meshoptions = [];
                
        % Mesh results
        meshgroup = [];
        meshphase = [];
        meshvolume = [];
        
    end

    
    methods (Access = private)
        
        function [array, outcome] = function_import_microstructure(app)
            app.LengthinvoxelEditField.Visible = 'off'; app.LengthinvoxelEditFieldLabel.Visible = 'off';
            [FileName,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)'},'Select volume');
            app.Select_Lamp.Color = 'r'; app.ActualizeresetdimensionsButton.Enable = 'off';
            app.TabisNOTsetupcorrectlyLabel_2.Text = 'Tab is NOT setup correctly';            
            if FileName==0
                % User clicked cancel button or closed the dialog box
                app.reset_microstructure
                array = []; outcome = false;
            else
                app.TiffileloadedLabel.Text = 'Loading... please wait'; app.TiffileloadedLabel.Visible = 'on'; pause(0.01);
                [array,outcome] = function_load_tif([PathName FileName],app.datatype_array);
                if outcome % Success
                    app.TiffileloadedLabel.Text = [PathName FileName];
                else
                    app.reset_microstructure
                end
            end
            
        end
        
        function [] = reset_microstructure(app)
            app.OneuniquemicrostructureLabel.BackgroundColor = [0.90 0.90 0.90];
            app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.90 0.90 0.90];
            app.CC_left.BackgroundColor = [0.90 0.90 0.90];
            app.LeftelectrodeLabel.BackgroundColor = [0.90 0.90 0.90];
            app.SeparatorLabel.BackgroundColor = [0.90 0.90 0.90];
            app.RightelectrodeLabel.BackgroundColor = [0.90 0.90 0.90];
            app.CC_right.BackgroundColor = [0.90 0.90 0.90];
            app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.90 0.90 0.90];
            app.micro_unique = false; app.micro_cell = false; app.micro_poly = false;
            app.imported_domain = [0 0 0 0 0];
            app.Select_Lamp.Color = 'r'; app.ActualizeresetdimensionsButton.Enable = 'off';
            app.TabisNOTsetupcorrectlyLabel_2.Text = 'Tab is NOT setup correctly';
            app.Domain(1).array = []; app.Domain(2).array = []; app.Domain(3).array = []; app.Domain(4).array = []; app.Domain(5).array = []; 
            app.TiffileloadedLabel.Visible = 'off'; 
            app.Select_VisualizemicrostructureButton.Enable = 'off'; 
            app.SelecttheregionofinterestROILabel.Visible = 'off'; 
            app.Select_ROI_UITable.Visible = 'off'; 
            app.SelecttheregionofinterestROIoptionalLabel_2.Visible = 'off'; 
            app.Select_Orientation_UITable.Visible = 'off'; 
            app.Select_Rescale_UITable.Visible = 'off'; 
            app.loaded_array = []; app.current_array = []; app.beforescaling_array = [];
            app.SelecttheregionofinterestROIoptionalLabel_3.Visible = 'off';
            app.Rows12unswapedButton.Visible = 'off'; app.Rows12unswapedButton.Enable = 'off'; app.Rows12unswapedButton.Value = 0; app.Rows12unswapedButton.Text = 'Rows 1-2 (un-swaped)';  
            app.Rows13unswapedButton.Visible = 'off'; app.Rows13unswapedButton.Enable = 'off'; app.Rows13unswapedButton.Value = 0; app.Rows13unswapedButton.Text = 'Rows 1-3 (un-swaped)'; 
            app.Rows23unswapedButton.Visible = 'off'; app.Rows23unswapedButton.Enable = 'off'; app.Rows23unswapedButton.Value = 0; app.Rows23unswapedButton.Text = 'Rows 2-3 (un-swaped)'; 
            app.Row1unflippedButton.Visible = 'off'; app.Row1unflippedButton.Enable = 'off'; app.Row1unflippedButton.Value = 0; app.Row1unflippedButton.Text = 'Rows 1 (un-flipped)'; 
            app.Row2unflippedButton.Visible = 'off'; app.Row2unflippedButton.Enable = 'off'; app.Row2unflippedButton.Value = 0; app.Row2unflippedButton.Text = 'Rows 2 (un-flipped)'; 
            app.Row3unflippedButton.Visible = 'off'; app.Row3unflippedButton.Enable = 'off'; app.Row3unflippedButton.Value = 0; app.Row3unflippedButton.Text = 'Rows 3 (un-flipped)'; 
            app.ScalingfactorEditField.Visible = 'off'; app.ScalingfactorEditField.Editable = 'off'; app.ScalingfactorEditFieldLabel.Visible = 'off';
            app.VoxelsizebeforeEditField.Visible = 'off'; app.VoxelsizebeforeEditField.Editable = 'off'; app.VoxelsizebeforeEditFieldLabel.Visible = 'off';
            app.VoxelsizeafterEditField.Visible = 'off'; app.VoxelsizeafterEditFieldLabel.Visible = 'off';
            app.ApplyscalingButton.Visible = 'off'; app.ApplyscalingButton.Enable = 'off';
            app.BackgroundidEditField.Visible = 'off'; app.BackgroundidEditField.Editable = 'off'; app.BackgroundidEditFieldLabel.Visible = 'off';
        end
        
        function [] = Microstructure_information(app,statut)
            app.VolumefractionsLabel.Visible = statut;
            app.Microstructure_vf_UITable.Visible = statut;
            app.Select_VisualizemicrostructureButton.Enable = statut;            
            app.SelecttheregionofinterestROILabel.Visible = statut; 
            app.Select_ROI_UITable.Visible = statut; 
            app.SelecttheregionofinterestROIoptionalLabel_2.Visible = statut; 
            app.Select_Orientation_UITable.Visible = statut;   
            app.Select_Rescale_UITable.Visible = statut;  
            app.SelecttheregionofinterestROIoptionalLabel_3.Visible = statut;
            app.Rows12unswapedButton.Visible = statut; app.Rows12unswapedButton.Enable = statut; app.Rows12unswapedButton.Value = 0; app.Rows12unswapedButton.Text = 'Rows 1-2 (un-swaped)';  
            app.Rows13unswapedButton.Visible = statut; app.Rows13unswapedButton.Enable = statut; app.Rows13unswapedButton.Value = 0; app.Rows13unswapedButton.Text = 'Rows 1-3 (un-swaped)'; 
            app.Rows23unswapedButton.Visible = statut; app.Rows23unswapedButton.Enable = statut; app.Rows23unswapedButton.Value = 0; app.Rows23unswapedButton.Text = 'Rows 2-3 (un-swaped)'; 
            app.Row1unflippedButton.Visible = statut; app.Row1unflippedButton.Enable = statut; app.Row1unflippedButton.Value = 0; app.Row1unflippedButton.Text = 'Rows 1 (un-flipped)'; 
            app.Row2unflippedButton.Visible = statut; app.Row2unflippedButton.Enable = statut; app.Row2unflippedButton.Value = 0; app.Row2unflippedButton.Text = 'Rows 2 (un-flipped)'; 
            app.Row3unflippedButton.Visible = statut; app.Row3unflippedButton.Enable = statut; app.Row3unflippedButton.Value = 0; app.Row3unflippedButton.Text = 'Rows 3 (un-flipped)'; 
            app.ScalingfactorEditField.Visible = statut; app.ScalingfactorEditField.Editable = statut; app.ScalingfactorEditFieldLabel.Visible = statut;
            app.VoxelsizebeforeEditField.Visible = statut; app.VoxelsizebeforeEditField.Editable = statut; app.VoxelsizebeforeEditFieldLabel.Visible = statut;
            app.VoxelsizeafterEditField.Visible = statut; app.VoxelsizeafterEditFieldLabel.Visible = statut;
            app.ApplyscalingButton.Visible = statut; app.ApplyscalingButton.Enable = statut;      
            app.BackgroundidEditField.Visible = statut; app.BackgroundidEditField.Editable = statut; app.BackgroundidEditFieldLabel.Visible = statut;            
            if strcmp(statut,'on')
                % Volume fraction
                phases_id = unique(app.loaded_array);
                n_phases = length(phases_id);
                vf = zeros(n_phases,4);
                for k=1:1:n_phases
                    vf(k,1) = phases_id(k);
                    vf(k,2) = sum(sum(sum( app.loaded_array==vf(k,1)  )))/numel(app.loaded_array);
                    vf(k,3) = vf(k,2); vf(k,4) = vf(k,2);
                end
                app.Microstructure_vf_UITable.Data = vf;
                % Roi
                domain_size = size(app.loaded_array);
                roi = zeros(3,5);
                roi(:,1) = [1;2;3];
                roi(:,2) = domain_size;
                roi(:,3) = [1;1;1];
                roi(:,4) = domain_size;
                roi(:,5) = domain_size;
                app.Select_ROI_UITable.Data = roi;
                app.roi_sav = roi;
                % Orientation
                app.Select_Orientation_UITable.Data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell([1;2;3]) num2cell(domain_size')];
                % Scaling
                app.Select_Rescale_UITable.Data =  [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(domain_size') num2cell(domain_size')];              
            end
        end
        

        
        function [] = Check_inplane_dimension(app)
            voxel_IP_1 = cell2mat(app.Dim_IP1_UITable.Data(:,4)) - cell2mat(app.Dim_IP1_UITable.Data(:,3)) +1;
            voxel_IP_2 = cell2mat(app.Dim_IP2_UITable.Data(:,4)) - cell2mat(app.Dim_IP2_UITable.Data(:,3)) +1;
            % Actualize table
            app.Dim_IP1_UITable.Data(:,5) = num2cell(voxel_IP_1.*app.VoxelsizeEditField.Value);
            app.Dim_IP2_UITable.Data(:,5) = num2cell(voxel_IP_2.*app.VoxelsizeEditField.Value);
            % Check compatibility
            compatibility_IP = [0 0];
            if length(unique(voxel_IP_1))==1
                app.Inplanedimension1arematchingLabel.Text = 'In-plane dimension 1 are matching';
                app.Inplanedimension1arematchingLabel.FontColor = 'k';
                compatibility_IP(1) = 1;
            else
                app.Inplanedimension1arematchingLabel.Text = 'In-plane dimension 1 are NOT matching';
                app.Inplanedimension1arematchingLabel.FontColor = 'r';
            end
            if length(unique(voxel_IP_2))==1
                app.Inplanedimension2arematchingLabel.Text = 'In-plane dimension 2 are matching';
                app.Inplanedimension2arematchingLabel.FontColor = 'k';
                compatibility_IP(2) = 1;
            else
                app.Inplanedimension2arematchingLabel.Text = 'In-plane dimension 2 are NOT matching';
                app.Inplanedimension2arematchingLabel.FontColor = 'r';
            end        
            if sum(compatibility_IP)==2
                app.SavecellButton.Enable = 'on'; app.Select_VisualizecellButton.Enable = 'on';
            else
                app.SavecellButton.Enable = 'off'; app.Select_VisualizecellButton.Enable = 'off';
            end
        end
        
        function [] = CropAssembleCell(app)
            n_domain = 0;
            cell_sz = zeros(1,3);
            
            % Through-plane dimension
            all_dimension_TP = cell2mat(app.Dim_TP_UITable.Data(:,2));
            cell_sz(1) = sum(all_dimension_TP);
            % In-plane dimension 1
            all_dimension_IP1 = cell2mat(app.Dim_IP1_UITable.Data(:,3:4));
            all_length_IP1 = all_dimension_IP1(:,2) - all_dimension_IP1(:,1) + 1;
            min_length = min(all_length_IP1);
            cell_sz(2) = all_length_IP1(1);
            % In-plane dimension 2
            all_dimension_IP2 = cell2mat(app.Dim_IP2_UITable.Data(:,3:4));
            all_length_IP2 = all_dimension_IP2(:,2) - all_dimension_IP2(:,1) + 1;
            min_length = min(all_length_IP2);           
            cell_sz(3) = all_length_IP2(1);
            
            % Initialize
            app.cellarray_firstassemble = zeros(cell_sz,'uint16');
            x0 = 1;
            for k=1:1:5
                if ~isempty(app.Domain(k).array)
                    n_domain = n_domain+1;
                    x1 = x0 + all_dimension_TP(n_domain) - 1;
                    app.domainbound(n_domain).x0 = x0;
                    app.domainbound(n_domain).x1 = x1;
                    if numel(app.Domain(k).array)>1 % Heterogeneous
                        sz = size(app.Domain(k).array);
                        y0 = all_dimension_IP1(n_domain,1);
                        y1 = all_dimension_IP1(n_domain,2);
                        z0 = all_dimension_IP2(n_domain,1);
                        z1 = all_dimension_IP2(n_domain,2);
                        app.Domain_in_cell(n_domain).array_firstassemble = app.Domain(k).array(:,y0:y1,z0:z1);
                        app.cellarray_firstassemble(x0:x1,:,:) = app.Domain_in_cell(n_domain).array_firstassemble;
                    else % Homogeneous
                        app.cellarray_firstassemble(x0:x1,:,:) = 1;
                        app.Domain_in_cell(n_domain).array_firstassemble = app.cellarray_firstassemble(x0:x1,:,:);
                    end
                    x0 = x1+1;
                    app.Domain_in_cell(n_domain).name = char(app.Dim_TP_UITable.Data(n_domain,1));
                end
            end
            
        end
        
        
        function [] = update_celltable_id(app)
            c1={}; c2={}; c3={}; c4={}; c5={};
            line=0;
            n_domain=length(app.Domain_in_cell);
            for kdomain=1:1:n_domain
                tmp = app.cell_before_reassignemnt(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:);
                ids = unique(tmp);
                for kid = 1:1:length(ids)
                    line=line+1;
                    c1(line,1) = {app.Domain_in_cell(kdomain).name};
                    c2(line,1) = {ids(kid)};
                    c3(line,1) = {line};
                    app.cell_after_reassignement( tmp==ids(kid) ) = line;
                    if app.micro_poly
                        if ids(kid)==0
                            c4(line,1) = {[app.Domain_in_cell(kdomain).name '_pore']};
                            c5(line,1) =  {0};
                        else
                            c4(line,1) = {[app.Domain_in_cell(kdomain).name '_grain_' num2str(ids(kid))]};
                            c5(line,1) =  {kdomain};
                        end
                        app.Assemble_UITable.Enable = 'off';
                        app.UpdateButton.Enable = 'off';
                    else
                        app.Assemble_UITable.Enable = 'on';
                        app.UpdateButton.Enable = 'on';
                        if length(ids)==1
                            c4(line,1) = {[app.Domain_in_cell(kdomain).name '_homogenized']};
                            c5(line,1) =  {0};
                        else
                            if ids(kid)==0
                                c4(line,1) = {[app.Domain_in_cell(kdomain).name '_pore']};
                                c5(line,1) =  {0};
                            elseif ids(kid)==1
                                c4(line,1) = {[app.Domain_in_cell(kdomain).name '_activematerial']};
                                c5(line,1) =  {kdomain};
                            elseif ids(kid)==2
                                c4(line,1) = {[app.Domain_in_cell(kdomain).name '_additive']};
                                c5(line,1) =  {kdomain};
                            else
                                c4(line,1) = {[app.Domain_in_cell(kdomain).name '_phase_' num2str(ids(kid))]};
                                c5(line,1) =  {kdomain};
                            end
                        end
                    end
                    
                end
            end
            allgroupid = cell2mat(c5);
            u = unique(allgroupid);
            counter=0;
            for kg = 1:1:length(u)
                if u(kg)>0
                    counter=counter+1;
                    idx = find(allgroupid==u(kg));
                    c5(idx) ={counter};
                end
            end
            % Update table
            app.Assemble_UITable.Data = [c1 c2 c3 c4 c5];  
            
            % Group name table
            allgroupid = cell2mat(c5);
            u = unique(allgroupid);            
            c1={}; c2={};
            for kg = 1:1:length(u)
                c1(kg,1) = {u(kg)};
                c2(kg,1) ={['Group_' num2str(u(kg))]};
            end
            app.Group_UITable.Data = [c1 c2]; 
            
            % Update array
            app.update_cellid
    
            % Update GUI
            app.Assemble_VisualizecellButton.Enable = 'on';
            app.Assemble_SavecellButton.Enable = 'on';
        end
        
        function [] = update_cellid(app)
            app.cell_after_reassignement = app.cell_before_reassignemnt; % Initialize
            app.cell_group_after_reassignement = app.cell_before_reassignemnt; % Initialize
            line=0;
            n_domain=length(app.Domain_in_cell);
            old_ids = cell2mat(app.Assemble_UITable.Data(:,2));
            new_ids = cell2mat(app.Assemble_UITable.Data(:,3));
            group_ids = cell2mat(app.Assemble_UITable.Data(:,5));
            for kdomain=1:1:n_domain
                tmp = app.cell_before_reassignemnt(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:);
                tmp2 = tmp; tmp3 = tmp;
                ids = unique(tmp);
                for kid = 1:1:length(ids)
                    line=line+1;
                    idx = tmp == old_ids(line);
                    tmp2( idx ) = new_ids(line);
                    tmp3( idx ) = group_ids(line);
                end
                app.cell_after_reassignement(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp2;
                app.cell_group_after_reassignement(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp3;
            end
            app.Assemble_VisualizecellButton.Enable = 'on';
            app.Assemble_SavecellButton.Enable = 'on';            
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
            app.reset_microstructure
            app.LengthinvoxelEditField.Visible = 'off'; app.LengthinvoxelEditField.Enable = 'off';
            app.LengthinvoxelEditFieldLabel.Visible = 'off';
            app.Setup_UITable.Data = table({'Parameters';'Loaded tif';'Modified Tif during pre-processing tasks' },{app.save_parameters;app.save_loaded_tif;app.save_modified_tif});
            app.Meshformat_UITable.Data = table({'MATLAB .mat';'Datasheet .csv';'GMSH .msh';'ABAQUS .inp';'STL .stl';'Binary STL .stl'},{app.save_mat;app.save_csv;app.save_msh;app.save_inp;app.save_stl;app.save_binarystl});
            app.Meshchoice_UITable.Data = table({'One mesh for each phase';'One mesh for each group of phase';'One mesh for the whole geometry'},{app.do_one_mesh_per_phase;app.do_one_mesh_per_group;app.do_one_mesh_wholegeometry});
            app.TiffileloadedLabel.Visible = 'off';     
            
            app.morphology_operation(1).name = 'Erosion + dilatation';
            app.morphology_operation(2).name = 'Clean voxel connection';
            app.morphology_operation(3).name = 'Convert to unique connected cluster';
            
            app.VolumefractionsLabel.Visible = 'off';
            app.Microstructure_vf_UITable.Visible = 'off';
            app.Select_VisualizemicrostructureButton.Enable = 'off'; 
            
            % Iso2mesh table
            str_radbound = 'Max radius of the Delaunay sphere (control the facet size at the surface)';
            str_distbound = 'Maximum deviation from the specified isosurfaces';
            data_ = [{'radbound';'distbound'} {str_radbound;str_distbound} num2cell([1;1])];
            app.CGAL_UITable.Data=data_;
            
            app.Surfacemesh_MethodDropDownValueChanged
            
            str_keepratio = 'Percentage of elements being kept after the simplification';
            str_maxvol = 'Maximum tetrahedra element volume';
            data_ = [{'keepratio';'maxvol'} {str_keepratio;str_maxvol} num2cell([1;10])];
            app.Meshdensity_UITable.Data=data_;
                        
            iso2mesh_is_installed = exist('iso2meshver');
            if iso2mesh_is_installed
                app.Iso2meshversionnumberLabel.Text = ['Iso2mesh version number: ' iso2meshver];
            else
                warning('Iso2mesh is not installed: you can do pre-processing tasks and create regular (structured) mesh, but you cannot create unstrucured mesh. Please look at documentation to know how to install Iso2mesh before continuing');
            end
            app.EachvoxelwillberepresentedwithDropDownValueChanged;
        end

        % Button pushed function: OpendocumentationButton
        function OpendocumentationButtonPushed(app, event)
            Find_file('NREL_MATBOX_Microstructure_analysis_toolbox_documentation.pdf','MATBOX_Microstructure_analysis_toolbox','Default location is \MATBOX_Microstructure_analysis_toolbox\Documentation\');            
        end

        % Image clicked function: About_Logo_NREL
        function About_Logo_NRELImageClicked(app, event)
            url = 'https://www.nrel.gov/transportation/';
            web(url)
        end

        % Button pushed function: Iso2MeshwebsiteButton
        function Iso2MeshwebsiteButtonPushed(app, event)
            url = 'http://iso2mesh.sourceforge.net/cgi-bin/index.cgi';
            web(url)
        end

        % Button pushed function: ClicktoselectsavefolderButton
        function ClicktoselectsavefolderButtonPushed(app, event)
            str_dialogbox = 'Select location where the save folder will be created';
            resultfolder_location = uigetdir(matlabroot,str_dialogbox); % Open dialog box
            if resultfolder_location==0
                % User clicked cancel button or closed the dialog box
                set(app.Save_folder_text,'Text','Save folder location: NOT DEFINED','FontColor','r');
                app.TabisNOTsetupcorrectlyLabel.Text = 'Tab is NOT setup correctly';
                app.Folder_Lamp.Color='r';
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
                app.Folder_Lamp.Color='g';
                app.TabisNOTsetupcorrectlyLabel.Text = 'Tab is setup correctly';
            end
        end

        % Button pushed function: ImporttifButton_Unique
        function ImporttifButton_UniquePushed(app, event)
            app.datatype_array = 'uint8';
            [array, outcome] = app.function_import_microstructure;
            if outcome
                app.loaded_array = array; app.current_array = array; app.beforescaling_array = array;
                app.micro_unique = true; app.micro_cell = false; app.micro_poly = false;
                app.SavemicrostructureButton.Enable = 'on';
                app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.9 0.9 0.9];
                app.CC_left.BackgroundColor = [0.9 0.9 0.9];
                app.LeftelectrodeLabel.BackgroundColor = [0.9 0.9 0.9];
                app.SeparatorLabel.BackgroundColor = [0.9 0.9 0.9];
                app.RightelectrodeLabel.BackgroundColor = [0.9 0.9 0.9];
                app.CC_right.BackgroundColor = [0.9 0.9 0.9];
                app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.9 0.9 0.9];
                app.Microstructure_information('on');
            else
                app.loaded_array = []; app.current_array = []; app.beforescaling_array = [];         
                app.Microstructure_information('off');
            end
        end

        % Button pushed function: ImporttifButton_Cell
        function ImporttifButton_CellPushed(app, event)
            app.datatype_array = 'uint8';
            [array, outcome] = app.function_import_microstructure;
            if outcome
                app.loaded_array = array;  app.current_array = array; app.beforescaling_array = array;
                app.micro_unique = false; app.micro_cell = true; app.micro_poly = false;
                app.SavemicrostructureButton.Enable = 'on';
                app.OneuniquemicrostructureLabel.BackgroundColor = [0.9 0.9 0.9];
                app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.9 0.9 0.9];
                app.Microstructure_information('on');
            else
                app.loaded_array = []; app.current_array = []; app.beforescaling_array = [];           
                app.Microstructure_information('off');
            end
            app.Choosedomain = app.ChoosedomainDropDown.Value;
            app.ChoosedomainDropDown.Value = '';
            app.ImporttifButton_Cell.Enable = 'off'; app.HomogenousButton_Cell.Enable = 'off';
        end

        % Button pushed function: ImporttifButton_Poly
        function ImporttifButton_PolyPushed(app, event)
            app.datatype_array = 'uint16';
            [array, outcome] = app.function_import_microstructure;
            if outcome
                background_code = 0; % By default
                % Add one layer of background (so that background connectivity will be correctly calculated later)
                tmp = zeros(size(array)+2)+background_code;
                tmp(2:end-1,2:end-1,2:end-1) = array;
                array = tmp; clear tmp;
                % Re-order grain, if required
                ids = unique(array);
                ids(ids==background_code)=[];
                tmp = zeros(size(array));
                for k=1:1:length(ids)
                    tmp(array==ids(k))=k;
                end
                array = tmp; clear tmp;
                
                app.loaded_array = array;  app.current_array = array; app.beforescaling_array = array;  
                app.micro_unique = false; app.micro_cell = false; app.micro_poly = true;
                app.SavemicrostructureButton.Enable = 'on';
                app.OneuniquemicrostructureLabel.BackgroundColor = [0.9 0.9 0.9];
                app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.9 0.9 0.9];
                app.CC_left.BackgroundColor = [0.9 0.9 0.9];
                app.LeftelectrodeLabel.BackgroundColor = [0.9 0.9 0.9];
                app.SeparatorLabel.BackgroundColor = [0.9 0.9 0.9];
                app.RightelectrodeLabel.BackgroundColor = [0.9 0.9 0.9];    
                app.CC_right.BackgroundColor = [0.9 0.9 0.9];
                app.Microstructure_information('on');
            else
                app.loaded_array = [];  app.current_array = []; app.beforescaling_array = [];
                app.Microstructure_information('off');
            end
        end

        % Value changed function: ChoosedomainDropDown
        function ChoosedomainDropDownValueChanged(app, event)
            app.ImporttifButton_Cell.Enable = 'on';
            app.HomogenousButton_Cell.Enable = 'on';
        end

        % Button pushed function: HomogenousButton_Cell
        function HomogenousButton_CellPushed(app, event)
            app.LengthinvoxelEditFieldLabel.Visible = 'on';
            app.LengthinvoxelEditField.Value = 5; app.LengthinvoxelEditField.Visible = 'on'; app.LengthinvoxelEditField.Enable = 'on';
            app.loaded_array = app.LengthinvoxelEditField.Value; app.current_array = app.LengthinvoxelEditField.Value; app.beforescaling_array = app.LengthinvoxelEditField.Value;
            app.micro_unique = false; app.micro_cell = true; app.micro_poly = false;
            app.SavemicrostructureButton.Enable = 'on';
            app.OneuniquemicrostructureLabel.BackgroundColor = [0.9 0.9 0.9];
            app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.9 0.9 0.9];
            app.Choosedomain = app.ChoosedomainDropDown.Value;
            app.ChoosedomainDropDown.Value = '';
            app.ImporttifButton_Cell.Enable = 'off'; app.HomogenousButton_Cell.Enable = 'off';
        end

        % Button pushed function: SavemicrostructureButton
        function SavemicrostructureButtonPushed(app, event)
            app.SavemicrostructureButton.Text = 'Saving... please wait'; pause(0.1);
            if app.micro_unique
                app.Domain(1).array = app.current_array;
                sz = size(app.Domain(1).array);
                app.domainbound(1).x0 = 1; app.domainbound(1).x1 = sz(1); 
                app.OneuniquemicrostructureLabel.BackgroundColor = [0.93 0.69 0.13];
                app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.9 0.9 0.9];
                app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.9 0.9 0.9];
                app.SavemicrostructureButton.Enable = 'off';
                app.Select_Lamp.Color = 'g'; app.ActualizeresetdimensionsButton.Enable = 'on';
                app.TabisNOTsetupcorrectlyLabel_2.Text = 'Tab is setup correctly';
                app.Dim_YoucanskipthistabLabel.Visible = 'on';
                app.Comp_Lamp.Color = 'g';
                app.TabisNOTsetupcorrectlyLabel_3.Text = 'Tab is setup correctly';
                app.Domain_in_cell(1).array_firstassemble = app.current_array;
                app.Domain_in_cell(1).name = 'Unique microstructure';
                app.cellarray_firstassemble = app.current_array;                
                app.cellarray_secondassemble = app.cellarray_firstassemble;
                app.MO_Microstructure = app.cellarray_secondassemble; app.MO_Microstructure_saved = app.MO_Microstructure;              
                domain_name = 'Unique microstructure';
                
                app.MO_vf_UITable.Visible = 'on';
                volume_name = {}; phases_id = {}; vf_before = {};
                k=1;
                volume_name(k,1) = {app.Domain_in_cell(k).name};
                tmp = app.Domain_in_cell(k).array_firstassemble;
                phases = unique(tmp); n=length(phases);
                if n>1
                    for p=1:1:n
                        if p==1
                            tmp2 = [num2str(phases(p),'%i') ' / '];
                        elseif p==n
                            tmp2 = [tmp2 num2str(phases(p),'%i')];
                        else
                            tmp2 = [tmp2 num2str(phases(p),'%i') ' / '];
                        end
                        
                        vf=sum(sum(sum( tmp==phases(p)  )))/numel(tmp);
                        if p==1
                            tmp3 = [num2str(vf,'%1.3f') ' / '];
                        elseif p==n
                            tmp3 = [tmp3 num2str(vf,'%1.3f')];
                        else
                            tmp3 = [tmp3 num2str(vf,'%1.3f') ' / '];
                        end
                    end
                else
                    tmp2 = [num2str(phases(1),'%i')];
                    vf=sum(sum(sum( tmp==phases(1)  )))/numel(tmp);
                    tmp3 = [num2str(vf,'%1.3f')];
                end
                phases_id(k,1)={tmp2};
                vf_before(k,1)={tmp3};
                vf_after=vf_before;
                app.MO_vf_UITable.Data = [char(volume_name) phases_id vf_before vf_after];
                
                c1 = {}; c2= {}; c3 ={};
                line=0;
                for kdomain=1:1:1
                    for kstep=1:1:app.n_mo_step
                        line=line+1;
                        c1(line,1) = {app.Domain_in_cell(kdomain).name};
                        c2(line,1) = {app.morphology_operation(kstep).name};
                        if length(unique(app.Domain_in_cell(kdomain).array_firstassemble))==1
                            c3(line,1) = {false}; % Homogeneous domain, there is no need to do morphology opening of any sort
                        else
                            c3(line,1) = {true};
                        end
                    end
                end
                app.MO_operations_UITable.Data = [c1 c2 c3];
                
                app.GrainvolumefilterinnumberofvoxelEditFieldLabel.Enable = 'off';
                app.GrainvolumefilterinnumberofvoxelEditField.Enable = 'off';
                app.MO_operations_UITable.Visible = 'on';    
                
                app.TabGroup.SelectedTab = app.MorphologyopeningTab; % Skip dimension tab
                app.MO_VisualizecellButton.Enable = 'on';  app.MO_SkipButton.Enable = 'on'; 
                app.MO_DoButton.Enable = 'on'; app.MO_UndoButton.Enable = 'off'; app.MO_SaveButton.Enable = 'off'; app.MO_statut_Label.Text = 'Idle';      
                app.MO_SavecellButton.Enable = 'on';    
                app.MO_UnionsoliduniqueclusterCheckBox.Enable = 'on';
                app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Enable = 'off';

            elseif app.micro_cell
                app.MO_vf_UITable.Visible = 'on';
                if strcmp(app.Choosedomain,'Left current collector')
                    app.Domain(1).array = app.current_array;
                    app.CC_left.BackgroundColor = [0.93 0.69 0.13];
                    app.imported_domain(1) = 1;
                    domain_name = 'Left_current_collector';
                elseif strcmp(app.Choosedomain,'Left electrode')
                    app.Domain(2).array = app.current_array;
                    app.LeftelectrodeLabel.BackgroundColor = [0.93 0.69 0.13];
                    app.imported_domain(2) = 1;
                    domain_name = 'Left_electrode';
                elseif strcmp(app.Choosedomain,'Separator')
                    app.Domain(3).array = app.current_array;
                    app.SeparatorLabel.BackgroundColor = [0.93 0.69 0.13];
                    app.imported_domain(3) = 1;
                    domain_name = 'Separator';
                elseif strcmp(app.Choosedomain,'Right electrode')
                    app.Domain(4).array = app.current_array;
                    app.RightelectrodeLabel.BackgroundColor = [0.93 0.69 0.13];
                    app.imported_domain(4) = 1;
                    domain_name = 'Right_electrode';
                elseif strcmp(app.Choosedomain,'Right current collector')
                    app.Domain(5).array = app.current_array;
                    app.CC_right.BackgroundColor = [0.93 0.69 0.13];
                    app.imported_domain(5) = 1;
                    domain_name = 'Right_current_collector';                    
                end
                if sum(app.imported_domain)>=2 && max(bwlabel(app.imported_domain))==1 % At least two domains, all domains in contact
                    app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.93 0.69 0.13];
                    app.Select_Lamp.Color = 'g'; app.ActualizeresetdimensionsButton.Enable = 'on';
                    app.TabisNOTsetupcorrectlyLabel_2.Text = 'Tab is setup correctly';
                end
                app.Dim_YoucanskipthistabLabel.Visible = 'off';
                app.Comp_Lamp.Color = 'r';
                app.TabisNOTsetupcorrectlyLabel_3.Text = 'Tab is NOT setup correctly';
                app.MO_UnionsoliduniqueclusterCheckBox.Enable = 'on';
                app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Enable = 'off';
                
            elseif app.micro_poly
                app.Domain(1).array = app.current_array;
                sz = size(app.Domain(1).array);
                app.domainbound(1).x0 = 1; app.domainbound(1).x1 = sz(1);                 
                app.OneuniquemicrostructureLabel.BackgroundColor = [0.9 0.9 0.9];
                app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.9 0.9 0.9];
                app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.93 0.69 0.13];
                app.SavemicrostructureButton.Enable = 'off';
                app.Select_Lamp.Color = 'g'; app.ActualizeresetdimensionsButton.Enable = 'on';
                app.TabisNOTsetupcorrectlyLabel_2.Text = 'Tab is setup correctly';
                app.Dim_YoucanskipthistabLabel.Visible = 'on';
                app.Comp_Lamp.Color = 'g';
                app.TabisNOTsetupcorrectlyLabel_3.Text = 'Tab is setup correctly';
                app.Domain_in_cell(1).array_firstassemble = app.current_array;
                app.Domain_in_cell(1).name = 'Polychrystalline';
                app.cellarray_firstassemble = app.current_array;        
                app.cellarray_secondassemble = app.cellarray_firstassemble;
                app.MO_Microstructure = app.cellarray_secondassemble; app.MO_Microstructure_saved = app.MO_Microstructure;
                app.MO_UnionsoliduniqueclusterCheckBox.Value = 0; app.MO_UnionsoliduniqueclusterCheckBox.Enable = 'off';
                app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Enable = 'on';
                domain_name = 'Polychrystalline';    
                
                c1 = {}; c2= {}; c3 ={};
                line=0;
                for kdomain=1:1:1
                    for kstep=1:1:app.n_mo_step
                        line=line+1;
                        c1(line,1) = {app.Domain_in_cell(kdomain).name};
                        c2(line,1) = {app.morphology_operation(kstep).name};
                        if length(unique(app.Domain_in_cell(kdomain).array_firstassemble))==1
                            c3(line,1) = {false}; % Homogeneous domain, there is no need to do morphology opening of any sort
                        else
                            c3(line,1) = {true};
                        end
                    end
                end
                
                app.GrainvolumefilterinnumberofvoxelEditFieldLabel.Enable = 'on';
                app.GrainvolumefilterinnumberofvoxelEditField.Enable = 'on';                
                app.MO_operations_UITable.Data = [c1 c2 c3];                

                tmp = app.Domain_in_cell(1).array_firstassemble;
                grain_id = unique(tmp); n=length(grain_id);
                vf_before = zeros(n,1); vf_after = zeros(n,1);
                volume_name ={}; gid = {}; vf_before={};
                for p=1:1:n
                    gid(p,1) = {grain_id(p)};
                    volume_name(p,1) = {app.Domain_in_cell(1).name};
                    vf_before(p,1) = {sum(sum(sum( tmp==grain_id(p) ))) / numel(tmp)};
                end
                vf_after=vf_before;
                app.MO_vf_UITable.Data = [volume_name gid vf_before vf_after];
                app.TabGroup.SelectedTab = app.MorphologyopeningTab; % Skip dimension tab
                app.MO_VisualizecellButton.Enable = 'on';  app.MO_SkipButton.Enable = 'on'; 
                app.MO_DoButton.Enable = 'on'; app.MO_UndoButton.Enable = 'off'; app.MO_SaveButton.Enable = 'off'; app.MO_statut_Label.Text = 'Idle';      
                app.MO_SavecellButton.Enable = 'on';                   
            end

            % Save setup
            if ~isempty(app.mainsavefoldersave)
                if app.save_parameters
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Selection\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Selection/'];
                    end
                    % Prepare the data
                    clear DATA_writetable
                    if numel(app.current_array)==1
                        DATA_writetable.sheet(1).name='Homogenous_medium';
                        DATA_writetable.sheet(1).table=table(app.current_array,'VariableNames',{'Length'});                        
                    else
                        DATA_writetable.sheet(1).name='Tiff';
                        DATA_writetable.sheet(1).table=table({app.TiffileloadedLabel.Text},'VariableNames',{'Path'});
                        DATA_writetable.sheet(2).name='Volume_fraction';
                        d=app.Microstructure_vf_UITable.Data;
                        DATA_writetable.sheet(2).table=table(d(:,1),d(:,2),d(:,3),d(:,4),'VariableNames',app.Microstructure_vf_UITable.ColumnName');
                        DATA_writetable.sheet(3).name='Region_of_interest';
                        d=app.Select_ROI_UITable.Data;
                        DATA_writetable.sheet(3).table=table(d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),'VariableNames',app.Select_ROI_UITable.ColumnName');
                        DATA_writetable.sheet(4).name='Orientation';
                        d=app.Select_Orientation_UITable.Data;
                        DATA_writetable.sheet(4).table=table(d(:,1),d(:,2),d(:,3),'VariableNames',app.Select_Orientation_UITable.ColumnName');
                        DATA_writetable.sheet(5).name='Scaling';
                        DATA_writetable.sheet(5).table=table(app.ScalingfactorEditField.Value,'VariableNames',{'Scaling factor'});
                    end
                    % Save function
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end
                    Function_Writetable(current_folder,domain_name,DATA_writetable)
                end
                if app.save_loaded_tif && numel(app.current_array)>1
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Selection\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Selection/'];
                    end
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end                    
                    function_save_tif(app.loaded_array,[current_folder domain_name '_loaded.tif']);
                    function_save_tif(app.current_array,[current_folder domain_name '_modified.tif']); 
                end                    
            end
            app.TiffileloadedLabel.Visible = 'off'; 
            app.LengthinvoxelEditFieldLabel.Visible = 'off'; 
            app.LengthinvoxelEditField.Visible = 'off'; 
            app.SavemicrostructureButton.Enable = 'off'; 
            app.Microstructure_information('off');         
            app.SavemicrostructureButton.Text = 'Save microstructure';
        end

        % Button pushed function: 
        % Select_VisualizemicrostructureButton
        function Select_VisualizemicrostructureButtonPushed(app, event)
            direction_name = {'Normal to Through-plane direction';'Normal to In-plane direction 1';'Normal to In-plane direction 2'};
            Microstructure_basic_visualization_interface(app.current_array, direction_name);            
        end

        % Cell edit callback: Select_ROI_UITable
        function Select_ROI_UITableCellEdit(app, event)
            newData = event.NewData;
            modified_table = app.Select_ROI_UITable.Data;
            domain_size = size(app.loaded_array);
            if isnumeric(newData)
                if round(newData)-newData == 0 % Check integer
                    if sum(modified_table(:,3) < modified_table(:,4))==3 && sum(modified_table(:,3) >= ones(3,1))==3 && sum(modified_table(:,4) <= domain_size')==3
                        roi = modified_table(:,3:4);
                        app.current_array = app.loaded_array(roi(1,1):roi(1,2), roi(2,1):roi(2,2), roi(3,1):roi(3,2));
                        app.beforescaling_array = app.current_array;
                        app.Select_ROI_UITable.Data(:,5) = roi(:,2)-roi(:,1)+ones(3,1);
                        app.roi_sav = app.Select_ROI_UITable.Data(:,3:end);
                    else
                        app.Select_ROI_UITable.Data(:,3:end) = app.roi_sav;
                    end
                else
                    app.Select_ROI_UITable.Data(:,3:end) = app.roi_sav;
                end
            else
                app.Select_ROI_UITable.Data(:,3:end) = app.roi_sav;
            end
            % Update volume fraction
            phases_id = app.Microstructure_vf_UITable.Data(:,1);
            n_phases = length(phases_id);
            vf = zeros(n_phases,1);
            for k=1:1:n_phases
                vf(k,1) = sum(sum(sum( app.current_array==phases_id(k)  )))/numel(app.current_array);
            end
            app.Microstructure_vf_UITable.Data(:,3) = vf;
            
            % Reset orientation
            domain_size = size(app.current_array);
            app.Select_Orientation_UITable.Data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell([1;2;3]) num2cell(domain_size')];
            app.Rows12unswapedButton.Text = 'Rows 1-2 (un-swaped)'; app.Rows12unswapedButton.Value = 0;
            app.Rows13unswapedButton.Text = 'Rows 1-3 (un-swaped)'; app.Rows13unswapedButton.Value = 0;
            app.Rows23unswapedButton.Text = 'Rows 2-3 (un-swaped)'; app.Rows23unswapedButton.Value = 0;
            app.Row1unflippedButton.Text = 'Row 1 (un-flipped)'; app.Row1unflippedButton.Value = 0;
            app.Row2unflippedButton.Text = 'Row 2 (un-flipped)'; app.Row2unflippedButton.Value = 0;
            app.Row3unflippedButton.Text = 'Row 3 (un-flipped)'; app.Row3unflippedButton.Value = 0;
            % Reset voxel size
            app.Select_Rescale_UITable.Data =  [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(domain_size') num2cell(domain_size')];              

            
            
        end

        % Value changed function: LengthinvoxelEditField
        function LengthinvoxelEditFieldValueChanged(app, event)
            app.loaded_array = app.LengthinvoxelEditField.Value; app.current_array = app.LengthinvoxelEditField.Value; app.beforescaling_array = app.LengthinvoxelEditField.Value;
        end

        % Value changed function: Rows12unswapedButton
        function Rows12unswapedButtonValueChanged(app, event)
            value = app.Rows12unswapedButton.Value;
            if ~value
                app.Rows12unswapedButton.Text = 'Rows 1-2 (un-swaped)';
            else
                app.Rows12unswapedButton.Text = 'Rows 1-2 (swaped)';
            end
            domain_size = size(app.current_array);
            New_domain_size(1)=domain_size(2); New_domain_size(2)=domain_size(1); New_domain_size(3)=domain_size(3);
            tmp=zeros(New_domain_size,app.datatype_array);
            for k=1:1:domain_size(2)
                slice=app.current_array(:,k,:);
                tmp(k,:,:)=slice;
            end
            app.current_array = tmp; app.beforescaling_array = tmp;
            domain_size = size(app.current_array);
            % Update table
            old_axis = app.Select_Orientation_UITable.Data(:,2);
            new_axis = old_axis;
            new_axis(2) = old_axis(1); new_axis(1) = old_axis(2);
            app.Select_Orientation_UITable.Data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} new_axis num2cell(domain_size')];
            % Reset voxel size
            app.Select_Rescale_UITable.Data =  [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(domain_size') num2cell(domain_size')];
        end

        % Value changed function: Rows13unswapedButton
        function Rows13unswapedButtonValueChanged(app, event)
            value = app.Rows13unswapedButton.Value;
            if ~value
                app.Rows13unswapedButton.Text = 'Rows 1-3 (un-swaped)';
            else
                app.Rows13unswapedButton.Text = 'Rows 1-3 (swaped)';
            end
            domain_size = size(app.current_array);
            New_domain_size(1)=domain_size(3); New_domain_size(2)=domain_size(2); New_domain_size(3)=domain_size(1);
            tmp=zeros(New_domain_size,app.datatype_array);
            for k=1:1:domain_size(3)
                slice=app.current_array(:,:,k)';
                tmp(k,:,:)=slice;
            end
            app.current_array = tmp;  app.beforescaling_array = tmp;
            domain_size = size(app.current_array);
            % Update table
            old_axis = app.Select_Orientation_UITable.Data(:,2);
            new_axis = old_axis;
            new_axis(1) = old_axis(3); new_axis(3) = old_axis(1);
            app.Select_Orientation_UITable.Data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} new_axis num2cell(domain_size')];
            % Reset voxel size
            app.Select_Rescale_UITable.Data =  [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(domain_size') num2cell(domain_size')];
        end

        % Value changed function: Rows23unswapedButton
        function Rows23unswapedButtonValueChanged(app, event)
            value = app.Rows23unswapedButton.Value;
            if ~value
                app.Rows23unswapedButton.Text = 'Rows 2-3 (un-swaped)';
            else
                app.Rows23unswapedButton.Text = 'Rows 2-3 (swaped)';
                
            end   
            domain_size = size(app.current_array);
            New_domain_size(1)=domain_size(1); New_domain_size(2)=domain_size(3); New_domain_size(3)=domain_size(2);
            tmp=zeros(New_domain_size,app.datatype_array);
            for k=1:1:domain_size(3)
                slice=app.current_array(:,:,k);
                tmp(:,k,:)=slice;
            end
            app.current_array = tmp;  app.beforescaling_array = tmp;
            domain_size = size(app.current_array);
            % Update table
            old_axis = app.Select_Orientation_UITable.Data(:,2);
            new_axis = old_axis;
            new_axis(2) = old_axis(3); new_axis(3) = old_axis(2);
            app.Select_Orientation_UITable.Data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} new_axis num2cell(domain_size')];
            % Reset voxel size
            app.Select_Rescale_UITable.Data =  [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(domain_size') num2cell(domain_size')];
        end

        % Value changed function: Row1unflippedButton
        function Row1unflippedButtonValueChanged(app, event)
            value = app.Row1unflippedButton.Value;
            if ~value
                app.Row1unflippedButton.Text = 'Row 1 (un-flipped)';
            else
                app.Row1unflippedButton.Text = 'Row 1 (flipped)';
            end
            domain_size = size(app.current_array);
            tmp=zeros(domain_size,app.datatype_array);
            for k=1:1:domain_size(1)
                slice=app.current_array(k,:,:);
                tmp(domain_size(1)-k+1,:,:)=slice;
            end
            app.current_array = tmp; app.beforescaling_array = tmp;
        end

        % Value changed function: Row2unflippedButton
        function Row2unflippedButtonValueChanged(app, event)
            value = app.Row2unflippedButton.Value;
            if ~value
                app.Row2unflippedButton.Text = 'Row 2 (un-flipped)';
            else
                app.Row2unflippedButton.Text = 'Row 2 (flipped)';
            end     
            domain_size = size(app.current_array);
            tmp=zeros(domain_size,app.datatype_array);
            for k=1:1:domain_size(2)
                slice=app.current_array(:,k,:);
                tmp(:,domain_size(2)-k+1,:)=slice;
            end
            app.current_array = tmp; app.beforescaling_array = tmp;
        end

        % Value changed function: Row3unflippedButton
        function Row3unflippedButtonValueChanged(app, event)
            value = app.Row3unflippedButton.Value;
            if ~value
                app.Row3unflippedButton.Text = 'Row 3 (un-flipped)';
            else
                app.Row3unflippedButton.Text = 'Row 3 (flipped)';
            end  
            domain_size = size(app.current_array);
            tmp=zeros(domain_size,app.datatype_array);
            for k=1:1:domain_size(3)
                slice=app.current_array(:,:,k);
                tmp(:,:,domain_size(3)-k+1)=slice;
            end
            app.current_array = tmp; app.beforescaling_array = tmp;
        end

        % Value changed function: ScalingfactorEditField
        function ScalingfactorEditFieldValueChanged(app, event)
            app.VoxelsizeafterEditField.Value = app.VoxelsizebeforeEditField.Value * app.ScalingfactorEditField.Value;
            app.Select_Rescale_UITable.Data(:,3) = num2cell (round( cell2mat(app.Select_Rescale_UITable.Data(:,2)) / app.ScalingfactorEditField.Value));
        end

        % Button pushed function: ApplyscalingButton
        function ApplyscalingButtonPushed(app, event)
            % Set parameters
            parameters_scaling.scaling_factor = app.ScalingfactorEditField.Value;
            parameters_scaling.label_or_greylevel = 'Label';
            parameters_scaling.background = app.BackgroundidEditField.Value;
            % Scale
            app.current_array = function_scaling(app.beforescaling_array,parameters_scaling);  
            % Update volume fraction
            phases_id = app.Microstructure_vf_UITable.Data(:,1);
            n_phases = length(phases_id);
            vf = zeros(n_phases,1);
            for k=1:1:n_phases
                vf(k,1) = sum(sum(sum( app.current_array==phases_id(k)  )))/numel(app.current_array);
            end
            app.Microstructure_vf_UITable.Data(:,4) = vf; 
        end

        % Value changed function: VoxelsizebeforeEditField
        function VoxelsizebeforeEditFieldValueChanged(app, event)
            app.VoxelsizeafterEditField.Value = app.VoxelsizebeforeEditField.Value * app.ScalingfactorEditField.Value;
        end

        % Cell edit callback: Setup_UITable
        function Setup_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if indices(1)==1
                app.save_parameters = newData;
            elseif indices(1)==2
                app.save_loaded_tif = newData;
            elseif indices(1)==3
                app.save_modified_tif = newData;             
            end
        end

        % Button pushed function: ActualizeresetdimensionsButton
        function ActualizeresetdimensionsButtonPushed(app, event)
            volume_name = {}; voxel_TP = []; voxel_IP_1 = []; voxel_IP_2 = []; 
            n_domain = 0;
            for k=1:1:5
                if ~isempty(app.Domain(k).array)
                    volume_name = [volume_name; app.ChoosedomainDropDown.Items(k+1)];
                    n_domain = n_domain+1;
                    if numel(app.Domain(k).array)>1 % Heterogeneous
                        sz = size(app.Domain(k).array);
                        voxel_TP = [voxel_TP; sz(1)];
                        voxel_IP_1 = [voxel_IP_1; sz(2)];
                        voxel_IP_2 = [voxel_IP_2; sz(3)];
                    else % Homogeneous
                        voxel_TP = [voxel_TP; app.Domain(k).array];
                        voxel_IP_1 = [voxel_IP_1; NaN];
                        voxel_IP_2 = [voxel_IP_2; NaN];                       
                    end
                end
            end
            % Replace NaN with minimum among the different
            voxel_IP_1(isnan(voxel_IP_1))=min(voxel_IP_1);
            voxel_IP_2(isnan(voxel_IP_2))=min(voxel_IP_2);
            % Actualize table
            length_IP_1 = voxel_IP_1.*app.VoxelsizeEditField.Value;
            app.Dim_IP1_UITable.Data = [volume_name num2cell(voxel_IP_1) num2cell(ones(n_domain,1)) num2cell(voxel_IP_1) num2cell(length_IP_1)];
            length_IP_2 = voxel_IP_2.*app.VoxelsizeEditField.Value;
            app.Dim_IP2_UITable.Data = [volume_name num2cell(voxel_IP_2) num2cell(ones(n_domain,1)) num2cell(voxel_IP_2) num2cell(length_IP_2)];            
            length_TP = voxel_TP.*app.VoxelsizeEditField.Value;
            app.Dim_TP_UITable.Data = [volume_name num2cell(voxel_TP) num2cell(length_TP)];
            % Check dimension
            app.Check_inplane_dimension;
        end

        % Button pushed function: Dim_Autocrop_1_Button
        function Dim_Autocrop_1_ButtonPushed(app, event)
            all_dimension = cell2mat(app.Dim_IP1_UITable.Data(:,3:4));
            all_length = all_dimension(:,2) - all_dimension(:,1) + 1;
            min_length = min(all_length);
            [n,~]=size(all_dimension);
            for k=1:1:n
                if all_length(k)>min_length
                    all_dimension(k,1) = round(all_length(k)/2 - min_length/2);
                    all_dimension(k,2) = round(all_length(k)/2 + min_length/2)-1;
                end
            end
            app.Dim_IP1_UITable.Data(:,3) = num2cell(all_dimension(:,1));
            app.Dim_IP1_UITable.Data(:,4) = num2cell(all_dimension(:,2));
            % Check dimension
            app.Check_inplane_dimension;
        end

        % Button pushed function: Dim_Autocrop_2_Button
        function Dim_Autocrop_2_ButtonPushed(app, event)
            all_dimension = cell2mat(app.Dim_IP2_UITable.Data(:,3:4));
            all_length = all_dimension(:,2) - all_dimension(:,1) + 1;
            min_length = min(all_length);
            [n,~]=size(all_dimension);
            for k=1:1:n
                if all_length(k)>min_length
                    all_dimension(k,1) = round(all_length(k)/2 - min_length/2);
                    all_dimension(k,2) = round(all_length(k)/2 + min_length/2)-1;
                end
            end
            app.Dim_IP2_UITable.Data(:,3) = num2cell(all_dimension(:,1));
            app.Dim_IP2_UITable.Data(:,4) = num2cell(all_dimension(:,2));
            % Check dimension
            app.Check_inplane_dimension;
        end

        % Value changed function: VoxelsizeEditField
        function VoxelsizeEditFieldValueChanged(app, event)
            value = app.VoxelsizeEditField.Value;
            
            % Actualize table
            tmp = app.Dim_IP1_UITable.Data;
            volume_name = tmp(:,1);
            voxel_IP_1 = cell2mat(tmp(:,2));
            n_domain = length(volume_name);
            length_IP_1 = voxel_IP_1.*value;
            app.Dim_IP1_UITable.Data = [volume_name num2cell(voxel_IP_1) num2cell(ones(n_domain,1)) num2cell(voxel_IP_1) num2cell(length_IP_1)];
            
            tmp = app.Dim_IP2_UITable.Data;
            voxel_IP_2 = cell2mat(tmp(:,2));
            length_IP_2 = voxel_IP_2.*value;
            app.Dim_IP2_UITable.Data = [volume_name num2cell(voxel_IP_2) num2cell(ones(n_domain,1)) num2cell(voxel_IP_2) num2cell(length_IP_2)];  
            
            tmp = app.Dim_TP_UITable.Data;
            voxel_TP = cell2mat(tmp(:,2));
            length_TP = voxel_TP.*value;
            app.Dim_TP_UITable.Data = [volume_name num2cell(voxel_TP) num2cell(length_TP)];      
        end

        % Button pushed function: Select_VisualizecellButton
        function Select_VisualizecellButtonPushed(app, event)
            app.CropAssembleCell
            direction_name = {'Normal to Through-plane direction';'Normal to In-plane direction 1';'Normal to In-plane direction 2'};
            Microstructure_basic_visualization_interface(app.cellarray_firstassemble, direction_name); 
        end

        % Button pushed function: SavecellButton
        function SavecellButtonPushed(app, event)
            app.CropAssembleCell
            app.SavecellButton.Text = 'Saving... please wait'; pause(0.1);
            n_domain=length(app.Domain_in_cell);

            % Save setup
            if ~isempty(app.mainsavefoldersave)
                if app.save_parameters
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Dimension_compatibility\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Dimension_compatibility/'];
                    end
                    % Prepare the data
                    clear DATA_writetable
                    d = app.Dim_TP_UITable.Data;
                    DATA_writetable.sheet(1).name='Along_thickness';
                    DATA_writetable.sheet(1).table=table(d(:,1),d(:,2),d(:,3),'VariableNames',app.Dim_TP_UITable.ColumnName');
                    d=app.Dim_IP1_UITable.Data;
                    DATA_writetable.sheet(2).name='In_plane_direction_1';
                    DATA_writetable.sheet(2).table=table(d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),'VariableNames',app.Dim_IP1_UITable.ColumnName');
                    d=app.Dim_IP2_UITable.Data;
                    DATA_writetable.sheet(3).name='In_plane_direction_2';
                    DATA_writetable.sheet(3).table=table(d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),'VariableNames',app.Dim_IP2_UITable.ColumnName');                    
                    % Save function
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end
                    Function_Writetable(current_folder,'Inplane_crop',DATA_writetable)
                end
                if app.save_modified_tif
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Dimension_compatibility\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Dimension_compatibility/'];
                    end
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end                    
                    function_save_tif(app.cellarray_firstassemble,[current_folder 'Cell.tif']);
                    
                    for k=1:1:n_domain
                        [ str_new ] = function_remove_emptyandspecialcharacter_string(app.Domain_in_cell(k).name);
                        function_save_tif(app.Domain_in_cell(k).array_firstassemble,[current_folder str_new '_cropped.tif']);
                    end
                end
            end
            
            volume_name = {}; phases_id = {}; vf_before = {}; vf_after = {};
            for k=1:1:n_domain
                volume_name(k,1) = {app.Domain_in_cell(k).name};
                tmp = app.Domain_in_cell(k).array_firstassemble;
                phases = unique(tmp); n=length(phases);
                if n>1
                    for p=1:1:n
                        if p==1
                            tmp2 = [num2str(phases(p),'%i') ' / '];
                        elseif p==n
                            tmp2 = [tmp2 num2str(phases(p),'%i')];
                        else
                            tmp2 = [tmp2 num2str(phases(p),'%i') ' / '];
                        end
                        
                        vf=sum(sum(sum( tmp==phases(p)  )))/numel(tmp);
                        if p==1
                            tmp3 = [num2str(vf,'%1.3f') ' / '];
                        elseif p==n
                            tmp3 = [tmp3 num2str(vf,'%1.3f')];
                        else
                            tmp3 = [tmp3 num2str(vf,'%1.3f') ' / '];
                        end
                    end
                else
                    tmp2 = [num2str(phases(1),'%i')];
                    vf=sum(sum(sum( tmp==phases(1)  )))/numel(tmp);
                    tmp3 = [num2str(vf,'%1.3f')];
                end
                
                phases_id(k,1)={tmp2};
                vf_before(k,1)={tmp3};
            end
            vf_after=vf_before;
            app.MO_vf_UITable.Data = [volume_name phases_id vf_before vf_after];
            
            c1 = {}; c2= {}; c3 ={};
            line=0;
            for kdomain=1:1:n_domain
                for kstep=1:1:app.n_mo_step
                    line=line+1;
                    c1(line,1) = {app.Domain_in_cell(kdomain).name};
                    c2(line,1) = {app.morphology_operation(kstep).name};
                    if length(unique(app.Domain_in_cell(kdomain).array_firstassemble))==1
                        c3(line,1) = {false}; % Homogeneous domain, there is no need to do morphology opening of any sort
                    else
                        c3(line,1) = {true};
                    end
                end
            end
            
            app.MO_operations_UITable.Data = [c1 c2 c3];
            app.MO_VisualizecellButton.Enable = 'on'; app.MO_SkipButton.Enable = 'on';
            app.MO_DoButton.Enable = 'on'; app.MO_UndoButton.Enable = 'off'; app.MO_SaveButton.Enable = 'off'; app.MO_statut_Label.Text = 'Idle';
            
            app.GrainvolumefilterinnumberofvoxelEditFieldLabel.Enable = 'off';
            app.GrainvolumefilterinnumberofvoxelEditField.Enable = 'off';            
            
            app.cellarray_secondassemble = app.cellarray_firstassemble;
            app.MO_Microstructure = app.cellarray_secondassemble; app.MO_Microstructure_saved = app.MO_Microstructure;  
            
            app.TabGroup.SelectedTab = app.MorphologyopeningTab; % Skip dimension tab
            
            app.SavecellButton.Text = 'Save cell';            
            app.Comp_Lamp.Color = 'g';
            app.TabisNOTsetupcorrectlyLabel_3.Text = 'Tab is setup correctly';            
        end

        % Cell edit callback: Dim_IP1_UITable
        function Dim_IP1_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if ~isnan(newData) && newData>=1 && round(newData)==newData
                if newData<=cell2mat(app.Dim_IP1_UITable.Data(indices(1),2)) && cell2mat(app.Dim_IP1_UITable.Data(indices(1),3))<=cell2mat(app.Dim_IP1_UITable.Data(indices(1),4))
                    foo=1; % Do nothing
                else % Error, re-initialize
                    app.Dim_IP1_UITable.Data(indices(1),3) = {1};
                    app.Dim_IP1_UITable.Data(indices(1),4) = app.Dim_IP1_UITable.Data(indices(1),2);
                end
            else % Error, re-initialize
                app.Dim_IP1_UITable.Data(indices(1),3) = {1};
                app.Dim_IP1_UITable.Data(indices(1),4) = app.Dim_IP1_UITable.Data(indices(1),2);
            end
            % Check dimension
            app.Check_inplane_dimension;            
        end

        % Cell edit callback: Dim_IP2_UITable
        function Dim_IP2_UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            if ~isnan(newData) && newData>=1 && round(newData)==newData
                if newData<=cell2mat(app.Dim_IP2_UITable.Data(indices(1),2)) && cell2mat(app.Dim_IP2_UITable.Data(indices(1),3))<=cell2mat(app.Dim_IP2_UITable.Data(indices(1),4))
                    foo=1; % Do nothing
                else % Error, re-initialize
                    app.Dim_IP2_UITable.Data(indices(1),3) = {1};
                    app.Dim_IP2_UITable.Data(indices(1),4) = app.Dim_IP2_UITable.Data(indices(1),2);
                end
            else % Error, re-initialize
                app.Dim_IP2_UITable.Data(indices(1),3) = {1};
                app.Dim_IP2_UITable.Data(indices(1),4) = app.Dim_IP2_UITable.Data(indices(1),2);
            end
            % Check dimension
            app.Check_inplane_dimension;             
        end

        % Button pushed function: MO_VisualizecellButton
        function MO_VisualizecellButtonPushed(app, event)
            direction_name = {'Normal to Through-plane direction';'Normal to In-plane direction 1';'Normal to In-plane direction 2'};
            Microstructure_comparison_visualization_interface(app.cellarray_secondassemble, app.MO_Microstructure_saved, app.MO_Microstructure, direction_name); 
        end

        % Button pushed function: MO_SavecellButton
        function MO_SavecellButtonPushed(app, event)
            % Remove background layer
            if app.micro_poly
                app.MO_Microstructure_saved(1,:,:) = [];
                app.MO_Microstructure_saved(end,:,:) = [];
                app.MO_Microstructure_saved(:,1,:) = [];
                app.MO_Microstructure_saved(:,end,:) = [];
                app.MO_Microstructure_saved(:,:,1) = [];
                app.MO_Microstructure_saved(:,:,end) = [];   
                app.domainbound(1).x1 = app.domainbound(1).x1-2;
            end
            
            % Save parameter
            if ~isempty(app.mainsavefoldersave)
                app.MO_SavecellButton.Text = 'Saving... please wait'; pause(0.1);
                if app.save_parameters
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Morphology_opening\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Morphology_opening/'];
                    end
                    % Prepare the data
                    clear DATA_writetable
                    d = app.MO_operations_UITable.Data;
                    DATA_writetable.sheet(1).name='Operations';
                    DATA_writetable.sheet(1).table=table(d(:,1),d(:,2),d(:,3),'VariableNames',{'Domain','Operation','Do'});
                    
                    DATA_writetable.sheet(2).name='Erosion_dilatation_parameters';
                    d={};
                    d(1,1)={'Erosion distance'}; d(1,2)={app.ErosiondistanceinvoxellengthEditField.Value};
                    d(2,1)={'Dilatation distance'}; d(2,2)={app.DilatationdistanceinvoxellengthEditField.Value};
                    d(3,1)={'Number of iteration'}; d(3,2)={app.ErosiondilatationNumberofiterationEditField.Value};
                    d(4,1)={'Perform on'}; d(4,2)={app.PerformErosionDropDown.Value};
                    d(5,1)={'Preserve background'}; d(5,2)={app.MO_PreservebackgroundCheckBox.Value};
                    DATA_writetable.sheet(2).table=table(d(:,1),d(:,2),'VariableNames',{'Parameters','Value'});
                    
                    DATA_writetable.sheet(3).name='Cleaning_voxel_connection_parameters';
                    d={};
                    d(1,1)={'Clean connection, perform on'}; d(1,2)={app.PerformCleanConnectionDropDown.Value};
                    d(2,1)={'Enforce pore background is a unique connected cluster'}; d(2,2)={app.MO_PoreuniqueclusterCheckBox.Value};
                    d(3,1)={'Enforce union of solid phase is a unique connected cluster'}; d(3,2)={app.MO_UnionsoliduniqueclusterCheckBox.Value};    
                    DATA_writetable.sheet(3).table=table(d(:,1),d(:,2),'VariableNames',{'Parameters','Value'});
                    
                    DATA_writetable.sheet(4).name='Volume_fractions';
                    DATA_writetable.sheet(4).table=table(app.MO_vf_UITable.Data(:,1),app.MO_vf_UITable.Data(:,2),app.MO_vf_UITable.Data(:,3),app.MO_vf_UITable.Data(:,4),'VariableNames',app.MO_vf_UITable.ColumnName);
                    
                    if app.micro_poly
                        DATA_writetable.sheet(5).name='Grains';
                        DATA_writetable.sheet(5).table=table(app.GrainvolumefilterinnumberofvoxelEditField.Value,app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Value,'VariableNames',{'Grain size filter','Grain contiguity'});
                    end
                    
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end
                    Function_Writetable(current_folder,'Morphology_opening',DATA_writetable)
                end
                if app.save_modified_tif
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Morphology_opening\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Morphology_opening/'];
                    end
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end                    
                    function_save_tif(uint8(app.MO_Microstructure_saved),[current_folder 'Cell.tif']);
                    n_domain=length(app.Domain_in_cell);
                    if app.micro_cell
                        for k=1:1:n_domain
                            [ str_new ] = function_remove_emptyandspecialcharacter_string(app.Domain_in_cell(k).name);
                            function_save_tif(uint8(app.MO_Microstructure_saved(app.domainbound(k).x0:app.domainbound(k).x1,:,:)),[current_folder str_new '.tif']);
                        end
                    end
                end
                app.MO_SavecellButton.Text = 'Save cell'; pause(0.1);
            end            
            
            % Save array
            app.cell_before_reassignemnt = app.MO_Microstructure_saved;
            
            % Update table
            app.update_celltable_id
                        
            % Next tab
            app.TabGroup.SelectedTab = app.AssemblecellTab;            

        end

        % Button pushed function: MO_DoButton
        function MO_DoButtonPushed(app, event)
            % Parameters
            n_domain=length(app.Domain_in_cell);
            erosion_distance = app.ErosiondistanceinvoxellengthEditField.Value;
            dilatation_distance = app.DilatationdistanceinvoxellengthEditField.Value;
            erositiondilatation_iteration = app.ErosiondilatationNumberofiterationEditField.Value;
            choice_erosiondilatation = app.PerformErosionDropDown.Value;
            preserve_background_erosiondilatation = app.MO_PreservebackgroundCheckBox.Value; % Not polychristalline architecture
            background_code = 0; % By default
            distance_method = 'chessboard';

            opt.background_code = background_code;
            opt.choice_cleansolidconnection = app.PerformCleanConnectionDropDown.Value;
            opt.choice_pore_uniquecluster = app.MO_PoreuniqueclusterCheckBox.Value;
            opt.choice_unionsolid_uniquecluster = app.MO_UnionsoliduniqueclusterCheckBox.Value;
            
            tmp = app.MO_operations_UITable.Data;
            [n_row, ~] = size(tmp);
            operation_choice = cell2mat(tmp(:,3));
            if app.micro_unique
                % Erosion dilatation
                if operation_choice(1)
                    app.MO_statut_Label.Text = 'Erosion + dilatation...'; pause(0.1);
                    [app.MO_Microstructure] = Function_morphologyopening_erosiondilatation(app.MO_Microstructure,background_code,erosion_distance,dilatation_distance,distance_method, erositiondilatation_iteration, choice_erosiondilatation, preserve_background_erosiondilatation);
                end
                % Correct microstructure
                if operation_choice(2)
                    if operation_choice(2) && operation_choice(3) % clean connection and convert to unique cluster
                        opt.todo_clean_connection = true; opt.todo_convert_to_uniquecluster = true;
                        app.MO_statut_Label.Text = 'Clean connection and convert to unique cluster...'; pause(0.1);
                        [app.MO_Microstructure] = Function_correct_microstructure(app.MO_Microstructure,opt);
                    elseif operation_choice(2) && ~operation_choice(3) % Only clean connection
                        opt.todo_clean_connection = true; opt.todo_convert_to_uniquecluster = false;
                        app.MO_statut_Label.Text = 'Clean connection...'; pause(0.1);
                        [app.MO_Microstructure] = Function_correct_microstructure(app.MO_Microstructure,opt);
                    elseif ~operation_choice(2) && operation_choice(3) % Only convert to unique cluster
                        opt.todo_clean_connection = false; opt.todo_convert_to_uniquecluster = true;
                        app.MO_statut_Label.Text = 'Convert to unique cluster...'; pause(0.1);
                        [app.MO_Microstructure] = Function_correct_microstructure(app.MO_Microstructure,opt);
                    end
                    if ~operation_choice(3) % At minimum, isolated cluster of size 1 are assigned to adjacent phases
                        sizefilter=1;
                        [app.MO_Microstructure] = Function_remove_isolated_clusters(app.MO_Microstructure, background_code, sizefilter);
                    end
                end
                
            elseif app.micro_cell
                % Erosion dilatation
                line=0;
                for kdomain=1:1:n_domain
                    for kstep=1:1:app.n_mo_step
                        line=line+1;
                        if strcmp(char(tmp(line,2)),'Erosion + dilatation') && operation_choice(line)
                            app.MO_statut_Label.Text = 'Erosion + dilatation...'; pause(0.1);
                            tmp2 = app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:);
                            [tmp2] = Function_morphologyopening_erosiondilatation(tmp2,background_code,erosion_distance,dilatation_distance,distance_method, erositiondilatation_iteration, choice_erosiondilatation, preserve_background_erosiondilatation);
                            app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp2;
                        end
                    end
                end
                % Correct microstructure
                line=0;
                for kdomain=1:1:n_domain
                    tmp2 = app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:);
                    for kstep=1:1:app.n_mo_step
                        line=line+1;
                        if line<n_row
                            if strcmp(char(tmp(line,2)),'Clean voxel connection')
                                if operation_choice(line) && operation_choice(line+1) % clean connection and convert to unique cluster
                                    opt.todo_clean_connection = true; opt.todo_convert_to_uniquecluster = true;
                                    app.MO_statut_Label.Text = 'Clean connection and convert to unique cluster...'; pause(0.1);
                                    [tmp2] = Function_correct_microstructure(tmp2,opt);
                                    app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp2;
                                elseif operation_choice(line) && ~operation_choice(line+1) % Only clean connection
                                    opt.todo_clean_connection = true; opt.todo_convert_to_uniquecluster = false;
                                    app.MO_statut_Label.Text = 'Clean connection...'; pause(0.1);
                                    [tmp2] = Function_correct_microstructure(tmp2,opt);
                                    app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp2;                                    
                                elseif ~operation_choice(line) && operation_choice(line+1) % Only convert to unique cluster
                                    opt.todo_clean_connection = false; opt.todo_convert_to_uniquecluster = true;
                                    app.MO_statut_Label.Text = 'Convert to unique cluster...'; pause(0.1);
                                    [tmp2] = Function_correct_microstructure(tmp2,opt);
                                    app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp2;
                                end
                                if ~operation_choice(line+1) % At minimum, isolated cluster of size 1 are assigned to adjacent phases
                                    sizefilter=1;
                                    [tmp2] = Function_remove_isolated_clusters(tmp2, background_code, sizefilter);
                                    app.MO_Microstructure(app.domainbound(kdomain).x0:app.domainbound(kdomain).x1,:,:) = tmp2;
                                end
                            end
                        end
                    end
                end
                

            elseif app.micro_poly
                check_contiguity = app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Value;
                remove_isolatedpore = app.MO_PoreuniqueclusterCheckBox.Value;
                grainsizefilter = app.GrainvolumefilterinnumberofvoxelEditField.Value;
                
                % Steps below could be re-order
                
                % Re-order grain
                ids = unique(app.MO_Microstructure);
                ids(ids==background_code)=[];
                tmp = zeros(size(app.MO_Microstructure));
                for k=1:1:length(ids)
                    tmp(app.MO_Microstructure==ids(k))=k;
                end
                app.MO_Microstructure = tmp; clear tmp;                
                % Grain connectivity
                if check_contiguity
                    % Grain dicontinuity is corrected (if some grains are not continguous, new id will be assinged to them)
                    connectivity = 6;
                    app.MO_statut_Label.Text = 'Check grain contiguity (1/3)...'; pause(0.1);
                    [app.MO_Microstructure] = Function_grain_connectivity(app.MO_Microstructure,background_code,connectivity);
                end
                if remove_isolatedpore && operation_choice(3)
                    app.MO_statut_Label.Text = 'Remove isolated pore/cracks...'; pause(0.1);
                    [app.MO_Microstructure] = Function_convert_to_unique_cluster(app.MO_Microstructure,background_code, true, false, true);
                end
                if grainsizefilter>0
                    app.MO_statut_Label.Text = 'Apply grain size filter (1/2)...'; pause(0.1);
                    [app.MO_Microstructure] = Function_grain_sizefilter(app.MO_Microstructure,background_code,grainsizefilter);
                end
                % Erosion dilatation
                if operation_choice(1)
                    app.MO_statut_Label.Text = 'Erosion + dilatation...'; pause(0.1);
                    [app.MO_Microstructure] = Function_morphologyopening_erosiondilatation(app.MO_Microstructure,background_code,erosion_distance,dilatation_distance,distance_method, erositiondilatation_iteration, choice_erosiondilatation, preserve_background_erosiondilatation);
                end
                 % Grain connectivity
                if check_contiguity
                    % Grain dicontinuity is corrected (if some grains are not continguous, new id will be assinged to them)
                    app.MO_statut_Label.Text = 'Check grain contiguity (2/3)...'; pause(0.1);
                    [app.MO_Microstructure] = Function_grain_connectivity(app.MO_Microstructure,background_code,connectivity);
                end               
                % Clean connection 
                if operation_choice(2)
                    app.MO_statut_Label.Text = 'Clean connection...'; pause(0.1);
                    opt.todo_clean_connection = true; opt.todo_convert_to_uniquecluster = false;
                    [app.MO_Microstructure] = Function_correct_microstructure(app.MO_Microstructure,opt);
                end
                % Grain connectivity
                if check_contiguity
                    % Grain dicontinuity is corrected (if some grains are not continguous, new id will be assinged to them)
                    app.MO_statut_Label.Text = 'Check grain contiguity (3/3)...'; pause(0.1);
                    [app.MO_Microstructure] = Function_grain_connectivity(app.MO_Microstructure,background_code,connectivity);
                end                    
                if grainsizefilter>0
                    app.MO_statut_Label.Text = 'Apply grain size filter (2/2)...'; pause(0.1);
                    [app.MO_Microstructure] = Function_grain_sizefilter(app.MO_Microstructure,background_code,grainsizefilter);
                end 
                if remove_isolatedpore && operation_choice(3)
                    app.MO_statut_Label.Text = 'Remove isolated pore/cracks...'; pause(0.1);
                    [app.MO_Microstructure] = Function_convert_to_unique_cluster(app.MO_Microstructure,background_code, true, false, true);
                end                
                
                % Re-order grain
                ids = unique(app.MO_Microstructure);
                ids(ids==background_code)=[];
                tmp = zeros(size(app.MO_Microstructure));
                for k=1:1:length(ids)
                    tmp(app.MO_Microstructure==ids(k))=k;
                end
                app.MO_Microstructure = tmp; clear tmp;                 
            end
            app.MO_statut_Label.Text = 'Done!';
            
            % Update volume fraction
            if ~app.micro_poly
                for k=1:1:n_domain
                    tmp = app.MO_Microstructure(app.domainbound(k).x0:app.domainbound(k).x1,:,:);
                    phases = unique(tmp); n=length(phases);
                    if n>1
                        for p=1:1:n
                            vf=sum(sum(sum( tmp==phases(p)  )))/numel(tmp);
                            if p==1
                                tmp3 = [num2str(vf,'%1.3f') ' / '];
                            elseif p==n
                                tmp3 = [tmp3 num2str(vf,'%1.3f')];
                            else
                                tmp3 = [tmp3 num2str(vf,'%1.3f') ' / '];
                            end
                        end
                    else
                        vf=sum(sum(sum( tmp==phases(1)  )))/numel(tmp);
                        tmp3 = [num2str(vf,'%1.3f')];
                    end
                    vf_after(k,1)={tmp3};
                end
                app.MO_vf_UITable.Data(:,4) = vf_after;
                
            else
                tmp = app.MO_Microstructure;
                grain_id = unique(tmp); n=length(grain_id);
                [nbefore,~] = size(app.MO_vf_UITable.Data(:,3));
                vf_after = cell(nbefore,1);
                for p=1:1:n
                    gid(p,1) = {grain_id(p)};
                    vf_after(p,1) = {sum(sum(sum( tmp==grain_id(p) ))) / numel(tmp)};
                end
                app.MO_vf_UITable.Data(:,4) = vf_after;
            end
            
            app.MO_DoButton.Enable = 'off'; app.MO_UndoButton.Enable = 'on'; app.MO_UndoAllButton.Enable = 'on'; app.MO_SaveButton.Enable = 'on';  app.MO_SavecellButton.Enable = 'off';
        end

        % Button pushed function: MO_UndoButton
        function MO_UndoButtonPushed(app, event)
            app.MO_statut_Label.Text = 'Reverse previous change...'; pause(0.1);
            app.MO_Microstructure = app.MO_Microstructure_saved;
            app.MO_DoButton.Enable = 'on'; app.MO_UndoButton.Enable = 'off'; app.MO_UndoAllButton.Enable = 'on'; app.MO_SaveButton.Enable = 'off'; app.MO_SavecellButton.Enable = 'on';
            % Update volume fraction
            if ~app.micro_poly
                n_domain=length(app.Domain_in_cell);
                for k=1:1:n_domain
                    tmp = app.MO_Microstructure(app.domainbound(k).x0:app.domainbound(k).x1,:,:);
                    phases = unique(tmp); n=length(phases);
                    if n>1
                        for p=1:1:n
                            vf=sum(sum(sum( tmp==phases(p)  )))/numel(tmp);
                            if p==1
                                tmp3 = [num2str(vf,'%1.3f') ' / '];
                            elseif p==n
                                tmp3 = [tmp3 num2str(vf,'%1.3f')];
                            else
                                tmp3 = [tmp3 num2str(vf,'%1.3f') ' / '];
                            end
                        end
                    else
                        vf=sum(sum(sum( tmp==phases(1)  )))/numel(tmp);
                        tmp3 = [num2str(vf,'%1.3f')];
                    end
                    vf_after(k,1)={tmp3};
                end
            else
                tmp = app.MO_Microstructure;
                grain_id = unique(tmp); n=length(grain_id);
                [nbefore,~] = size(app.MO_vf_UITable.Data(:,3));
                vf_after = cell(nbefore,1);
                for p=1:1:n
                    gid(p,1) = {grain_id(p)};
                    vf_after(p,1) = {sum(sum(sum( tmp==grain_id(p) ))) / numel(tmp)};
                end
            end
            app.MO_vf_UITable.Data(:,4) = vf_after;
            app.MO_statut_Label.Text = 'Done!';
        end

        % Button pushed function: MO_SaveButton
        function MO_SaveButtonPushed(app, event)
            app.MO_statut_Label.Text = 'Saving for next tab...';  pause(0.1);
            app.MO_Microstructure_saved = app.MO_Microstructure;
            app.MO_DoButton.Enable = 'on'; app.MO_UndoButton.Enable = 'off'; app.MO_UndoAllButton.Enable = 'on'; app.MO_SaveButton.Enable = 'off'; app.MO_SavecellButton.Enable = 'on';           
            app.MO_statut_Label.Text = 'Done!';
        end

        % Button pushed function: MO_SkipButton
        function MO_SkipButtonPushed(app, event)
            if app.micro_poly
                % Re-order grain at least
                background_code = 0;
                ids = unique(app.MO_Microstructure);
                ids(ids==background_code)=[];
                tmp = zeros(size(app.MO_Microstructure));
                for k=1:1:length(ids)
                    tmp(app.MO_Microstructure==ids(k))=k;
                end
                app.MO_Microstructure = tmp; clear tmp;
                
                % Remove background layer
                app.MO_Microstructure_saved(1,:,:) = [];
                app.MO_Microstructure_saved(end,:,:) = [];
                app.MO_Microstructure_saved(:,1,:) = [];
                app.MO_Microstructure_saved(:,end,:) = [];
                app.MO_Microstructure_saved(:,:,1) = [];
                app.MO_Microstructure_saved(:,:,end) = [];
                app.domainbound(1).x1 = app.domainbound(1).x1-2;
            end
            app.MO_statut_Label.Text = 'Saving for next tab...';  pause(0.1);
            app.cell_before_reassignemnt = app.cellarray_secondassemble;
            app.MO_statut_Label.Text = 'Done!';
            % Next tab
            if app.micro_cell || app.micro_unique
                % Update table
                app.update_celltable_id
                app.TabGroup.SelectedTab = app.AssemblecellTab;
            else
                app.TabGroup.SelectedTab = app.Meshgenerationchoice; % Skip assemble tab
            end            
        end

        % Button pushed function: MO_UndoAllButton
        function MO_UndoAllButtonPushed(app, event)
            app.MO_statut_Label.Text = 'Reverse all changes...'; pause(0.1);
            app.MO_Microstructure = app.cellarray_secondassemble; app.MO_Microstructure_saved = app.cellarray_secondassemble;
            app.MO_DoButton.Enable = 'on'; app.MO_UndoButton.Enable = 'off'; app.MO_UndoAllButton.Enable = 'off'; app.MO_SaveButton.Enable = 'off'; app.MO_SavecellButton.Enable = 'on';
            % Update volume fraction
            if ~app.micro_poly
            n_domain=length(app.Domain_in_cell);
            for k=1:1:n_domain
                tmp = app.MO_Microstructure(app.domainbound(k).x0:app.domainbound(k).x1,:,:);
                phases = unique(tmp); n=length(phases);
                if n>1
                    for p=1:1:n
                        vf=sum(sum(sum( tmp==phases(p)  )))/numel(tmp);
                        if p==1
                            tmp3 = [num2str(vf,'%1.3f') ' / '];
                        elseif p==n
                            tmp3 = [tmp3 num2str(vf,'%1.3f')];
                        else
                            tmp3 = [tmp3 num2str(vf,'%1.3f') ' / '];
                        end
                    end
                else
                    vf=sum(sum(sum( tmp==phases(1)  )))/numel(tmp);
                    tmp3 = [num2str(vf,'%1.3f')];
                end
                vf_after(k,1)={tmp3};
            end
            else
                tmp = app.MO_Microstructure;
                grain_id = unique(tmp); n=length(grain_id);
                [nbefore,~] = size(app.MO_vf_UITable.Data(:,3));
                vf_after = cell(nbefore,1);
                for p=1:1:n
                    gid(p,1) = {grain_id(p)};
                    vf_after(p,1) = {sum(sum(sum( tmp==grain_id(p) ))) / numel(tmp)};
                end
            end
            app.MO_vf_UITable.Data(:,4) = vf_after;
            app.MO_statut_Label.Text = 'Done!';
        end

        % Cell edit callback: Assemble_UITable
        function Assemble_UITableCellEdit(app, event)
            app.Assemble_VisualizecellButton.Enable = 'off';
            app.Assemble_SavecellButton.Enable = 'off';
            app.Assemble_Lamp.Color = 'r';
            app.TabisNOTsetupcorrectlyLabel_4.Text = 'Tab is NOT setup correctly';
        end

        % Button pushed function: Assemble_VisualizecellButton
        function Assemble_VisualizecellButtonPushed(app, event)
            direction_name = {'Normal to Through-plane direction';'Normal to In-plane direction 1';'Normal to In-plane direction 2'};
            Microstructure_comparison_visualization_interface(app.cell_before_reassignemnt, app.cell_after_reassignement, app.cell_group_after_reassignement, direction_name)
        end

        % Button pushed function: Assemble_SavecellButton
        function Assemble_SavecellButtonPushed(app, event)
            % Save parameter
            if ~isempty(app.mainsavefoldersave)
                app.Assemble_SavecellButton.Text = 'Saving... please wait'; pause(0.1);
                if app.save_parameters || app.save_modified_tif
                    if ispc
                        current_folder= [app.mainsavefoldersave 'Parameters\Assemble\'];
                    else
                        current_folder = [app.mainsavefoldersave 'Parameters/Assemble/'];
                    end
                    if exist(current_folder,'dir')==0 % Folder existence is checked, and created if necessary
                        mkdir(current_folder);
                    end                    
                end
                
                if app.save_parameters
                    % Prepare the data
                    clear DATA_writetable
                    d = app.Assemble_UITable.Data;
                    DATA_writetable.sheet(1).name='Assignment';
                    DATA_writetable.sheet(1).table=table(d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),'VariableNames',{'Domain','Phase id','Assigned id in cell','Name','Group id'});
                    Function_Writetable(current_folder,'Assemble',DATA_writetable)
                end
                if app.save_modified_tif
                    function_save_tif(uint8(app.cell_after_reassignement),[current_folder 'Cell.tif']);
                    function_save_tif(uint8(app.cell_group_after_reassignement),[current_folder 'Group.tif']);
                end
                app.Assemble_SavecellButton.Text = 'Save cell'; pause(0.1);
            end
            if app.micro_cell % Optional
                % Microstructure dimension: [voxel size; number of voxel along axe 1; number of voxel along axe 2; number of voxel along axe 3];
                app.Microstructure_dimension = [app.VoxelsizeEditField.Value; size(app.cell_after_reassignement)'];
                % Separator bounds
                tmp = app.Assemble_UITable.Data;
                cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
                logical_cells = cellfun(cellfind('Separator'),tmp(:,1));
                id_separator = find(logical_cells);
                if ~isempty(id_separator)
                    xminmax=[1e9;-1e9]; % Initialization
                    for k=1:1:length(id_separator)
                        idx = find(app.cell_after_reassignement==cell2mat(tmp(id_separator(k),3)));
                        [IX, ~, ~]=ind2sub(size(app.cell_after_reassignement),idx);
                        xminmax = [min(xminmax(1),min(IX)-1);max(xminmax(2),max(IX))];
                    end
                    app.Separator_bounds = xminmax;
                    %idx = find(app.cell_after_reassignement==cell2mat(tmp(id_separator,3)));
                    %[IX, ~, ~]=ind2sub(size(app.cell_after_reassignement),idx);
                    %app.Separator_bounds = [min(IX)-1; max(IX)];
                end
            end
            
            app.Assemble_Lamp.Color = 'g';
            app.TabisNOTsetupcorrectlyLabel_4.Text = 'Tab is setup correctly';             
                        
            % Next tab
            app.Meshgen_Iso2mesh.Enable = 'on';
            app.Meshgen_Cuboidmesh.Enable = 'on';
            app.TabGroup.SelectedTab = app.Meshgenerationchoice;
        end

        % Button pushed function: UpdateButton
        function UpdateButtonPushed(app, event)
            % Check entry are correct
            new_id = cell2mat(app.Assemble_UITable.Data(:,3));
            group_id = cell2mat(app.Assemble_UITable.Data(:,5));
            
            id_are_correct = true;
            if max(isnan(new_id)) || min(new_id)<1 || length(unique(new_id))~=numel(new_id)
                id_are_correct = false;
            end
            if max(isnan(group_id)) || min(group_id)<0
                id_are_correct = false;
            end     
      
            % Update array
            if id_are_correct
                % Group name table
                u = unique(group_id);
                c1={}; c2={};
                for kg = 1:1:length(u)
                    c1(kg,1) = {u(kg)};
                    c2(kg,1) ={['Group_' num2str(u(kg))]};
                end
                app.Group_UITable.Data = [c1 c2];
                app.UpdateButton_label.Text = 'Once table is modified, please click on the update button';
                app.UpdateButton_label.FontColor = 'k';
                app.update_cellid
            else
                app.UpdateButton_label.Text = 'Id are NOT correct!';  
                app.UpdateButton_label.FontColor = 'r';
                app.Assemble_VisualizecellButton.Enable = 'off';
                app.Assemble_SavecellButton.Enable = 'off';                
            end            
        end

        % Button pushed function: Meshgen_Iso2mesh
        function Meshgen_Iso2meshButtonPushed(app, event)
            app.TabGroup.SelectedTab = app.Iso2meshoptionsTab;
            app.Createmesh_button.Enable = 'on';
            app.Createmesh_choicemesh.Text = 'Mesh method is: Iso2mesh';
            app.Createmesh_choicemesh.FontColor = 'k';
            app.meshmethod = 'Iso2mesh';
            app.Face_CheckBox.Value = true; app.Face_CheckBox.Enable = false;
        end

        % Button pushed function: Meshgen_Cuboidmesh
        function Meshgen_CuboidmeshButtonPushed(app, event)
            app.TabGroup.SelectedTab = app.StructuredmeshoptionsTab;
            app.Createmesh_button.Enable = 'on';
            app.Createmesh_choicemesh.Text = 'Mesh method is: Cuboid mesh';
            app.Createmesh_choicemesh.FontColor = 'k';        
            app.meshmethod = 'Cuboid';
            app.Face_CheckBox.Value = true; app.Face_CheckBox.Enable = true;
        end

        % Value changed function: 
        % EachvoxelwillberepresentedwithDropDown
        function EachvoxelwillberepresentedwithDropDownValueChanged(app, event)
            value = app.EachvoxelwillberepresentedwithDropDown.Value;
            if strcmp(value,'5 tetrahedrons') % Can generate the mesh but FEM calculations bug for this one, so there is sth to fix
                app.Regularmesh_Image6.Visible = 'off'; app.Regularmesh_Image24.Visible = 'off';
                app.Regularmesh_Image5.Visible = 'on';
            elseif strcmp(value,'6 tetrahedrons')
                app.Regularmesh_Image5.Visible = 'off'; app.Regularmesh_Image24.Visible = 'off';
                app.Regularmesh_Image6.Visible = 'on';                
            elseif strcmp(value,'24 tetrahedrons')
                app.Regularmesh_Image5.Visible = 'off'; app.Regularmesh_Image6.Visible = 'off';
                app.Regularmesh_Image24.Visible = 'on';
            end            
        end

        % Button pushed function: Createmesh_button
        function Createmesh_buttonButtonPushed(app, event)
            app.Mesh_VisualizeButton.Enable = 'off';
            app.Mesh_SaveButton.Enable = 'off';         
            app.Meshquality_JoeLiu_button.Enable = 'off';
            tmp=app.Meshchoice_UITable.Data; tmp=tmp.Variables;
            app.do_one_mesh_per_phase = cell2mat(tmp(1,2));
            app.do_one_mesh_per_group = cell2mat(tmp(2,2));
            app.do_one_mesh_wholegeometry = cell2mat(tmp(3,2));
            tmp=app.Assemble_UITable.Data;
            phase_id = cell2mat(tmp(:,3)); number_phase = length(phase_id);
            group_id = cell2mat(tmp(:,5)); unique_group = unique(group_id); number_group = length(unique_group);
            phase_name = tmp(:,4); 
            tmp=app.Group_UITable.Data;
            group_name = tmp(:,2); 
            app.meshgroup(1).number_group = number_group;
            app.meshphase(1).number_phase = number_phase;
            for kgroup = 1:1:number_group
                app.meshgroup(kgroup).name = char(group_name(kgroup));
            end
            for kphase = 1:1:number_phase
                app.meshphase(kphase).name = char(phase_name(kphase));
            end
            app.meshvolume.name = 'Full volume';
            do_face = app.Face_CheckBox.Value;
            app.Createmesh_StatutTextArea.Value = 'Generating meshes...'; pause(0.1)
            if strcmp(app.meshmethod,'Cuboid')
                app.meshoptions.method = 'Cuboid';
                cellchoice = app.EachvoxelwillberepresentedwithDropDown.Value;
                app.meshoptions.cellchoice = cellchoice;
                if app.do_one_mesh_per_phase
                    for kphase = 1:1:number_phase
                        app.Createmesh_StatutTextArea.Value = ['Generating mesh of phase #' num2str(kphase) ': ' char(phase_name(kphase))]; pause(0.1);
                        current_phase = phase_id(kphase);
                        [node_, face_, elem_, subdomain_] = function_regularmesh_from_array(app.cell_after_reassignement,current_phase,cellchoice,do_face);
                        app.meshphase(kphase).node_ = node_;
                        app.meshphase(kphase).face_ = face_;
                        app.meshphase(kphase).elem_ = elem_;
                        app.meshphase(kphase).subdomain_ = subdomain_;
                    end                       
                end
                if app.do_one_mesh_per_group
                    for kgroup = 1:1:number_group
                        app.Createmesh_StatutTextArea.Value = ['Generating mesh of group #' num2str(kgroup) ': ' char(group_name(kgroup))]; pause(0.1);
                        current_group = unique_group(kgroup);
                        idx = find(group_id==current_group);
                        ids = phase_id(idx);
                        [node_, face_, elem_, subdomain_] = function_regularmesh_from_array(app.cell_after_reassignement,ids,cellchoice,do_face);
                        app.meshgroup(kgroup).node_ = node_;
                        app.meshgroup(kgroup).face_ = face_;
                        app.meshgroup(kgroup).elem_ = elem_;
                        app.meshgroup(kgroup).subdomain_ = subdomain_;
                    end                    
                end
                if app.do_one_mesh_wholegeometry
                    app.Createmesh_StatutTextArea.Value = 'Generating mesh of full volume'; pause(0.1);
                    [node_, face_, elem_, subdomain_] = function_regularmesh_from_array(app.cell_after_reassignement,phase_id,cellchoice,do_face);
                    app.meshvolume.node_ = node_;
                    app.meshvolume.face_ = face_;
                    app.meshvolume.elem_ = elem_;
                    app.meshvolume.subdomain_ = subdomain_;
                end
                
            elseif strcmp(app.meshmethod,'Iso2mesh')
                % Read parameters from GUI
                tmp = app.CGAL_UITable.Data;
                app.meshoptions.method = 'Iso2mesh';
                app.meshoptions.radbound  = cell2mat(tmp(1,3));
                app.meshoptions.distbound = cell2mat(tmp(2,3));
                
                app.meshoptions.method_surfacemesh = app.Surfacemesh_MethodDropDown.Value;
                tmp = app.Smoothing_UITable.Data;
                app.meshoptions.iteration_smoothing  = cell2mat(tmp(1,3));
                app.meshoptions.useralpha = cell2mat(tmp(2,3));
                if strcmp(app.meshoptions.method_surfacemesh,'laplacianhc')
                    app.meshoptions.userbeta = cell2mat(tmp(3,3));
                end
                
                tmp = app.Meshdensity_UITable.Data;
                app.meshoptions.keepratio  = cell2mat(tmp(1,3));
                app.meshoptions.maxvol  = cell2mat(tmp(2,3));
                
                % Create mesh (full volume)
                [node_, face_, elem_, subdomain_] = function_iso2mesh_from_array(app.cell_after_reassignement,app.meshoptions);
                app.meshvolume.node_ = node_;
                app.meshvolume.face_ = face_;
                app.meshvolume.elem_ = elem_;
                app.meshvolume.subdomain_ = subdomain_;
                
                if app.do_one_mesh_per_phase
                    for kphase = 1:1:number_phase
                        current_phase = phase_id(kphase);
                        idx = find(subdomain_==current_phase);
                        vertices_idx_phase = unique(elem_(idx,:));
                        vertices_phase = node_(vertices_idx_phase,:);
                        check_face = ismember(face_,vertices_idx_phase);
                        faces_idx_phase = find( sum(check_face,2)==3 );
                        
                        %[Lia,Locb] = ismember(round(node_,5), round(vertices_group,5),'rows'); % Precision with 5 digits (tolerance)
                        [Lia,Locb] = ismembertol(node_, vertices_phase,1e-6,'ByRows',true);
                        old_id = find(Lia==1);
                        new_id = Locb(old_id);
                        
                        app.meshphase(kphase).node_ = vertices_phase;
                        app.meshphase(kphase).face_ = my_changem(face_(faces_idx_phase,:), new_id, old_id);
                        app.meshphase(kphase).elem_ = my_changem(elem_(idx,:), new_id, old_id);
                        app.meshphase(kphase).subdomain_ = subdomain_(idx);
                    end
                    
                end
                if app.do_one_mesh_per_group
                    for kgroup = 1:1:number_group
                        current_group = unique_group(kgroup);
                        idx = find(group_id==current_group);
                        ids = phase_id(idx);
                        idx=[];
                        for k=1:1:length(ids)
                            idx = [idx; find(subdomain_==ids(k))];
                        end
                        vertices_idx_group = unique(elem_(idx,:));
                        vertices_group = node_(vertices_idx_group,:);
                        check_face = ismember(face_,vertices_idx_group);
                        faces_idx_group = find( sum(check_face,2)==3 );
                        
                        %[Lia,Locb] = ismember(round(node_,5), round(vertices_group,5),'rows'); % Precision with 5 digits (tolerance)
                        [Lia,Locb] = ismembertol(node_, vertices_group,1e-6,'ByRows',true);
                        old_id = find(Lia==1);
                        new_id = Locb(old_id);
                        
                        app.meshgroup(kgroup).node_ = vertices_group;
                        app.meshgroup(kgroup).face_ = my_changem(face_(faces_idx_group,:), new_id, old_id);
                        app.meshgroup(kgroup).elem_ = my_changem(elem_(idx,:), new_id, old_id);
                        app.meshgroup(kgroup).subdomain_ = subdomain_(idx);
                    end
                end
                
                
            end
            app.Createmesh_StatutTextArea.Value = 'Done';
            
            if app.Face_CheckBox.Value
                app.Mesh_VisualizeButton.Enable = 'on';
            end
            app.Mesh_SaveButton.Enable = 'on';
            app.Meshquality_JoeLiu_button.Enable = 'on';
        end

        % Button pushed function: Mesh_VisualizeButton
        function Mesh_VisualizeButtonPushed(app, event)
            options.show = app.ShowDropDown.Value;
            options.lighting = app.LightingDropDown.Value;
            options.edges = app.MeshedgesDropDown.Value;
            options.coloris = app.ColorisDropDown.Value;
            options.colormap = app.ColormapDropDown.Value;
            options.transparency = app.Transparency01EditField.Value;
            options.grid = app.GridCheckBox.Value;
            options.axislabel = app.AxislabelCheckBox.Value;
            options.savefig = app.SavefigCheckBox.Value;
            options.savepng = app.SavepngCheckBox.Value;
            options.indexstart_zero = app.indexstart_zero;
            if ispc
                options.folder= [app.mainsavefoldersave 'Results\Pictures\'];
            else
                options.folder = [app.mainsavefoldersave 'Results/Pictures/'];
            end
            if (options.savefig || options.savepng) && exist(options.folder,'dir')==0 % Folder existence is checked, and created if necessary
                mkdir(options.folder);
            end
            if app.do_one_mesh_per_phase
                for kphase = 1:1:app.meshphase(1).number_phase
                    function_showmesh(app.meshphase(kphase).node_, app.meshphase(kphase).elem_, app.meshphase(kphase).face_, app.meshphase(kphase).subdomain_, app.meshphase(kphase).name, options);
                end
            end
            if app.do_one_mesh_per_group
                for kgroup = 1:1:app.meshgroup(1).number_group
                    function_showmesh(app.meshgroup(kgroup).node_, app.meshgroup(kgroup).elem_, app.meshgroup(kgroup).face_, app.meshgroup(kgroup).subdomain_, app.meshgroup(kgroup).name, options);
                end
            end
            if app.do_one_mesh_wholegeometry
                function_showmesh(app.meshvolume.node_, app.meshvolume.elem_, app.meshvolume.face_, app.meshvolume.subdomain_, app.meshvolume.name, options);
            end
        end

        % Button pushed function: Mesh_SaveButton
        function Mesh_SaveButtonPushed(app, event)
            if app.savefolder_is_defined
                if ispc
                    options.folder= [app.mainsavefoldersave 'Results\Meshes\'];
                else
                    options.folder = [app.mainsavefoldersave 'Results/Meshes/'];
                end
                if exist(options.folder,'dir')==0 % Folder existence is checked, and created if necessary
                    mkdir(options.folder);
                end
                
                tmp = app.Meshformat_UITable.Data.Variables;
                choice = cell2mat(tmp(:,2));
                options.save_mat = choice(1);
                options.save_csv = choice(2);
                options.save_msh = choice(3);
                options.save_inp = choice(4);
                options.save_stl = choice(5);
                options.save_binarystl = choice(6);
                options.indexstart_zero = app.Cellindexstartsat0CheckBox.Value;
                
                % Save parameters
                if app.save_parameters
                    % Prepare the data
                    clear DATA_writetable
                    if strcmp(app.meshoptions.method,'Cuboid')
                        DATA_writetable.sheet(1).name='Cuboid';
                        DATA_writetable.sheet(1).table=table({app.meshoptions.cellchoice},'VariableNames',{'Cell choice'});
                        Function_Writetable(options.folder,'Meshoptions',DATA_writetable)                        
                    elseif strcmp(app.meshoptions.method,'Iso2mesh')
                        DATA_writetable.sheet(1).name='CGAL';
                        DATA_writetable.sheet(1).table=table(app.meshoptions.radbound,app.meshoptions.distbound,'VariableNames',{'radbound','distbound'});
                        DATA_writetable.sheet(2).name='Smoothing';
                        if strcmp(app.meshoptions.method_surfacemesh,'laplacianhc')
                            DATA_writetable.sheet(2).table=table({app.meshoptions.method_surfacemesh},app.meshoptions.iteration_smoothing,app.meshoptions.useralpha,app.meshoptions.userbeta,'VariableNames',{'Method','Iteration number','useralpha','userbeta'});  
                        else
                            DATA_writetable.sheet(2).table=table({app.meshoptions.method_surfacemesh},app.meshoptions.iteration_smoothing,app.meshoptions.useralpha,'VariableNames',{'Method','Iteration number','useralpha'});  
                        end
                        DATA_writetable.sheet(3).name='Mesh density';
                        DATA_writetable.sheet(3).table=table(app.meshoptions.keepratio,app.meshoptions.maxvol,'VariableNames',{'keepratio','maxvol'});
                        Function_Writetable(options.folder,'Meshoptions',DATA_writetable)                           
                    end
                end
                
                app.Createmesh_StatutTextArea.Value = 'Saving...'; pause(0.1);
                if app.do_one_mesh_per_phase
                    for kphase = 1:1:app.meshphase(1).number_phase
                        function_savemesh(app.meshphase(kphase).node_, app.meshphase(kphase).face_, app.meshphase(kphase).elem_, app.meshphase(kphase).subdomain_, app.meshphase(kphase).name, options)
                    end                    
                end
                if app.do_one_mesh_per_group
                    for kgroup = 1:1:app.meshgroup(1).number_group
                        function_savemesh(app.meshgroup(kgroup).node_, app.meshgroup(kgroup).face_, app.meshgroup(kgroup).elem_, app.meshgroup(kgroup).subdomain_, app.meshgroup(kgroup).name, options)
                    end
                end
                if app.do_one_mesh_wholegeometry
                     function_savemesh(app.meshvolume.node_, app.meshvolume.face_, app.meshvolume.elem_, app.meshvolume.subdomain_, app.meshvolume.name, options)
                end
                
                if app.micro_cell
                    % Microstructure dimension: [voxel size in micrometers; number of voxel along axe 1; number of voxel along axe 2; number of voxel along axe 3]
                    tmp = app.Microstructure_dimension; save([options.folder 'Microstructure_dimension.txt'],"tmp",'-ascii');
                    % dlmwrite([options.folder 'Microstructure_dimension.txt'],app.Microstructure_dimension)
                    % Separator bound
                    if ~isempty(app.Separator_bounds)
                        tmp = app.Separator_bounds; save([options.folder 'Separator_bounds.txt'],"tmp",'-ascii');
                        % dlmwrite([options.folder 'Separator_bounds.txt'],app.Separator_bounds);
                    end
                end
                app.Createmesh_StatutTextArea.Value = 'Done';
            else
                app.TabGroup.SelectedTab = app.FolderandsaveoptionsTab; % Set saving folder
            end
            
        end

        % Button pushed function: Meshquality_JoeLiu_button
        function Meshquality_JoeLiu_buttonPushed(app, event)
            parameters_distributionfigure.figureposition = [100 100 1500 800];
            parameters_distributionfigure.fontname = 'Times New Roman';
            parameters_distributionfigure.grid = 'on';
            parameters_distributionfigure.minorgrid = 'on';
            parameters_distributionfigure.fullpath = 1;
            parameters_distributionfigure.save = app.savefolder_is_defined;
            parameters_distributionfigure.closefig = false;
            if parameters_distributionfigure.save
                if ispc
                    parameters_distributionfigure.fullpath= [app.mainsavefoldersave 'Results\Mesh_quality\'];
                else
                    parameters_distributionfigure.fullpath = [app.mainsavefoldersave 'Results/Mesh_quality/'];
                end
                if exist(parameters_distributionfigure.fullpath,'dir')==0 % Folder existence is checked, and created if necessary
                    mkdir(parameters_distributionfigure.fullpath);
                end
            end
            parameters_smoothingmovingaverage.smooth_cumulative_fct = true;
            parameters_smoothingmovingaverage.minimum_array_length = 5;
            parameters_smoothingmovingaverage.number_point = 250;
            parameters_smoothingmovingaverage.moving_rangeratio = 0.05;
            parameters_smoothingmovingaverage.moving_range = 0;
            parameters_smoothingmovingaverage.enforce_samelength = false;
            parameters_smoothingmovingaverage.origin = 'symmetrical';
            parameters_smoothingmovingaverage.boundary_behavior = 'keep origin';
            parameters_smoothingmovingaverage.round_value = 2;

            if app.do_one_mesh_per_phase
                number_phase = app.meshphase(1).number_phase;
                clear DATA_writetable
                for kphase = 1:1:number_phase
                    % Cell volume
                    [cell_volume, ~] = Function_CalculateCellVolume(app.meshphase(kphase).node_,app.meshphase(kphase).elem_); % Mesh volume information
                    % % Mesh quality
                    mesh_quality=meshquality(app.meshphase(kphase).node_,app.meshphase(kphase).elem_);
                    % Calculate cumulative and distribution function
                    [results, ~] = Function_probability_density(mesh_quality,[],parameters_smoothingmovingaverage);
                    parameters_distributionfigure.figurename = ['Tetrahedron quality of ' app.meshphase(kphase).name];
                    parameters_distributionfigure.filename = ['Tetrahedron_quality_distribution_' app.meshphase(kphase).name];
                    parameters_distributionfigure.subaxe1_title = ['Tetrahedron quality: cumulative fct., ' app.meshphase(kphase).name];
                    parameters_distributionfigure.subaxe2_title = ['Tetrahedron quality: distribution fct., ' app.meshphase(kphase).name];
                    parameters_distributionfigure.xlabel = {'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'};
                    data_distributionfigure.cumulative_fct = results.cumulative_fct;
                    data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
                    data_distributionfigure.probability_density_fct = results.probability_density_fct;
                    data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
                    data_distributionfigure.x50 = results.x50;
                    data_distributionfigure.smoothed_x50 = results.smoothed_x50;
                    data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
                    data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
                    data_distributionfigure.unit = '';
                    parameters_distributionfigure.title = ['Tetrahedron quality (Joe-Liu) of ' app.meshphase(kphase).name '. Minimum = ' num2str(min(mesh_quality),5)];
                    function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);
                    % Cell quality as function of volume
                    Fig_= figure; % Create figure
                    Fig_.Name= ['Tetrahedron quality with cell volume of ' app.meshphase(kphase).name];
                    Fig_.Color='white'; % Background colour
                    set(Fig_, 'Position', [460 428 660 520]);
                    axes_ = axes('Parent',Fig_); % Create axes
                    hold(axes_,'on');
                    t_=title (' ','FontName','Times New Roman','FontSize',16); % Title
                    t_.String= {'Tetrahedron quality with cell volume',app.meshphase(kphase).name};
                    h_=plot(cell_volume,mesh_quality); % Curves
                    set(h_,'LineStyle','none','Marker','x'); % Colors
                    xlabel('Cell volume (normalized with initial voxel cube size)'); % Axis label
                    ylabel({'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'});
                    grid(axes_,'on'); % Display grid
                    set(axes_,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                    set(axes_,'FontName','Times New Roman','FontSize',14); % - Fontname and fontsize
                    hold(axes_,'off');
                    % Save figures
                    if app.savefolder_is_defined
                        function_savefig(Fig_, parameters_distributionfigure.fullpath, ['Tetrahedron_quality_withcellvolume_' app.meshphase(kphase).name]);
                        DATA_writetable.sheet(kphase).name = app.meshphase(kphase).name;
                        DATA_writetable.sheet(kphase).table=table(cell_volume,mesh_quality,'VariableNames',{'Cell volume','Mesh quality'});
                    end
                end
                if app.savefolder_is_defined
                    Function_Writetable(parameters_distributionfigure.fullpath,'Tetrahedron_quality_withcellvolume_phases',DATA_writetable);
                end
            end
            if app.do_one_mesh_per_group
                number_group = app.meshgroup(1).number_group;
                clear DATA_writetable
                for kgroup = 1:1:number_group
                    % Cell volume
                    [cell_volume, ~] = Function_CalculateCellVolume(app.meshgroup(kgroup).node_,app.meshgroup(kgroup).elem_); % Mesh volume information
                    % % Mesh quality
                    mesh_quality=meshquality(app.meshgroup(kgroup).node_,app.meshgroup(kgroup).elem_);
                    % Calculate cumulative and distribution function
                    [results, ~] = Function_probability_density(mesh_quality,[],parameters_smoothingmovingaverage);
                    parameters_distributionfigure.figurename = ['Tetrahedron quality of ' app.meshgroup(kgroup).name];
                    parameters_distributionfigure.filename = ['Tetrahedron_quality_distribution_' app.meshgroup(kgroup).name];
                    parameters_distributionfigure.subaxe1_title = ['Tetrahedron quality: cumulative fct., ' app.meshgroup(kgroup).name];
                    parameters_distributionfigure.subaxe2_title = ['Tetrahedron quality: distribution fct., ' app.meshgroup(kgroup).name];
                    parameters_distributionfigure.xlabel = {'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'};
                    data_distributionfigure.cumulative_fct = results.cumulative_fct;
                    data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
                    data_distributionfigure.probability_density_fct = results.probability_density_fct;
                    data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
                    data_distributionfigure.x50 = results.x50;
                    data_distributionfigure.smoothed_x50 = results.smoothed_x50;
                    data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
                    data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
                    data_distributionfigure.unit = '';
                    parameters_distributionfigure.title = ['Tetrahedron quality (Joe-Liu) of ' app.meshgroup(kgroup).name '. Minimum = ' num2str(min(mesh_quality),5)];
                    function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);
                    % Cell quality as function of volume
                    Fig_= figure; % Create figure
                    Fig_.Name= ['Tetrahedron quality with cell volume of ' app.meshgroup(kgroup).name];
                    Fig_.Color='white'; % Background colour
                    set(Fig_, 'Position', [460 428 660 520]);
                    axes_ = axes('Parent',Fig_); % Create axes
                    hold(axes_,'on');
                    t_=title (' ','FontName','Times New Roman','FontSize',16); % Title
                    t_.String= {'Tetrahedron quality with cell volume',app.meshgroup(kgroup).name};
                    h_=plot(cell_volume,mesh_quality); % Curves
                    set(h_,'LineStyle','none','Marker','x'); % Colors
                    xlabel('Cell volume (normalized with initial voxel cube size)'); % Axis label
                    ylabel({'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'});
                    grid(axes_,'on'); % Display grid
                    set(axes_,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                    set(axes_,'FontName','Times New Roman','FontSize',14); % - Fontname and fontsize
                    hold(axes_,'off');
                    % Save figures
                    if app.savefolder_is_defined
                        function_savefig(Fig_, parameters_distributionfigure.fullpath, ['Tetrahedron_quality_withcellvolume_' app.meshgroup(kgroup).name]);
                        DATA_writetable.sheet(kgroup).name = app.meshgroup(kgroup).name;
                        DATA_writetable.sheet(kgroup).table=table(cell_volume(1:10),mesh_quality(1:10),'VariableNames',{'Cell volume','Mesh quality'});
                    end
                end
                if app.savefolder_is_defined
                    Function_Writetable(parameters_distributionfigure.fullpath,'Tetrahedron_quality_withcellvolume_group',DATA_writetable);
                end
            end
            if app.do_one_mesh_wholegeometry
                clear DATA_writetable
                % Cell volume
                [cell_volume, ~] = Function_CalculateCellVolume(app.meshvolume.node_,app.meshvolume.elem_); % Mesh volume information
                % % Mesh quality
                mesh_quality=meshquality(app.meshvolume.node_,app.meshvolume.elem_);
                % Calculate cumulative and distribution function
                [results, ~] = Function_probability_density(mesh_quality,[],parameters_smoothingmovingaverage);
                parameters_distributionfigure.figurename = ['Tetrahedron quality of ' app.meshvolume.name];
                parameters_distributionfigure.filename = ['Tetrahedron_quality_distribution_' app.meshvolume.name];
                parameters_distributionfigure.subaxe1_title = ['Tetrahedron quality: cumulative fct., ' app.meshvolume.name];
                parameters_distributionfigure.subaxe2_title = ['Tetrahedron quality: distribution fct., ' app.meshvolume.name];
                parameters_distributionfigure.xlabel = {'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'};
                data_distributionfigure.cumulative_fct = results.cumulative_fct;
                data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
                data_distributionfigure.probability_density_fct = results.probability_density_fct;
                data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
                data_distributionfigure.x50 = results.x50;
                data_distributionfigure.smoothed_x50 = results.smoothed_x50;
                data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
                data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
                data_distributionfigure.unit = '';
                parameters_distributionfigure.title = ['Tetrahedron quality (Joe-Liu) of ' app.meshvolume.name '. Minimum = ' num2str(min(mesh_quality),5)];
                function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);
                % Cell quality as function of volume
                Fig_= figure; % Create figure
                Fig_.Name= ['Tetrahedron quality with cell volume of ' app.meshvolume.name];
                Fig_.Color='white'; % Background colour
                set(Fig_, 'Position', [460 428 660 520]);
                axes_ = axes('Parent',Fig_); % Create axes
                hold(axes_,'on');
                t_=title (' ','FontName','Times New Roman','FontSize',16); % Title
                t_.String= {'Tetrahedron quality with cell volume',app.meshvolume.name};
                h_=plot(cell_volume,mesh_quality); % Curves
                set(h_,'LineStyle','none','Marker','x'); % Colors
                xlabel('Cell volume (normalized with initial voxel cube size)'); % Axis label
                ylabel({'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'});
                grid(axes_,'on'); % Display grid
                set(axes_,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                set(axes_,'FontName','Times New Roman','FontSize',14); % - Fontname and fontsize
                hold(axes_,'off');
                % Save figures
                if app.savefolder_is_defined
                    function_savefig(Fig_, parameters_distributionfigure.fullpath, ['Tetrahedron_quality_withcellvolume_' app.meshvolume.name]);
                    DATA_writetable.sheet(1).name = app.meshvolume.name;
                    DATA_writetable.sheet(1).table=table(cell_volume,mesh_quality,'VariableNames',{'Cell volume','Mesh quality'});
                end
            end
            if app.savefolder_is_defined
                Function_Writetable(parameters_distributionfigure.fullpath,['Tetrahedron_quality_withcellvolume_' app.meshvolume.name],DATA_writetable);
            end
        end

        % Value changed function: Surfacemesh_MethodDropDown
        function Surfacemesh_MethodDropDownValueChanged(app, event)
            choice = app.Surfacemesh_MethodDropDown.Value;
            str_iteration = 'Number of smoothing iteration';
            str_useralpha = 'Scaler, smoothing parameter, v(k+1)=(1-alpha)*v(k)+alpha*mean(neighbors)';
            str_userbeta = 'Scaler, smoothing parameter, for laplacianhc';
            if strcmp(choice,'lowpass') || strcmp(choice,'laplacian')
                data_ = [{'Iteration number';'useralpha'} {str_iteration;str_useralpha} num2cell([2;0.5])];
            elseif strcmp(choice,'laplacianhc') 
                data_ = [{'Iteration number';'useralpha';'userbeta'} {str_iteration;str_useralpha;str_userbeta} num2cell([2;0.5;1])];
            end
            app.Smoothing_UITable.Data=data_;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create MeshingmoduleUIFigure and hide until all components are created
            app.MeshingmoduleUIFigure = uifigure('Visible', 'off');
            app.MeshingmoduleUIFigure.Position = [100 100 1121 667];
            app.MeshingmoduleUIFigure.Name = 'Meshing module';
            app.MeshingmoduleUIFigure.Icon = 'Icon_Meshing.png';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.MeshingmoduleUIFigure);
            app.TabGroup.TabLocation = 'left';
            app.TabGroup.Position = [1 1 1121 667];

            % Create InstructionsTab
            app.InstructionsTab = uitab(app.TabGroup);
            app.InstructionsTab.Title = 'Instructions';

            % Create Instructions_title
            app.Instructions_title = uilabel(app.InstructionsTab);
            app.Instructions_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_title.HorizontalAlignment = 'center';
            app.Instructions_title.FontWeight = 'bold';
            app.Instructions_title.Position = [10 634 944 22];
            app.Instructions_title.Text = 'Main instructions';

            % Create Instructions_Label_3
            app.Instructions_Label_3 = uitextarea(app.InstructionsTab);
            app.Instructions_Label_3.Position = [11 287 944 280];
            app.Instructions_Label_3.Value = {'- Select your save folder in the next tab. Parameters, intermediate tif files, and meshes will be saved in this folder.'; ''; '- You can either choose a unique microstructure (i.e. loading a unique tif file), a half-cell or full cell (i.e. loading from 2 to 5 different tif files), or a polychristalline microstructure (i.e. loading a unique tif file with dozens-hundreds of different integer values) in the third tab. Note that in both cases, the tif files must represent segmented/binarize image(s). You will be able to apply some basic operations (e.g. cropping, scaling) to select a Region Of Interest. For more complex operations (e.g., rotations) please use before the "ROI, Filtering, and segmentation" module. '; 'For the half-cell or full cell cases, the different volumes may not share the same in-plane dimensions. If you did not crop them accordingly in the third tab, the next one "Dimension compatibility" will care of this. Note that this tab will be skipped if you have selected either a unique microstructure or a polychristalline microstructure.'; ''; '- Of high importance for the mesh generation to be successful and for the to-be-generated mesh to be high quality is the "Morphology opening" tab. Depending on the morphology complexity, the "erosion+dilation" step can be skipped (a complex geometry should use it typically). The "unique connected cluster" options depend on your model - if it requires unique connected cluster for each phase, or if it allows disconnected clusters. Laslty, the "clean voxel connections" is a must-go for most geometries (at the exception of simple test geometries for which all voxels are well-connected, such as periodically aligned spheres) as it will significanly improve the mesh quality. These operations will affect the morphology and volume fractions, although if the volume is nearly fully connected and enough refined then the modifications will be negligible since restricted only to the phase interfaces.'; ''; '- "The assemble cell" tab is where you will decide to assign each phase to a distinct id (and group id for segragated models).'; ''; 'At each step you will be able to visualize the microstructures to verify your choice.'};

            % Create Instructions_Label_1
            app.Instructions_Label_1 = uilabel(app.InstructionsTab);
            app.Instructions_Label_1.Position = [10 599 938 22];
            app.Instructions_Label_1.Text = 'The meshing module breaks down the microstructure volume meshing into three main steps: pre-processing, setting meshing options, and lastly the meshing creation itself.';

            % Create Instructions_Label_2
            app.Instructions_Label_2 = uilabel(app.InstructionsTab);
            app.Instructions_Label_2.FontWeight = 'bold';
            app.Instructions_Label_2.FontColor = [0 0 1];
            app.Instructions_Label_2.Position = [11 572 717 22];
            app.Instructions_Label_2.Text = '- Pre-processing tasks (blue tabs) are all the modifications applied to the 3D array in prevision of its subsequent meshing.';

            % Create Instructions_Label_4
            app.Instructions_Label_4 = uilabel(app.InstructionsTab);
            app.Instructions_Label_4.FontWeight = 'bold';
            app.Instructions_Label_4.FontColor = [0.851 0.3255 0.098];
            app.Instructions_Label_4.Position = [11 255 611 22];
            app.Instructions_Label_4.Text = '- Meshing options (orange tabs) allows you to select the meshing methods and their respective options.';

            % Create Instructions_Label_5
            app.Instructions_Label_5 = uitextarea(app.InstructionsTab);
            app.Instructions_Label_5.Position = [11 206 944 44];
            app.Instructions_Label_5.Value = {'- The "Mesh generation choice" let you decide between structured and unstructured meshes. You will be bring to the associated tab to set up the parameters. Note that if you select the Iso2mesh option, you will need to quote also Iso2mesh article in addition to this toolbox (see about tab for quotation instructions).'};

            % Create Instructions_Label_6
            app.Instructions_Label_6 = uilabel(app.InstructionsTab);
            app.Instructions_Label_6.FontWeight = 'bold';
            app.Instructions_Label_6.FontColor = [1 0 0];
            app.Instructions_Label_6.Position = [11 174 376 22];
            app.Instructions_Label_6.Text = '- The create mesh (red tab) is where you will generate the mesh.';

            % Create Instructions_Label_7
            app.Instructions_Label_7 = uitextarea(app.InstructionsTab);
            app.Instructions_Label_7.Position = [11 55 944 114];
            app.Instructions_Label_7.Value = {'- You can decide to generate a unique mesh containing the whole geometry (with cells assigned to an id that correspond to the associated phases) and/or meshes for each group (group being a union of different phases, for instance a group for the electrolyte and the separator and another group for the electrode and its current collector), and/or meshed for each phases. Id selected in the "assemble cell" tab are used to inform the cell assignements. '; ''; '- Mesh generation can be very time-expensive and RAM-expensive. It is recommended to test the module on a test small volume before running a large geometries.'; ''; 'Once the mesh generated, you can visualize it and save it in the same tab. In addition you can calculate mesh quality in the dedicated magenta tab.'};

            % Create FolderandsaveoptionsTab
            app.FolderandsaveoptionsTab = uitab(app.TabGroup);
            app.FolderandsaveoptionsTab.Title = 'Folder and save options';
            app.FolderandsaveoptionsTab.ForegroundColor = [0 0 1];

            % Create Saveoptions_title
            app.Saveoptions_title = uilabel(app.FolderandsaveoptionsTab);
            app.Saveoptions_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Saveoptions_title.HorizontalAlignment = 'center';
            app.Saveoptions_title.FontWeight = 'bold';
            app.Saveoptions_title.Position = [10 634 944 22];
            app.Saveoptions_title.Text = 'Folder and save options';

            % Create Save_Label
            app.Save_Label = uilabel(app.FolderandsaveoptionsTab);
            app.Save_Label.Position = [10 599 944 22];
            app.Save_Label.Text = 'Choose a save folder and then select the data in the table below you would like to save.';

            % Create ClicktoselectsavefolderButton
            app.ClicktoselectsavefolderButton = uibutton(app.FolderandsaveoptionsTab, 'push');
            app.ClicktoselectsavefolderButton.ButtonPushedFcn = createCallbackFcn(app, @ClicktoselectsavefolderButtonPushed, true);
            app.ClicktoselectsavefolderButton.BackgroundColor = [0.8 0.8 0.8];
            app.ClicktoselectsavefolderButton.Position = [10 567 166 26];
            app.ClicktoselectsavefolderButton.Text = 'Click to select save folder';

            % Create Save_folder_text
            app.Save_folder_text = uilabel(app.FolderandsaveoptionsTab);
            app.Save_folder_text.FontColor = [1 0 0];
            app.Save_folder_text.Position = [186 569 768 22];
            app.Save_folder_text.Text = 'Save folder location: NOT DEFINED';

            % Create TabisNOTsetupcorrectlyLabel
            app.TabisNOTsetupcorrectlyLabel = uilabel(app.FolderandsaveoptionsTab);
            app.TabisNOTsetupcorrectlyLabel.HorizontalAlignment = 'right';
            app.TabisNOTsetupcorrectlyLabel.FontSize = 14;
            app.TabisNOTsetupcorrectlyLabel.FontWeight = 'bold';
            app.TabisNOTsetupcorrectlyLabel.Position = [703 25 183 22];
            app.TabisNOTsetupcorrectlyLabel.Text = 'Tab is NOT setup correctly';

            % Create Folder_Lamp
            app.Folder_Lamp = uilamp(app.FolderandsaveoptionsTab);
            app.Folder_Lamp.Position = [901 9 53 53];
            app.Folder_Lamp.Color = [1 0 0];

            % Create Setup_UITable
            app.Setup_UITable = uitable(app.FolderandsaveoptionsTab);
            app.Setup_UITable.ColumnName = {'Parameters and tif'; 'Choice'};
            app.Setup_UITable.ColumnWidth = {'1x'};
            app.Setup_UITable.RowName = {};
            app.Setup_UITable.ColumnEditable = [false true];
            app.Setup_UITable.CellEditCallback = createCallbackFcn(app, @Setup_UITableCellEdit, true);
            app.Setup_UITable.Position = [10 445 473 105];

            % Create SelectmicrostructuresTab
            app.SelectmicrostructuresTab = uitab(app.TabGroup);
            app.SelectmicrostructuresTab.Title = 'Select microstructures';
            app.SelectmicrostructuresTab.ForegroundColor = [0 0 1];

            % Create SelectMicrostructure_Title
            app.SelectMicrostructure_Title = uilabel(app.SelectmicrostructuresTab);
            app.SelectMicrostructure_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.SelectMicrostructure_Title.HorizontalAlignment = 'center';
            app.SelectMicrostructure_Title.FontWeight = 'bold';
            app.SelectMicrostructure_Title.Position = [10 634 944 22];
            app.SelectMicrostructure_Title.Text = 'Microstructure choice (One unique microstructure OR Half-cell or full-cell OR Polycrystalline architecture). Import segmented volumes only.';

            % Create OneuniquemicrostructureLabel
            app.OneuniquemicrostructureLabel = uilabel(app.SelectmicrostructuresTab);
            app.OneuniquemicrostructureLabel.BackgroundColor = [0.902 0.902 0.902];
            app.OneuniquemicrostructureLabel.HorizontalAlignment = 'center';
            app.OneuniquemicrostructureLabel.FontWeight = 'bold';
            app.OneuniquemicrostructureLabel.Position = [10 575 288 45];
            app.OneuniquemicrostructureLabel.Text = 'One unique microstructure';

            % Create Halfcellorfullcellselectatleast2volumesLabel
            app.Halfcellorfullcellselectatleast2volumesLabel = uilabel(app.SelectmicrostructuresTab);
            app.Halfcellorfullcellselectatleast2volumesLabel.BackgroundColor = [0.902 0.902 0.902];
            app.Halfcellorfullcellselectatleast2volumesLabel.HorizontalAlignment = 'center';
            app.Halfcellorfullcellselectatleast2volumesLabel.FontWeight = 'bold';
            app.Halfcellorfullcellselectatleast2volumesLabel.Position = [325 598 315 22];
            app.Halfcellorfullcellselectatleast2volumesLabel.Text = 'Half-cell or full-cell (select at least 2 volumes)';

            % Create PolycrystallinearchitectureBetaLabel
            app.PolycrystallinearchitectureBetaLabel = uilabel(app.SelectmicrostructuresTab);
            app.PolycrystallinearchitectureBetaLabel.BackgroundColor = [0.902 0.902 0.902];
            app.PolycrystallinearchitectureBetaLabel.HorizontalAlignment = 'center';
            app.PolycrystallinearchitectureBetaLabel.FontWeight = 'bold';
            app.PolycrystallinearchitectureBetaLabel.Position = [666 575 288 45];
            app.PolycrystallinearchitectureBetaLabel.Text = 'Polycrystalline architecture (Beta)';

            % Create ImporttifButton_Unique
            app.ImporttifButton_Unique = uibutton(app.SelectmicrostructuresTab, 'push');
            app.ImporttifButton_Unique.ButtonPushedFcn = createCallbackFcn(app, @ImporttifButton_UniquePushed, true);
            app.ImporttifButton_Unique.BackgroundColor = [0.8 0.8 0.8];
            app.ImporttifButton_Unique.Position = [10 514 288 51];
            app.ImporttifButton_Unique.Text = 'Import .tif';

            % Create ImporttifButton_Cell
            app.ImporttifButton_Cell = uibutton(app.SelectmicrostructuresTab, 'push');
            app.ImporttifButton_Cell.ButtonPushedFcn = createCallbackFcn(app, @ImporttifButton_CellPushed, true);
            app.ImporttifButton_Cell.BackgroundColor = [0.8 0.8 0.8];
            app.ImporttifButton_Cell.Enable = 'off';
            app.ImporttifButton_Cell.Position = [338 514 106 22];
            app.ImporttifButton_Cell.Text = 'Import .tif';

            % Create ImporttifButton_Poly
            app.ImporttifButton_Poly = uibutton(app.SelectmicrostructuresTab, 'push');
            app.ImporttifButton_Poly.ButtonPushedFcn = createCallbackFcn(app, @ImporttifButton_PolyPushed, true);
            app.ImporttifButton_Poly.BackgroundColor = [0.8 0.8 0.8];
            app.ImporttifButton_Poly.Position = [666 514 288 51];
            app.ImporttifButton_Poly.Text = 'Import .tif';

            % Create LeftelectrodeLabel
            app.LeftelectrodeLabel = uilabel(app.SelectmicrostructuresTab);
            app.LeftelectrodeLabel.BackgroundColor = [0.902 0.902 0.902];
            app.LeftelectrodeLabel.HorizontalAlignment = 'center';
            app.LeftelectrodeLabel.Position = [354 575 90 24];
            app.LeftelectrodeLabel.Text = 'Left electrode';

            % Create SeparatorLabel
            app.SeparatorLabel = uilabel(app.SelectmicrostructuresTab);
            app.SeparatorLabel.BackgroundColor = [0.902 0.902 0.902];
            app.SeparatorLabel.HorizontalAlignment = 'center';
            app.SeparatorLabel.Position = [449 575 67 24];
            app.SeparatorLabel.Text = 'Separator';

            % Create RightelectrodeLabel
            app.RightelectrodeLabel = uilabel(app.SelectmicrostructuresTab);
            app.RightelectrodeLabel.BackgroundColor = [0.902 0.902 0.902];
            app.RightelectrodeLabel.HorizontalAlignment = 'center';
            app.RightelectrodeLabel.Position = [521 575 90 24];
            app.RightelectrodeLabel.Text = 'Right electrode';

            % Create ChoosedomainDropDownLabel
            app.ChoosedomainDropDownLabel = uilabel(app.SelectmicrostructuresTab);
            app.ChoosedomainDropDownLabel.HorizontalAlignment = 'right';
            app.ChoosedomainDropDownLabel.Position = [338 543 90 22];
            app.ChoosedomainDropDownLabel.Text = 'Choose domain';

            % Create ChoosedomainDropDown
            app.ChoosedomainDropDown = uidropdown(app.SelectmicrostructuresTab);
            app.ChoosedomainDropDown.Items = {'', 'Left current collector', 'Left electrode', 'Separator', 'Right electrode', 'Right current collector'};
            app.ChoosedomainDropDown.ValueChangedFcn = createCallbackFcn(app, @ChoosedomainDropDownValueChanged, true);
            app.ChoosedomainDropDown.Position = [443 543 183 22];
            app.ChoosedomainDropDown.Value = '';

            % Create HomogenousButton_Cell
            app.HomogenousButton_Cell = uibutton(app.SelectmicrostructuresTab, 'push');
            app.HomogenousButton_Cell.ButtonPushedFcn = createCallbackFcn(app, @HomogenousButton_CellPushed, true);
            app.HomogenousButton_Cell.BackgroundColor = [0.8 0.8 0.8];
            app.HomogenousButton_Cell.Enable = 'off';
            app.HomogenousButton_Cell.Position = [482 514 144 22];
            app.HomogenousButton_Cell.Text = 'Homogenous medium';

            % Create orLabel
            app.orLabel = uilabel(app.SelectmicrostructuresTab);
            app.orLabel.HorizontalAlignment = 'center';
            app.orLabel.Position = [451 514 25 22];
            app.orLabel.Text = 'or';

            % Create LengthinvoxelEditFieldLabel
            app.LengthinvoxelEditFieldLabel = uilabel(app.SelectmicrostructuresTab);
            app.LengthinvoxelEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthinvoxelEditFieldLabel.Position = [479 488 86 22];
            app.LengthinvoxelEditFieldLabel.Text = 'Length in voxel';

            % Create LengthinvoxelEditField
            app.LengthinvoxelEditField = uieditfield(app.SelectmicrostructuresTab, 'numeric');
            app.LengthinvoxelEditField.Limits = [1 Inf];
            app.LengthinvoxelEditField.RoundFractionalValues = 'on';
            app.LengthinvoxelEditField.ValueChangedFcn = createCallbackFcn(app, @LengthinvoxelEditFieldValueChanged, true);
            app.LengthinvoxelEditField.Enable = 'off';
            app.LengthinvoxelEditField.Position = [580 488 46 22];
            app.LengthinvoxelEditField.Value = 5;

            % Create SavemicrostructureButton
            app.SavemicrostructureButton = uibutton(app.SelectmicrostructuresTab, 'push');
            app.SavemicrostructureButton.ButtonPushedFcn = createCallbackFcn(app, @SavemicrostructureButtonPushed, true);
            app.SavemicrostructureButton.BackgroundColor = [0.8 0.8 0.8];
            app.SavemicrostructureButton.Enable = 'off';
            app.SavemicrostructureButton.Position = [460 9 225 22];
            app.SavemicrostructureButton.Text = 'Save microstructure';

            % Create Label
            app.Label = uilabel(app.SelectmicrostructuresTab);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontAngle = 'italic';
            app.Label.Position = [11 9 452 51];
            app.Label.Text = {'Once finished, you must click on ''Save microstructure''.'; 'You must repeat the import process at least 2 times for the half-cell of full-cell case'};

            % Create TiffileloadedLabel
            app.TiffileloadedLabel = uilabel(app.SelectmicrostructuresTab);
            app.TiffileloadedLabel.FontAngle = 'italic';
            app.TiffileloadedLabel.Position = [10 457 944 22];
            app.TiffileloadedLabel.Text = 'Tif file loaded';

            % Create Microstructure_vf_UITable
            app.Microstructure_vf_UITable = uitable(app.SelectmicrostructuresTab);
            app.Microstructure_vf_UITable.ColumnName = {'Id'; 'Initial'; 'ROI'; 'Scaling'};
            app.Microstructure_vf_UITable.ColumnWidth = {50, 'auto', 'auto', 'auto'};
            app.Microstructure_vf_UITable.RowName = {};
            app.Microstructure_vf_UITable.Position = [11 72 288 355];

            % Create VolumefractionsLabel
            app.VolumefractionsLabel = uilabel(app.SelectmicrostructuresTab);
            app.VolumefractionsLabel.Position = [12 427 285 22];
            app.VolumefractionsLabel.Text = 'Volume fractions';

            % Create Select_VisualizemicrostructureButton
            app.Select_VisualizemicrostructureButton = uibutton(app.SelectmicrostructuresTab, 'push');
            app.Select_VisualizemicrostructureButton.ButtonPushedFcn = createCallbackFcn(app, @Select_VisualizemicrostructureButtonPushed, true);
            app.Select_VisualizemicrostructureButton.BackgroundColor = [0.8 0.8 0.8];
            app.Select_VisualizemicrostructureButton.Enable = 'off';
            app.Select_VisualizemicrostructureButton.Position = [460 38 225 22];
            app.Select_VisualizemicrostructureButton.Text = 'Visualize microstructure';

            % Create Select_Lamp
            app.Select_Lamp = uilamp(app.SelectmicrostructuresTab);
            app.Select_Lamp.Position = [901 9 53 53];
            app.Select_Lamp.Color = [1 0 0];

            % Create TabisNOTsetupcorrectlyLabel_2
            app.TabisNOTsetupcorrectlyLabel_2 = uilabel(app.SelectmicrostructuresTab);
            app.TabisNOTsetupcorrectlyLabel_2.HorizontalAlignment = 'right';
            app.TabisNOTsetupcorrectlyLabel_2.FontSize = 14;
            app.TabisNOTsetupcorrectlyLabel_2.FontWeight = 'bold';
            app.TabisNOTsetupcorrectlyLabel_2.Position = [703 25 183 22];
            app.TabisNOTsetupcorrectlyLabel_2.Text = 'Tab is NOT setup correctly';

            % Create Select_ROI_UITable
            app.Select_ROI_UITable = uitable(app.SelectmicrostructuresTab);
            app.Select_ROI_UITable.ColumnName = {'Axis'; 'Length (before)'; 'Start'; 'End'; 'Length (after)'};
            app.Select_ROI_UITable.ColumnWidth = {50, 'auto', 'auto', 'auto', 'auto'};
            app.Select_ROI_UITable.RowName = {};
            app.Select_ROI_UITable.ColumnEditable = [false false true true false];
            app.Select_ROI_UITable.CellEditCallback = createCallbackFcn(app, @Select_ROI_UITableCellEdit, true);
            app.Select_ROI_UITable.Tooltip = {'If you need to perform a more advanced ROI selection (such as rotation), please use first the filtering and segmentation module.'};
            app.Select_ROI_UITable.Position = [339 329 615 98];

            % Create SelecttheregionofinterestROILabel
            app.SelecttheregionofinterestROILabel = uilabel(app.SelectmicrostructuresTab);
            app.SelecttheregionofinterestROILabel.FontWeight = 'bold';
            app.SelecttheregionofinterestROILabel.Position = [340 427 614 22];
            app.SelecttheregionofinterestROILabel.Text = '1) Select the region of interest (ROI). Lengths are expressed in number of voxels';

            % Create SelecttheregionofinterestROIoptionalLabel_2
            app.SelecttheregionofinterestROIoptionalLabel_2 = uilabel(app.SelectmicrostructuresTab);
            app.SelecttheregionofinterestROIoptionalLabel_2.FontWeight = 'bold';
            app.SelecttheregionofinterestROIoptionalLabel_2.Tooltip = {'Slices along the through-plane are ordered so that first slices are left-side, while last slices are right-side.'};
            app.SelecttheregionofinterestROIoptionalLabel_2.Position = [340 297 608 22];
            app.SelecttheregionofinterestROIoptionalLabel_2.Text = '2) Set orientation so that for half-cell or full-cell, the thickness (through-plane direction) is the first axis.';

            % Create Select_Orientation_UITable
            app.Select_Orientation_UITable = uitable(app.SelectmicrostructuresTab);
            app.Select_Orientation_UITable.ColumnName = {'Orientation'; 'Axis'; 'Length'};
            app.Select_Orientation_UITable.RowName = {};
            app.Select_Orientation_UITable.Position = [341 199 302 98];

            % Create Select_Rescale_UITable
            app.Select_Rescale_UITable = uitable(app.SelectmicrostructuresTab);
            app.Select_Rescale_UITable.ColumnName = {'Orientation'; 'Length (before)'; 'Length (after)'};
            app.Select_Rescale_UITable.RowName = {};
            app.Select_Rescale_UITable.Position = [338 69 347 98];

            % Create SelecttheregionofinterestROIoptionalLabel_3
            app.SelecttheregionofinterestROIoptionalLabel_3 = uilabel(app.SelectmicrostructuresTab);
            app.SelecttheregionofinterestROIoptionalLabel_3.FontWeight = 'bold';
            app.SelecttheregionofinterestROIoptionalLabel_3.Position = [340 167 587 22];
            app.SelecttheregionofinterestROIoptionalLabel_3.Text = '3) Apply up/down scaling. For a full-cell, you must rescale so that each tif share the same voxel size';

            % Create Rows12unswapedButton
            app.Rows12unswapedButton = uibutton(app.SelectmicrostructuresTab, 'state');
            app.Rows12unswapedButton.ValueChangedFcn = createCallbackFcn(app, @Rows12unswapedButtonValueChanged, true);
            app.Rows12unswapedButton.Text = 'Rows 1-2 (un-swaped)';
            app.Rows12unswapedButton.Position = [656 253 136 22];

            % Create Rows13unswapedButton
            app.Rows13unswapedButton = uibutton(app.SelectmicrostructuresTab, 'state');
            app.Rows13unswapedButton.ValueChangedFcn = createCallbackFcn(app, @Rows13unswapedButtonValueChanged, true);
            app.Rows13unswapedButton.Text = 'Rows 1-3 (un-swaped)';
            app.Rows13unswapedButton.Position = [656 226 136 22];

            % Create Rows23unswapedButton
            app.Rows23unswapedButton = uibutton(app.SelectmicrostructuresTab, 'state');
            app.Rows23unswapedButton.ValueChangedFcn = createCallbackFcn(app, @Rows23unswapedButtonValueChanged, true);
            app.Rows23unswapedButton.Text = 'Rows 2-3 (un-swaped)';
            app.Rows23unswapedButton.Position = [656 199 136 22];

            % Create Row1unflippedButton
            app.Row1unflippedButton = uibutton(app.SelectmicrostructuresTab, 'state');
            app.Row1unflippedButton.ValueChangedFcn = createCallbackFcn(app, @Row1unflippedButtonValueChanged, true);
            app.Row1unflippedButton.Text = 'Row 1 (un-flipped)';
            app.Row1unflippedButton.Position = [807 253 114 22];

            % Create Row2unflippedButton
            app.Row2unflippedButton = uibutton(app.SelectmicrostructuresTab, 'state');
            app.Row2unflippedButton.ValueChangedFcn = createCallbackFcn(app, @Row2unflippedButtonValueChanged, true);
            app.Row2unflippedButton.Text = 'Row 2 (un-flipped)';
            app.Row2unflippedButton.Position = [807 226 114 22];

            % Create Row3unflippedButton
            app.Row3unflippedButton = uibutton(app.SelectmicrostructuresTab, 'state');
            app.Row3unflippedButton.ValueChangedFcn = createCallbackFcn(app, @Row3unflippedButtonValueChanged, true);
            app.Row3unflippedButton.Text = 'Row 3 (un-flipped)';
            app.Row3unflippedButton.Position = [807 199 114 22];

            % Create ScalingfactorEditFieldLabel
            app.ScalingfactorEditFieldLabel = uilabel(app.SelectmicrostructuresTab);
            app.ScalingfactorEditFieldLabel.HorizontalAlignment = 'right';
            app.ScalingfactorEditFieldLabel.Position = [703 145 78 22];
            app.ScalingfactorEditFieldLabel.Text = 'Scaling factor';

            % Create ScalingfactorEditField
            app.ScalingfactorEditField = uieditfield(app.SelectmicrostructuresTab, 'numeric');
            app.ScalingfactorEditField.Limits = [1e-06 Inf];
            app.ScalingfactorEditField.ValueChangedFcn = createCallbackFcn(app, @ScalingfactorEditFieldValueChanged, true);
            app.ScalingfactorEditField.Position = [823 145 42 22];
            app.ScalingfactorEditField.Value = 1;

            % Create VoxelsizebeforeEditFieldLabel
            app.VoxelsizebeforeEditFieldLabel = uilabel(app.SelectmicrostructuresTab);
            app.VoxelsizebeforeEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelsizebeforeEditFieldLabel.Position = [703 121 105 22];
            app.VoxelsizebeforeEditFieldLabel.Text = 'Voxel size (before)';

            % Create VoxelsizebeforeEditField
            app.VoxelsizebeforeEditField = uieditfield(app.SelectmicrostructuresTab, 'numeric');
            app.VoxelsizebeforeEditField.Limits = [0 Inf];
            app.VoxelsizebeforeEditField.ValueChangedFcn = createCallbackFcn(app, @VoxelsizebeforeEditFieldValueChanged, true);
            app.VoxelsizebeforeEditField.Position = [823 118 42 22];
            app.VoxelsizebeforeEditField.Value = 1;

            % Create VoxelsizeafterEditFieldLabel
            app.VoxelsizeafterEditFieldLabel = uilabel(app.SelectmicrostructuresTab);
            app.VoxelsizeafterEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelsizeafterEditFieldLabel.Position = [703 94 95 22];
            app.VoxelsizeafterEditFieldLabel.Text = 'Voxel size (after)';

            % Create VoxelsizeafterEditField
            app.VoxelsizeafterEditField = uieditfield(app.SelectmicrostructuresTab, 'numeric');
            app.VoxelsizeafterEditField.Limits = [0 Inf];
            app.VoxelsizeafterEditField.Editable = 'off';
            app.VoxelsizeafterEditField.Position = [823 94 42 22];
            app.VoxelsizeafterEditField.Value = 1;

            % Create ApplyscalingButton
            app.ApplyscalingButton = uibutton(app.SelectmicrostructuresTab, 'push');
            app.ApplyscalingButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyscalingButtonPushed, true);
            app.ApplyscalingButton.BackgroundColor = [0.8 0.8 0.8];
            app.ApplyscalingButton.Position = [879 72 55 95];
            app.ApplyscalingButton.Text = {'Apply'; 'scaling'};

            % Create BackgroundidEditFieldLabel
            app.BackgroundidEditFieldLabel = uilabel(app.SelectmicrostructuresTab);
            app.BackgroundidEditFieldLabel.HorizontalAlignment = 'right';
            app.BackgroundidEditFieldLabel.Position = [703 70 82 22];
            app.BackgroundidEditFieldLabel.Text = 'Background id';

            % Create BackgroundidEditField
            app.BackgroundidEditField = uieditfield(app.SelectmicrostructuresTab, 'numeric');
            app.BackgroundidEditField.Limits = [0 Inf];
            app.BackgroundidEditField.RoundFractionalValues = 'on';
            app.BackgroundidEditField.Position = [823 70 42 22];

            % Create CC_left
            app.CC_left = uilabel(app.SelectmicrostructuresTab);
            app.CC_left.BackgroundColor = [0.902 0.902 0.902];
            app.CC_left.HorizontalAlignment = 'center';
            app.CC_left.Position = [325 575 24 24];
            app.CC_left.Text = 'CC';

            % Create CC_right
            app.CC_right = uilabel(app.SelectmicrostructuresTab);
            app.CC_right.BackgroundColor = [0.902 0.902 0.902];
            app.CC_right.HorizontalAlignment = 'center';
            app.CC_right.Position = [616 575 24 24];
            app.CC_right.Text = 'CC';

            % Create DimensioncompatibilityTab
            app.DimensioncompatibilityTab = uitab(app.TabGroup);
            app.DimensioncompatibilityTab.Title = 'Dimension compatibility';
            app.DimensioncompatibilityTab.ForegroundColor = [0 0 1];

            % Create Compatibility_Title
            app.Compatibility_Title = uilabel(app.DimensioncompatibilityTab);
            app.Compatibility_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Compatibility_Title.HorizontalAlignment = 'center';
            app.Compatibility_Title.FontWeight = 'bold';
            app.Compatibility_Title.Position = [10 634 944 22];
            app.Compatibility_Title.Text = 'Dimensions. Only required for half-cell and full-cell case';

            % Create Comp_voxelsize_text
            app.Comp_voxelsize_text = uilabel(app.DimensioncompatibilityTab);
            app.Comp_voxelsize_text.FontWeight = 'bold';
            app.Comp_voxelsize_text.Position = [10 534 902 22];
            app.Comp_voxelsize_text.Text = '2) Enter voxel size. At this stage all volumes must share the same voxel size (for half-cell and full-cell case). If not, you must rescale it in the previous tab.';

            % Create Comp_Lamp
            app.Comp_Lamp = uilamp(app.DimensioncompatibilityTab);
            app.Comp_Lamp.Position = [901 9 53 53];
            app.Comp_Lamp.Color = [1 0 0];

            % Create TabisNOTsetupcorrectlyLabel_3
            app.TabisNOTsetupcorrectlyLabel_3 = uilabel(app.DimensioncompatibilityTab);
            app.TabisNOTsetupcorrectlyLabel_3.HorizontalAlignment = 'right';
            app.TabisNOTsetupcorrectlyLabel_3.FontSize = 14;
            app.TabisNOTsetupcorrectlyLabel_3.FontWeight = 'bold';
            app.TabisNOTsetupcorrectlyLabel_3.Position = [703 25 183 22];
            app.TabisNOTsetupcorrectlyLabel_3.Text = 'Tab is NOT setup correctly';

            % Create Label_2
            app.Label_2 = uilabel(app.DimensioncompatibilityTab);
            app.Label_2.FontAngle = 'italic';
            app.Label_2.Position = [10 512 655 22];
            app.Label_2.Text = 'Note that voxel size is optional (meshing is adimensional). It is here to help you get the actual size of the microstructure.';

            % Create Comp_voxelsize_text_2
            app.Comp_voxelsize_text_2 = uilabel(app.DimensioncompatibilityTab);
            app.Comp_voxelsize_text_2.Position = [11 210 418 22];
            app.Comp_voxelsize_text_2.Text = 'Along thickness (no action required)';

            % Create Dim_TP_UITable
            app.Dim_TP_UITable = uitable(app.DimensioncompatibilityTab);
            app.Dim_TP_UITable.ColumnName = {'Volume'; 'Voxel'; 'Length'};
            app.Dim_TP_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.Dim_TP_UITable.RowName = {};
            app.Dim_TP_UITable.Position = [11 63 418 147];

            % Create Comp_voxelsize_text_3
            app.Comp_voxelsize_text_3 = uilabel(app.DimensioncompatibilityTab);
            app.Comp_voxelsize_text_3.FontWeight = 'bold';
            app.Comp_voxelsize_text_3.Position = [10 448 538 22];
            app.Comp_voxelsize_text_3.Text = '3) All volumes must share the same in-plane dimensions to be able to assemble the full cell';

            % Create VoxelsizeEditFieldLabel
            app.VoxelsizeEditFieldLabel = uilabel(app.DimensioncompatibilityTab);
            app.VoxelsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelsizeEditFieldLabel.Position = [10 485 60 22];
            app.VoxelsizeEditFieldLabel.Text = 'Voxel size';

            % Create VoxelsizeEditField
            app.VoxelsizeEditField = uieditfield(app.DimensioncompatibilityTab, 'numeric');
            app.VoxelsizeEditField.Limits = [0 Inf];
            app.VoxelsizeEditField.ValueChangedFcn = createCallbackFcn(app, @VoxelsizeEditFieldValueChanged, true);
            app.VoxelsizeEditField.Position = [85 485 62 22];
            app.VoxelsizeEditField.Value = 1;

            % Create ActualizeresetdimensionsButton
            app.ActualizeresetdimensionsButton = uibutton(app.DimensioncompatibilityTab, 'push');
            app.ActualizeresetdimensionsButton.ButtonPushedFcn = createCallbackFcn(app, @ActualizeresetdimensionsButtonPushed, true);
            app.ActualizeresetdimensionsButton.BackgroundColor = [0.8 0.8 0.8];
            app.ActualizeresetdimensionsButton.Enable = 'off';
            app.ActualizeresetdimensionsButton.Position = [10 571 158 22];
            app.ActualizeresetdimensionsButton.Text = 'Actualize/reset dimensions';

            % Create Actualizedim_text
            app.Actualizedim_text = uilabel(app.DimensioncompatibilityTab);
            app.Actualizedim_text.FontWeight = 'bold';
            app.Actualizedim_text.Position = [10 598 520 22];
            app.Actualizedim_text.Text = '1) Once previous tab finished, actualize dimensions and check dimensions compatibility';

            % Create Comp_voxelsize_text_4
            app.Comp_voxelsize_text_4 = uilabel(app.DimensioncompatibilityTab);
            app.Comp_voxelsize_text_4.Position = [11 421 441 22];
            app.Comp_voxelsize_text_4.Text = 'Along in-plane direction 1 (dimensions must match)';

            % Create Dim_IP1_UITable
            app.Dim_IP1_UITable = uitable(app.DimensioncompatibilityTab);
            app.Dim_IP1_UITable.ColumnName = {'Volume'; 'Voxel'; 'Start'; 'End'; 'Length'};
            app.Dim_IP1_UITable.ColumnWidth = {'1x', 'auto', 'auto', 'auto', 'auto'};
            app.Dim_IP1_UITable.RowName = {};
            app.Dim_IP1_UITable.ColumnEditable = [false false true true false];
            app.Dim_IP1_UITable.CellEditCallback = createCallbackFcn(app, @Dim_IP1_UITableCellEdit, true);
            app.Dim_IP1_UITable.Position = [11 274 441 147];

            % Create Comp_voxelsize_text_5
            app.Comp_voxelsize_text_5 = uilabel(app.DimensioncompatibilityTab);
            app.Comp_voxelsize_text_5.Position = [513 421 441 22];
            app.Comp_voxelsize_text_5.Text = 'Along in-plane direction 2 (dimensions must match)';

            % Create Dim_IP2_UITable
            app.Dim_IP2_UITable = uitable(app.DimensioncompatibilityTab);
            app.Dim_IP2_UITable.ColumnName = {'Volume'; 'Voxel'; 'Start'; 'End'; 'Length'};
            app.Dim_IP2_UITable.ColumnWidth = {'1x', 'auto', 'auto', 'auto', 'auto'};
            app.Dim_IP2_UITable.RowName = {};
            app.Dim_IP2_UITable.ColumnEditable = [false false true true false];
            app.Dim_IP2_UITable.CellEditCallback = createCallbackFcn(app, @Dim_IP2_UITableCellEdit, true);
            app.Dim_IP2_UITable.Position = [513 274 441 147];

            % Create Dim_Autocrop_1_Button
            app.Dim_Autocrop_1_Button = uibutton(app.DimensioncompatibilityTab, 'push');
            app.Dim_Autocrop_1_Button.ButtonPushedFcn = createCallbackFcn(app, @Dim_Autocrop_1_ButtonPushed, true);
            app.Dim_Autocrop_1_Button.BackgroundColor = [0.8 0.8 0.8];
            app.Dim_Autocrop_1_Button.Position = [11 247 100 22];
            app.Dim_Autocrop_1_Button.Text = 'Autocrop';

            % Create Dim_Autocrop_2_Button
            app.Dim_Autocrop_2_Button = uibutton(app.DimensioncompatibilityTab, 'push');
            app.Dim_Autocrop_2_Button.ButtonPushedFcn = createCallbackFcn(app, @Dim_Autocrop_2_ButtonPushed, true);
            app.Dim_Autocrop_2_Button.BackgroundColor = [0.8 0.8 0.8];
            app.Dim_Autocrop_2_Button.Position = [513 247 100 22];
            app.Dim_Autocrop_2_Button.Text = 'Autocrop';

            % Create Label_3
            app.Label_3 = uilabel(app.DimensioncompatibilityTab);
            app.Label_3.HorizontalAlignment = 'center';
            app.Label_3.FontAngle = 'italic';
            app.Label_3.Position = [11 9 452 51];
            app.Label_3.Text = 'When in-plane dimensions are identical for all volumes, save the cell to move on';

            % Create SavecellButton
            app.SavecellButton = uibutton(app.DimensioncompatibilityTab, 'push');
            app.SavecellButton.ButtonPushedFcn = createCallbackFcn(app, @SavecellButtonPushed, true);
            app.SavecellButton.BackgroundColor = [0.8 0.8 0.8];
            app.SavecellButton.Enable = 'off';
            app.SavecellButton.Position = [460 9 225 22];
            app.SavecellButton.Text = 'Save cell';

            % Create Select_VisualizecellButton
            app.Select_VisualizecellButton = uibutton(app.DimensioncompatibilityTab, 'push');
            app.Select_VisualizecellButton.ButtonPushedFcn = createCallbackFcn(app, @Select_VisualizecellButtonPushed, true);
            app.Select_VisualizecellButton.BackgroundColor = [0.8 0.8 0.8];
            app.Select_VisualizecellButton.Enable = 'off';
            app.Select_VisualizecellButton.Position = [460 38 225 22];
            app.Select_VisualizecellButton.Text = 'Visualize cell';

            % Create Inplanedimension1arematchingLabel
            app.Inplanedimension1arematchingLabel = uilabel(app.DimensioncompatibilityTab);
            app.Inplanedimension1arematchingLabel.FontWeight = 'bold';
            app.Inplanedimension1arematchingLabel.Position = [121 247 331 22];
            app.Inplanedimension1arematchingLabel.Text = 'In-plane dimension 1 are matching';

            % Create Inplanedimension2arematchingLabel
            app.Inplanedimension2arematchingLabel = uilabel(app.DimensioncompatibilityTab);
            app.Inplanedimension2arematchingLabel.FontWeight = 'bold';
            app.Inplanedimension2arematchingLabel.Position = [623 247 331 22];
            app.Inplanedimension2arematchingLabel.Text = 'In-plane dimension 2 are matching';

            % Create Dim_YoucanskipthistabLabel
            app.Dim_YoucanskipthistabLabel = uilabel(app.DimensioncompatibilityTab);
            app.Dim_YoucanskipthistabLabel.BackgroundColor = [0.302 0.7451 0.9333];
            app.Dim_YoucanskipthistabLabel.HorizontalAlignment = 'center';
            app.Dim_YoucanskipthistabLabel.FontSize = 18;
            app.Dim_YoucanskipthistabLabel.FontWeight = 'bold';
            app.Dim_YoucanskipthistabLabel.Visible = 'off';
            app.Dim_YoucanskipthistabLabel.Position = [503 79 451 53];
            app.Dim_YoucanskipthistabLabel.Text = 'Only one domain: you can ignore and skip this tab.';

            % Create micrometersLabel
            app.micrometersLabel = uilabel(app.DimensioncompatibilityTab);
            app.micrometersLabel.FontAngle = 'italic';
            app.micrometersLabel.Position = [151 485 72 22];
            app.micrometersLabel.Text = 'micrometers';

            % Create MorphologyopeningTab
            app.MorphologyopeningTab = uitab(app.TabGroup);
            app.MorphologyopeningTab.Title = 'Morphology opening';
            app.MorphologyopeningTab.ForegroundColor = [0 0 1];

            % Create PerformalloperationssequentiallyLabel
            app.PerformalloperationssequentiallyLabel = uilabel(app.MorphologyopeningTab);
            app.PerformalloperationssequentiallyLabel.BackgroundColor = [1 1 1];
            app.PerformalloperationssequentiallyLabel.HorizontalAlignment = 'center';
            app.PerformalloperationssequentiallyLabel.VerticalAlignment = 'top';
            app.PerformalloperationssequentiallyLabel.FontWeight = 'bold';
            app.PerformalloperationssequentiallyLabel.Position = [721 54 233 135];
            app.PerformalloperationssequentiallyLabel.Text = 'Perform all operations sequentially';

            % Create MO_Title
            app.MO_Title = uilabel(app.MorphologyopeningTab);
            app.MO_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.MO_Title.HorizontalAlignment = 'center';
            app.MO_Title.FontWeight = 'bold';
            app.MO_Title.Position = [10 634 944 22];
            app.MO_Title.Text = 'Morphology opening: simplify microstrucutre to improve mesh quality and ease mesh generation (optional)';

            % Create MO_text_1
            app.MO_text_1 = uilabel(app.MorphologyopeningTab);
            app.MO_text_1.Position = [10 598 953 22];
            app.MO_text_1.Text = 'Choose the operation(s) you want to apply for each domain, set operation parameters, and then apply changes clicking on the Do button. If satisfied, click the save cell button.';

            % Create MO_SavecellButton
            app.MO_SavecellButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_SavecellButton.ButtonPushedFcn = createCallbackFcn(app, @MO_SavecellButtonPushed, true);
            app.MO_SavecellButton.BackgroundColor = [0.8 0.8 0.8];
            app.MO_SavecellButton.Enable = 'off';
            app.MO_SavecellButton.Position = [460 9 103 22];
            app.MO_SavecellButton.Text = 'Save cell';

            % Create MO_text_2
            app.MO_text_2 = uilabel(app.MorphologyopeningTab);
            app.MO_text_2.HorizontalAlignment = 'center';
            app.MO_text_2.FontAngle = 'italic';
            app.MO_text_2.Position = [8 9 439 51];
            app.MO_text_2.Text = {' Modfications are saved only if you click on the green button ''Save''.'; 'Once finished, you must click on ''Save cell''.  Alternatively you can skip this step.'};

            % Create MO_VisualizecellButton
            app.MO_VisualizecellButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_VisualizecellButton.ButtonPushedFcn = createCallbackFcn(app, @MO_VisualizecellButtonPushed, true);
            app.MO_VisualizecellButton.BackgroundColor = [0.8 0.8 0.8];
            app.MO_VisualizecellButton.Enable = 'off';
            app.MO_VisualizecellButton.Tooltip = {'Three volumes are compared:'; 'Left: initial'; 'Center: saved'; 'Right: current'};
            app.MO_VisualizecellButton.Position = [460 38 225 22];
            app.MO_VisualizecellButton.Text = 'Visualize cell';

            % Create MO_vf_UITable
            app.MO_vf_UITable = uitable(app.MorphologyopeningTab);
            app.MO_vf_UITable.ColumnName = {'Domain'; 'Phases'; 'Before'; 'After'};
            app.MO_vf_UITable.ColumnWidth = {'1x', 'auto', 'auto', 'auto'};
            app.MO_vf_UITable.RowName = {};
            app.MO_vf_UITable.Position = [482 297 472 264];

            % Create MO_VolumefractionsLabel
            app.MO_VolumefractionsLabel = uilabel(app.MorphologyopeningTab);
            app.MO_VolumefractionsLabel.HorizontalAlignment = 'center';
            app.MO_VolumefractionsLabel.FontWeight = 'bold';
            app.MO_VolumefractionsLabel.Position = [482 566 472 22];
            app.MO_VolumefractionsLabel.Text = 'Volume fractions';

            % Create MO_operations_UITable
            app.MO_operations_UITable = uitable(app.MorphologyopeningTab);
            app.MO_operations_UITable.ColumnName = {'Domain'; 'Operation'; 'Do'};
            app.MO_operations_UITable.ColumnWidth = {'1x', 'auto', 'auto'};
            app.MO_operations_UITable.RowName = {};
            app.MO_operations_UITable.ColumnEditable = [false false true];
            app.MO_operations_UITable.Position = [10 297 440 264];

            % Create MO_VolumefractionsLabel_2
            app.MO_VolumefractionsLabel_2 = uilabel(app.MorphologyopeningTab);
            app.MO_VolumefractionsLabel_2.HorizontalAlignment = 'center';
            app.MO_VolumefractionsLabel_2.FontWeight = 'bold';
            app.MO_VolumefractionsLabel_2.Position = [10 566 440 22];
            app.MO_VolumefractionsLabel_2.Text = 'Operations';

            % Create MO_text_4
            app.MO_text_4 = uilabel(app.MorphologyopeningTab);
            app.MO_text_4.FontSize = 14;
            app.MO_text_4.FontWeight = 'bold';
            app.MO_text_4.Position = [10 262 177 22];
            app.MO_text_4.Text = 'Parameters';

            % Create MO_text_3
            app.MO_text_3 = uilabel(app.MorphologyopeningTab);
            app.MO_text_3.FontWeight = 'bold';
            app.MO_text_3.Position = [10 234 177 22];
            app.MO_text_3.Text = 'Step 1: Erosion + dilatation';

            % Create ErosiondistanceinvoxellengthEditFieldLabel
            app.ErosiondistanceinvoxellengthEditFieldLabel = uilabel(app.MorphologyopeningTab);
            app.ErosiondistanceinvoxellengthEditFieldLabel.HorizontalAlignment = 'right';
            app.ErosiondistanceinvoxellengthEditFieldLabel.Position = [10 181 182 22];
            app.ErosiondistanceinvoxellengthEditFieldLabel.Text = 'Erosion distance (in voxel length)';

            % Create ErosiondistanceinvoxellengthEditField
            app.ErosiondistanceinvoxellengthEditField = uieditfield(app.MorphologyopeningTab, 'numeric');
            app.ErosiondistanceinvoxellengthEditField.Limits = [0 Inf];
            app.ErosiondistanceinvoxellengthEditField.RoundFractionalValues = 'on';
            app.ErosiondistanceinvoxellengthEditField.Position = [218 181 42 22];
            app.ErosiondistanceinvoxellengthEditField.Value = 1;

            % Create DilatationdistanceinvoxellengthEditFieldLabel
            app.DilatationdistanceinvoxellengthEditFieldLabel = uilabel(app.MorphologyopeningTab);
            app.DilatationdistanceinvoxellengthEditFieldLabel.HorizontalAlignment = 'right';
            app.DilatationdistanceinvoxellengthEditFieldLabel.Position = [10 155 193 22];
            app.DilatationdistanceinvoxellengthEditFieldLabel.Text = 'Dilatation distance (in voxel length)';

            % Create DilatationdistanceinvoxellengthEditField
            app.DilatationdistanceinvoxellengthEditField = uieditfield(app.MorphologyopeningTab, 'numeric');
            app.DilatationdistanceinvoxellengthEditField.Limits = [0 Inf];
            app.DilatationdistanceinvoxellengthEditField.RoundFractionalValues = 'on';
            app.DilatationdistanceinvoxellengthEditField.Position = [218 155 42 22];
            app.DilatationdistanceinvoxellengthEditField.Value = 1;

            % Create NumberofiterationEditFieldLabel
            app.NumberofiterationEditFieldLabel = uilabel(app.MorphologyopeningTab);
            app.NumberofiterationEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofiterationEditFieldLabel.Position = [10 129 108 22];
            app.NumberofiterationEditFieldLabel.Text = 'Number of iteration';

            % Create ErosiondilatationNumberofiterationEditField
            app.ErosiondilatationNumberofiterationEditField = uieditfield(app.MorphologyopeningTab, 'numeric');
            app.ErosiondilatationNumberofiterationEditField.Limits = [0 Inf];
            app.ErosiondilatationNumberofiterationEditField.RoundFractionalValues = 'on';
            app.ErosiondilatationNumberofiterationEditField.Position = [218 129 42 22];
            app.ErosiondilatationNumberofiterationEditField.Value = 1;

            % Create PerformDropDownLabel
            app.PerformDropDownLabel = uilabel(app.MorphologyopeningTab);
            app.PerformDropDownLabel.HorizontalAlignment = 'right';
            app.PerformDropDownLabel.Position = [10 103 48 22];
            app.PerformDropDownLabel.Text = 'Perform';

            % Create PerformErosionDropDown
            app.PerformErosionDropDown = uidropdown(app.MorphologyopeningTab);
            app.PerformErosionDropDown.Items = {'All solid phases in one step', 'Solid phase per solid phase'};
            app.PerformErosionDropDown.Tooltip = {'All solid phases in one step: erosion dilatation will be applied on the solid union phase. Then solid voxels will be assigned to the nearest phase from the initial volume.'; 'Solid phase per solid phase: erosion dilatation will be applied sequentially (the latter phases overwritting the first phases).'};
            app.PerformErosionDropDown.Position = [73 103 187 22];
            app.PerformErosionDropDown.Value = 'All solid phases in one step';

            % Create ErositiondilatationisalwaysperformedfirstLabel
            app.ErositiondilatationisalwaysperformedfirstLabel = uilabel(app.MorphologyopeningTab);
            app.ErositiondilatationisalwaysperformedfirstLabel.FontAngle = 'italic';
            app.ErositiondilatationisalwaysperformedfirstLabel.Position = [12 207 260 23];
            app.ErositiondilatationisalwaysperformedfirstLabel.Text = 'Erosition dilatation is always performed first';

            % Create MO_text_5
            app.MO_text_5 = uilabel(app.MorphologyopeningTab);
            app.MO_text_5.FontWeight = 'bold';
            app.MO_text_5.Position = [284 207 177 22];
            app.MO_text_5.Text = 'Clean voxel connections';

            % Create DetectandremoveverticeverticeandedgeedgeconnectionsLabel
            app.DetectandremoveverticeverticeandedgeedgeconnectionsLabel = uilabel(app.MorphologyopeningTab);
            app.DetectandremoveverticeverticeandedgeedgeconnectionsLabel.Position = [284 181 341 23];
            app.DetectandremoveverticeverticeandedgeedgeconnectionsLabel.Text = 'Detect and remove vertice-vertice and edge-edge connections';

            % Create PerformDropDown_2Label
            app.PerformDropDown_2Label = uilabel(app.MorphologyopeningTab);
            app.PerformDropDown_2Label.HorizontalAlignment = 'right';
            app.PerformDropDown_2Label.Position = [284 155 48 22];
            app.PerformDropDown_2Label.Text = 'Perform';

            % Create PerformCleanConnectionDropDown
            app.PerformCleanConnectionDropDown = uidropdown(app.MorphologyopeningTab);
            app.PerformCleanConnectionDropDown.Items = {'All solid phases in one step', 'Solid phase per solid phase'};
            app.PerformCleanConnectionDropDown.Position = [347 155 187 22];
            app.PerformCleanConnectionDropDown.Value = 'All solid phases in one step';

            % Create MO_text_6
            app.MO_text_6 = uilabel(app.MorphologyopeningTab);
            app.MO_text_6.FontWeight = 'bold';
            app.MO_text_6.Position = [284 129 222 22];
            app.MO_text_6.Text = 'Convert to unique connected clusters';

            % Create GrainvolumefilterinnumberofvoxelEditFieldLabel
            app.GrainvolumefilterinnumberofvoxelEditFieldLabel = uilabel(app.MorphologyopeningTab);
            app.GrainvolumefilterinnumberofvoxelEditFieldLabel.FontWeight = 'bold';
            app.GrainvolumefilterinnumberofvoxelEditFieldLabel.Position = [647 234 232 22];
            app.GrainvolumefilterinnumberofvoxelEditFieldLabel.Text = 'Grain volume filter (in number of voxel)';

            % Create GrainvolumefilterinnumberofvoxelEditField
            app.GrainvolumefilterinnumberofvoxelEditField = uieditfield(app.MorphologyopeningTab, 'numeric');
            app.GrainvolumefilterinnumberofvoxelEditField.Limits = [0 Inf];
            app.GrainvolumefilterinnumberofvoxelEditField.RoundFractionalValues = 'on';
            app.GrainvolumefilterinnumberofvoxelEditField.Position = [890 234 42 22];
            app.GrainvolumefilterinnumberofvoxelEditField.Value = 1000;

            % Create MO_DoButton
            app.MO_DoButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_DoButton.ButtonPushedFcn = createCallbackFcn(app, @MO_DoButtonPushed, true);
            app.MO_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.MO_DoButton.FontSize = 14;
            app.MO_DoButton.FontWeight = 'bold';
            app.MO_DoButton.FontColor = [1 1 1];
            app.MO_DoButton.Enable = 'off';
            app.MO_DoButton.Position = [721 139 71 25];
            app.MO_DoButton.Text = 'Do';

            % Create MO_UndoButton
            app.MO_UndoButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @MO_UndoButtonPushed, true);
            app.MO_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.MO_UndoButton.FontSize = 14;
            app.MO_UndoButton.FontWeight = 'bold';
            app.MO_UndoButton.FontColor = [1 1 1];
            app.MO_UndoButton.Enable = 'off';
            app.MO_UndoButton.Position = [802 139 71 25];
            app.MO_UndoButton.Text = 'Undo';

            % Create MO_SaveButton
            app.MO_SaveButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @MO_SaveButtonPushed, true);
            app.MO_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.MO_SaveButton.FontSize = 14;
            app.MO_SaveButton.FontWeight = 'bold';
            app.MO_SaveButton.FontColor = [1 1 1];
            app.MO_SaveButton.Enable = 'off';
            app.MO_SaveButton.Position = [883 139 71 25];
            app.MO_SaveButton.Text = 'Save';

            % Create MO_statut_Label
            app.MO_statut_Label = uilabel(app.MorphologyopeningTab);
            app.MO_statut_Label.HorizontalAlignment = 'center';
            app.MO_statut_Label.FontAngle = 'italic';
            app.MO_statut_Label.Position = [721 54 233 51];
            app.MO_statut_Label.Text = 'Idle';

            % Create MO_text_7
            app.MO_text_7 = uilabel(app.MorphologyopeningTab);
            app.MO_text_7.FontWeight = 'bold';
            app.MO_text_7.Position = [284 234 279 22];
            app.MO_text_7.Text = 'Subsequent steps (in a loop until convergence)';

            % Create MO_PreservebackgroundCheckBox
            app.MO_PreservebackgroundCheckBox = uicheckbox(app.MorphologyopeningTab);
            app.MO_PreservebackgroundCheckBox.Tooltip = {'If true, after erosion dilatation, voxels will be reassigned so that background is unchanged.'};
            app.MO_PreservebackgroundCheckBox.Text = 'Preserve background';
            app.MO_PreservebackgroundCheckBox.Position = [10 77 246 22];

            % Create MO_UnionsoliduniqueclusterCheckBox
            app.MO_UnionsoliduniqueclusterCheckBox = uicheckbox(app.MorphologyopeningTab);
            app.MO_UnionsoliduniqueclusterCheckBox.Tooltip = {'If true, isolated solid clusters are removed and assigned to pore phase'};
            app.MO_UnionsoliduniqueclusterCheckBox.Text = 'Enforce union of solid phase is a unique connected cluster';
            app.MO_UnionsoliduniqueclusterCheckBox.Position = [288 79 337 22];
            app.MO_UnionsoliduniqueclusterCheckBox.Value = true;

            % Create MO_PoreuniqueclusterCheckBox
            app.MO_PoreuniqueclusterCheckBox = uicheckbox(app.MorphologyopeningTab);
            app.MO_PoreuniqueclusterCheckBox.Tooltip = {'If true, isolated pore clusters are removed and assigned to solid phase'};
            app.MO_PoreuniqueclusterCheckBox.Text = 'Enforce pore background is a unique connected cluster';
            app.MO_PoreuniqueclusterCheckBox.Position = [288 104 320 22];
            app.MO_PoreuniqueclusterCheckBox.Value = true;

            % Create MO_SkipButton
            app.MO_SkipButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_SkipButton.ButtonPushedFcn = createCallbackFcn(app, @MO_SkipButtonPushed, true);
            app.MO_SkipButton.BackgroundColor = [0.8 0.8 0.8];
            app.MO_SkipButton.Enable = 'off';
            app.MO_SkipButton.Position = [582 9 103 22];
            app.MO_SkipButton.Text = 'Skip step';

            % Create MO_UndoAllButton
            app.MO_UndoAllButton = uibutton(app.MorphologyopeningTab, 'push');
            app.MO_UndoAllButton.ButtonPushedFcn = createCallbackFcn(app, @MO_UndoAllButtonPushed, true);
            app.MO_UndoAllButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.MO_UndoAllButton.FontSize = 14;
            app.MO_UndoAllButton.FontWeight = 'bold';
            app.MO_UndoAllButton.FontColor = [1 1 1];
            app.MO_UndoAllButton.Enable = 'off';
            app.MO_UndoAllButton.Position = [803 108 71 25];
            app.MO_UndoAllButton.Text = 'Undo all';

            % Create CheckgraincontiguityandfixitifdiscontinousslowCheckBox
            app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox = uicheckbox(app.MorphologyopeningTab);
            app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Text = 'Check grain contiguity and fix it if discontinous (slow)';
            app.CheckgraincontiguityandfixitifdiscontinousslowCheckBox.Position = [647 210 307 22];

            % Create AssemblecellTab
            app.AssemblecellTab = uitab(app.TabGroup);
            app.AssemblecellTab.Title = 'Assemble cell';
            app.AssemblecellTab.ForegroundColor = [0 0 1];

            % Create Assemble_Title
            app.Assemble_Title = uilabel(app.AssemblecellTab);
            app.Assemble_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Assemble_Title.HorizontalAlignment = 'center';
            app.Assemble_Title.FontWeight = 'bold';
            app.Assemble_Title.Position = [10 634 944 22];
            app.Assemble_Title.Text = 'Build 3D array that will be used for the meshing generation';

            % Create Assemble_Instructions
            app.Assemble_Instructions = uilabel(app.AssemblecellTab);
            app.Assemble_Instructions.Position = [10 598 856 22];
            app.Assemble_Instructions.Text = 'Assign each phase to a unique id and name it. Phase id must be positive non zero integer. Group Id must be positive integer. Default values are pre-entered.';

            % Create Assemble_Lamp
            app.Assemble_Lamp = uilamp(app.AssemblecellTab);
            app.Assemble_Lamp.Position = [901 9 53 53];
            app.Assemble_Lamp.Color = [1 0 0];

            % Create TabisNOTsetupcorrectlyLabel_4
            app.TabisNOTsetupcorrectlyLabel_4 = uilabel(app.AssemblecellTab);
            app.TabisNOTsetupcorrectlyLabel_4.HorizontalAlignment = 'right';
            app.TabisNOTsetupcorrectlyLabel_4.FontSize = 14;
            app.TabisNOTsetupcorrectlyLabel_4.FontWeight = 'bold';
            app.TabisNOTsetupcorrectlyLabel_4.Position = [703 25 183 22];
            app.TabisNOTsetupcorrectlyLabel_4.Text = 'Tab is NOT setup correctly';

            % Create Assemble_UITable
            app.Assemble_UITable = uitable(app.AssemblecellTab);
            app.Assemble_UITable.ColumnName = {'Domain'; 'Phase Id'; 'Assigned Id in cell'; 'Phase name'; 'Group Id'};
            app.Assemble_UITable.ColumnWidth = {'1x', 'auto', 'auto', 'auto', 'auto'};
            app.Assemble_UITable.RowName = {};
            app.Assemble_UITable.ColumnEditable = [false false true true true];
            app.Assemble_UITable.CellEditCallback = createCallbackFcn(app, @Assemble_UITableCellEdit, true);
            app.Assemble_UITable.Position = [10 268 520 310];

            % Create Assemble_Label_4
            app.Assemble_Label_4 = uilabel(app.AssemblecellTab);
            app.Assemble_Label_4.HorizontalAlignment = 'center';
            app.Assemble_Label_4.FontAngle = 'italic';
            app.Assemble_Label_4.Position = [11 9 452 51];
            app.Assemble_Label_4.Text = 'Once done, please click on ''Save cell'' to go on.';

            % Create Assemble_SavecellButton
            app.Assemble_SavecellButton = uibutton(app.AssemblecellTab, 'push');
            app.Assemble_SavecellButton.ButtonPushedFcn = createCallbackFcn(app, @Assemble_SavecellButtonPushed, true);
            app.Assemble_SavecellButton.BackgroundColor = [0.8 0.8 0.8];
            app.Assemble_SavecellButton.Enable = 'off';
            app.Assemble_SavecellButton.Position = [460 9 225 22];
            app.Assemble_SavecellButton.Text = 'Save cell';

            % Create Assemble_VisualizecellButton
            app.Assemble_VisualizecellButton = uibutton(app.AssemblecellTab, 'push');
            app.Assemble_VisualizecellButton.ButtonPushedFcn = createCallbackFcn(app, @Assemble_VisualizecellButtonPushed, true);
            app.Assemble_VisualizecellButton.BackgroundColor = [0.8 0.8 0.8];
            app.Assemble_VisualizecellButton.Enable = 'off';
            app.Assemble_VisualizecellButton.Tooltip = {'Three volumes are compared:'; 'Left: initial phase id'; 'Center: assigned id in cell'; 'Right: group id'};
            app.Assemble_VisualizecellButton.Position = [460 38 225 22];
            app.Assemble_VisualizecellButton.Text = 'Visualize cell';

            % Create Label_4
            app.Label_4 = uilabel(app.AssemblecellTab);
            app.Label_4.VerticalAlignment = 'top';
            app.Label_4.FontAngle = 'italic';
            app.Label_4.Position = [547 329 393 249];
            app.Label_4.Text = {'Groups are used for segregated models, for which domains are not'; 'solved at once (monolithic model) but sequentially, usually embedded'; 'in a self-converging iterative numerical scheme.'; ''; 'Phase with the same group id will be meshed together.'; 'Meshes from different group have a conforming interface.'; ''; 'If your model is monolithic and/or you will select the ''One mesh for'; 'the whole geometry'' option in the create mesh tab, then you can'; 'ignore this fifth coloum.'; ''; 'Phase name will be used for mesh filename'};

            % Create UpdateButton
            app.UpdateButton = uibutton(app.AssemblecellTab, 'push');
            app.UpdateButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateButtonPushed, true);
            app.UpdateButton.Position = [547 268 100 22];
            app.UpdateButton.Text = 'Update';

            % Create UpdateButton_label
            app.UpdateButton_label = uilabel(app.AssemblecellTab);
            app.UpdateButton_label.FontAngle = 'italic';
            app.UpdateButton_label.Position = [547 296 310 22];
            app.UpdateButton_label.Text = 'Once table is modified, please click on the update button';

            % Create Group_UITable
            app.Group_UITable = uitable(app.AssemblecellTab);
            app.Group_UITable.ColumnName = {'Group Id'; 'Group Name'};
            app.Group_UITable.ColumnWidth = {'1x', 'auto'};
            app.Group_UITable.RowName = {};
            app.Group_UITable.ColumnEditable = [false true];
            app.Group_UITable.Position = [10 100 520 155];

            % Create UpdateButton_label_2
            app.UpdateButton_label_2 = uilabel(app.AssemblecellTab);
            app.UpdateButton_label_2.FontAngle = 'italic';
            app.UpdateButton_label_2.Position = [547 199 239 22];
            app.UpdateButton_label_2.Text = 'Group name will be used for mesh filename';

            % Create Meshgenerationchoice
            app.Meshgenerationchoice = uitab(app.TabGroup);
            app.Meshgenerationchoice.Title = 'Mesh generation choice';
            app.Meshgenerationchoice.ForegroundColor = [0.851 0.3255 0.098];

            % Create Meshgen_Title
            app.Meshgen_Title = uilabel(app.Meshgenerationchoice);
            app.Meshgen_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Meshgen_Title.HorizontalAlignment = 'center';
            app.Meshgen_Title.FontWeight = 'bold';
            app.Meshgen_Title.Position = [10 634 944 22];
            app.Meshgen_Title.Text = 'Choose which method you want to use';

            % Create Meshgen_Iso2mesh
            app.Meshgen_Iso2mesh = uibutton(app.Meshgenerationchoice, 'push');
            app.Meshgen_Iso2mesh.ButtonPushedFcn = createCallbackFcn(app, @Meshgen_Iso2meshButtonPushed, true);
            app.Meshgen_Iso2mesh.BackgroundColor = [0.8 0.8 0.8];
            app.Meshgen_Iso2mesh.Enable = 'off';
            app.Meshgen_Iso2mesh.Position = [10 566 288 51];
            app.Meshgen_Iso2mesh.Text = 'Iso2mesh';

            % Create Meshgen_Cuboidmesh
            app.Meshgen_Cuboidmesh = uibutton(app.Meshgenerationchoice, 'push');
            app.Meshgen_Cuboidmesh.ButtonPushedFcn = createCallbackFcn(app, @Meshgen_CuboidmeshButtonPushed, true);
            app.Meshgen_Cuboidmesh.BackgroundColor = [0.8 0.8 0.8];
            app.Meshgen_Cuboidmesh.Enable = 'off';
            app.Meshgen_Cuboidmesh.Position = [10 495 288 51];
            app.Meshgen_Cuboidmesh.Text = 'Cuboid mesh';

            % Create Label_5
            app.Label_5 = uilabel(app.Meshgenerationchoice);
            app.Label_5.Position = [307 566 645 51];
            app.Label_5.Text = {'Tetrahedron-based unstructured mesh with mesh density control and smooth surface. Time and RAM expensive and'; 'may fail is volume too large and/or too complex.'};

            % Create Label_6
            app.Label_6 = uilabel(app.Meshgenerationchoice);
            app.Label_6.Position = [307 495 644 51];
            app.Label_6.Text = {'Tetrahedron-based structured mesh with cuboid representation. Fast and robust but vertices expensive method and'; 'no surface smoothing and no mesh density control.'};

            % Create Label_7
            app.Label_7 = uilabel(app.Meshgenerationchoice);
            app.Label_7.FontAngle = 'italic';
            app.Label_7.Position = [10 443 824 32];
            app.Label_7.Text = {'You will be bring to the orange coloured tab associated with your choice. You can ignore the other orange coloured tab.'; 'If you select the ''Iso2mesh'' choice, you will have to quote an additional article in addition to quote this toolbox (see About tab for quotation instructions)'};

            % Create Iso2meshoptionsTab
            app.Iso2meshoptionsTab = uitab(app.TabGroup);
            app.Iso2meshoptionsTab.Title = 'Iso2mesh options';
            app.Iso2meshoptionsTab.ForegroundColor = [0.851 0.3255 0.098];

            % Create Iso2mesh_Title
            app.Iso2mesh_Title = uilabel(app.Iso2meshoptionsTab);
            app.Iso2mesh_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Iso2mesh_Title.HorizontalAlignment = 'center';
            app.Iso2mesh_Title.FontWeight = 'bold';
            app.Iso2mesh_Title.Position = [10 634 944 22];
            app.Iso2mesh_Title.Text = 'Choose Iso2mesh options';

            % Create CGAL_UITable
            app.CGAL_UITable = uitable(app.Iso2meshoptionsTab);
            app.CGAL_UITable.ColumnName = {'Parameters'; 'Description'; 'Value'};
            app.CGAL_UITable.ColumnWidth = {'auto', '1x', 'auto'};
            app.CGAL_UITable.RowName = {};
            app.CGAL_UITable.ColumnEditable = [false false true];
            app.CGAL_UITable.Position = [10 493 938 100];

            % Create Iso2mesh_label
            app.Iso2mesh_label = uilabel(app.Iso2meshoptionsTab);
            app.Iso2mesh_label.FontWeight = 'bold';
            app.Iso2mesh_label.Position = [10 598 121 22];
            app.Iso2mesh_label.Text = 'CGAL mesh options';

            % Create Iso2mesh_label_2
            app.Iso2mesh_label_2 = uilabel(app.Iso2meshoptionsTab);
            app.Iso2mesh_label_2.FontWeight = 'bold';
            app.Iso2mesh_label_2.Position = [10 451 149 22];
            app.Iso2mesh_label_2.Text = 'Smoothing mesh options';

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.Iso2meshoptionsTab);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.Position = [10 424 46 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create Surfacemesh_MethodDropDown
            app.Surfacemesh_MethodDropDown = uidropdown(app.Iso2meshoptionsTab);
            app.Surfacemesh_MethodDropDown.Items = {'lowpass', 'laplacian', 'laplacianhc'};
            app.Surfacemesh_MethodDropDown.ValueChangedFcn = createCallbackFcn(app, @Surfacemesh_MethodDropDownValueChanged, true);
            app.Surfacemesh_MethodDropDown.Position = [71 424 100 22];
            app.Surfacemesh_MethodDropDown.Value = 'lowpass';

            % Create Smoothing_UITable
            app.Smoothing_UITable = uitable(app.Iso2meshoptionsTab);
            app.Smoothing_UITable.ColumnName = {'Parameters'; 'Description'; 'Value'};
            app.Smoothing_UITable.ColumnWidth = {'auto', '1x', 'auto'};
            app.Smoothing_UITable.RowName = {};
            app.Smoothing_UITable.ColumnEditable = [false false true];
            app.Smoothing_UITable.Position = [10 286 938 100];

            % Create Iso2mesh_label_3
            app.Iso2mesh_label_3 = uilabel(app.Iso2meshoptionsTab);
            app.Iso2mesh_label_3.FontWeight = 'bold';
            app.Iso2mesh_label_3.Position = [13 244 128 22];
            app.Iso2mesh_label_3.Text = 'Mesh density options';

            % Create Meshdensity_UITable
            app.Meshdensity_UITable = uitable(app.Iso2meshoptionsTab);
            app.Meshdensity_UITable.ColumnName = {'Parameters'; 'Description'; 'Value'};
            app.Meshdensity_UITable.ColumnWidth = {'auto', '1x', 'auto'};
            app.Meshdensity_UITable.RowName = {};
            app.Meshdensity_UITable.ColumnEditable = [false false true];
            app.Meshdensity_UITable.Position = [13 139 938 100];

            % Create Label_12
            app.Label_12 = uilabel(app.Iso2meshoptionsTab);
            app.Label_12.FontAngle = 'italic';
            app.Label_12.Position = [10 391 937 28];
            app.Label_12.Text = {'''Lowpass'' method outperforms ''Laplacian-HC'' in volume preserving and both are significantly better than the standard Laplacian method. See R. Bade, H. Haase, B. Preim,'; '"Comparison of Fundamental Mesh Smoothing Algorithms for Medical Surface Models," Simulation and Visualization, pp. 289-304, 2006. '};

            % Create Iso2meshversionnumberLabel
            app.Iso2meshversionnumberLabel = uilabel(app.Iso2meshoptionsTab);
            app.Iso2meshversionnumberLabel.FontAngle = 'italic';
            app.Iso2meshversionnumberLabel.Position = [10 12 284 22];
            app.Iso2meshversionnumberLabel.Text = 'Iso2mesh version number';

            % Create StructuredmeshoptionsTab
            app.StructuredmeshoptionsTab = uitab(app.TabGroup);
            app.StructuredmeshoptionsTab.Title = 'Structured mesh options';
            app.StructuredmeshoptionsTab.ForegroundColor = [0.851 0.3255 0.098];

            % Create Regularmesh_Title
            app.Regularmesh_Title = uilabel(app.StructuredmeshoptionsTab);
            app.Regularmesh_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Regularmesh_Title.HorizontalAlignment = 'center';
            app.Regularmesh_Title.FontWeight = 'bold';
            app.Regularmesh_Title.Position = [10 634 944 22];
            app.Regularmesh_Title.Text = 'Choose the voxel mesh representation';

            % Create EachvoxelwillberepresentedwithDropDownLabel
            app.EachvoxelwillberepresentedwithDropDownLabel = uilabel(app.StructuredmeshoptionsTab);
            app.EachvoxelwillberepresentedwithDropDownLabel.HorizontalAlignment = 'right';
            app.EachvoxelwillberepresentedwithDropDownLabel.Position = [10 599 194 22];
            app.EachvoxelwillberepresentedwithDropDownLabel.Text = 'Each voxel will be represented with';

            % Create EachvoxelwillberepresentedwithDropDown
            app.EachvoxelwillberepresentedwithDropDown = uidropdown(app.StructuredmeshoptionsTab);
            app.EachvoxelwillberepresentedwithDropDown.Items = {'6 tetrahedrons', '24 tetrahedrons'};
            app.EachvoxelwillberepresentedwithDropDown.ValueChangedFcn = createCallbackFcn(app, @EachvoxelwillberepresentedwithDropDownValueChanged, true);
            app.EachvoxelwillberepresentedwithDropDown.Position = [219 599 136 22];
            app.EachvoxelwillberepresentedwithDropDown.Value = '6 tetrahedrons';

            % Create Regularmesh_Image5
            app.Regularmesh_Image5 = uiimage(app.StructuredmeshoptionsTab);
            app.Regularmesh_Image5.Position = [10 232 944 346];
            app.Regularmesh_Image5.ImageSource = 'Voxel_5tets.png';

            % Create Regularmesh_Image24
            app.Regularmesh_Image24 = uiimage(app.StructuredmeshoptionsTab);
            app.Regularmesh_Image24.Visible = 'off';
            app.Regularmesh_Image24.Position = [10 228 944 346];
            app.Regularmesh_Image24.ImageSource = 'Voxel_24tets.png';

            % Create Regularmesh_Image6
            app.Regularmesh_Image6 = uiimage(app.StructuredmeshoptionsTab);
            app.Regularmesh_Image6.Visible = 'off';
            app.Regularmesh_Image6.Position = [10 228 944 346];
            app.Regularmesh_Image6.ImageSource = 'Voxel_6tets.png';

            % Create CreatemeshTab
            app.CreatemeshTab = uitab(app.TabGroup);
            app.CreatemeshTab.Title = 'Create mesh';
            app.CreatemeshTab.ForegroundColor = [1 0 0];

            % Create Createmesh_Title
            app.Createmesh_Title = uilabel(app.CreatemeshTab);
            app.Createmesh_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Createmesh_Title.HorizontalAlignment = 'center';
            app.Createmesh_Title.FontWeight = 'bold';
            app.Createmesh_Title.Position = [10 634 944 22];
            app.Createmesh_Title.Text = 'Create, visualize and save mesh';

            % Create Createmesh_button
            app.Createmesh_button = uibutton(app.CreatemeshTab, 'push');
            app.Createmesh_button.ButtonPushedFcn = createCallbackFcn(app, @Createmesh_buttonButtonPushed, true);
            app.Createmesh_button.BackgroundColor = [0.8 0.8 0.8];
            app.Createmesh_button.Enable = 'off';
            app.Createmesh_button.Position = [10 427 329 51];
            app.Createmesh_button.Text = 'Create mesh';

            % Create Meshchoice_UITable
            app.Meshchoice_UITable = uitable(app.CreatemeshTab);
            app.Meshchoice_UITable.ColumnName = {'What do you want to mesh ?'; 'Choice'};
            app.Meshchoice_UITable.ColumnWidth = {'1x'};
            app.Meshchoice_UITable.RowName = {};
            app.Meshchoice_UITable.ColumnEditable = [false true];
            app.Meshchoice_UITable.Position = [10 488 329 100];

            % Create Label_8
            app.Label_8 = uilabel(app.CreatemeshTab);
            app.Label_8.Position = [10 598 665 22];
            app.Label_8.Text = 'Select the mesh(es) you want to generate and click on the button ''create mesh''. You can then visualize and save meshes.';

            % Create Mesh_VisualizeButton
            app.Mesh_VisualizeButton = uibutton(app.CreatemeshTab, 'push');
            app.Mesh_VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @Mesh_VisualizeButtonPushed, true);
            app.Mesh_VisualizeButton.BackgroundColor = [0.8 0.8 0.8];
            app.Mesh_VisualizeButton.Enable = 'off';
            app.Mesh_VisualizeButton.Position = [10 106 329 51];
            app.Mesh_VisualizeButton.Text = 'Visualize';

            % Create Cellindexstartsat0CheckBox
            app.Cellindexstartsat0CheckBox = uicheckbox(app.CreatemeshTab);
            app.Cellindexstartsat0CheckBox.Text = 'Cell index starts at 0 (otherwise at 1) for .mat and .csv files.';
            app.Cellindexstartsat0CheckBox.Position = [385 171 342 22];
            app.Cellindexstartsat0CheckBox.Value = true;

            % Create Createmesh_choicemesh
            app.Createmesh_choicemesh = uilabel(app.CreatemeshTab);
            app.Createmesh_choicemesh.FontWeight = 'bold';
            app.Createmesh_choicemesh.FontColor = [1 0 0];
            app.Createmesh_choicemesh.Position = [348 456 280 22];
            app.Createmesh_choicemesh.Text = 'Mesh method not selected';

            % Create Mesh_SaveButton
            app.Mesh_SaveButton = uibutton(app.CreatemeshTab, 'push');
            app.Mesh_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Mesh_SaveButtonPushed, true);
            app.Mesh_SaveButton.BackgroundColor = [0.8 0.8 0.8];
            app.Mesh_SaveButton.Enable = 'off';
            app.Mesh_SaveButton.Position = [385 106 329 51];
            app.Mesh_SaveButton.Text = 'Save mesh';

            % Create StatutTextAreaLabel
            app.StatutTextAreaLabel = uilabel(app.CreatemeshTab);
            app.StatutTextAreaLabel.HorizontalAlignment = 'right';
            app.StatutTextAreaLabel.Position = [626 32 40 22];
            app.StatutTextAreaLabel.Text = 'Statut:';

            % Create Createmesh_StatutTextArea
            app.Createmesh_StatutTextArea = uitextarea(app.CreatemeshTab);
            app.Createmesh_StatutTextArea.Editable = 'off';
            app.Createmesh_StatutTextArea.Position = [681 13 273 43];
            app.Createmesh_StatutTextArea.Value = {'Idle'};

            % Create CreateMesh_text_1
            app.CreateMesh_text_1 = uilabel(app.CreatemeshTab);
            app.CreateMesh_text_1.HorizontalAlignment = 'center';
            app.CreateMesh_text_1.FontSize = 14;
            app.CreateMesh_text_1.FontWeight = 'bold';
            app.CreateMesh_text_1.Position = [8 379 331 22];
            app.CreateMesh_text_1.Text = 'Mesh visualization options';

            % Create LightingDropDownLabel
            app.LightingDropDownLabel = uilabel(app.CreatemeshTab);
            app.LightingDropDownLabel.HorizontalAlignment = 'right';
            app.LightingDropDownLabel.Position = [8 293 48 22];
            app.LightingDropDownLabel.Text = 'Lighting';

            % Create LightingDropDown
            app.LightingDropDown = uidropdown(app.CreatemeshTab);
            app.LightingDropDown.Items = {'none', 'flat', 'gouraud'};
            app.LightingDropDown.Position = [94 293 134 22];
            app.LightingDropDown.Value = 'none';

            % Create MeshedgesDropDownLabel
            app.MeshedgesDropDownLabel = uilabel(app.CreatemeshTab);
            app.MeshedgesDropDownLabel.HorizontalAlignment = 'right';
            app.MeshedgesDropDownLabel.Position = [8 320 71 22];
            app.MeshedgesDropDownLabel.Text = 'Mesh edges';

            % Create MeshedgesDropDown
            app.MeshedgesDropDown = uidropdown(app.CreatemeshTab);
            app.MeshedgesDropDown.Items = {'Visible', 'Not visible'};
            app.MeshedgesDropDown.Position = [94 320 134 22];
            app.MeshedgesDropDown.Value = 'Visible';

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.CreatemeshTab);
            app.ColormapDropDownLabel.HorizontalAlignment = 'right';
            app.ColormapDropDownLabel.Position = [8 239 61 22];
            app.ColormapDropDownLabel.Text = 'Color map';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.CreatemeshTab);
            app.ColormapDropDown.Items = {'Default', 'Jet', 'Turbo', 'gray', 'bone', 'copper'};
            app.ColormapDropDown.Position = [94 239 134 22];
            app.ColormapDropDown.Value = 'Default';

            % Create Transparency01EditFieldLabel
            app.Transparency01EditFieldLabel = uilabel(app.CreatemeshTab);
            app.Transparency01EditFieldLabel.HorizontalAlignment = 'right';
            app.Transparency01EditFieldLabel.Position = [8 212 105 22];
            app.Transparency01EditFieldLabel.Text = 'Transparency [0,1]';

            % Create Transparency01EditField
            app.Transparency01EditField = uieditfield(app.CreatemeshTab, 'numeric');
            app.Transparency01EditField.Limits = [0 1];
            app.Transparency01EditField.Position = [128 212 100 22];

            % Create ShowDropDownLabel
            app.ShowDropDownLabel = uilabel(app.CreatemeshTab);
            app.ShowDropDownLabel.HorizontalAlignment = 'right';
            app.ShowDropDownLabel.Position = [8 347 36 22];
            app.ShowDropDownLabel.Text = 'Show';

            % Create ShowDropDown
            app.ShowDropDown = uidropdown(app.CreatemeshTab);
            app.ShowDropDown.Items = {'Volume', 'Slice'};
            app.ShowDropDown.Position = [94 347 134 22];
            app.ShowDropDown.Value = 'Volume';

            % Create SavepngCheckBox
            app.SavepngCheckBox = uicheckbox(app.CreatemeshTab);
            app.SavepngCheckBox.Text = 'Save .png';
            app.SavepngCheckBox.Position = [8 160 76 22];
            app.SavepngCheckBox.Value = true;

            % Create SavefigCheckBox
            app.SavefigCheckBox = uicheckbox(app.CreatemeshTab);
            app.SavefigCheckBox.Text = 'Save .fig (Warning: it will create large files)';
            app.SavefigCheckBox.Position = [94 160 250 22];

            % Create Label_9
            app.Label_9 = uilabel(app.CreatemeshTab);
            app.Label_9.FontColor = [0.851 0.3255 0.098];
            app.Label_9.Position = [10 53 329 45];
            app.Label_9.Text = {'WARNING: for very large meshes, MATLAB may crash'; 'trying showing the figure. You may consider saving your'; 'meshes first to avoid losing your work.'};

            % Create GridCheckBox
            app.GridCheckBox = uicheckbox(app.CreatemeshTab);
            app.GridCheckBox.Text = 'Grid';
            app.GridCheckBox.Position = [8 185 45 22];
            app.GridCheckBox.Value = true;

            % Create AxislabelCheckBox
            app.AxislabelCheckBox = uicheckbox(app.CreatemeshTab);
            app.AxislabelCheckBox.Text = 'Axis label';
            app.AxislabelCheckBox.Position = [63 185 73 22];
            app.AxislabelCheckBox.Value = true;

            % Create Meshformat_UITable
            app.Meshformat_UITable = uitable(app.CreatemeshTab);
            app.Meshformat_UITable.ColumnName = {'Mesh format'; 'Choice'};
            app.Meshformat_UITable.ColumnWidth = {'1x'};
            app.Meshformat_UITable.RowName = {};
            app.Meshformat_UITable.ColumnEditable = [false true];
            app.Meshformat_UITable.Position = [385 202 329 167];

            % Create CreateMesh_text_2
            app.CreateMesh_text_2 = uilabel(app.CreatemeshTab);
            app.CreateMesh_text_2.HorizontalAlignment = 'center';
            app.CreateMesh_text_2.FontSize = 14;
            app.CreateMesh_text_2.FontWeight = 'bold';
            app.CreateMesh_text_2.Position = [384 379 331 22];
            app.CreateMesh_text_2.Text = 'Mesh saving options';

            % Create csvandmatLabel
            app.csvandmatLabel = uilabel(app.CreatemeshTab);
            app.csvandmatLabel.FontAngle = 'italic';
            app.csvandmatLabel.Position = [721 201 233 168];
            app.csvandmatLabel.Text = {'.mat and .csv save mesh in three files per'; ' mesh:'; '- a node file that contains a n x 3 arrays'; 'of the n vertices coordinates.'; '- a tetrahedron file that contains a'; 'm x 4 arrays of the m cells connectivities.'; '- a subomain file that contains a'; 'm x 4 arrrays of the m cells domain id.'; ''; 'You can later use these three files to'; 'recreate a mesh (see documentation'; 'for an example).'};

            % Create Face_CheckBox
            app.Face_CheckBox = uicheckbox(app.CreatemeshTab);
            app.Face_CheckBox.Text = 'Create face array (optional, and only required for visualization). Disable to save RAM for large geometry.';
            app.Face_CheckBox.Position = [348 427 587 22];
            app.Face_CheckBox.Value = true;

            % Create ColorisDropDownLabel
            app.ColorisDropDownLabel = uilabel(app.CreatemeshTab);
            app.ColorisDropDownLabel.HorizontalAlignment = 'right';
            app.ColorisDropDownLabel.Position = [8 266 46 22];
            app.ColorisDropDownLabel.Text = 'Color is';

            % Create ColorisDropDown
            app.ColorisDropDown = uidropdown(app.CreatemeshTab);
            app.ColorisDropDown.Items = {'z-axis', 'Phase Id'};
            app.ColorisDropDown.Position = [94 266 134 22];
            app.ColorisDropDown.Value = 'z-axis';

            % Create CalculatemeshqualityTab
            app.CalculatemeshqualityTab = uitab(app.TabGroup);
            app.CalculatemeshqualityTab.Title = 'Calculate mesh quality';
            app.CalculatemeshqualityTab.ForegroundColor = [0.4941 0.1843 0.5569];

            % Create Meshquality_Title
            app.Meshquality_Title = uilabel(app.CalculatemeshqualityTab);
            app.Meshquality_Title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Meshquality_Title.HorizontalAlignment = 'center';
            app.Meshquality_Title.FontWeight = 'bold';
            app.Meshquality_Title.Position = [10 634 944 22];
            app.Meshquality_Title.Text = 'Calculate mesh quality';

            % Create Label_11
            app.Label_11 = uilabel(app.CalculatemeshqualityTab);
            app.Label_11.Position = [10 592 962 28];
            app.Label_11.Text = {'Calculate the Joe-Liu mesh quality measure. For each cell, a value close to 1 represents a higher mesh quality (1 means equilateral tetrahedron) while a value close to 0 means'; 'nearly degenerated element. Reference: A. Liu, B. Joe, Relationship between tetrahedron shape measures, BIT 34 (2) (1994) 268-287.'};

            % Create Meshquality_JoeLiu_button
            app.Meshquality_JoeLiu_button = uibutton(app.CalculatemeshqualityTab, 'push');
            app.Meshquality_JoeLiu_button.ButtonPushedFcn = createCallbackFcn(app, @Meshquality_JoeLiu_buttonPushed, true);
            app.Meshquality_JoeLiu_button.BackgroundColor = [0.8 0.8 0.8];
            app.Meshquality_JoeLiu_button.Enable = 'off';
            app.Meshquality_JoeLiu_button.Position = [12 560 179 22];
            app.Meshquality_JoeLiu_button.Text = 'Calculate Joe-Liu mesh quality';

            % Create AboutTab
            app.AboutTab = uitab(app.TabGroup);
            app.AboutTab.Title = 'About';

            % Create About_title
            app.About_title = uilabel(app.AboutTab);
            app.About_title.BackgroundColor = [0.4667 0.6745 0.1882];
            app.About_title.HorizontalAlignment = 'center';
            app.About_title.FontWeight = 'bold';
            app.About_title.Position = [10 634 944 22];
            app.About_title.Text = 'About the meshing module';

            % Create About_TextArea
            app.About_TextArea = uitextarea(app.AboutTab);
            app.About_TextArea.Editable = 'off';
            app.About_TextArea.Position = [10 406 944 214];
            app.About_TextArea.Value = {'This module enables you to create volumetric tetrahedron-based mesh for a polycrystalline architecture, a single electrode, a half-cell or a full-cell.'; ''; 'The meshing can be performed through two different approaches:'; ''; '1) With Iso2mesh. Iso2mesh is using constrained Delaunay tetrahedralization (CGAL) for surface mesh extraction, Laplacian, Laplacian-HC, and Low-pass filters for surface mesh smoothing, and Tetgen for volumetric mesh generation and adaptative mesh resolution. Meshes are un-structured with a mesh density control. You will have to download Iso2mesh and add it to the Matlab path before using it.'; ''; '2) Alternatively, you can also create a structured grid mesh but without mesh density control. This approach is not relying on Iso2mesh.'; ''; 'The meshing can be performed:'; '- phase per phase (grain per grain for the polycrystalline architecture) with conforming meshes at the interfaces, i.e. one mesh per phase (per grain). Relevant for a segregated model'; '- for the whole geometry, i.e. with a unique mesh. Relevant for a monolithic model.'};

            % Create QuotationinstructionsTextAreaLabel
            app.QuotationinstructionsTextAreaLabel = uilabel(app.AboutTab);
            app.QuotationinstructionsTextAreaLabel.HorizontalAlignment = 'right';
            app.QuotationinstructionsTextAreaLabel.FontWeight = 'bold';
            app.QuotationinstructionsTextAreaLabel.Position = [12 217 134 22];
            app.QuotationinstructionsTextAreaLabel.Text = 'Quotation instructions';

            % Create About_Quotationinstructions
            app.About_Quotationinstructions = uitextarea(app.AboutTab);
            app.About_Quotationinstructions.Editable = 'off';
            app.About_Quotationinstructions.Position = [165 128 788 114];
            app.About_Quotationinstructions.Value = {'- For any mesh created with this module:'; 'F. L. E. Usseglio-Viretta et al., MATBOX: An Open-source Microstructure Analysis Toolbox for microstructure generation, segmentation, characterization, visualization, correlation, and meshing, SoftwareX, in preparation'; ''; '- If you create a mesh using Iso2mesh, please ALSO quote:'; 'Q. Fang and D. A. Boas, Tetrahedral Mesh Generation From Volumetric Binary and Gray-scale Images, Proceedings of IEEE International Symposium on Biomedical Imaging 2009, 2009, Pages 1142-1145'};

            % Create About_Logo_NREL
            app.About_Logo_NREL = uiimage(app.AboutTab);
            app.About_Logo_NREL.ImageClickedFcn = createCallbackFcn(app, @About_Logo_NRELImageClicked, true);
            app.About_Logo_NREL.Position = [10 16 264 100];
            app.About_Logo_NREL.ImageSource = 'logo_NREL.png';

            % Create OpendocumentationButton
            app.OpendocumentationButton = uibutton(app.AboutTab, 'push');
            app.OpendocumentationButton.ButtonPushedFcn = createCallbackFcn(app, @OpendocumentationButtonPushed, true);
            app.OpendocumentationButton.BackgroundColor = [0.8 0.8 0.8];
            app.OpendocumentationButton.Position = [807 78 147 38];
            app.OpendocumentationButton.Text = 'Open documentation';

            % Create Iso2MeshwebsiteButton
            app.Iso2MeshwebsiteButton = uibutton(app.AboutTab, 'push');
            app.Iso2MeshwebsiteButton.ButtonPushedFcn = createCallbackFcn(app, @Iso2MeshwebsiteButtonPushed, true);
            app.Iso2MeshwebsiteButton.BackgroundColor = [0.8 0.8 0.8];
            app.Iso2MeshwebsiteButton.FontSize = 20;
            app.Iso2MeshwebsiteButton.Position = [296 16 264 100];
            app.Iso2MeshwebsiteButton.Text = 'Iso2Mesh website';

            % Create Label_10
            app.Label_10 = uilabel(app.AboutTab);
            app.Label_10.FontSize = 11;
            app.Label_10.FontAngle = 'italic';
            app.Label_10.Position = [12 329 942 60];
            app.Label_10.Text = {'You can find examples of meshes generated with this module in the following articles and conferences:'; '- F. Usseglio-Viretta, A. Colclasure, K. Smith, J. Allen, J. Chang, P. Graf, S. Santhanagopalan, M. Keyser, D. Abraham, T. M. M. Heenan, Microstructure Modeling of Lithium-ion Graphite'; 'Electrodes Under Fast Charging, 36th International Battery Seminar, Fort Lauderdale, FL, USA, March 2019.'; '- J. Allen, J. Chang, F. Usseglio-Viretta, P. Graf, K. Smith, A Segregated Approach for Modeling the Electrochemistry in the 3-D Microstructure of Li-Ion Batteries and its Acceleration Using'; 'Block Preconditioners, submitted to Journal Scientific Computing, 86:42, 2021'};

            % Show the figure after all components are created
            app.MeshingmoduleUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = microstructure_meshing_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.MeshingmoduleUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MeshingmoduleUIFigure)
        end
    end
end