classdef Segmentation_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ROIfilteringandsegmentationmoduleUIFigure  matlab.ui.Figure
        VolumeMenu                      matlab.ui.container.Menu
        LoadvolumeMenu                  matlab.ui.container.Menu
        SavevolumeMenu                  matlab.ui.container.Menu
        DefaultlocationMenu             matlab.ui.container.Menu
        CustomlocationMenu              matlab.ui.container.Menu
        MenuCalculate                   matlab.ui.container.Menu
        PlothistogramMenu               matlab.ui.container.Menu
        forthewholevolumeMenu           matlab.ui.container.Menu
        histogram_standard_wholevolume_auto  matlab.ui.container.Menu
        histogram_standard_wholevolume_1binpervalue  matlab.ui.container.Menu
        histogram_distribution          matlab.ui.container.Menu
        foreachsliceMenu                matlab.ui.container.Menu
        PlotgreylevelMenu               matlab.ui.container.Menu
        asafunctionofpositionMenu       matlab.ui.container.Menu
        inaveraged2DmapsMenu            matlab.ui.container.Menu
        PlotnoiselevelMenu              matlab.ui.container.Menu
        foreachsliceMenu_2              matlab.ui.container.Menu
        HelpMenu                        matlab.ui.container.Menu
        DocumentationMenu               matlab.ui.container.Menu
        TabGroup                        matlab.ui.container.TabGroup
        InstructionsTab                 matlab.ui.container.Tab
        Loading_message                 matlab.ui.control.Label
        Fileloaded_text                 matlab.ui.control.Label
        Fileloaded_min_max              matlab.ui.control.Label
        Fileloaded_memory_datatype      matlab.ui.control.Label
        Fileloaded_numberofvoxel        matlab.ui.control.Label
        Fileloaded_dimension            matlab.ui.control.Label
        Fileloaded_Filename             matlab.ui.control.Label
        Fileloaded_Path                 matlab.ui.control.Label
        VolumeLoadedStatus              matlab.ui.control.Label
        Instructions7                   matlab.ui.control.Label
        Instructions6                   matlab.ui.control.Label
        Instructions5                   matlab.ui.control.Label
        Instructions4                   matlab.ui.control.Label
        Instructions3                   matlab.ui.control.Label
        Instructions2                   matlab.ui.control.Label
        Instructions_Instructions       matlab.ui.control.Label
        Instructions1                   matlab.ui.control.Label
        Logo                            matlab.ui.control.Image
        SaveOptions                     matlab.ui.container.Tab
        Save_choiceFont_DropDown        matlab.ui.control.DropDown
        FontnameusedforfiguresLabel     matlab.ui.control.Label
        Save_checkbox                   matlab.ui.control.CheckBox
        FolderandfilenameusedifyouselectMenuSavevolumeLabel  matlab.ui.control.Label
        DefaultfilenamenoneLabel        matlab.ui.control.Label
        DefaultsavefoldernoneselectedLabel  matlab.ui.control.Label
        ClicktoselectsavefolderButton   matlab.ui.control.Button
        Save_instructions               matlab.ui.control.Label
        RegionofInterestViewTab         matlab.ui.container.Tab
        ROI_background_OptionsDropDown  matlab.ui.control.DropDown
        AutocropforLabel                matlab.ui.control.Label
        SetidentifiedbackgroundvaluetoEditField  matlab.ui.control.NumericEditField
        SetidentifiedbackgroundvaluetoEditFieldLabel  matlab.ui.control.Label
        ROI_background_SaveButton       matlab.ui.control.Button
        ROI_background_UndoButton       matlab.ui.control.Button
        ROI_background_DoButton         matlab.ui.control.Button
        ROI_backgroundetection_DropDown  matlab.ui.control.DropDown
        Label                           matlab.ui.control.Label
        ROI_ModifysomeslicesalongLabel_2  matlab.ui.control.Label
        Savecurrent2DviewButton         matlab.ui.control.Button
        VieworthogonalslicesButton      matlab.ui.control.Button
        StdthresholdEditField           matlab.ui.control.NumericEditField
        StdthresholdEditFieldLabel      matlab.ui.control.Label
        MedianfilterrangeEditField      matlab.ui.control.NumericEditField
        MedianfilterrangeEditFieldLabel  matlab.ui.control.Label
        BackgroundvalueEditField        matlab.ui.control.NumericEditField
        BackgroundvalueEditFieldLabel   matlab.ui.control.Label
        CropbackgroundButton            matlab.ui.control.Button
        ROI_modify_SaveButton           matlab.ui.control.Button
        ROI_modify_UndoButton           matlab.ui.control.Button
        ROI_modify_DoButton             matlab.ui.control.Button
        ROI_modify_ActionDropDown       matlab.ui.control.DropDown
        ActionforSlicesfromAtoBDropDownLabel  matlab.ui.control.Label
        ROI_modify_TosliceB             matlab.ui.control.NumericEditField
        TosliceBEditFieldLabel          matlab.ui.control.Label
        ROI_modify_FromsliceA           matlab.ui.control.NumericEditField
        FromsliceAEditFieldLabel        matlab.ui.control.Label
        ROI_modify_axeAxeDropDown       matlab.ui.control.DropDown
        ROI_ModifysomeslicesalongLabel  matlab.ui.control.Label
        ROI_flipswap_SaveButton         matlab.ui.control.Button
        ROI_flipswap_UndoButton         matlab.ui.control.Button
        ROI_flipswap_DoButton           matlab.ui.control.Button
        FliporswapaxisLabel             matlab.ui.control.Label
        ROI_flipswapDropDown            matlab.ui.control.DropDown
        ROI_rotation_SaveButton         matlab.ui.control.Button
        ROI_rotation_UndoButton         matlab.ui.control.Button
        ROI_rotation_DoButton           matlab.ui.control.Button
        ROI_angle_choice_DropDown       matlab.ui.control.DropDown
        AngleDropDownLabel              matlab.ui.control.Label
        ROI_angle_spinner               matlab.ui.control.Spinner
        ROI_rotation_text               matlab.ui.control.Label
        ROI_voxelsize                   matlab.ui.control.NumericEditField
        VoxelsizeoptionalLabel          matlab.ui.control.Label
        ROI_crop_SaveButton             matlab.ui.control.Button
        ROI_crop_UndoButton             matlab.ui.control.Button
        ROI_crop_DoButton               matlab.ui.control.Button
        ROI_crop_text                   matlab.ui.control.Label
        ROI_table_crop                  matlab.ui.control.Table
        ROI_SliceSpinner                matlab.ui.control.Spinner
        ROI_View3D                      matlab.ui.control.Button
        ROI_checkbox_updateslicer       matlab.ui.control.CheckBox
        ROI_instructions_tab            matlab.ui.control.Label
        ROI_colormap                    matlab.ui.control.DropDown
        ColormapLabel                   matlab.ui.control.Label
        ROI_SliceSlider                 matlab.ui.control.Slider
        SliceselectionLabel             matlab.ui.control.Label
        ROI_ViewnormaltoDropDown        matlab.ui.control.DropDown
        ViewnormaltoDropDown_2Label     matlab.ui.control.Label
        ROI_2Dslice                     matlab.ui.control.UIAxes
        FormatconversionTab             matlab.ui.container.Tab
        Convert_UpdatetableButton       matlab.ui.control.Button
        Convert_Utilization_checkbox    matlab.ui.control.CheckBox
        Convert_table                   matlab.ui.control.Table
        Convert_SaveButton              matlab.ui.control.Button
        Convert_UndoButton              matlab.ui.control.Button
        Convert_DoButton                matlab.ui.control.Button
        Convert_instructions            matlab.ui.control.Label
        Convert_Choice_DropDown         matlab.ui.control.DropDown
        ConverttoLabel                  matlab.ui.control.Label
        TovisualizeimpactonhistogramCalculateHistogramLabel  matlab.ui.control.Label
        UpdownscalingTab                matlab.ui.container.Tab
        Scale_BackgroundvalueEditField  matlab.ui.control.NumericEditField
        BackgroundvalueEditField_2Label  matlab.ui.control.Label
        Scale_SaveButton                matlab.ui.control.Button
        Scale_UndoButton                matlab.ui.control.Button
        Scale_DoButton                  matlab.ui.control.Button
        Scale_instructions              matlab.ui.control.Label
        Scale_SliceselectionSlider      matlab.ui.control.Slider
        SliceselectionSliderLabel       matlab.ui.control.Label
        Scale_DatatypeButtonGroup       matlab.ui.container.ButtonGroup
        LabelButton                     matlab.ui.control.ToggleButton
        GreylevelButton                 matlab.ui.control.ToggleButton
        Scale_text_4                    matlab.ui.control.Label
        Scale_newvoxelsizeunitLabel     matlab.ui.control.Label
        Scale_NewvoxelsizeEditField     matlab.ui.control.NumericEditField
        NewvoxelsizeEditFieldLabel      matlab.ui.control.Label
        Scale_table                     matlab.ui.control.Table
        Scale_VoxelunitDropDown         matlab.ui.control.DropDown
        VoxelunitDropDownLabel          matlab.ui.control.Label
        Scale_InitialvoxelsizeEditField  matlab.ui.control.NumericEditField
        InitialvoxelsizeEditFieldLabel  matlab.ui.control.Label
        Scale_text_1                    matlab.ui.control.Label
        Scale_scalingfactorEditField    matlab.ui.control.NumericEditField
        ScalingfactorEditFieldLabel     matlab.ui.control.Label
        Scale_text_3                    matlab.ui.control.Label
        Scale_text_2                    matlab.ui.control.Label
        Scaling_UIAxes_after            matlab.ui.control.UIAxes
        Scaling_UIAxes_before           matlab.ui.control.UIAxes
        ImagequalityTab                 matlab.ui.container.Tab
        Quality_CalculateSeparabilityButton  matlab.ui.control.Button
        Quality_Text1                   matlab.ui.control.Label
        Quality_numberofphase           matlab.ui.control.NumericEditField
        NumberofphaseincludingporebackgroundLabel  matlab.ui.control.Label
        Quality_instructions            matlab.ui.control.Label
        ContrastcorrectionTab           matlab.ui.container.Tab
        Contrast_Method_DropDown        matlab.ui.control.DropDown
        OperationLabel                  matlab.ui.control.Label
        Contrast_Advanced_SaveButton    matlab.ui.control.Button
        Contrast_Advanced_UndoButton    matlab.ui.control.Button
        Contrast_Advanced_DoButton      matlab.ui.control.Button
        Contrast_Text8                  matlab.ui.control.Label
        Contrast_Text7                  matlab.ui.control.Label
        Contrast_TableParameters        matlab.ui.control.Table
        Contrast_Text6                  matlab.ui.control.Label
        Contrast_Text4                  matlab.ui.control.Label
        Contrast_Text5                  matlab.ui.control.Label
        Contrast_TableWhere             matlab.ui.control.Table
        Contrast_Text3                  matlab.ui.control.Label
        Contrast_SliceselectionSlider   matlab.ui.control.Slider
        SliceselectionSliderLabel_2     matlab.ui.control.Label
        Contrast_Customrange_SaveButton  matlab.ui.control.Button
        Contrast_Customrange_UndoButton  matlab.ui.control.Button
        Contrast_Customrange_DoButton   matlab.ui.control.Button
        Contrast_Text2                  matlab.ui.control.Label
        Contrast_Burn_SaveButton        matlab.ui.control.Button
        Contrast_Burn_UndoButton        matlab.ui.control.Button
        Contrast_Burn_DoButton          matlab.ui.control.Button
        Contrast_instructions           matlab.ui.control.Label
        Contrast_Text1                  matlab.ui.control.Label
        Contrast_CustomrangerescaleusingthisnumberofvalueEditField  matlab.ui.control.NumericEditField
        RescaleimageusingthisnumberofvalueLabel  matlab.ui.control.Label
        Contrast_VolumepercentthresholdforthelowvaluesEditField  matlab.ui.control.NumericEditField
        VolumepercentthresholdforthelowvaluesEditFieldLabel  matlab.ui.control.Label
        Contrast_VolumepercentthresholdforthehighvaluesEditField  matlab.ui.control.NumericEditField
        VolumepercentthresholdforthehighvaluesEditFieldLabel  matlab.ui.control.Label
        Contrast_UIAxes_after           matlab.ui.control.UIAxes
        Contrast_UIAxes_before          matlab.ui.control.UIAxes
        FilteringTab                    matlab.ui.container.Tab
        Filtering_Anisotropy_EstimateCheckBox  matlab.ui.control.CheckBox
        Filtering_ProgressnotrunningLabel  matlab.ui.control.Label
        Filtering_previewCheckBox       matlab.ui.control.CheckBox
        Filtering_Text4                 matlab.ui.control.Label
        Filtering_NLMF_SaveButton       matlab.ui.control.Button
        Filtering_NLMF_UndoButton       matlab.ui.control.Button
        Filtering_NLMF_DoButton         matlab.ui.control.Button
        Filetring_NLMF_ComparisonwindowsizeEditField  matlab.ui.control.NumericEditField
        ComparisonwindowsizeEditFieldLabel  matlab.ui.control.Label
        Filetring_NLMF_SearchwindowsizeEditField  matlab.ui.control.NumericEditField
        SearchwindowsizeEditFieldLabel  matlab.ui.control.Label
        Filetring_NLMF_DegreeofsmoothingEditField  matlab.ui.control.NumericEditField
        DegreeofsmoothingEditFieldLabel  matlab.ui.control.Label
        Filetring_NLMF_EstimateCheckBox  matlab.ui.control.CheckBox
        Filtering_Text2                 matlab.ui.control.Label
        Filtering_SliceselectionSlider  matlab.ui.control.Slider
        SliceselectionSliderLabel_3     matlab.ui.control.Label
        Filtering_Text3                 matlab.ui.control.Label
        Filtering_Anisotropic_SaveButton  matlab.ui.control.Button
        Filtering_Anisotropic_UndoButton  matlab.ui.control.Button
        Filtering_Anisotropic_DoButton  matlab.ui.control.Button
        Filtering_AnisotropicGradient_Table  matlab.ui.control.Table
        Filtering_AnisotropicFilter_IterEditField  matlab.ui.control.NumericEditField
        NumberofiterationEditFieldLabel  matlab.ui.control.Label
        Filtering_AnisotropicFilter_ConnDropDown  matlab.ui.control.DropDown
        ConnectivityDropDownLabel       matlab.ui.control.Label
        Filtering_AnisotropicFilter_MethodDropDown  matlab.ui.control.DropDown
        ConductionmethodDropDownLabel   matlab.ui.control.Label
        Filtering_AnisotropicFilter_2D3DDropDown  matlab.ui.control.DropDown
        ApplyDropDownLabel              matlab.ui.control.Label
        Filtering_Text1                 matlab.ui.control.Label
        Filtering_instructions          matlab.ui.control.Label
        Filetring_UIAxes_after          matlab.ui.control.UIAxes
        Filtering_UIAxes_before         matlab.ui.control.UIAxes
        SegmentationTab                 matlab.ui.container.Tab
        Segmentation_MovingaveragefilterrangeEditField  matlab.ui.control.NumericEditField
        MovingaveragefilterrangeEditFieldLabel  matlab.ui.control.Label
        Segmentation_SmoothSliceThreshold_CheckBox  matlab.ui.control.CheckBox
        Segmentation_OpacityoverlayEditField  matlab.ui.control.NumericEditField
        OpacityoverlayEditFieldLabel    matlab.ui.control.Label
        Segmentation_Highthreshold_SliceSpinner  matlab.ui.control.Spinner
        Segmentation_Lowthreshold_SliceSpinner  matlab.ui.control.Spinner
        Segmentation_update_range       matlab.ui.control.Button
        Segmentation_Highthreshold_SliceSlider  matlab.ui.control.Slider
        ThresholdhighLabel              matlab.ui.control.Label
        Segmentation_Lowthreshold_SliceSlider  matlab.ui.control.Slider
        ThresholdlowLabel               matlab.ui.control.Label
        Segmentation_SliceSlider        matlab.ui.control.Slider
        SliceselectionLabel_2           matlab.ui.control.Label
        Segmentation_Local_SaveButton   matlab.ui.control.Button
        Segmentation_Local_UndoButton   matlab.ui.control.Button
        Segmentation_Local_DoButton     matlab.ui.control.Button
        Segmentation_Local_MethodDropDown  matlab.ui.control.DropDown
        MethodDropDown_2Label           matlab.ui.control.Label
        Segmentation_Text2              matlab.ui.control.Label
        Segmentation_Global_SaveButton  matlab.ui.control.Button
        Segmentation_Global_UndoButton  matlab.ui.control.Button
        Segmentation_Global_DoButton    matlab.ui.control.Button
        Segmentation_CompareOtsuandexpectedvolumefractionsButton  matlab.ui.control.Button
        Segmentation_Global_UITable     matlab.ui.control.Table
        Segmentation_NumberofphaseEditField  matlab.ui.control.NumericEditField
        NumberofphaseEditFieldLabel     matlab.ui.control.Label
        Segmentation_Phase_UITable      matlab.ui.control.Table
        Segmentation_VoxelsizemicrometersEditField  matlab.ui.control.NumericEditField
        VoxelsizemicrometersEditFieldLabel  matlab.ui.control.Label
        Segmentation_Global_MethodDropDown  matlab.ui.control.DropDown
        MethodDropDownLabel             matlab.ui.control.Label
        Segmentation_Text1              matlab.ui.control.Label
        Segmentation_instructions       matlab.ui.control.Label
        Segmentation_UIAxes             matlab.ui.control.UIAxes
        SegmentationsensitivityTab      matlab.ui.container.Tab
        Sensitivity_ProgressnotrunningLabel  matlab.ui.control.Label
        Sensitivity_CalculateparametersensitivityButton  matlab.ui.control.Button
        Sensitivity_Text_2              matlab.ui.control.Label
        Sensitivity_Text_1              matlab.ui.control.Label
        Sensitivity_instructions        matlab.ui.control.Label
        Sensitivity_VoxelsizenanometersEditField  matlab.ui.control.NumericEditField
        VoxelsizemicrometersLabel       matlab.ui.control.Label
        Sensitivity_MaximumnumberofthresholdevaluatedEditField  matlab.ui.control.NumericEditField
        MaximumnumberofthresholdevaluatedEditFieldLabel  matlab.ui.control.Label
        Sensitivity_ExpectedporosityEditField  matlab.ui.control.NumericEditField
        ExpectedporosityEditFieldLabel  matlab.ui.control.Label
        Sensitivity_HigherboundEditField  matlab.ui.control.NumericEditField
        HigherboundEditFieldLabel       matlab.ui.control.Label
        Sensitivity_LowerboundEditField  matlab.ui.control.NumericEditField
        LowerboundEditFieldLabel        matlab.ui.control.Label
        Sensitivity_PoretortuosityCheckBox  matlab.ui.control.CheckBox
        Sensitivity_SpecificsurfaceareaCheckBox  matlab.ui.control.CheckBox
        Sensitivity_ParticlediameterCheckBox  matlab.ui.control.CheckBox
        Sensitivity_PorosityCheckBox    matlab.ui.control.CheckBox
        PhasereassignmentTab            matlab.ui.container.Tab
        PhaseAssignment_UpdatetableButton  matlab.ui.control.Button
        PhaseAssignment_SaveButton      matlab.ui.control.Button
        PhaseAssignment_UndoButton      matlab.ui.control.Button
        PhaseAssignment_DoButton        matlab.ui.control.Button
        PhaseAssignment_UITable         matlab.ui.control.Table
        PhaseAssignment_Instructions    matlab.ui.control.Label
        HistoryTab                      matlab.ui.container.Tab
        UITable_History                 matlab.ui.control.Table
        Instructions_History            matlab.ui.control.Label
    end

    
    properties (Access = public)
        % Set global variable
        Microstructure=[]; % Microstructure array, initialized empty
        Microstructure_undo=[]; % Microstructure array, initialized empty
        Do_parameters=[]; % Operation parameters
        Do_elapsedtime=[]; % Operation time
        Filename=[]; % Loaded file
        Savefolder=[]; % Default save folder
        defaultfilesave=[]; % Default file savename
        Dimension=[]; % Microstructure dimension, initialized empty
        History_step=[]; History_operations=[]; History_parameters=[]; History_elapsedtime=[]; % History log, initialized empty
        Operation_step=0;
        fontisze_oneplot_axe = 12;
        fontisze_oneplot_title = 14;
        fontisze_subplot4_axe = 12;
        fontisze_subplot4_title = 14;       
        sav_thresholds = []; sav_thresholds_smoothed=[]; sav_translated_thresholds=[];
    end
    
    
    methods (Access = private)
        
        function [] = update_history_log(app)
            History_table = table(app.History_step,app.History_operations,app.History_parameters,app.History_elapsedtime);
            History_table.Properties.VariableNames = {'Step','Operations','Parameters','Elapsed time'};
            app.UITable_History.Data=History_table;
            [~,filesave,~] = fileparts(app.Filename);
            app.defaultfilesave= [filesave '_step_' num2str(app.Operation_step) '.tif'];
            app.DefaultfilenamenoneLabel.Text = ['Current file name: ' app.defaultfilesave];
        end
        
        function [] = update_GUI_TabInstruction_status(app,status_GUI)
            % Instructions tab
            app.Fileloaded_Path.Visible =status_GUI;
            app.Fileloaded_Filename.Visible =status_GUI;
            app.Fileloaded_dimension.Visible =status_GUI;
            app.Fileloaded_numberofvoxel.Visible =status_GUI;
            app.Fileloaded_memory_datatype.Visible =status_GUI;
            app.Fileloaded_min_max.Visible =status_GUI;
            app.Fileloaded_text.Visible = status_GUI;
        end
        
        function [] = Turn_on_off(app,status_GUI,tab)
            child_handles = allchild(tab); % Get all children
            [n_children,~]=size(child_handles);
            for k=1:1:n_children
                child_handles(k).Visible=status_GUI;
                if isprop(child_handles(k),'Enable')
                    child_handles(k).Enable=status_GUI;
                end
            end
            app.MenuCalculate.Enable=status_GUI;
            
            if strcmp(tab.Title,'Region of Interest / View')
                if strcmp(status_GUI,'off')
                    cla(app.ROI_2Dslice);
                    colorbar(app.ROI_2Dslice,'off')
                else
                    app.ROI_crop_UndoButton.Enable='off'; app.ROI_crop_SaveButton.Enable='off';
                    app.ROI_rotation_UndoButton.Enable='off'; app.ROI_rotation_SaveButton.Enable='off';
                    app.ROI_modify_UndoButton.Enable='off'; app.ROI_modify_SaveButton.Enable='off';
                    app.ROI_flipswap_UndoButton.Enable='off'; app.ROI_flipswap_SaveButton.Enable='off';
                    app.ROI_background_UndoButton.Enable='off'; app.ROI_background_SaveButton.Enable='off';
                end
            elseif strcmp(tab.Title,'Format conversion')
                if strcmp(status_GUI,'on')
                    app.Convert_UndoButton.Enable='off'; app.Convert_SaveButton.Enable='off';
                end
            elseif strcmp(tab.Title,'Up/down scaling')
                if strcmp(status_GUI,'off')
                    cla(app.Scaling_UIAxes_before); cla(app.Scaling_UIAxes_after);
                    colorbar(app.Scaling_UIAxes_before,'off'); colorbar(app.Scaling_UIAxes_after,'off')
                else
                    app.Scale_UndoButton.Enable='off'; app.Scale_SaveButton.Enable='off';
                    app.Scale_SliceselectionSlider.Enable='off';
                    app.Scale_BackgroundvalueEditField.Enable='off';
                end
                
            elseif strcmp(tab.Title,'Phase re-assignment')
                if strcmp(status_GUI,'on')
                    app.PhaseAssignment_UndoButton.Enable='off'; app.PhaseAssignment_SaveButton.Enable='off';
                end   
            elseif strcmp(tab.Title,'Contrast correction')
                if strcmp(status_GUI,'off')
                    cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after);
                    colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');
                else
                    app.Contrast_Burn_UndoButton.Enable='off'; app.Contrast_Burn_SaveButton.Enable='off';
                    app.Contrast_Customrange_UndoButton.Enable='off'; app.Contrast_Customrange_SaveButton.Enable='off';
                    app.Contrast_Advanced_UndoButton.Enable='off'; app.Contrast_Advanced_SaveButton.Enable='off';
                    app.Contrast_SliceselectionSlider.Enable='off';
                end                 
            elseif strcmp(tab.Title,'Segmentation (sensitivity)')
                if strcmp(status_GUI,'on')
                    app.Sensitivity_PorosityCheckBox.Enable='off';
                end
            elseif strcmp(tab.Title,'Filtering')
                if strcmp(status_GUI,'off')
                    cla(app.Filtering_UIAxes_before); cla(app.Filetring_UIAxes_after);
                    colorbar(app.Filtering_UIAxes_before,'off'); colorbar(app.Filetring_UIAxes_after,'off');
                else
                    app.Filtering_Anisotropic_UndoButton.Enable='off'; app.Filtering_Anisotropic_SaveButton.Enable='off';
                    app.Filtering_NLMF_UndoButton.Enable='off'; app.Filtering_NLMF_SaveButton.Enable='off';
                    app.Filtering_SliceselectionSlider.Enable='off';
                    app.Filetring_NLMF_DegreeofsmoothingEditField.Enable='off';
                    app.Filtering_AnisotropicFilter_IterEditField.Enable='off'; app.Filtering_AnisotropicGradient_Table.Enable='off';
                end
            elseif strcmp(tab.Title,'Segmentation')
                if strcmp(status_GUI,'off')
                    cla(app.Segmentation_UIAxes);
                    colorbar(app.Segmentation_UIAxes,'off');
                else
                    app.Segmentation_Local_UndoButton.Enable='off'; app.Segmentation_Local_SaveButton.Enable='off';
                    app.Segmentation_Global_UndoButton.Enable='off'; app.Segmentation_Global_SaveButton.Enable='off';
                end
            end
        end
        
        function [] = DoUndoSave_button(app,tab,operation,action)
            if strcmp(action,'undo') || strcmp(action,'save')
                % % All Do button 'on', All Undo and Save button 'off'
                % Menu calculate
                app.MenuCalculate.Enable='on';
                % ROI tab
                app.ROI_crop_DoButton.Enable='on'; app.ROI_crop_UndoButton.Enable='off'; app.ROI_crop_SaveButton.Enable='off';
                app.ROI_rotation_DoButton.Enable='on'; app.ROI_rotation_UndoButton.Enable='off'; app.ROI_rotation_SaveButton.Enable='off';
                app.ROI_modify_DoButton.Enable='on'; app.ROI_modify_UndoButton.Enable='off'; app.ROI_modify_SaveButton.Enable='off';
                app.ROI_flipswap_DoButton.Enable='on'; app.ROI_flipswap_UndoButton.Enable='off'; app.ROI_flipswap_SaveButton.Enable='off';
                app.ROI_background_DoButton.Enable='on'; app.ROI_background_UndoButton.Enable='off'; app.ROI_background_SaveButton.Enable='off';
                % Convert tab
                app.Convert_DoButton.Enable='on'; app.Convert_UndoButton.Enable='off'; app.Convert_SaveButton.Enable='off';
                % Scale tab
                app.Scale_DoButton.Enable='on'; app.Scale_UndoButton.Enable='off'; app.Scale_SaveButton.Enable='off';
                % Reassign tab
                app.PhaseAssignment_DoButton.Enable='on'; app.PhaseAssignment_UndoButton.Enable='off'; app.PhaseAssignment_SaveButton.Enable='off';
                % Contrast tab
                app.Contrast_Burn_DoButton.Enable='on'; app.Contrast_Burn_UndoButton.Enable='off'; app.Contrast_Burn_SaveButton.Enable='off';
                app.Contrast_Customrange_DoButton.Enable='on'; app.Contrast_Customrange_UndoButton.Enable='off'; app.Contrast_Customrange_SaveButton.Enable='off';
                app.Contrast_Advanced_DoButton.Enable='on'; app.Contrast_Advanced_UndoButton.Enable='off'; app.Contrast_Advanced_SaveButton.Enable='off';
                % Sensitivity tab
                app.Sensitivity_CalculateparametersensitivityButton.Enable='on';
                % Filtering tab
                app.Filtering_Anisotropic_DoButton.Enable='on'; app.Filtering_Anisotropic_UndoButton.Enable='off'; app.Filtering_Anisotropic_SaveButton.Enable='off';
                app.Filtering_NLMF_DoButton.Enable='on'; app.Filtering_NLMF_UndoButton.Enable='off'; app.Filtering_NLMF_SaveButton.Enable='off';
                % Segmentation tab
                app.Segmentation_Global_DoButton.Enable='on'; app.Segmentation_Global_UndoButton.Enable='off'; app.Segmentation_Global_SaveButton.Enable='off';
                app.Segmentation_Local_DoButton.Enable='on'; app.Segmentation_Local_UndoButton.Enable='off'; app.Segmentation_Local_SaveButton.Enable='off';
                                                                
            elseif strcmp(action,'do')
                % % All Do button 'off'
                % Menu calculate
                app.MenuCalculate.Enable='off';                
                % ROI tab
                app.ROI_crop_DoButton.Enable='off'; app.ROI_rotation_DoButton.Enable='off'; app.ROI_modify_DoButton.Enable='off'; app.ROI_flipswap_DoButton.Enable='off'; app.ROI_background_DoButton.Enable='off';
                % Convert tab
                app.Convert_DoButton.Enable='off';
                % Scale tab
                app.Scale_DoButton.Enable='off';
                % Reassign tab
                app.PhaseAssignment_DoButton.Enable='off';  
                % Contrast tab
                app.Contrast_Burn_DoButton.Enable='off';  app.Contrast_Customrange_DoButton.Enable='off';  app.Contrast_Advanced_DoButton.Enable='off';  
                % Sensitivity tab
                app.Sensitivity_CalculateparametersensitivityButton.Enable='off';
                % Filtering tab
                app.Filtering_Anisotropic_DoButton.Enable='off';  app.Filtering_NLMF_DoButton.Enable='off';
                % Segmentation tab
                app.Segmentation_Global_DoButton.Enable='off';  app.Segmentation_Local_DoButton.Enable='off';
                
                % % Undo and Save button specific to this operation are 'on'
                if strcmp(tab,'ROI')
                    if strcmp(operation,'crop')
                        app.ROI_crop_UndoButton.Enable='on'; app.ROI_crop_SaveButton.Enable='on';
                    elseif strcmp(operation,'rotation')
                        app.ROI_rotation_UndoButton.Enable='on'; app.ROI_rotation_SaveButton.Enable='on';
                    elseif strcmp(operation,'modify')
                        app.ROI_modify_UndoButton.Enable='on'; app.ROI_modify_SaveButton.Enable='on';
                    elseif strcmp(operation,'flipswap')
                        app.ROI_flipswap_UndoButton.Enable='on'; app.ROI_flipswap_SaveButton.Enable='on';
                    elseif strcmp(operation,'background')
                        app.ROI_background_UndoButton.Enable='on'; app.ROI_background_SaveButton.Enable='on';
                    end
                elseif strcmp(tab,'Convert')
                    app.Convert_UndoButton.Enable='on'; app.Convert_SaveButton.Enable='on';
                elseif strcmp(tab,'Background')
                    app.ROI_background_UndoButton.Enable='on'; app.ROI_background_SaveButton.Enable='on';
                elseif strcmp(tab,'Scale')
                    app.Scale_UndoButton.Enable='on'; app.Scale_SaveButton.Enable='on';
                elseif strcmp(tab,'Reassign')
                    app.PhaseAssignment_UndoButton.Enable='on'; app.PhaseAssignment_SaveButton.Enable='on';                    
                elseif strcmp(tab,'Burn')
                    app.Contrast_Burn_UndoButton.Enable='on'; app.Contrast_Burn_SaveButton.Enable='on';                 
                elseif strcmp(tab,'CustomRange')
                    app.Contrast_Customrange_UndoButton.Enable='on'; app.Contrast_Customrange_SaveButton.Enable='on';    
                elseif strcmp(tab,'AdvContrast')
                    app.Contrast_Advanced_UndoButton.Enable='on'; app.Contrast_Advanced_SaveButton.Enable='on';   
                elseif strcmp(tab,'NLMF')
                    app.Filtering_NLMF_UndoButton.Enable='on'; app.Filtering_NLMF_SaveButton.Enable='on';    
                elseif strcmp(tab,'Anisotropicfilter')
                    app.Filtering_Anisotropic_UndoButton.Enable='on'; app.Filtering_Anisotropic_SaveButton.Enable='on'; 
                elseif strcmp(tab,'SegmentationGlobal')
                    app.Segmentation_Global_UndoButton.Enable='on'; app.Segmentation_Global_SaveButton.Enable='on';    
                elseif strcmp(tab,'SegmentationLocal')
                    app.Segmentation_Local_UndoButton.Enable='on'; app.Segmentation_Local_SaveButton.Enable='on'; 
                end                
                
            end
        end
        
        function [] = update_GUI_dimension(app)
            domain_size = size(app.Microstructure);
            app.Dimension = length(domain_size);
            if app.Dimension==2
                statut_GUI_3D='off';
                % Region of Interest tab
                app.ROI_angle_choice_DropDown.Items = {'Normal to axe 3'};
                app.ROI_angle_choice_DropDown.Value = {'Normal to axe 3'};
                app.ROI_flipswapDropDown.Items = {'Flip axis 1','Flip axis 2','Swap axis 1 with axis 2'};
                app.ROI_flipswapDropDown.Value = {'Flip axis 1'};
                app.ROI_ViewnormaltoDropDown.Items = {'Axe 3'};
                app.ROI_ViewnormaltoDropDown.Value = 'Axe 3';
                xlabel(app.ROI_2Dslice, '2nd Axis');
                ylabel(app.ROI_2Dslice, '1st Axis');
                app.ROI_modify_axeAxeDropDown.Items = {'Axe 1', 'Axe 2'};
                app.ROI_modify_axeAxeDropDown.Value = {'Axe 1'};
                app.ROI_modify_FromsliceA.Limits = [1 domain_size(1)]; app.ROI_modify_TosliceB.Limits = [1 domain_size(1)];
            else
                statut_GUI_3D='on';
                % Region of Interest tab
                app.ROI_angle_choice_DropDown.Items = {'Normal to axe 1', 'Normal to axe 2', 'Normal to axe 3'};
                app.ROI_flipswapDropDown.Items = {'Flip axis 1','Flip axis 2','Flip axis 3','Swap axis 1 with axis 2','Swap axis 1 with axis 3','Swap axis 2 with axis 3'};
                app.ROI_flipswapDropDown.Value = {'Flip axis 1'};
                app.ROI_ViewnormaltoDropDown.Items = {'Axe 1', 'Axe 2', 'Axe 3'};
                if isempty(app.ROI_ViewnormaltoDropDown.UserData)
                    app.ROI_ViewnormaltoDropDown.UserData = 3;
                    app.ROI_ViewnormaltoDropDown.Value = 'Axe 3';
                end
                if app.ROI_ViewnormaltoDropDown.UserData==1
                    str_x = '3rd Axis'; str_y = '2nd Axis';
                elseif app.ROI_ViewnormaltoDropDown.UserData==2
                    str_x = '3rd Axis'; str_y = '1st Axis';
                elseif app.ROI_ViewnormaltoDropDown.UserData==3
                    str_x = '2nd Axis'; str_y = '1st Axis';
                end
                xlabel(app.ROI_2Dslice, str_x);
                ylabel(app.ROI_2Dslice, str_y);
                app.ROI_SliceSlider.Limits = [1 domain_size(app.ROI_ViewnormaltoDropDown.UserData)];
                app.ROI_SliceSlider.Value = round(domain_size(app.ROI_ViewnormaltoDropDown.UserData)/2);
                app.ROI_SliceSpinner.Limits = app.ROI_SliceSlider.Limits;
                app.ROI_SliceSpinner.Value = app.ROI_SliceSlider.Value;
                app.ROI_SliceSlider.MajorTicks = round(linspace(1,domain_size(app.ROI_ViewnormaltoDropDown.UserData),5));
                app.ROI_SliceSlider.MinorTicks = round(linspace(1,domain_size(app.ROI_ViewnormaltoDropDown.UserData),9));
                app.ROI_modify_axeAxeDropDown.Items = {'Axe 1', 'Axe 2', 'Axe 3'};
                app.ROI_modify_axeAxeDropDown.Value = {'Axe 1'};
                app.ROI_modify_FromsliceA.Limits = [1 domain_size(1)]; app.ROI_modify_TosliceB.Limits = [1 domain_size(1)];
                app.Scale_SliceselectionSlider.Limits = [1 domain_size(3)];
                app.Scale_SliceselectionSlider.Value = round(domain_size(3)/2);
                app.Scale_SliceselectionSlider.MajorTicks = round(linspace(1,domain_size(3),10));
                app.Scale_SliceselectionSlider.MinorTicks = round(linspace(1,domain_size(3),19));                 
                app.Contrast_SliceselectionSlider.Limits = [1 domain_size(3)];
                app.Contrast_SliceselectionSlider.Value = round(domain_size(3)/2);
                app.Contrast_SliceselectionSlider.MajorTicks = round(linspace(1,domain_size(3),10));
                app.Contrast_SliceselectionSlider.MinorTicks = round(linspace(1,domain_size(3),19));    
                app.Filtering_SliceselectionSlider.Limits = [1 domain_size(3)];
                app.Filtering_SliceselectionSlider.Value = round(domain_size(3)/2);
                app.Filtering_SliceselectionSlider.MajorTicks = round(linspace(1,domain_size(3),10));
                app.Filtering_SliceselectionSlider.MinorTicks = round(linspace(1,domain_size(3),19));  
                app.Segmentation_SliceSlider.Limits = [1 domain_size(3)];
                app.Segmentation_SliceSlider.Value = round(domain_size(3)/2);
                app.Segmentation_SliceSlider.MajorTicks = round(linspace(1,domain_size(3),10));
                app.Segmentation_SliceSlider.MinorTicks = round(linspace(1,domain_size(3),19));
            end
            % Region of Interest tab
            app.ROI_ViewnormaltoDropDown.Visible = statut_GUI_3D; app.ROI_ViewnormaltoDropDown.Enable = statut_GUI_3D;
            app.ViewnormaltoDropDown_2Label.Visible = statut_GUI_3D;
            app.ROI_SliceSlider.Visible = statut_GUI_3D; app.ROI_SliceSlider.Enable = statut_GUI_3D;
            app.SliceselectionLabel.Visible = statut_GUI_3D;
            app.ROI_checkbox_updateslicer.Visible = statut_GUI_3D; app.ROI_checkbox_updateslicer.Enable = statut_GUI_3D;
            app.ROI_View3D.Visible = statut_GUI_3D; app.ROI_View3D.Enable = statut_GUI_3D;
            app.ROI_SliceSpinner.Visible = statut_GUI_3D; app.ROI_SliceSpinner.Enable = statut_GUI_3D;
            app.VieworthogonalslicesButton.Visible = statut_GUI_3D; app.VieworthogonalslicesButton.Enable = statut_GUI_3D;
            % Scaling tab
            app.Scale_SliceselectionSlider.Visible = statut_GUI_3D; app.SliceselectionSliderLabel.Visible = statut_GUI_3D;
            % Contrast tab
            app.Contrast_SliceselectionSlider.Visible = statut_GUI_3D; app.SliceselectionSliderLabel_2.Visible = statut_GUI_3D;
            % Filterint tab
            app.Filtering_SliceselectionSlider.Visible = statut_GUI_3D; app.SliceselectionSliderLabel_3.Visible = statut_GUI_3D;
            % Segmentation tab
            app.Segmentation_SliceSlider.Visible = statut_GUI_3D; app.SliceselectionLabel_2.Visible = statut_GUI_3D;
        end
        
        function [] = update_visualization_slice_ROI(app,slicer_value)
            if app.Dimension==3
                direction_ROI = app.ROI_ViewnormaltoDropDown.UserData;
                if direction_ROI==1
                    slice_ = squeeze(app.Microstructure(slicer_value,:,:));
                elseif direction_ROI==2
                    slice_ = squeeze(app.Microstructure(:,slicer_value,:));
                elseif direction_ROI==3
                    slice_ = app.Microstructure(:,:,slicer_value);
                end
            else
                slice_ = app.Microstructure(:,:,1);
            end
            image(slice_(:,:,1),'CDataMapping','scaled','Parent', app.ROI_2Dslice);
            set(app.ROI_2Dslice,'YDir','normal');
        end
        
        function [] = update_all_slicer(app)
            domain_size = size(app.Microstructure);
            app.Dimension = length(domain_size);
            % ROI
            if app.Dimension==2
                slicer_value=1;
            elseif app.Dimension==3
                %slicer_value = round(app.ROI_SliceSlider.Value);
                slicer_value=1;
            end
            % Format ROI axis
            col_ = eval(app.ROI_colormap.Value);
            colorbar(app.ROI_2Dslice)
            axis(app.ROI_2Dslice,'equal');
            axis(app.ROI_2Dslice,'tight');
            set(app.ROI_2Dslice,'YDir','normal','Colormap',col_);
            title(app.ROI_2Dslice, '')
            % ROI plot
            app.update_visualization_slice_ROI(slicer_value);
        end
        
        function [] = update_ROI_crop_table(app)
            domain_size = size(app.Microstructure);
            tmp = table([1:1:app.Dimension]',ones(app.Dimension,1),ones(app.Dimension,1).*(domain_size'),[domain_size*app.ROI_voxelsize.Value]');
            app.ROI_table_crop.Data=tmp;
            % Scale tab
            tmp = table([1:1:app.Dimension]',domain_size',[round(domain_size/app.Scale_scalingfactorEditField.Value)]');
            app.Scale_table.Data=tmp;
            % Contrast tab
            tmp = table([1:1:app.Dimension]',ones(app.Dimension,1),domain_size');
            app.Contrast_TableWhere.Data=tmp;
        end
        
        function [] = update_Convert_table(app,column)
            if column>0
                min_ = num2str( min(min(min(app.Microstructure))), '%1.1d');
                max_ = num2str( max(max(max(app.Microstructure))), '%1.1d');
                a=app.Microstructure;
                tmp=whos('a'); clear a;
                data_type = tmp.class;
                data_MB   = num2str(tmp.bytes*9.53674e-7,'%1.1f');
                % Data utilization
                if app.Convert_Utilization_checkbox.Value
                    n_used = length(unique(app.Microstructure));
                    if strcmp(data_type,'uint8')
                        n_possible = 2^8;
                    elseif strcmp(data_type,'uint16')
                        n_possible = 2^16;
                    elseif strcmp(data_type,'uint32')
                        n_possible = 2^32;
                    elseif strcmp(data_type,'uint64')
                        n_possible = 2^64;
                    elseif strcmp(data_type,'double')
                        n_possible = 2^64;
                    end
                    data_type_utilization =  num2str(n_used/n_possible,'%1.3e');
                else
                    data_type_utilization ='-';
                end
                new_column = {min_; max_; data_type; data_type_utilization;data_MB};
            else
                new_column = {'-'; '-'; '-'; '-';'-'};
                column=-column;
            end
            current_table = app.Convert_table.Data;
            current_table(:,column)=new_column;
            app.Convert_table.Data = current_table;
        end
        
        function [] = create_histogram(app,options)
            Fig_= figure;
            Fig_.Name= ['Histogram, step ' num2str(app.Operation_step)];
            Fig_.Color='white'; % Background colour
            ax_ = axes('Parent',Fig_);
            hold(ax_,'on');
            if strcmp(options.bin,'one_per_uniquevalue')
                nbins = length(unique(app.Microstructure));
                h=histogram(app.Microstructure,nbins);
                str2 = 'one bin per unique value';
            elseif strcmp(options.bin,'auto')
                h=histogram(app.Microstructure);
                str2 = 'auto bin';
            end
            title_str = {Fig_.Name,str2};
            set(h,'LineStyle','none','Normalization','probability');
            xlabel('Value');
            ylabel('Probability');
            grid(ax_,'on');
            set(ax_,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_axe);
            title (title_str,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_title);
            hold(ax_,'off');
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string([Fig_.Name '_' str2]);
                function_savefig(Fig_, app.Savefolder, filename);
            end 
        end
        
        function [] = initialize_image_contrastcorrection(app)
            col_ = eval(app.ROI_colormap.Value);
            
            % Axe before
            set(app.Contrast_UIAxes_before,'Colormap',col_);
            axis(app.Contrast_UIAxes_before,'equal');
            axis(app.Contrast_UIAxes_before,'tight');
            colorbar('peer',app.Contrast_UIAxes_before);
            min_=min(min(min(app.Microstructure_undo)));
            max_=max(max(max(app.Microstructure_undo)));
            caxis(app.Contrast_UIAxes_before,[min_ max_]);
  
            % Axe after
            set(app.Contrast_UIAxes_after,'Colormap',col_);
            axis(app.Contrast_UIAxes_after,'equal');
            axis(app.Contrast_UIAxes_after,'tight');
            colorbar('peer',app.Contrast_UIAxes_after); 
            min_=min(min(min(app.Microstructure)));
            max_=max(max(max(app.Microstructure)));
            if min_<max_
                caxis(app.Contrast_UIAxes_after,[min_ max_]);
            end
            
            % Update image
            app.update_image_contrastcorrection
            
        end
        
        function [] = update_image_contrastcorrection(app)
            if app.Dimension==3
                slicer_value = round(app.Contrast_SliceselectionSlider.Value);
                slice_before = app.Microstructure_undo(:,:,slicer_value);
                slice_after = app.Microstructure(:,:,slicer_value);
            else
                slice_before = app.Microstructure_undo(:,:,1);
                slice_after = app.Microstructure(:,:,1);
            end            
            image(slice_before(:,:,1),'CDataMapping','scaled','Parent', app.Contrast_UIAxes_before);
            image(slice_after(:,:,1),'CDataMapping','scaled','Parent', app.Contrast_UIAxes_after);
            set(app.Contrast_UIAxes_before,'YDir','normal');
            set(app.Contrast_UIAxes_after,'YDir','normal');            
        end
        
        function [] = initialize_image_filter(app)
            col_ = eval(app.ROI_colormap.Value);
            
            % Axe before
            set(app.Filtering_UIAxes_before,'Colormap',col_);
            axis(app.Filtering_UIAxes_before,'equal');
            axis(app.Filtering_UIAxes_before,'tight');
            colorbar('peer',app.Filtering_UIAxes_before);
            min_=min(min(min(app.Microstructure_undo)));
            max_=max(max(max(app.Microstructure_undo)));
            caxis(app.Filtering_UIAxes_before,[min_ max_]);
  
            % Axe after
            set(app.Filetring_UIAxes_after,'Colormap',col_);
            axis(app.Filetring_UIAxes_after,'equal');
            axis(app.Filetring_UIAxes_after,'tight');
            colorbar('peer',app.Filetring_UIAxes_after); 
            min_=min(min(min(app.Microstructure)));
            max_=max(max(max(app.Microstructure)));
            caxis(app.Filetring_UIAxes_after,[min_ max_]);       
            
            % Update image
            app.update_image_filter         
            
        end
        
        function [] = update_image_filter(app)
            if app.Dimension==3
                slicer_value = round(app.Filtering_SliceselectionSlider.Value);
                slice_before = app.Microstructure_undo(:,:,slicer_value);
                slice_after = app.Microstructure(:,:,slicer_value);
            else
                slice_before = app.Microstructure_undo(:,:,1);
                slice_after = app.Microstructure(:,:,1);
            end
            image(slice_before(:,:,1),'CDataMapping','scaled','Parent', app.Filtering_UIAxes_before);
            image(slice_after(:,:,1),'CDataMapping','scaled','Parent', app.Filetring_UIAxes_after);
            set(app.Filtering_UIAxes_before,'YDir','normal');
            set(app.Filetring_UIAxes_after,'YDir','normal');
            
        end
        
        function [] = figure_segmentation(app,localthreshold)
            if app.Dimension==3
                number_phase = app.Segmentation_NumberofphaseEditField.Value;
                tmp = app.Segmentation_Phase_UITable.Data.Variables;
                expected_vf = tmp(:,3); obtained_vf = tmp(:,4); phase_id = tmp(:,2);
                c = [get(0, 'DefaultAxesColorOrder')];
                domain_size=size(app.Microstructure);
                for direction=1:1:3
                    if direction==1
                        str_view = '1st Axis'; str_x='3rd Axis'; str_y ='2nd Axis';
                        area=domain_size(2)*domain_size(3);
                        x_=([1:1:domain_size(1)]-0.5) * app.Segmentation_VoxelsizemicrometersEditField.Value;
                        vol_=zeros(domain_size(direction),number_phase);
                        for pos=1:1:domain_size(direction)
                            for current_phase=1:1:number_phase
                                slice_segmented = squeeze(app.Microstructure(pos,:,:));
                                vol_(pos,current_phase) = sum(sum( slice_segmented == phase_id(current_phase) ))/area;
                            end
                        end

                    elseif direction==2
                        str_view = '2nd Axis'; str_x='3rd Axis'; str_y ='1st Axis';
                        area=domain_size(1)*domain_size(3);
                        x_=([1:1:domain_size(direction)]-0.5) * app.Segmentation_VoxelsizemicrometersEditField.Value;
                        vol_=zeros(domain_size(direction),number_phase);
                        for pos=1:1:domain_size(direction)
                            for current_phase=1:1:number_phase
                                slice_segmented = squeeze(app.Microstructure(:,pos,:));
                                vol_(pos,current_phase) = sum(sum( slice_segmented == phase_id(current_phase) ))/area;
                            end
                        end
                        
                    elseif direction==3
                        str_view = '3rd Axis'; str_x='2nd Axis'; str_y ='1st Axis';
                        area=domain_size(1)*domain_size(2);
                        x_=([1:1:domain_size(direction)]-0.5) * app.Segmentation_VoxelsizemicrometersEditField.Value;
                        vol_=zeros(domain_size(direction),number_phase);
                        for pos=1:1:domain_size(direction)
                            for current_phase=1:1:number_phase
                                slice_segmented = app.Microstructure(:,:,pos);
                                vol_(pos,current_phase) = sum(sum( slice_segmented == phase_id(current_phase) ))/area;
                            end
                        end                        
                    end

                    Fig_= figure;
                    Fig_.Name= ['Segmented volume, view normal to ' str_view ', step ' num2str(app.Operation_step)];
                    Fig_.Color='white'; % Background colour
                    scrsz = get(0,'ScreenSize'); % Screen resolution
                    set(Fig_, 'Position', scrsz);
                    for id_axe=1:1:4 % Iterate over axe
                        sub_axes=subplot(2,2,id_axe,'Parent',Fig_);
                        hold(sub_axes,'on');
                        if id_axe==1
                            str_legend = cell(1,number_phase*3);
                            for current_phase=1:1:number_phase
                                plot(x_,vol_(:,current_phase),'LineWidth',2,'Color',c(current_phase,:));
                                plot([x_(1) x_(end)],[expected_vf(current_phase) expected_vf(current_phase)],'LineWidth',1.5,'LineStyle','--','Color',c(current_phase,:));
                                plot([x_(1) x_(end)],[obtained_vf(current_phase) obtained_vf(current_phase)],'LineWidth',1.5,'LineStyle',':','Color',c(current_phase,:));
                                str_legend(1,3*(current_phase-1)+1)={sprintf('Phase %i, assigned to %i',current_phase,phase_id(current_phase))};
                                str_legend(1,3*(current_phase-1)+2)={sprintf('Expected avg. volume fraction %1.3f',expected_vf(current_phase))};
                                str_legend(1,3*(current_phase-1)+3)={sprintf('Obtained avg. volume fraction %1.3f',obtained_vf(current_phase))};
                            end
                            legend(sub_axes,str_legend,'Location','best');
                            title (['Volume fractions along ' str_view],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            
                        elseif id_axe==2
                            if direction==1
                                slice_segmented = squeeze(app.Microstructure(1,:,:));
                                slice_unsegmented = squeeze(app.Microstructure_undo(1,:,:));
                            elseif direction==2
                                slice_segmented = squeeze(app.Microstructure(:,1,:));
                                slice_unsegmented = squeeze(app.Microstructure_undo(:,1,:));
                            elseif direction==3
                                slice_segmented = app.Microstructure(:,:,1);
                                slice_unsegmented = app.Microstructure_undo(:,:,1);
                            end
                            title ('First slice','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                                                        
                        elseif id_axe==3
                            if direction==1
                                slice_segmented = squeeze(app.Microstructure(round(domain_size(1)/2),:,:));
                                slice_unsegmented = squeeze(app.Microstructure_undo(round(domain_size(1)/2),:,:));
                            elseif direction==2
                                slice_segmented = squeeze(app.Microstructure(:,round(domain_size(2)/2),:));
                                slice_unsegmented = squeeze(app.Microstructure_undo(:,round(domain_size(2)/2),:));
                            elseif direction==3
                                slice_segmented = app.Microstructure(:,:,round(domain_size(3)/2));
                                slice_unsegmented = app.Microstructure_undo(:,:,round(domain_size(3)/2));
                            end 
                            title ('Slice in the middle','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            
                        elseif id_axe==4
                            if direction==1
                                slice_segmented = squeeze(app.Microstructure(domain_size(1),:,:));
                                slice_unsegmented = squeeze(app.Microstructure_undo(domain_size(1),:,:));
                            elseif direction==2
                                slice_segmented = squeeze(app.Microstructure(:,domain_size(2),:));
                                slice_unsegmented = squeeze(app.Microstructure_undo(:,domain_size(2),:));
                            elseif direction==3
                                slice_segmented = app.Microstructure(:,:,domain_size(3));
                                slice_unsegmented = app.Microstructure_undo(:,:,domain_size(3));
                            end
                            title ('Last slice','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            
                        end
                        if id_axe==1
                            grid(sub_axes,'on');
                            set(sub_axes,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
                            xlabel(sub_axes, 'Position along axis (\mum)');
                            ylabel(sub_axes, 'Volume fraction');
                        else
                            imshowpair(slice_unsegmented,slice_segmented,'montage','Parent', sub_axes);
                            xlabel(sub_axes, str_x); ylabel(sub_axes, str_y);
                            set(sub_axes,'YDir','normal','Colormap',gray);
                            set(sub_axes,'Colormap',gray);
                            axis(sub_axes,'equal');
                        end
                        axis(sub_axes,'tight');
                        set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        hold(sub_axes,'off');
                    end
                    sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                    if app.Save_checkbox.Value
                        filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                        function_savefig(Fig_, app.Savefolder, filename);
                    end
                end
                
                if exist('Method_local')
                    col_ = eval(app.ROI_colormap.Value);
                    Fig_= figure;
                    Fig_.Color='white'; % Background colour
                    ax_ = axes('Parent',Fig_);
                    hold(ax_,'on');
                    Fig_.Name= 'Local threshold slice per slice, along 3rd Axis';
                    Method_local = app.Segmentation_Local_MethodDropDown.Value;
                    if strcmp(Method_local,'Local Otsu')
                        if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                            str_legend = cell(1,number_phase*2);
                        else
                            str_legend = cell(1,number_phase*1);
                        end
                    elseif strcmp(Method_local,'Local Otsu translated to match expected volume fractions')
                        if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                            str_legend = cell(1,number_phase*3);
                        else
                            str_legend = cell(1,number_phase*2);
                        end
                    end
                    
                    x_=([1:1:domain_size(direction)]-0.5) * app.Segmentation_VoxelsizemicrometersEditField.Value;
                    iter=0;
                    for current_phase=1:1:number_phase
                        iter=iter+1;
                        plot(x_,app.sav_thresholds(:,current_phase+1),'LineWidth',2,'Color',c(current_phase,:));
                        str_legend(1,iter)={sprintf('Phase %i, assigned to %i, threshold per slice',current_phase,phase_id(current_phase))};
                        if strcmp(Method_local,'Local Otsu')
                            if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                                iter=iter+1;
                                plot(x_,app.sav_thresholds_smoothed(:,current_phase+1),'LineWidth',2,'LineStyle','--','Color',c(current_phase,:));
                                str_legend(1,iter)={sprintf('Phase %i, assigned to %i, threshold per slice, smoothed (used for segmentation)',current_phase,phase_id(current_phase))};
                            end
                        elseif strcmp(Method_local,'Local Otsu translated to match expected volume fractions')
                            if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                                iter=iter+1;
                                plot(x_,app.sav_thresholds_smoothed(:,current_phase+1),'LineWidth',2,'LineStyle','--','Color',c(current_phase,:));
                                str_legend(1,iter)={sprintf('Phase %i, assigned to %i, threshold per slice, smoothed',current_phase,phase_id(current_phase))};
                                iter=iter+1;
                                plot(x_,app.sav_translated_thresholds(:,current_phase+1),'LineWidth',2,'LineStyle',':','Color',c(current_phase,:));
                                str_legend(1,iter)={sprintf('Phase %i, assigned to %i, threshold per slice, translated to match expected volume fraction',current_phase,phase_id(current_phase))};
                            else
                                iter=iter+1;
                                plot(x_,app.sav_translated_thresholds(:,current_phase+1),'LineWidth',2,'LineStyle',':','Color',c(current_phase,:));
                                str_legend(1,iter)={sprintf('Phase %i, assigned to %i, threshold per slice, translated to match expected volume fraction',current_phase,phase_id(current_phase))};
                            end
                        end
                    end
                    
                    legend(ax_,str_legend,'Location','best');
                    grid(ax_,'on');
                    set(ax_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
                    xlabel(ax_, 'Position along axis (\mum)');
                    ylabel(ax_, 'Local threshold');
                    axis(ax_,'tight');
                    set(ax_,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_axe);
                    title (Fig_.Name,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_title);
                    set(ax_,'YDir','normal');
                    hold(ax_,'off');
                    if app.Save_checkbox.Value
                        filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                        function_savefig(Fig_, app.Savefolder, filename);
                    end
                end
                
                
            end
        end
        
        function [] = initialize_image_scaling(app)
            col_ = eval(app.ROI_colormap.Value);
            
            % Axe before
            set(app.Scaling_UIAxes_before,'Colormap',col_);
            axis(app.Scaling_UIAxes_before,'equal');
            axis(app.Scaling_UIAxes_before,'tight');
            colorbar('peer',app.Scaling_UIAxes_before);
            min_=min(min(min(app.Microstructure_undo)));
            max_=max(max(max(app.Microstructure_undo)));
            caxis(app.Scaling_UIAxes_before,[min_ max_]);
  
            % Axe after
            set(app.Scaling_UIAxes_after,'Colormap',col_);
            axis(app.Scaling_UIAxes_after,'equal');
            axis(app.Scaling_UIAxes_after,'tight');
            colorbar('peer',app.Scaling_UIAxes_after); 
            min_=min(min(min(app.Microstructure)));
            max_=max(max(max(app.Microstructure)));
            caxis(app.Scaling_UIAxes_after,[min_ max_]);       
            
            % Update image
            app.update_image_scaling
            
        end

        function [] = update_image_scaling(app)
            if app.Dimension==3
                slicer_value = round(app.Scale_SliceselectionSlider.Value);
                slice_before = app.Microstructure_undo(:,:,slicer_value);
                
                slicer_value_after = round(slicer_value/app.Scale_scalingfactorEditField.Value);
                slicer_value_after = max([1 slicer_value_after]);
                new_domain_size=size(app.Microstructure);
                slicer_value_after = min([new_domain_size(3) slicer_value_after]);
                slice_after = app.Microstructure(:,:,slicer_value_after);   
            else
                slice_before = app.Microstructure_undo(:,:,1);
                slice_after = app.Microstructure(:,:,1);
            end            
            image(slice_before(:,:,1),'CDataMapping','scaled','Parent', app.Scaling_UIAxes_before);
            image(slice_after(:,:,1),'CDataMapping','scaled','Parent', app.Scaling_UIAxes_after);
            set(app.Scaling_UIAxes_before,'YDir','normal');
            set(app.Scaling_UIAxes_after,'YDir','normal');            
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
            app.Turn_on_off('off',app.RegionofInterestViewTab)
            app.Turn_on_off('off',app.FormatconversionTab)
            app.Turn_on_off('off',app.UpdownscalingTab)
            app.Turn_on_off('off',app.ImagequalityTab)
            app.Turn_on_off('off',app.PhasereassignmentTab)
            app.Turn_on_off('off',app.ContrastcorrectionTab)
            app.Turn_on_off('off',app.SegmentationsensitivityTab)
            app.Turn_on_off('off',app.FilteringTab)
            app.Turn_on_off('off',app.SegmentationTab)
            app.Save_choiceFont_DropDown.Items = listfonts;
            app.Save_choiceFont_DropDown.Value = 'Times New Roman';
            app.Scale_table.RowName = {};
            numberiteration = 5;
            tmp=zeros(numberiteration,2);
            tmp(:,1)=1:1:numberiteration;
            tmp(:,2)=ones(numberiteration,1)*25;
            app.Filtering_AnisotropicGradient_Table.Data = table(tmp(:,1),tmp(:,2));
        end

        % Menu selected function: LoadvolumeMenu
        function LoadvolumeMenuSelected(app, event)
            str_dialogbox = 'Select a .tif file';
            
            % Open dialog box to choose file path
            [FileName,PathName,~] = uigetfile({'*.tif;*.tiff','Tif image (*.tif, *.tiff)'},str_dialogbox);
            if FileName==0
                % User clicked cancel button or closed the dialog box
                app.VolumeLoadedStatus.Text ='No volume loaded. Please load a volume: Menu/Load volume';
                app.VolumeLoadedStatus.FontColor = 'r';
                app.update_GUI_TabInstruction_status('off')
                app.Microstructure=[]; % Reset Microstructure
                app.Microstructure_undo=[]; % Reset Microstructure
                app.History_step=[]; app.History_operations=[]; app.History_parameters=[]; app.History_elapsedtime=[]; % Reset log
                app.Operation_step=0;
                app.DefaultsavefoldernoneselectedLabel.Text = 'Current save folder: none selected';
                app.DefaultfilenamenoneLabel.Text = 'Current file name: none';
                % Turn off GUI
                app.Turn_on_off('off',app.RegionofInterestViewTab)
                app.Turn_on_off('off',app.FormatconversionTab)
                app.Turn_on_off('off',app.UpdownscalingTab)
                app.Turn_on_off('off',app.ImagequalityTab)
                app.Turn_on_off('off',app.PhasereassignmentTab)
                app.Turn_on_off('off',app.ContrastcorrectionTab)
                app.Turn_on_off('off',app.SegmentationsensitivityTab)
                app.Turn_on_off('off',app.FilteringTab)
                app.Turn_on_off('off',app.SegmentationTab)
                
            else
                tic
                app.Loading_message.Visible = 'on'; pause(0.1);
                [Microstructure_local,outcome] = function_load_tif([PathName FileName]);
                app.Loading_message.Visible = 'off';
                elapsed_time = toc;
                if outcome % Success to import
                    app.VolumeLoadedStatus.Text ='Volume successfully loaded';
                    app.VolumeLoadedStatus.FontColor = 'k';
                    domain_size = size(Microstructure_local);
                    app.Dimension = length(domain_size);
                    tmp=whos('Microstructure_local');
                    
                    app.Fileloaded_Path.Text = PathName;
                    app.Fileloaded_Filename.Text = FileName;
                    app.Fileloaded_dimension.Text = ['Image dimension: ' num2str(app.Dimension)];
                    if app.Dimension==2
                        app.Fileloaded_numberofvoxel.Text = sprintf('Number of voxels: %i x %i', domain_size(1), domain_size(2));
                    elseif app.Dimension==3
                        app.Fileloaded_numberofvoxel.Text = sprintf('Number of voxels: %i x %i x %i', domain_size(1), domain_size(2), domain_size(3));
                    end
                    app.Fileloaded_memory_datatype.Text = sprintf('Datatype and memory: %s and %1.1f MB', tmp.class, tmp.bytes*9.53674e-7);
                    app.Fileloaded_min_max.Text = sprintf('Minimum and maximum values: %1.3f and %1.3f', min(min(min(Microstructure_local))), max(max(max(Microstructure_local))));
                    app.update_GUI_TabInstruction_status('on')
                    
                    % Update global variables
                    app.Microstructure=Microstructure_local; % Update Microstructure
                    app.Microstructure_undo = app.Microstructure;
                    % Update save folder
                    app.Savefolder = PathName;
                    app.DefaultsavefoldernoneselectedLabel.Text = app.Savefolder;
                    % Update history log
                    app.Operation_step=1;
                    [~,filesave,~] = fileparts(FileName); % Remove .tif
                    app.defaultfilesave= [filesave '_step_' num2str(app.Operation_step) '.tif'];
                    app.DefaultfilenamenoneLabel.Text = ['Current file name: ' app.defaultfilesave];
                    app.History_step=["-";"-";"-";"-";"-";"1"];
                    app.History_operations=["Start date";"User name";"Computer name";"Operating system";"MATLAB version";"Loading file"];
                    time_ = convertCharsToStrings(char(datetime('now','TimeZone','local','Format','eeee, HH:mm:ss Z, MMMM d, y')));
                    username_ = convertCharsToStrings(getenv('USERNAME'));
                    computername_ = convertCharsToStrings(getenv('COMPUTERNAME'));
                    osname_ = convertCharsToStrings(getenv('OS'));
                    Matlabversion_ = convertCharsToStrings(version);
                    filepath = [PathName FileName];
                    app.Filename = FileName;
                    
                    app.History_parameters=[time_;username_;computername_;osname_;Matlabversion_;filepath];
                    app.History_elapsedtime=["n/a";"n/a";"n/a";"n/a";"n/a";[num2str(elapsed_time,'%1.1f') 's']];
                    app.update_history_log
                    % Update all slicer
                    app.update_GUI_dimension
                    app.update_all_slicer
                    % update tables
                    app.update_ROI_crop_table
                    row_names = {'Minimum';'Maximum';'Format';'Utilization ratio';'Size (MB)'};
                    val_current = {'-';'-';'-';'-';'-'};
                    val_new = {'-';'-';'-';'-';'-'};
                    app.Convert_table.Data = table(row_names,val_current,val_new);
                    app.update_Convert_table(2)
                    % Initialize segmentation tab
                    app.Segmentation_NumberofphaseEditFieldValueChanged                    
                    % Turn on GUI
                    app.Turn_on_off('on',app.RegionofInterestViewTab)
                    app.Turn_on_off('on',app.FormatconversionTab)
                    app.Turn_on_off('on',app.UpdownscalingTab)
                    app.Turn_on_off('on',app.ImagequalityTab)
                    app.Turn_on_off('on',app.PhasereassignmentTab)
                    app.Turn_on_off('on',app.ContrastcorrectionTab)
                    app.Turn_on_off('on',app.SegmentationsensitivityTab)
                    app.Turn_on_off('on',app.FilteringTab)
                    app.Turn_on_off('on',app.SegmentationTab)
                    
                else
                    app.VolumeLoadedStatus.Text ='Error: volume cannot be loaded! File is not a .tif file';
                    app.VolumeLoadedStatus.FontColor = 'r';
                    app.update_GUI_TabInstruction_status('off')
                    app.Microstructure=[]; % Reset Microstructure
                    app.Microstructure_undo=[]; % Reset Microstructure
                    app.Dimension=[];
                    app.History_step=[]; app.History_operations=[]; app.History_parameters=[]; app.History_elapsedtime=[]; % Reset log
                    app.Operation_step=0;
                    app.DefaultsavefoldernoneselectedLabel.Text = 'Current save folder: none selected';
                    app.DefaultfilenamenoneLabel.text = 'Current file name: none';
                    % Turn off GUI
                    app.Turn_on_off('off',app.RegionofInterestViewTab)
                    app.Turn_on_off('off',app.FormatconversionTab)
                    app.Turn_on_off('off',app.UpdownscalingTab)
                    app.Turn_on_off('off',app.ImagequalityTab)
                    app.Turn_on_off('off',app.PhasereassignmentTab)
                    app.Turn_on_off('off',app.ContrastcorrectionTab)
                    app.Turn_on_off('off',app.SegmentationsensitivityTab)
                    app.Turn_on_off('off',app.FilteringTab)
                    app.Turn_on_off('off',app.SegmentationTab)
                end
            end
            
        end

        % Menu selected function: DefaultlocationMenu
        function SavevolumeMenuSelected(app, event)
            if ~isempty(app.Microstructure) % Only if microstructure exists
                % Save volume
                function_save_tif(app.Microstructure, [app.Savefolder app.defaultfilesave]);
                % Save history log
                [~,filesave,~] = fileparts(app.defaultfilesave); % Remove .tif
                DATA_writetable.sheet(1).name='History log';
                DATA_writetable.sheet(1).table=app.UITable_History.Data;
                Function_Writetable(app.Savefolder,filesave,DATA_writetable)
            end
        end

        % Value changed function: ROI_SliceSlider
        function ROI_SliceSliderValueChanged(app, event)
            slicer_value = round(app.ROI_SliceSlider.Value);
            app.ROI_SliceSpinner.Value = slicer_value;
            app.update_visualization_slice_ROI(slicer_value);
            
        end

        % Value changing function: ROI_SliceSlider
        function ROI_SliceSliderValueChanging(app, event)
            slicer_value = round(event.Value);
            if app.ROI_checkbox_updateslicer.Value
                app.ROI_SliceSpinner.Value = slicer_value;
                app.update_visualization_slice_ROI(slicer_value);
            end
            
        end

        % Value changed function: ROI_ViewnormaltoDropDown
        function ROI_ViewnormaltoDropDownValueChanged(app, event)
            domain_size = size(app.Microstructure);
            tmp = app.ROI_ViewnormaltoDropDown.Value;
            if strcmp(tmp,'Axe 1')
                direction_ROI=1;
                str_x = '3rd Axis'; str_y = '2nd Axis';
            elseif strcmp(tmp,'Axe 2')
                direction_ROI=2;
                str_x = '3rd Axis'; str_y = '1st Axis';
            elseif strcmp(tmp,'Axe 3')
                direction_ROI=3;
                str_x = '2nd Axis'; str_y = '1st Axis';
            end
            
            app.ROI_ViewnormaltoDropDown.UserData = direction_ROI;
            
            app.ROI_SliceSlider.Limits = [1 domain_size(direction_ROI)];
            app.ROI_SliceSlider.Value = round(domain_size(direction_ROI)/2);
            app.ROI_SliceSlider.MajorTicks = round(linspace(1,domain_size(direction_ROI),5));
            app.ROI_SliceSlider.MinorTicks = round(linspace(1,domain_size(direction_ROI),9));
            
            app.ROI_SliceSpinner.Limits = [1 domain_size(direction_ROI)];
            app.ROI_SliceSpinner.Value = app.ROI_SliceSlider.Value;
            
            slicer_value = round(app.ROI_SliceSlider.Value);
            app.update_visualization_slice_ROI(slicer_value);
            
            if app.Dimension==2
                xlabel(app.ROI_2Dslice, '2nd Axis')
                ylabel(app.ROI_2Dslice, '1st Axis')
            else
                xlabel(app.ROI_2Dslice, str_x)
                ylabel(app.ROI_2Dslice, str_y)
            end
            col_ = eval(app.ROI_colormap.Value);
            colorbar(app.ROI_2Dslice)
            set(app.ROI_2Dslice,'YDir','normal','Colormap',col_);
            axis(app.ROI_2Dslice,'equal');
            axis(app.ROI_2Dslice,'tight');
            title(app.ROI_2Dslice, '')
        end

        % Button pushed function: ROI_View3D
        function ROI_View3DButtonPushed(app, event)
            Fig_3D=figure;
            col_ = eval(app.ROI_colormap.Value);
            volshow(app.Microstructure,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_);
        end

        % Value changed function: ROI_colormap
        function ROI_colormapValueChanged(app, event)
            col_ = eval(app.ROI_colormap.Value);
            colorbar(app.ROI_2Dslice)
            set(app.ROI_2Dslice,'YDir','normal','Colormap',col_);
        end

        % Value changed function: ROI_SliceSpinner
        function ROI_SliceSpinnerValueChanged(app, event)
            slicer_value = app.ROI_SliceSpinner.Value;
            app.ROI_SliceSlider.Value = slicer_value;
            app.update_visualization_slice_ROI(slicer_value);
            
        end

        % Value changing function: ROI_SliceSpinner
        function ROI_SliceSpinnerValueChanging(app, event)
            slicer_value = event.Value;
            if app.ROI_checkbox_updateslicer.Value
                app.ROI_SliceSlider.Value = slicer_value;
                app.update_visualization_slice_ROI(slicer_value);
            end
        end

        % Value changed function: ROI_voxelsize
        function ROI_voxelsizeValueChanged(app, event)
            app.update_ROI_crop_table;
            
        end

        % Callback function: ROI_table_crop, ROI_table_crop
        function ROI_table_cropCellEdit(app, event)
            tmp = app.ROI_table_crop.Data.Variables; % Extract array
            tmp(isnan(tmp))=-1e9; % NaN -> -1e9
            tmp(:,2) = max( ones(app.Dimension,1),tmp(:,2) ); % > 1
            tmp(:,3) = max( ones(app.Dimension,1),tmp(:,3) );
            tmp(:,2) = min( size(app.Microstructure)',tmp(:,2) ); % < Domain size
            tmp(:,3) = min( size(app.Microstructure)',tmp(:,3) );
            tmp(:,2:3) = sort(tmp(:,2:3),2); % Start < End
            app.ROI_table_crop.Data=table(tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4)); % Insert table
            
            dimension_equality = sum(tmp(:,2) == tmp(:,3));
            if strcmp(app.ROI_crop_UndoButton.Enable,'off')
                if app.Dimension==3 && dimension_equality>1 || app.Dimension==2 && dimension_equality>0
                    app.ROI_crop_DoButton.Enable='off';
                else
                    app.ROI_crop_DoButton.Enable='on';
                end
            end
            
            % Plot ROI line
            if app.Dimension==2 || strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 3')
                slicer_value = round(app.ROI_SliceSlider.Value);
                app.update_visualization_slice_ROI(slicer_value);
                hold(app.ROI_2Dslice,'on');
                plot(app.ROI_2Dslice,[tmp(2,2) tmp(2,3)],[tmp(1,2) tmp(1,2)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(2,2) tmp(2,3)],[tmp(1,3) tmp(1,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(2,2) tmp(2,2)],[tmp(1,2) tmp(1,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(2,3) tmp(2,3)],[tmp(1,2) tmp(1,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                hold(app.ROI_2Dslice,'off');
            elseif strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 2')
                slicer_value = round(app.ROI_SliceSlider.Value);
                app.update_visualization_slice_ROI(slicer_value);
                hold(app.ROI_2Dslice,'on');
                plot(app.ROI_2Dslice,[tmp(3,2) tmp(3,3)],[tmp(1,2) tmp(1,2)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(3,2) tmp(3,3)],[tmp(1,3) tmp(1,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(3,2) tmp(3,2)],[tmp(1,2) tmp(1,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(3,3) tmp(3,3)],[tmp(1,2) tmp(1,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                hold(app.ROI_2Dslice,'off');
            elseif strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 1')
                slicer_value = round(app.ROI_SliceSlider.Value);
                app.update_visualization_slice_ROI(slicer_value);
                hold(app.ROI_2Dslice,'on');
                plot(app.ROI_2Dslice,[tmp(3,2) tmp(3,3)],[tmp(2,2) tmp(2,2)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(3,2) tmp(3,3)],[tmp(2,3) tmp(2,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(3,2) tmp(3,2)],[tmp(2,2) tmp(2,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                plot(app.ROI_2Dslice,[tmp(3,3) tmp(3,3)],[tmp(2,2) tmp(2,3)],'LineWidth',2,'LineStyle',"--",'Color','r');
                hold(app.ROI_2Dslice,'off');
            end
            
            
        end

        % Button pushed function: ROI_crop_DoButton
        function ROI_crop_DoButtonPushed(app, event)
            tmp = app.ROI_table_crop.Data.Variables;
            
            % Format parameter string
            [n,~]=size(tmp);
            app.Do_parameters = [];
            for k=1:1:n
                app.Do_parameters = [app.Do_parameters sprintf('Axe %i: %i-%i ',k,tmp(k,2),tmp(k,3))];
            end
            
            tic
            app.Microstructure_undo = app.Microstructure;
            if app.Dimension==3
                app.Microstructure = app.Microstructure( tmp(1,2):tmp(1,3), tmp(2,2):tmp(2,3), tmp(3,2):tmp(3,3) );
            else
                app.Microstructure = app.Microstructure( tmp(1,2):tmp(1,3), tmp(2,2):tmp(2,3));
            end
            app.Do_elapsedtime = toc;
            
            new_size = size(app.Microstructure);
            if sum( new_size==1)==1
                if new_size(1)==1
                    app.Microstructure = squeeze(app.Microstructure(1,:,:));
                elseif new_size(2)==1
                    app.Microstructure = squeeze(app.Microstructure(:,1,:));
                elseif new_size(3)==1
                    app.Microstructure = app.Microstructure(:,:,1);
                end
            end
            
            app.update_GUI_dimension
            app.update_all_slicer
            app.update_ROI_crop_table
            app.DoUndoSave_button('ROI','crop','do'); % Update button
        end

        % Button pushed function: ROI_crop_UndoButton
        function ROI_crop_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.update_GUI_dimension
            app.update_all_slicer
            app.update_ROI_crop_table
            app.DoUndoSave_button('ROI','crop','undo'); % Update button
        end

        % Button pushed function: ROI_crop_SaveButton
        function ROI_crop_SaveButtonPushed(app, event)
            app.DoUndoSave_button('ROI','crop','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; "Crop volume"];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
        end

        % Value changed function: ROI_angle_spinner
        function ROI_angle_spinnerValueChanged(app, event)
            angle_deg = app.ROI_angle_spinner.Value;
            domain_size = size(app.Microstructure);
            nline=32;
            if app.Dimension==3
                slicer_value = round(app.ROI_SliceSlider.Value);
                app.update_visualization_slice_ROI(slicer_value);
                hold(app.ROI_2Dslice,'on');
                if strcmp(app.ROI_angle_choice_DropDown.Value,'Normal to axe 1') && strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 1')
                    xlimsav=[1 domain_size(3)]; ylimsav=[1 domain_size(2)];
                    y0s=unique([linspace(-3*domain_size(2),1,nline+2) linspace(1,3*domain_size(2),nline+2)]);
                    for kline=1:1:length(y0s)
                        y0 = round(y0s(kline));
                        y1 = y0 - tan(deg2rad(angle_deg))*(domain_size(3)-1);
                        plot(app.ROI_2Dslice,[1 domain_size(3)],[y0 y1],'LineWidth',2,'LineStyle',"--",'Color','r');
                    end
                    hold(app.ROI_2Dslice,'off');
                    xlim(app.ROI_2Dslice,xlimsav);
                    ylim(app.ROI_2Dslice,ylimsav);
                    
                elseif strcmp(app.ROI_angle_choice_DropDown.Value,'Normal to axe 2') && strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 2')
                    xlimsav=[1 domain_size(3)]; ylimsav=[1 domain_size(1)];
                    y0s=unique([linspace(-3*domain_size(1),1,nline+2) linspace(1,3*domain_size(1),nline+2)]);
                    for kline=1:1:length(y0s)
                        y0 = round(y0s(kline));
                        y1 = y0 - tan(deg2rad(angle_deg))*(domain_size(3)-1);
                        plot(app.ROI_2Dslice,[1 domain_size(3)],[y0 y1],'LineWidth',2,'LineStyle',"--",'Color','r');
                    end
                    hold(app.ROI_2Dslice,'off');
                    xlim(app.ROI_2Dslice,xlimsav);
                    ylim(app.ROI_2Dslice,ylimsav);
                    
                elseif strcmp(app.ROI_angle_choice_DropDown.Value,'Normal to axe 3') && strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 3')
                    xlimsav=[1 domain_size(2)]; ylimsav=[1 domain_size(1)];
                    y0s=unique([linspace(-3*domain_size(1),1,nline+2) linspace(1,3*domain_size(1),nline+2)]);
                    for kline=1:1:length(y0s)
                        y0 = round(y0s(kline));
                        y1 = y0 - tan(deg2rad(angle_deg))*(domain_size(2)-1);
                        plot(app.ROI_2Dslice,[1 domain_size(2)],[y0 y1],'LineWidth',2,'LineStyle',"--",'Color','r');
                    end
                    hold(app.ROI_2Dslice,'off');
                    xlim(app.ROI_2Dslice,xlimsav);
                    ylim(app.ROI_2Dslice,ylimsav);
                end
            else
                slicer_value = 1;
                app.update_visualization_slice_ROI(slicer_value);
                hold(app.ROI_2Dslice,'on');
                xlimsav=[1 domain_size(2)]; ylimsav=[1 domain_size(1)];
                y0s=unique([linspace(-3*domain_size(1),1,nline+2) linspace(1,3*domain_size(1),nline+2)]);
                for kline=1:1:length(y0s)
                    y0 = round(y0s(kline));
                    y1 = y0 - tan(deg2rad(angle_deg))*(domain_size(2)-1);
                    plot(app.ROI_2Dslice,[1 domain_size(2)],[y0 y1],'LineWidth',2,'LineStyle',"--",'Color','r');
                end
                hold(app.ROI_2Dslice,'off');
                xlim(app.ROI_2Dslice,xlimsav);
                ylim(app.ROI_2Dslice,ylimsav);
            end

        end

        % Button pushed function: ROI_rotation_DoButton
        function ROI_rotation_DoButtonPushed(app, event)
            angle_deg = app.ROI_angle_spinner.Value;
            t=[]; % Initialize empty
            % Format parameter string
            app.Do_parameters = [app.ROI_angle_choice_DropDown.Value ': ' num2str(angle_deg,'%1.1f') ' degrees'];
            
            tic
            app.Microstructure_undo = app.Microstructure;
            % https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
            %Transformation matrix that rotates the image around the x-axis
            angle_rad = deg2rad(angle_deg);
            if app.Dimension==3
                if strcmp(app.ROI_angle_choice_DropDown.Value,'Normal to axe 1') && strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 1')
                    t = [cos(angle_rad)  0      -sin(angle_rad)   0
                        0             1              0     0
                        sin(angle_rad)    0       cos(angle_rad)   0
                        0             0              0     1];
                elseif strcmp(app.ROI_angle_choice_DropDown.Value,'Normal to axe 2') && strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 2')
                    t = [1 0 0 0
                        0 cos(-angle_rad)  sin(-angle_rad) 0
                        0 -sin(-angle_rad) cos(-angle_rad) 0
                        0             0              0       1];
                elseif strcmp(app.ROI_angle_choice_DropDown.Value,'Normal to axe 3') && strcmp(app.ROI_ViewnormaltoDropDown.Value,'Axe 3')
                    %Transformation matrix that rotates the image around the z-axis
                    t = [cos(angle_rad) sin(angle_rad) 0 0
                        -sin(angle_rad)  cos(angle_rad) 0 0
                        0 0 1 0
                        0 0 0 1];
                end
            else
                %Transformation matrix that rotates the image around the z-axis
                t = [cos(angle_rad) sin(angle_rad) 0
                    -sin(angle_rad)  cos(angle_rad) 0
                    0 0 1];
            end
            if ~isempty(t)
                % Then pass the matrix to the affine3d object constructor.
                if app.Dimension==3
                    t_form = affine3d(t);
                else
                    t_form = affine2d(t);
                end
                app.Microstructure = imwarp(app.Microstructure,t_form);
                app.Do_elapsedtime = toc;
                
                app.update_all_slicer
                app.update_ROI_crop_table
                app.DoUndoSave_button('ROI','rotation','do'); % Update button
            end
        end

        % Button pushed function: ROI_rotation_UndoButton
        function ROI_rotation_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.update_all_slicer
            app.update_ROI_crop_table
            app.DoUndoSave_button('ROI','rotation','undo'); % Update button
        end

        % Button pushed function: ROI_rotation_SaveButton
        function ROI_rotation_SaveButtonPushed(app, event)
            app.DoUndoSave_button('ROI','rotation','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; "Rotation volume"];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
        end

        % Button pushed function: ROI_flipswap_DoButton
        function ROI_flipswap_DoButtonPushed(app, event)
            
            % Format parameter string
            app.Do_parameters ='No parameters';
            
            tic
            app.Microstructure_undo = app.Microstructure;
            domain_size = size(app.Microstructure); New_domain_size=zeros(1,3);
            
            if app.Dimension==3
                if strcmp(app.ROI_flipswapDropDown.Value,'Flip axis 1')
                    tmp=zeros(domain_size);
                    for k=1:1:domain_size(1)
                        slice=app.Microstructure(k,:,:);
                        tmp(domain_size(1)-k+1,:,:)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Flip axis 2')
                    tmp=zeros(domain_size);
                    for k=1:1:domain_size(2)
                        slice=app.Microstructure(:,k,:);
                        tmp(:,domain_size(2)-k+1,:)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Flip axis 3')
                    tmp=zeros(domain_size);
                    for k=1:1:domain_size(3)
                        slice=app.Microstructure(:,:,k);
                        tmp(:,:,domain_size(3)-k+1)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Swap axis 1 with axis 2')
                    New_domain_size(1)=domain_size(2); New_domain_size(2)=domain_size(1); New_domain_size(3)=domain_size(3);
                    tmp=zeros(New_domain_size);
                    for k=1:1:domain_size(2)
                        slice=app.Microstructure(:,k,:);
                        tmp(k,:,:)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Swap axis 1 with axis 3')
                    New_domain_size(1)=domain_size(3); New_domain_size(2)=domain_size(2); New_domain_size(3)=domain_size(1);
                    tmp=zeros(New_domain_size);
                    for k=1:1:domain_size(3)
                        slice=app.Microstructure(:,:,k)';
                        tmp(k,:,:)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Swap axis 2 with axis 3')
                    New_domain_size(1)=domain_size(1); New_domain_size(2)=domain_size(3); New_domain_size(3)=domain_size(2);
                    tmp=zeros(New_domain_size);
                    for k=1:1:domain_size(3)
                        slice=app.Microstructure(:,:,k);
                        tmp(:,k,:)=slice;
                    end
                end
                
            else
                if strcmp(app.ROI_flipswapDropDown.Value,'Flip axis 1')
                    tmp=zeros(domain_size);
                    for k=1:1:domain_size(1)
                        slice=app.Microstructure(k,:);
                        tmp(domain_size(1)-k+1,:)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Flip axis 2')
                    tmp=zeros(domain_size);
                    for k=1:1:domain_size(2)
                        slice=app.Microstructure(:,k);
                        tmp(:,domain_size(2)-k+1)=slice;
                    end
                elseif strcmp(app.ROI_flipswapDropDown.Value,'Swap axis 1 with axis 2')
                    New_domain_size=zeros(1,2); New_domain_size(1)=domain_size(2); New_domain_size(2)=domain_size(1);
                    tmp=zeros(New_domain_size);
                    for k=1:1:domain_size(2)
                        slice=app.Microstructure(:,k);
                        tmp(k,:)=slice;
                    end
                end
                
            end
            app.Microstructure=tmp;
            
            app.Do_elapsedtime = toc;
            app.update_all_slicer
            app.update_ROI_crop_table
            app.DoUndoSave_button('ROI','flipswap','do'); % Update button
        end

        % Button pushed function: ROI_flipswap_UndoButton
        function ROI_flipswap_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.update_all_slicer
            app.update_ROI_crop_table
            app.DoUndoSave_button('ROI','flipswap','undo'); % Update button
        end

        % Button pushed function: ROI_flipswap_SaveButton
        function ROI_flipswap_SaveButtonPushed(app, event)
            app.DoUndoSave_button('ROI','flipswap','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; app.ROI_flipswapDropDown.Value];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
        end

        % Value changed function: ROI_modify_axeAxeDropDown
        function ROI_modify_axeAxeDropDownValueChanged(app, event)
            value = app.ROI_modify_axeAxeDropDown.Value;
            domain_size = size(app.Microstructure);
            if strcmp(value,'Axe 1')
                app.ROI_modify_FromsliceA.Limits = [1 domain_size(1)]; app.ROI_modify_TosliceB.Limits = [1 domain_size(1)];
            elseif strcmp(value,'Axe 2')
                app.ROI_modify_FromsliceA.Limits = [1 domain_size(2)]; app.ROI_modify_TosliceB.Limits = [1 domain_size(2)];
            elseif strcmp(value,'Axe 3')
                app.ROI_modify_FromsliceA.Limits = [1 domain_size(3)]; app.ROI_modify_TosliceB.Limits = [1 domain_size(3)];
            end
        end

        % Value changed function: ROI_modify_FromsliceA, 
        % ROI_modify_TosliceB
        function ROI_modify_FromsliceAValueChanged(app, event)
            fromA = app.ROI_modify_FromsliceA.Value;
            toB = app.ROI_modify_TosliceB.Value;
            if fromA>toB
                app.ROI_modify_FromsliceA.Value=toB;
                app.ROI_modify_TosliceB.Value=fromA;
            end
        end

        % Button pushed function: ROI_modify_DoButton
        function ROI_modify_DoButtonPushed(app, event)
            fromA = app.ROI_modify_FromsliceA.Value;
            toB = app.ROI_modify_TosliceB.Value;
            if fromA<toB && (toB-fromA)>1
                % Format parameter string
                app.Do_parameters =sprintf('Along %s, from slice %i (not included) to %i (not included)',app.ROI_modify_axeAxeDropDown.Value,fromA,toB);
                
                tic
                if strcmp(app.ROI_modify_axeAxeDropDown.Value,'Axe 1')
                    direction = 1;
                elseif strcmp(app.ROI_modify_axeAxeDropDown.Value,'Axe 2')
                    direction = 2;
                elseif strcmp(app.ROI_modify_axeAxeDropDown.Value,'Axe 3')
                    direction = 3;
                end
                app.Microstructure_undo = app.Microstructure;
                domain_size = size(app.Microstructure); New_domain_size=zeros(1,3);
                if strcmp(app.ROI_modify_ActionDropDown.Value,'Remove slices')
                    if app.Dimension==3
                        if direction==1
                            app.Microstructure(fromA+1:toB-1,:,:)=[];
                        elseif direction==2
                            app.Microstructure(:,fromA+1:toB-1,:)=[];
                        elseif direction==3
                            app.Microstructure(:,:,fromA+1:toB-1)=[];
                        end
                    else
                        if direction==1
                            app.Microstructure(fromA+1:toB-1,:)=[];
                        elseif direction==2
                            app.Microstructure(:,fromA+1:toB-1)=[];
                        end
                    end
                elseif strcmp(app.ROI_modify_ActionDropDown.Value,'Interpolate slices (grey level)')
                    a = (0-1)/(toB-fromA);
                    b = 0-a*toB;
                    for k=fromA+1:1:toB-1
                        ya = a*k+b;
                        if app.Dimension==3
                            if direction==1
                                app.Microstructure(k,:,:) = round(ya*double(app.Microstructure(fromA,:,:)) + (1-ya)*double(app.Microstructure(toB,:,:)));
                            elseif direction==2
                                app.Microstructure(:,k,:) = round(ya*double(app.Microstructure(:,fromA,:)) + (1-ya)*double(app.Microstructure(:,toB,:)));
                            elseif direction==3
                                app.Microstructure(:,:,k) = round(ya*double(app.Microstructure(:,:,fromA)) + (1-ya)*double(app.Microstructure(:,:,toB)));
                            end
                        else
                            if direction==1
                                app.Microstructure(k,:) = round(ya*double(app.Microstructure(fromA,:)) + (1-ya)*double(app.Microstructure(toB,:)));
                            elseif direction==2
                                app.Microstructure(:,k) = round(ya*double(app.Microstructure(:,fromA)) + (1-ya)*double(app.Microstructure(:,toB)));
                            end
                        end
                    end
                elseif strcmp(app.ROI_modify_ActionDropDown.Value,'Interpolate slices (label)')
                    A0 = fromA+1; A1 = round((fromA+toB)/2);
                    B0 = A1+1; B1 = toB-1;
                    if app.Dimension==3
                        if direction==1
                            for k=A0:A1
                                app.Microstructure(k,:,:) = app.Microstructure(fromA,:,:);
                            end
                            if B0<=B1
                                for k=B0:B1
                                    app.Microstructure(k,:,:) = app.Microstructure(toB,:,:);
                                end
                            end
                        elseif direction==2
                            for k=A0:A1
                                app.Microstructure(:,k,:) = app.Microstructure(:,fromA,:);
                            end
                            if B0<=B1
                                for k=B0:B1
                                    app.Microstructure(:,k,:) = app.Microstructure(:,toB,:);
                                end
                            end
                        elseif direction==3
                            for k=A0:A1
                                app.Microstructure(:,:,k) = app.Microstructure(:,:,fromA);
                            end
                            if B0<=B1
                                for k=B0:B1
                                    app.Microstructure(:,:,k) = app.Microstructure(:,:,toB);
                                end
                            end
                        end
                    else
                        if direction==1
                            for k=A0:A1
                                app.Microstructure(k,:) = app.Microstructure(fromA,:);
                            end
                            if B0<=B1
                                for k=B0:B1
                                    app.Microstructure(k,:) = app.Microstructure(toB,:);
                                end
                            end
                        elseif direction==2
                            for k=A0:A1
                                app.Microstructure(:,k) = app.Microstructure(:,fromA);
                            end
                            if B0<=B1
                                for k=B0:B1
                                    app.Microstructure(:,k) = app.Microstructure(:,toB);
                                end
                            end
                        end
                    end
                end
                app.Do_elapsedtime = toc;
                app.update_all_slicer
                app.update_ROI_crop_table
                app.DoUndoSave_button('ROI','modify','do'); % Update button
            end
            
        end

        % Button pushed function: ROI_modify_UndoButton
        function ROI_modify_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.update_all_slicer
            app.update_ROI_crop_table
            app.DoUndoSave_button('ROI','modify','undo'); % Update button
        end

        % Button pushed function: ROI_modify_SaveButton
        function ROI_modify_SaveButtonPushed(app, event)
            app.DoUndoSave_button('ROI','modify','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; app.ROI_modify_ActionDropDown.Value];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
        end

        % Menu selected function: CustomlocationMenu
        function CustomlocationMenuSelected(app, event)
            if ~isempty(app.Microstructure) % Only if microstructure exists
                [filesave,pathsave,indx] = uiputfile('*.tif','Select the folder where to save the current state of the volume, and write its file name',app.Filename);
                if indx~=0 % User clicked save button
                    % Save volume
                    function_save_tif(app.Microstructure, [pathsave filesave]);
                    % Save history log
                    [~,filesave,~] = fileparts(filesave); % Remove .tif
                    DATA_writetable.sheet(1).name='History log';
                    DATA_writetable.sheet(1).table=app.UITable_History.Data;
                    Function_Writetable(pathsave,filesave,DATA_writetable)
                end
            end
        end

        % Button pushed function: ClicktoselectsavefolderButton
        function ClicktoselectsavefolderButtonPushed(app, event)
            if ~isempty(app.Savefolder)
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
                app.DefaultsavefoldernoneselectedLabel.Text = app.Savefolder;
            end
        end

        % Image clicked function: Logo
        function LogoImageClicked(app, event)
            url = 'https://www.nrel.gov/transportation/';
            web(url)
        end

        % Menu selected function: histogram_standard_wholevolume_auto
        function histogram_standard_wholevolume_autoSelected(app, event)
            options.bin = 'auto';
            app.create_histogram(options);
        end

        % Button pushed function: Convert_DoButton
        function Convert_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            new_data_type = app.Convert_Choice_DropDown.Value;
            if strcmp(new_data_type,'8-bit unsigned integer arrays')
                app.Microstructure =uint8(app.Microstructure);
            elseif strcmp(new_data_type,'8-bit unsigned integer arrays (rescale)')
                app.Microstructure =im2uint8(app.Microstructure); % rescaling or offsetting the data as necessary.
            elseif strcmp(new_data_type,'16-bit unsigned integer arrays')
                app.Microstructure = uint16(app.Microstructure);
            elseif strcmp(new_data_type,'16-bit unsigned integer arrays (rescale)')
                app.Microstructure = im2uint16(app.Microstructure); % rescaling or offsetting the data as necessary.
            elseif strcmp(new_data_type,'32-bit unsigned integer arrays')
                app.Microstructure = uint32(app.Microstructure);
            elseif strcmp(new_data_type,'64-bit unsigned integer arrays')
                app.Microstructure = uint64(app.Microstructure);
            elseif strcmp(new_data_type,'double precision array')
                app.Microstructure = im2double(app.Microstructure); % rescaling or offsetting the data as necessary.
            elseif strcmp(new_data_type,'double precision array (rescale)')
                app.Microstructure = im2double(app.Microstructure); % rescaling or offsetting the data as necessary.
            end
            app.update_Convert_table(3);
            app.Do_elapsedtime = toc;
            app.DoUndoSave_button('Convert','-','do'); % Update button
            
        end

        % Button pushed function: Convert_UndoButton
        function Convert_UndoButtonPushed(app, event)
            app.update_Convert_table(-3);
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('Convert','-','undo'); % Update button
        end

        % Button pushed function: Convert_SaveButton
        function Convert_SaveButtonPushed(app, event)
            app.update_Convert_table(2); app.update_Convert_table(-3);
            app.DoUndoSave_button('Convert','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Data type conversion'];
            app.History_parameters = [app.History_parameters; app.Convert_Choice_DropDown.Value];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
        end

        % Menu selected function: 
        % histogram_standard_wholevolume_1binpervalue
        function histogram_standard_wholevolume_1binpervalueMenuSelected(app, event)
            options.bin = 'one_per_uniquevalue';
            app.create_histogram(options);
            
        end

        % Button pushed function: CropbackgroundButton
        function CropbackgroundButtonPushed(app, event)
            domain_size=size(app.Microstructure);
            if strcmp(app.ROI_background_OptionsDropDown.Value,'Cylindrical FOV')
                % Binary slice of background
                slice_=app.Microstructure(:,:,1);
                binary_slice=zeros(domain_size(1),domain_size(2));
                binary_slice(slice_==app.BackgroundvalueEditField.Value)=1;
                
                % identify background edge
                distance_transform=bwdist(binary_slice,'euclidean');
                distance_transform(distance_transform>2)=0;
                distance_transform(binary_slice==1)=0;
                distance_transform(distance_transform~=0)=1;
                
                % Find largest rectangle within this edge
                max_area=0;
                fov_x_min=[]; fov_x_max=[]; fov_y_min=[]; fov_y_max=[];
                for min_x=round(0.05*domain_size(1)):1:round(0.45*domain_size(1))
                    % line y
                    line_y=distance_transform(min_x,:);
                    index_line_y=find(line_y==1);
                    if ~isempty(index_line_y) && length(index_line_y>1)
                        min_y=min(index_line_y);
                        max_y=max(index_line_y);
                        
                        % Lines x
                        line_x_1=distance_transform(:,min_y);
                        line_x_2=distance_transform(:,max_y);
                        index_line_x_1=find(line_x_1==1);
                        index_line_x_2=find(line_x_2==1);
                        if ~isempty(index_line_x_1) && length(index_line_x_1>1) && ~isempty(index_line_x_2>1) && length(index_line_x_2>1)
                            max_x=min(max(index_line_x_1),max(index_line_x_2));
                            % Area
                            area=(max_x-min_x)*(max_y-min_y);
                            % Keep only maximum in memory
                            if area>max_area
                                fov_x_min=min_x; fov_x_max=max_x;
                                fov_y_min=min_y; fov_y_max=max_y;
                                max_area=area;
                            end
                        end
                    end
                end
                if ~isempty(fov_x_min) && ~isempty(fov_x_max) && ~isempty(fov_y_min) && ~isempty(fov_y_max)
                    % Slight crop
                    fov_x_min=fov_x_min+2; fov_x_max=fov_x_max-2;
                    fov_y_min=fov_y_min+2; fov_y_max=fov_y_max-2;
                    
                    % Update table
                    tmp = app.ROI_table_crop.Data.Variables;
                    tmp(1,2) = fov_x_min;  tmp(1,3) = fov_x_max;
                    tmp(2,2) = fov_y_min;  tmp(2,3) = fov_y_max;
                    app.ROI_table_crop.Data.Variables = tmp;
                    % 'Click' on table
                    app.ROI_table_cropCellEdit
                end
                
            elseif strcmp(app.ROI_background_OptionsDropDown.Value,'After rotation')
                % Loop over all slices
                x_min=1;
                for x=1:1:domain_size(1)
                    % Extract the current slice
                    slice_=app.Microstructure(x,round(domain_size(2)/4):round(3*domain_size(2)/4),:);
                    % Check background values
                    isanybackground = any(slice_(:)==app.BackgroundvalueEditField.Value);
                    if isanybackground
                        % Update x_min
                        x_min=x+1;
                    else
                        break % Exit loop
                    end
                end
                x_max=domain_size(1);
                for x=domain_size(1):-1:1
                    % Extract the current slice
                    slice_=app.Microstructure(x,round(domain_size(2)/4):round(3*domain_size(2)/4),:);
                    % Check background values
                    isanybackground = any(slice_(:)==app.BackgroundvalueEditField.Value);
                    if isanybackground
                        % Update x_max
                        x_max=x-1;
                    else
                        break % Exit loop
                    end
                end
                y_min=1;
                for y=1:1:domain_size(2)
                    % Extract the current slice
                    slice_=app.Microstructure(round(domain_size(1)/4):round(3*domain_size(1)/4),y,:);
                    % Check background values
                    isanybackground = any(slice_(:)==app.BackgroundvalueEditField.Value);
                    if isanybackground
                        % Update y_min
                        y_min=y+1;
                    else
                        break % Exit loop
                    end
                end
                y_max=domain_size(2);
                for y=domain_size(2):-1:1
                    % Extract the current slice
                    slice_=app.Microstructure(round(domain_size(1)/4):round(3*domain_size(1)/4),y,:);
                    % Check background values
                    isanybackground = any(slice_(:)==app.BackgroundvalueEditField.Value);
                    if isanybackground
                        % Update y_max
                        y_max=y-1;
                    else
                        break % Exit loop
                    end
                end
                % Update table
                tmp = app.ROI_table_crop.Data.Variables;
                tmp(1,2) = x_min;  tmp(1,3) = x_max;
                tmp(2,2) = y_min;  tmp(2,3) = y_max;
                app.ROI_table_crop.Data.Variables = tmp;
                % 'Click' on table
                app.ROI_table_cropCellEdit
                
            end
            
            
        end

        % Button pushed function: VieworthogonalslicesButton
        function VieworthogonalslicesButtonPushed(app, event)
            domain_size=size(app.Microstructure);
            col_ = eval(app.ROI_colormap.Value);
            c = [get(0, 'DefaultAxesColorOrder')];
            
            Fig_= figure;
            Fig_.Name= ['Orthogonal slice view, step ' num2str(app.Operation_step)];
            Fig_.Color='white'; % Background colour
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig_, 'Position', scrsz);
            for id_axe=1:1:4 % Iterate over axe
                sub_axes=subplot(2,2,id_axe,'Parent',Fig_);
                hold(sub_axes,'on');
                if id_axe==1
                    slice_ = squeeze(app.Microstructure(round(domain_size(1)/2),:,:));
                    image(slice_(:,:,1),'CDataMapping','scaled','Parent', sub_axes);
                    set(sub_axes,'YDir','normal');
                    title (['Middle slice, view nornal to axe 1, step ' num2str(app.Operation_step)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    xlabel('3rd Axis');
                    ylabel('2nd Axis');
                elseif id_axe==2
                    slice_ = squeeze(app.Microstructure(:,round(domain_size(2)/2),:));
                    image(slice_(:,:,1),'CDataMapping','scaled','Parent', sub_axes);
                    set(sub_axes,'YDir','normal');
                    title (['Middle slice, view nornal to axe 2, step ' num2str(app.Operation_step)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    xlabel('3rd Axis');
                    ylabel('1st Axis');
                elseif id_axe==3
                    slice_ = app.Microstructure(:,:,round(domain_size(3)/2));
                    image(slice_(:,:,1),'CDataMapping','scaled','Parent', sub_axes);
                    set(sub_axes,'YDir','normal');
                    title (['Middle slice, view nornal to axe 3, step ' num2str(app.Operation_step)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    xlabel('2nd Axis');
                    ylabel('1st Axis');
                elseif id_axe==4
                    slice(double(app.Microstructure),round(domain_size(2)/2),round(domain_size(1)/2),round(domain_size(3)/2));
                    plot3([round(domain_size(2)/2) round(domain_size(2)/2)],[0 domain_size(1)],[round(domain_size(3)/2) round(domain_size(3)/2)],'--','Color',c(1,:),'LineWidth',4,'Parent', sub_axes)
                    plot3([0 domain_size(2)],[round(domain_size(1)/2) round(domain_size(1)/2)],[round(domain_size(3)/2) round(domain_size(3)/2)],'--','Color',c(2,:),'LineWidth',4,'Parent', sub_axes)
                    plot3([round(domain_size(2)/2) round(domain_size(2)/2)],[round(domain_size(1)/2) round(domain_size(1)/2)],[0 domain_size(3)],'--','Color',c(3,:),'LineWidth',4,'Parent', sub_axes)
                    xlabel('2nd Axis');
                    ylabel('1st Axis');
                    zlabel('3rd Axis');
                    shading interp
                    grid(sub_axes,'on');
                    set(sub_axes,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
                    title (['Orthogonal slice 3D view, step ' num2str(app.Operation_step)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    view(sub_axes,30,30)
                end
                set(sub_axes,'Colormap',col_);
                axis(sub_axes,'equal');
                axis(sub_axes,'tight');
                colorbar('peer',sub_axes);
                set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                hold(sub_axes,'off');
            end
            sgtitle(Fig_,['Microstructure slice view, step ' num2str(app.Operation_step)],'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(['Microstructure slice view, step ' num2str(app.Operation_step)]);
                function_savefig(Fig_, app.Savefolder, filename);
            end
        end

        % Button pushed function: Savecurrent2DviewButton
        function Savecurrent2DviewButtonPushed(app, event)
            col_ = eval(app.ROI_colormap.Value);
            Fig_= figure;
            Fig_.Color='white'; % Background colour
            ax_ = axes('Parent',Fig_);
            hold(ax_,'on');
            if app.Dimension==3
                slicer_value = round(app.ROI_SliceSlider.Value);
                direction = app.ROI_ViewnormaltoDropDown.UserData;
                Fig_.Name= ['Slice ' num2str(slicer_value) ' normal to ' app.ROI_ViewnormaltoDropDown.Value ', step ' num2str(app.Operation_step)];
                if direction==1
                    slice_ = squeeze(app.Microstructure(slicer_value,:,:));
                    xlabel('3rd Axis');
                    ylabel('2nd Axis');
                elseif direction==2
                    slice_ = squeeze(app.Microstructure(:,slicer_value,:));
                    xlabel('3rd Axis');
                    ylabel('1st Axis');
                elseif direction==3
                    slice_ = app.Microstructure(:,:,slicer_value);
                    xlabel('2nd Axis');
                    ylabel('1st Axis');
                end
            else
                Fig_.Name= ['Slice, step ' num2str(app.Operation_step)];
                slice_ = app.Microstructure(:,:,1);
            end
            set(ax_,'Colormap',col_);
            axis(ax_,'equal');
            axis(ax_,'tight');
            colorbar('peer',ax_);
            set(ax_,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_axe);
            title (Fig_.Name,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_title);
            image(slice_(:,:,1),'CDataMapping','scaled','Parent', ax_);
            set(ax_,'YDir','normal');
            hold(ax_,'off');
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
        end

        % Button pushed function: Convert_UpdatetableButton
        function Convert_UpdatetableButtonPushed(app, event)
            app.update_Convert_table(2)
        end

        % Button pushed function: ROI_background_DoButton
        function ROI_background_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('Background','-','do'); % Update button
            domain_size=size(app.Microstructure_undo);
            if app.Dimension==2
                domain_size(3)=1;
            end
            if strcmp(app.ROI_backgroundetection_DropDown.Value,'Search only first slice of axe 3, then apply to whole volume')
                app.Do_parameters = sprintf('For first slice of axe3, median filtering range %i, standard deviation threshold %1.3f',app.MedianfilterrangeEditField.Value,app.StdthresholdEditField.Value);
                sdImage = stdfilt(app.Microstructure(:,:,1)); % Calculate standard deviation
                
                figure
                imagesc(sdImage)
                
                sdImage_filtered = medfilt2(sdImage,[app.MedianfilterrangeEditField.Value app.MedianfilterrangeEditField.Value]); % Apply mean filter
                index_ = find(sdImage_filtered<app.StdthresholdEditField.Value); % Threshold it and find coordinates
                if ~isempty(index_)
                    background_id = zeros(length(index_),1) + app.SetidentifiedbackgroundvaluetoEditField.Value;
                    for k=1:1:domain_size(3) % Replace slice per slice
                        slice_tmp = app.Microstructure(:,:,k);
                        slice_tmp(index_) = background_id; % Assign to background
                        app.Microstructure(:,:,k) = slice_tmp;
                    end
                end
            elseif strcmp(app.ROI_backgroundetection_DropDown.Value,'Search slice per slice along axe 3')
                app.Do_parameters = sprintf('Slice per slice, median filtering range %i, standard deviation threshold %1.3f',app.MedianfilterrangeEditField.Value,app.StdthresholdEditField.Value);
                for k=1:1:domain_size(3) % Slice per slice
                    sdImage = stdfilt(app.Microstructure(:,:,k)); % Calculate standard deviation
                    sdImage_filtered = medfilt2(sdImage,[app.MedianfilterrangeEditField.Value app.MedianfilterrangeEditField.Value]); % Apply mean filter
                    index_ = find(sdImage_filtered<app.StdthresholdEditField.Value); % Threshold it and find coordinates
                    if ~isempty(index_)
                        background_id = zeros(length(index_),1) + app.SetidentifiedbackgroundvaluetoEditField.Value;
                        slice_tmp = app.Microstructure(:,:,k);
                        slice_tmp(index_) = background_id; % Assign to background
                        app.Microstructure(:,:,k) = slice_tmp;
                    end
                end
            end
            % Update visualization
            slicer_value = round(app.ROI_SliceSlider.Value);
            app.update_visualization_slice_ROI(slicer_value);
            app.Do_elapsedtime = toc;
        end

        % Button pushed function: ROI_background_UndoButton
        function ROI_background_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('Background','-','undo'); % Update button
            slicer_value = round(app.ROI_SliceSlider.Value);
            app.update_visualization_slice_ROI(slicer_value);
        end

        % Button pushed function: ROI_background_SaveButton
        function ROI_background_SaveButtonPushed(app, event)
            app.DoUndoSave_button('Background','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Background detection'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
            
        end

        % Selection changed function: Scale_DatatypeButtonGroup
        function Scale_DatatypeButtonGroupSelectionChanged(app, event)
            if strcmp(app.Scale_DatatypeButtonGroup.SelectedObject.Text,'Grey level')
                app.Scale_BackgroundvalueEditField.Enable='off';
            elseif strcmp(app.Scale_DatatypeButtonGroup.SelectedObject.Text,'Label')
                app.Scale_BackgroundvalueEditField.Enable='on';
            end
            
        end

        % Value changed function: Scale_VoxelunitDropDown
        function Scale_VoxelunitDropDownValueChanged(app, event)
            app.Scale_newvoxelsizeunitLabel.Text = app.Scale_VoxelunitDropDown.Value;
            
        end

        % Value changed function: Scale_InitialvoxelsizeEditField
        function Scale_InitialvoxelsizeEditFieldValueChanged(app, event)
            app.Scale_NewvoxelsizeEditField.Value = app.Scale_InitialvoxelsizeEditField.Value * app.Scale_scalingfactorEditField.Value;
        end

        % Value changed function: Scale_scalingfactorEditField
        function Scale_scalingfactorEditFieldValueChanged(app, event)
            app.Scale_NewvoxelsizeEditField.Value = app.Scale_InitialvoxelsizeEditField.Value * app.Scale_scalingfactorEditField.Value;
            app.update_ROI_crop_table
        end

        % Button pushed function: Scale_DoButton
        function Scale_DoButtonPushed(app, event)
            if app.Scale_scalingfactorEditField.Value~=1
                tic
                app.Microstructure_undo = app.Microstructure;
                app.DoUndoSave_button('Scale','-','do'); % Update button
                
                % Set parameters
                parameters_scaling.scaling_factor = app.Scale_scalingfactorEditField.Value;
                parameters_scaling.label_or_greylevel = app.Scale_DatatypeButtonGroup.SelectedObject.Text;
                parameters_scaling.background = app.Scale_BackgroundvalueEditField.Value;
                % Scale
                before_domain_size=size(app.Microstructure);
                app.Microstructure = function_scaling(app.Microstructure,parameters_scaling);
                app.Do_elapsedtime = toc;
                
                app.initialize_image_scaling; % Update visualization
                app.Scale_SliceselectionSlider.Enable='on';
                app.update_GUI_dimension
                if app.Dimension==3 % Configure slider
                    app.Scale_SliceselectionSlider.Limits = [1 before_domain_size(3)];
                    app.Scale_SliceselectionSlider.Value = round(before_domain_size(3)/2);
                    app.Scale_SliceselectionSlider.MajorTicks = round(linspace(1,before_domain_size(3),10));
                    app.Scale_SliceselectionSlider.MinorTicks = round(linspace(1,before_domain_size(3),19));
                end
                app.update_all_slicer
                app.update_ROI_crop_table
                app.Do_parameters = sprintf('Scaling facor: %1.3f, option %s',parameters_scaling.scaling_factor,parameters_scaling.label_or_greylevel);
            end
        end

        % Button pushed function: Scale_UndoButton
        function Scale_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('Scale','-','undo'); % Update button
            app.update_GUI_dimension
            app.update_all_slicer
            app.update_ROI_crop_table
            % Clean axis
            cla(app.Scaling_UIAxes_before); colorbar(app.Scaling_UIAxes_before,'off');
            cla(app.Scaling_UIAxes_after); colorbar(app.Scaling_UIAxes_after,'off');
            app.Scale_SliceselectionSlider.Enable='off';
        end

        % Button pushed function: Scale_SaveButton
        function Scale_SaveButtonPushed(app, event)
            app.DoUndoSave_button('Scale','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Rescaling'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
            % Clean axis
            cla(app.Scaling_UIAxes_before); colorbar(app.Scaling_UIAxes_before,'off');
            cla(app.Scaling_UIAxes_after); colorbar(app.Scaling_UIAxes_after,'off');
            app.Scale_SliceselectionSlider.Enable='off';
        end

        % Value changed function: Scale_SliceselectionSlider
        function Scale_SliceselectionSliderValueChanged(app, event)
            app.update_image_scaling;           
            %slicer_value = round(app.Scale_SliceselectionSlider.Value);
            %slice_before = app.Microstructure_undo(:,:,slicer_value);
            %slice_new = app.Microstructure(:,:,round(slicer_value/app.Scale_scalingfactorEditField.Value));
            %imshowpair(slice_before,slice_new,'montage','Parent', app.Scale_UIAxes)
            %montage({double(slice_before), double(slice_new)},'Parent', app.Scale_UIAxes)
            %col_ = eval(app.ROI_colormap.Value);
            %set(app.Scale_UIAxes,'YDir','normal','Colormap',col_);
            %axis(app.Scale_UIAxes,'equal');
            %axis(app.Scale_UIAxes,'tight');
        end

        % Menu selected function: asafunctionofpositionMenu
        function asafunctionofpositionMenuSelected(app, event)
            [Fig_] = function_evolution_along_direction(app.Microstructure);
            % Re-format
            Fig_.Name= ['Grey level along directions, step ' num2str(app.Operation_step)];
            allaxes = findobj( get(Fig_,'Children'), '-depth', 1, 'type', 'axes');
            for k=1:1:length(allaxes)
                set(allaxes(k),'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                allaxes(k).Title.FontName = app.Save_choiceFont_DropDown.Value;
                allaxes(k).Title.FontSize = app.fontisze_subplot4_title;
            end
            sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
        end

        % Menu selected function: inaveraged2DmapsMenu
        function TwoDmapMenuSelected(app, event)
            [Fig_] = function_2Dmap_3Darray(app.Microstructure);
            % Re-format
            Fig_.Name= ['Averaged grey level, view normal to axes, step ' num2str(app.Operation_step)];
            allaxes = findobj( get(Fig_,'Children'), '-depth', 1, 'type', 'axes');
            for k=1:1:length(allaxes)
                set(allaxes(k),'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                allaxes(k).Title.FontName = app.Save_choiceFont_DropDown.Value;
                allaxes(k).Title.FontSize = app.fontisze_subplot4_title;
            end
            sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
            
            
        end

        % Menu selected function: histogram_distribution
        function histogram_distributionMenuSelected(app, event)
            % Calculate histogram
            [histogram_microstructure] = function_calculate_histogram(app.Microstructure);
            
            % Calculate probability density function
            parameters.round_value=0;
            parameters.smooth_cumulative_fct=true;
            parameters.minimum_array_length=5;
            parameters.number_point=length(histogram_microstructure(:,1));
            parameters.moving_range = 0;
            parameters.moving_rangeratio = 0.025;
            parameters.enforce_samelength = false;
            parameters.origin = 'symmetrical';
            parameters.boundary_behavior = 'keep origin';
            [results, ~] = Function_probability_density(histogram_microstructure(:,1),histogram_microstructure(:,2),parameters);
            results.unit ='grey level';
            
            % Figure
            parameters_distributionfigure.title=['Grey level distribution, step ' num2str(app.Operation_step)];
            parameters_distributionfigure.figurename = ['Grey level distribution, step ' num2str(app.Operation_step)];
            parameters_distributionfigure.subaxe1_title = ['Grey level: cummulative fucntion'];
            parameters_distributionfigure.subaxe2_title = ['Grey level: probability density function'];
            parameters_distributionfigure.xlabel = 'Value';
            parameters_distributionfigure.figureposition = [100 100 1500 800];
            parameters_distributionfigure.fontname = app.Save_choiceFont_DropDown.Value;
            parameters_distributionfigure.grid = 'on';
            parameters_distributionfigure.minorgrid = 'on';
            parameters_distributionfigure.save=false;
            parameters_distributionfigure.filename = [];
            parameters_distributionfigure.fullpath = [];
            Fig_ = function_probability_distribution_figure(results,parameters_distributionfigure);
            allaxes = findobj( get(Fig_,'Children'), '-depth', 1, 'type', 'axes');
            for k=1:1:length(allaxes)
                set(allaxes(k),'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                allaxes(k).Title.FontName = app.Save_choiceFont_DropDown.Value;
                allaxes(k).Title.FontSize = app.fontisze_subplot4_title;
            end
            sgtitle(Fig_,['Grey level distribution, step ' num2str(app.Operation_step)],'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
        end

        % Button pushed function: Quality_CalculateSeparabilityButton
        function Quality_CalculateSeparabilityButtonPushed(app, event)
            n_class = app.Quality_numberofphase.Value;
            
            % % Whole volume
            % Calculate histogram
            [histogram_microstructure] = function_calculate_histogram(app.Microstructure);
            % Otsu's algorithm applied to the whole volume
            Otsu_result_wholevolume=Function_otsu_algorithm(histogram_microstructure,n_class);
            
            all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:n_class);
            all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
            best_ratio_wholevolume = max(all_ratio_wholevolume);
            idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
            best_threshold_wholevolume = all_threshold_wholevolume(idx,:);

            % % Slice per slice
            domain_size=size(app.Microstructure);
            app.Dimension = length(domain_size);
            if app.Dimension==3
                for direction=1:1:3
                    val(direction).separability = zeros(domain_size(direction),1);
                    val(direction).threshold = zeros(domain_size(direction),n_class-1);
                    for k=1:1:domain_size(direction)
                        if direction==1
                            slice_ = squeeze(app.Microstructure(k,:,:));
                        elseif direction==2
                            slice_ = squeeze(app.Microstructure(:,k,:));
                        elseif direction==3
                            slice_ = app.Microstructure(:,:,k);
                        end
                        [current_histogram] = function_calculate_histogram(slice_);
                        % Otsu's algorithm applied to the whole volume
                        Otsu_result_perslice=Function_otsu_algorithm(current_histogram,n_class);
                        val(direction).separability(k,1) = Otsu_result_perslice.sigmab2_sigmat2;
                        val(direction).threshold(k,:) = Otsu_result_perslice.threshold(2:n_class);
                    end
                end
            end

            Fig_= figure;
            Fig_.Name= ['Separability criterion (Otsu), number of class:' num2str(n_class) ', step ' num2str(app.Operation_step)];
            Fig_.Color='white'; % Background colour
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig_, 'Position', scrsz);
            if app.Dimension==3
                max_id=4; nplot=2;
            else
                max_id=1; nplot=1;
            end
            col_ = [get(0, 'DefaultAxesColorOrder')];
            for id_axe=1:1:max_id % Iterate over axe
                sub_axes=subplot(nplot,nplot,id_axe,'Parent',Fig_);
                hold(sub_axes,'on');
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also                
                if id_axe==1
                    title ('Whole volume','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    if n_class==2
                        plot(all_threshold_wholevolume,all_ratio_wholevolume,'LineWidth',2,'LineStyle','-','Marker','none');
                        xlabel('Global threshold');
                        ylabel('Separability criterion');
                        legend(sub_axes,['Maximum ratio: ' num2str(best_ratio_wholevolume,'%1.3f') ', for threshold ' num2str(best_threshold_wholevolume)],'Location','best');
                    elseif n_class==3
                        x = all_threshold_wholevolume(:,1); y= all_threshold_wholevolume(:,2); z = all_ratio_wholevolume;
                        % https://www.mathworks.com/matlabcentral/fileexchange/5105-making-surface-plots-from-scatter-data
                        tri = delaunay(x,y);
                        h1=trisurf(tri, x, y, z,'LineStyle','none');
                        % Find optimal                       
                        [optimal_threshold_1,optimal_threshold_2,optimal_threshold] = function_optimal_XYZ(x,y,z);
                        %plot3(optimal_threshold_1(:,1),optimal_threshold_1(:,2),optimal_threshold_1(:,3),'Linestyle','--','Linewidth',2,'Color','k');
                        %plot3(optimal_threshold_2(:,2),optimal_threshold_2(:,1),optimal_threshold_2(:,3),'Linestyle','--','Linewidth',2,'Color','k');
                        h2=plot3(optimal_threshold(:,1),optimal_threshold(:,2),optimal_threshold(:,3),'Linestyle','-','Linewidth',2,'Color','r');
                        h3=plot3(best_threshold_wholevolume(1), best_threshold_wholevolume(2), best_ratio_wholevolume,'Linestyle','none','Linewidth',4,'Marker','x','MarkerSize',14,'Color','k');
                        xlabel('First global threshold');
                        ylabel('Second global threshold');
                        zlabel('Separability criterion');
                        axis(sub_axes,'equal');  
                        colorbar
                        h=colorbar(sub_axes);
                        ylabel(h, 'Separability criterion');
                        set(h,'FontSize',app.fontisze_subplot4_axe,'FontName',app.Save_choiceFont_DropDown.Value);
                        legend(sub_axes,[h2,h3],{'Local maximum',['Maximum ratio: ' num2str(best_ratio_wholevolume,'%1.3f') ', for ' num2str(best_threshold_wholevolume(1)) ' and ' num2str(best_threshold_wholevolume(2))]},'Location','best');
                    else
                        str=cell(1,n_class-1);
                        for k=1:1:n_class
                            if k==1
                               str(1,1) ={['Maximum ratio: ' num2str(best_ratio_wholevolume,'%1.3f')]};
                            else
                               str(1,k)={['Threshold ' num2str(k-1,'%i') ' is ' num2str(Otsu_result_wholevolume.threshold(k))]};
                            end
                        end
                        t = text(0.1,0.5,str);
                        set(t,'FontSize',app.fontisze_subplot4_axe,'FontName',app.Save_choiceFont_DropDown.Value,'EdgeColor','k','BackgroundColor','w');
                        grid(sub_axes,'off'); % Display grid
                        set(sub_axes,'XMinorGrid','off','YMinorGrid','off'); % Display grid for minor thicks also
                        set(sub_axes,'XColor',[1 1 1],'YColor',[1 1 1]);
                    end
                elseif id_axe==2
                    title (['Slice per slice, along Axis ' num2str(id_axe-1)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    str_legend=cell(1,n_class);
                    yyaxis left
                    for k=1:1:n_class-1
                        plot( [1:1:domain_size(1)], val(1).threshold(:,k),'LineStyle','-','Linewidth',1.5,'Color',col_(k,:));
                        str_legend(1,k)={['Threshold ' num2str(k)]};
                    end
                    ylabel('Threshold');
                    xlabel('Slice');
                    yyaxis right
                    plot( [1:1:domain_size(1)], val(1).separability ,'Color','k','Linewidth',2);
                    set(sub_axes,'YColor','k');
                    str_legend(1,n_class) = {'Separability criterion';};
                    legend(sub_axes,str_legend,'Location','best');
                    ylabel('Separability criterion');
                    xlabel('Slice');
                    
                elseif id_axe==3
                    title (['Slice per slice, along Axis ' num2str(id_axe-1)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    yyaxis left
                    for k=1:1:n_class-1
                        plot( [1:1:domain_size(2)], val(2).threshold(:,k),'LineStyle','-','Linewidth',1.5,'Color',col_(k,:));
                    end
                    ylabel('Threshold');
                    xlabel('Slice');
                    yyaxis right
                    plot( [1:1:domain_size(2)], val(2).separability ,'Color','k','Linewidth',2);
                    set(sub_axes,'YColor','k');
                    legend(sub_axes,str_legend,'Location','best');
                    ylabel('Separability criterion');
                    xlabel('Slice');                    
                    
                elseif id_axe==4
                    title (['Slice per slice, along Axis ' num2str(id_axe-1)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    yyaxis left
                    for k=1:1:n_class-1
                        plot( [1:1:domain_size(3)], val(3).threshold(:,k),'LineStyle','-','Linewidth',1.5,'Color',col_(k,:));
                    end
                    ylabel('Threshold');
                    xlabel('Slice');
                    yyaxis right
                    plot( [1:1:domain_size(3)], val(3).separability ,'Color','k','Linewidth',2);
                    set(sub_axes,'YColor','k');
                    legend(sub_axes,str_legend,'Location','best');
                    ylabel('Separability criterion');
                    xlabel('Slice');                      
                    
                end
                axis(sub_axes,'tight');                

                set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                hold(sub_axes,'off');
            end
            sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
        end

        % Menu selected function: foreachsliceMenu
        function foreachsliceMenuSelected(app, event)
            domain_size=size(app.Microstructure);
            app.Dimension = length(domain_size);
            if app.Dimension==3
                for direction=1:1:3
                    for k=1:1:domain_size(direction)
                        if direction==1
                            slice_ = squeeze(app.Microstructure(k,:,:));
                        elseif direction==2
                            slice_ = squeeze(app.Microstructure(:,k,:));
                        elseif direction==3
                            slice_ = app.Microstructure(:,:,k);
                        end
                        [Allhistogram(direction).position(k).val] = function_calculate_histogram(slice_);
                    end
                end

                Fig_= figure;
                Fig_.Name= ['Grey level histogram per slice, step ' num2str(app.Operation_step)];
                Fig_.Color='white'; % Background colour
                scrsz = get(0,'ScreenSize'); % Screen resolution
                set(Fig_, 'Position', [scrsz(1) scrsz(2) scrsz(3) 3/5*scrsz(4)]);
                for id_axe=1:1:3 % Iterate over axe
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig_);
                    hold(sub_axes,'on');
                    grid(sub_axes,'on'); % Display grid
                    set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                    title (['Slice per slice, along Axis ' num2str(id_axe)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    if id_axe==1
                        for pos=1:1:domain_size(1)
                            x_=Allhistogram(1).position(pos).val(:,1);
                            y_=Allhistogram(1).position(pos).val(:,3);
                            z_=ones(length(x_),1)*pos;
                            h_slice=plot3(x_,y_,z_);
                            % Set colors according to depth
                            red_=pos/domain_size(1);
                            blue_=abs(red_-1);
                            green_=0;
                            set(h_slice,'LineWidth',1,'LineStyle','-','Marker','none','color',[red_ green_ blue_]);
                        end
                    elseif id_axe==2
                        for pos=1:1:domain_size(2)
                            x_=Allhistogram(2).position(pos).val(:,1);
                            y_=Allhistogram(2).position(pos).val(:,3);
                            z_=ones(length(x_),1)*pos;
                            h_slice=plot3(x_,y_,z_);
                            % Set colors according to depth
                            red_=pos/domain_size(2);
                            blue_=abs(red_-1);
                            green_=0;
                            set(h_slice,'LineWidth',1,'LineStyle','-','Marker','none','color',[red_ green_ blue_]);
                        end                        
                    elseif id_axe==3
                        for pos=1:1:domain_size(3)
                            x_=Allhistogram(3).position(pos).val(:,1);
                            y_=Allhistogram(3).position(pos).val(:,3);
                            z_=ones(length(x_),1)*pos;
                            h_slice=plot3(x_,y_,z_);
                            % Set colors according to depth
                            red_=pos/domain_size(3);
                            blue_=abs(red_-1);
                            green_=0;
                            set(h_slice,'LineWidth',1,'LineStyle','-','Marker','none','color',[red_ green_ blue_]);
                        end                        
                    end
                    xlabel('Grey level'); zlabel('Slice');
                    view(180,-25)
                    axis(sub_axes,'tight');
                    set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                    hold(sub_axes,'off');
                end
                sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                if app.Save_checkbox.Value
                    filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                    function_savefig(Fig_, app.Savefolder, filename);
                end
            end
        end

        % Menu selected function: foreachsliceMenu_2
        function foreachsliceMenu_2Selected(app, event)
            % modified_img = img + randn(size(img)) * level;
            % true noise = level;
            % [calculated noise,~] = NoiseLevel(modified_img);
            domain_size=size(app.Microstructure);
            app.Dimension = length(domain_size);
            if app.Dimension==3
                array = im2double(app.Microstructure); % Noise function not supported for integer
                for direction=1:1:3
                    Allnoise(direction).nlevel=zeros(1,domain_size(direction));
                    for k=1:1:domain_size(direction)
                        if direction==1
                            slice_ = squeeze(array(k,:,:));
                        elseif direction==2
                            slice_ = squeeze(array(:,k,:));
                        elseif direction==3
                            slice_ = array(:,:,k);
                        end
                        [nlevel, th] = NoiseLevel(slice_);
                        Allnoise(direction).nlevel(1,k) = nlevel / mean(mean(slice_)); % Normalized
                    end
                end
                
                Fig_= figure;
                Fig_.Name= ['Noise estimation per slice, step ' num2str(app.Operation_step)];
                Fig_.Color='white'; % Background colour
                scrsz = get(0,'ScreenSize'); % Screen resolution
                set(Fig_, 'Position', [scrsz(1) scrsz(2) scrsz(3) 3/5*scrsz(4)]);
                for id_axe=1:1:3 % Iterate over axe
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig_);
                    hold(sub_axes,'on');
                    grid(sub_axes,'on'); % Display grid
                    set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                    title (['Slice per slice, along Axis ' num2str(id_axe)],'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    x = [1:1:domain_size(id_axe)];
                    y = Allnoise(id_axe).nlevel(1,:);
                    plot(x,y,'LineWidth',2);
                    legend(sub_axes,['Mean noise: ' num2str(sum(y)/length(x))],'Location','best');
                    xlabel('Slice'); ylabel('Normalized noise');
                    axis(sub_axes,'tight');
                    set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                    hold(sub_axes,'off');
                end
                sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                if app.Save_checkbox.Value
                    filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                    function_savefig(Fig_, app.Savefolder, filename);
                end
            else
                array = double(app.Microstructure); % Noise function not supported for integer
                [nlevel, th] = NoiseLevel(array);
                noise = nlevel / mean(mean(array)); % Normalized
                fprintf('Image noise is %e',noise);
            end
        end

        % Button pushed function: PhaseAssignment_UpdatetableButton
        function PhaseAssignment_UpdatetableButtonPushed(app, event)
            [Phases] = function_calculate_histogram(app.Microstructure);
            Phases(:,3) = Phases(:,2)/numel(app.Microstructure); % Normalize to obtain volume fraction
            tmp = table(Phases(:,1),Phases(:,3),Phases(:,1));
            app.PhaseAssignment_UITable.Data=tmp;
        end

        % Button pushed function: PhaseAssignment_DoButton
        function PhaseAssignment_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('Reassign','-','do'); % Update button
            tmp = app.PhaseAssignment_UITable.Data;
            old_id = tmp(:,1).Variables;
            new_id = tmp(:,3).Variables;
            app.Do_parameters = [];
            for k=1:1:length(old_id)
                if old_id(k)~=new_id(k)
                    app.Microstructure( app.Microstructure_undo==old_id(k) ) = new_id(k);
                    app.Do_parameters=[app.Do_parameters ' ' num2str(old_id(k)) ' to ' num2str(new_id(k))];
                end
            end  
            app.Do_elapsedtime = toc;
        end

        % Button pushed function: PhaseAssignment_UndoButton
        function PhaseAssignment_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('Reassign','-','undo'); % Update button
        end

        % Button pushed function: PhaseAssignment_SaveButton
        function PhaseAssignment_SaveButtonPushed(app, event)
            app.PhaseAssignment_UpdatetableButtonPushed
            app.DoUndoSave_button('Reassign','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Phase re-assignment'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log
        end

        % Cell edit callback: PhaseAssignment_UITable
        function PhaseAssignment_UITableCellEdit(app, event)
            tmp = app.PhaseAssignment_UITable.Data;
            new_Id = tmp(:,3).Variables;
            idx = find(isnan(new_Id)); % Prevent NaN
            if ~isempty(idx)
                tmp(idx,3) = tmp(idx,1);
                app.PhaseAssignment_UITable.Data = tmp;
            end
        end

        % Button pushed function: Contrast_Burn_DoButton
        function Contrast_Burn_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('Burn','-','do'); % Update button
            
            Threshold_H = app.Contrast_VolumepercentthresholdforthehighvaluesEditField.Value;
            Threshold_L = app.Contrast_VolumepercentthresholdforthelowvaluesEditField.Value;

            % % Threshold will be identified using a Dichotomy approach.
            % Faster than compute histogram.
            max_=max(max(max(app.Microstructure)));
            min_=min(min(min(app.Microstructure)));
            number_voxel = numel(app.Microstructure);            
            max_iteration=100;
            target_error=0.1; % in percent            
            if Threshold_H~=0
                threshold_highervalue=max_;
                error_=1e9; iteration_no=0;
                while error_>target_error && iteration_no<max_iteration % Dichotomy loop
                    iteration_no=iteration_no+1; % Update iteration number
                    current_threshold=min_+round((max_-min_)/2); % Current threshold
                    % Calculate volume
                    volume_ratio = 100*sum(sum(sum(app.Microstructure>current_threshold)))/number_voxel;
                    % Check error
                    rel_error=volume_ratio-Threshold_H;
                    error_=abs(rel_error);
                    % Keep going if error is still too high
                    if error_<target_error
                        % Exit the while loop
                        break
                    else
                        % Update the interval
                        if rel_error>0
                            min_= current_threshold;
                        else
                            max_=current_threshold;
                        end
                    end
                end
                threshold_highervalue=current_threshold;
            end
            
            if Threshold_L~=0
                max_=max(max(max(app.Microstructure)));
                min_=min(min(min(app.Microstructure)));
                threshold_lowervalue=min_;
                error_=1e9; iteration_no=0;
                while error_>target_error && iteration_no<max_iteration % Dichotomy loop
                    iteration_no=iteration_no+1; % Update iteration number
                    current_threshold=min_+round((max_-min_)/2); % Current threshold
                    % Calculate volume
                    volume_ratio = 100*sum(sum(sum(app.Microstructure<current_threshold)))/number_voxel;
                    % Check error
                    rel_error=volume_ratio-Threshold_L;
                    error_=abs(rel_error);
                    % Keep going if error is still too high
                    if error_<target_error
                        % Exit the while loop
                        break
                    else
                        % Update the interval
                        if rel_error>0
                            max_= current_threshold;
                        else
                            min_=current_threshold;
                        end
                    end
                end
                threshold_lowervalue=current_threshold;
            end
            
            app.Do_parameters=[];
            % Apply changes
            if Threshold_H~=0
                app.Microstructure(app.Microstructure>=threshold_highervalue)=threshold_highervalue;
                app.Do_parameters=[app.Do_parameters num2str(Threshold_H,'%1.2f') ' percent of volume with higher values burnt '];
            end
            if Threshold_L~=0
                app.Microstructure(app.Microstructure<=threshold_lowervalue)=threshold_lowervalue;
                app.Do_parameters=[app.Do_parameters num2str(Threshold_H,'%1.2f') ' percent of volume with lower values burnt'];
            end
            if Threshold_H~=0 || Threshold_L~=0
                app.initialize_image_contrastcorrection; % Update visualization
                app.Contrast_SliceselectionSlider.Enable='on';
            end

            app.Do_elapsedtime = toc;
        end

        % Button pushed function: Contrast_Burn_UndoButton
        function Contrast_Burn_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('Burn','-','undo'); % Update button
            cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after); 
            colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');            
            app.Contrast_SliceselectionSlider.Enable='off';
        end

        % Button pushed function: Contrast_Burn_SaveButton
        function Contrast_Burn_SaveButtonPushed(app, event)
            app.DoUndoSave_button('Burn','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Burn image extreme values'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.Contrast_SliceselectionSlider.Enable='off';
            cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after); 
            colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');               
            app.update_history_log 
        end

        % Button pushed function: Contrast_Customrange_DoButton
        function Contrast_Customrange_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('CustomRange','-','do'); % Update button
            
            Custom_greyscale=app.Contrast_CustomrangerescaleusingthisnumberofvalueEditField.Value;
            max_=max(max(max(app.Microstructure)));
            min_=min(min(min(app.Microstructure)));
            initial_delta=max_-min_;
            final_delta=Custom_greyscale-1;
            % Initialise
            domain_size=size(app.Microstructure);
            app.Microstructure=zeros(domain_size);
            % Do the calculation
            app.Microstructure=  round( ((double((app.Microstructure_undo-min_)) ./ double(initial_delta)) .* final_delta)+1 );
            app.Do_parameters=['Number of value: ' num2str(Custom_greyscale,'%i')];
            app.initialize_image_contrastcorrection; % Update visualization
            app.Contrast_SliceselectionSlider.Enable='on';
                        
            app.Do_elapsedtime = toc;
        end

        % Button pushed function: Contrast_Customrange_UndoButton
        function Contrast_Customrange_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('CustomRange','-','undo'); % Update button
            app.Contrast_SliceselectionSlider.Enable='off';
            cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after); 
            colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');               
        end

        % Button pushed function: Contrast_Customrange_SaveButton
        function Contrast_Customrange_SaveButtonPushed(app, event)
            app.DoUndoSave_button('CustomRange','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Rescale image with custom range'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.Contrast_SliceselectionSlider.Enable='off';
            cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after); 
            colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');               
            app.update_history_log             
        end

        % Value changed function: Contrast_SliceselectionSlider
        function Contrast_SliceselectionSliderValueChanged(app, event)
            app.update_image_contrastcorrection;
            
        end

        % Cell edit callback: Contrast_TableWhere
        function Contrast_TableWhereCellEdit(app, event)
            tmp = app.Contrast_TableWhere.Data.Variables;
            domain_size=size(app.Microstructure);
            for k=1:1:app.Dimension
                if isnan(tmp(k,2)) || tmp(k,2)<1 || tmp(k,2)~=round(tmp(k,2))
                    tmp(k,2)=1;
                end
                if isnan(tmp(k,3)) || tmp(k,3)>domain_size(k) || tmp(k,3)~=round(tmp(k,3))
                    tmp(k,3)=domain_size(k);
                end    
                if tmp(k,2)>tmp(k,3)
                    s=tmp(k,3); tmp(k,3)=tmp(k,2); tmp(k,2)=s;
                end
                if tmp(k,2)<1 || tmp(k,2)~=round(tmp(k,2))
                    tmp(k,2)=1;
                end
                if tmp(k,3)>domain_size(k) || tmp(k,3)~=round(tmp(k,3))
                    tmp(k,3)=domain_size(k);
                end
            end
            app.Contrast_TableWhere.Data.Variables = tmp;
        end

        % Value changed function: Contrast_Method_DropDown
        function Contrast_Method_DropDownValueChanged(app, event)
            if strcmp(app.Contrast_Method_DropDown.Value,'Adjust contrast')
                parameter_names = {'low_in';'high_in';'low_out';'high_out';'gamma';'Default'};
                parameter_values = {0.0;1.0;0.0;1.0;1.0;1.0};
                app.Contrast_TableParameters.Tooltip ='J = imadjust(I,[low_in high_in],[low_out high_out],gamma) maps intensity values in I to new values in J such that values between low_in and high_in map to values between low_out and high_out, where gamma specifies the shape of the curve describing the relationship between the values in I and J. If you set the last parameter 1, then J = imadjust(I) saturates the bottom 1% and the top 1% of all pixel values.';
                                
            elseif strcmp(app.Contrast_Method_DropDown.Value,'Histogram equalization')
                parameter_names = {'Number of discrete gray level n'};
                parameter_values = {64};
                app.Contrast_TableParameters.Tooltip ='J=histeq(I,n) transforms the grayscale image I, returning in J an grayscale image with n discrete gray levels. A roughly equal number of pixels is mapped to each of the n levels in J, so that the histogram of J is approximately flat. The histogram of J is flatter when n is much smaller than the number of discrete levels in I.';
                
            elseif strcmp(app.Contrast_Method_DropDown.Value,'Contrast-limited adaptive histogram equalization')
                parameter_names = {'NumTiles A';'NumTiles B';'ClipLimit';'NBins';'Range: 0=full, 1=original';'Distribution: 0=uniform, 1=rayleigh, 2=exponential';'Alpha (only for rayleigh and exponential)'};
                parameter_values = {8;8;0.01;256;0;0;0.4};
                app.Contrast_TableParameters.Tooltip ='J = adapthisteq(I) enhances the contrast of the grayscale image I by transforming the values using contrast-limited adaptive histogram equalization (CLAHE).NumTiles: Number of rectangular contextual regions (tiles) into which adapthisteq divides the image, specified as a 2-element vector of positive integers. ClipLimit: Contrast enhancement limit, specified as a number in the range [0, 1]. Higher limits result in more contrast. NBins: Number of histogram bins used to build a contrast enhancing transformation, specified as a positive integer. Higher values result in greater dynamic range at the cost of slower processing speed. Range: range of the output image data, full format range or min-max of input. Distribution: desired histogram shape (uniform=flat, rayliegh=bell-shaped, exponential=curved). Alpha: Distribution parameter, specified as a nonnegative number';
            end
            app.Contrast_TableParameters.Data = table(parameter_names,parameter_values);
        end

        % Button pushed function: Contrast_Advanced_DoButton
        function Contrast_Advanced_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('AdvContrast','-','do'); % Update button
            
            tmp = app.Contrast_TableParameters.Data.Variables;
            p=cell2mat(tmp(:,2));
            if strcmp(app.Contrast_Method_DropDown.Value,'Adjust contrast')
                app.Do_parameters=sprintf('low_in=%1.3f, high_in=%1.3f, low_out=%1.3f, high_out=%1.3f, gamma=%1.3f, default=%i',p(1),p(2),p(3),p(4),p(5),p(6));
            elseif strcmp(app.Contrast_Method_DropDown.Value,'Histogram equalization')
                app.Do_parameters=['n=' num2str(p)];
            elseif strcmp(app.Contrast_Method_DropDown.Value,'Contrast-limited adaptive histogram equalization')
                if p(5)==0
                    str_range = 'full';
                elseif p(5)==1
                    str_range = 'original';
                end
                if p(6)==0
                    str_distribution = 'uniform';
                elseif p(6)==1
                    str_distribution = 'rayleigh';
                elseif p(6)==2
                    str_distribution = 'exponential';
                end
                app.Do_parameters=sprintf('NumTiles=[%i %i], ClipLimit=%1.3f, Nbins=%i, Range=%s, Distribution=%s, Alpha=%1.3f',p(1),p(2),p(3),p(4),str_range,str_distribution,p(7));
            end

            domain_size = size(app.Microstructure);
            tmp = app.Contrast_TableWhere.Data.Variables;
            if app.Dimension==3
                z0=tmp(3,2); z1=tmp(3,3);
            else
                z0=1; z1=1;
            end
            x0 = tmp(1,2); x1 = tmp(1,3); 
            y0 = tmp(2,2); y1 = tmp(2,3); 
            app.Do_parameters = [app.Do_parameters ' from slices (3rd axis) ' num2str(z0) ' to ' num2str(z1) ', within x:' num2str(x0) '-' num2str(x1) ' and y:' num2str(y0) '-' num2str(y1)];
            
            for z=z0:1:z1
                slice_ = app.Microstructure(x0:x1,y0:y1,z);
                if strcmp(app.Contrast_Method_DropDown.Value,'Adjust contrast')
                    if p(6)==1
                        enhanced_slice = imadjust(slice_);
                    else
                        enhanced_slice = imadjust(slice_,[p(1) p(2)],[p(3) p(4)],p(5));
                    end
                elseif strcmp(app.Contrast_Method_DropDown.Value,'Histogram equalization')
                    enhanced_slice = histeq(slice_,p(1));
                elseif strcmp(app.Contrast_Method_DropDown.Value,'Contrast-limited adaptive histogram equalization')
                    if strcmp(str_distribution,'uniform')
                        enhanced_slice = adapthisteq(slice_,'NumTiles',[p(1) p(2)],'ClipLimit',p(3),'NBins',p(4),'Range',str_range,'Distribution',str_distribution);
                    else
                        enhanced_slice = adapthisteq(slice_,'NumTiles',[p(1) p(2)],'ClipLimit',p(3),'NBins',p(4),'Range',str_range,'Distribution',str_distribution,'Alpha',p(7));
                    end
                end
                app.Microstructure(x0:x1,y0:y1,z) = enhanced_slice(:,:,1);
            end
            app.initialize_image_contrastcorrection; % Update visualization
            app.Contrast_SliceselectionSlider.Enable='on';            
            app.Contrast_Advanced_DoButton.UserData = app.Contrast_Method_DropDown.Value;
            app.Do_elapsedtime = toc;
        end

        % Button pushed function: Contrast_Advanced_UndoButton
        function Contrast_Advanced_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('AdvContrast','-','undo'); % Update button
            app.Contrast_SliceselectionSlider.Enable='off';
            cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after); 
            colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');          
        end

        % Button pushed function: Contrast_Advanced_SaveButton
        function Contrast_Advanced_SaveButtonPushed(app, event)
            app.DoUndoSave_button('AdvContrast','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; app.Contrast_Advanced_DoButton.UserData];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.Contrast_SliceselectionSlider.Enable='off';
            cla(app.Contrast_UIAxes_before); cla(app.Contrast_UIAxes_after); 
            colorbar(app.Contrast_UIAxes_before,'off'); colorbar(app.Contrast_UIAxes_after,'off');               
            app.update_history_log         
        end

        % Cell edit callback: Contrast_TableParameters
        function Contrast_TableParametersCellEdit(app, event)
            tmp = app.Contrast_TableParameters.Data.Variables;
            parameters=cell2mat(tmp(:,2));
            if strcmp(app.Contrast_Method_DropDown.Value,'Adjust contrast')
                if isnan(parameters(1)) || parameters(1)<0 || parameters(1)>1
                    parameters(1)=0.0; % Default value
                end
                if isnan(parameters(2)) || parameters(2)<0 || parameters(2)>1
                    parameters(2)=1.0; % Default value
                end
                if parameters(1)>parameters(2)
                    s=parameters(1); parameters(1)=parameters(2); parameters(1)=s;
                end
                if isnan(parameters(3)) || parameters(3)<0 || parameters(3)>1
                    parameters(3)=0.0; % Default value
                end
                if isnan(parameters(4)) || parameters(4)<0 || parameters(4)>1
                    parameters(4)=1.0; % Default value
                end
                if parameters(3)>parameters(4)
                    s=parameters(3); parameters(3)=parameters(4); parameters(4)=s;
                end
                if isnan(parameters(5)) || parameters(5)<0
                    parameters(5)=1.0; % Default value
                end
            elseif strcmp(app.Contrast_Method_DropDown.Value,'Histogram equalization')
                if isnan(parameters(1)) || parameters(1)<1
                    parameters(1)=64; % Default value
                end
            elseif strcmp(app.Contrast_Method_DropDown.Value,'Contrast-limited adaptive histogram equalization')
               
                
            end
            for k=1:1:length(parameters)
                tmp(k,2)={parameters(k)};
            end
            app.Contrast_TableParameters.Data.Variables = tmp;
        end

        % Value changed function: Sensitivity_ExpectedporosityEditField, 
        % Sensitivity_HigherboundEditField, Sensitivity_LowerboundEditField
        function Sensitivity_LowerboundEditFieldValueChanged(app, event)
            expected_porosity = app.Sensitivity_ExpectedporosityEditField.Value;
            lower_bound = app.Sensitivity_LowerboundEditField.Value;
            higher_bound = app.Sensitivity_HigherboundEditField.Value;
            if lower_bound>higher_bound
                app.Sensitivity_LowerboundEditField.Value=higher_bound;
                app.Sensitivity_HigherboundEditField.Value=lower_bound;
            end
            lower_bound = app.Sensitivity_LowerboundEditField.Value;
            higher_bound = app.Sensitivity_HigherboundEditField.Value;            
            if lower_bound>=expected_porosity
                app.Sensitivity_LowerboundEditField.Value=max([expected_porosity-0.05 0]);
            end
            if higher_bound<=expected_porosity            
                app.Sensitivity_HigherboundEditField.Value=min([expected_porosity+0.05 1]);
            end
            
            if app.Sensitivity_LowerboundEditField.Value==app.Sensitivity_ExpectedporosityEditField.Value || app.Sensitivity_LowerboundEditField.Value==app.Sensitivity_HigherboundEditField.Value || app.Sensitivity_HigherboundEditField.Value==app.Sensitivity_ExpectedporosityEditField.Value
                app.Sensitivity_CalculateparametersensitivityButton.Enable='off';
            else
                app.Sensitivity_CalculateparametersensitivityButton.Enable='on';
            end
            
        end

        % Button pushed function: 
        % Sensitivity_CalculateparametersensitivityButton
        function Sensitivity_CalculateparametersensitivityButtonPushed(app, event)
            
            app.Sensitivity_ProgressnotrunningLabel.Text = ['Progress: initializing...'];
            pause(0.1);
            number_voxel=numel(app.Microstructure);
            unique_values=double(unique(app.Microstructure));
            n_=length(unique_values); % Number of unique value
            range_greyscale = max(unique_values)-min(unique_values)+1;
            range_greylevel=100*unique_values/range_greyscale;
            expected_porosity = app.Sensitivity_ExpectedporosityEditField.Value;
            
            % % Porosity
            Volume_fraction_perthreshold=zeros(n_,3); % initialisation
            Volume_fraction_perthreshold(:,1)=unique_values;
            for k=1:1:n_ % Loop over all values
                progress_percent = k*100/n_;
                app.Sensitivity_ProgressnotrunningLabel.Text = ['Progress, porosity step: ' num2str(progress_percent,'%1.3f') '%'];
                current_threshold=unique_values(k); % Select threshold
                Volume_fraction_perthreshold(k,2)=sum(sum(sum(app.Microstructure<=current_threshold)))/number_voxel;
                Volume_fraction_perthreshold(k,3)=1-Volume_fraction_perthreshold(k,2);
            end

            a=abs(Volume_fraction_perthreshold(:,2)-expected_porosity);
            index_reference_ = find( a==min(a)); index_reference_=index_reference_(1);
            Fig_= figure;
            Fig_.Name= ['Porosity sensitivity with threshold, step ' num2str(app.Operation_step)];
            pause(0.1);
            Fig_.Color='white'; % Background colour
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig_, 'Position', [scrsz(1) scrsz(2) scrsz(3) 3/5*scrsz(4)]);
            for id_axe=1:1:3 % Iterate over axe
                sub_axes=subplot(1,3,id_axe,'Parent',Fig_);
                hold(sub_axes,'on');
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                title ('Volume fraction as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                
                if id_axe==1 || id_axe==2
                    x_=Volume_fraction_perthreshold(:,1);
                    xlabel(sub_axes,'Global threshold');
                else
                    x_=range_greylevel;
                    xlabel(sub_axes,'Grey-scale range (%)');
                end
                ylabel(sub_axes,'Volume fraction');
                % Volume fraction
                plot(x_,Volume_fraction_perthreshold(:,2),'LineWidth',2,'LineStyle','-','Marker','none');
                plot(x_,Volume_fraction_perthreshold(:,3),'LineWidth',2,'LineStyle','-','Marker','none');
                % Reference
                ylim_=ylim;
                plot([x_(index_reference_,1) x_(index_reference_,1)],[ylim_(1) ylim_(2)],'LineWidth',2,'LineStyle','--','Marker','none','color',[0.5 0.5 0.5]);
                ylim([ylim_(1) ylim_(2)]);
                if id_axe==1
                    % Set secondary x-axis
                    ax2 = axes('Position',get(sub_axes,'Position'),...
                        'XAxisLocation','top',...
                        'YAxisLocation','right',...
                        'Color','none',...
                        'XColor','k','YColor','w');
                    % Set coincident axe
                    xlimits1 = get(sub_axes,'XLim');
                    xlim(ax2,[100*xlimits1(1)/range_greyscale 100*xlimits1(2)/range_greyscale])
                    % Label
                    xlabel(ax2,'Grey-scale range (%)');
                    % - Fontname and fontsize
                    set(ax2,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                end
                legend(sub_axes,{'Dark phase','Bright phase'},'Location','best');
                axis(sub_axes,'tight');
                set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                hold(sub_axes,'off');
            end
            sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
            
            % % Threshold bounds
            lower_bound_porosity = app.Sensitivity_LowerboundEditField.Value;
            upper_bound_porosity = app.Sensitivity_HigherboundEditField.Value;
            kk=0;
            for k=1:1:n_
                current_porosity=Volume_fraction_perthreshold(k,2);
                current_threshold=unique_values(k);
                if current_porosity>=lower_bound_porosity && current_porosity<=upper_bound_porosity
                    kk=kk+1;
                    all_threshold(kk)=current_threshold;
                    all_porosity(kk)=current_porosity;
                end
            end
            number_threshold=length(all_threshold);
            
            if number_threshold>2
                % % Select only few thresholds
                max_number_thresholds = app.Sensitivity_MaximumnumberofthresholdevaluatedEditField.Value+1;
                n_=length(all_threshold); % number of thresholds
                if n_>max_number_thresholds % Too many?
                    % Find index reference
                    index_reference = find ( abs(all_porosity-expected_porosity)==min( abs(all_porosity-expected_porosity) ) );
                    % distance in ratio
                    distance_lower=round((index_reference/n_)*max_number_thresholds);
                    distance_higher=round(((n_-index_reference)/n_)*max_number_thresholds);
                    % Get threshold
                    threshold_lower = round(linspace(double(all_threshold(1)),double(all_threshold(index_reference)),distance_lower));
                    threshold_higher = round(linspace(double(all_threshold(index_reference)),double(all_threshold(end)),distance_higher));
                    % restriected threshold
                    restrict_all_threshold=[threshold_lower threshold_higher(2:end)];
                    
                    % Porosity calculation
                    nn_=length(restrict_all_threshold);
                    Volume_fraction_restricted_perthreshold=zeros(nn_,1);
                    for k=1:1:nn_
                        current_threshold=restrict_all_threshold(k);
                        % Count voxels
                        Volume_fraction_restricted_perthreshold(k,1)=sum(sum(sum(app.Microstructure<=current_threshold)))/number_voxel;
                    end
                    % Find index reference
                    restrict_index_reference = find ( abs(Volume_fraction_restricted_perthreshold(:,1)-expected_porosity)==min( abs(Volume_fraction_restricted_perthreshold(:,1)-expected_porosity) ) );
                    restrict_number_threshold=length(restrict_all_threshold);
                else
                    % Assign same values
                    restrict_all_threshold=all_threshold;
                    index_reference = find ( abs(all_porosity-expected_porosity)==min( abs(all_porosity-expected_porosity) ) );
                    restrict_index_reference=index_reference;
                    restrict_number_threshold=length(all_threshold);
                end
                % All threshold distance from referenc
                restrict_all_threshold_distancefromref = double(restrict_all_threshold)-double(restrict_all_threshold(restrict_index_reference));
                % Expressed in grey scale percent
                restrict_all_threshold_distancefromref_scale = 100*restrict_all_threshold_distancefromref/range_greyscale;
                
                                
                % % Calculate error for porosity
                reference_value=all_porosity(index_reference);
                error_porosity=zeros(number_threshold,1);
                for k=1:1:number_threshold
                    error_porosity(k,1)=abs( all_porosity(k) - reference_value); % Absolute error
                    error_porosity(k,2)= error_porosity(k,1)*100/reference_value; % Relative error
                end
                % All threshold distance from reference
                all_threshold_distancefromref = double(all_threshold)-double(all_threshold(index_reference));
                % Expressed in grey scale percent
                all_threshold_distancefromref_scale = 100*all_threshold_distancefromref/range_greyscale;
                
                Fig_= figure;
                Fig_.Name= ['Porosity sensitivity with threshold, error, step ' num2str(app.Operation_step)];
                Fig_.Color='white'; % Background colour
                scrsz = get(0,'ScreenSize'); % Screen resolution
                set(Fig_, 'Position', [scrsz(1) scrsz(2) scrsz(3) 3/5*scrsz(4)]);
                for id_axe=1:1:3 % Iterate over axe
                    sub_axes=subplot(1,3,id_axe,'Parent',Fig_);
                    hold(sub_axes,'on');
                    grid(sub_axes,'on'); % Display grid
                    set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                    title ('Volume fraction error as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                    
                    if id_axe==1 || id_axe==2
                        x_=all_threshold_distancefromref;
                        xlabel(sub_axes,'Global threshold variation from the reference');
                    else
                        x_=all_threshold_distancefromref_scale;
                        xlabel(sub_axes,'Grey-scale range variation from the reference (%)');
                    end
                    ylabel(sub_axes,'Porosity error');
                    % Volume fraction
                    plot(x_,error_porosity(:,1),'LineWidth',2,'LineStyle','-','Marker','none','Color','k');
                    if id_axe==1
                        % Set secondary x-axis
                        ax2 = axes('Position',get(sub_axes,'Position'),...
                            'XAxisLocation','top',...
                            'YAxisLocation','right',...
                            'Color','none',...
                            'XColor','k','YColor','w');
                        % Set coincident axe
                        xlimits1 = get(sub_axes,'XLim');
                        xlim(ax2,[100*xlimits1(1)/range_greyscale 100*xlimits1(2)/range_greyscale])
                        % Label
                        xlabel(ax2,'Grey-scale range variation from the reference (%)');
                        % - Fontname and fontsize
                        set(ax2,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                    end
                    legend(sub_axes,{'Porosity'},'Location','best');
                    axis(sub_axes,'tight');
                    set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                    hold(sub_axes,'off');
                end
                sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                if app.Save_checkbox.Value
                    filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                    function_savefig(Fig_, app.Savefolder, filename);
                end
 
                % % Specific surface area
                voxel_size = app.Sensitivity_VoxelsizenanometersEditField.Value;
                if app.Sensitivity_SpecificsurfaceareaCheckBox.Value
                    all_specific_surface_area=zeros(restrict_number_threshold,1); % Initialise
                    for k=1:1:restrict_number_threshold
                        progress_percent = k*100/restrict_number_threshold;
                        app.Sensitivity_ProgressnotrunningLabel.Text = ['Progress, specific surface area step: ' num2str(progress_percent,'%1.3f') '%'];
                        pause(0.01)
                        current_threshold=restrict_all_threshold(k); % Select threshold
                        % binary_phase
                        domain_size=size(app.Microstructure);
                        binary_phase=zeros(domain_size); 
                        binary_phase(app.Microstructure>current_threshold)=1;
                        volume_domain=numel(binary_phase); % Domain volume
                        volume_domain=volume_domain*(voxel_size^3); % Domain volume in cubic micrometer
                        [~,surface_] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
                        surface_=surface_*(voxel_size^2); % Surface area in square micrometer
                        surface_=surface_/volume_domain; % Specific surface area in micrometer-1
                        surface_=surface_*2/3; % Specific surface area in micrometer-1 modified with the corrective factor
                        all_specific_surface_area(k,1)=surface_; % Save it
                    end
                    % Calculate error
                    reference_value=all_specific_surface_area(restrict_index_reference,1);
                    error_specific_surface_area=zeros(restrict_number_threshold,1);
                    for k=1:1:restrict_number_threshold
                        error_specific_surface_area(k,1)=abs( all_specific_surface_area(k,1) - reference_value); % Absolute error
                        error_specific_surface_area(k,2)= error_specific_surface_area(k,1)*100/reference_value; % Relative error
                    end
                    
                    % Figure
                    Fig_= figure;
                    Fig_.Name= ['Specific surface area sensitivity with threshold, step ' num2str(app.Operation_step)];
                    pause(0.1);
                    Fig_.Color='white'; % Background colour
                    scrsz = get(0,'ScreenSize'); % Screen resolution
                    set(Fig_, 'Position', scrsz);
                    for id_axe=1:1:6 % Iterate over axe
                        sub_axes=subplot(2,3,id_axe,'Parent',Fig_);
                        hold(sub_axes,'on');
                        grid(sub_axes,'on'); % Display grid
                        set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                        if id_axe==1 || id_axe==2 || id_axe==3
                            title ('Specific surface area as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            ylabel(sub_axes,'Specifc surface area (\mum^{-1})');
                            y_ = all_specific_surface_area;
                        else
                            title ('Specific surface area error as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            ylabel(sub_axes,'Specifc surface area error (\mum^{-1})');
                            y_ =error_specific_surface_area(:,1);
                        end
                        if id_axe==1 || id_axe==2
                            x_=restrict_all_threshold;
                            xlabel(sub_axes,'Global threshold');
                        elseif id_axe==3
                            x_=100*restrict_all_threshold/range_greyscale;
                            xlabel(sub_axes,'Grey-scale range (%)');
                        elseif id_axe==4 || id_axe==5
                            x_=restrict_all_threshold_distancefromref;
                            xlabel(sub_axes,'Global threshold variation from the reference');
                        elseif id_axe==6
                            x_=restrict_all_threshold_distancefromref_scale;
                            xlabel(sub_axes,'Grey-scale range variation from the reference (%)');
                        end
                        plot(x_,y_,'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','r');
                        % Reference
                        ylim_=ylim;
                        plot([x_(restrict_index_reference) x_(restrict_index_reference)],[ylim_(1) ylim_(2)],'LineWidth',2,'LineStyle','--','Marker','none','color',[0.5 0.5 0.5]);
                        ylim([ylim_(1) ylim_(2)]);
                        if id_axe==1 || id_axe==4
                            % Set secondary x-axis
                            ax2 = axes('Position',get(sub_axes,'Position'),...
                                'XAxisLocation','top',...
                                'YAxisLocation','right',...
                                'Color','none',...
                                'XColor','k','YColor','w');
                            % Set coincident axe
                            xlimits1 = get(sub_axes,'XLim');
                            xlim(ax2,[100*xlimits1(1)/range_greyscale 100*xlimits1(2)/range_greyscale])
                            % Label
                            if id_axe==1
                                xlabel(ax2,'Grey-scale range (%)');
                            else
                                xlabel(ax2,'Grey-scale range variation from the reference (%)');
                            end
                            % - Fontname and fontsize
                            set(ax2,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        end
                        legend(sub_axes,{'Specific surface area'},'Location','best');
                        axis(sub_axes,'tight');
                        set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        hold(sub_axes,'off');
                    end
                    sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                    if app.Save_checkbox.Value
                        filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                        function_savefig(Fig_, app.Savefolder, filename);
                    end
                end
                
                % % Particle diameter
                if app.Sensitivity_ParticlediameterCheckBox.Value
                    all_mean_diameter=zeros(restrict_number_threshold,1); % Initialise
                    for k=1:1:restrict_number_threshold
                        progress_percent = k*100/restrict_number_threshold;
                        app.Sensitivity_ProgressnotrunningLabel.Text = ['Progress, solid particle diameter, step: ' num2str(progress_percent,'%1.3f') '%'];
                        pause(0.01)
                        current_threshold=restrict_all_threshold(k); % Select threshold
                        % binary_phase
                        domain_size=size(app.Microstructure);
                        binary_phase=zeros(domain_size); 
                        binary_phase(app.Microstructure>current_threshold)=1;                        
                        % Call function
                        [Particle_size] = Function_particle_size_CPSD_Algorithm(binary_phase);
                        Particle_size = Particle_size*voxel_size; % um
                        all_diameters = Particle_size(binary_phase==1);                        
                        % Save it
                        all_mean_diameter(k)=mean(all_diameters);
                        
                    end
                    % Calculate error
                    reference_value=all_mean_diameter(restrict_index_reference,1);
                    error_mean_diameter=zeros(restrict_number_threshold,1);
                    for k=1:1:restrict_number_threshold
                        error_mean_diameter(k,1)=abs( all_mean_diameter(k,1) - reference_value);% Absolute error
                        error_mean_diameter(k,2)= error_mean_diameter(k,1)*100/reference_value; % Relative error
                    end
                    % Figure
                    Fig_= figure;
                    Fig_.Name= ['Particle diameter sensitivity with threshold, step ' num2str(app.Operation_step)];
                    pause(0.1);
                    Fig_.Color='white'; % Background colour
                    scrsz = get(0,'ScreenSize'); % Screen resolution
                    set(Fig_, 'Position', scrsz);
                    for id_axe=1:1:6 % Iterate over axe
                        sub_axes=subplot(2,3,id_axe,'Parent',Fig_);
                        hold(sub_axes,'on');
                        grid(sub_axes,'on'); % Display grid
                        set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                        if id_axe==1 || id_axe==2 || id_axe==3
                            title ('Particle diameter as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            ylabel(sub_axes,'Mean equivalent diameter (\mum)');
                            y_ = all_mean_diameter;
                        else
                            title ('Particle diameter error as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            ylabel(sub_axes,'Equivalent diameter mean error (\mum^{-1})');
                            y_ =error_mean_diameter(:,1);
                        end
                        if id_axe==1 || id_axe==2
                            x_=restrict_all_threshold;
                            xlabel(sub_axes,'Global threshold');
                        elseif id_axe==3
                            x_=100*restrict_all_threshold/range_greyscale;
                            xlabel(sub_axes,'Grey-scale range (%)');
                        elseif id_axe==4 || id_axe==5
                            x_=restrict_all_threshold_distancefromref;
                            xlabel(sub_axes,'Global threshold variation from the reference');
                        elseif id_axe==6
                            x_=restrict_all_threshold_distancefromref_scale;
                            xlabel(sub_axes,'Grey-scale range variation from the reference (%)');
                        end
                        plot(x_,y_,'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','b');
                        % Reference
                        ylim_=ylim;
                        plot([x_(restrict_index_reference) x_(restrict_index_reference)],[ylim_(1) ylim_(2)],'LineWidth',2,'LineStyle','--','Marker','none','color',[0.5 0.5 0.5]);
                        ylim([ylim_(1) ylim_(2)]);
                        if id_axe==1 || id_axe==4
                            % Set secondary x-axis
                            ax2 = axes('Position',get(sub_axes,'Position'),...
                                'XAxisLocation','top',...
                                'YAxisLocation','right',...
                                'Color','none',...
                                'XColor','k','YColor','w');
                            % Set coincident axe
                            xlimits1 = get(sub_axes,'XLim');
                            xlim(ax2,[100*xlimits1(1)/range_greyscale 100*xlimits1(2)/range_greyscale])
                            % Label
                            if id_axe==1
                                xlabel(ax2,'Grey-scale range (%)');
                            else
                                xlabel(ax2,'Grey-scale range variation from the reference (%)');
                            end
                            % - Fontname and fontsize
                            set(ax2,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        end
                        legend(sub_axes,{'Mean particle diameter'},'Location','best');
                        axis(sub_axes,'tight');
                        set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        hold(sub_axes,'off');
                    end
                    sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                    if app.Save_checkbox.Value
                        filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                        function_savefig(Fig_, app.Savefolder, filename);
                    end
                end                    
                    
                    
                if app.Sensitivity_PoretortuosityCheckBox.Value
                    all_tortuosity=zeros(restrict_number_threshold,1); % Initialise
                    for k=1:1:restrict_number_threshold
                        progress_percent = k*100/restrict_number_threshold;
                        app.Sensitivity_ProgressnotrunningLabel.Text = ['Progress, pore tortuosity, step: ' num2str(progress_percent,'%1.3f') '%'];
                        pause(0.01)
                        current_threshold=restrict_all_threshold(k); % Select threshold
                        % binary_phase
                        domain_size=size(app.Microstructure);
                        binary_phase=zeros(domain_size); 
                        binary_phase(app.Microstructure<=current_threshold)=1;                        
                        % Call function
                        Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                        % Save it
                        all_tortuosity(k)=Tau_factor_result.Tau_W3.Tau;
                    end
                    % Calculate error
                    reference_value=all_tortuosity(restrict_index_reference,1);
                    error_tortuosity=zeros(restrict_number_threshold,1);
                    for k=1:1:restrict_number_threshold
                        error_tortuosity(k,1)=abs( all_tortuosity(k,1) - reference_value);% Absolute error
                        error_tortuosity(k,2)= error_tortuosity(k,1)*100/reference_value; % Relative error
                    end
                    % Figure
                    Fig_= figure;
                    Fig_.Name= ['Pore tortuosity factor sensitivity with threshold, step ' num2str(app.Operation_step)];
                    pause(0.1);
                    Fig_.Color='white'; % Background colour
                    scrsz = get(0,'ScreenSize'); % Screen resolution
                    set(Fig_, 'Position', scrsz);
                    for id_axe=1:1:6 % Iterate over axe
                        sub_axes=subplot(2,3,id_axe,'Parent',Fig_);
                        hold(sub_axes,'on');
                        grid(sub_axes,'on'); % Display grid
                        set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                        if id_axe==1 || id_axe==2 || id_axe==3
                            title ('Pore tortuosity as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            ylabel(sub_axes,'Pore tortuosity factor \tau_{3rd axis}');
                            y_ = all_tortuosity;
                        else
                            title ('Pore tortuosity error as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                            ylabel(sub_axes,'Pore tortuosity factor \tau_{3rd axis} error');
                            y_ =error_tortuosity(:,1);
                        end
                        if id_axe==1 || id_axe==2
                            x_=restrict_all_threshold;
                            xlabel(sub_axes,'Global threshold');
                        elseif id_axe==3
                            x_=100*restrict_all_threshold/range_greyscale;
                            xlabel(sub_axes,'Grey-scale range (%)');
                        elseif id_axe==4 || id_axe==5
                            x_=restrict_all_threshold_distancefromref;
                            xlabel(sub_axes,'Global threshold variation from the reference');
                        elseif id_axe==6
                            x_=restrict_all_threshold_distancefromref_scale;
                            xlabel(sub_axes,'Grey-scale range variation from the reference (%)');
                        end
                        plot(x_,y_,'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color',[0.929 0.690 0.129]);
                        % Reference
                        ylim_=ylim;
                        plot([x_(restrict_index_reference) x_(restrict_index_reference)],[ylim_(1) ylim_(2)],'LineWidth',2,'LineStyle','--','Marker','none','color',[0.5 0.5 0.5]);
                        ylim([ylim_(1) ylim_(2)]);
                        if id_axe==1 || id_axe==4
                            % Set secondary x-axis
                            ax2 = axes('Position',get(sub_axes,'Position'),...
                                'XAxisLocation','top',...
                                'YAxisLocation','right',...
                                'Color','none',...
                                'XColor','k','YColor','w');
                            % Set coincident axe
                            xlimits1 = get(sub_axes,'XLim');
                            xlim(ax2,[100*xlimits1(1)/range_greyscale 100*xlimits1(2)/range_greyscale])
                            % Label
                            if id_axe==1
                                xlabel(ax2,'Grey-scale range (%)');
                            else
                                xlabel(ax2,'Grey-scale range variation from the reference (%)');
                            end
                            % - Fontname and fontsize
                            set(ax2,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        end
                        legend(sub_axes,{'Pore tortuosity factor'},'Location','best');
                        axis(sub_axes,'tight');
                        set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                        hold(sub_axes,'off');
                    end
                    sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                    if app.Save_checkbox.Value
                        filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                        function_savefig(Fig_, app.Savefolder, filename);
                    end                    
                end
    
                % % Relative error per grey-level
                % Initialise
                slope_x=double(all_threshold)+0.5;
                slope_x(end)=[];
                slope_x_restricted=double(restrict_all_threshold)+0.5;
                slope_x_restricted(end)=[];
                slope_x=slope_x-double(restrict_all_threshold(restrict_index_reference));
                slope_x_restricted=slope_x_restricted-double(restrict_all_threshold(restrict_index_reference));
                slope_x_scale=100*slope_x/range_greyscale;
                slope_x_restricted_scale=100*slope_x_restricted/range_greyscale;
                
                slope_porosity=zeros(number_threshold-1,2);
                slope_specificsurfacearea=zeros(restrict_number_threshold-1,2);
                slope_meandiameter=zeros(restrict_number_threshold-1,2);
                slope_tortuosity=zeros(restrict_number_threshold-1,2);
                
                for k=1:1:number_threshold-1
                    % Abs slope of porosity
                    % Property Error in percent per grey level
                    slope_porosity(k,1)=abs((error_porosity(k+1,2)-error_porosity(k,2))/(double(all_threshold(k+1))-double(all_threshold(k))));
                    % Property error in percent per percentage of grey range
                    slope_porosity(k,2)=slope_porosity(k,1)*range_greyscale/100;
                    
                end
                
                for k=1:1:restrict_number_threshold-1
                    % specific surface area
                    if app.Sensitivity_SpecificsurfaceareaCheckBox.Value
                        slope_specificsurfacearea(k,1)=abs((error_specific_surface_area(k+1,2)-error_specific_surface_area(k,2))/(double(restrict_all_threshold(k+1))-double(restrict_all_threshold(k))));
                        slope_specificsurfacearea(k,2)= slope_specificsurfacearea(k,1)*range_greyscale/100;
                    end
                    % mean diameter
                    if app.Sensitivity_ParticlediameterCheckBox.Value
                        slope_meandiameter(k,1)=abs((error_mean_diameter(k+1,2)-error_mean_diameter(k,2))/(double(restrict_all_threshold(k+1))-double(restrict_all_threshold(k))));
                        slope_meandiameter(k,2)= slope_meandiameter(k,1)*range_greyscale/100;
                    end
                    % tortuosity
                    if app.Sensitivity_PoretortuosityCheckBox.Value
                        slope_tortuosity(k,1)=abs((error_tortuosity(k+1,2)-error_tortuosity(k,2))/(double(restrict_all_threshold(k+1))-double(restrict_all_threshold(k))));
                        slope_tortuosity(k,2)= slope_tortuosity(k,1)*range_greyscale/100;
                    end
                end
                
                % % Reltaive error figure
                Fig_= figure;
                Fig_.Name= ['Relative error sensitivity with threshold, step ' num2str(app.Operation_step)];
                pause(0.1);
                Fig_.Color='white'; % Background colour
                scrsz = get(0,'ScreenSize'); % Screen resolution
                set(Fig_, 'Position', scrsz);
                str_legend={'Porosity'};
                if app.Sensitivity_SpecificsurfaceareaCheckBox.Value
                    str_legend = [str_legend {'Specific surface area'}];
                end
                if app.Sensitivity_ParticlediameterCheckBox.Value
                    str_legend = [str_legend {'Mean particle diameter'}];
                end
                if app.Sensitivity_PoretortuosityCheckBox.Value
                    str_legend = [str_legend {'Pore tortuosity'}];
                end
                for id_axe=1:1:6 % Iterate over axe
                    sub_axes=subplot(2,3,id_axe,'Parent',Fig_);
                    hold(sub_axes,'on');
                    grid(sub_axes,'on'); % Display grid
                    set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                    if id_axe==1 || id_axe==2 || id_axe==3
                        title ('Relative error as a function of global threshold','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                        ylabel(sub_axes,'Relative error (%)');
                    else
                        title ('Relative error per grey-level','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                        ylabel(sub_axes,'Relative error (%) per grey-level change');
                    end
                    if id_axe==1 || id_axe==2
                        x_=restrict_all_threshold_distancefromref;
                        plot(all_threshold_distancefromref,error_porosity(:,2),'LineWidth',2,'LineStyle','-','Marker','none','color','k');
                        xlabel(sub_axes,'Global threshold variation from the reference');
                    elseif id_axe==3
                        x_=restrict_all_threshold_distancefromref_scale;
                        plot(all_threshold_distancefromref_scale,error_porosity(:,2),'LineWidth',2,'LineStyle','-','Marker','none','color','k');
                        xlabel(sub_axes,'Grey-scale range variation from the reference (%)');
                    elseif id_axe==4 || id_axe==5
                        x_=slope_x_restricted;
                        plot(slope_x,slope_porosity(:,1),'LineWidth',2,'LineStyle','-','Marker','none','color','k');
                        xlabel(sub_axes,'Global threshold variation from the reference');
                    elseif id_axe==6
                        x_=slope_x_restricted_scale;
                        plot(slope_x_scale,slope_porosity(:,1),'LineWidth',2,'LineStyle','-','Marker','none','color','k');
                        xlabel(sub_axes,'Grey-scale range variation from the reference (%)');
                    end
                    if id_axe==1 || id_axe==2 || id_axe==3
                        if app.Sensitivity_SpecificsurfaceareaCheckBox.Value
                            plot(x_,error_specific_surface_area(:,2),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','r');
                        end
                        if app.Sensitivity_ParticlediameterCheckBox.Value
                            plot(x_,error_mean_diameter(:,2),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','b');
                        end
                        if app.Sensitivity_PoretortuosityCheckBox.Value
                            plot(x_,error_tortuosity(:,2),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color',[0.929 0.690 0.129]);
                        end
                    else
                        if app.Sensitivity_SpecificsurfaceareaCheckBox.Value
                            plot(x_,slope_specificsurfacearea(:,1),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','r');
                        end
                        if app.Sensitivity_ParticlediameterCheckBox.Value
                            plot(x_,slope_meandiameter(:,1),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','b');
                        end
                        if app.Sensitivity_PoretortuosityCheckBox.Value
                            plot(x_,slope_tortuosity(:,1),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color',[0.929 0.690 0.129]);
                        end
                    end
                    % Reference
                    ylim_=ylim;
                    plot([x_(restrict_index_reference) x_(restrict_index_reference)],[ylim_(1) ylim_(2)],'LineWidth',2,'LineStyle','--','Marker','none','color',[0.5 0.5 0.5]);
                    ylim([ylim_(1) ylim_(2)]);
                    if id_axe==1 || id_axe==4
                        % Set secondary x-axis
                        ax2 = axes('Position',get(sub_axes,'Position'),...
                            'XAxisLocation','top',...
                            'YAxisLocation','right',...
                            'Color','none',...
                            'XColor','k','YColor','w');
                        % Set coincident axe
                        xlimits1 = get(sub_axes,'XLim');
                        xlim(ax2,[100*xlimits1(1)/range_greyscale 100*xlimits1(2)/range_greyscale])
                        % Label
                        xlabel(ax2,'Grey-scale range variation from the reference (%)');
                        % - Fontname and fontsize
                        set(ax2,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                    end
                    legend(sub_axes,str_legend,'Location','best');
                    axis(sub_axes,'tight');
                    set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                    hold(sub_axes,'off');
                end
                sgtitle(Fig_,Fig_.Name,'FontWeight','bold','FontSize',app.fontisze_subplot4_title+2,'FontName',app.Save_choiceFont_DropDown.Value);
                if app.Save_checkbox.Value
                    filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                    function_savefig(Fig_, app.Savefolder, filename);
                end
                
                
                % - Create figure
                Fig_= figure;
                Fig_.Name= 'Relative error per percentage of grey-level range';
                Fig_.Color='white'; % Background colour
                set(Fig_, 'Position', [50 50 660 600]);
                axes_ = axes('Parent',Fig_,'Position',[0.13 0.11 0.775 0.736153846153846]);
                hold(axes_,'on');
                title ('Relative error per percentage of grey-level range','FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_title);
                x_=slope_x_restricted_scale;
                plot(slope_x_scale,slope_porosity(:,2),'LineWidth',2,'LineStyle','-','Marker','none','color','k');
                if app.Sensitivity_SpecificsurfaceareaCheckBox.Value
                    plot(x_,slope_specificsurfacearea(:,2),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','r');
                end
                if app.Sensitivity_ParticlediameterCheckBox.Value
                    plot(x_,slope_meandiameter(:,2),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color','b');
                end
                if app.Sensitivity_PoretortuosityCheckBox.Value
                    plot(x_,slope_tortuosity(:,2),'LineWidth',2,'LineStyle','-','Marker','x','Markersize',14,'color',[0.929 0.690 0.129]);
                end
                % Reference
                ylim_=ylim;
                h_reference=plot([x_(restrict_index_reference) x_(restrict_index_reference)],[ylim_(1) ylim_(2)],'LineWidth',2,'LineStyle','--','Marker','none','color',[0.5 0.5 0.5]);
                ylim([ylim_(1) ylim_(2)]);
                legend(axes_,str_legend,'Location','best');
                xlabel(axes_,'Grey-scale range variation from the reference (%)');
                ylabel(axes_,'Relative error (%) per percentage of grey-level change');
                grid(axes_,'on'); % Display grid
                set(axes_,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
                set(axes_,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_oneplot_axe);
                hold(axes_,'off');
                if app.Save_checkbox.Value
                    filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                    function_savefig(Fig_, app.Savefolder, filename);
                end                
            end
            app.Sensitivity_ProgressnotrunningLabel.Text = 'Progress: not running'; 
        end

        % Value changed function: Filtering_AnisotropicFilter_IterEditField
        function Filtering_AnisotropicFilter_IterEditFieldValueChanged(app, event)
            numberiteration = app.Filtering_AnisotropicFilter_IterEditField.Value;
            tmp=zeros(numberiteration,2);
            tmp(:,1)=1:1:numberiteration;
            tmp(:,2)=ones(numberiteration,1)*25;
            app.Filtering_AnisotropicGradient_Table.Data = table(tmp(:,1),tmp(:,2));            
        end

        % Cell edit callback: Filtering_AnisotropicGradient_Table
        function Filtering_AnisotropicGradient_TableCellEdit(app, event)
            tmp = app.Filtering_AnisotropicGradient_Table.Data.Variables;
            tmp(isnan(tmp))=25;
            tmp(tmp<=0)=25;
            app.Filtering_AnisotropicGradient_Table.Data.Variables = tmp;            
        end

        % Button pushed function: Filtering_Anisotropic_DoButton
        function Filtering_Anisotropic_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('Anisotropicfilter','-','do'); % Update button
            
            Choice_2D3D = app.Filtering_AnisotropicFilter_2D3DDropDown.Value;
            Conduction_Method = app.Filtering_AnisotropicFilter_MethodDropDown.Value;
            Connectivity_ = app.Filtering_AnisotropicFilter_ConnDropDown.Value;
            number_iteration = app.Filtering_AnisotropicFilter_IterEditField.Value;
            tmp = app.Filtering_AnisotropicGradient_Table.Data.Variables; Gradient_threshold=tmp(:,2);
            if app.Filtering_Anisotropy_EstimateCheckBox.Value
                app.Do_parameters = ['Number of iteration and gradient threshold: auto estimate'];
            else
                str=sprintf(',%1.3f',Gradient_threshold);
                app.Do_parameters = sprintf('Number of iteration: %i, gradient threshold %s',number_iteration,str(2:end));
            end
            app.Do_parameters = [app.Do_parameters '. Apply ' Choice_2D3D '. Conduction method: ' Conduction_Method '. Connectivity: ' Connectivity_];

            domain_size=size(app.Microstructure);
            if app.Dimension==3
                if app.Filtering_previewCheckBox.Value
                    n=round(app.Filtering_SliceselectionSlider.Value);
                    z0=n; z1=n;
                else
                    z0=1; z1=domain_size(3);
                end
            else
                z0=1; z1=1;
            end
                        
            if strcmp(Choice_2D3D,'Whole volume') && (~app.Filtering_previewCheckBox.Value || app.Dimension==2)
                if app.Filtering_Anisotropy_EstimateCheckBox.Value
                    [Gradient_threshold,number_iteration] = imdiffuseest(app.Microstructure,'Connectivity',Connectivity_,'ConductionMethod',Conduction_Method);
                end
                app.Filtering_ProgressnotrunningLabel.Text = 'Progress: whole volume in one step...';
                pause(0.01)
                app.Microstructure = imdiffusefilt(app.Microstructure,'Connectivity',Connectivity_,'ConductionMethod',Conduction_Method,'GradientThreshold',Gradient_threshold,'NumberOfIterations',number_iteration);
            else
                for k=z0:1:z1
                    app.Filtering_ProgressnotrunningLabel.Text = sprintf('Progress, slice: %i / %i',k,z1);
                    pause(0.01)
                    slice_ = app.Microstructure(:,:,k);
                    if app.Filtering_Anisotropy_EstimateCheckBox.Value
                        [Gradient_threshold,number_iteration] = imdiffuseest(slice_,'Connectivity',Connectivity_,'ConductionMethod',Conduction_Method);
                        %Gradient_threshold
                        %number_iteration
                    end
                    app.Microstructure(:,:,k) = imdiffusefilt(slice_,'Connectivity',Connectivity_,'ConductionMethod',Conduction_Method,'GradientThreshold',Gradient_threshold,'NumberOfIterations',number_iteration);
                end
            end
            app.Filtering_ProgressnotrunningLabel.Text = 'Progress: not running';
            
            app.initialize_image_filter; % Update visualization
            app.Filtering_SliceselectionSlider.Enable='on';            
            app.Do_elapsedtime = toc;            
        end

        % Button pushed function: Filtering_Anisotropic_UndoButton
        function Filtering_Anisotropic_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('Anisotropicfilter','-','undo'); % Update button
            app.Filtering_SliceselectionSlider.Enable='off';
            cla(app.Filtering_UIAxes_before); cla(app.Filetring_UIAxes_after); 
            colorbar(app.Filtering_UIAxes_before,'off'); colorbar(app.Filetring_UIAxes_after,'off');             
        end

        % Button pushed function: Filtering_Anisotropic_SaveButton
        function Filtering_Anisotropic_SaveButtonPushed(app, event)
            app.DoUndoSave_button('Anisotropicfilter','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Image filtering: Anisotropic diffusion filter'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.Filtering_SliceselectionSlider.Enable='off';
            cla(app.Filtering_UIAxes_before); cla(app.Filetring_UIAxes_after); 
            colorbar(app.Filtering_UIAxes_before,'off'); colorbar(app.Filetring_UIAxes_after,'off');               
            app.update_history_log                   
        end

        % Value changed function: Filtering_SliceselectionSlider
        function Filtering_SliceselectionSliderValueChanged(app, event)
            app.update_image_filter;
            
        end

        % Value changed function: 
        % Filetring_NLMF_ComparisonwindowsizeEditField, 
        % Filetring_NLMF_SearchwindowsizeEditField
        function Filetring_NLMF_SearchwindowsizeEditFieldValueChanged(app, event)
            Searchwindowssize = app.Filetring_NLMF_SearchwindowsizeEditField.Value;
            ComparisonWindowSize = app.Filetring_NLMF_ComparisonwindowsizeEditField.Value;
            if rem(Searchwindowssize, 2)==0
                Searchwindowssize=Searchwindowssize+1;
            end
            if rem(ComparisonWindowSize, 2)==0
                ComparisonWindowSize=ComparisonWindowSize+1;
            end            
            if Searchwindowssize<ComparisonWindowSize
                s=Searchwindowssize; Searchwindowssize=ComparisonWindowSize; ComparisonWindowSize=s;
            end
            app.Filetring_NLMF_SearchwindowsizeEditField.Value = Searchwindowssize;
            app.Filetring_NLMF_ComparisonwindowsizeEditField.Value = ComparisonWindowSize;            
        end

        % Value changed function: Filetring_NLMF_EstimateCheckBox
        function Filetring_NLMF_EstimateCheckBoxValueChanged(app, event)
            app.Filetring_NLMF_DegreeofsmoothingEditField.Enable = ~app.Filetring_NLMF_EstimateCheckBox.Value;           
        end

        % Button pushed function: Filtering_NLMF_DoButton
        function Filtering_NLMF_DoButtonPushed(app, event)
            tic
            app.Microstructure_undo = app.Microstructure;
            app.DoUndoSave_button('NLMF','-','do'); % Update button
            
            degreeofsmoothing = app.Filetring_NLMF_DegreeofsmoothingEditField.Value;
            searchwindowssize = app.Filetring_NLMF_SearchwindowsizeEditField.Value;
            comparizonwindowssize = app.Filetring_NLMF_ComparisonwindowsizeEditField.Value;
            
            if app.Filetring_NLMF_EstimateCheckBox.Value
                app.Do_parameters = ['Degree of smoothing: auto estimate'];
            else
                app.Do_parameters = ['Degree of smoothing: ' num2str(app.Filetring_NLMF_DegreeofsmoothingEditField.Value)];
            end
            app.Do_parameters = [app.Do_parameters ', search windows size:' num2str(searchwindowssize,'%i') ', comparison windows size:' num2str(comparizonwindowssize,'%i')];
            
            domain_size=size(app.Microstructure);
            if app.Dimension==3
                if app.Filtering_previewCheckBox.Value
                    n=round(app.Filtering_SliceselectionSlider.Value);
                    z0=n; z1=n;
                else
                    z0=1; z1=domain_size(3);
                end
            else
                z0=1; z1=1;
            end
            for k=z0:1:z1
                app.Filtering_ProgressnotrunningLabel.Text = sprintf('Progress, slice: %i / %i',k,z1);
                pause(0.01)
                slice_ = app.Microstructure(:,:,k);
                if app.Filetring_NLMF_EstimateCheckBox.Value
                    [app.Microstructure(:,:,k), estDos] = imnlmfilt(slice_,'SearchWindowSize',searchwindowssize,'ComparisonWindowSize',comparizonwindowssize);
                    %estDos % Estimated degree of smoothing
                else
                    app.Microstructure(:,:,k) = imnlmfilt(slice_,'SearchWindowSize',searchwindowssize,'ComparisonWindowSize',comparizonwindowssize,'DegreeOfSmoothing',degreeofsmoothing);
                end
            end
            app.Filtering_ProgressnotrunningLabel.Text = 'Progress: not running';
            
            app.initialize_image_filter; % Update visualization
            app.Filtering_SliceselectionSlider.Enable='on';            
            app.Do_elapsedtime = toc;
        end

        % Button pushed function: Filtering_NLMF_UndoButton
        function Filtering_NLMF_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('NLMF','-','undo'); % Update button
            app.Filtering_SliceselectionSlider.Enable='off';
            cla(app.Filtering_UIAxes_before); cla(app.Filetring_UIAxes_after); 
            colorbar(app.Filtering_UIAxes_before,'off'); colorbar(app.Filetring_UIAxes_after,'off');  
        end

        % Button pushed function: Filtering_NLMF_SaveButton
        function Filtering_NLMF_SaveButtonPushed(app, event)
            app.DoUndoSave_button('NLMF','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Image filtering: non local mean filter'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.Filtering_SliceselectionSlider.Enable='off';
            cla(app.Filtering_UIAxes_before); cla(app.Filetring_UIAxes_after); 
            colorbar(app.Filtering_UIAxes_before,'off'); colorbar(app.Filetring_UIAxes_after,'off');               
            app.update_history_log                 
        end

        % Value changed function: Filtering_Anisotropy_EstimateCheckBox
        function Filtering_Anisotropy_EstimateCheckBoxValueChanged(app, event)
            app.Filtering_AnisotropicFilter_IterEditField.Enable = ~app.Filtering_Anisotropy_EstimateCheckBox.Value;  
            if app.Filtering_Anisotropy_EstimateCheckBox.Value
                app.Filtering_AnisotropicGradient_Table.Enable ='off';
            else
                app.Filtering_AnisotropicGradient_Table.Enable ='on';
            end
        end

        % Value changed function: Segmentation_NumberofphaseEditField
        function Segmentation_NumberofphaseEditFieldValueChanged(app, event)
            numberphase = app.Segmentation_NumberofphaseEditField.Value;
            tmp=zeros(numberphase,4);
            tmp(:,1) = [1:1:numberphase]';
            tmp(:,2) = [1:1:numberphase]'-1;
            tmp(:,3) = ones(numberphase,1)*1/numberphase;
            tmp(:,4) = zeros(numberphase,1);
            app.Segmentation_Phase_UITable.Data = table(tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4));     
            
            min_ = double( min(min(min( app.Microstructure ))));
            max_ = double( max(max(max( app.Microstructure ))));
            tmpthrehsold = linspace(min_,max_,numberphase+1);
            tmp=zeros(numberphase,3);
            tmp(:,1) = [1:1:numberphase]';
            tmp(:,2) = [1:1:numberphase]'-1;
            tmp(:,3) = tmpthrehsold(2:end)';
            app.Segmentation_Global_UITable.Data = table(tmp(:,1),tmp(:,2),tmp(:,3));                 
            
        end

        % Cell edit callback: Segmentation_Phase_UITable
        function Segmentation_Phase_UITableCellEdit(app, event)
            numberphase = app.Segmentation_NumberofphaseEditField.Value;
            tmp = app.Segmentation_Phase_UITable.Data.Variables;
            tmp(isnan(tmp))=0;
            if length(tmp(:,2))~=length(unique(tmp(:,2)))
                tmp(:,2) = [1:1:numberphase]'-1;
            end
            if (sum(tmp(:,3)>1) + sum(tmp(:,3)<0))>0
                tmp(:,3) = ones(numberphase,1)*1/numberphase;
            end
            if sum( abs(sum(tmp(:,3))-1) > 0.001  )
                app.Segmentation_Phase_UITable.ForegroundColor = 'r';
            else
                app.Segmentation_Phase_UITable.ForegroundColor = 'k';
            end
            app.Segmentation_Phase_UITable.Data.Variables = tmp;
            
            tmp2 = app.Segmentation_Global_UITable.Data.Variables;
            tmp2(:,2) = tmp(:,2);
            app.Segmentation_Global_UITable.Data.Variables = tmp2;
            
        end

        % Value changed function: Segmentation_Global_MethodDropDown
        function Segmentation_Global_MethodDropDownValueChanged(app, event)
            Method_global = app.Segmentation_Global_MethodDropDown.Value;
            if strcmp(Method_global,'Manual')
                app.Segmentation_Global_UITable.Enable='on';
            else
                app.Segmentation_Global_UITable.Enable='off';
            end
        end

        % Cell edit callback: Segmentation_Global_UITable
        function Segmentation_Global_UITableCellEdit(app, event)
            numberphase = app.Segmentation_NumberofphaseEditField.Value;
            tmp = app.Segmentation_Global_UITable.Data.Variables;
            min_ = double( min(min(min( app.Microstructure ))));
            max_ = double( max(max(max( app.Microstructure ))));
            tmpthrehsold = linspace(min_,max_,numberphase+1);
            if length(tmp(:,3))~=length(unique(tmp(:,3)))
                tmp(:,3) = tmpthrehsold(2:end)';
            end       
            if sum( tmp(:,3)==sort(tmp(:,3)) ) < numberphase
                tmp(:,3) = tmpthrehsold(2:end)';
            end
            if sum( tmp(:,3)>=min_ ) < numberphase
                tmp(:,3) = tmpthrehsold(2:end)';
            end     
            if sum( tmp(:,3)<=max_ ) < numberphase
                tmp(:,3) = tmpthrehsold(2:end)';
            end 
            app.Segmentation_Global_UITable.Data.Variables = tmp;
        end

        % Button pushed function: 
        % Segmentation_CompareOtsuandexpectedvolumefractionsButton
        function Segmentation_CompareOtsuandexpectedvolumefractionsButtonPushed(app, event)
            numberphase = app.Segmentation_NumberofphaseEditField.Value;
            numbervoxel = numel(app.Microstructure);
            min_ = min(min(min( app.Microstructure )));
            max_ = max(min(max( app.Microstructure )));

            % Calculate histogram
            histogram_microstructure = function_calculate_histogram(app.Microstructure);
            
            % % Determined thresholds: expected volume fraction
            ref_expected_porosity = zeros(numberphase-1,1);
            histogram_microstructure(:,3) = histogram_microstructure(:,2)/numbervoxel;
            % Cumulative function
            n=length(histogram_microstructure(:,3));
            cumulative_fct = zeros(n,1);
            cumulative_fct(end,1)=histogram_microstructure(end,3);
            for current_value=n-1:-1:1
                cumulative_fct(current_value,1)=cumulative_fct(current_value+1,1)+histogram_microstructure(current_value,3);
            end
            % Determine thresholds
            tmp = app.Segmentation_Phase_UITable.Data.Variables;
            expected_volume_fraction = tmp(:,3);
            target = 1-expected_volume_fraction(1);
            for k=1:1:numberphase-1
                idx = find( abs(cumulative_fct-target) == min(abs(cumulative_fct-target)) );
                ref_expected_porosity(k)=histogram_microstructure(idx(1)-1,1);
                target = target-expected_volume_fraction(k+1);
            end
            ref_expected_porosity = [min_; ref_expected_porosity; max_];
            
            % % Determine thresholds: Otsu
            ref_Otsu = zeros(numberphase-1,1);
            % Otsu's algorithm applied to the whole volume
            Otsu_result_wholevolume=Function_otsu_algorithm(histogram_microstructure,numberphase);
            % Best threshold
            all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:numberphase);
            all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
            best_ratio_wholevolume = max(all_ratio_wholevolume);
            idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
            ref_Otsu = all_threshold_wholevolume(idx,:);
            ref_Otsu = [min_; ref_Otsu'; max_];

            for current_phase=1:1:numberphase-1
                LB_expected = ref_expected_porosity(current_phase);
                HB_expected = ref_expected_porosity(current_phase+2);
                thresholds_expected = histogram_microstructure(:,1);
                thresholds_expected(thresholds_expected<LB_expected) = [];
                thresholds_expected(thresholds_expected>HB_expected) = [];
                n_expected = length(thresholds_expected);
                vf_ = zeros(n_expected,1);
                for k=1:1:n_expected
                    condition_1 = app.Microstructure > LB_expected;
                    condition_2 = app.Microstructure <= thresholds_expected(k);
                    vf_(k) = sum(sum(sum( condition_1+condition_2==2 ))) / numbervoxel;
                end
                phase(current_phase).expected.thresholds = thresholds_expected;
                phase(current_phase).expected.vf_abs_error = abs(vf_ - expected_volume_fraction(current_phase));
                phase(current_phase).expected.vf_rel_error = 100 * phase(current_phase).expected.vf_abs_error / expected_volume_fraction(current_phase);

                val_ = zeros(size(all_threshold_wholevolume));
                for k=1:1:numberphase-1
                    if k~=current_phase
                        idx = find(all_threshold_wholevolume(:,k)==ref_Otsu(k+1));
                        val_(idx,k)=1;
                    end
                end
                idx = find ( sum(val_,2)==numberphase-2);
                thresholds_Otsu = all_threshold_wholevolume(idx,current_phase);
                ratio_Otsu = all_ratio_wholevolume(idx);
                phase(current_phase).Otsu.thresholds = thresholds_Otsu;
                phase(current_phase).Otsu.abs_error = abs(ratio_Otsu - Otsu_result_wholevolume.sigmab2_sigmat2);
                phase(current_phase).Otsu.rel_error = 100 * phase(current_phase).Otsu.abs_error / Otsu_result_wholevolume.sigmab2_sigmat2;              
            end            
            
            % Figure
            c = [get(0, 'DefaultAxesColorOrder')];
            Fig_= figure;
            Fig_.Name= ['Global thresholding volume fractions: expected volume fractions vs Otsu, step ' num2str(app.Operation_step)];
            Fig_.Color='white'; % Background colour
            scrsz = get(0,'ScreenSize'); % Screen resolution
            set(Fig_, 'Position', scrsz);
            str_legend = cell(1,(numberphase-1)*2); 
            for current_phase=1:1:numberphase-1
                str_legend(1,2*current_phase-1)={sprintf('Phase %i assgined to %i: error with expected volume fraction',current_phase,tmp(current_phase,2))};
                str_legend(1,2*current_phase)={sprintf('Phase %i assgined to %i: error with Otsu''s separability criterion',current_phase,tmp(current_phase,2))};
            end
            for id_axe=1:1:2
                sub_axes=subplot(2,1,id_axe,'Parent',Fig_);
                hold(sub_axes,'on');
                title (Fig_.Name,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_title);
                if id_axe==1
                    ylabel(sub_axes, 'Absolute error');
                    ylim(sub_axes,[0 1]);
                else
                    ylabel(sub_axes, 'Relative error (%)');
                    ylim(sub_axes,[0 100]);
                end
                for current_phase=1:1:numberphase-1
                    if id_axe==1 % Absolute error
                        plot(phase(current_phase).expected.thresholds, phase(current_phase).expected.vf_abs_error,'LineWidth',2,'Color',c(current_phase,:));
                        plot(phase(current_phase).Otsu.thresholds, phase(current_phase).Otsu.abs_error,'LineWidth',2,'LineStyle','--','Color',c(current_phase,:));
                    else % Relative error
                        plot(phase(current_phase).expected.thresholds, phase(current_phase).expected.vf_rel_error,'LineWidth',2,'Color',c(current_phase,:));
                        plot(phase(current_phase).Otsu.thresholds, phase(current_phase).Otsu.rel_error,'LineWidth',2,'LineStyle','--','Color',c(current_phase,:));
                    end
                end   
                legend(sub_axes,str_legend,'Location','best');
                grid(sub_axes,'on');
                set(sub_axes,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
                xlabel(sub_axes, 'Global threshold');
                set(sub_axes,'FontName',app.Save_choiceFont_DropDown.Value,'FontSize',app.fontisze_subplot4_axe);
                hold(sub_axes,'off');
            end
            if app.Save_checkbox.Value
                filename = function_remove_emptyandspecialcharacter_string(Fig_.Name);
                function_savefig(Fig_, app.Savefolder, filename);
            end
            
        end

        % Button pushed function: Segmentation_Global_DoButton
        function Segmentation_Global_DoButtonPushed(app, event)
            tmp = app.Segmentation_Phase_UITable.Data.Variables;
            if ~sum( abs(sum(tmp(:,3))-1) > 0.001  )
                
                tic
                domain_size = size(app.Microstructure);
                app.Microstructure_undo = app.Microstructure;
                app.DoUndoSave_button('SegmentationGlobal','-','do'); % Update button
                
                numberphase = app.Segmentation_NumberofphaseEditField.Value;
                Method_global = app.Segmentation_Global_MethodDropDown.Value;
                app.Do_parameters = Method_global;
                
                % Find thresholds
                if strcmp(Method_global,'Manual')
                    foo=1; %
                    
                elseif strcmp(Method_global,'Match expected volume fraction')
                    histogram_microstructure = function_calculate_histogram(app.Microstructure);
                    histogram_microstructure(:,3) = histogram_microstructure(:,2)/numel(app.Microstructure);
                    % Cumulative function
                    n=length(histogram_microstructure(:,3));
                    cumulative_fct = zeros(n,1);
                    cumulative_fct(end,1)=histogram_microstructure(end,3);
                    for current_value=n-1:-1:1
                        cumulative_fct(current_value,1)=cumulative_fct(current_value+1,1)+histogram_microstructure(current_value,3);
                    end
                    % Determine thresholds
                    tmp = app.Segmentation_Phase_UITable.Data.Variables;
                    tmp2 = app.Segmentation_Global_UITable.Data.Variables;
                    expected_volume_fraction = tmp(:,3);
                    target = 1-expected_volume_fraction(1);
                    for k=1:1:numberphase-1
                        idx = find( abs(cumulative_fct-target) == min(abs(cumulative_fct-target)) );
                        tmp2(k,3)=histogram_microstructure(idx(1)-1,1);
                        target = target-expected_volume_fraction(k+1);
                    end
                    tmp2(end,3) = max(max(max( app.Microstructure )));
                    app.Segmentation_Global_UITable.Data.Variables = tmp2;
                    
                elseif strcmp(Method_global,'Otsu')
                    [histogram_microstructure] = function_calculate_histogram(app.Microstructure);
                    % Otsu's algorithm applied to the whole volume
                    Otsu_result_wholevolume=Function_otsu_algorithm(histogram_microstructure,numberphase);
                    % Best threshold
                    all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:numberphase);
                    all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
                    best_ratio_wholevolume = max(all_ratio_wholevolume);
                    idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
                    best_threshold_wholevolume = all_threshold_wholevolume(idx,:);
                    % Update table
                    tmp = app.Segmentation_Global_UITable.Data.Variables;
                    tmp(:,3) = [best_threshold_wholevolume';  max(max(max( app.Microstructure)))];
                    %tmp(end,3) = max(max(max( app.Microstructure )));
                    app.Segmentation_Global_UITable.Data.Variables = tmp;
                end
                
                % Assign voxels to phase based on threshold
                tmp = app.Segmentation_Global_UITable.Data.Variables;
                thresholds = [-1e9; tmp(:,3)];
                app.Microstructure = zeros(domain_size);
                for k=1:1:numberphase
                    condition_1 = app.Microstructure_undo > thresholds(k);
                    condition_2 = app.Microstructure_undo <= thresholds(k+1);
                    app.Microstructure( condition_1+condition_2==2 ) = tmp(k,2);
                end
                
                % Update volume fraction
                tmp = app.Segmentation_Phase_UITable.Data.Variables;
                number_voxel=numel(app.Microstructure);
                for k=1:1:numberphase
                    tmp(k,4) = sum(sum(sum( app.Microstructure==tmp(k,2) ))) / number_voxel;
                end
                app.Segmentation_Phase_UITable.Data.Variables = tmp;
                
                % Parameters
                tmp = app.Segmentation_Global_UITable.Data.Variables;
                for k=1:1:numberphase
                    app.Do_parameters = [app.Do_parameters ', Phase ' num2str(k) ' assigned to ' num2str(tmp(k,2)) ', threshold <=' num2str(tmp(k,3))];
                end
                app.Do_elapsedtime = toc;
            end
        end

        % Button pushed function: Segmentation_Global_UndoButton
        function Segmentation_Global_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('SegmentationGlobal','-','undo'); % Update button
        end

        % Button pushed function: Segmentation_Global_SaveButton
        function Segmentation_Global_SaveButtonPushed(app, event)
            app.DoUndoSave_button('SegmentationGlobal','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Segmentation: global threshold'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log 
            app.figure_segmentation(0) % Segmentation result
        end

        % Button pushed function: Segmentation_Local_DoButton
        function Segmentation_Local_DoButtonPushed(app, event)
            tmp = app.Segmentation_Phase_UITable.Data.Variables;
            if ~sum( abs(sum(tmp(:,3))-1) > 0.001  )
                tic
                domain_size = size(app.Microstructure);
                app.Microstructure_undo = app.Microstructure;
                app.DoUndoSave_button('SegmentationLocal','-','do'); % Update button
                numberphase = app.Segmentation_NumberofphaseEditField.Value;
                Method_local = app.Segmentation_Local_MethodDropDown.Value;
                
                if strcmp(Method_local,'Local Otsu') || strcmp(Method_local,'Local Otsu translated to match expected volume fractions')
                    % Parameters
                    app.Do_parameters = Method_local;
                    max_ = max(max(max( app.Microstructure )));
                    
                    % Determine thresholds
                    thresholds=zeros(domain_size(3),numberphase+1);
                    for k=1:1:domain_size(3)
                        slice_ = app.Microstructure(:,:,k);
                        [current_histogram] = function_calculate_histogram(slice_);
                        Otsu_result=Function_otsu_algorithm(current_histogram,numberphase);
                        % Best threshold
                        all_threshold=Otsu_result.allpermuation(:,2:numberphase);
                        all_ratio=Otsu_result.allratio;
                        best_ratio = max(all_ratio);
                        idx = find(all_ratio==best_ratio);
                        best_threshold = all_threshold(idx,:);
                        thresholds(k,:) = [-1e9; best_threshold'; max_]';
                    end
                    app.sav_thresholds = thresholds;
                    
                    % Smooth threshold
                    if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                        mean_range = app.Segmentation_MovingaveragefilterrangeEditField.Value;
                        app.Do_parameters = [app.Do_parameters ', threshold smoothed along axis with moving average range:' num2str(mean_range,'%i')];
                        thresholds_smoothed = thresholds;
                        for current_phase=1:1:numberphase
                            t_=thresholds(:,current_phase+1);
                            for z=1:1:domain_size(3)
                                min_z=max(1,z-mean_range);
                                max_z=min(domain_size(3),z+mean_range);
                                thresholds_smoothed(z,current_phase+1)= mean(t_(min_z:max_z));
                            end
                        end
                        app.sav_thresholds_smoothed = thresholds_smoothed;
                    end
                    
                    % Replace slice
                    if strcmp(Method_local,'Local Otsu')
                        for k=1:1:domain_size(3)
                            slice_ = app.Microstructure(:,:,k);
                            slice_segmented = zeros(size(slice_));
                            for current_phase=1:1:numberphase
                                if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                                    condition_1 = double(slice_) > double(thresholds_smoothed(k,current_phase));
                                    condition_2 = double(slice_) <= double(thresholds_smoothed(k,current_phase+1));
                                else
                                    condition_1 = slice_ > thresholds(k,current_phase);
                                    condition_2 = slice_ <= thresholds(k,current_phase+1);
                                end
                                slice_segmented( condition_1+condition_2==2 ) = tmp(current_phase,2);
                            end
                            app.Microstructure(:,:,k) = slice_segmented;
                        end
                        
                    elseif strcmp(Method_local,'Local Otsu translated to match expected volume fractions')
                        % Threshold will be translated to match expected volume fraction
                        if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                            initial_thresholds = app.sav_thresholds_smoothed; % Initialize
                        else
                            initial_thresholds = app.sav_thresholds_smoothed;
                        end
                        translated_thresholds = initial_thresholds;
                        
                        number_voxel=numel(app.Microstructure);
                        tmp = app.Segmentation_Phase_UITable.Data.Variables;
                        min_ = double(min(min(min( app.Microstructure_undo))));
                        max_ = double(max(max(max( app.Microstructure_undo))));
                        for current_phase=1:1:numberphase-1 % No need to check last phase.
                            expected_vf = tmp(current_phase,3);
                            
                            % Dichotomy algorithm
                            error_vf=1e9;
                            max_iteration = 50;
                            min_error = 1e-3;
                            n_iteration=1;
                            L = -(max_-min_); % Lower bound
                            H = max_-min_; % Higher bound
                            threshold_translation = 0;
                            while abs(error_vf)>min_error && n_iteration<max_iteration
                                for pos=1:1:domain_size(3)
                                    slice_ = app.Microstructure_undo(:,:,pos);
                                    slice_segmented = zeros(size(slice_));
                                    condition_1 = double(slice_) > double(initial_thresholds(pos,current_phase));
                                    condition_2 = double(slice_) <= double(initial_thresholds(pos,current_phase+1)) + threshold_translation;
                                    slice_segmented( condition_1+condition_2==2 ) = 1;
                                    app.Microstructure(:,:,pos) = slice_segmented;
                                end
                                obtained_vf = sum(sum(sum( app.Microstructure==1 ))) / number_voxel;
                                error_vf = expected_vf-obtained_vf;
                                if abs(error_vf)>min_error && n_iteration<max_iteration
                                    if error_vf>0 % High threshold must increase
                                        L=threshold_translation; H=H;
                                        threshold_translation = (threshold_translation+H)/2;
                                    else % High threshold must decrease
                                        L=L; H=threshold_translation;
                                        threshold_translation = (threshold_translation+L)/2;
                                    end
                                end
                                n_iteration=n_iteration+1;
                            end
                            translated_thresholds(:,current_phase+1) = double(initial_thresholds(:,current_phase+1)) + threshold_translation;
                        end
                        for k=1:1:domain_size(3)
                            slice_ = app.Microstructure_undo(:,:,k);
                            slice_segmented = zeros(size(slice_));
                            for current_phase=1:1:numberphase
                                if app.Segmentation_SmoothSliceThreshold_CheckBox.Value
                                    condition_1 = double(slice_) > double(translated_thresholds(k,current_phase));
                                    condition_2 = double(slice_) <= double(translated_thresholds(k,current_phase+1));
                                else
                                    condition_1 = slice_ > thresholds(k,current_phase);
                                    condition_2 = slice_ <= thresholds(k,current_phase+1);
                                end
                                slice_segmented( condition_1+condition_2==2 ) = tmp(current_phase,2);
                            end
                            app.Microstructure(:,:,k) = slice_segmented;
                        end
                        app.sav_translated_thresholds = translated_thresholds; % Save
                    end
                end
                
                % Update volume fraction
                tmp = app.Segmentation_Phase_UITable.Data.Variables;
                number_voxel=numel(app.Microstructure);
                for k=1:1:numberphase
                    tmp(k,4) = sum(sum(sum( app.Microstructure==tmp(k,2) ))) / number_voxel;
                end
                app.Segmentation_Phase_UITable.Data.Variables = tmp;                
                
                app.Do_elapsedtime = toc;
            end
        end

        % Button pushed function: Segmentation_Local_UndoButton
        function Segmentation_Local_UndoButtonPushed(app, event)
            app.Microstructure = app.Microstructure_undo;
            app.DoUndoSave_button('SegmentationLocal','-','undo'); % Update button
        end

        % Button pushed function: Segmentation_Local_SaveButton
        function Segmentation_Local_SaveButtonPushed(app, event)
            app.DoUndoSave_button('SegmentationLocal','-','save'); % Update button
            app.Operation_step=app.Operation_step+1;
            app.History_step = [app.History_step; num2str(app.Operation_step)];
            app.History_operations = [app.History_operations; 'Segmentation: local threshold'];
            app.History_parameters = [app.History_parameters; app.Do_parameters];
            app.History_elapsedtime = [app.History_elapsedtime; [num2str(app.Do_elapsedtime,'%1.1f') 's']];
            app.update_history_log             
            app.figure_segmentation(1) % Segmentation result
        end

        % Value changed function: Segmentation_SliceSlider
        function Segmentation_SliceSliderValueChanged(app, event)
            if app.Dimension==3
                slicer_value = round(app.Segmentation_SliceSlider.Value);
            else
                slicer_value=1;
            end
            slice_ = app.Microstructure(:,:,slicer_value);
            
            if strcmp(app.Segmentation_Global_DoButton.Enable,'on') && strcmp(app.Segmentation_Local_DoButton.Enable,'on')
                overlay_=app.Segmentation_OpacityoverlayEditField.Value;
                
                % Between bounds
                min_ = app.Segmentation_Lowthreshold_SliceSlider.Limits(1);
                max_ = app.Segmentation_Lowthreshold_SliceSlider.Limits(2);
                a=(1-0)/(max_-min_); b=1-a*max_;
                low = app.Segmentation_Lowthreshold_SliceSlider.Value;
                high = app.Segmentation_Highthreshold_SliceSlider.Value;
                d=size(slice_);
                between_bounds = zeros(d(1),d(2),3); %RGB
                tmp = zeros(size(slice_));
                condition_1 = slice_ >= low;
                condition_2 = slice_ <= high;
                tmp( condition_1+condition_2==2 ) = 1;
                between_bounds(:,:,1) = tmp;
                
                % Convert grey image in RGB
                slice_in_RGB = cat(3, slice_, slice_, slice_);
                slice_in_RGB = a*double(slice_in_RGB)+b;
                
                % Linear combination (overlay)
                C = overlay_*between_bounds + (1-overlay_)*slice_in_RGB;
                image(C,'CDataMapping','scaled','Parent', app.Segmentation_UIAxes);
                                
                set(app.Segmentation_UIAxes,'YDir','normal','Colormap',gray);
                axis(app.Segmentation_UIAxes,'equal');
                axis(app.Segmentation_UIAxes,'tight');
                
            else
                slice_unsegmented = app.Microstructure_undo(:,:,slicer_value);
                imshowpair(slice_unsegmented,slice_,'montage','Parent', app.Segmentation_UIAxes);
                set(app.Segmentation_UIAxes,'YDir','normal','Colormap',gray);
                axis(app.Segmentation_UIAxes,'equal');
                axis(app.Segmentation_UIAxes,'tight');                
            end
            
            
        end

        % Value changed function: Segmentation_Highthreshold_SliceSlider, 
        % Segmentation_Lowthreshold_SliceSlider
        function Segmentation_Lowthreshold_SliceSliderValueChanged(app, event)
            min_ = double(min(min(min( app.Microstructure ))));
            max_ = double(max(max(max( app.Microstructure ))));            
            low = app.Segmentation_Lowthreshold_SliceSlider.Value;
            high = app.Segmentation_Highthreshold_SliceSlider.Value;
            if low==app.Segmentation_Lowthreshold_SliceSlider.Limits(2)
                low= min_ + 0.99*(max_-min_);
            end
            if high==app.Segmentation_Lowthreshold_SliceSlider.Limits(1)
                high= min_ + 0.01*(max_-min_);
            end            
            if low>=high
                low=high-0.01*(max_-min_);
            end
            app.Segmentation_Lowthreshold_SliceSlider.Value=low;
            app.Segmentation_Highthreshold_SliceSlider.Value=high;
            app.Segmentation_Lowthreshold_SliceSpinner.Value = low;
            app.Segmentation_Highthreshold_SliceSpinner.Value = high;
            app.Segmentation_SliceSliderValueChanged % Update visualization
        end

        % Button pushed function: Segmentation_update_range
        function Segmentation_update_rangeButtonPushed(app, event)
            min_ = double(min(min(min( app.Microstructure ))));
            max_ = double(max(max(max( app.Microstructure ))));
            app.Segmentation_Lowthreshold_SliceSlider.Limits = [min_ max_];
            app.Segmentation_Lowthreshold_SliceSlider.Value = min_;
            app.Segmentation_Lowthreshold_SliceSlider.MajorTicks = linspace(min_,max_,7);
            app.Segmentation_Lowthreshold_SliceSlider.MinorTicks = linspace(min_,max_,15);
            app.Segmentation_Highthreshold_SliceSlider.Limits = [min_ max_];
            app.Segmentation_Highthreshold_SliceSlider.Value = max_;
            app.Segmentation_Highthreshold_SliceSlider.MajorTicks = linspace(min_,max_,7);
            app.Segmentation_Highthreshold_SliceSlider.MinorTicks = linspace(min_,max_,15);
            app.Segmentation_Lowthreshold_SliceSpinner.Limits = [min_ max_];
            app.Segmentation_Lowthreshold_SliceSpinner.Value = min_;
            app.Segmentation_Highthreshold_SliceSpinner.Limits = [min_ max_];
            app.Segmentation_Highthreshold_SliceSpinner.Value = max_;
        end

        % Value changed function: Segmentation_Lowthreshold_SliceSpinner
        function Segmentation_Lowthreshold_SliceSpinnerValueChanged(app, event)
            app.Segmentation_Lowthreshold_SliceSlider.Value = app.Segmentation_Lowthreshold_SliceSpinner.Value;
            app.Segmentation_Lowthreshold_SliceSliderValueChanged
        end

        % Value changed function: Segmentation_Highthreshold_SliceSpinner
        function Segmentation_Highthreshold_SliceSpinnerValueChanged(app, event)
            app.Segmentation_Highthreshold_SliceSlider.Value = app.Segmentation_Highthreshold_SliceSpinner.Value;
            app.Segmentation_Lowthreshold_SliceSliderValueChanged
            
        end

        % Value changed function: 
        % Segmentation_SmoothSliceThreshold_CheckBox
        function Segmentation_SmoothSliceThreshold_CheckBoxValueChanged(app, event)
             app.Segmentation_MovingaveragefilterrangeEditField.Enable = app.Segmentation_SmoothSliceThreshold_CheckBox.Value;            
        end

        % Menu selected function: DocumentationMenu
        function DocumentationMenuSelected(app, event)
            Find_file('NREL_MATBOX_Microstructure_analysis_toolbox_documentation.pdf','MATBOX_Microstructure_analysis_toolbox','Default location is \MATBOX_Microstructure_analysis_toolbox\Documentation\');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ROIfilteringandsegmentationmoduleUIFigure and hide until all components are created
            app.ROIfilteringandsegmentationmoduleUIFigure = uifigure('Visible', 'off');
            app.ROIfilteringandsegmentationmoduleUIFigure.Position = [100 100 930 776];
            app.ROIfilteringandsegmentationmoduleUIFigure.Name = 'ROI, filtering and segmentation module';
            app.ROIfilteringandsegmentationmoduleUIFigure.Icon = 'Icon_ROI.png';

            % Create VolumeMenu
            app.VolumeMenu = uimenu(app.ROIfilteringandsegmentationmoduleUIFigure);
            app.VolumeMenu.Text = 'Volume';

            % Create LoadvolumeMenu
            app.LoadvolumeMenu = uimenu(app.VolumeMenu);
            app.LoadvolumeMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadvolumeMenuSelected, true);
            app.LoadvolumeMenu.Text = 'Load volume';

            % Create SavevolumeMenu
            app.SavevolumeMenu = uimenu(app.VolumeMenu);
            app.SavevolumeMenu.Text = 'Save volume';

            % Create DefaultlocationMenu
            app.DefaultlocationMenu = uimenu(app.SavevolumeMenu);
            app.DefaultlocationMenu.MenuSelectedFcn = createCallbackFcn(app, @SavevolumeMenuSelected, true);
            app.DefaultlocationMenu.Text = 'Default location';

            % Create CustomlocationMenu
            app.CustomlocationMenu = uimenu(app.SavevolumeMenu);
            app.CustomlocationMenu.MenuSelectedFcn = createCallbackFcn(app, @CustomlocationMenuSelected, true);
            app.CustomlocationMenu.Text = 'Custom location';

            % Create MenuCalculate
            app.MenuCalculate = uimenu(app.ROIfilteringandsegmentationmoduleUIFigure);
            app.MenuCalculate.Enable = 'off';
            app.MenuCalculate.Text = 'Calculate';

            % Create PlothistogramMenu
            app.PlothistogramMenu = uimenu(app.MenuCalculate);
            app.PlothistogramMenu.Text = 'Plot histogram';

            % Create forthewholevolumeMenu
            app.forthewholevolumeMenu = uimenu(app.PlothistogramMenu);
            app.forthewholevolumeMenu.Text = 'for the whole volume';

            % Create histogram_standard_wholevolume_auto
            app.histogram_standard_wholevolume_auto = uimenu(app.forthewholevolumeMenu);
            app.histogram_standard_wholevolume_auto.MenuSelectedFcn = createCallbackFcn(app, @histogram_standard_wholevolume_autoSelected, true);
            app.histogram_standard_wholevolume_auto.Text = 'histogram(x) - bin auto';

            % Create histogram_standard_wholevolume_1binpervalue
            app.histogram_standard_wholevolume_1binpervalue = uimenu(app.forthewholevolumeMenu);
            app.histogram_standard_wholevolume_1binpervalue.MenuSelectedFcn = createCallbackFcn(app, @histogram_standard_wholevolume_1binpervalueMenuSelected, true);
            app.histogram_standard_wholevolume_1binpervalue.Text = 'histogram(x) - one bin per value';

            % Create histogram_distribution
            app.histogram_distribution = uimenu(app.forthewholevolumeMenu);
            app.histogram_distribution.MenuSelectedFcn = createCallbackFcn(app, @histogram_distributionMenuSelected, true);
            app.histogram_distribution.Text = 'histogram: probability density function';

            % Create foreachsliceMenu
            app.foreachsliceMenu = uimenu(app.PlothistogramMenu);
            app.foreachsliceMenu.MenuSelectedFcn = createCallbackFcn(app, @foreachsliceMenuSelected, true);
            app.foreachsliceMenu.Text = 'for each slice';

            % Create PlotgreylevelMenu
            app.PlotgreylevelMenu = uimenu(app.MenuCalculate);
            app.PlotgreylevelMenu.Text = 'Plot grey level';

            % Create asafunctionofpositionMenu
            app.asafunctionofpositionMenu = uimenu(app.PlotgreylevelMenu);
            app.asafunctionofpositionMenu.MenuSelectedFcn = createCallbackFcn(app, @asafunctionofpositionMenuSelected, true);
            app.asafunctionofpositionMenu.Text = 'as a function of position';

            % Create inaveraged2DmapsMenu
            app.inaveraged2DmapsMenu = uimenu(app.PlotgreylevelMenu);
            app.inaveraged2DmapsMenu.MenuSelectedFcn = createCallbackFcn(app, @TwoDmapMenuSelected, true);
            app.inaveraged2DmapsMenu.Text = 'in averaged 2D maps';

            % Create PlotnoiselevelMenu
            app.PlotnoiselevelMenu = uimenu(app.MenuCalculate);
            app.PlotnoiselevelMenu.Text = 'Plot noise level';

            % Create foreachsliceMenu_2
            app.foreachsliceMenu_2 = uimenu(app.PlotnoiselevelMenu);
            app.foreachsliceMenu_2.MenuSelectedFcn = createCallbackFcn(app, @foreachsliceMenu_2Selected, true);
            app.foreachsliceMenu_2.Text = 'for each slice';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.ROIfilteringandsegmentationmoduleUIFigure);
            app.HelpMenu.Text = 'Help';

            % Create DocumentationMenu
            app.DocumentationMenu = uimenu(app.HelpMenu);
            app.DocumentationMenu.MenuSelectedFcn = createCallbackFcn(app, @DocumentationMenuSelected, true);
            app.DocumentationMenu.Text = 'Documentation';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.ROIfilteringandsegmentationmoduleUIFigure);
            app.TabGroup.TabLocation = 'left';
            app.TabGroup.Position = [1 1 930 776];

            % Create InstructionsTab
            app.InstructionsTab = uitab(app.TabGroup);
            app.InstructionsTab.Title = 'Instructions';

            % Create Logo
            app.Logo = uiimage(app.InstructionsTab);
            app.Logo.ImageClickedFcn = createCallbackFcn(app, @LogoImageClicked, true);
            app.Logo.Position = [17 14 802 100];
            app.Logo.ImageSource = 'logo_NREL.png';

            % Create Instructions1
            app.Instructions1 = uilabel(app.InstructionsTab);
            app.Instructions1.Position = [17 710 802 22];
            app.Instructions1.Text = '- Start by loading a volume data file (.tif or .tiff) with menu selection: Volume/Load volume. Then choose save/display options in the next tab.';

            % Create Instructions_Instructions
            app.Instructions_Instructions = uilabel(app.InstructionsTab);
            app.Instructions_Instructions.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_Instructions.HorizontalAlignment = 'center';
            app.Instructions_Instructions.FontWeight = 'bold';
            app.Instructions_Instructions.Position = [17 742 802 22];
            app.Instructions_Instructions.Text = 'Main instructions';

            % Create Instructions2
            app.Instructions2 = uilabel(app.InstructionsTab);
            app.Instructions2.Position = [17 678 802 22];
            app.Instructions2.Text = '- You can then perform operations in the other tabs.';

            % Create Instructions3
            app.Instructions3 = uilabel(app.InstructionsTab);
            app.Instructions3.Position = [17 646 802 22];
            app.Instructions3.Text = '- Each time you perform an operation (e.g. cropping, downscaling, rotation etc.) the history log is updated to keep track of all the sequential';

            % Create Instructions4
            app.Instructions4 = uilabel(app.InstructionsTab);
            app.Instructions4.Position = [17 625 254 22];
            app.Instructions4.Text = 'operations. The log is displayed in the last tab.';

            % Create Instructions5
            app.Instructions5 = uilabel(app.InstructionsTab);
            app.Instructions5.Position = [17 593 802 22];
            app.Instructions5.Text = '- Once modifcations are done, you can save your changes in a new volume file with Volume/Save volume. The history log is also saved in';

            % Create Instructions6
            app.Instructions6 = uilabel(app.InstructionsTab);
            app.Instructions6.Position = [17 566 802 22];
            app.Instructions6.Text = 'a excel file with the same filename of the new .tif file. You can save the volume at any step, for instance to save intermediate states.';

            % Create Instructions7
            app.Instructions7 = uilabel(app.InstructionsTab);
            app.Instructions7.Position = [17 534 802 22];
            app.Instructions7.Text = '- The algorithms used in this module are explained in the documentation Help/Documentation.';

            % Create VolumeLoadedStatus
            app.VolumeLoadedStatus = uilabel(app.InstructionsTab);
            app.VolumeLoadedStatus.FontWeight = 'bold';
            app.VolumeLoadedStatus.FontColor = [1 0 0];
            app.VolumeLoadedStatus.Position = [17 469 462 22];
            app.VolumeLoadedStatus.Text = 'No volume loaded. Please load a volume: menu selection Volume/Load volume';

            % Create Fileloaded_Path
            app.Fileloaded_Path = uilabel(app.InstructionsTab);
            app.Fileloaded_Path.Visible = 'off';
            app.Fileloaded_Path.Position = [30 412 731 22];
            app.Fileloaded_Path.Text = 'Path: *';

            % Create Fileloaded_Filename
            app.Fileloaded_Filename = uilabel(app.InstructionsTab);
            app.Fileloaded_Filename.Visible = 'off';
            app.Fileloaded_Filename.Position = [30 385 731 22];
            app.Fileloaded_Filename.Text = 'Filename: *.tif';

            % Create Fileloaded_dimension
            app.Fileloaded_dimension = uilabel(app.InstructionsTab);
            app.Fileloaded_dimension.Visible = 'off';
            app.Fileloaded_dimension.Position = [30 358 326 22];
            app.Fileloaded_dimension.Text = 'Image dimension: 3';

            % Create Fileloaded_numberofvoxel
            app.Fileloaded_numberofvoxel = uilabel(app.InstructionsTab);
            app.Fileloaded_numberofvoxel.Visible = 'off';
            app.Fileloaded_numberofvoxel.Position = [30 331 326 22];
            app.Fileloaded_numberofvoxel.Text = 'Number of voxels: * x * x * x';

            % Create Fileloaded_memory_datatype
            app.Fileloaded_memory_datatype = uilabel(app.InstructionsTab);
            app.Fileloaded_memory_datatype.Visible = 'off';
            app.Fileloaded_memory_datatype.Position = [30 304 326 22];
            app.Fileloaded_memory_datatype.Text = 'Format and memory: uint8 and * Mb';

            % Create Fileloaded_min_max
            app.Fileloaded_min_max = uilabel(app.InstructionsTab);
            app.Fileloaded_min_max.Visible = 'off';
            app.Fileloaded_min_max.Position = [30 277 597 22];
            app.Fileloaded_min_max.Text = 'Minimum and maximum values: * and *';

            % Create Fileloaded_text
            app.Fileloaded_text = uilabel(app.InstructionsTab);
            app.Fileloaded_text.FontAngle = 'italic';
            app.Fileloaded_text.Visible = 'off';
            app.Fileloaded_text.Position = [17 439 610 22];
            app.Fileloaded_text.Text = 'The information below corresponds to the loaded volume, not the current state of the volume after modifications';

            % Create Loading_message
            app.Loading_message = uilabel(app.InstructionsTab);
            app.Loading_message.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Loading_message.HorizontalAlignment = 'center';
            app.Loading_message.FontSize = 18;
            app.Loading_message.Visible = 'off';
            app.Loading_message.Position = [307 123 223 38];
            app.Loading_message.Text = 'Loading... please wait';

            % Create SaveOptions
            app.SaveOptions = uitab(app.TabGroup);
            app.SaveOptions.Title = 'Save/display options';

            % Create Save_instructions
            app.Save_instructions = uilabel(app.SaveOptions);
            app.Save_instructions.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Save_instructions.HorizontalAlignment = 'center';
            app.Save_instructions.FontWeight = 'bold';
            app.Save_instructions.Position = [17 742 802 22];
            app.Save_instructions.Text = 'Instructions: select your save folder and format options';

            % Create ClicktoselectsavefolderButton
            app.ClicktoselectsavefolderButton = uibutton(app.SaveOptions, 'push');
            app.ClicktoselectsavefolderButton.ButtonPushedFcn = createCallbackFcn(app, @ClicktoselectsavefolderButtonPushed, true);
            app.ClicktoselectsavefolderButton.Position = [31 678 152 22];
            app.ClicktoselectsavefolderButton.Text = 'Click to select save folder';

            % Create DefaultsavefoldernoneselectedLabel
            app.DefaultsavefoldernoneselectedLabel = uilabel(app.SaveOptions);
            app.DefaultsavefoldernoneselectedLabel.Position = [31 646 781 22];
            app.DefaultsavefoldernoneselectedLabel.Text = 'Default save folder: none selected';

            % Create DefaultfilenamenoneLabel
            app.DefaultfilenamenoneLabel = uilabel(app.SaveOptions);
            app.DefaultfilenamenoneLabel.Position = [31 614 781 22];
            app.DefaultfilenamenoneLabel.Text = 'Default file name: none';

            % Create FolderandfilenameusedifyouselectMenuSavevolumeLabel
            app.FolderandfilenameusedifyouselectMenuSavevolumeLabel = uilabel(app.SaveOptions);
            app.FolderandfilenameusedifyouselectMenuSavevolumeLabel.Position = [17 710 802 22];
            app.FolderandfilenameusedifyouselectMenuSavevolumeLabel.Text = 'Folder and filename used if you select Menu/Save volume/default location.';

            % Create Save_checkbox
            app.Save_checkbox = uicheckbox(app.SaveOptions);
            app.Save_checkbox.Text = 'Save figures (e.g. microstructure views, histogram, etc.).';
            app.Save_checkbox.Position = [17 558 348 22];
            app.Save_checkbox.Value = true;

            % Create FontnameusedforfiguresLabel
            app.FontnameusedforfiguresLabel = uilabel(app.SaveOptions);
            app.FontnameusedforfiguresLabel.HorizontalAlignment = 'right';
            app.FontnameusedforfiguresLabel.Position = [17 529 146 22];
            app.FontnameusedforfiguresLabel.Text = 'Fontname used for figures';

            % Create Save_choiceFont_DropDown
            app.Save_choiceFont_DropDown = uidropdown(app.SaveOptions);
            app.Save_choiceFont_DropDown.Position = [178 529 178 22];

            % Create RegionofInterestViewTab
            app.RegionofInterestViewTab = uitab(app.TabGroup);
            app.RegionofInterestViewTab.Title = 'Region of Interest / View';

            % Create ROI_2Dslice
            app.ROI_2Dslice = uiaxes(app.RegionofInterestViewTab);
            xlabel(app.ROI_2Dslice, 'X')
            ylabel(app.ROI_2Dslice, 'Y')
            app.ROI_2Dslice.Position = [368 202 451 498];

            % Create ViewnormaltoDropDown_2Label
            app.ViewnormaltoDropDown_2Label = uilabel(app.RegionofInterestViewTab);
            app.ViewnormaltoDropDown_2Label.HorizontalAlignment = 'center';
            app.ViewnormaltoDropDown_2Label.Position = [369 156 84 22];
            app.ViewnormaltoDropDown_2Label.Text = 'View normal to ';

            % Create ROI_ViewnormaltoDropDown
            app.ROI_ViewnormaltoDropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_ViewnormaltoDropDown.Items = {'Axe 1', 'Axe 2', 'Axe 3'};
            app.ROI_ViewnormaltoDropDown.ValueChangedFcn = createCallbackFcn(app, @ROI_ViewnormaltoDropDownValueChanged, true);
            app.ROI_ViewnormaltoDropDown.Position = [467 156 84 22];
            app.ROI_ViewnormaltoDropDown.Value = 'Axe 1';

            % Create SliceselectionLabel
            app.SliceselectionLabel = uilabel(app.RegionofInterestViewTab);
            app.SliceselectionLabel.HorizontalAlignment = 'right';
            app.SliceselectionLabel.Position = [369 114 86 22];
            app.SliceselectionLabel.Text = 'Slice  selection';

            % Create ROI_SliceSlider
            app.ROI_SliceSlider = uislider(app.RegionofInterestViewTab);
            app.ROI_SliceSlider.ValueChangedFcn = createCallbackFcn(app, @ROI_SliceSliderValueChanged, true);
            app.ROI_SliceSlider.ValueChangingFcn = createCallbackFcn(app, @ROI_SliceSliderValueChanging, true);
            app.ROI_SliceSlider.Position = [467 123 342 3];

            % Create ColormapLabel
            app.ColormapLabel = uilabel(app.RegionofInterestViewTab);
            app.ColormapLabel.HorizontalAlignment = 'center';
            app.ColormapLabel.Position = [572 58 79 22];
            app.ColormapLabel.Text = 'Colormap';

            % Create ROI_colormap
            app.ROI_colormap = uidropdown(app.RegionofInterestViewTab);
            app.ROI_colormap.Items = {'gray', 'bone', 'copper', 'jet', 'parula', 'cool'};
            app.ROI_colormap.ValueChangedFcn = createCallbackFcn(app, @ROI_colormapValueChanged, true);
            app.ROI_colormap.Position = [650 58 79 22];
            app.ROI_colormap.Value = 'copper';

            % Create ROI_instructions_tab
            app.ROI_instructions_tab = uilabel(app.RegionofInterestViewTab);
            app.ROI_instructions_tab.BackgroundColor = [0.4706 0.6706 0.1882];
            app.ROI_instructions_tab.HorizontalAlignment = 'center';
            app.ROI_instructions_tab.FontWeight = 'bold';
            app.ROI_instructions_tab.Position = [17 720 802 44];
            app.ROI_instructions_tab.Text = {'Instructions: You can perform operations (cropping, rotation, etc.) in any order. To visualize changes, click on the ''Do'' button.'; 'To cancel, click on the ''Undo'' button. To keep changes and update the history log, click on the ''Save'' button.'};

            % Create ROI_checkbox_updateslicer
            app.ROI_checkbox_updateslicer = uicheckbox(app.RegionofInterestViewTab);
            app.ROI_checkbox_updateslicer.Text = 'Update slice while moving slicer';
            app.ROI_checkbox_updateslicer.Position = [627 156 193 22];

            % Create ROI_View3D
            app.ROI_View3D = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_View3D.ButtonPushedFcn = createCallbackFcn(app, @ROI_View3DButtonPushed, true);
            app.ROI_View3D.Position = [744 58 76 22];
            app.ROI_View3D.Text = '3D View';

            % Create ROI_SliceSpinner
            app.ROI_SliceSpinner = uispinner(app.RegionofInterestViewTab);
            app.ROI_SliceSpinner.ValueChangingFcn = createCallbackFcn(app, @ROI_SliceSpinnerValueChanging, true);
            app.ROI_SliceSpinner.Limits = [0 100];
            app.ROI_SliceSpinner.ValueChangedFcn = createCallbackFcn(app, @ROI_SliceSpinnerValueChanged, true);
            app.ROI_SliceSpinner.Position = [467 58 94 22];

            % Create ROI_table_crop
            app.ROI_table_crop = uitable(app.RegionofInterestViewTab);
            app.ROI_table_crop.ColumnName = {'Axe'; 'Start'; 'End'; 'Length'};
            app.ROI_table_crop.RowName = {};
            app.ROI_table_crop.ColumnEditable = [false true true false];
            app.ROI_table_crop.CellEditCallback = createCallbackFcn(app, @ROI_table_cropCellEdit, true);
            app.ROI_table_crop.CellSelectionCallback = createCallbackFcn(app, @ROI_table_cropCellEdit, true);
            app.ROI_table_crop.Position = [18 99 304 110];

            % Create ROI_crop_text
            app.ROI_crop_text = uilabel(app.RegionofInterestViewTab);
            app.ROI_crop_text.FontSize = 14;
            app.ROI_crop_text.FontWeight = 'bold';
            app.ROI_crop_text.Position = [18 214 91 22];
            app.ROI_crop_text.Text = 'Crop volume';

            % Create ROI_crop_DoButton
            app.ROI_crop_DoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_crop_DoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_crop_DoButtonPushed, true);
            app.ROI_crop_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.ROI_crop_DoButton.FontSize = 14;
            app.ROI_crop_DoButton.FontWeight = 'bold';
            app.ROI_crop_DoButton.FontColor = [1 1 1];
            app.ROI_crop_DoButton.Position = [100 15 71 25];
            app.ROI_crop_DoButton.Text = 'Do';

            % Create ROI_crop_UndoButton
            app.ROI_crop_UndoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_crop_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_crop_UndoButtonPushed, true);
            app.ROI_crop_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.ROI_crop_UndoButton.FontSize = 14;
            app.ROI_crop_UndoButton.FontWeight = 'bold';
            app.ROI_crop_UndoButton.FontColor = [1 1 1];
            app.ROI_crop_UndoButton.Enable = 'off';
            app.ROI_crop_UndoButton.Position = [181 15 71 25];
            app.ROI_crop_UndoButton.Text = 'Undo';

            % Create ROI_crop_SaveButton
            app.ROI_crop_SaveButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_crop_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_crop_SaveButtonPushed, true);
            app.ROI_crop_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.ROI_crop_SaveButton.FontSize = 14;
            app.ROI_crop_SaveButton.FontWeight = 'bold';
            app.ROI_crop_SaveButton.FontColor = [1 1 1];
            app.ROI_crop_SaveButton.Enable = 'off';
            app.ROI_crop_SaveButton.Position = [262 15 71 25];
            app.ROI_crop_SaveButton.Text = 'Save';

            % Create VoxelsizeoptionalLabel
            app.VoxelsizeoptionalLabel = uilabel(app.RegionofInterestViewTab);
            app.VoxelsizeoptionalLabel.HorizontalAlignment = 'right';
            app.VoxelsizeoptionalLabel.Position = [156 214 113 22];
            app.VoxelsizeoptionalLabel.Text = 'Voxel size (optional)';

            % Create ROI_voxelsize
            app.ROI_voxelsize = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.ROI_voxelsize.Limits = [1e-09 Inf];
            app.ROI_voxelsize.ValueChangedFcn = createCallbackFcn(app, @ROI_voxelsizeValueChanged, true);
            app.ROI_voxelsize.Position = [275 214 47 22];
            app.ROI_voxelsize.Value = 1;

            % Create ROI_rotation_text
            app.ROI_rotation_text = uilabel(app.RegionofInterestViewTab);
            app.ROI_rotation_text.FontSize = 14;
            app.ROI_rotation_text.FontWeight = 'bold';
            app.ROI_rotation_text.Position = [17 684 62 22];
            app.ROI_rotation_text.Text = 'Rotation';

            % Create ROI_angle_spinner
            app.ROI_angle_spinner = uispinner(app.RegionofInterestViewTab);
            app.ROI_angle_spinner.Step = 0.5;
            app.ROI_angle_spinner.Limits = [-180 180];
            app.ROI_angle_spinner.ValueChangedFcn = createCallbackFcn(app, @ROI_angle_spinnerValueChanged, true);
            app.ROI_angle_spinner.Position = [250 657 73 22];

            % Create AngleDropDownLabel
            app.AngleDropDownLabel = uilabel(app.RegionofInterestViewTab);
            app.AngleDropDownLabel.HorizontalAlignment = 'right';
            app.AngleDropDownLabel.Position = [17 657 36 22];
            app.AngleDropDownLabel.Text = 'Angle';

            % Create ROI_angle_choice_DropDown
            app.ROI_angle_choice_DropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_angle_choice_DropDown.Items = {'Normal to axe 1', 'Normal to axe 2', 'Normal to axe 3'};
            app.ROI_angle_choice_DropDown.Position = [68 657 174 22];
            app.ROI_angle_choice_DropDown.Value = 'Normal to axe 1';

            % Create ROI_rotation_DoButton
            app.ROI_rotation_DoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_rotation_DoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_rotation_DoButtonPushed, true);
            app.ROI_rotation_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.ROI_rotation_DoButton.FontSize = 14;
            app.ROI_rotation_DoButton.FontWeight = 'bold';
            app.ROI_rotation_DoButton.FontColor = [1 1 1];
            app.ROI_rotation_DoButton.Position = [90 627 71 25];
            app.ROI_rotation_DoButton.Text = 'Do';

            % Create ROI_rotation_UndoButton
            app.ROI_rotation_UndoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_rotation_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_rotation_UndoButtonPushed, true);
            app.ROI_rotation_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.ROI_rotation_UndoButton.FontSize = 14;
            app.ROI_rotation_UndoButton.FontWeight = 'bold';
            app.ROI_rotation_UndoButton.FontColor = [1 1 1];
            app.ROI_rotation_UndoButton.Enable = 'off';
            app.ROI_rotation_UndoButton.Position = [171 627 71 25];
            app.ROI_rotation_UndoButton.Text = 'Undo';

            % Create ROI_rotation_SaveButton
            app.ROI_rotation_SaveButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_rotation_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_rotation_SaveButtonPushed, true);
            app.ROI_rotation_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.ROI_rotation_SaveButton.FontSize = 14;
            app.ROI_rotation_SaveButton.FontWeight = 'bold';
            app.ROI_rotation_SaveButton.FontColor = [1 1 1];
            app.ROI_rotation_SaveButton.Enable = 'off';
            app.ROI_rotation_SaveButton.Position = [252 627 71 25];
            app.ROI_rotation_SaveButton.Text = 'Save';

            % Create ROI_flipswapDropDown
            app.ROI_flipswapDropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_flipswapDropDown.Items = {'Flip axis 1', 'Flip axis 2', 'Flip axis 3', 'Swap axis 1 with axis 2', 'Swap axis 1 with axis 3', 'Swap axis 2 with axis 3'};
            app.ROI_flipswapDropDown.Position = [17 553 306 22];
            app.ROI_flipswapDropDown.Value = 'Flip axis 1';

            % Create FliporswapaxisLabel
            app.FliporswapaxisLabel = uilabel(app.RegionofInterestViewTab);
            app.FliporswapaxisLabel.FontSize = 14;
            app.FliporswapaxisLabel.FontWeight = 'bold';
            app.FliporswapaxisLabel.Position = [19 580 118 22];
            app.FliporswapaxisLabel.Text = 'Flip or swap axis';

            % Create ROI_flipswap_DoButton
            app.ROI_flipswap_DoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_flipswap_DoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_flipswap_DoButtonPushed, true);
            app.ROI_flipswap_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.ROI_flipswap_DoButton.FontSize = 14;
            app.ROI_flipswap_DoButton.FontWeight = 'bold';
            app.ROI_flipswap_DoButton.FontColor = [1 1 1];
            app.ROI_flipswap_DoButton.Position = [90 523 71 25];
            app.ROI_flipswap_DoButton.Text = 'Do';

            % Create ROI_flipswap_UndoButton
            app.ROI_flipswap_UndoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_flipswap_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_flipswap_UndoButtonPushed, true);
            app.ROI_flipswap_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.ROI_flipswap_UndoButton.FontSize = 14;
            app.ROI_flipswap_UndoButton.FontWeight = 'bold';
            app.ROI_flipswap_UndoButton.FontColor = [1 1 1];
            app.ROI_flipswap_UndoButton.Enable = 'off';
            app.ROI_flipswap_UndoButton.Position = [171 523 71 25];
            app.ROI_flipswap_UndoButton.Text = 'Undo';

            % Create ROI_flipswap_SaveButton
            app.ROI_flipswap_SaveButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_flipswap_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_flipswap_SaveButtonPushed, true);
            app.ROI_flipswap_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.ROI_flipswap_SaveButton.FontSize = 14;
            app.ROI_flipswap_SaveButton.FontWeight = 'bold';
            app.ROI_flipswap_SaveButton.FontColor = [1 1 1];
            app.ROI_flipswap_SaveButton.Enable = 'off';
            app.ROI_flipswap_SaveButton.Position = [252 523 71 25];
            app.ROI_flipswap_SaveButton.Text = 'Save';

            % Create ROI_ModifysomeslicesalongLabel
            app.ROI_ModifysomeslicesalongLabel = uilabel(app.RegionofInterestViewTab);
            app.ROI_ModifysomeslicesalongLabel.FontSize = 14;
            app.ROI_ModifysomeslicesalongLabel.FontWeight = 'bold';
            app.ROI_ModifysomeslicesalongLabel.Tooltip = {'From slice A (not included) to slice B (not included). Interpolate with ''Label'': no intermediate values between slice A and B. Interpolate with ''grey level'': intermediate values between slice A and B'};
            app.ROI_ModifysomeslicesalongLabel.Position = [17 476 175 22];
            app.ROI_ModifysomeslicesalongLabel.Text = 'Modify some slices along';

            % Create ROI_modify_axeAxeDropDown
            app.ROI_modify_axeAxeDropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_modify_axeAxeDropDown.Items = {'Axe 1', 'Axe 2', 'Axe 3'};
            app.ROI_modify_axeAxeDropDown.ValueChangedFcn = createCallbackFcn(app, @ROI_modify_axeAxeDropDownValueChanged, true);
            app.ROI_modify_axeAxeDropDown.Tooltip = {''};
            app.ROI_modify_axeAxeDropDown.Position = [241 476 82 22];
            app.ROI_modify_axeAxeDropDown.Value = 'Axe 1';

            % Create FromsliceAEditFieldLabel
            app.FromsliceAEditFieldLabel = uilabel(app.RegionofInterestViewTab);
            app.FromsliceAEditFieldLabel.HorizontalAlignment = 'right';
            app.FromsliceAEditFieldLabel.Tooltip = {''};
            app.FromsliceAEditFieldLabel.Position = [17 449 72 22];
            app.FromsliceAEditFieldLabel.Text = 'From slice A';

            % Create ROI_modify_FromsliceA
            app.ROI_modify_FromsliceA = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.ROI_modify_FromsliceA.Limits = [0 Inf];
            app.ROI_modify_FromsliceA.ValueChangedFcn = createCallbackFcn(app, @ROI_modify_FromsliceAValueChanged, true);
            app.ROI_modify_FromsliceA.Tooltip = {''};
            app.ROI_modify_FromsliceA.Position = [100 449 54 22];

            % Create TosliceBEditFieldLabel
            app.TosliceBEditFieldLabel = uilabel(app.RegionofInterestViewTab);
            app.TosliceBEditFieldLabel.HorizontalAlignment = 'right';
            app.TosliceBEditFieldLabel.Tooltip = {''};
            app.TosliceBEditFieldLabel.Position = [196 449 57 22];
            app.TosliceBEditFieldLabel.Text = 'To slice B';

            % Create ROI_modify_TosliceB
            app.ROI_modify_TosliceB = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.ROI_modify_TosliceB.Limits = [0 Inf];
            app.ROI_modify_TosliceB.ValueChangedFcn = createCallbackFcn(app, @ROI_modify_FromsliceAValueChanged, true);
            app.ROI_modify_TosliceB.Tooltip = {''};
            app.ROI_modify_TosliceB.Position = [269 449 54 22];

            % Create ActionforSlicesfromAtoBDropDownLabel
            app.ActionforSlicesfromAtoBDropDownLabel = uilabel(app.RegionofInterestViewTab);
            app.ActionforSlicesfromAtoBDropDownLabel.HorizontalAlignment = 'right';
            app.ActionforSlicesfromAtoBDropDownLabel.Tooltip = {''};
            app.ActionforSlicesfromAtoBDropDownLabel.Position = [17 422 155 22];
            app.ActionforSlicesfromAtoBDropDownLabel.Text = 'Action for Slices from A to B';

            % Create ROI_modify_ActionDropDown
            app.ROI_modify_ActionDropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_modify_ActionDropDown.Items = {'Remove slices', 'Interpolate slices (grey level)', 'Interpolate slices (label)'};
            app.ROI_modify_ActionDropDown.Tooltip = {''};
            app.ROI_modify_ActionDropDown.Position = [187 422 136 22];
            app.ROI_modify_ActionDropDown.Value = 'Remove slices';

            % Create ROI_modify_DoButton
            app.ROI_modify_DoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_modify_DoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_modify_DoButtonPushed, true);
            app.ROI_modify_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.ROI_modify_DoButton.FontSize = 14;
            app.ROI_modify_DoButton.FontWeight = 'bold';
            app.ROI_modify_DoButton.FontColor = [1 1 1];
            app.ROI_modify_DoButton.Tooltip = {''};
            app.ROI_modify_DoButton.Position = [90 392 71 25];
            app.ROI_modify_DoButton.Text = 'Do';

            % Create ROI_modify_UndoButton
            app.ROI_modify_UndoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_modify_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_modify_UndoButtonPushed, true);
            app.ROI_modify_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.ROI_modify_UndoButton.FontSize = 14;
            app.ROI_modify_UndoButton.FontWeight = 'bold';
            app.ROI_modify_UndoButton.FontColor = [1 1 1];
            app.ROI_modify_UndoButton.Enable = 'off';
            app.ROI_modify_UndoButton.Tooltip = {''};
            app.ROI_modify_UndoButton.Position = [171 392 71 25];
            app.ROI_modify_UndoButton.Text = 'Undo';

            % Create ROI_modify_SaveButton
            app.ROI_modify_SaveButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_modify_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_modify_SaveButtonPushed, true);
            app.ROI_modify_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.ROI_modify_SaveButton.FontSize = 14;
            app.ROI_modify_SaveButton.FontWeight = 'bold';
            app.ROI_modify_SaveButton.FontColor = [1 1 1];
            app.ROI_modify_SaveButton.Enable = 'off';
            app.ROI_modify_SaveButton.Tooltip = {''};
            app.ROI_modify_SaveButton.Position = [252 392 71 25];
            app.ROI_modify_SaveButton.Text = 'Save';

            % Create CropbackgroundButton
            app.CropbackgroundButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.CropbackgroundButton.ButtonPushedFcn = createCallbackFcn(app, @CropbackgroundButtonPushed, true);
            app.CropbackgroundButton.Tooltip = {'1) set the value of the background. 2) choose for what case the crop value will be estimated (''Auto-crop for'' dropdown selection). 3) Click to estimate crop values.'};
            app.CropbackgroundButton.Position = [214 72 108 22];
            app.CropbackgroundButton.Text = 'Crop background';

            % Create BackgroundvalueEditFieldLabel
            app.BackgroundvalueEditFieldLabel = uilabel(app.RegionofInterestViewTab);
            app.BackgroundvalueEditFieldLabel.HorizontalAlignment = 'right';
            app.BackgroundvalueEditFieldLabel.Position = [18 72 102 22];
            app.BackgroundvalueEditFieldLabel.Text = 'Background value';

            % Create BackgroundvalueEditField
            app.BackgroundvalueEditField = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.BackgroundvalueEditField.Position = [129 72 40 22];

            % Create MedianfilterrangeEditFieldLabel
            app.MedianfilterrangeEditFieldLabel = uilabel(app.RegionofInterestViewTab);
            app.MedianfilterrangeEditFieldLabel.HorizontalAlignment = 'right';
            app.MedianfilterrangeEditFieldLabel.Position = [17 318 105 22];
            app.MedianfilterrangeEditFieldLabel.Text = 'Median filter range';

            % Create MedianfilterrangeEditField
            app.MedianfilterrangeEditField = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.MedianfilterrangeEditField.Limits = [0 Inf];
            app.MedianfilterrangeEditField.Position = [127 318 40 22];
            app.MedianfilterrangeEditField.Value = 7;

            % Create StdthresholdEditFieldLabel
            app.StdthresholdEditFieldLabel = uilabel(app.RegionofInterestViewTab);
            app.StdthresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.StdthresholdEditFieldLabel.Position = [207 318 76 22];
            app.StdthresholdEditFieldLabel.Text = 'Std threshold';

            % Create StdthresholdEditField
            app.StdthresholdEditField = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.StdthresholdEditField.Limits = [0 Inf];
            app.StdthresholdEditField.Position = [288 318 40 22];
            app.StdthresholdEditField.Value = 0.01;

            % Create VieworthogonalslicesButton
            app.VieworthogonalslicesButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.VieworthogonalslicesButton.ButtonPushedFcn = createCallbackFcn(app, @VieworthogonalslicesButtonPushed, true);
            app.VieworthogonalslicesButton.Tooltip = {'Warning: 3D slice requires a lor of RAM and CPU. You may have to crop your volume first.'};
            app.VieworthogonalslicesButton.Position = [582 21 238 22];
            app.VieworthogonalslicesButton.Text = 'Plot and save 3D view (orthogonal slices)';

            % Create Savecurrent2DviewButton
            app.Savecurrent2DviewButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.Savecurrent2DviewButton.ButtonPushedFcn = createCallbackFcn(app, @Savecurrent2DviewButtonPushed, true);
            app.Savecurrent2DviewButton.Position = [431 21 130 22];
            app.Savecurrent2DviewButton.Text = 'Save current 2D view';

            % Create ROI_ModifysomeslicesalongLabel_2
            app.ROI_ModifysomeslicesalongLabel_2 = uilabel(app.RegionofInterestViewTab);
            app.ROI_ModifysomeslicesalongLabel_2.FontSize = 14;
            app.ROI_ModifysomeslicesalongLabel_2.FontWeight = 'bold';
            app.ROI_ModifysomeslicesalongLabel_2.Tooltip = {'Median filter is applied to the image standard deviation and then thresholded.'};
            app.ROI_ModifysomeslicesalongLabel_2.Position = [17 345 153 22];
            app.ROI_ModifysomeslicesalongLabel_2.Text = 'Background detection';

            % Create Label
            app.Label = uilabel(app.RegionofInterestViewTab);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [188 345 25 22];
            app.Label.Text = '';

            % Create ROI_backgroundetection_DropDown
            app.ROI_backgroundetection_DropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_backgroundetection_DropDown.Items = {'Search only first slice of axe 3, then apply to whole volume', 'Search slice per slice along axe 3'};
            app.ROI_backgroundetection_DropDown.Position = [192 345 136 22];
            app.ROI_backgroundetection_DropDown.Value = 'Search only first slice of axe 3, then apply to whole volume';

            % Create ROI_background_DoButton
            app.ROI_background_DoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_background_DoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_background_DoButtonPushed, true);
            app.ROI_background_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.ROI_background_DoButton.FontSize = 14;
            app.ROI_background_DoButton.FontWeight = 'bold';
            app.ROI_background_DoButton.FontColor = [1 1 1];
            app.ROI_background_DoButton.Tooltip = {''};
            app.ROI_background_DoButton.Position = [100 261 71 25];
            app.ROI_background_DoButton.Text = 'Do';

            % Create ROI_background_UndoButton
            app.ROI_background_UndoButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_background_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_background_UndoButtonPushed, true);
            app.ROI_background_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.ROI_background_UndoButton.FontSize = 14;
            app.ROI_background_UndoButton.FontWeight = 'bold';
            app.ROI_background_UndoButton.FontColor = [1 1 1];
            app.ROI_background_UndoButton.Enable = 'off';
            app.ROI_background_UndoButton.Tooltip = {''};
            app.ROI_background_UndoButton.Position = [181 261 71 25];
            app.ROI_background_UndoButton.Text = 'Undo';

            % Create ROI_background_SaveButton
            app.ROI_background_SaveButton = uibutton(app.RegionofInterestViewTab, 'push');
            app.ROI_background_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @ROI_background_SaveButtonPushed, true);
            app.ROI_background_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.ROI_background_SaveButton.FontSize = 14;
            app.ROI_background_SaveButton.FontWeight = 'bold';
            app.ROI_background_SaveButton.FontColor = [1 1 1];
            app.ROI_background_SaveButton.Enable = 'off';
            app.ROI_background_SaveButton.Tooltip = {''};
            app.ROI_background_SaveButton.Position = [262 261 71 25];
            app.ROI_background_SaveButton.Text = 'Save';

            % Create SetidentifiedbackgroundvaluetoEditFieldLabel
            app.SetidentifiedbackgroundvaluetoEditFieldLabel = uilabel(app.RegionofInterestViewTab);
            app.SetidentifiedbackgroundvaluetoEditFieldLabel.HorizontalAlignment = 'right';
            app.SetidentifiedbackgroundvaluetoEditFieldLabel.Position = [32 291 186 22];
            app.SetidentifiedbackgroundvaluetoEditFieldLabel.Text = 'Set identified background value to';

            % Create SetidentifiedbackgroundvaluetoEditField
            app.SetidentifiedbackgroundvaluetoEditField = uieditfield(app.RegionofInterestViewTab, 'numeric');
            app.SetidentifiedbackgroundvaluetoEditField.Position = [223 291 40 22];

            % Create AutocropforLabel
            app.AutocropforLabel = uilabel(app.RegionofInterestViewTab);
            app.AutocropforLabel.HorizontalAlignment = 'right';
            app.AutocropforLabel.Position = [18 45 75 22];
            app.AutocropforLabel.Text = 'Auto-crop for';

            % Create ROI_background_OptionsDropDown
            app.ROI_background_OptionsDropDown = uidropdown(app.RegionofInterestViewTab);
            app.ROI_background_OptionsDropDown.Items = {'Cylindrical FOV', 'After rotation'};
            app.ROI_background_OptionsDropDown.Tooltip = {'Cylindrical FOV: the largest rectangle that does not contain background will be found, based on analysis of first slice of Axe 3.'; 'After rotation: rotation creates background. Axe 1 and 2 bounds will be detected.'};
            app.ROI_background_OptionsDropDown.Position = [108 45 100 22];
            app.ROI_background_OptionsDropDown.Value = 'Cylindrical FOV';

            % Create FormatconversionTab
            app.FormatconversionTab = uitab(app.TabGroup);
            app.FormatconversionTab.Title = 'Format conversion';

            % Create TovisualizeimpactonhistogramCalculateHistogramLabel
            app.TovisualizeimpactonhistogramCalculateHistogramLabel = uilabel(app.FormatconversionTab);
            app.TovisualizeimpactonhistogramCalculateHistogramLabel.Position = [17 436 377 22];
            app.TovisualizeimpactonhistogramCalculateHistogramLabel.Text = 'To visualize impact on histogram: Calculate/Histogram';

            % Create ConverttoLabel
            app.ConverttoLabel = uilabel(app.FormatconversionTab);
            app.ConverttoLabel.HorizontalAlignment = 'right';
            app.ConverttoLabel.Position = [17 677 61 22];
            app.ConverttoLabel.Text = 'Convert to';

            % Create Convert_Choice_DropDown
            app.Convert_Choice_DropDown = uidropdown(app.FormatconversionTab);
            app.Convert_Choice_DropDown.Items = {'8-bit unsigned integer arrays', '8-bit unsigned integer arrays (rescale)', '16-bit unsigned integer arrays', '16-bit unsigned integer arrays (rescale)', '32-bit unsigned integer arrays', '64-bit unsigned integer arrays', 'double precision array (rescale)'};
            app.Convert_Choice_DropDown.Position = [93 677 207 22];
            app.Convert_Choice_DropDown.Value = '8-bit unsigned integer arrays (rescale)';

            % Create Convert_instructions
            app.Convert_instructions = uilabel(app.FormatconversionTab);
            app.Convert_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Convert_instructions.HorizontalAlignment = 'center';
            app.Convert_instructions.FontWeight = 'bold';
            app.Convert_instructions.Position = [17 720 802 44];
            app.Convert_instructions.Text = {'Instructions: Select new data type. To apply, click on the ''Do'' button.'; 'To cancel, click on the ''Undo'' button. To keep changes and update the history log, click on the ''Save'' button.'};

            % Create Convert_DoButton
            app.Convert_DoButton = uibutton(app.FormatconversionTab, 'push');
            app.Convert_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Convert_DoButtonPushed, true);
            app.Convert_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Convert_DoButton.FontSize = 14;
            app.Convert_DoButton.FontWeight = 'bold';
            app.Convert_DoButton.FontColor = [1 1 1];
            app.Convert_DoButton.Position = [67 647 71 25];
            app.Convert_DoButton.Text = 'Do';

            % Create Convert_UndoButton
            app.Convert_UndoButton = uibutton(app.FormatconversionTab, 'push');
            app.Convert_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Convert_UndoButtonPushed, true);
            app.Convert_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Convert_UndoButton.FontSize = 14;
            app.Convert_UndoButton.FontWeight = 'bold';
            app.Convert_UndoButton.FontColor = [1 1 1];
            app.Convert_UndoButton.Enable = 'off';
            app.Convert_UndoButton.Position = [148 647 71 25];
            app.Convert_UndoButton.Text = 'Undo';

            % Create Convert_SaveButton
            app.Convert_SaveButton = uibutton(app.FormatconversionTab, 'push');
            app.Convert_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Convert_SaveButtonPushed, true);
            app.Convert_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Convert_SaveButton.FontSize = 14;
            app.Convert_SaveButton.FontWeight = 'bold';
            app.Convert_SaveButton.FontColor = [1 1 1];
            app.Convert_SaveButton.Enable = 'off';
            app.Convert_SaveButton.Position = [229 647 71 25];
            app.Convert_SaveButton.Text = 'Save';

            % Create Convert_table
            app.Convert_table = uitable(app.FormatconversionTab);
            app.Convert_table.ColumnName = {''; 'Current'; 'New'};
            app.Convert_table.RowName = {};
            app.Convert_table.Position = [17 470 377 137];

            % Create Convert_Utilization_checkbox
            app.Convert_Utilization_checkbox = uicheckbox(app.FormatconversionTab);
            app.Convert_Utilization_checkbox.Tooltip = {'Number of unique values / number of possible values. For instance, with uint8, the number of possible values is 2^8=256.'};
            app.Convert_Utilization_checkbox.Text = 'Calculate utilization ratio of data type';
            app.Convert_Utilization_checkbox.Position = [17 407 219 22];

            % Create Convert_UpdatetableButton
            app.Convert_UpdatetableButton = uibutton(app.FormatconversionTab, 'push');
            app.Convert_UpdatetableButton.ButtonPushedFcn = createCallbackFcn(app, @Convert_UpdatetableButtonPushed, true);
            app.Convert_UpdatetableButton.Position = [294 407 100 22];
            app.Convert_UpdatetableButton.Text = 'Update table';

            % Create UpdownscalingTab
            app.UpdownscalingTab = uitab(app.TabGroup);
            app.UpdownscalingTab.Title = 'Up/down scaling';

            % Create Scaling_UIAxes_before
            app.Scaling_UIAxes_before = uiaxes(app.UpdownscalingTab);
            title(app.Scaling_UIAxes_before, 'Before scaling. View normal to Axe 3')
            xlabel(app.Scaling_UIAxes_before, '2nd Axis')
            ylabel(app.Scaling_UIAxes_before, '1st Axis')
            app.Scaling_UIAxes_before.Position = [17 103 400 385];

            % Create Scaling_UIAxes_after
            app.Scaling_UIAxes_after = uiaxes(app.UpdownscalingTab);
            title(app.Scaling_UIAxes_after, 'After scaling. View normal to Axe 3')
            xlabel(app.Scaling_UIAxes_after, '2nd Axis')
            ylabel(app.Scaling_UIAxes_after, '1st Axis')
            app.Scaling_UIAxes_after.Position = [419 103 400 385];

            % Create Scale_text_2
            app.Scale_text_2 = uilabel(app.UpdownscalingTab);
            app.Scale_text_2.Position = [19 584 397 27];
            app.Scale_text_2.Text = 'size times scaling factor. A scaling factor value above 1 is for downscaling,';

            % Create Scale_text_3
            app.Scale_text_3 = uilabel(app.UpdownscalingTab);
            app.Scale_text_3.Position = [19 558 397 27];
            app.Scale_text_3.Text = 'while a scaling factor value below 1 is upscaling.';

            % Create ScalingfactorEditFieldLabel
            app.ScalingfactorEditFieldLabel = uilabel(app.UpdownscalingTab);
            app.ScalingfactorEditFieldLabel.HorizontalAlignment = 'right';
            app.ScalingfactorEditFieldLabel.Position = [19 531 78 22];
            app.ScalingfactorEditFieldLabel.Text = 'Scaling factor';

            % Create Scale_scalingfactorEditField
            app.Scale_scalingfactorEditField = uieditfield(app.UpdownscalingTab, 'numeric');
            app.Scale_scalingfactorEditField.Limits = [0 Inf];
            app.Scale_scalingfactorEditField.ValueChangedFcn = createCallbackFcn(app, @Scale_scalingfactorEditFieldValueChanged, true);
            app.Scale_scalingfactorEditField.Position = [112 531 50 22];
            app.Scale_scalingfactorEditField.Value = 1;

            % Create Scale_text_1
            app.Scale_text_1 = uilabel(app.UpdownscalingTab);
            app.Scale_text_1.Position = [19 670 411 27];
            app.Scale_text_1.Text = 'Please enter the initial voxel size (optional, you can let the unit value).';

            % Create InitialvoxelsizeEditFieldLabel
            app.InitialvoxelsizeEditFieldLabel = uilabel(app.UpdownscalingTab);
            app.InitialvoxelsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.InitialvoxelsizeEditFieldLabel.Position = [20 647 89 22];
            app.InitialvoxelsizeEditFieldLabel.Text = 'Initial voxel size';

            % Create Scale_InitialvoxelsizeEditField
            app.Scale_InitialvoxelsizeEditField = uieditfield(app.UpdownscalingTab, 'numeric');
            app.Scale_InitialvoxelsizeEditField.Limits = [0 Inf];
            app.Scale_InitialvoxelsizeEditField.ValueChangedFcn = createCallbackFcn(app, @Scale_InitialvoxelsizeEditFieldValueChanged, true);
            app.Scale_InitialvoxelsizeEditField.Position = [124 647 50 22];
            app.Scale_InitialvoxelsizeEditField.Value = 1;

            % Create VoxelunitDropDownLabel
            app.VoxelunitDropDownLabel = uilabel(app.UpdownscalingTab);
            app.VoxelunitDropDownLabel.HorizontalAlignment = 'right';
            app.VoxelunitDropDownLabel.Position = [307 647 57 22];
            app.VoxelunitDropDownLabel.Text = 'Voxel unit';

            % Create Scale_VoxelunitDropDown
            app.Scale_VoxelunitDropDown = uidropdown(app.UpdownscalingTab);
            app.Scale_VoxelunitDropDown.Items = {'nanometers', 'micrometers', 'millimeters', 'meters'};
            app.Scale_VoxelunitDropDown.ValueChangedFcn = createCallbackFcn(app, @Scale_VoxelunitDropDownValueChanged, true);
            app.Scale_VoxelunitDropDown.Position = [187 647 120 22];
            app.Scale_VoxelunitDropDown.Value = 'nanometers';

            % Create Scale_table
            app.Scale_table = uitable(app.UpdownscalingTab);
            app.Scale_table.ColumnName = {'Axe'; 'Initial number of voxel'; 'New number of voxel'};
            app.Scale_table.RowName = {''};
            app.Scale_table.ColumnEditable = [false false false];
            app.Scale_table.Position = [442 595 377 102];

            % Create NewvoxelsizeEditFieldLabel
            app.NewvoxelsizeEditFieldLabel = uilabel(app.UpdownscalingTab);
            app.NewvoxelsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.NewvoxelsizeEditFieldLabel.Position = [186 531 86 22];
            app.NewvoxelsizeEditFieldLabel.Text = 'New voxel size';

            % Create Scale_NewvoxelsizeEditField
            app.Scale_NewvoxelsizeEditField = uieditfield(app.UpdownscalingTab, 'numeric');
            app.Scale_NewvoxelsizeEditField.Limits = [0 Inf];
            app.Scale_NewvoxelsizeEditField.Editable = 'off';
            app.Scale_NewvoxelsizeEditField.Position = [286 531 50 22];
            app.Scale_NewvoxelsizeEditField.Value = 1;

            % Create Scale_newvoxelsizeunitLabel
            app.Scale_newvoxelsizeunitLabel = uilabel(app.UpdownscalingTab);
            app.Scale_newvoxelsizeunitLabel.Position = [339 531 77 22];
            app.Scale_newvoxelsizeunitLabel.Text = 'nanometers';

            % Create Scale_text_4
            app.Scale_text_4 = uilabel(app.UpdownscalingTab);
            app.Scale_text_4.Position = [19 610 397 27];
            app.Scale_text_4.Text = 'Please enter below the scaling factor. New voxel size equals initial voxel';

            % Create Scale_DatatypeButtonGroup
            app.Scale_DatatypeButtonGroup = uibuttongroup(app.UpdownscalingTab);
            app.Scale_DatatypeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @Scale_DatatypeButtonGroupSelectionChanged, true);
            app.Scale_DatatypeButtonGroup.Tooltip = {'Grey level: intermediate values can be created. Label: intermediate values cannot be created.'};
            app.Scale_DatatypeButtonGroup.TitlePosition = 'centertop';
            app.Scale_DatatypeButtonGroup.Title = 'Data type';
            app.Scale_DatatypeButtonGroup.Position = [442 531 194 53];

            % Create GreylevelButton
            app.GreylevelButton = uitogglebutton(app.Scale_DatatypeButtonGroup);
            app.GreylevelButton.Text = 'Grey level';
            app.GreylevelButton.Position = [13 6 81 22];
            app.GreylevelButton.Value = true;

            % Create LabelButton
            app.LabelButton = uitogglebutton(app.Scale_DatatypeButtonGroup);
            app.LabelButton.Text = 'Label';
            app.LabelButton.Position = [104 6 81 22];

            % Create SliceselectionSliderLabel
            app.SliceselectionSliderLabel = uilabel(app.UpdownscalingTab);
            app.SliceselectionSliderLabel.HorizontalAlignment = 'right';
            app.SliceselectionSliderLabel.Position = [37 72 82 22];
            app.SliceselectionSliderLabel.Text = 'Slice selection';

            % Create Scale_SliceselectionSlider
            app.Scale_SliceselectionSlider = uislider(app.UpdownscalingTab);
            app.Scale_SliceselectionSlider.ValueChangedFcn = createCallbackFcn(app, @Scale_SliceselectionSliderValueChanged, true);
            app.Scale_SliceselectionSlider.Position = [131 81 658 3];

            % Create Scale_instructions
            app.Scale_instructions = uilabel(app.UpdownscalingTab);
            app.Scale_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Scale_instructions.HorizontalAlignment = 'center';
            app.Scale_instructions.FontWeight = 'bold';
            app.Scale_instructions.Position = [17 720 802 44];
            app.Scale_instructions.Text = {'Instructions: 1) Enter scaling factor, 2) Select data type, 3) To visualize changes, click on the ''Do'' button.'; 'To cancel, click on the ''Undo'' button. To keep changes and update the history log, click on the ''Save'' button.'};

            % Create Scale_DoButton
            app.Scale_DoButton = uibutton(app.UpdownscalingTab, 'push');
            app.Scale_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Scale_DoButtonPushed, true);
            app.Scale_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Scale_DoButton.FontSize = 14;
            app.Scale_DoButton.FontWeight = 'bold';
            app.Scale_DoButton.FontColor = [1 1 1];
            app.Scale_DoButton.Position = [302 14 71 25];
            app.Scale_DoButton.Text = 'Do';

            % Create Scale_UndoButton
            app.Scale_UndoButton = uibutton(app.UpdownscalingTab, 'push');
            app.Scale_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Scale_UndoButtonPushed, true);
            app.Scale_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Scale_UndoButton.FontSize = 14;
            app.Scale_UndoButton.FontWeight = 'bold';
            app.Scale_UndoButton.FontColor = [1 1 1];
            app.Scale_UndoButton.Enable = 'off';
            app.Scale_UndoButton.Position = [383 14 71 25];
            app.Scale_UndoButton.Text = 'Undo';

            % Create Scale_SaveButton
            app.Scale_SaveButton = uibutton(app.UpdownscalingTab, 'push');
            app.Scale_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Scale_SaveButtonPushed, true);
            app.Scale_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Scale_SaveButton.FontSize = 14;
            app.Scale_SaveButton.FontWeight = 'bold';
            app.Scale_SaveButton.FontColor = [1 1 1];
            app.Scale_SaveButton.Enable = 'off';
            app.Scale_SaveButton.Position = [464 14 71 25];
            app.Scale_SaveButton.Text = 'Save';

            % Create BackgroundvalueEditField_2Label
            app.BackgroundvalueEditField_2Label = uilabel(app.UpdownscalingTab);
            app.BackgroundvalueEditField_2Label.HorizontalAlignment = 'right';
            app.BackgroundvalueEditField_2Label.Enable = 'off';
            app.BackgroundvalueEditField_2Label.Position = [662 547 102 22];
            app.BackgroundvalueEditField_2Label.Text = 'Background value';

            % Create Scale_BackgroundvalueEditField
            app.Scale_BackgroundvalueEditField = uieditfield(app.UpdownscalingTab, 'numeric');
            app.Scale_BackgroundvalueEditField.Enable = 'off';
            app.Scale_BackgroundvalueEditField.Tooltip = {'Only required for Data type being ''Label'''};
            app.Scale_BackgroundvalueEditField.Position = [779 547 40 22];

            % Create ImagequalityTab
            app.ImagequalityTab = uitab(app.TabGroup);
            app.ImagequalityTab.Title = 'Image quality';

            % Create Quality_instructions
            app.Quality_instructions = uilabel(app.ImagequalityTab);
            app.Quality_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Quality_instructions.HorizontalAlignment = 'center';
            app.Quality_instructions.FontWeight = 'bold';
            app.Quality_instructions.Position = [17 720 802 44];
            app.Quality_instructions.Text = {'Instructions: You can perform image quality calculations from the menu Calculate an any step.'; 'In this tab, you can calculate other metrics that require parameters you can define here.'};

            % Create NumberofphaseincludingporebackgroundLabel
            app.NumberofphaseincludingporebackgroundLabel = uilabel(app.ImagequalityTab);
            app.NumberofphaseincludingporebackgroundLabel.HorizontalAlignment = 'right';
            app.NumberofphaseincludingporebackgroundLabel.Position = [17 657 252 22];
            app.NumberofphaseincludingporebackgroundLabel.Text = 'Number of phase (including pore/background)';

            % Create Quality_numberofphase
            app.Quality_numberofphase = uieditfield(app.ImagequalityTab, 'numeric');
            app.Quality_numberofphase.Limits = [2 Inf];
            app.Quality_numberofphase.RoundFractionalValues = 'on';
            app.Quality_numberofphase.Position = [284 657 50 22];
            app.Quality_numberofphase.Value = 2;

            % Create Quality_Text1
            app.Quality_Text1 = uilabel(app.ImagequalityTab);
            app.Quality_Text1.FontSize = 14;
            app.Quality_Text1.FontWeight = 'bold';
            app.Quality_Text1.Position = [17 684 334 22];
            app.Quality_Text1.Text = 'Separability criterion (for non-segmented image)';

            % Create Quality_CalculateSeparabilityButton
            app.Quality_CalculateSeparabilityButton = uibutton(app.ImagequalityTab, 'push');
            app.Quality_CalculateSeparabilityButton.ButtonPushedFcn = createCallbackFcn(app, @Quality_CalculateSeparabilityButtonPushed, true);
            app.Quality_CalculateSeparabilityButton.Position = [17 630 146 22];
            app.Quality_CalculateSeparabilityButton.Text = 'Plot separability criterion';

            % Create ContrastcorrectionTab
            app.ContrastcorrectionTab = uitab(app.TabGroup);
            app.ContrastcorrectionTab.Title = 'Contrast correction';

            % Create Contrast_UIAxes_before
            app.Contrast_UIAxes_before = uiaxes(app.ContrastcorrectionTab);
            title(app.Contrast_UIAxes_before, 'Before contrast correction. View normal to Axe 3')
            xlabel(app.Contrast_UIAxes_before, '2nd Axis')
            ylabel(app.Contrast_UIAxes_before, '1st Axis')
            app.Contrast_UIAxes_before.Position = [17 64 400 250];

            % Create Contrast_UIAxes_after
            app.Contrast_UIAxes_after = uiaxes(app.ContrastcorrectionTab);
            title(app.Contrast_UIAxes_after, 'After contrast correction. View normal to Axe 3')
            xlabel(app.Contrast_UIAxes_after, '2nd Axis')
            ylabel(app.Contrast_UIAxes_after, '1st Axis')
            app.Contrast_UIAxes_after.Position = [419 64 400 250];

            % Create VolumepercentthresholdforthehighvaluesEditFieldLabel
            app.VolumepercentthresholdforthehighvaluesEditFieldLabel = uilabel(app.ContrastcorrectionTab);
            app.VolumepercentthresholdforthehighvaluesEditFieldLabel.HorizontalAlignment = 'right';
            app.VolumepercentthresholdforthehighvaluesEditFieldLabel.Position = [17 657 245 22];
            app.VolumepercentthresholdforthehighvaluesEditFieldLabel.Text = 'Volume percent threshold for the high values';

            % Create Contrast_VolumepercentthresholdforthehighvaluesEditField
            app.Contrast_VolumepercentthresholdforthehighvaluesEditField = uieditfield(app.ContrastcorrectionTab, 'numeric');
            app.Contrast_VolumepercentthresholdforthehighvaluesEditField.Limits = [0 100];
            app.Contrast_VolumepercentthresholdforthehighvaluesEditField.Position = [277 657 40 22];
            app.Contrast_VolumepercentthresholdforthehighvaluesEditField.Value = 1;

            % Create VolumepercentthresholdforthelowvaluesEditFieldLabel
            app.VolumepercentthresholdforthelowvaluesEditFieldLabel = uilabel(app.ContrastcorrectionTab);
            app.VolumepercentthresholdforthelowvaluesEditFieldLabel.HorizontalAlignment = 'right';
            app.VolumepercentthresholdforthelowvaluesEditFieldLabel.Position = [17 630 240 22];
            app.VolumepercentthresholdforthelowvaluesEditFieldLabel.Text = 'Volume percent threshold for the low values';

            % Create Contrast_VolumepercentthresholdforthelowvaluesEditField
            app.Contrast_VolumepercentthresholdforthelowvaluesEditField = uieditfield(app.ContrastcorrectionTab, 'numeric');
            app.Contrast_VolumepercentthresholdforthelowvaluesEditField.Limits = [0 100];
            app.Contrast_VolumepercentthresholdforthelowvaluesEditField.Position = [277 630 40 22];
            app.Contrast_VolumepercentthresholdforthelowvaluesEditField.Value = 1;

            % Create RescaleimageusingthisnumberofvalueLabel
            app.RescaleimageusingthisnumberofvalueLabel = uilabel(app.ContrastcorrectionTab);
            app.RescaleimageusingthisnumberofvalueLabel.HorizontalAlignment = 'right';
            app.RescaleimageusingthisnumberofvalueLabel.Position = [17 526 230 22];
            app.RescaleimageusingthisnumberofvalueLabel.Text = 'Rescale image using this number of value';

            % Create Contrast_CustomrangerescaleusingthisnumberofvalueEditField
            app.Contrast_CustomrangerescaleusingthisnumberofvalueEditField = uieditfield(app.ContrastcorrectionTab, 'numeric');
            app.Contrast_CustomrangerescaleusingthisnumberofvalueEditField.Limits = [1 Inf];
            app.Contrast_CustomrangerescaleusingthisnumberofvalueEditField.Position = [267 526 50 22];
            app.Contrast_CustomrangerescaleusingthisnumberofvalueEditField.Value = 1024;

            % Create Contrast_Text1
            app.Contrast_Text1 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text1.FontSize = 14;
            app.Contrast_Text1.FontWeight = 'bold';
            app.Contrast_Text1.Position = [17 684 127 22];
            app.Contrast_Text1.Text = 'Saturate extremes';

            % Create Contrast_instructions
            app.Contrast_instructions = uilabel(app.ContrastcorrectionTab);
            app.Contrast_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Contrast_instructions.HorizontalAlignment = 'center';
            app.Contrast_instructions.FontWeight = 'bold';
            app.Contrast_instructions.Position = [17 720 802 44];
            app.Contrast_instructions.Text = {'Instructions: You can perform operations in any order. To visualize changes, click on the ''Do'' button.'; 'To cancel, click on the ''Undo'' button. To keep changes and update the history log, click on the ''Save'' button.'};

            % Create Contrast_Burn_DoButton
            app.Contrast_Burn_DoButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Burn_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Burn_DoButtonPushed, true);
            app.Contrast_Burn_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Contrast_Burn_DoButton.FontSize = 14;
            app.Contrast_Burn_DoButton.FontWeight = 'bold';
            app.Contrast_Burn_DoButton.FontColor = [1 1 1];
            app.Contrast_Burn_DoButton.Position = [84 600 71 25];
            app.Contrast_Burn_DoButton.Text = 'Do';

            % Create Contrast_Burn_UndoButton
            app.Contrast_Burn_UndoButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Burn_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Burn_UndoButtonPushed, true);
            app.Contrast_Burn_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Contrast_Burn_UndoButton.FontSize = 14;
            app.Contrast_Burn_UndoButton.FontWeight = 'bold';
            app.Contrast_Burn_UndoButton.FontColor = [1 1 1];
            app.Contrast_Burn_UndoButton.Enable = 'off';
            app.Contrast_Burn_UndoButton.Position = [165 600 71 25];
            app.Contrast_Burn_UndoButton.Text = 'Undo';

            % Create Contrast_Burn_SaveButton
            app.Contrast_Burn_SaveButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Burn_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Burn_SaveButtonPushed, true);
            app.Contrast_Burn_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Contrast_Burn_SaveButton.FontSize = 14;
            app.Contrast_Burn_SaveButton.FontWeight = 'bold';
            app.Contrast_Burn_SaveButton.FontColor = [1 1 1];
            app.Contrast_Burn_SaveButton.Enable = 'off';
            app.Contrast_Burn_SaveButton.Position = [246 600 71 25];
            app.Contrast_Burn_SaveButton.Text = 'Save';

            % Create Contrast_Text2
            app.Contrast_Text2 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text2.FontSize = 14;
            app.Contrast_Text2.FontWeight = 'bold';
            app.Contrast_Text2.Position = [17 553 100 22];
            app.Contrast_Text2.Text = 'Custom range';

            % Create Contrast_Customrange_DoButton
            app.Contrast_Customrange_DoButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Customrange_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Customrange_DoButtonPushed, true);
            app.Contrast_Customrange_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Contrast_Customrange_DoButton.FontSize = 14;
            app.Contrast_Customrange_DoButton.FontWeight = 'bold';
            app.Contrast_Customrange_DoButton.FontColor = [1 1 1];
            app.Contrast_Customrange_DoButton.Position = [84 496 71 25];
            app.Contrast_Customrange_DoButton.Text = 'Do';

            % Create Contrast_Customrange_UndoButton
            app.Contrast_Customrange_UndoButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Customrange_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Customrange_UndoButtonPushed, true);
            app.Contrast_Customrange_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Contrast_Customrange_UndoButton.FontSize = 14;
            app.Contrast_Customrange_UndoButton.FontWeight = 'bold';
            app.Contrast_Customrange_UndoButton.FontColor = [1 1 1];
            app.Contrast_Customrange_UndoButton.Enable = 'off';
            app.Contrast_Customrange_UndoButton.Position = [165 496 71 25];
            app.Contrast_Customrange_UndoButton.Text = 'Undo';

            % Create Contrast_Customrange_SaveButton
            app.Contrast_Customrange_SaveButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Customrange_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Customrange_SaveButtonPushed, true);
            app.Contrast_Customrange_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Contrast_Customrange_SaveButton.FontSize = 14;
            app.Contrast_Customrange_SaveButton.FontWeight = 'bold';
            app.Contrast_Customrange_SaveButton.FontColor = [1 1 1];
            app.Contrast_Customrange_SaveButton.Enable = 'off';
            app.Contrast_Customrange_SaveButton.Position = [246 496 71 25];
            app.Contrast_Customrange_SaveButton.Text = 'Save';

            % Create SliceselectionSliderLabel_2
            app.SliceselectionSliderLabel_2 = uilabel(app.ContrastcorrectionTab);
            app.SliceselectionSliderLabel_2.HorizontalAlignment = 'right';
            app.SliceselectionSliderLabel_2.Position = [17 30 82 22];
            app.SliceselectionSliderLabel_2.Text = 'Slice selection';

            % Create Contrast_SliceselectionSlider
            app.Contrast_SliceselectionSlider = uislider(app.ContrastcorrectionTab);
            app.Contrast_SliceselectionSlider.ValueChangedFcn = createCallbackFcn(app, @Contrast_SliceselectionSliderValueChanged, true);
            app.Contrast_SliceselectionSlider.Position = [120 39 688 3];

            % Create Contrast_Text3
            app.Contrast_Text3 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text3.FontAngle = 'italic';
            app.Contrast_Text3.Position = [18 324 514 22];
            app.Contrast_Text3.Text = 'Colormap of axes below is set to [global minimum, global maximum] of their respective image.';

            % Create Contrast_TableWhere
            app.Contrast_TableWhere = uitable(app.ContrastcorrectionTab);
            app.Contrast_TableWhere.ColumnName = {'Axis'; 'From'; 'To'};
            app.Contrast_TableWhere.RowName = {};
            app.Contrast_TableWhere.ColumnEditable = [false true true];
            app.Contrast_TableWhere.CellEditCallback = createCallbackFcn(app, @Contrast_TableWhereCellEdit, true);
            app.Contrast_TableWhere.Position = [369 462 204 115];

            % Create Contrast_Text5
            app.Contrast_Text5 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text5.FontSize = 14;
            app.Contrast_Text5.FontWeight = 'bold';
            app.Contrast_Text5.Position = [369 684 228 22];
            app.Contrast_Text5.Text = 'Advanced contrast enhancement';

            % Create Contrast_Text4
            app.Contrast_Text4 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text4.Position = [369 657 421 22];
            app.Contrast_Text4.Text = 'Choose where to apply contrast correction. It will be performed slice per slice';

            % Create Contrast_Text6
            app.Contrast_Text6 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text6.Position = [369 632 410 22];
            app.Contrast_Text6.Text = 'along Axe 3, within the bounds defined in the left table below. Then choose';

            % Create Contrast_TableParameters
            app.Contrast_TableParameters = uitable(app.ContrastcorrectionTab);
            app.Contrast_TableParameters.ColumnName = {'Parameters'; 'Value'};
            app.Contrast_TableParameters.RowName = {};
            app.Contrast_TableParameters.ColumnEditable = [false true];
            app.Contrast_TableParameters.CellEditCallback = createCallbackFcn(app, @Contrast_TableParametersCellEdit, true);
            app.Contrast_TableParameters.Position = [585 424 234 153];

            % Create Contrast_Text7
            app.Contrast_Text7 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text7.Position = [369 607 447 22];
            app.Contrast_Text7.Text = 'what operation you want to apply in the dropdown button, and then its parameters';

            % Create Contrast_Text8
            app.Contrast_Text8 = uilabel(app.ContrastcorrectionTab);
            app.Contrast_Text8.Position = [369 582 441 22];
            app.Contrast_Text8.Text = 'using the right table below. Put the mouse cursor over the right table to see help.';

            % Create Contrast_Advanced_DoButton
            app.Contrast_Advanced_DoButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Advanced_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Advanced_DoButtonPushed, true);
            app.Contrast_Advanced_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Contrast_Advanced_DoButton.FontSize = 14;
            app.Contrast_Advanced_DoButton.FontWeight = 'bold';
            app.Contrast_Advanced_DoButton.FontColor = [1 1 1];
            app.Contrast_Advanced_DoButton.Position = [586 391 71 25];
            app.Contrast_Advanced_DoButton.Text = 'Do';

            % Create Contrast_Advanced_UndoButton
            app.Contrast_Advanced_UndoButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Advanced_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Advanced_UndoButtonPushed, true);
            app.Contrast_Advanced_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Contrast_Advanced_UndoButton.FontSize = 14;
            app.Contrast_Advanced_UndoButton.FontWeight = 'bold';
            app.Contrast_Advanced_UndoButton.FontColor = [1 1 1];
            app.Contrast_Advanced_UndoButton.Enable = 'off';
            app.Contrast_Advanced_UndoButton.Position = [667 391 71 25];
            app.Contrast_Advanced_UndoButton.Text = 'Undo';

            % Create Contrast_Advanced_SaveButton
            app.Contrast_Advanced_SaveButton = uibutton(app.ContrastcorrectionTab, 'push');
            app.Contrast_Advanced_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Contrast_Advanced_SaveButtonPushed, true);
            app.Contrast_Advanced_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Contrast_Advanced_SaveButton.FontSize = 14;
            app.Contrast_Advanced_SaveButton.FontWeight = 'bold';
            app.Contrast_Advanced_SaveButton.FontColor = [1 1 1];
            app.Contrast_Advanced_SaveButton.Enable = 'off';
            app.Contrast_Advanced_SaveButton.Position = [748 391 71 25];
            app.Contrast_Advanced_SaveButton.Text = 'Save';

            % Create OperationLabel
            app.OperationLabel = uilabel(app.ContrastcorrectionTab);
            app.OperationLabel.HorizontalAlignment = 'right';
            app.OperationLabel.Position = [385 424 58 22];
            app.OperationLabel.Text = 'Operation';

            % Create Contrast_Method_DropDown
            app.Contrast_Method_DropDown = uidropdown(app.ContrastcorrectionTab);
            app.Contrast_Method_DropDown.Items = {'Adjust contrast', 'Histogram equalization', 'Contrast-limited adaptive histogram equalization'};
            app.Contrast_Method_DropDown.ValueChangedFcn = createCallbackFcn(app, @Contrast_Method_DropDownValueChanged, true);
            app.Contrast_Method_DropDown.Position = [458 424 100 22];
            app.Contrast_Method_DropDown.Value = 'Adjust contrast';

            % Create FilteringTab
            app.FilteringTab = uitab(app.TabGroup);
            app.FilteringTab.Tooltip = {''};
            app.FilteringTab.Title = 'Filtering';

            % Create Filtering_UIAxes_before
            app.Filtering_UIAxes_before = uiaxes(app.FilteringTab);
            title(app.Filtering_UIAxes_before, 'Before filtering. View normal to Axe 3')
            xlabel(app.Filtering_UIAxes_before, '2nd Axis')
            ylabel(app.Filtering_UIAxes_before, '1st Axis')
            app.Filtering_UIAxes_before.Position = [17 64 400 250];

            % Create Filetring_UIAxes_after
            app.Filetring_UIAxes_after = uiaxes(app.FilteringTab);
            title(app.Filetring_UIAxes_after, 'After filtering. View normal to Axe 3')
            xlabel(app.Filetring_UIAxes_after, '2nd Axis')
            ylabel(app.Filetring_UIAxes_after, '1st Axis')
            app.Filetring_UIAxes_after.Position = [419 64 400 250];

            % Create Filtering_instructions
            app.Filtering_instructions = uilabel(app.FilteringTab);
            app.Filtering_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Filtering_instructions.HorizontalAlignment = 'center';
            app.Filtering_instructions.FontWeight = 'bold';
            app.Filtering_instructions.Position = [17 720 802 44];
            app.Filtering_instructions.Text = {'Instructions: You can perform operations in any order. To visualize changes, click on the ''Do'' button.'; 'To cancel, click on the ''Undo'' button. To keep changes and update the history log, click on the ''Save'' button.'};

            % Create Filtering_Text1
            app.Filtering_Text1 = uilabel(app.FilteringTab);
            app.Filtering_Text1.FontSize = 14;
            app.Filtering_Text1.FontWeight = 'bold';
            app.Filtering_Text1.Tooltip = {'See: https://www.mathworks.com/help/images/ref/imdiffusefilt.html and  Perona, P., and J. Malik. "Scale-space and edge detection using anisotropic diffusion." IEEE® Transactions on Pattern Analysis and Machine Intelligence. Vol. 12, No. 7, July 1990, pp. 629–639.'};
            app.Filtering_Text1.Position = [17 684 180 22];
            app.Filtering_Text1.Text = 'Anisotropic diffusion filter';

            % Create ApplyDropDownLabel
            app.ApplyDropDownLabel = uilabel(app.FilteringTab);
            app.ApplyDropDownLabel.HorizontalAlignment = 'right';
            app.ApplyDropDownLabel.Position = [17 657 36 22];
            app.ApplyDropDownLabel.Text = 'Apply';

            % Create Filtering_AnisotropicFilter_2D3DDropDown
            app.Filtering_AnisotropicFilter_2D3DDropDown = uidropdown(app.FilteringTab);
            app.Filtering_AnisotropicFilter_2D3DDropDown.Items = {'Slice per slice, along Axis 3', 'Whole volume'};
            app.Filtering_AnisotropicFilter_2D3DDropDown.Position = [68 657 242 22];
            app.Filtering_AnisotropicFilter_2D3DDropDown.Value = 'Slice per slice, along Axis 3';

            % Create ConductionmethodDropDownLabel
            app.ConductionmethodDropDownLabel = uilabel(app.FilteringTab);
            app.ConductionmethodDropDownLabel.HorizontalAlignment = 'right';
            app.ConductionmethodDropDownLabel.Position = [17 630 110 22];
            app.ConductionmethodDropDownLabel.Text = 'Conduction method';

            % Create Filtering_AnisotropicFilter_MethodDropDown
            app.Filtering_AnisotropicFilter_MethodDropDown = uidropdown(app.FilteringTab);
            app.Filtering_AnisotropicFilter_MethodDropDown.Items = {'exponential', 'quadratic'};
            app.Filtering_AnisotropicFilter_MethodDropDown.Tooltip = {'Exponential diffusion favors high-contrast edges over low-contrast edges. Quadratic diffusion favors wide regions over smaller regions.'};
            app.Filtering_AnisotropicFilter_MethodDropDown.Position = [142 630 168 22];
            app.Filtering_AnisotropicFilter_MethodDropDown.Value = 'exponential';

            % Create ConnectivityDropDownLabel
            app.ConnectivityDropDownLabel = uilabel(app.FilteringTab);
            app.ConnectivityDropDownLabel.HorizontalAlignment = 'right';
            app.ConnectivityDropDownLabel.Position = [17 603 71 22];
            app.ConnectivityDropDownLabel.Text = 'Connectivity';

            % Create Filtering_AnisotropicFilter_ConnDropDown
            app.Filtering_AnisotropicFilter_ConnDropDown = uidropdown(app.FilteringTab);
            app.Filtering_AnisotropicFilter_ConnDropDown.Items = {'maximal', 'minimal'};
            app.Filtering_AnisotropicFilter_ConnDropDown.Tooltip = {'''maximal'' — Considers 8 nearest neighbors for 2-D images, and 26 nearest neighbors for 3-D images. ''minimal'' — Considers 4 nearest neighbors for 2-D images, and 6 nearest neighbors for 3-D images.'};
            app.Filtering_AnisotropicFilter_ConnDropDown.Position = [103 603 207 22];
            app.Filtering_AnisotropicFilter_ConnDropDown.Value = 'maximal';

            % Create NumberofiterationEditFieldLabel
            app.NumberofiterationEditFieldLabel = uilabel(app.FilteringTab);
            app.NumberofiterationEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofiterationEditFieldLabel.Enable = 'off';
            app.NumberofiterationEditFieldLabel.Position = [17 576 108 22];
            app.NumberofiterationEditFieldLabel.Text = 'Number of iteration';

            % Create Filtering_AnisotropicFilter_IterEditField
            app.Filtering_AnisotropicFilter_IterEditField = uieditfield(app.FilteringTab, 'numeric');
            app.Filtering_AnisotropicFilter_IterEditField.Limits = [1 Inf];
            app.Filtering_AnisotropicFilter_IterEditField.ValueChangedFcn = createCallbackFcn(app, @Filtering_AnisotropicFilter_IterEditFieldValueChanged, true);
            app.Filtering_AnisotropicFilter_IterEditField.Enable = 'off';
            app.Filtering_AnisotropicFilter_IterEditField.Position = [186 576 124 22];
            app.Filtering_AnisotropicFilter_IterEditField.Value = 5;

            % Create Filtering_AnisotropicGradient_Table
            app.Filtering_AnisotropicGradient_Table = uitable(app.FilteringTab);
            app.Filtering_AnisotropicGradient_Table.ColumnName = {'Iteration'; 'Gradient threshold'};
            app.Filtering_AnisotropicGradient_Table.RowName = {};
            app.Filtering_AnisotropicGradient_Table.ColumnEditable = [false true];
            app.Filtering_AnisotropicGradient_Table.CellEditCallback = createCallbackFcn(app, @Filtering_AnisotropicGradient_TableCellEdit, true);
            app.Filtering_AnisotropicGradient_Table.Tooltip = {'The value of Gradient threshold controls the conduction process by classifying gradient values as an actual edge or as noise. Increasing the value of Gradient threshold smooths the image more. The default value is 10% of the dynamic range of the image.'};
            app.Filtering_AnisotropicGradient_Table.Enable = 'off';
            app.Filtering_AnisotropicGradient_Table.Position = [17 447 293 124];

            % Create Filtering_Anisotropic_DoButton
            app.Filtering_Anisotropic_DoButton = uibutton(app.FilteringTab, 'push');
            app.Filtering_Anisotropic_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Filtering_Anisotropic_DoButtonPushed, true);
            app.Filtering_Anisotropic_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Filtering_Anisotropic_DoButton.FontSize = 14;
            app.Filtering_Anisotropic_DoButton.FontWeight = 'bold';
            app.Filtering_Anisotropic_DoButton.FontColor = [1 1 1];
            app.Filtering_Anisotropic_DoButton.Position = [77 390 71 25];
            app.Filtering_Anisotropic_DoButton.Text = 'Do';

            % Create Filtering_Anisotropic_UndoButton
            app.Filtering_Anisotropic_UndoButton = uibutton(app.FilteringTab, 'push');
            app.Filtering_Anisotropic_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Filtering_Anisotropic_UndoButtonPushed, true);
            app.Filtering_Anisotropic_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Filtering_Anisotropic_UndoButton.FontSize = 14;
            app.Filtering_Anisotropic_UndoButton.FontWeight = 'bold';
            app.Filtering_Anisotropic_UndoButton.FontColor = [1 1 1];
            app.Filtering_Anisotropic_UndoButton.Enable = 'off';
            app.Filtering_Anisotropic_UndoButton.Position = [158 390 71 25];
            app.Filtering_Anisotropic_UndoButton.Text = 'Undo';

            % Create Filtering_Anisotropic_SaveButton
            app.Filtering_Anisotropic_SaveButton = uibutton(app.FilteringTab, 'push');
            app.Filtering_Anisotropic_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Filtering_Anisotropic_SaveButtonPushed, true);
            app.Filtering_Anisotropic_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Filtering_Anisotropic_SaveButton.FontSize = 14;
            app.Filtering_Anisotropic_SaveButton.FontWeight = 'bold';
            app.Filtering_Anisotropic_SaveButton.FontColor = [1 1 1];
            app.Filtering_Anisotropic_SaveButton.Enable = 'off';
            app.Filtering_Anisotropic_SaveButton.Position = [239 390 71 25];
            app.Filtering_Anisotropic_SaveButton.Text = 'Save';

            % Create Filtering_Text3
            app.Filtering_Text3 = uilabel(app.FilteringTab);
            app.Filtering_Text3.FontAngle = 'italic';
            app.Filtering_Text3.Position = [18 324 514 22];
            app.Filtering_Text3.Text = 'Colormap of axes below is set to [global minimum, global maximum] of their respective image.';

            % Create SliceselectionSliderLabel_3
            app.SliceselectionSliderLabel_3 = uilabel(app.FilteringTab);
            app.SliceselectionSliderLabel_3.HorizontalAlignment = 'right';
            app.SliceselectionSliderLabel_3.Position = [17 30 82 22];
            app.SliceselectionSliderLabel_3.Text = 'Slice selection';

            % Create Filtering_SliceselectionSlider
            app.Filtering_SliceselectionSlider = uislider(app.FilteringTab);
            app.Filtering_SliceselectionSlider.ValueChangedFcn = createCallbackFcn(app, @Filtering_SliceselectionSliderValueChanged, true);
            app.Filtering_SliceselectionSlider.Position = [120 39 688 3];

            % Create Filtering_Text2
            app.Filtering_Text2 = uilabel(app.FilteringTab);
            app.Filtering_Text2.FontSize = 14;
            app.Filtering_Text2.FontWeight = 'bold';
            app.Filtering_Text2.Tooltip = {'See: https://www.mathworks.com/help/images/ref/imnlmfilt.html and Buades, A., B. Coll, and J.-M. Morel. "A Non-Local Algorithm for Image Denoising." 2005 IEEE® Computer Society Conference on Computer Vision and Pattern Recognition. Vol. 2, June 2005, pp. 60–65.'};
            app.Filtering_Text2.Position = [340 684 152 22];
            app.Filtering_Text2.Text = 'Non-local means filter';

            % Create Filetring_NLMF_EstimateCheckBox
            app.Filetring_NLMF_EstimateCheckBox = uicheckbox(app.FilteringTab);
            app.Filetring_NLMF_EstimateCheckBox.ValueChangedFcn = createCallbackFcn(app, @Filetring_NLMF_EstimateCheckBoxValueChanged, true);
            app.Filetring_NLMF_EstimateCheckBox.Text = 'Estimate degree of smoothing';
            app.Filetring_NLMF_EstimateCheckBox.Position = [340 657 181 22];
            app.Filetring_NLMF_EstimateCheckBox.Value = true;

            % Create DegreeofsmoothingEditFieldLabel
            app.DegreeofsmoothingEditFieldLabel = uilabel(app.FilteringTab);
            app.DegreeofsmoothingEditFieldLabel.HorizontalAlignment = 'right';
            app.DegreeofsmoothingEditFieldLabel.Position = [340 630 117 22];
            app.DegreeofsmoothingEditFieldLabel.Text = 'Degree of smoothing';

            % Create Filetring_NLMF_DegreeofsmoothingEditField
            app.Filetring_NLMF_DegreeofsmoothingEditField = uieditfield(app.FilteringTab, 'numeric');
            app.Filetring_NLMF_DegreeofsmoothingEditField.Limits = [0 Inf];
            app.Filetring_NLMF_DegreeofsmoothingEditField.Enable = 'off';
            app.Filetring_NLMF_DegreeofsmoothingEditField.Tooltip = {'As this value increases, the smoothing in the resulting image increases.'};
            app.Filetring_NLMF_DegreeofsmoothingEditField.Position = [472 630 69 22];
            app.Filetring_NLMF_DegreeofsmoothingEditField.Value = 5;

            % Create SearchwindowsizeEditFieldLabel
            app.SearchwindowsizeEditFieldLabel = uilabel(app.FilteringTab);
            app.SearchwindowsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.SearchwindowsizeEditFieldLabel.Position = [340 603 112 22];
            app.SearchwindowsizeEditFieldLabel.Text = 'Search window size';

            % Create Filetring_NLMF_SearchwindowsizeEditField
            app.Filetring_NLMF_SearchwindowsizeEditField = uieditfield(app.FilteringTab, 'numeric');
            app.Filetring_NLMF_SearchwindowsizeEditField.Limits = [1 Inf];
            app.Filetring_NLMF_SearchwindowsizeEditField.RoundFractionalValues = 'on';
            app.Filetring_NLMF_SearchwindowsizeEditField.ValueChangedFcn = createCallbackFcn(app, @Filetring_NLMF_SearchwindowsizeEditFieldValueChanged, true);
            app.Filetring_NLMF_SearchwindowsizeEditField.Tooltip = {'Odd-valued positive integer s. The search for similar neighborhoods to a pixel is limited to the s-by-s region surrounding that pixel.'};
            app.Filetring_NLMF_SearchwindowsizeEditField.Position = [493 603 48 22];
            app.Filetring_NLMF_SearchwindowsizeEditField.Value = 21;

            % Create ComparisonwindowsizeEditFieldLabel
            app.ComparisonwindowsizeEditFieldLabel = uilabel(app.FilteringTab);
            app.ComparisonwindowsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.ComparisonwindowsizeEditFieldLabel.Position = [340 576 138 22];
            app.ComparisonwindowsizeEditFieldLabel.Text = 'Comparison window size';

            % Create Filetring_NLMF_ComparisonwindowsizeEditField
            app.Filetring_NLMF_ComparisonwindowsizeEditField = uieditfield(app.FilteringTab, 'numeric');
            app.Filetring_NLMF_ComparisonwindowsizeEditField.Limits = [1 Inf];
            app.Filetring_NLMF_ComparisonwindowsizeEditField.RoundFractionalValues = 'on';
            app.Filetring_NLMF_ComparisonwindowsizeEditField.ValueChangedFcn = createCallbackFcn(app, @Filetring_NLMF_SearchwindowsizeEditFieldValueChanged, true);
            app.Filetring_NLMF_ComparisonwindowsizeEditField.Tooltip = {'Odd-valued positive integer c. The non-local means filter computes similarity weights using the c-by-c neighborhood surrounding pixels. Must be less than or equal to Search window size.'};
            app.Filetring_NLMF_ComparisonwindowsizeEditField.Position = [493 576 48 22];
            app.Filetring_NLMF_ComparisonwindowsizeEditField.Value = 5;

            % Create Filtering_NLMF_DoButton
            app.Filtering_NLMF_DoButton = uibutton(app.FilteringTab, 'push');
            app.Filtering_NLMF_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Filtering_NLMF_DoButtonPushed, true);
            app.Filtering_NLMF_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Filtering_NLMF_DoButton.FontSize = 14;
            app.Filtering_NLMF_DoButton.FontWeight = 'bold';
            app.Filtering_NLMF_DoButton.FontColor = [1 1 1];
            app.Filtering_NLMF_DoButton.Position = [340 546 71 25];
            app.Filtering_NLMF_DoButton.Text = 'Do';

            % Create Filtering_NLMF_UndoButton
            app.Filtering_NLMF_UndoButton = uibutton(app.FilteringTab, 'push');
            app.Filtering_NLMF_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Filtering_NLMF_UndoButtonPushed, true);
            app.Filtering_NLMF_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Filtering_NLMF_UndoButton.FontSize = 14;
            app.Filtering_NLMF_UndoButton.FontWeight = 'bold';
            app.Filtering_NLMF_UndoButton.FontColor = [1 1 1];
            app.Filtering_NLMF_UndoButton.Enable = 'off';
            app.Filtering_NLMF_UndoButton.Position = [340 516 71 25];
            app.Filtering_NLMF_UndoButton.Text = 'Undo';

            % Create Filtering_NLMF_SaveButton
            app.Filtering_NLMF_SaveButton = uibutton(app.FilteringTab, 'push');
            app.Filtering_NLMF_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Filtering_NLMF_SaveButtonPushed, true);
            app.Filtering_NLMF_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Filtering_NLMF_SaveButton.FontSize = 14;
            app.Filtering_NLMF_SaveButton.FontWeight = 'bold';
            app.Filtering_NLMF_SaveButton.FontColor = [1 1 1];
            app.Filtering_NLMF_SaveButton.Enable = 'off';
            app.Filtering_NLMF_SaveButton.Position = [340 486 71 25];
            app.Filtering_NLMF_SaveButton.Text = 'Save';

            % Create Filtering_Text4
            app.Filtering_Text4 = uilabel(app.FilteringTab);
            app.Filtering_Text4.FontAngle = 'italic';
            app.Filtering_Text4.Position = [340 459 183 22];
            app.Filtering_Text4.Text = 'Apply slice per slice, along Axis 3';

            % Create Filtering_previewCheckBox
            app.Filtering_previewCheckBox = uicheckbox(app.FilteringTab);
            app.Filtering_previewCheckBox.Text = 'Preview (apply filter only at the current slice)';
            app.Filtering_previewCheckBox.Position = [559 324 260 22];

            % Create Filtering_ProgressnotrunningLabel
            app.Filtering_ProgressnotrunningLabel = uilabel(app.FilteringTab);
            app.Filtering_ProgressnotrunningLabel.FontAngle = 'italic';
            app.Filtering_ProgressnotrunningLabel.Position = [559 345 233 22];
            app.Filtering_ProgressnotrunningLabel.Text = 'Progress: not running';

            % Create Filtering_Anisotropy_EstimateCheckBox
            app.Filtering_Anisotropy_EstimateCheckBox = uicheckbox(app.FilteringTab);
            app.Filtering_Anisotropy_EstimateCheckBox.ValueChangedFcn = createCallbackFcn(app, @Filtering_Anisotropy_EstimateCheckBoxValueChanged, true);
            app.Filtering_Anisotropy_EstimateCheckBox.Text = 'Estimate gradient threshold and number of iteration';
            app.Filtering_Anisotropy_EstimateCheckBox.Position = [18 426 299 22];
            app.Filtering_Anisotropy_EstimateCheckBox.Value = true;

            % Create SegmentationTab
            app.SegmentationTab = uitab(app.TabGroup);
            app.SegmentationTab.Title = 'Segmentation';

            % Create Segmentation_UIAxes
            app.Segmentation_UIAxes = uiaxes(app.SegmentationTab);
            title(app.Segmentation_UIAxes, 'View normal to Axe 3')
            xlabel(app.Segmentation_UIAxes, '2nd Axis')
            ylabel(app.Segmentation_UIAxes, '1st Axis')
            app.Segmentation_UIAxes.Position = [353 179 466 527];

            % Create Segmentation_instructions
            app.Segmentation_instructions = uilabel(app.SegmentationTab);
            app.Segmentation_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Segmentation_instructions.HorizontalAlignment = 'center';
            app.Segmentation_instructions.FontWeight = 'bold';
            app.Segmentation_instructions.Position = [17 720 802 44];
            app.Segmentation_instructions.Text = {'Instructions: Enter voxel size, number of phase and phase table. To apply segmentation, click on the ''Do'' button. '; 'To cancel, click on the ''Undo'' button. To keep changes  and update the history log, click on the ''Save'' button.'; 'To visualize, open ''Microstructure and results visualization'' module, ''Segmented volume'' tab.'};

            % Create Segmentation_Text1
            app.Segmentation_Text1 = uilabel(app.SegmentationTab);
            app.Segmentation_Text1.FontSize = 14;
            app.Segmentation_Text1.FontWeight = 'bold';
            app.Segmentation_Text1.Tooltip = {'See: https://www.mathworks.com/help/images/ref/imdiffusefilt.html and  Perona, P., and J. Malik. "Scale-space and edge detection using anisotropic diffusion." IEEE® Transactions on Pattern Analysis and Machine Intelligence. Vol. 12, No. 7, July 1990, pp. 629–639.'};
            app.Segmentation_Text1.Position = [17 497 117 22];
            app.Segmentation_Text1.Text = 'Global threshold';

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.SegmentationTab);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.Position = [17 470 46 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create Segmentation_Global_MethodDropDown
            app.Segmentation_Global_MethodDropDown = uidropdown(app.SegmentationTab);
            app.Segmentation_Global_MethodDropDown.Items = {'Manual', 'Match expected volume fraction', 'Otsu'};
            app.Segmentation_Global_MethodDropDown.ValueChangedFcn = createCallbackFcn(app, @Segmentation_Global_MethodDropDownValueChanged, true);
            app.Segmentation_Global_MethodDropDown.Position = [78 470 241 22];
            app.Segmentation_Global_MethodDropDown.Value = 'Manual';

            % Create VoxelsizemicrometersEditFieldLabel
            app.VoxelsizemicrometersEditFieldLabel = uilabel(app.SegmentationTab);
            app.VoxelsizemicrometersEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelsizemicrometersEditFieldLabel.Position = [17 684 137 22];
            app.VoxelsizemicrometersEditFieldLabel.Text = 'Voxel size (micrometers)';

            % Create Segmentation_VoxelsizemicrometersEditField
            app.Segmentation_VoxelsizemicrometersEditField = uieditfield(app.SegmentationTab, 'numeric');
            app.Segmentation_VoxelsizemicrometersEditField.Limits = [0 Inf];
            app.Segmentation_VoxelsizemicrometersEditField.Position = [172 684 45 22];
            app.Segmentation_VoxelsizemicrometersEditField.Value = 1;

            % Create Segmentation_Phase_UITable
            app.Segmentation_Phase_UITable = uitable(app.SegmentationTab);
            app.Segmentation_Phase_UITable.ColumnName = {'Phase'; 'Assign value'; 'Expected volume fraction'; 'Obtained'};
            app.Segmentation_Phase_UITable.RowName = {};
            app.Segmentation_Phase_UITable.ColumnEditable = [false true true false];
            app.Segmentation_Phase_UITable.CellEditCallback = createCallbackFcn(app, @Segmentation_Phase_UITableCellEdit, true);
            app.Segmentation_Phase_UITable.Position = [17 533 302 119];

            % Create NumberofphaseEditFieldLabel
            app.NumberofphaseEditFieldLabel = uilabel(app.SegmentationTab);
            app.NumberofphaseEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofphaseEditFieldLabel.Position = [17 657 98 22];
            app.NumberofphaseEditFieldLabel.Text = 'Number of phase';

            % Create Segmentation_NumberofphaseEditField
            app.Segmentation_NumberofphaseEditField = uieditfield(app.SegmentationTab, 'numeric');
            app.Segmentation_NumberofphaseEditField.Limits = [2 Inf];
            app.Segmentation_NumberofphaseEditField.RoundFractionalValues = 'on';
            app.Segmentation_NumberofphaseEditField.ValueChangedFcn = createCallbackFcn(app, @Segmentation_NumberofphaseEditFieldValueChanged, true);
            app.Segmentation_NumberofphaseEditField.Position = [130 657 45 22];
            app.Segmentation_NumberofphaseEditField.Value = 2;

            % Create Segmentation_Global_UITable
            app.Segmentation_Global_UITable = uitable(app.SegmentationTab);
            app.Segmentation_Global_UITable.ColumnName = {'Phase'; 'Assign value'; 'Threshold <='};
            app.Segmentation_Global_UITable.RowName = {};
            app.Segmentation_Global_UITable.ColumnEditable = [false false true];
            app.Segmentation_Global_UITable.CellEditCallback = createCallbackFcn(app, @Segmentation_Global_UITableCellEdit, true);
            app.Segmentation_Global_UITable.Position = [17 346 302 119];

            % Create Segmentation_CompareOtsuandexpectedvolumefractionsButton
            app.Segmentation_CompareOtsuandexpectedvolumefractionsButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_CompareOtsuandexpectedvolumefractionsButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_CompareOtsuandexpectedvolumefractionsButtonPushed, true);
            app.Segmentation_CompareOtsuandexpectedvolumefractionsButton.Position = [17 319 263 22];
            app.Segmentation_CompareOtsuandexpectedvolumefractionsButton.Text = 'Compare Otsu and expected volume fractions';

            % Create Segmentation_Global_DoButton
            app.Segmentation_Global_DoButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_Global_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_Global_DoButtonPushed, true);
            app.Segmentation_Global_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Segmentation_Global_DoButton.FontSize = 14;
            app.Segmentation_Global_DoButton.FontWeight = 'bold';
            app.Segmentation_Global_DoButton.FontColor = [1 1 1];
            app.Segmentation_Global_DoButton.Position = [86 289 71 25];
            app.Segmentation_Global_DoButton.Text = 'Do';

            % Create Segmentation_Global_UndoButton
            app.Segmentation_Global_UndoButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_Global_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_Global_UndoButtonPushed, true);
            app.Segmentation_Global_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Segmentation_Global_UndoButton.FontSize = 14;
            app.Segmentation_Global_UndoButton.FontWeight = 'bold';
            app.Segmentation_Global_UndoButton.FontColor = [1 1 1];
            app.Segmentation_Global_UndoButton.Enable = 'off';
            app.Segmentation_Global_UndoButton.Position = [167 289 71 25];
            app.Segmentation_Global_UndoButton.Text = 'Undo';

            % Create Segmentation_Global_SaveButton
            app.Segmentation_Global_SaveButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_Global_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_Global_SaveButtonPushed, true);
            app.Segmentation_Global_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Segmentation_Global_SaveButton.FontSize = 14;
            app.Segmentation_Global_SaveButton.FontWeight = 'bold';
            app.Segmentation_Global_SaveButton.FontColor = [1 1 1];
            app.Segmentation_Global_SaveButton.Enable = 'off';
            app.Segmentation_Global_SaveButton.Position = [248 289 71 25];
            app.Segmentation_Global_SaveButton.Text = 'Save';

            % Create Segmentation_Text2
            app.Segmentation_Text2 = uilabel(app.SegmentationTab);
            app.Segmentation_Text2.FontSize = 14;
            app.Segmentation_Text2.FontWeight = 'bold';
            app.Segmentation_Text2.Tooltip = {'See: https://www.mathworks.com/help/images/ref/imdiffusefilt.html and  Perona, P., and J. Malik. "Scale-space and edge detection using anisotropic diffusion." IEEE® Transactions on Pattern Analysis and Machine Intelligence. Vol. 12, No. 7, July 1990, pp. 629–639.'};
            app.Segmentation_Text2.Position = [17 247 314 22];
            app.Segmentation_Text2.Text = 'Local threshold, slice per slice along 3rd Axis';

            % Create MethodDropDown_2Label
            app.MethodDropDown_2Label = uilabel(app.SegmentationTab);
            app.MethodDropDown_2Label.HorizontalAlignment = 'right';
            app.MethodDropDown_2Label.Position = [17 220 46 22];
            app.MethodDropDown_2Label.Text = 'Method';

            % Create Segmentation_Local_MethodDropDown
            app.Segmentation_Local_MethodDropDown = uidropdown(app.SegmentationTab);
            app.Segmentation_Local_MethodDropDown.Items = {'Local Otsu', 'Local Otsu translated to match expected volume fractions'};
            app.Segmentation_Local_MethodDropDown.Position = [78 220 241 22];
            app.Segmentation_Local_MethodDropDown.Value = 'Local Otsu';

            % Create Segmentation_Local_DoButton
            app.Segmentation_Local_DoButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_Local_DoButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_Local_DoButtonPushed, true);
            app.Segmentation_Local_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.Segmentation_Local_DoButton.FontSize = 14;
            app.Segmentation_Local_DoButton.FontWeight = 'bold';
            app.Segmentation_Local_DoButton.FontColor = [1 1 1];
            app.Segmentation_Local_DoButton.Position = [86 136 71 25];
            app.Segmentation_Local_DoButton.Text = 'Do';

            % Create Segmentation_Local_UndoButton
            app.Segmentation_Local_UndoButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_Local_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_Local_UndoButtonPushed, true);
            app.Segmentation_Local_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.Segmentation_Local_UndoButton.FontSize = 14;
            app.Segmentation_Local_UndoButton.FontWeight = 'bold';
            app.Segmentation_Local_UndoButton.FontColor = [1 1 1];
            app.Segmentation_Local_UndoButton.Enable = 'off';
            app.Segmentation_Local_UndoButton.Position = [167 136 71 25];
            app.Segmentation_Local_UndoButton.Text = 'Undo';

            % Create Segmentation_Local_SaveButton
            app.Segmentation_Local_SaveButton = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_Local_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_Local_SaveButtonPushed, true);
            app.Segmentation_Local_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Segmentation_Local_SaveButton.FontSize = 14;
            app.Segmentation_Local_SaveButton.FontWeight = 'bold';
            app.Segmentation_Local_SaveButton.FontColor = [1 1 1];
            app.Segmentation_Local_SaveButton.Enable = 'off';
            app.Segmentation_Local_SaveButton.Position = [248 136 71 25];
            app.Segmentation_Local_SaveButton.Text = 'Save';

            % Create SliceselectionLabel_2
            app.SliceselectionLabel_2 = uilabel(app.SegmentationTab);
            app.SliceselectionLabel_2.HorizontalAlignment = 'right';
            app.SliceselectionLabel_2.Position = [353 152 86 22];
            app.SliceselectionLabel_2.Text = 'Slice  selection';

            % Create Segmentation_SliceSlider
            app.Segmentation_SliceSlider = uislider(app.SegmentationTab);
            app.Segmentation_SliceSlider.ValueChangedFcn = createCallbackFcn(app, @Segmentation_SliceSliderValueChanged, true);
            app.Segmentation_SliceSlider.Position = [451 161 357 3];

            % Create ThresholdlowLabel
            app.ThresholdlowLabel = uilabel(app.SegmentationTab);
            app.ThresholdlowLabel.HorizontalAlignment = 'right';
            app.ThresholdlowLabel.Position = [473 104 80 22];
            app.ThresholdlowLabel.Text = 'Threshold low';

            % Create Segmentation_Lowthreshold_SliceSlider
            app.Segmentation_Lowthreshold_SliceSlider = uislider(app.SegmentationTab);
            app.Segmentation_Lowthreshold_SliceSlider.ValueChangedFcn = createCallbackFcn(app, @Segmentation_Lowthreshold_SliceSliderValueChanged, true);
            app.Segmentation_Lowthreshold_SliceSlider.Position = [565 113 243 3];

            % Create ThresholdhighLabel
            app.ThresholdhighLabel = uilabel(app.SegmentationTab);
            app.ThresholdhighLabel.HorizontalAlignment = 'right';
            app.ThresholdhighLabel.Position = [473 56 85 22];
            app.ThresholdhighLabel.Text = 'Threshold high';

            % Create Segmentation_Highthreshold_SliceSlider
            app.Segmentation_Highthreshold_SliceSlider = uislider(app.SegmentationTab);
            app.Segmentation_Highthreshold_SliceSlider.ValueChangedFcn = createCallbackFcn(app, @Segmentation_Lowthreshold_SliceSliderValueChanged, true);
            app.Segmentation_Highthreshold_SliceSlider.Position = [570 65 238 3];

            % Create Segmentation_update_range
            app.Segmentation_update_range = uibutton(app.SegmentationTab, 'push');
            app.Segmentation_update_range.ButtonPushedFcn = createCallbackFcn(app, @Segmentation_update_rangeButtonPushed, true);
            app.Segmentation_update_range.Position = [359 8 146 22];
            app.Segmentation_update_range.Text = 'Update bounds for slider';

            % Create Segmentation_Lowthreshold_SliceSpinner
            app.Segmentation_Lowthreshold_SliceSpinner = uispinner(app.SegmentationTab);
            app.Segmentation_Lowthreshold_SliceSpinner.Limits = [0 100];
            app.Segmentation_Lowthreshold_SliceSpinner.ValueChangedFcn = createCallbackFcn(app, @Segmentation_Lowthreshold_SliceSpinnerValueChanged, true);
            app.Segmentation_Lowthreshold_SliceSpinner.Position = [380 104 94 22];

            % Create Segmentation_Highthreshold_SliceSpinner
            app.Segmentation_Highthreshold_SliceSpinner = uispinner(app.SegmentationTab);
            app.Segmentation_Highthreshold_SliceSpinner.Limits = [0 100];
            app.Segmentation_Highthreshold_SliceSpinner.ValueChangedFcn = createCallbackFcn(app, @Segmentation_Highthreshold_SliceSpinnerValueChanged, true);
            app.Segmentation_Highthreshold_SliceSpinner.Position = [380 56 94 22];

            % Create OpacityoverlayEditFieldLabel
            app.OpacityoverlayEditFieldLabel = uilabel(app.SegmentationTab);
            app.OpacityoverlayEditFieldLabel.HorizontalAlignment = 'right';
            app.OpacityoverlayEditFieldLabel.Position = [515 8 96 22];
            app.OpacityoverlayEditFieldLabel.Text = 'Opacity (overlay)';

            % Create Segmentation_OpacityoverlayEditField
            app.Segmentation_OpacityoverlayEditField = uieditfield(app.SegmentationTab, 'numeric');
            app.Segmentation_OpacityoverlayEditField.Limits = [0 1];
            app.Segmentation_OpacityoverlayEditField.Position = [621 8 45 22];
            app.Segmentation_OpacityoverlayEditField.Value = 0.5;

            % Create Segmentation_SmoothSliceThreshold_CheckBox
            app.Segmentation_SmoothSliceThreshold_CheckBox = uicheckbox(app.SegmentationTab);
            app.Segmentation_SmoothSliceThreshold_CheckBox.ValueChangedFcn = createCallbackFcn(app, @Segmentation_SmoothSliceThreshold_CheckBoxValueChanged, true);
            app.Segmentation_SmoothSliceThreshold_CheckBox.Text = 'Smooth threshold along axis';
            app.Segmentation_SmoothSliceThreshold_CheckBox.Position = [17 193 173 22];
            app.Segmentation_SmoothSliceThreshold_CheckBox.Value = true;

            % Create MovingaveragefilterrangeEditFieldLabel
            app.MovingaveragefilterrangeEditFieldLabel = uilabel(app.SegmentationTab);
            app.MovingaveragefilterrangeEditFieldLabel.HorizontalAlignment = 'right';
            app.MovingaveragefilterrangeEditFieldLabel.Position = [17 166 151 22];
            app.MovingaveragefilterrangeEditFieldLabel.Text = 'Moving average filter range';

            % Create Segmentation_MovingaveragefilterrangeEditField
            app.Segmentation_MovingaveragefilterrangeEditField = uieditfield(app.SegmentationTab, 'numeric');
            app.Segmentation_MovingaveragefilterrangeEditField.Limits = [1 Inf];
            app.Segmentation_MovingaveragefilterrangeEditField.RoundFractionalValues = 'on';
            app.Segmentation_MovingaveragefilterrangeEditField.Position = [183 166 52 22];
            app.Segmentation_MovingaveragefilterrangeEditField.Value = 10;

            % Create SegmentationsensitivityTab
            app.SegmentationsensitivityTab = uitab(app.TabGroup);
            app.SegmentationsensitivityTab.Title = 'Segmentation (sensitivity)';

            % Create Sensitivity_PorosityCheckBox
            app.Sensitivity_PorosityCheckBox = uicheckbox(app.SegmentationsensitivityTab);
            app.Sensitivity_PorosityCheckBox.Enable = 'off';
            app.Sensitivity_PorosityCheckBox.Tooltip = {'Volume fraction of phase with low grey level'};
            app.Sensitivity_PorosityCheckBox.Text = 'Porosity';
            app.Sensitivity_PorosityCheckBox.Position = [17 502 65 22];
            app.Sensitivity_PorosityCheckBox.Value = true;

            % Create Sensitivity_ParticlediameterCheckBox
            app.Sensitivity_ParticlediameterCheckBox = uicheckbox(app.SegmentationsensitivityTab);
            app.Sensitivity_ParticlediameterCheckBox.Tooltip = {'Particle diameter calculatef for phase with high grey level, using a spherical assumption (continuum particle size distribution, C-PSD).'};
            app.Sensitivity_ParticlediameterCheckBox.Text = 'Particle diameter';
            app.Sensitivity_ParticlediameterCheckBox.Position = [17 475 112 22];

            % Create Sensitivity_SpecificsurfaceareaCheckBox
            app.Sensitivity_SpecificsurfaceareaCheckBox = uicheckbox(app.SegmentationsensitivityTab);
            app.Sensitivity_SpecificsurfaceareaCheckBox.Tooltip = {'Surface of phase with high grey level / volume domain * 2/3. Corrective factor 2/3 to correct the surface overestimation due to cuboid representation.'};
            app.Sensitivity_SpecificsurfaceareaCheckBox.Text = 'Specific surface area';
            app.Sensitivity_SpecificsurfaceareaCheckBox.Position = [17 448 134 22];

            % Create Sensitivity_PoretortuosityCheckBox
            app.Sensitivity_PoretortuosityCheckBox = uicheckbox(app.SegmentationsensitivityTab);
            app.Sensitivity_PoretortuosityCheckBox.Tooltip = {'Tortuosity factor of phase with low grey level'};
            app.Sensitivity_PoretortuosityCheckBox.Text = 'Pore tortuosity';
            app.Sensitivity_PoretortuosityCheckBox.Position = [17 421 99 22];

            % Create LowerboundEditFieldLabel
            app.LowerboundEditFieldLabel = uilabel(app.SegmentationsensitivityTab);
            app.LowerboundEditFieldLabel.HorizontalAlignment = 'right';
            app.LowerboundEditFieldLabel.Position = [17 630 75 22];
            app.LowerboundEditFieldLabel.Text = 'Lower bound';

            % Create Sensitivity_LowerboundEditField
            app.Sensitivity_LowerboundEditField = uieditfield(app.SegmentationsensitivityTab, 'numeric');
            app.Sensitivity_LowerboundEditField.Limits = [0 1];
            app.Sensitivity_LowerboundEditField.ValueChangedFcn = createCallbackFcn(app, @Sensitivity_LowerboundEditFieldValueChanged, true);
            app.Sensitivity_LowerboundEditField.Position = [133 630 60 22];
            app.Sensitivity_LowerboundEditField.Value = 0.4;

            % Create HigherboundEditFieldLabel
            app.HigherboundEditFieldLabel = uilabel(app.SegmentationsensitivityTab);
            app.HigherboundEditFieldLabel.HorizontalAlignment = 'right';
            app.HigherboundEditFieldLabel.Position = [17 603 78 22];
            app.HigherboundEditFieldLabel.Text = 'Higher bound';

            % Create Sensitivity_HigherboundEditField
            app.Sensitivity_HigherboundEditField = uieditfield(app.SegmentationsensitivityTab, 'numeric');
            app.Sensitivity_HigherboundEditField.Limits = [0 1];
            app.Sensitivity_HigherboundEditField.ValueChangedFcn = createCallbackFcn(app, @Sensitivity_LowerboundEditFieldValueChanged, true);
            app.Sensitivity_HigherboundEditField.Position = [133 603 60 22];
            app.Sensitivity_HigherboundEditField.Value = 0.6;

            % Create ExpectedporosityEditFieldLabel
            app.ExpectedporosityEditFieldLabel = uilabel(app.SegmentationsensitivityTab);
            app.ExpectedporosityEditFieldLabel.HorizontalAlignment = 'right';
            app.ExpectedporosityEditFieldLabel.Position = [17 657 101 22];
            app.ExpectedporosityEditFieldLabel.Text = 'Expected porosity';

            % Create Sensitivity_ExpectedporosityEditField
            app.Sensitivity_ExpectedporosityEditField = uieditfield(app.SegmentationsensitivityTab, 'numeric');
            app.Sensitivity_ExpectedporosityEditField.Limits = [0 1];
            app.Sensitivity_ExpectedporosityEditField.ValueChangedFcn = createCallbackFcn(app, @Sensitivity_LowerboundEditFieldValueChanged, true);
            app.Sensitivity_ExpectedporosityEditField.Position = [133 657 60 22];
            app.Sensitivity_ExpectedporosityEditField.Value = 0.5;

            % Create MaximumnumberofthresholdevaluatedEditFieldLabel
            app.MaximumnumberofthresholdevaluatedEditFieldLabel = uilabel(app.SegmentationsensitivityTab);
            app.MaximumnumberofthresholdevaluatedEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumnumberofthresholdevaluatedEditFieldLabel.Position = [17 576 225 22];
            app.MaximumnumberofthresholdevaluatedEditFieldLabel.Text = 'Maximum number of threshold evaluated';

            % Create Sensitivity_MaximumnumberofthresholdevaluatedEditField
            app.Sensitivity_MaximumnumberofthresholdevaluatedEditField = uieditfield(app.SegmentationsensitivityTab, 'numeric');
            app.Sensitivity_MaximumnumberofthresholdevaluatedEditField.Limits = [0 Inf];
            app.Sensitivity_MaximumnumberofthresholdevaluatedEditField.Position = [257 576 60 22];
            app.Sensitivity_MaximumnumberofthresholdevaluatedEditField.Value = 40;

            % Create VoxelsizemicrometersLabel
            app.VoxelsizemicrometersLabel = uilabel(app.SegmentationsensitivityTab);
            app.VoxelsizemicrometersLabel.HorizontalAlignment = 'right';
            app.VoxelsizemicrometersLabel.Position = [168 475 137 22];
            app.VoxelsizemicrometersLabel.Text = 'Voxel size (micrometers)';

            % Create Sensitivity_VoxelsizenanometersEditField
            app.Sensitivity_VoxelsizenanometersEditField = uieditfield(app.SegmentationsensitivityTab, 'numeric');
            app.Sensitivity_VoxelsizenanometersEditField.Limits = [0 Inf];
            app.Sensitivity_VoxelsizenanometersEditField.Position = [320 475 60 22];
            app.Sensitivity_VoxelsizenanometersEditField.Value = 0.5;

            % Create Sensitivity_instructions
            app.Sensitivity_instructions = uilabel(app.SegmentationsensitivityTab);
            app.Sensitivity_instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.Sensitivity_instructions.HorizontalAlignment = 'center';
            app.Sensitivity_instructions.FontWeight = 'bold';
            app.Sensitivity_instructions.Position = [17 720 802 44];
            app.Sensitivity_instructions.Text = {'Instructions: ONLY for non-segmented volumes with two phases. 1) Choose the porosity range you want to investigate. 2) Select'; 'the parameters to calculate for each threshold. 3) Click on ''Calculate parameter sensitivity with threshold'' button.'};

            % Create Sensitivity_Text_1
            app.Sensitivity_Text_1 = uilabel(app.SegmentationsensitivityTab);
            app.Sensitivity_Text_1.FontSize = 14;
            app.Sensitivity_Text_1.FontWeight = 'bold';
            app.Sensitivity_Text_1.Position = [17 684 104 22];
            app.Sensitivity_Text_1.Text = 'Porosity range';

            % Create Sensitivity_Text_2
            app.Sensitivity_Text_2 = uilabel(app.SegmentationsensitivityTab);
            app.Sensitivity_Text_2.FontSize = 14;
            app.Sensitivity_Text_2.FontWeight = 'bold';
            app.Sensitivity_Text_2.Position = [17 529 267 22];
            app.Sensitivity_Text_2.Text = 'Microstructure parameters to calculate';

            % Create Sensitivity_CalculateparametersensitivityButton
            app.Sensitivity_CalculateparametersensitivityButton = uibutton(app.SegmentationsensitivityTab, 'push');
            app.Sensitivity_CalculateparametersensitivityButton.ButtonPushedFcn = createCallbackFcn(app, @Sensitivity_CalculateparametersensitivityButtonPushed, true);
            app.Sensitivity_CalculateparametersensitivityButton.Position = [17 353 259 22];
            app.Sensitivity_CalculateparametersensitivityButton.Text = 'Calculate parameter sensitivity with threshold';

            % Create Sensitivity_ProgressnotrunningLabel
            app.Sensitivity_ProgressnotrunningLabel = uilabel(app.SegmentationsensitivityTab);
            app.Sensitivity_ProgressnotrunningLabel.FontAngle = 'italic';
            app.Sensitivity_ProgressnotrunningLabel.Position = [17 326 363 22];
            app.Sensitivity_ProgressnotrunningLabel.Text = 'Progress: not running';

            % Create PhasereassignmentTab
            app.PhasereassignmentTab = uitab(app.TabGroup);
            app.PhasereassignmentTab.Title = 'Phase re-assignment';

            % Create PhaseAssignment_Instructions
            app.PhaseAssignment_Instructions = uilabel(app.PhasereassignmentTab);
            app.PhaseAssignment_Instructions.BackgroundColor = [0.4706 0.6706 0.1882];
            app.PhaseAssignment_Instructions.HorizontalAlignment = 'center';
            app.PhaseAssignment_Instructions.FontWeight = 'bold';
            app.PhaseAssignment_Instructions.Position = [17 720 802 44];
            app.PhaseAssignment_Instructions.Text = {'Instructions: Modify ''New id'' column. To apply, click on the ''Do'' button.'; 'To cancel, click on the ''Undo'' button. To keep changes and update the history log, click on the ''Save'' button.'};

            % Create PhaseAssignment_UITable
            app.PhaseAssignment_UITable = uitable(app.PhasereassignmentTab);
            app.PhaseAssignment_UITable.ColumnName = {'Id'; 'Volume fraction'; 'New Id'};
            app.PhaseAssignment_UITable.RowName = {};
            app.PhaseAssignment_UITable.ColumnEditable = [false false true];
            app.PhaseAssignment_UITable.CellEditCallback = createCallbackFcn(app, @PhaseAssignment_UITableCellEdit, true);
            app.PhaseAssignment_UITable.Position = [17 349 302 353];

            % Create PhaseAssignment_DoButton
            app.PhaseAssignment_DoButton = uibutton(app.PhasereassignmentTab, 'push');
            app.PhaseAssignment_DoButton.ButtonPushedFcn = createCallbackFcn(app, @PhaseAssignment_DoButtonPushed, true);
            app.PhaseAssignment_DoButton.BackgroundColor = [0.0745 0.6235 1];
            app.PhaseAssignment_DoButton.FontSize = 14;
            app.PhaseAssignment_DoButton.FontWeight = 'bold';
            app.PhaseAssignment_DoButton.FontColor = [1 1 1];
            app.PhaseAssignment_DoButton.Position = [343 677 71 25];
            app.PhaseAssignment_DoButton.Text = 'Do';

            % Create PhaseAssignment_UndoButton
            app.PhaseAssignment_UndoButton = uibutton(app.PhasereassignmentTab, 'push');
            app.PhaseAssignment_UndoButton.ButtonPushedFcn = createCallbackFcn(app, @PhaseAssignment_UndoButtonPushed, true);
            app.PhaseAssignment_UndoButton.BackgroundColor = [0.9294 0.6941 0.1255];
            app.PhaseAssignment_UndoButton.FontSize = 14;
            app.PhaseAssignment_UndoButton.FontWeight = 'bold';
            app.PhaseAssignment_UndoButton.FontColor = [1 1 1];
            app.PhaseAssignment_UndoButton.Enable = 'off';
            app.PhaseAssignment_UndoButton.Position = [424 677 71 25];
            app.PhaseAssignment_UndoButton.Text = 'Undo';

            % Create PhaseAssignment_SaveButton
            app.PhaseAssignment_SaveButton = uibutton(app.PhasereassignmentTab, 'push');
            app.PhaseAssignment_SaveButton.ButtonPushedFcn = createCallbackFcn(app, @PhaseAssignment_SaveButtonPushed, true);
            app.PhaseAssignment_SaveButton.BackgroundColor = [0.4667 0.6745 0.1882];
            app.PhaseAssignment_SaveButton.FontSize = 14;
            app.PhaseAssignment_SaveButton.FontWeight = 'bold';
            app.PhaseAssignment_SaveButton.FontColor = [1 1 1];
            app.PhaseAssignment_SaveButton.Enable = 'off';
            app.PhaseAssignment_SaveButton.Position = [505 677 71 25];
            app.PhaseAssignment_SaveButton.Text = 'Save';

            % Create PhaseAssignment_UpdatetableButton
            app.PhaseAssignment_UpdatetableButton = uibutton(app.PhasereassignmentTab, 'push');
            app.PhaseAssignment_UpdatetableButton.ButtonPushedFcn = createCallbackFcn(app, @PhaseAssignment_UpdatetableButtonPushed, true);
            app.PhaseAssignment_UpdatetableButton.Position = [17 322 302 22];
            app.PhaseAssignment_UpdatetableButton.Text = 'Update table';

            % Create HistoryTab
            app.HistoryTab = uitab(app.TabGroup);
            app.HistoryTab.Title = 'History';

            % Create Instructions_History
            app.Instructions_History = uilabel(app.HistoryTab);
            app.Instructions_History.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Instructions_History.HorizontalAlignment = 'center';
            app.Instructions_History.FontWeight = 'bold';
            app.Instructions_History.Position = [17 742 802 22];
            app.Instructions_History.Text = 'Each operations performed on the volume are recorded in the table below.';

            % Create UITable_History
            app.UITable_History = uitable(app.HistoryTab);
            app.UITable_History.ColumnName = {'Step'; 'Operations'; 'Parameters'; 'Elapsed time'};
            app.UITable_History.RowName = {};
            app.UITable_History.Position = [17 14 802 709];

            % Show the figure after all components are created
            app.ROIfilteringandsegmentationmoduleUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Segmentation_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ROIfilteringandsegmentationmoduleUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ROIfilteringandsegmentationmoduleUIFigure)
        end
    end
end