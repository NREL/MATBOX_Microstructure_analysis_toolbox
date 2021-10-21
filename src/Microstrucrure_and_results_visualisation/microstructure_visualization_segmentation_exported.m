classdef microstructure_visualization_segmentation_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GreylevelsegmentedimagecomparisonUIFigure  matlab.ui.Figure
        SaveMenu                        matlab.ui.container.Menu
        SelectsavefolderMenu            matlab.ui.container.Menu
        SavefigureMenu                  matlab.ui.container.Menu
        SavevideoMenu                   matlab.ui.container.Menu
        VideoformatDropDown             matlab.ui.control.DropDown
        VideoformatDropDownLabel        matlab.ui.control.Label
        Phase7volumefractionxLabel      matlab.ui.control.Label
        colphase7                       matlab.ui.control.Label
        Phase6volumefractionxLabel      matlab.ui.control.Label
        colphase6                       matlab.ui.control.Label
        Phase5volumefractionxLabel      matlab.ui.control.Label
        colphase5                       matlab.ui.control.Label
        Phase4volumefractionxLabel      matlab.ui.control.Label
        colphase4                       matlab.ui.control.Label
        Phase3volumefractionxLabel      matlab.ui.control.Label
        colphase3                       matlab.ui.control.Label
        Phase2volumefractionxLabel      matlab.ui.control.Label
        colphase2                       matlab.ui.control.Label
        Phase1volumefractionxLabel      matlab.ui.control.Label
        colphase1                       matlab.ui.control.Label
        line2                           matlab.ui.control.Label
        PhasesLabel                     matlab.ui.control.Label
        Domainssizeaxbxcum3voxelsizeumLabel  matlab.ui.control.Label
        SlicexyabLabel_2                matlab.ui.control.Label
        ViewnornaltoaxexLabel           matlab.ui.control.Label
        line1                           matlab.ui.control.Label
        FieldofviewLabel                matlab.ui.control.Label
        Button_minus_2                  matlab.ui.control.Button
        Button_plus_2                   matlab.ui.control.Button
        Button_minus                    matlab.ui.control.Button
        Button_plus                     matlab.ui.control.Button
        SaveoptionsLabel                matlab.ui.control.Label
        FramepersEditField              matlab.ui.control.NumericEditField
        FramepersEditFieldLabel         matlab.ui.control.Label
        FilenameEditField               matlab.ui.control.EditField
        FilenameLabel                   matlab.ui.control.Label
        ComparisonandcoloroptionsLabel  matlab.ui.control.Label
        SlicexyabLabel                  matlab.ui.control.Label
        SelectdirectionandsliceLabel    matlab.ui.control.Label
        UITable                         matlab.ui.control.Table
        MethodDropDown                  matlab.ui.control.DropDown
        MethodDropDownLabel             matlab.ui.control.Label
        Slider                          matlab.ui.control.Slider
        SliderLabel                     matlab.ui.control.Label
        ViewnormaltoDropDown            matlab.ui.control.DropDown
        ViewnormaltoDropDownLabel       matlab.ui.control.Label
        GridnumberEditField             matlab.ui.control.NumericEditField
        GridnumberEditFieldLabel        matlab.ui.control.Label
        OverlayEditField                matlab.ui.control.NumericEditField
        OverlayEditFieldLabel           matlab.ui.control.Label
        UIAxes_comparison               matlab.ui.control.UIAxes
        UIAxes_seg                      matlab.ui.control.UIAxes
        UIAxes_grey                     matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        greyvolume = [];
        segmentedvolume = [];
        greyname = [];
        segmentedname = [];
        voxel_size = [];
        voxel_unit = [];
        slice_segmented = [];
        slice_grey = [];
        domain_size = []; dimension = [];
        number_phase= [];
        phases_id = [];
        position_slice = [];
        direction = [];
        color_phase = [];
        RGB_phase = [];
        max_grey = [];
        x1_start = []; x1_end=[];
        x2_start = []; x2_end=[];
        savefolder = [];
    end
    
    methods (Access = private)
        
        function [] = func_update_segmented_figure(app)
            if app.direction==1
                % Initializaion
                app.slice_segmented = zeros(app.domain_size(2),app.domain_size(3),3); % RGB color map
                slice_r = zeros(app.domain_size(2),app.domain_size(3)); % Red color map
                slice_g = zeros(app.domain_size(2),app.domain_size(3)); % Green color map
                slice_b = zeros(app.domain_size(2),app.domain_size(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:app.number_phase
                    id = app.phases_id(current_phase); % Current phase code
                    slice_r(app.segmentedvolume(app.position_slice,:,:)==id)= app.RGB_phase(current_phase,1); % Attribute red color
                    slice_g(app.segmentedvolume(app.position_slice,:,:)==id)= app.RGB_phase(current_phase,2); % Attribute green color
                    slice_b(app.segmentedvolume(app.position_slice,:,:)==id)= app.RGB_phase(current_phase,3); % Attribute blue color
                end
            elseif app.direction==2
                % Initializaion
                app.slice_segmented = zeros(app.domain_size(1),app.domain_size(3),3); % RGB color map
                slice_r = zeros(app.domain_size(1),app.domain_size(3)); % Red color map
                slice_g = zeros(app.domain_size(1),app.domain_size(3)); % Green color map
                slice_b = zeros(app.domain_size(1),app.domain_size(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:app.number_phase
                    id = app.phases_id(current_phase); % Current phase code
                    slice_r(app.segmentedvolume(:,app.position_slice,:)==id)= app.RGB_phase(current_phase,1); % Attribute red color
                    slice_g(app.segmentedvolume(:,app.position_slice,:)==id)= app.RGB_phase(current_phase,2); % Attribute green color
                    slice_b(app.segmentedvolume(:,app.position_slice,:)==id)= app.RGB_phase(current_phase,3); % Attribute blue color
                end
            elseif app.direction==3
                % Initializaion
                app.slice_segmented = zeros(app.domain_size(1),app.domain_size(2),3); % RGB color map
                slice_r = zeros(app.domain_size(1),app.domain_size(2)); % Red color map
                slice_g = zeros(app.domain_size(1),app.domain_size(2)); % Green color map
                slice_b = zeros(app.domain_size(1),app.domain_size(2)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:app.number_phase
                    id =app.phases_id(current_phase); % Current phase code
                    slice_r(app.segmentedvolume(:,:,app.position_slice)==id)= app.RGB_phase(current_phase,1); % Attribute red color
                    slice_g(app.segmentedvolume(:,:,app.position_slice)==id)= app.RGB_phase(current_phase,2); % Attribute green color
                    slice_b(app.segmentedvolume(:,:,app.position_slice)==id)= app.RGB_phase(current_phase,3); % Attribute blue color
                end
            end
            app.slice_segmented(:,:,1)=slice_r; % Attribute RGB color
            app.slice_segmented(:,:,2)=slice_g;
            app.slice_segmented(:,:,3)=slice_b;
            % Display the slice
            image(app.slice_segmented,'parent',app.UIAxes_seg);
            [t,s] = title(app.UIAxes_seg,'Segmented volume',app.segmentedname,'Color','k');
            t.FontSize = 16; t.FontName = 'Helvetica'; s.FontAngle = 'italic'; s.FontSize = 14;
            % Get axis position
            axe_position =app.UIAxes_seg.Position;
            set(app.UIAxes_seg ,'Position', axe_position);
            % Remove tick and label
            set(app.UIAxes_seg,'xtick',[],'ytick',[]);
            set(app.UIAxes_seg,'xticklabel',[],'yticklabel',[]);
            set(app.UIAxes_seg,'Xlabel',[],'Ylabel',[]);
            box(app.UIAxes_seg,'on');
            % Fit the axes box
            axis(app.UIAxes_seg,'tight');
            % Aspect ratio is 1:1
            axis(app.UIAxes_seg,'equal');
            app.UIAxes_seg.YDir='normal';
        end
        
        function [] = func_update_grey_figure(app)
            if app.direction==1
                app.slice_grey = squeeze(app.greyvolume(app.position_slice,:,:));
            elseif app.direction==2
                app.slice_grey = squeeze(app.greyvolume(:,app.position_slice,:));
            elseif app.direction==3
                app.slice_grey = app.greyvolume(:,:,app.position_slice);
            end
            [t,s] = title(app.UIAxes_grey,'Grey level volume',app.greyname,'Color','k');
            t.FontSize = 16; t.FontName = 'Helvetica'; s.FontAngle = 'italic'; s.FontSize = 14;
            image(app.slice_grey(:,:,1),'CDataMapping','scaled','Parent', app.UIAxes_grey);
            colormap(app.UIAxes_grey,gray)
            % Get axis position
            axe_position =app.UIAxes_grey.Position;
            set(app.UIAxes_grey ,'Position', axe_position);
            % Remove tick and label
            set(app.UIAxes_grey,'xtick',[],'ytick',[]);
            set(app.UIAxes_grey,'xticklabel',[],'yticklabel',[]);
            set(app.UIAxes_grey,'Xlabel',[],'Ylabel',[]);
            box(app.UIAxes_grey,'on');
            % Fit the axes box
            axis(app.UIAxes_grey,'tight');
            % Aspect ratio is 1:1
            axis(app.UIAxes_grey,'equal');
            app.UIAxes_grey.YDir='normal';
        end
        
        function [] = func_update_overlay(app)
            % Convert grey image in RGB
            grey_in_RGB = cat(3, app.slice_grey, app.slice_grey, app.slice_grey);
            grey_in_RGB=double(grey_in_RGB)/app.max_grey; % Scale
            % Linear combination
            C = app.OverlayEditField.Value*app.slice_segmented + (1-app.OverlayEditField.Value)*grey_in_RGB;
            % Display
            image(C,'Parent', app.UIAxes_comparison);
            [t] = title(app.UIAxes_comparison,'Overlay','Color','k');
            t.FontSize = 16; t.FontName = 'Helvetica';
            % Get axis position
            axe_position =app.UIAxes_comparison.Position;
            set(app.UIAxes_comparison ,'Position', axe_position);
            % Remove tick and label
            set(app.UIAxes_comparison,'xtick',[],'ytick',[]);
            set(app.UIAxes_comparison,'xticklabel',[],'yticklabel',[]);
            set(app.UIAxes_comparison,'Xlabel',[],'Ylabel',[]);
            box(app.UIAxes_comparison,'on');
            % Fit the axes box
            axis(app.UIAxes_comparison,'tight');
            % Aspect ratio is 1:1
            axis(app.UIAxes_comparison,'equal');
            app.UIAxes_comparison.YDir='normal';
        end
        
        function [] = func_update_checkerboard(app)
            % Convert grey image in RGB
            grey_in_RGB = cat(3, app.slice_grey, app.slice_grey, app.slice_grey);
            grey_in_RGB=double(grey_in_RGB)/app.max_grey; % Scale
            % Create grid
            checkerboard_slice = app.slice_segmented;
            check=1;
            for k1=1:1:length(app.x1_start)
                check=abs(check-1);
                for k2=1+check:2:length(app.x2_start)
                    checkerboard_slice(app.x1_start(k1):app.x1_end(k1),app.x2_start(k2):app.x2_end(k2),1) = grey_in_RGB(app.x1_start(k1):app.x1_end(k1),app.x2_start(k2):app.x2_end(k2),1);
                    checkerboard_slice(app.x1_start(k1):app.x1_end(k1),app.x2_start(k2):app.x2_end(k2),2) = grey_in_RGB(app.x1_start(k1):app.x1_end(k1),app.x2_start(k2):app.x2_end(k2),2);
                    checkerboard_slice(app.x1_start(k1):app.x1_end(k1),app.x2_start(k2):app.x2_end(k2),3) = grey_in_RGB(app.x1_start(k1):app.x1_end(k1),app.x2_start(k2):app.x2_end(k2),3);
                end
            end
            % Display
            image(checkerboard_slice,'Parent', app.UIAxes_comparison);
            [t] = title(app.UIAxes_comparison,'Checkerboard','Color','k');
            t.FontSize = 16; t.FontName = 'Helvetica';
            % Get axis position
            axe_position =app.UIAxes_comparison.Position;
            set(app.UIAxes_comparison ,'Position', axe_position);
            % Remove tick and label
            set(app.UIAxes_comparison,'xtick',[],'ytick',[]);
            set(app.UIAxes_comparison,'xticklabel',[],'yticklabel',[]);
            set(app.UIAxes_comparison,'Xlabel',[],'Ylabel',[]);
            box(app.UIAxes_comparison,'on');
            % Fit the axes box
            axis(app.UIAxes_comparison,'tight');
            % Aspect ratio is 1:1
            axis(app.UIAxes_comparison,'equal');
            app.UIAxes_comparison.YDir='normal';
        end
        
        function [] = func_update_allaxes(app)
            app.func_update_segmented_figure;
            app.func_update_grey_figure;
            if strcmp(app.MethodDropDown.Value,'Overlay')
                app.func_update_overlay;
            elseif strcmp(app.MethodDropDown.Value,'Checkerboard')
                app.func_update_checkerboard;
            end
            app.SlicexyabLabel.Text = ['Slice: ' num2str(app.position_slice) ' / ' num2str(app.domain_size(app.direction)) ', ' num2str(app.position_slice*app.voxel_size,'%1.1f') ' / ' num2str(app.domain_size(app.direction)*app.voxel_size,'%1.1f') app.voxel_unit];
        end
        
        function [] = editmode(app,bool)
            if bool
                display_edit = 'on';
                display_save = 'off';
            else
                display_edit = 'off';
                display_save = 'on';
            end
            
            app.ViewnormaltoDropDown.Visible = display_edit;
            app.Button_minus_2.Visible = display_edit;
            app.Button_minus.Visible = display_edit;
            app.Button_plus_2.Visible = display_edit;
            app.Button_plus.Visible = display_edit;
            app.Slider.Visible = display_edit;
            app.SliderLabel.Visible = display_edit;
            app.SlicexyabLabel.Visible = display_edit;
            app.ComparisonandcoloroptionsLabel.Visible = display_edit;
            app.UITable.Visible = display_edit;
            app.MethodDropDown.Visible = display_edit;
            app.MethodDropDownLabel.Visible = display_edit;
            app.GridnumberEditField.Visible = display_edit;
            app.OverlayEditField.Visible = display_edit;
            app.GridnumberEditFieldLabel.Visible = display_edit;
            app.OverlayEditFieldLabel.Visible = display_edit;
            app.SelectdirectionandsliceLabel.Visible = display_edit;
            app.ViewnormaltoDropDownLabel.Visible = display_edit;
            app.SaveoptionsLabel.Visible = display_edit;
            app.SaveoptionsLabel.Visible = display_edit;
            app.SaveoptionsLabel.Visible = display_edit;
            app.FilenameEditField.Visible = display_edit;
            app.FilenameLabel.Visible = display_edit;
            app.FramepersEditField.Visible = display_edit;
            app.FramepersEditFieldLabel.Visible = display_edit;
            app.VideoformatDropDown.Visible = display_edit;
            app.VideoformatDropDownLabel.Visible = display_edit;
            
            if bool
                app.Phase1volumefractionxLabel.Visible = display_save;
                app.Phase2volumefractionxLabel.Visible = display_save;
                app.Phase3volumefractionxLabel.Visible = display_save;
                app.Phase4volumefractionxLabel.Visible = display_save;
                app.Phase5volumefractionxLabel.Visible = display_save;
                app.Phase6volumefractionxLabel.Visible = display_save;
                app.Phase7volumefractionxLabel.Visible = display_save;
                app.colphase1.Visible = display_save;
                app.colphase2.Visible = display_save;
                app.colphase3.Visible = display_save;
                app.colphase4.Visible = display_save;
                app.colphase5.Visible = display_save;
                app.colphase6.Visible = display_save;
                app.colphase7.Visible = display_save;
            end
            
            app.FieldofviewLabel.Visible = display_save;
            app.PhasesLabel.Visible = display_save;
            app.Domainssizeaxbxcum3voxelsizeumLabel.Visible = display_save;
            app.ViewnornaltoaxexLabel.Visible = display_save;
            app.SlicexyabLabel_2.Visible = display_save;
            app.line1.Visible = display_save;
            app.line2.Visible = display_save;
            
            if bool
                app.MethodDropDownValueChanged;
            else
                app.ViewnornaltoaxexLabel.Text = ['View nornal to axe ' num2str(app.direction)];
                app.SlicexyabLabel_2.Text = app.SlicexyabLabel.Text;
                s=app.domain_size*app.voxel_size;
                app.Domainssizeaxbxcum3voxelsizeumLabel.Text = ['Domain''s size: ' num2str(s(1),'%1.1f') ' x ' num2str(s(2),'%1.1f') ' x ' num2str(s(3),'%1.1f') ' ' app.voxel_unit char(179) ', voxel size: ' num2str(app.voxel_size,'%1.3f')  app.voxel_unit];
                datatable = app.UITable.Data;
                name = datatable(:,6);
                vf = datatable(:,2);
                rgb = datatable(:,3:5);
                for phase=1:1:app.number_phase % tedious...
                    if phase==1
                        app.Phase1volumefractionxLabel.Visible = display_save;
                        app.colphase1.Visible = display_save;
                        app.Phase1volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase1.BackgroundColor = cell2mat(rgb(phase,:));
                    elseif phase==2
                        app.Phase2volumefractionxLabel.Visible = display_save;
                        app.colphase2.Visible = display_save;
                        app.Phase2volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase2.BackgroundColor = cell2mat(rgb(phase,:));
                    elseif phase==3
                        app.Phase3volumefractionxLabel.Visible = display_save;
                        app.colphase3.Visible = display_save;
                        app.Phase3volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase3.BackgroundColor = cell2mat(rgb(phase,:));
                    elseif phase==4
                        app.Phase4volumefractionxLabel.Visible = display_save;
                        app.colphase4.Visible = display_save;
                        app.Phase4volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase4.BackgroundColor = cell2mat(rgb(phase,:));
                    elseif phase==5
                        app.Phase5volumefractionxLabel.Visible = display_save;
                        app.colphase5.Visible = display_save;
                        app.Phase5volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase5.BackgroundColor = cell2mat(rgb(phase,:));
                    elseif phase==6
                        app.Phase6volumefractionxLabel.Visible = display_save;
                        app.colphase6.Visible = display_save;
                        app.Phase6volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase6.BackgroundColor = cell2mat(rgb(phase,:));
                    elseif phase==7
                        app.Phase7volumefractionxLabel.Visible = display_save;
                        app.colphase7.Visible = display_save;
                        app.Phase7volumefractionxLabel.Text = [char(name(phase)) ', volume fraction ' num2str(cell2mat(vf(phase)),'%1.3f')];
                        app.colphase7.BackgroundColor = cell2mat(rgb(phase,:));
                    end
                end
            end
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, greyvolume, segmentedvolume, greyname, segmentedname, voxel_size, voxel_unit)
            app.greyvolume = greyvolume;
            app.segmentedvolume = segmentedvolume;
            app.greyname = greyname;
            app.segmentedname = segmentedname;
            app.voxel_size = voxel_size;
            app.voxel_unit = voxel_unit;
            app.domain_size = size(app.greyvolume);
            app.dimension = length(app.domain_size);
            app.phases_id = unique(app.segmentedvolume);
            app.number_phase = length(app.phases_id);
            app.direction = 3;
            app.color_phase=[colororder;rand(1000,3)];
            app.max_grey = double(max(max(max(app.greyvolume))));
            
            if strcmp(app.voxel_unit,'um') || strcmp(app.voxel_unit,'micrometer') || strcmp(app.voxel_unit,'micrometers')
                app.voxel_unit = [char(956) 'm'];
            end
            
            if app.dimension==2
                app.ViewnormaltoDropDown.Enable = 'off';
                app.Slider.Visible = 'off'; app.Slider.Enable = 'off';
                app.Button_plus.Visible = 'off'; app.Button_plus.Enable = 'off';
                app.Button_minus.Visible = 'off'; app.Button_minus.Enable = 'off';
                app.SlicexyabLabel.Visible = 'off'; app.SlicexyabLabel.Enable = 'off';
                app.FramepersEditField.Visible = 'off'; app.FramepersEditField.Enable = 'off';
                app.FramepersEditFieldLabel.Visible = 'off';
                app.SavevideoMenu.Enable = 'off';
                app.position_slice = 1;
            else
                app.Slider.Limits = [1 app.domain_size(app.direction)];
                app.Slider.Value = round(app.domain_size(app.direction)/2);
                app.position_slice = round(app.domain_size(app.direction)/2);
                app.SlicexyabLabel.Text = ['Slice: ' num2str(app.position_slice) ' / ' num2str(app.domain_size(app.direction)) ', ' num2str(app.position_slice*app.voxel_size,'%1.1f') ' / ' num2str(app.domain_size(app.direction)*app.voxel_size,'%1.1f') app.voxel_unit];
            end
            
            % Set default colors
            app.RGB_phase = app.color_phase(1:app.number_phase,:);
            % Table
            volumefraction=zeros(app.number_phase,1);
            for n=1:1:app.number_phase
                volumefraction(n,1)=sum(sum(sum(app.segmentedvolume==app.phases_id(n))))/numel(app.segmentedvolume);
                Phase(n).name=['Phase ' num2str(app.phases_id(n))];
            end
            table_data = [num2cell(app.phases_id) num2cell(volumefraction) num2cell(app.RGB_phase(:,1)) num2cell(app.RGB_phase(:,2)) num2cell(app.RGB_phase(:,3)) {Phase.name}'];
            set(app.UITable,'Data',table_data);
            
            app.editmode(true);
            app.MethodDropDownValueChanged;
            % Update axes
            app.func_update_allaxes;
        end

        % Value changed function: ViewnormaltoDropDown
        function ViewnormaltoDropDownValueChanged(app, event)
            value = app.ViewnormaltoDropDown.Value;
            if strcmp(value,'Axe 1')
                app.direction=1;
            elseif strcmp(value,'Axe 2')
                app.direction=2;
            elseif strcmp(value,'Axe 3')
                app.direction=3;
            end
            app.Slider.Limits = [1 app.domain_size(app.direction)];
            app.Slider.Value = round(app.domain_size(app.direction)/2);
            app.func_update_allaxes;
            
        end

        % Value changed function: GridnumberEditField
        function GridnumberEditFieldValueChanged(app, event)
            g1 = app.GridnumberEditField.Value;
            s = size(app.slice_segmented);
            tmp = round(linspace(1,s(1),g1+1));
            app.x1_start = tmp(1:end-1);
            app.x1_end = [tmp(2:end-1)-1 s(1)];
            app.x2_start = 1; app.x2_end = [];
            endnotreach = true;
            while endnotreach
                app.x2_end = [app.x2_end min(app.x2_start(end)+app.x1_end(1)-1,s(2))];
                if app.x2_end(end)==s(2)
                    endnotreach=false;
                else
                    app.x2_start=[app.x2_start app.x2_end(end)+1];
                end
            end
            app.func_update_allaxes;
        end

        % Value changed function: OverlayEditField
        function OverlayEditFieldValueChanged(app, event)
            app.func_update_allaxes;
        end

        % Callback function: Slider, Slider
        function SliderValueChanged(app, event)
            value = round(app.Slider.Value);
            app.Slider.Value = value;
            app.position_slice = app.Slider.Value;
            app.func_update_allaxes;
            
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            % Get back data of the table
            table_data = app.UITable.Data;
            % Check color
            Color_=table_data(:,3:5);
            Color_ = cell2mat(Color_);
            % Only number
            detect_nan = isnan(Color_);
            tmp=sum(sum(detect_nan));
            if tmp==0
                Color_hasnot_nan=1;
            else
                Color_hasnot_nan=0;
            end
            % Only between 0 and 1
            correct_color=1;
            if Color_hasnot_nan==1
                tmp1=Color_; tmp2=Color_;
                tmp1(Color_>=0)=1; tmp2(Color_<=1)=1;
                tmp=tmp1.*tmp2;
                if (length(unique(tmp))==1 && unique(tmp)==1)
                    correct_color=1;
                else
                    correct_color=0;
                end
            else
                correct_color=0;
            end
            if ~correct_color
                Color_=app.color_phase(1:app.number_phase,:);
                table_data(:,3:5) = [num2cell(Color_(:,1)) num2cell(Color_(:,2)) num2cell(Color_(:,3))];
                set(app.UITable,'Data',table_data);
            end
            app.RGB_phase = Color_; % Update color
            app.func_update_allaxes; % Update figure
        end

        % Value changed function: MethodDropDown
        function MethodDropDownValueChanged(app, event)
            if strcmp(app.MethodDropDown.Value,'Overlay')
                app.GridnumberEditField.Visible = 'off'; app.GridnumberEditField.Enable = 'off';
                app.GridnumberEditFieldLabel.Visible = 'off';
                app.OverlayEditField.Visible = 'on'; app.OverlayEditField.Enable = 'on';
                app.OverlayEditFieldLabel.Visible = 'on';
                app.func_update_overlay;
            elseif strcmp(app.MethodDropDown.Value,'Checkerboard')
                app.OverlayEditField.Visible = 'off'; app.OverlayEditField.Enable = 'off';
                app.OverlayEditFieldLabel.Visible = 'off';
                app.GridnumberEditField.Visible = 'on'; app.GridnumberEditField.Enable = 'on';
                app.GridnumberEditFieldLabel.Visible = 'on';
                app.GridnumberEditFieldValueChanged;
                app.func_update_checkerboard;
            end
        end

        % Menu selected function: SelectsavefolderMenu
        function SelectsavefolderMenuSelected(app, event)
            str_dialogbox = 'Select save folder';
            res = uigetdir(matlabroot,str_dialogbox); % Open dialog box
            if res==0
                % User clicked cancel button or closed the dialog box
                app.SavefigureMenu.Enable = 'off';
                if app.dimension==3
                    app.SavevideoMenu.Enable = 'off';
                end
            else
                if ispc
                    app.savefolder = [res '\'];
                else
                    app.savefolder = [res '/'];
                end
                app.SavefigureMenu.Enable = 'on';
                app.SavevideoMenu.Enable = 'on';
            end
        end

        % Menu selected function: SavefigureMenu
        function SavefigureMenuSelected(app, event)
            app.editmode(false);
            figure_filename = [app.FilenameEditField.Value '_' app.MethodDropDown.Value '_axe' num2str(app.direction) '_slice' num2str(app.position_slice) '.png'];
            % Can't save an app design figure yet (below won't work in 2020b)
            % function_savefig(app.GreylevelsegmentedimagecomparisonUIFigure, app.savefolder, figure_filename, OPTIONS);
            % Workaround proposed by Matlab will only capture the axis
            % exportgraphics(app.GreylevelsegmentedimagecomparisonUIFigure,[app.savefolder figure_filename])
            % Instead use screencapture
            % Yair Altman (2021). ScreenCapture - get a screen-capture of a figure frame or component
            % https://www.mathworks.com/matlabcentral/fileexchange/24323-screencapture-get-a-screen-capture-of-a-figure-frame-or-component
            % MATLAB Central File Exchange. Retrieved March 4, 2021.
            imageData = screencapture(app.GreylevelsegmentedimagecomparisonUIFigure);
            imwrite(imageData,[app.savefolder figure_filename]);
            app.editmode(true);
        end

        % Menu selected function: SavevideoMenu
        function SavevideoMenuSelected(app, event)
            app.editmode(false);
            video_filename = [app.FilenameEditField.Value '_' app.MethodDropDown.Value '_axe' num2str(app.direction)];
            video_format = app.VideoformatDropDown.Value;
            video_handle = VideoWriter([app.savefolder video_filename],video_format);
            if strcmp(video_format,'mpeg-4')
                set(video_handle,'Quality',100); % Set video quality
            end
            % Set video framerate
            set(video_handle,'FrameRate',app.FramepersEditField.Value);
            % Open video
            open(video_handle)
            for frame_number=1:1:app.domain_size(app.direction)
                % Update figure
                app.position_slice = frame_number;
                app.SlicexyabLabel_2.Text = ['Slice: ' num2str(app.position_slice) ' / ' num2str(app.domain_size(app.direction)) ', ' num2str(app.position_slice*app.voxel_size,'%1.1f') ' / ' num2str(app.domain_size(app.direction)*app.voxel_size,'%1.1f') app.voxel_unit];
                app.func_update_allaxes;
                stored_frame(frame_number) = getframe(app.GreylevelsegmentedimagecomparisonUIFigure);
                writeVideo(video_handle,stored_frame(frame_number))
                pause(0.01);
            end
            % Close video
            close(video_handle)
            app.editmode(true);
        end

        % Button pushed function: Button_minus
        function Button_minusPushed(app, event)
            app.Slider.Value = round(app.Slider.Value-1);
            app.Slider.Value = max(app.Slider.Limits(1), round(app.Slider.Value));
            app.position_slice = app.Slider.Value;
            app.func_update_allaxes;
        end

        % Button pushed function: Button_plus
        function Button_plusPushed(app, event)
            app.Slider.Value = round(app.Slider.Value+1);
            app.Slider.Value = min(app.Slider.Limits(2), round(app.Slider.Value));
            app.position_slice = app.Slider.Value;
            app.func_update_allaxes;
        end

        % Button pushed function: Button_plus_2
        function Button_plus_2Pushed(app, event)
            major_step = 1/10 * app.domain_size(app.direction);
            app.Slider.Value = min(app.Slider.Limits(2), round(app.Slider.Value+major_step));
            app.position_slice = app.Slider.Value;
            app.func_update_allaxes;
        end

        % Button pushed function: Button_minus_2
        function Button_minus_2Pushed(app, event)
            major_step = 1/10 * app.domain_size(app.direction);
            app.Slider.Value = max(app.Slider.Limits(1), round(app.Slider.Value-major_step));
            app.position_slice = app.Slider.Value;
            app.func_update_allaxes;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GreylevelsegmentedimagecomparisonUIFigure and hide until all components are created
            app.GreylevelsegmentedimagecomparisonUIFigure = uifigure('Visible', 'off');
            app.GreylevelsegmentedimagecomparisonUIFigure.Position = [100 100 1047 782];
            app.GreylevelsegmentedimagecomparisonUIFigure.Name = 'Grey level - segmented image comparison';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.SaveMenu.Text = 'Save';

            % Create SelectsavefolderMenu
            app.SelectsavefolderMenu = uimenu(app.SaveMenu);
            app.SelectsavefolderMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectsavefolderMenuSelected, true);
            app.SelectsavefolderMenu.Text = 'Select save folder';

            % Create SavefigureMenu
            app.SavefigureMenu = uimenu(app.SaveMenu);
            app.SavefigureMenu.MenuSelectedFcn = createCallbackFcn(app, @SavefigureMenuSelected, true);
            app.SavefigureMenu.Enable = 'off';
            app.SavefigureMenu.Text = 'Save figure';

            % Create SavevideoMenu
            app.SavevideoMenu = uimenu(app.SaveMenu);
            app.SavevideoMenu.MenuSelectedFcn = createCallbackFcn(app, @SavevideoMenuSelected, true);
            app.SavevideoMenu.Enable = 'off';
            app.SavevideoMenu.Text = 'Save video';

            % Create UIAxes_grey
            app.UIAxes_grey = uiaxes(app.GreylevelsegmentedimagecomparisonUIFigure);
            title(app.UIAxes_grey, 'Title')
            xlabel(app.UIAxes_grey, 'X')
            ylabel(app.UIAxes_grey, 'Y')
            zlabel(app.UIAxes_grey, 'Z')
            app.UIAxes_grey.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_grey.Position = [58 410 408 362];

            % Create UIAxes_seg
            app.UIAxes_seg = uiaxes(app.GreylevelsegmentedimagecomparisonUIFigure);
            title(app.UIAxes_seg, 'Title')
            xlabel(app.UIAxes_seg, 'X')
            ylabel(app.UIAxes_seg, 'Y')
            zlabel(app.UIAxes_seg, 'Z')
            app.UIAxes_seg.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_seg.Position = [573 410 408 362];

            % Create UIAxes_comparison
            app.UIAxes_comparison = uiaxes(app.GreylevelsegmentedimagecomparisonUIFigure);
            title(app.UIAxes_comparison, 'Title')
            xlabel(app.UIAxes_comparison, 'X')
            ylabel(app.UIAxes_comparison, 'Y')
            zlabel(app.UIAxes_comparison, 'Z')
            app.UIAxes_comparison.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_comparison.Position = [58 25 408 362];

            % Create OverlayEditFieldLabel
            app.OverlayEditFieldLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.OverlayEditFieldLabel.HorizontalAlignment = 'right';
            app.OverlayEditFieldLabel.Position = [808 198 47 22];
            app.OverlayEditFieldLabel.Text = 'Overlay';

            % Create OverlayEditField
            app.OverlayEditField = uieditfield(app.GreylevelsegmentedimagecomparisonUIFigure, 'numeric');
            app.OverlayEditField.Limits = [0 1];
            app.OverlayEditField.ValueChangedFcn = createCallbackFcn(app, @OverlayEditFieldValueChanged, true);
            app.OverlayEditField.Position = [870 198 52 22];
            app.OverlayEditField.Value = 0.2;

            % Create GridnumberEditFieldLabel
            app.GridnumberEditFieldLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.GridnumberEditFieldLabel.HorizontalAlignment = 'right';
            app.GridnumberEditFieldLabel.Position = [780 198 72 22];
            app.GridnumberEditFieldLabel.Text = 'Grid number';

            % Create GridnumberEditField
            app.GridnumberEditField = uieditfield(app.GreylevelsegmentedimagecomparisonUIFigure, 'numeric');
            app.GridnumberEditField.Limits = [2 Inf];
            app.GridnumberEditField.RoundFractionalValues = 'on';
            app.GridnumberEditField.ValueChangedFcn = createCallbackFcn(app, @GridnumberEditFieldValueChanged, true);
            app.GridnumberEditField.Position = [867 198 55 22];
            app.GridnumberEditField.Value = 8;

            % Create ViewnormaltoDropDownLabel
            app.ViewnormaltoDropDownLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.ViewnormaltoDropDownLabel.HorizontalAlignment = 'right';
            app.ViewnormaltoDropDownLabel.Position = [677 345 85 22];
            app.ViewnormaltoDropDownLabel.Text = 'View normal to';

            % Create ViewnormaltoDropDown
            app.ViewnormaltoDropDown = uidropdown(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.ViewnormaltoDropDown.Items = {'Axe 1', 'Axe 2', 'Axe 3'};
            app.ViewnormaltoDropDown.ValueChangedFcn = createCallbackFcn(app, @ViewnormaltoDropDownValueChanged, true);
            app.ViewnormaltoDropDown.Position = [777 345 100 22];
            app.ViewnormaltoDropDown.Value = 'Axe 3';

            % Create SliderLabel
            app.SliderLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.SliderLabel.HorizontalAlignment = 'right';
            app.SliderLabel.Position = [583 312 36 22];
            app.SliderLabel.Text = 'Slider';

            % Create Slider
            app.Slider = uislider(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Slider.Position = [640 321 295 3];

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.Position = [573 198 46 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create MethodDropDown
            app.MethodDropDown = uidropdown(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.MethodDropDown.Items = {'Overlay', 'Checkerboard'};
            app.MethodDropDown.ValueChangedFcn = createCallbackFcn(app, @MethodDropDownValueChanged, true);
            app.MethodDropDown.Position = [634 198 128 22];
            app.MethodDropDown.Value = 'Overlay';

            % Create UITable
            app.UITable = uitable(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.UITable.ColumnName = {'Phase'; 'Volume fraction'; 'Red'; 'Green'; 'Blue'; 'Name'};
            app.UITable.ColumnWidth = {75, 'auto', 'auto', 'auto', 'auto', '1x'};
            app.UITable.RowName = {};
            app.UITable.ColumnEditable = [false false true true true true];
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.Position = [527 74 499 112];

            % Create SelectdirectionandsliceLabel
            app.SelectdirectionandsliceLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.SelectdirectionandsliceLabel.FontSize = 14;
            app.SelectdirectionandsliceLabel.FontWeight = 'bold';
            app.SelectdirectionandsliceLabel.Position = [690 375 174 22];
            app.SelectdirectionandsliceLabel.Text = 'Select direction and slice';

            % Create SlicexyabLabel
            app.SlicexyabLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.SlicexyabLabel.HorizontalAlignment = 'center';
            app.SlicexyabLabel.Position = [583 261 388 22];
            app.SlicexyabLabel.Text = 'Slice: x / y, a / b';

            % Create ComparisonandcoloroptionsLabel
            app.ComparisonandcoloroptionsLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.ComparisonandcoloroptionsLabel.FontSize = 14;
            app.ComparisonandcoloroptionsLabel.FontWeight = 'bold';
            app.ComparisonandcoloroptionsLabel.Position = [672 230 209 22];
            app.ComparisonandcoloroptionsLabel.Text = 'Comparison and color options';

            % Create FilenameLabel
            app.FilenameLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.FilenameLabel.HorizontalAlignment = 'right';
            app.FilenameLabel.Position = [492 13 58 22];
            app.FilenameLabel.Text = 'Filename ';

            % Create FilenameEditField
            app.FilenameEditField = uieditfield(app.GreylevelsegmentedimagecomparisonUIFigure, 'text');
            app.FilenameEditField.Tooltip = {'Axis, slice, and comparison will be automatically added to the filename.'};
            app.FilenameEditField.Position = [565 13 156 22];

            % Create FramepersEditFieldLabel
            app.FramepersEditFieldLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.FramepersEditFieldLabel.HorizontalAlignment = 'right';
            app.FramepersEditFieldLabel.Position = [730 13 74 22];
            app.FramepersEditFieldLabel.Text = 'Frame per s.';

            % Create FramepersEditField
            app.FramepersEditField = uieditfield(app.GreylevelsegmentedimagecomparisonUIFigure, 'numeric');
            app.FramepersEditField.Limits = [1 Inf];
            app.FramepersEditField.RoundFractionalValues = 'on';
            app.FramepersEditField.Position = [811 13 44 22];
            app.FramepersEditField.Value = 25;

            % Create SaveoptionsLabel
            app.SaveoptionsLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.SaveoptionsLabel.FontSize = 14;
            app.SaveoptionsLabel.FontWeight = 'bold';
            app.SaveoptionsLabel.Position = [730 42 93 22];
            app.SaveoptionsLabel.Text = 'Save options';

            % Create Button_plus
            app.Button_plus = uibutton(app.GreylevelsegmentedimagecomparisonUIFigure, 'push');
            app.Button_plus.ButtonPushedFcn = createCallbackFcn(app, @Button_plusPushed, true);
            app.Button_plus.Position = [958 325 28 22];
            app.Button_plus.Text = '+';

            % Create Button_minus
            app.Button_minus = uibutton(app.GreylevelsegmentedimagecomparisonUIFigure, 'push');
            app.Button_minus.ButtonPushedFcn = createCallbackFcn(app, @Button_minusPushed, true);
            app.Button_minus.Position = [958 299 28 22];
            app.Button_minus.Text = '-';

            % Create Button_plus_2
            app.Button_plus_2 = uibutton(app.GreylevelsegmentedimagecomparisonUIFigure, 'push');
            app.Button_plus_2.ButtonPushedFcn = createCallbackFcn(app, @Button_plus_2Pushed, true);
            app.Button_plus_2.Position = [988 325 30 22];
            app.Button_plus_2.Text = '++';

            % Create Button_minus_2
            app.Button_minus_2 = uibutton(app.GreylevelsegmentedimagecomparisonUIFigure, 'push');
            app.Button_minus_2.ButtonPushedFcn = createCallbackFcn(app, @Button_minus_2Pushed, true);
            app.Button_minus_2.Position = [988 299 30 22];
            app.Button_minus_2.Text = '--';

            % Create FieldofviewLabel
            app.FieldofviewLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.FieldofviewLabel.FontSize = 14;
            app.FieldofviewLabel.FontWeight = 'bold';
            app.FieldofviewLabel.Position = [543 355 445 22];
            app.FieldofviewLabel.Text = 'Field of view';

            % Create line1
            app.line1 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.line1.Position = [543 340 468 22];
            app.line1.Text = '____________________________________________________________________';

            % Create ViewnornaltoaxexLabel
            app.ViewnornaltoaxexLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.ViewnornaltoaxexLabel.Position = [554 286 457 22];
            app.ViewnornaltoaxexLabel.Text = 'View nornal to axe x';

            % Create SlicexyabLabel_2
            app.SlicexyabLabel_2 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.SlicexyabLabel_2.Position = [554 259 457 22];
            app.SlicexyabLabel_2.Text = 'Slice: x / y, a / b';

            % Create Domainssizeaxbxcum3voxelsizeumLabel
            app.Domainssizeaxbxcum3voxelsizeumLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Domainssizeaxbxcum3voxelsizeumLabel.Position = [554 313 457 22];
            app.Domainssizeaxbxcum3voxelsizeumLabel.Text = 'Domain''s size: a x b x c um3, voxel size: um';

            % Create PhasesLabel
            app.PhasesLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.PhasesLabel.FontSize = 14;
            app.PhasesLabel.FontWeight = 'bold';
            app.PhasesLabel.Position = [543 229 445 22];
            app.PhasesLabel.Text = 'Phases';

            % Create line2
            app.line2 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.line2.Position = [543 214 468 22];
            app.line2.Text = '____________________________________________________________________';

            % Create colphase1
            app.colphase1 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase1.BackgroundColor = [0 0.4471 0.7412];
            app.colphase1.Position = [554 184 66 22];
            app.colphase1.Text = '';

            % Create Phase1volumefractionxLabel
            app.Phase1volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase1volumefractionxLabel.Position = [632 184 379 22];
            app.Phase1volumefractionxLabel.Text = 'Phase 1, volume fraction x';

            % Create colphase2
            app.colphase2 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase2.BackgroundColor = [0 0.4471 0.7412];
            app.colphase2.Position = [554 156 66 22];
            app.colphase2.Text = '';

            % Create Phase2volumefractionxLabel
            app.Phase2volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase2volumefractionxLabel.Position = [632 156 379 22];
            app.Phase2volumefractionxLabel.Text = 'Phase 2, volume fraction x';

            % Create colphase3
            app.colphase3 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase3.BackgroundColor = [0 0.4471 0.7412];
            app.colphase3.Position = [555 128 66 22];
            app.colphase3.Text = '';

            % Create Phase3volumefractionxLabel
            app.Phase3volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase3volumefractionxLabel.Position = [633 128 378 22];
            app.Phase3volumefractionxLabel.Text = 'Phase 3, volume fraction x';

            % Create colphase4
            app.colphase4 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase4.BackgroundColor = [0 0.4471 0.7412];
            app.colphase4.Position = [554 100 66 22];
            app.colphase4.Text = '';

            % Create Phase4volumefractionxLabel
            app.Phase4volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase4volumefractionxLabel.Position = [632 100 379 22];
            app.Phase4volumefractionxLabel.Text = 'Phase 4, volume fraction x';

            % Create colphase5
            app.colphase5 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase5.BackgroundColor = [0 0.4471 0.7412];
            app.colphase5.Position = [554 72 66 22];
            app.colphase5.Text = '';

            % Create Phase5volumefractionxLabel
            app.Phase5volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase5volumefractionxLabel.Position = [632 72 379 22];
            app.Phase5volumefractionxLabel.Text = 'Phase 5, volume fraction x';

            % Create colphase6
            app.colphase6 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase6.BackgroundColor = [0 0.4471 0.7412];
            app.colphase6.Position = [555 44 66 22];
            app.colphase6.Text = '';

            % Create Phase6volumefractionxLabel
            app.Phase6volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase6volumefractionxLabel.Position = [633 44 378 22];
            app.Phase6volumefractionxLabel.Text = 'Phase 6, volume fraction x';

            % Create colphase7
            app.colphase7 = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.colphase7.BackgroundColor = [0 0.4471 0.7412];
            app.colphase7.Position = [554 16 66 22];
            app.colphase7.Text = '';

            % Create Phase7volumefractionxLabel
            app.Phase7volumefractionxLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.Phase7volumefractionxLabel.Position = [632 16 379 22];
            app.Phase7volumefractionxLabel.Text = 'Phase 7, volume fraction x';

            % Create VideoformatDropDownLabel
            app.VideoformatDropDownLabel = uilabel(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.VideoformatDropDownLabel.HorizontalAlignment = 'right';
            app.VideoformatDropDownLabel.Position = [863 13 72 22];
            app.VideoformatDropDownLabel.Text = 'Video format';

            % Create VideoformatDropDown
            app.VideoformatDropDown = uidropdown(app.GreylevelsegmentedimagecomparisonUIFigure);
            app.VideoformatDropDown.Items = {'mpeg-4', 'Uncompressed AVI'};
            app.VideoformatDropDown.Position = [944 13 82 22];
            app.VideoformatDropDown.Value = 'mpeg-4';

            % Show the figure after all components are created
            app.GreylevelsegmentedimagecomparisonUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = microstructure_visualization_segmentation_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.GreylevelsegmentedimagecomparisonUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GreylevelsegmentedimagecomparisonUIFigure)
        end
    end
end