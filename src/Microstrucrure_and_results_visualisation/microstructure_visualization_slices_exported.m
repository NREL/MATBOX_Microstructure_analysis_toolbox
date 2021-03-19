classdef microstructure_visualization_slices_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        ColoroptionsPanel              matlab.ui.container.Panel
        ColormapDropDownLabel          matlab.ui.control.Label
        ColormapDropDown               matlab.ui.control.DropDown
        CustomchangeDropDownLabel      matlab.ui.control.Label
        CustomchangeDropDown           matlab.ui.control.DropDown
        ScaleDropDownLabel             matlab.ui.control.Label
        ScaleDropDown                  matlab.ui.control.DropDown
        MaxEditFieldLabel              matlab.ui.control.Label
        MaxEditField                   matlab.ui.control.NumericEditField
        MinEditField_2Label            matlab.ui.control.Label
        MinEditField                   matlab.ui.control.NumericEditField
        AddcolorbarCheckBox            matlab.ui.control.CheckBox
        TakeoppositevaluesButton       matlab.ui.control.StateButton
        FigurebackgroundDropDownLabel  matlab.ui.control.Label
        FigurebackgroundDropDown       matlab.ui.control.DropDown
        AxisoptionsPanel               matlab.ui.container.Panel
        Viewnormaltoaxe1CheckBox       matlab.ui.control.CheckBox
        Viewnormaltoaxe2CheckBox       matlab.ui.control.CheckBox
        Viewnormaltoaxe3CheckBox       matlab.ui.control.CheckBox
        ShowingslicesisCPUandRAMexpensiveLabel  matlab.ui.control.Label
        DslicesDropDownLabel           matlab.ui.control.Label
        DslicesDropDown                matlab.ui.control.DropDown
        Button_p_1                     matlab.ui.control.Button
        Button_pp_1                    matlab.ui.control.Button
        Slider_1                       matlab.ui.control.Slider
        Button_m_1                     matlab.ui.control.Button
        Button_mm_1                    matlab.ui.control.Button
        Positionalongaxe1Label         matlab.ui.control.Label
        Positionalongaxe2Label         matlab.ui.control.Label
        Slider_2                       matlab.ui.control.Slider
        Button_m_2                     matlab.ui.control.Button
        Button_mm_2                    matlab.ui.control.Button
        Button_p_2                     matlab.ui.control.Button
        Button_pp_2                    matlab.ui.control.Button
        Positionalongaxe3Label         matlab.ui.control.Label
        Slider_3                       matlab.ui.control.Slider
        Button_m_3                     matlab.ui.control.Button
        Button_mm_3                    matlab.ui.control.Button
        Button_p_3                     matlab.ui.control.Button
        Button_pp_3                    matlab.ui.control.Button
        AdjustviewasmovingsliderCheckBox  matlab.ui.control.CheckBox
        FigureandvideooptionsPanel     matlab.ui.container.Panel
        ClicktoselectsavefolderButton  matlab.ui.control.Button
        FilenameEditFieldLabel         matlab.ui.control.Label
        FilenameEditField              matlab.ui.control.EditField
        NosavefolderselectedLabel      matlab.ui.control.Label
        VideoformatDropDownLabel       matlab.ui.control.Label
        VideoformatDropDown            matlab.ui.control.DropDown
        FramepersEditFieldLabel        matlab.ui.control.Label
        FramepersEditField             matlab.ui.control.NumericEditField
        Savefigure_button              matlab.ui.control.Button
        Savevideo_button               matlab.ui.control.Button
        VolumenameLabel                matlab.ui.control.Label
        UIAxes_1                       matlab.ui.control.UIAxes
        UIAxes_2                       matlab.ui.control.UIAxes
        UIAxes_3                       matlab.ui.control.UIAxes
        UIAxes_3D                      matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        volume = [];
        volume_name = [];
        voxelsize = [];
        voxelsize_unitname = [];
        position_slice = [];
        domain_size = [];
        dimension = [];
        color_map = [];
        arrayunit = [];
        global_max = [];
        global_min= [];
        font_name_GUI = [];
        font_size_axe = [];
        axesposition = [];
        major_tick = [];
        savefolder = [];        
    end
    
    methods (Access = private)
        
        function [] = update_normalview_axis(app, direction, axe_)
            if app.dimension==3
                if direction==1
                    slice_ = squeeze(app.volume(app.position_slice(1),:,:));
                    str_xaxis = ['Third axis ' num2str(app.domain_size(3)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                    str_yaxis = ['Second axis ' num2str(app.domain_size(2)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                    str_title = 'View normal to first axis';
                elseif direction==2
                    slice_ = squeeze(app.volume(:,app.position_slice(2),:));
                    str_xaxis = ['Third axis ' num2str(app.domain_size(3)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                    str_yaxis = ['First axis ' num2str(app.domain_size(1)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                    str_title = 'View normal to second axis';
                elseif direction==3
                    slice_ = app.volume(:,:,app.position_slice(3));
                    str_xaxis = ['Second axis ' num2str(app.domain_size(2)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                    str_yaxis = ['First axis ' num2str(app.domain_size(1)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                    str_title = 'View normal to third axis';
                end
                str_subtitle = ['Slice: ' num2str(app.position_slice(direction)) ' / ' num2str(app.domain_size(direction)) ', ' num2str(app.position_slice(direction)*app.voxelsize,'%1.1f') ' / ' num2str(app.domain_size(direction)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
            else
                slice_ = app.volume(:,:,app.position_slice(3));
                str_xaxis = ['Second axis ' num2str(app.domain_size(2)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                str_yaxis = ['First axis ' num2str(app.domain_size(1)*app.voxelsize,'%1.1f') app.voxelsize_unitname];
                str_title = 'View normal to third axis';
            end
            h=pcolor(slice_(:,:,1),'Parent', axe_); % Pcolor: NaN value are transparent
            set(h,'EdgeColor','none');
            colormap(axe_,app.color_map)
            if app.AddcolorbarCheckBox.Value
                h=colorbar(axe_);
                ylabel(h, app.arrayunit);
                set(h,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe);
            else
               colorbar(axe_,'off');
            end
            if strcmp(app.ScaleDropDown.Value,'Global')
                caxis(axe_,[app.global_min app.global_max]);
            elseif strcmp(app.ScaleDropDown.Value,'Local')
                local_min = nanmin(nanmin(slice_));
                local_max = nanmax(nanmax(slice_));
                caxis(axe_,[local_min local_max])
            elseif strcmp(app.ScaleDropDown.Value,'Custom')
                caxis(axe_,[app.MinEditField.Value app.MaxEditField.Value]);
            end
            xlabel(axe_,str_xaxis);
            ylabel(axe_,str_yaxis);
            set(axe_,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe); % - Fontname and fontsize
            if app.dimension==3
                [t,s] = title(axe_,str_title,str_subtitle,'Color','k');
                t.FontSize = app.font_size_axe+2; t.FontName = app.font_name_GUI; s.FontAngle = 'italic'; s.FontSize = app.font_size_axe;
            else
                [t] = title(axe_,str_title,'Color','k');
                t.FontSize = app.font_size_axe+2; t.FontName = app.font_name_GUI;                
            end
            axe_position =axe_.Position; % Get axis position
            set(axe_ ,'Position', axe_position);
            set(axe_,'xtick',[],'ytick',[]); % Remove tick and label
            set(axe_,'xticklabel',[],'yticklabel',[]);
            axis(axe_,'tight'); % Fit the axes box
            axis(axe_,'equal'); % Aspect ratio is 1:1
        end
        
        function [] = Update_axes(app, selection)
            % Plot axis
            if app.dimension==3
                if app.Viewnormaltoaxe1CheckBox.Value && (strcmp(selection,'all') || strcmp(selection,'axe1'))
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if app.Viewnormaltoaxe2CheckBox.Value && (strcmp(selection,'all') || strcmp(selection,'axe2'))
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if app.Viewnormaltoaxe3CheckBox.Value && (strcmp(selection,'all') || strcmp(selection,'axe3'))
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show') || strcmp(selection,'axe3D')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            else
                app.update_normalview_axis(3,app.UIAxes_3);
            end
        end
        
        function [] = update_3d_slices(app,axe_)
            cla(axe_);
            hold(axe_,'on');
            view(axe_,3)
            if strcmp(app.DslicesDropDown.Value,'Show slices')
                slice( double(app.volume), app.position_slice(2),app.position_slice(1),app.position_slice(3),'Parent', axe_);
                if strcmp(app.CustomchangeDropDown.Value,'1st color is white')
                    line_color1='k';
                    line_color2='k';
                    line_color3='k';
                elseif strcmp(app.CustomchangeDropDown.Value,'1st color is black')
                    line_color1='w';
                    line_color2='w';
                    line_color3='w';
                else
                    c = [get(0, 'DefaultAxesColorOrder')];
                    line_color1=c(1,:);
                    line_color2=c(2,:);
                    line_color3=c(3,:);
                end
            elseif strcmp(app.DslicesDropDown.Value,'Show only locations')
                X = [0 app.domain_size(2) app.domain_size(2) 0];
                Y = [app.position_slice(1) app.position_slice(1) app.position_slice(1) app.position_slice(1)];
                Z = [0 0 app.domain_size(3) app.domain_size(3)];
                fill3(X,Y,Z,[0.7 0.7 0.7 0.7],'FaceAlpha',.5,'Parent', axe_)
                
                Y = [0 app.domain_size(1) app.domain_size(1) 0];
                X = [app.position_slice(2) app.position_slice(2) app.position_slice(2) app.position_slice(2)];
                Z = [0 0 app.domain_size(3) app.domain_size(3)];
                fill3(X,Y,Z,[0.6 0.6 0.6 0.6],'FaceAlpha',.5,'Parent', axe_)
                
                X = [0 app.domain_size(2) app.domain_size(2) 0];
                Y = [0 0 app.domain_size(1) app.domain_size(1)];
                Z = [app.position_slice(3) app.position_slice(3) app.position_slice(3) app.position_slice(3)];
                fill3(X,Y,Z,[0.5 0.5 0.5 0.5],'FaceAlpha',0.5,'Parent', axe_)
                
                line_color1='k';
                line_color2='k';
                line_color3='k';
            end
            xlabel(axe_,['Second axis (' app.voxelsize_unitname ')']);
            ylabel(axe_,['First axis (' app.voxelsize_unitname ')']);
            zlabel(axe_,['Third axis (' app.voxelsize_unitname ')']);
            X=[app.position_slice(2) app.position_slice(2)];
            Y=[0 app.domain_size(1)];
            Z=[app.position_slice(3) app.position_slice(3)];
            plot3(X,Y,Z,'--','Color',line_color1,'LineWidth',2,'Parent', axe_)
            X=[0 app.domain_size(2)];
            Y=[app.position_slice(1) app.position_slice(1)];
            Z=[app.position_slice(3) app.position_slice(3)];
            plot3(X,Y,Z,'--','Color',line_color2','LineWidth',2,'Parent', axe_)
            X=[app.position_slice(2) app.position_slice(2)];
            Y=[app.position_slice(1) app.position_slice(1)];
            Z=[0 app.domain_size(3)];
            plot3(X,Y,Z,'--','Color',line_color3,'LineWidth',2,'Parent', axe_)
            grid(axe_,'on'); % Display grid
            set(axe_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
            set(axe_,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe); % - Fontname and fontsize
            axis(axe_,'tight'); % Fit the axes box
            axis(axe_,'equal'); % Aspect ratio is 1:1
            if strcmp(app.DslicesDropDown.Value,'Show slices')
                colormap(axe_,app.color_map)
                if app.AddcolorbarCheckBox.Value
                    h=colorbar(axe_);
                    ylabel(h, app.arrayunit);
                    set(h,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe);
                else
                   colorbar(axe_,'off'); 
                end
                if strcmp(app.ScaleDropDown.Value,'Global') || strcmp(app.ScaleDropDown.Value,'Local')
                    caxis(axe_,[app.global_min app.global_max])
                elseif strcmp(app.ScaleDropDown.Value,'Custom')
                    caxis(axe_,[app.MinEditField.Value cusapp.MaxEditField.Valuetom_max])
                end
            else
                colormap(axe_,gray)
            end
            shading(axe_,'flat');
            set(axe_,'XTickMode','auto')
            set(axe_,'YTickMode','auto')
            set(axe_,'ZTickMode','auto')
            x_value = get(axe_,'XTick');
            set(axe_,'XtickLabel',round(x_value*app.voxelsize,1));
            y_value = get(axe_,'YTick');
            set(axe_,'YtickLabel',round(y_value*app.voxelsize,1));
            z_value = get(axe_,'ZTick');
            set(axe_,'ZtickLabel',round(z_value*app.voxelsize,1));
            hold(axe_,'off');
        end
        
        function [] = update_position(app)
            num_axis = app.Viewnormaltoaxe1CheckBox.Value + app.Viewnormaltoaxe2CheckBox.Value + app.Viewnormaltoaxe3CheckBox.Value;
            if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                num_axis=num_axis+1;
            end
            if num_axis==4 || num_axis==3 % Reset
                app.UIAxes_1.Position = app.axesposition(1).pos;
                app.UIAxes_2.Position = app.axesposition(2).pos;
                app.UIAxes_3.Position = app.axesposition(3).pos;
                app.UIAxes_3D.Position = app.axesposition(4).pos;
            elseif num_axis==2
                if app.Viewnormaltoaxe1CheckBox.Value && app.Viewnormaltoaxe2CheckBox.Value
                    app.UIAxes_1.Position = [app.axesposition(1).x0_a app.axesposition(1).y0_a app.axesposition(1).dx_a app.axesposition(1).dy_a];
                    app.UIAxes_2.Position = [app.axesposition(1).x0_b app.axesposition(1).y0_b app.axesposition(1).dx_b app.axesposition(1).dy_b];
                elseif app.Viewnormaltoaxe1CheckBox.Value && app.Viewnormaltoaxe3CheckBox.Value
                    app.UIAxes_1.Position = [app.axesposition(1).x0_a app.axesposition(1).y0_a app.axesposition(1).dx_a app.axesposition(1).dy_a];
                    app.UIAxes_3.Position = [app.axesposition(1).x0_b app.axesposition(1).y0_b app.axesposition(1).dx_b app.axesposition(1).dy_b];
                elseif app.Viewnormaltoaxe1CheckBox.Value && ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.UIAxes_1.Position = [app.axesposition(1).x0_a app.axesposition(1).y0_a app.axesposition(1).dx_a app.axesposition(1).dy_a];
                    app.UIAxes_3D.Position = [app.axesposition(1).x0_b app.axesposition(1).y0_b app.axesposition(1).dx_b app.axesposition(1).dy_b];
                elseif app.Viewnormaltoaxe2CheckBox.Value && app.Viewnormaltoaxe3CheckBox.Value
                    app.UIAxes_2.Position = [app.axesposition(1).x0_a app.axesposition(1).y0_a app.axesposition(1).dx_a app.axesposition(1).dy_a];
                    app.UIAxes_3.Position = [app.axesposition(1).x0_b app.axesposition(1).y0_b app.axesposition(1).dx_b app.axesposition(1).dy_b];
                elseif app.Viewnormaltoaxe2CheckBox.Value && ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.UIAxes_2.Position = [app.axesposition(1).x0_a app.axesposition(1).y0_a app.axesposition(1).dx_a app.axesposition(1).dy_a];
                    app.UIAxes_3D.Position = [app.axesposition(1).x0_b app.axesposition(1).y0_b app.axesposition(1).dx_b app.axesposition(1).dy_b];
                elseif app.Viewnormaltoaxe3CheckBox.Value && ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.UIAxes_3.Position = [app.axesposition(1).x0_a app.axesposition(1).y0_a app.axesposition(1).dx_a app.axesposition(1).dy_a];
                    app.UIAxes_3D.Position = [app.axesposition(1).x0_b app.axesposition(1).y0_b app.axesposition(1).dx_b app.axesposition(1).dy_b];
                end
            elseif num_axis==1
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.UIAxes_1.Position = [app.axesposition(1).x0 app.axesposition(1).y0 app.axesposition(1).dx app.axesposition(1).dy];
                elseif app.Viewnormaltoaxe2CheckBox.Value
                    app.UIAxes_2.Position = [app.axesposition(1).x0 app.axesposition(1).y0 app.axesposition(1).dx app.axesposition(1).dy];
                elseif app.Viewnormaltoaxe3CheckBox.Value                    
                    app.UIAxes_3.Position = [app.axesposition(1).x0 app.axesposition(1).y0 app.axesposition(1).dx app.axesposition(1).dy];
                else
                    app.UIAxes_3D.Position = [app.axesposition(1).x0 app.axesposition(1).y0 app.axesposition(1).dx app.axesposition(1).dy];
                end                
            end
        end
        
        function [] = savemode(app,bool)
            if bool
                app.ColoroptionsPanel.Visible = 'off';
                app.AxisoptionsPanel.Visible = 'off';
                app.FigureandvideooptionsPanel.Visible = 'off';
            else
                app.ColoroptionsPanel.Visible = 'on';
                app.AxisoptionsPanel.Visible = 'on';
                app.FigureandvideooptionsPanel.Visible = 'on';                
            end

        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, volume, volume_name, voxelsize, voxelsize_unitname, arrayunit)
            app.volume = volume;
            app.volume_name = volume_name;
            app.voxelsize = voxelsize;
            app.voxelsize_unitname = voxelsize_unitname;
            app.arrayunit = arrayunit;
            [~,name,~] = fileparts(app.volume_name);
            app.FilenameEditField.Value = name;
            if strcmp(app.voxelsize_unitname,'um') || strcmp(app.voxelsize_unitname,'micrometer') || strcmp(app.voxelsize_unitname,'micrometers')
                app.voxelsize_unitname = [char(956) 'm'];
            end
            if strcmp(app.arrayunit,'um') || strcmp(app.arrayunit,'micrometer') || strcmp(app.arrayunit,'micrometers')
                app.arrayunit = [char(956) 'm'];
            end
            app.VolumenameLabel.Text = app.volume_name;
            app.ColormapDropDown.Items = {'gray','copper','bone','jet','turbo','parula','default','hsv','cool','random'};
            app.global_min = double( min(min(min(app.volume))) );
            app.global_max = double( max(max(max(app.volume))) );
            app.MinEditField.Value = app.global_min;
            app.MaxEditField.Value = app.global_max;
            app.domain_size = size(app.volume);
            app.dimension = length(app.domain_size);
            if app.dimension==2
                app.UIAxes_1.Visible = 'off'; app.UIAxes_2.Visible = 'off'; app.UIAxes_3D.Visible = 'off';
                app.AxisoptionsPanel.Enable = 'off';
                app.position_slice=ones(1,3);
            else
                app.position_slice = round(app.domain_size/2);
                app.Slider_1.Limits = [1 app.domain_size(1)];
                app.Slider_2.Limits = [1 app.domain_size(2)];
                app.Slider_3.Limits = [1 app.domain_size(3)];
                app.Slider_1.Value = app.position_slice(1);
                app.Slider_2.Value = app.position_slice(2);
                app.Slider_3.Value = app.position_slice(3);
                app.major_tick = app.domain_size/10;
            end
            app.color_map = eval(app.ColormapDropDown.Value);
            app.font_name_GUI = 'Helvetica';
            app.font_size_axe = 12;

            % Plot axis
            app.UIAxes_1.Units = 'normalized';
            app.UIAxes_2.Units = 'normalized';
            app.UIAxes_3.Units = 'normalized';
            app.UIAxes_3D.Units = 'normalized';
            if app.dimension==3
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if app.Viewnormaltoaxe2CheckBox.Value
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if app.Viewnormaltoaxe3CheckBox.Value
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
                app.axesposition(1).pos = app.UIAxes_1.Position;
                app.axesposition(2).pos =app.UIAxes_2.Position;
                app.axesposition(3).pos =app.UIAxes_3.Position;
                app.axesposition(4).pos =app.UIAxes_3D.Position;
                % Position for 1 axe displayed
                app.axesposition(1).x0 = app.axesposition(3).pos(1)+0.02;
                app.axesposition(1).y0 = app.axesposition(3).pos(2)+0.02;
                app.axesposition(1).dx = app.axesposition(1).pos(3) + app.axesposition(2).pos(3) + (app.axesposition(2).pos(1) - (app.axesposition(1).pos(1)+app.axesposition(1).pos(3) ) ) - 0.04;
                app.axesposition(1).dy = app.axesposition(3).pos(2) + app.axesposition(1).pos(2) + (app.axesposition(1).pos(2) - (app.axesposition(3).pos(2)+app.axesposition(3).pos(2) ) ) - 0.04;
                % Position for 2 axes displayed
                app.axesposition(1).x0_a = app.axesposition(3).pos(1);
                app.axesposition(1).y0_a = app.axesposition(3).pos(2)+0.02;                
                app.axesposition(1).dx_a = app.axesposition(1).pos(3);
                app.axesposition(1).dy_a = app.axesposition(3).pos(2) + app.axesposition(1).pos(2) + (app.axesposition(1).pos(2) - (app.axesposition(3).pos(2)+app.axesposition(3).pos(2) ) ) - 0.04;
                app.axesposition(1).x0_b = app.axesposition(2).pos(1);
                app.axesposition(1).y0_b = app.axesposition(3).pos(2)+0.02;                
                app.axesposition(1).dx_b = app.axesposition(2).pos(3);
                app.axesposition(1).dy_b = app.axesposition(3).pos(2) + app.axesposition(1).pos(2) + (app.axesposition(1).pos(2) - (app.axesposition(3).pos(2)+app.axesposition(3).pos(2) ) ) - 0.04;;
            else
                app.update_normalview_axis(3,app.UIAxes_3);
                app.axesposition(1).pos =app.UIAxes_1.Position;
                app.axesposition(2).pos =app.UIAxes_2.Position;
                app.axesposition(3).pos =app.UIAxes_3.Position;
                app.axesposition(4).pos =app.UIAxes_3D.Position;
                app.axesposition(1).x0 = app.axesposition(3).pos(1)+0.02;
                app.axesposition(1).y0 = app.axesposition(3).pos(2)+0.02;
                app.axesposition(1).dx = app.axesposition(1).pos(3) + app.axesposition(2).pos(3) + (app.axesposition(2).pos(1) - (app.axesposition(1).pos(1)+app.axesposition(1).pos(3) ) ) - 0.04;
                app.axesposition(1).dy = app.axesposition(3).pos(2) + app.axesposition(1).pos(2) + (app.axesposition(1).pos(2) - (app.axesposition(3).pos(2)+app.axesposition(3).pos(2) ) ) - 0.04;
                app.UIAxes_3.Position = [app.axesposition(1).x0 app.axesposition(1).y0 app.axesposition(1).dx app.axesposition(1).dy];
            end 
        end

        % Value changed function: ColormapDropDown
        function ColormapDropDownValueChanged(app, event)
            if strcmp(app.ColormapDropDown.Value,'default')
                app.color_map = [colororder; rand(1000,3)];
            elseif strcmp(app.ColormapDropDown.Value,'random')
                app.color_map = rand(1000,3);
            else
                app.color_map = eval(app.ColormapDropDown.Value);
            end
            app.Update_axes('all');            
        end

        % Value changed function: CustomchangeDropDown
        function CustomchangeDropDownValueChanged(app, event)
            if strcmp(app.CustomchangeDropDown.Value,'Default colormap')
                app.ColormapDropDownValueChanged;
            elseif strcmp(app.CustomchangeDropDown.Value,'1st color is white')
                app.color_map(1,:) = 1; % background color is white
            elseif strcmp(app.CustomchangeDropDown.Value,'1st color is black')
                app.color_map(1,:) = 0; % background color is dark
            end
            app.Update_axes('all'); 
        end

        % Value changed function: ScaleDropDown
        function ScaleDropDownValueChanged(app, event)
            if strcmp(app.ScaleDropDown.Value,'Custom')
                app.MinEditField.Enable = 'on'; app.MaxEditField.Enable = 'on';
                app.MinEditField_2Label.Enable = 'on'; app.MaxEditFieldLabel.Enable = 'on';
            else
                app.MinEditField.Enable = 'off'; app.MaxEditField.Enable = 'off';
                app.MinEditField_2Label.Enable = 'off'; app.MaxEditFieldLabel.Enable = 'off';
            end
            app.Update_axes('all');
        end

        % Value changed function: AddcolorbarCheckBox
        function AddcolorbarCheckBoxValueChanged(app, event)
            if app.AddcolorbarCheckBox.Value
                if strcmp(app.UIAxes_1.Visible,'on')
                    h=colorbar(app.UIAxes_1);
                    ylabel(h, app.arrayunit);
                    set(h,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe);
                end
                if strcmp(app.UIAxes_2.Visible,'on')
                    h=colorbar(app.UIAxes_2);
                    ylabel(h, app.arrayunit);
                    set(h,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe);
                end
                if strcmp(app.UIAxes_3.Visible,'on')
                    h=colorbar(app.UIAxes_3);
                    ylabel(h, app.arrayunit);
                    set(h,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe);
                end
                if strcmp(app.UIAxes_3D.Visible,'on') && strcmp(app.DslicesDropDown.Value,'Show slices')
                    h=colorbar(app.UIAxes_3D);
                    ylabel(h, app.arrayunit);
                    set(h,'FontName',app.font_name_GUI,'FontSize',app.font_size_axe);
                end                
            else
                colorbar(app.UIAxes_1,'off');
                colorbar(app.UIAxes_2,'off');
                colorbar(app.UIAxes_3,'off');
                colorbar(app.UIAxes_3D,'off');
            end
        end

        % Value changed function: Viewnormaltoaxe1CheckBox
        function Viewnormaltoaxe1CheckBoxValueChanged(app, event)
            if app.Viewnormaltoaxe1CheckBox.Value
                app.UIAxes_1.Visible = 'on'; app.update_position; app.update_normalview_axis(1,app.UIAxes_1);
            else
                app.UIAxes_1.Visible = 'off'; cla(app.UIAxes_1); app.update_position;
            end
        end

        % Value changed function: MaxEditField, MinEditField
        function MinEditField_2ValueChanged(app, event)
            min_edit = app.MinEditField.Value;
            max_edit = app.MaxEditField.Value;
            if max_edit<=min_edit || min_edit<app.global_min || max_edit>app.global_max % Reset to global min-max
                app.MinEditField.Value = app.global_min;
                app.MaxEditField.Value = app.global_max;
            end
            app.Update_axes('all'); 
        end

        % Button pushed function: ClicktoselectsavefolderButton
        function ClicktoselectsavefolderButtonPushed(app, event)
            str_dialogbox = 'Select save folder';
            res = uigetdir(matlabroot,str_dialogbox); % Open dialog box
            if res==0
                % User clicked cancel button or closed the dialog box
                app.NosavefolderselectedLabel.Text = 'No save folder selected';
                app.Savefigure_button.Enable = 'off';
                if app.dimension==3
                    app.Savevideo_button.Enable = 'off';
                end
            else
                if ispc
                    app.savefolder = [res '\'];
                else
                    app.savefolder = [res '/'];
                end
                app.NosavefolderselectedLabel.Text = app.savefolder;
                app.Savefigure_button.Enable = 'on';
                app.Savevideo_button.Enable = 'on';
            end            
            
        end

        % Button pushed function: Savefigure_button
        function Savefigure_buttonButtonPushed(app, event)
            if app.dimension==3
                figure_filename = [app.FilenameEditField.Value '_axe1_slice' num2str(app.position_slice(1)) '_axe2_slice' num2str(app.position_slice(2)) '_axe3_slice' num2str(app.position_slice(3)) '.png'];
            else
                figure_filename = app.FilenameEditField.Value;
            end
            % Can't save an app design figure yet (below won't work in 2020b)
            % function_savefig(app.GreylevelsegmentedimagecomparisonUIFigure, app.savefolder, figure_filename, OPTIONS);
            % Workaround proposed by Matlab will only capture the axis
            exportgraphics(app.UIFigure,[app.savefolder figure_filename]); % Better quality than screencapture
            % Instead use screencapture
            % Yair Altman (2021). ScreenCapture - get a screen-capture of a figure frame or component
            % https://www.mathworks.com/matlabcentral/fileexchange/24323-screencapture-get-a-screen-capture-of-a-figure-frame-or-component
            % MATLAB Central File Exchange. Retrieved March 4, 2021.
            %imageData = screencapture(app.UIFigure);
            %imwrite(imageData,[app.savefolder figure_filename]);
        end

        % Button pushed function: Savevideo_button
        function Savevideo_buttonButtonPushed(app, event)
            app.savemode(true); pause(0.1);
            video_filename = app.FilenameEditField.Value;
            video_format = app.VideoformatDropDown.Value;
            video_handle = VideoWriter([app.savefolder video_filename],video_format);
            if strcmp(video_format,'mpeg-4')
                set(video_handle,'Quality',100); % Set video quality
            end
            % Set video framerate
            set(video_handle,'FrameRate',app.FramepersEditField.Value);
            % Open video
            open(video_handle)
            frame_number = 0;
            if app.Viewnormaltoaxe1CheckBox.Value
                for value=1:1:app.domain_size(1)
                    frame_number = frame_number+1;
                    % Update figure
                    app.position_slice(1) = value;
                    app.Slider_1.Value = value;
                    app.update_normalview_axis(1,app.UIAxes_1);
                    if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                        app.update_3d_slices(app.UIAxes_3D);
                    end
                    stored_frame(frame_number) = getframe(app.UIFigure);
                    writeVideo(video_handle,stored_frame(frame_number))
                    pause(0.01);
                end
            end
            if app.Viewnormaltoaxe2CheckBox.Value
                for value=1:1:app.domain_size(2)
                    frame_number = frame_number+1;
                    % Update figure
                    app.position_slice(2) = value;
                    app.Slider_2.Value = value;
                    app.update_normalview_axis(2,app.UIAxes_2);
                    if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                        app.update_3d_slices(app.UIAxes_3D);
                    end
                    stored_frame(frame_number) = getframe(app.UIFigure);
                    writeVideo(video_handle,stored_frame(frame_number))
                    pause(0.01);
                end
            end       
            if app.Viewnormaltoaxe3CheckBox.Value
                for value=1:1:app.domain_size(3)
                    frame_number = frame_number+1;
                    % Update figure
                    app.position_slice(3) = value;
                    app.Slider_3.Value = value;
                    app.update_normalview_axis(3,app.UIAxes_3);
                    if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                        app.update_3d_slices(app.UIAxes_3D);
                    end
                    stored_frame(frame_number) = getframe(app.UIFigure);
                    writeVideo(video_handle,stored_frame(frame_number))
                    pause(0.01);
                end
            end                
            % Close video
            close(video_handle)
            app.savemode(false);
        end

        % Value changed function: DslicesDropDown
        function DslicesDropDownValueChanged(app, event)
            if strcmp(app.DslicesDropDown.Value,'Do not show')
                app.UIAxes_3D.Visible = 'off'; cla(app.UIAxes_3D); app.update_position;
            else
                app.UIAxes_3D.Visible = 'on'; app.update_position; app.update_3d_slices(app.UIAxes_3D);
            end     
        end

        % Value changed function: Slider_1
        function Slider_1ValueChanged(app, event)
            value = round(app.Slider_1.Value);
            app.position_slice(1) = value;
            app.Slider_1.Value = value;
            if app.Viewnormaltoaxe1CheckBox.Value
                app.update_normalview_axis(1,app.UIAxes_1);
            end
            if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                app.update_3d_slices(app.UIAxes_3D);
            end
        end

        % Value changing function: Slider_1
        function Slider_1ValueChanging(app, event)
            if app.AdjustviewasmovingsliderCheckBox.Value
                value = round(event.Value);
                app.position_slice(1) = value;
                app.Slider_1.Value = value;
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Value changed function: Slider_2
        function Slider_2ValueChanged(app, event)
            value = round(app.Slider_2.Value);
            app.position_slice(2) = value;
            app.Slider_2.Value = value;
            if app.Viewnormaltoaxe2CheckBox.Value
                app.update_normalview_axis(2,app.UIAxes_2);
            end
            if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                app.update_3d_slices(app.UIAxes_3D);
            end            
        end

        % Value changing function: Slider_2
        function Slider_2ValueChanging(app, event)
            if app.AdjustviewasmovingsliderCheckBox.Value
                value = round(event.Value);
                app.position_slice(2) = value;
                app.Slider_2.Value = value;
                if app.Viewnormaltoaxe2CheckBox.Value
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Value changed function: Slider_3
        function Slider_3ValueChanged(app, event)
            value = round(app.Slider_3.Value);
            app.position_slice(3) = value;
            app.Slider_3.Value = value;
            if app.Viewnormaltoaxe3CheckBox.Value
                app.update_normalview_axis(3,app.UIAxes_3);
            end
            if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                app.update_3d_slices(app.UIAxes_3D);
            end
        end

        % Value changing function: Slider_3
        function Slider_3ValueChanging(app, event)
            if app.AdjustviewasmovingsliderCheckBox.Value
                value = round(event.Value);
                app.position_slice(3) = value;
                app.Slider_3.Value = value;
                if app.Viewnormaltoaxe3CheckBox.Value
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_p_1
        function Button_p_1Pushed(app, event)
            value = round(app.Slider_1.Value);
            new_value = min(value+1, app.domain_size(1));
            app.position_slice(1) = new_value;
            app.Slider_1.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_pp_1
        function Button_pp_1Pushed(app, event)
            value = round(app.Slider_1.Value);
            new_value = min(round(value+app.major_tick(1)), app.domain_size(1));
            app.position_slice(1) = new_value;
            app.Slider_1.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_m_1
        function Button_m_1Pushed(app, event)
            value = round(app.Slider_1.Value);
            new_value = max(value-1, 1);
            app.position_slice(1) = new_value;
            app.Slider_1.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end            
        end

        % Button pushed function: Button_mm_1
        function Button_mm_1Pushed(app, event)
            value = round(app.Slider_1.Value);
            new_value = max(round(value-app.major_tick(1)), 1);
            app.position_slice(1) = new_value;
            app.Slider_1.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe1CheckBox.Value
                    app.update_normalview_axis(1,app.UIAxes_1);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end               
        end

        % Button pushed function: Button_p_2
        function Button_p_2Pushed(app, event)
            value = round(app.Slider_2.Value);
            new_value = min(value+1, app.domain_size(2));
            app.position_slice(2) = new_value;
            app.Slider_2.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe2CheckBox.Value
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end            
        end

        % Button pushed function: Button_pp_2
        function Button_pp_2Pushed(app, event)
            value = round(app.Slider_2.Value);
            new_value = min(round(value+app.major_tick(2)), app.domain_size(2));
            app.position_slice(2) = new_value;
            app.Slider_2.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe2CheckBox.Value
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end            
        end

        % Button pushed function: Button_m_2
        function Button_m_2Pushed(app, event)
            value = round(app.Slider_2.Value);
            new_value = max(value-1, 1);
            app.position_slice(2) = new_value;
            app.Slider_2.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe2CheckBox.Value
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_mm_2
        function Button_mm_2Pushed(app, event)
            value = round(app.Slider_2.Value);
            new_value = max(round(value-app.major_tick(2)), 1);
            app.position_slice(2) = new_value;
            app.Slider_2.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe2CheckBox.Value
                    app.update_normalview_axis(2,app.UIAxes_2);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end               
        end

        % Button pushed function: Button_p_3
        function Button_p_3Pushed(app, event)
            value = round(app.Slider_3.Value);
            new_value = min(value+1, app.domain_size(3));
            app.position_slice(3) = new_value;
            app.Slider_3.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe3CheckBox.Value
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_pp_3
        function Button_pp_3Pushed(app, event)
            value = round(app.Slider_3.Value);
            new_value = min(round(value+app.major_tick(3)), app.domain_size(3));
            app.position_slice(3) = new_value;
            app.Slider_3.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe3CheckBox.Value
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_m_3
        function Button_m_3Pushed(app, event)
            value = round(app.Slider_3.Value);
            new_value = max(value-1, 1);
            app.position_slice(3) = new_value;
            app.Slider_3.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe3CheckBox.Value
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end
        end

        % Button pushed function: Button_mm_3
        function Button_mm_3Pushed(app, event)
            value = round(app.Slider_3.Value);
            new_value = max(round(value-app.major_tick(3)), 1);
            app.position_slice(3) = new_value;
            app.Slider_3.Value = new_value;
            if value~=new_value
                if app.Viewnormaltoaxe3CheckBox.Value
                    app.update_normalview_axis(3,app.UIAxes_3);
                end
                if ~strcmp(app.DslicesDropDown.Value,'Do not show')
                    app.update_3d_slices(app.UIAxes_3D);
                end
            end             
        end

        % Value changed function: TakeoppositevaluesButton
        function TakeoppositevaluesButtonValueChanged(app, event)
            if app.TakeoppositevaluesButton.Value
                app.volume = -double(app.volume);
            else
                app.volume = -double(app.volume); % reverse
            end
            app.global_min = double( min(min(min(app.volume))) );
            app.global_max = double( max(max(max(app.volume))) );
            app.MinEditField.Value = app.global_min;
            app.MaxEditField.Value = app.global_max;      
            app.Update_axes('all');
        end

        % Value changed function: Viewnormaltoaxe2CheckBox
        function Viewnormaltoaxe2CheckBoxValueChanged(app, event)
            if app.Viewnormaltoaxe2CheckBox.Value
                app.UIAxes_2.Visible = 'on'; app.update_position; app.update_normalview_axis(2,app.UIAxes_2);
            else
                app.UIAxes_2.Visible = 'off'; cla(app.UIAxes_2); app.update_position;
            end
        end

        % Value changed function: Viewnormaltoaxe3CheckBox
        function Viewnormaltoaxe3CheckBoxValueChanged(app, event)
            if app.Viewnormaltoaxe3CheckBox.Value
                app.UIAxes_3.Visible = 'on'; app.update_position; app.update_normalview_axis(3,app.UIAxes_3);
            else
                app.UIAxes_3.Visible = 'off'; cla(app.UIAxes_3); app.update_position;
            end
        end

        % Value changed function: FigurebackgroundDropDown
        function FigurebackgroundDropDownValueChanged(app, event)
            value = app.FigurebackgroundDropDown.Value;
            if strcmp(value,'White')
                app.UIFigure.Color = [1.0 1.0 1.0];
            elseif strcmp(value,'Light grey')
                app.UIFigure.Color = [0.94 0.94 0.94];
            elseif strcmp(value,'Yellow')
                app.UIFigure.Color = [0.94 0.94 0];                
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1556 939];
            app.UIFigure.Name = 'MATLAB App';

            % Create ColoroptionsPanel
            app.ColoroptionsPanel = uipanel(app.UIFigure);
            app.ColoroptionsPanel.Title = 'Color options';
            app.ColoroptionsPanel.BackgroundColor = [1 1 1];
            app.ColoroptionsPanel.FontWeight = 'bold';
            app.ColoroptionsPanel.Position = [1301 741 256 199];

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.ColoroptionsPanel);
            app.ColormapDropDownLabel.Position = [7 152 65 23];
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.ColoroptionsPanel);
            app.ColormapDropDown.ValueChangedFcn = createCallbackFcn(app, @ColormapDropDownValueChanged, true);
            app.ColormapDropDown.Position = [120 153 128 22];

            % Create CustomchangeDropDownLabel
            app.CustomchangeDropDownLabel = uilabel(app.ColoroptionsPanel);
            app.CustomchangeDropDownLabel.Position = [7 124 89 23];
            app.CustomchangeDropDownLabel.Text = 'Custom change';

            % Create CustomchangeDropDown
            app.CustomchangeDropDown = uidropdown(app.ColoroptionsPanel);
            app.CustomchangeDropDown.Items = {'Default colormap', '1st color is white', '1st color is black'};
            app.CustomchangeDropDown.ValueChangedFcn = createCallbackFcn(app, @CustomchangeDropDownValueChanged, true);
            app.CustomchangeDropDown.Position = [120 125 128 22];
            app.CustomchangeDropDown.Value = 'Default colormap';

            % Create ScaleDropDownLabel
            app.ScaleDropDownLabel = uilabel(app.ColoroptionsPanel);
            app.ScaleDropDownLabel.Position = [7 93 36 23];
            app.ScaleDropDownLabel.Text = 'Scale';

            % Create ScaleDropDown
            app.ScaleDropDown = uidropdown(app.ColoroptionsPanel);
            app.ScaleDropDown.Items = {'Global', 'Local', 'Custom'};
            app.ScaleDropDown.ValueChangedFcn = createCallbackFcn(app, @ScaleDropDownValueChanged, true);
            app.ScaleDropDown.Position = [72 94 105 22];
            app.ScaleDropDown.Value = 'Global';

            % Create MaxEditFieldLabel
            app.MaxEditFieldLabel = uilabel(app.ColoroptionsPanel);
            app.MaxEditFieldLabel.Enable = 'off';
            app.MaxEditFieldLabel.Position = [98 65 36 23];
            app.MaxEditFieldLabel.Text = 'Max';

            % Create MaxEditField
            app.MaxEditField = uieditfield(app.ColoroptionsPanel, 'numeric');
            app.MaxEditField.ValueChangedFcn = createCallbackFcn(app, @MinEditField_2ValueChanged, true);
            app.MaxEditField.Enable = 'off';
            app.MaxEditField.Position = [133 66 44 22];

            % Create MinEditField_2Label
            app.MinEditField_2Label = uilabel(app.ColoroptionsPanel);
            app.MinEditField_2Label.Enable = 'off';
            app.MinEditField_2Label.Position = [7 65 36 23];
            app.MinEditField_2Label.Text = 'Min';

            % Create MinEditField
            app.MinEditField = uieditfield(app.ColoroptionsPanel, 'numeric');
            app.MinEditField.ValueChangedFcn = createCallbackFcn(app, @MinEditField_2ValueChanged, true);
            app.MinEditField.Enable = 'off';
            app.MinEditField.Position = [42 66 44 22];

            % Create AddcolorbarCheckBox
            app.AddcolorbarCheckBox = uicheckbox(app.ColoroptionsPanel);
            app.AddcolorbarCheckBox.ValueChangedFcn = createCallbackFcn(app, @AddcolorbarCheckBoxValueChanged, true);
            app.AddcolorbarCheckBox.Text = 'Add colorbar';
            app.AddcolorbarCheckBox.Position = [7 37 90 23];

            % Create TakeoppositevaluesButton
            app.TakeoppositevaluesButton = uibutton(app.ColoroptionsPanel, 'state');
            app.TakeoppositevaluesButton.ValueChangedFcn = createCallbackFcn(app, @TakeoppositevaluesButtonValueChanged, true);
            app.TakeoppositevaluesButton.Text = 'Take opposite values';
            app.TakeoppositevaluesButton.Position = [119 38 128 22];

            % Create FigurebackgroundDropDownLabel
            app.FigurebackgroundDropDownLabel = uilabel(app.ColoroptionsPanel);
            app.FigurebackgroundDropDownLabel.HorizontalAlignment = 'right';
            app.FigurebackgroundDropDownLabel.Position = [7 10 106 22];
            app.FigurebackgroundDropDownLabel.Text = 'Figure background';

            % Create FigurebackgroundDropDown
            app.FigurebackgroundDropDown = uidropdown(app.ColoroptionsPanel);
            app.FigurebackgroundDropDown.Items = {'Light grey', 'White', 'Yellow'};
            app.FigurebackgroundDropDown.ValueChangedFcn = createCallbackFcn(app, @FigurebackgroundDropDownValueChanged, true);
            app.FigurebackgroundDropDown.Position = [128 10 119 22];
            app.FigurebackgroundDropDown.Value = 'Light grey';

            % Create AxisoptionsPanel
            app.AxisoptionsPanel = uipanel(app.UIFigure);
            app.AxisoptionsPanel.Title = 'Axis options';
            app.AxisoptionsPanel.BackgroundColor = [1 1 1];
            app.AxisoptionsPanel.FontWeight = 'bold';
            app.AxisoptionsPanel.Position = [1301 277 256 465];

            % Create Viewnormaltoaxe1CheckBox
            app.Viewnormaltoaxe1CheckBox = uicheckbox(app.AxisoptionsPanel);
            app.Viewnormaltoaxe1CheckBox.ValueChangedFcn = createCallbackFcn(app, @Viewnormaltoaxe1CheckBoxValueChanged, true);
            app.Viewnormaltoaxe1CheckBox.Text = 'View normal to axe 1';
            app.Viewnormaltoaxe1CheckBox.Position = [5 423 134 22];
            app.Viewnormaltoaxe1CheckBox.Value = true;

            % Create Viewnormaltoaxe2CheckBox
            app.Viewnormaltoaxe2CheckBox = uicheckbox(app.AxisoptionsPanel);
            app.Viewnormaltoaxe2CheckBox.ValueChangedFcn = createCallbackFcn(app, @Viewnormaltoaxe2CheckBoxValueChanged, true);
            app.Viewnormaltoaxe2CheckBox.Text = 'View normal to axe 2';
            app.Viewnormaltoaxe2CheckBox.Position = [5 396 134 22];
            app.Viewnormaltoaxe2CheckBox.Value = true;

            % Create Viewnormaltoaxe3CheckBox
            app.Viewnormaltoaxe3CheckBox = uicheckbox(app.AxisoptionsPanel);
            app.Viewnormaltoaxe3CheckBox.ValueChangedFcn = createCallbackFcn(app, @Viewnormaltoaxe3CheckBoxValueChanged, true);
            app.Viewnormaltoaxe3CheckBox.Text = 'View normal to axe 3';
            app.Viewnormaltoaxe3CheckBox.Position = [5 369 134 22];
            app.Viewnormaltoaxe3CheckBox.Value = true;

            % Create ShowingslicesisCPUandRAMexpensiveLabel
            app.ShowingslicesisCPUandRAMexpensiveLabel = uilabel(app.AxisoptionsPanel);
            app.ShowingslicesisCPUandRAMexpensiveLabel.FontSize = 10;
            app.ShowingslicesisCPUandRAMexpensiveLabel.FontAngle = 'italic';
            app.ShowingslicesisCPUandRAMexpensiveLabel.Position = [5 310 199 22];
            app.ShowingslicesisCPUandRAMexpensiveLabel.Text = 'Showing slices is CPU and RAM expensive';

            % Create DslicesDropDownLabel
            app.DslicesDropDownLabel = uilabel(app.AxisoptionsPanel);
            app.DslicesDropDownLabel.HorizontalAlignment = 'right';
            app.DslicesDropDownLabel.Position = [5 337 54 22];
            app.DslicesDropDownLabel.Text = '3D slices';

            % Create DslicesDropDown
            app.DslicesDropDown = uidropdown(app.AxisoptionsPanel);
            app.DslicesDropDown.Items = {'Show only locations', 'Show slices', 'Do not show'};
            app.DslicesDropDown.ValueChangedFcn = createCallbackFcn(app, @DslicesDropDownValueChanged, true);
            app.DslicesDropDown.Position = [75 337 172 22];
            app.DslicesDropDown.Value = 'Show only locations';

            % Create Button_p_1
            app.Button_p_1 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_p_1.ButtonPushedFcn = createCallbackFcn(app, @Button_p_1Pushed, true);
            app.Button_p_1.Position = [219 248 29 22];
            app.Button_p_1.Text = '+';

            % Create Button_pp_1
            app.Button_pp_1 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_pp_1.ButtonPushedFcn = createCallbackFcn(app, @Button_pp_1Pushed, true);
            app.Button_pp_1.Position = [219 224 29 22];
            app.Button_pp_1.Text = '++';

            % Create Slider_1
            app.Slider_1 = uislider(app.AxisoptionsPanel);
            app.Slider_1.ValueChangedFcn = createCallbackFcn(app, @Slider_1ValueChanged, true);
            app.Slider_1.ValueChangingFcn = createCallbackFcn(app, @Slider_1ValueChanging, true);
            app.Slider_1.Position = [53 256 151 3];

            % Create Button_m_1
            app.Button_m_1 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_m_1.ButtonPushedFcn = createCallbackFcn(app, @Button_m_1Pushed, true);
            app.Button_m_1.Position = [10 248 28 22];
            app.Button_m_1.Text = '-';

            % Create Button_mm_1
            app.Button_mm_1 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_mm_1.ButtonPushedFcn = createCallbackFcn(app, @Button_mm_1Pushed, true);
            app.Button_mm_1.Position = [10 224 28 22];
            app.Button_mm_1.Text = '--';

            % Create Positionalongaxe1Label
            app.Positionalongaxe1Label = uilabel(app.AxisoptionsPanel);
            app.Positionalongaxe1Label.HorizontalAlignment = 'center';
            app.Positionalongaxe1Label.FontAngle = 'italic';
            app.Positionalongaxe1Label.Position = [1 269 256 22];
            app.Positionalongaxe1Label.Text = 'Position along axe 1';

            % Create Positionalongaxe2Label
            app.Positionalongaxe2Label = uilabel(app.AxisoptionsPanel);
            app.Positionalongaxe2Label.HorizontalAlignment = 'center';
            app.Positionalongaxe2Label.FontAngle = 'italic';
            app.Positionalongaxe2Label.Position = [62 189 132 22];
            app.Positionalongaxe2Label.Text = 'Position along axe 2';

            % Create Slider_2
            app.Slider_2 = uislider(app.AxisoptionsPanel);
            app.Slider_2.ValueChangedFcn = createCallbackFcn(app, @Slider_2ValueChanged, true);
            app.Slider_2.ValueChangingFcn = createCallbackFcn(app, @Slider_2ValueChanging, true);
            app.Slider_2.Position = [53 176 151 3];

            % Create Button_m_2
            app.Button_m_2 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_m_2.ButtonPushedFcn = createCallbackFcn(app, @Button_m_2Pushed, true);
            app.Button_m_2.Position = [11 168 28 22];
            app.Button_m_2.Text = '-';

            % Create Button_mm_2
            app.Button_mm_2 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_mm_2.ButtonPushedFcn = createCallbackFcn(app, @Button_mm_2Pushed, true);
            app.Button_mm_2.Position = [11 144 28 22];
            app.Button_mm_2.Text = '--';

            % Create Button_p_2
            app.Button_p_2 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_p_2.ButtonPushedFcn = createCallbackFcn(app, @Button_p_2Pushed, true);
            app.Button_p_2.Position = [218 168 29 22];
            app.Button_p_2.Text = '+';

            % Create Button_pp_2
            app.Button_pp_2 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_pp_2.ButtonPushedFcn = createCallbackFcn(app, @Button_pp_2Pushed, true);
            app.Button_pp_2.Position = [218 144 29 22];
            app.Button_pp_2.Text = '++';

            % Create Positionalongaxe3Label
            app.Positionalongaxe3Label = uilabel(app.AxisoptionsPanel);
            app.Positionalongaxe3Label.HorizontalAlignment = 'center';
            app.Positionalongaxe3Label.FontAngle = 'italic';
            app.Positionalongaxe3Label.Position = [62 109 132 22];
            app.Positionalongaxe3Label.Text = 'Position along axe 3';

            % Create Slider_3
            app.Slider_3 = uislider(app.AxisoptionsPanel);
            app.Slider_3.ValueChangedFcn = createCallbackFcn(app, @Slider_3ValueChanged, true);
            app.Slider_3.ValueChangingFcn = createCallbackFcn(app, @Slider_3ValueChanging, true);
            app.Slider_3.Position = [53 96 151 3];

            % Create Button_m_3
            app.Button_m_3 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_m_3.ButtonPushedFcn = createCallbackFcn(app, @Button_m_3Pushed, true);
            app.Button_m_3.Position = [11 88 28 22];
            app.Button_m_3.Text = '-';

            % Create Button_mm_3
            app.Button_mm_3 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_mm_3.ButtonPushedFcn = createCallbackFcn(app, @Button_mm_3Pushed, true);
            app.Button_mm_3.Position = [11 64 28 22];
            app.Button_mm_3.Text = '--';

            % Create Button_p_3
            app.Button_p_3 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_p_3.ButtonPushedFcn = createCallbackFcn(app, @Button_p_3Pushed, true);
            app.Button_p_3.Position = [219 88 29 22];
            app.Button_p_3.Text = '+';

            % Create Button_pp_3
            app.Button_pp_3 = uibutton(app.AxisoptionsPanel, 'push');
            app.Button_pp_3.ButtonPushedFcn = createCallbackFcn(app, @Button_pp_3Pushed, true);
            app.Button_pp_3.Position = [219 64 29 22];
            app.Button_pp_3.Text = '++';

            % Create AdjustviewasmovingsliderCheckBox
            app.AdjustviewasmovingsliderCheckBox = uicheckbox(app.AxisoptionsPanel);
            app.AdjustviewasmovingsliderCheckBox.Text = 'Adjust view as moving slider';
            app.AdjustviewasmovingsliderCheckBox.Position = [5 25 173 22];

            % Create FigureandvideooptionsPanel
            app.FigureandvideooptionsPanel = uipanel(app.UIFigure);
            app.FigureandvideooptionsPanel.Title = 'Figure and video options';
            app.FigureandvideooptionsPanel.BackgroundColor = [1 1 1];
            app.FigureandvideooptionsPanel.FontWeight = 'bold';
            app.FigureandvideooptionsPanel.Position = [1301 1 256 277];

            % Create ClicktoselectsavefolderButton
            app.ClicktoselectsavefolderButton = uibutton(app.FigureandvideooptionsPanel, 'push');
            app.ClicktoselectsavefolderButton.ButtonPushedFcn = createCallbackFcn(app, @ClicktoselectsavefolderButtonPushed, true);
            app.ClicktoselectsavefolderButton.Position = [10 231 152 22];
            app.ClicktoselectsavefolderButton.Text = 'Click to select save folder';

            % Create FilenameEditFieldLabel
            app.FilenameEditFieldLabel = uilabel(app.FigureandvideooptionsPanel);
            app.FilenameEditFieldLabel.HorizontalAlignment = 'right';
            app.FilenameEditFieldLabel.Position = [11 177 55 22];
            app.FilenameEditFieldLabel.Text = 'Filename';

            % Create FilenameEditField
            app.FilenameEditField = uieditfield(app.FigureandvideooptionsPanel, 'text');
            app.FilenameEditField.Position = [81 177 171 22];

            % Create NosavefolderselectedLabel
            app.NosavefolderselectedLabel = uilabel(app.FigureandvideooptionsPanel);
            app.NosavefolderselectedLabel.Position = [10 204 242 22];
            app.NosavefolderselectedLabel.Text = 'No save folder selected';

            % Create VideoformatDropDownLabel
            app.VideoformatDropDownLabel = uilabel(app.FigureandvideooptionsPanel);
            app.VideoformatDropDownLabel.HorizontalAlignment = 'right';
            app.VideoformatDropDownLabel.Position = [10 145 73 22];
            app.VideoformatDropDownLabel.Text = 'Video format';

            % Create VideoformatDropDown
            app.VideoformatDropDown = uidropdown(app.FigureandvideooptionsPanel);
            app.VideoformatDropDown.Items = {'mpeg-4', 'Uncompressed AVI'};
            app.VideoformatDropDown.Position = [98 145 100 22];
            app.VideoformatDropDown.Value = 'mpeg-4';

            % Create FramepersEditFieldLabel
            app.FramepersEditFieldLabel = uilabel(app.FigureandvideooptionsPanel);
            app.FramepersEditFieldLabel.HorizontalAlignment = 'right';
            app.FramepersEditFieldLabel.Position = [10 118 74 22];
            app.FramepersEditFieldLabel.Text = 'Frame per s.';

            % Create FramepersEditField
            app.FramepersEditField = uieditfield(app.FigureandvideooptionsPanel, 'numeric');
            app.FramepersEditField.Limits = [1 Inf];
            app.FramepersEditField.RoundFractionalValues = 'on';
            app.FramepersEditField.Position = [99 118 38 22];
            app.FramepersEditField.Value = 25;

            % Create Savefigure_button
            app.Savefigure_button = uibutton(app.FigureandvideooptionsPanel, 'push');
            app.Savefigure_button.ButtonPushedFcn = createCallbackFcn(app, @Savefigure_buttonButtonPushed, true);
            app.Savefigure_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Savefigure_button.FontSize = 14;
            app.Savefigure_button.FontWeight = 'bold';
            app.Savefigure_button.FontColor = [1 1 1];
            app.Savefigure_button.Enable = 'off';
            app.Savefigure_button.Position = [64 70 122 33];
            app.Savefigure_button.Text = 'Save figure';

            % Create Savevideo_button
            app.Savevideo_button = uibutton(app.FigureandvideooptionsPanel, 'push');
            app.Savevideo_button.ButtonPushedFcn = createCallbackFcn(app, @Savevideo_buttonButtonPushed, true);
            app.Savevideo_button.BackgroundColor = [0.4667 0.6745 0.1882];
            app.Savevideo_button.FontSize = 14;
            app.Savevideo_button.FontWeight = 'bold';
            app.Savevideo_button.FontColor = [1 1 1];
            app.Savevideo_button.Enable = 'off';
            app.Savevideo_button.Position = [65 32 121 33];
            app.Savevideo_button.Text = 'Save video';

            % Create VolumenameLabel
            app.VolumenameLabel = uilabel(app.UIFigure);
            app.VolumenameLabel.HorizontalAlignment = 'center';
            app.VolumenameLabel.FontSize = 16;
            app.VolumenameLabel.FontWeight = 'bold';
            app.VolumenameLabel.Position = [127 911 1076 22];
            app.VolumenameLabel.Text = 'Volume name';

            % Create UIAxes_1
            app.UIAxes_1 = uiaxes(app.UIFigure);
            title(app.UIAxes_1, 'View normal to first axis')
            xlabel(app.UIAxes_1, 'Third axis x um')
            ylabel(app.UIAxes_1, 'Second axis x um')
            zlabel(app.UIAxes_1, 'Z')
            subtitle(app.UIAxes_1, 'Slice: x / y, a / b')
            app.UIAxes_1.DataAspectRatio = [1 1 1];
            app.UIAxes_1.XColor = [0 0 0];
            app.UIAxes_1.BoxStyle = 'full';
            app.UIAxes_1.Box = 'on';
            app.UIAxes_1.Position = [127 466 503 437];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.UIFigure);
            title(app.UIAxes_2, 'View normal to second axis')
            xlabel(app.UIAxes_2, 'Third axis x um')
            ylabel(app.UIAxes_2, 'First axis x um')
            zlabel(app.UIAxes_2, 'Z')
            subtitle(app.UIAxes_2, 'Slice: x / y, a / b')
            app.UIAxes_2.DataAspectRatio = [1 1 1];
            app.UIAxes_2.XColor = [0 0 0];
            app.UIAxes_2.BoxStyle = 'full';
            app.UIAxes_2.Box = 'on';
            app.UIAxes_2.Position = [700 466 503 437];

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.UIFigure);
            title(app.UIAxes_3, 'View normal to third axis')
            xlabel(app.UIAxes_3, 'First axis x um')
            ylabel(app.UIAxes_3, 'Second axis x um')
            zlabel(app.UIAxes_3, 'Z')
            subtitle(app.UIAxes_3, 'Slice: x / y, a / b')
            app.UIAxes_3.DataAspectRatio = [1 1 1];
            app.UIAxes_3.XColor = [0 0 0];
            app.UIAxes_3.BoxStyle = 'full';
            app.UIAxes_3.Box = 'on';
            app.UIAxes_3.Position = [127 12 503 437];

            % Create UIAxes_3D
            app.UIAxes_3D = uiaxes(app.UIFigure);
            title(app.UIAxes_3D, '3D view')
            xlabel(app.UIAxes_3D, 'First axis x um')
            ylabel(app.UIAxes_3D, 'Second axis x um')
            zlabel(app.UIAxes_3D, 'Z')
            app.UIAxes_3D.DataAspectRatio = [1 1 1];
            app.UIAxes_3D.XColor = [0 0 0];
            app.UIAxes_3D.BoxStyle = 'full';
            app.UIAxes_3D.Box = 'on';
            app.UIAxes_3D.Position = [700 12 503 437];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = microstructure_visualization_slices_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end