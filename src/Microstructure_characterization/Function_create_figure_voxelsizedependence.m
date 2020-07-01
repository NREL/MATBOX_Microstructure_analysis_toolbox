function [] = Function_create_figure_voxelsizedependence(p)

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter
Fig = figure; % Create figure

if isfield(p,'property_voxelsizedependence2') && isfield(p,'property_voxelsizedependence3')
    if isempty(p.method)
        Fig.Name= [p.propertyname ', Voxel size dependence analysis']; % Figure name
        str_title = [p.propertyname ', Voxel size dependence analysis'];
    else
        Fig.Name= [p.propertyname ', ' p.method ', Voxel size dependence analysis']; % Figure name
        str_title = {[p.propertyname ', ' p.method],'Voxel size dependence analysis'};
    end    
    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',scrsz); % Full screen figure
    % - Create axes as a subplot
    for id_axe=1:1:3
        sub_axes=subplot(1,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        h_title=title (p.INFO.direction(id_axe).name); % Title
        if id_axe==1
            pp = p.property_voxelsizedependence;
            pp2 = p.interpolation_voxelsize;
        elseif id_axe==2
            pp = p.property_voxelsizedependence2;
            pp2 = p.interpolation_voxelsize2;
        elseif id_axe==3
            pp = p.property_voxelsizedependence3;
            pp2 = p.interpolation_voxelsize3;
        end
        x_=pp(:,1); % Horizontal axis
        if isempty(pp2)
            for current_phase=1:1:p.number_phase % Loop over phases
                h=plot(x_,pp(:,current_phase+1));
                set(h, 'Color', p.INFO.phase(current_phase).color,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o','LineWidth',p.OPTIONS.Linewidth,'Linestyle','none');
            end
        else
            for current_phase=1:1:p.number_phase % Loop over phases
                h=plot(x_(2:end),pp(2:end,current_phase+1));
                set(h, 'Color', p.INFO.phase(current_phase).color,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o','LineWidth',p.OPTIONS.Linewidth,'Linestyle','none');
            end
            for current_phase=1:1:p.number_phase % Loop over phases
                x_fit=linspace(min(x_),max(x_),100);
                y_fit = polyval(pp2(current_phase).p,x_fit);
                h_extrapolated = plot(x_fit,y_fit);
                %h_extrapolated = plot([min(x_),max(x_)],[polyval(pp2(current_phase).p,min(x_)), polyval(pp2(current_phase).p,max(x_))]);
                set(h_extrapolated, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'Linestyle','--');
            end
        end
        yl=ylim;
        h=plot([p.INFO.voxel_size p.INFO.voxel_size],[yl(1) yl(2)]);
        set(h, 'Color','k','LineWidth',p.OPTIONS.Linewidth,'Linestyle','-.');
        ylim([yl(1) yl(2)]); % Reset y axis limit
        % Axis label
        xlabel('Voxel size (\mum)');
        t_ = ylabel(' ');
        t_.String= p.str_ylabel; % Sprintf does not accept greek characters
        % Legend
        idx_initial = find(pp(:,1) == p.INFO.voxel_size);
        if ~ isempty(pp2)
            idx_0 = find(pp(:,1) == 0);
        end
        for current_phase=1:1:p.number_phase
            if isempty(p.propertynameunit)
                str_legend(current_phase).name = [p.INFO.phase(current_phase).name ', ' num2str(pp(idx_initial,current_phase+1),'%1.3f') ' (' num2str(pp(idx_0,current_phase+1),'%1.3f') ')'];
            else
                str_legend(current_phase).name = [p.INFO.phase(current_phase).name ', ' num2str(pp(idx_initial,current_phase+1),'%1.3f') p.propertynameunit ' (' num2str(pp(idx_0,current_phase+1),'%1.3f') p.propertynameunit  ')'];
            end
        end
        h_legend = legend(sub_axes,str_legend.name,'Location','best');
        % - Grid
        if strcmp(p.OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
        h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,str_title ,'FontWeight','bold','FontSize',p.OPTIONS.Fontsize_title+2,'FontName',p.OPTIONS.fontname);

else

    if isempty(p.method)
        Fig.Name= [p.propertyname ', Voxel size dependence analysis']; % Figure name
        str_title = [p.propertyname ', Voxel size dependence analysis'];
    else
        Fig.Name= [p.propertyname ', ' p.method ', Voxel size dependence analysis']; % Figure name
        str_title = {[p.propertyname ', ' p.method],'Voxel size dependence analysis'};
    end
    Fig.Color='white'; % Background colour
    axe_ = axes('Parent',Fig); % Create axes
    hold(axe_,'on');
    h_title=title (str_title); % Set title
    % Plot graphs
    x_=p.property_voxelsizedependence(:,1); % Horizontal axis
    if isempty(p.interpolation_voxelsize)
        for current_phase=1:1:p.number_phase % Loop over phases
            h=plot(x_,p.property_voxelsizedependence(:,current_phase+1));
            set(h, 'Color', p.INFO.phase(current_phase).color,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o','LineWidth',p.OPTIONS.Linewidth,'Linestyle','none');
        end
    else
        for current_phase=1:1:p.number_phase % Loop over phases
            h=plot(x_(2:end),p.property_voxelsizedependence(2:end,current_phase+1));
            set(h, 'Color', p.INFO.phase(current_phase).color,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o','LineWidth',p.OPTIONS.Linewidth,'Linestyle','none');
        end
        for current_phase=1:1:p.number_phase % Loop over phases
            x_fit=linspace(min(x_),max(x_),100);
            y_fit = polyval(p.interpolation_voxelsize(current_phase).p,x_fit);
            h_extrapolated = plot(x_fit,y_fit);
            %h_extrapolated = plot([min(x_),max(x_)],[polyval(p.interpolation_voxelsize(current_phase).p,min(x_)), polyval(p.interpolation_voxelsize(current_phase).p,max(x_))]);
            set(h_extrapolated, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'Linestyle','--');
        end
    end
    yl=ylim;
    h=plot([p.INFO.voxel_size p.INFO.voxel_size],[yl(1) yl(2)]);
    set(h, 'Color','k','LineWidth',p.OPTIONS.Linewidth,'Linestyle','-.');
    ylim([yl(1) yl(2)]); % Reset y axis limit
    % Axis label
    xlabel('Voxel size (\mum)');
    t_ = ylabel(' ');
    t_.String= p.str_ylabel; % Sprintf does not accept greek characters
    % Legend
    idx_initial = find(p.property_voxelsizedependence(:,1) == p.INFO.voxel_size);
    if ~ isempty(p.interpolation_voxelsize)
        idx_0 = find(p.property_voxelsizedependence(:,1) == 0);
    end
    for current_phase=1:1:p.number_phase
        if isempty(p.propertynameunit)
            str_legend(current_phase).name = [p.INFO.phase(current_phase).name ', ' num2str(p.property_voxelsizedependence(idx_initial,current_phase+1),'%1.3f') ' (' num2str(p.property_voxelsizedependence(idx_0,current_phase+1),'%1.3f') ')'];
        else
            str_legend(current_phase).name = [p.INFO.phase(current_phase).name ', ' num2str(p.property_voxelsizedependence(idx_initial,current_phase+1),'%1.3f') p.propertynameunit ' (' num2str(p.property_voxelsizedependence(idx_0,current_phase+1),'%1.3f') p.propertynameunit  ')'];
        end
    end
    h_legend = legend(axe_,str_legend.name,'Location','best');
    % - Grid
    if strcmp(p.OPTIONS.grid,'on')
        grid(axe_,'on'); % Display grid
        set(axe_,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
    end
    set(axe_,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
    h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
    h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set title fontsize
    hold(axe_,'off'); % Relase figure

end

if p.OPTIONS.save_fig == true % Save figure
    function_savefig(Fig, p.Current_folder, p.filename, p.OPTIONS); % Call function
end
if p.OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end

end

