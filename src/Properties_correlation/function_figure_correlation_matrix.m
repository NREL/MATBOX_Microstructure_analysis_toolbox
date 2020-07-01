function function_figure_correlation_matrix(p)

% Number of volume and parameter
[number_volume, number_paramater] = size(p.matrix);
% Initialize axis (parameter) name
for k=1:1:number_paramater
    if ~isempty(char(p.parameterunitnames(k)))
        parameter_str(k).name = [char(p.parameternames(k)) ' ' char(p.parameterunitnames(k))];
    else
        parameter_str(k).name = char(p.parameternames(k));
    end
    diagonale(k).name = char(p.parameternames(k));
end

for k_vol=1:1:number_volume
    str_vol_legend(k_vol).name = char(p.volumenames(k_vol,1));
end

plot_group = find(cell2mat(p.group(:,end)));
number_group = length(plot_group);
for k_group=1:1:number_group
    current_group = cell2mat (p.group(k_group,1));
    group(k_group).idx = find( cell2mat(p.volumegroup) == current_group );
    str_group_legend(k_group).name = char(p.group(k_group,2));
end

if p.plot_per_volume
    name_ext=' per volume';
else
    name_ext=' per group';
end

%% FIGURE
Fig_ = figure; % Create figure
Fig_.Name= [p.figurename name_ext];
Fig_.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig_,'position',scrsz); % Full screen figure
id_axe=0;
for row=1:1:number_paramater
    for column=1:1:number_paramater
        id_axe=id_axe+1; % Update axis selection
        
        sub_axes=subplot(number_paramater,number_paramater,id_axe,'Parent',Fig_); % Create new axis
        hold(sub_axes,'on'); % Active subplot
        
        if row>column % Plot only the lower triangle
            if p.plot_per_volume
                for k_vol=1:1:number_volume
                    x=p.matrix(k_vol,column);
                    y=p.matrix(k_vol,row);
                    plot(x,y,'MarkerSize',p.MarkerSize,'Marker',char(p.marker_order(k_vol,1)),'LineWidth',p.Linewidth,'Linestyle','none');
                end
                h_legend = legend(sub_axes,str_vol_legend.name,'Location','best');
            
            elseif p.plot_per_group
                for k_group=1:1:number_group
                    x=p.matrix(group(k_group).idx,column);
                    y=p.matrix(group(k_group).idx,row);
                    color_group = cell2mat(p.group(k_group,5:7));
                    if logical(cell2mat(p.group(k_group,4)))
                        plot(x,y,'MarkerSize',p.MarkerSize,'color',color_group,'Marker',char(p.group(k_group,3)), 'MarkerFaceColor',color_group,'LineWidth',p.Linewidth,'Linestyle','none');
                    else
                        plot(x,y,'MarkerSize',p.MarkerSize,'color',color_group,'Marker',char(p.group(k_group,3)),'LineWidth',p.Linewidth,'Linestyle','none');
                    end
                end
                h_legend = legend(sub_axes,str_group_legend.name,'Location','best');
            end
            xlabel(parameter_str(column).name); % X Axis label
            ylabel(parameter_str(row).name); % Y Axis label            
            % - Grid
            if strcmp(p.grid,'on')
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid',p.minorgrid,'YMinorGrid',p.minorgrid); % Display grid for minor thicks
            end
            set(sub_axes,'FontName',p.fontname,'FontSize',p.Fontsize_axe); % Fontname and fontsize
            h_legend.FontSize = p.Fontsize_legend; % Set title fontsize
        end
        
        % Parameter name on the diagonale
        if row==column
            % Remove the tick
            set(sub_axes,'xtick',[])
            set(sub_axes,'ytick',[])
            % White axe
            set(sub_axes,'Xcolor','w')
            set(sub_axes,'Ycolor','w')
            % - Box
            box(sub_axes,'off'); % Display box
            % - Grid
            grid(sub_axes,'off'); % Remove grid
            % Text
            xt = 0.5; yt = 0.5;
            text(xt,yt,diagonale(row).name,'Color','k','FontSize',p.font_diagonale,'FontName',p.fontname,'HorizontalAlignment','center');
        end
        
        % Correlation on the upper triangle
        if row<column
            % Remove the tick
            set(sub_axes,'xtick',[])
            set(sub_axes,'ytick',[])
            % - Box
            box(sub_axes,'on'); % Display box
            % - Grid
            grid(sub_axes,'off'); % Remove grid
            % Text
            xt = 0.5; yt = 0.5;
            if p.plot_per_volume || number_group==1
                str_={'Kendal''s Tau ranking',num2str(p.correlation(row,column),'%1.3f')};
            else
                str_={'Kendal''s Tau ranking',[num2str(p.correlation(row,column),'%1.3f') ' (' num2str(p.correlation_group(row,column),'%1.3f') ')']};
            end
            text(xt,yt,str_,'Color','k','FontSize',p.font_diagonale,'FontName',p.fontname,'HorizontalAlignment','center');
        end
        
        hold(sub_axes,'off'); % Relase figure
    end
end
sgtitle(Fig_,p.figurename,'FontWeight','bold','FontSize',p.Fontsize_title,'FontName',p.fontname);
if p.save_fig == true % Save figure
    filename= function_remove_emptyandspecialcharacter_string([p.figurename name_ext]);
    function_savefig(Fig_, p.Save_folder, filename); % Call function
end

%% CORRELATION HEAT MAP

% Correlation map figure
Fig_correlationmap = figure;
Fig_correlationmap.Name= ['Correlation map: ' p.figurename];
Fig_correlationmap.Color='white'; % Background colour
% Generate axis labels
clear xvalues yvalues yvalues2
for k=1:1:number_paramater
    xvalues(k)={parameter_str(k).name};
    yvalues(k)={parameter_str(k).name};
    yvalues2(k)={num2str(k)};
end
if p.plot_per_volume || number_group==1
    n_axe=2; row_axe=1;
    set(Fig_correlationmap,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)/2]); % Full screen figure
else
    n_axe=4; row_axe=2;
    set(Fig_correlationmap,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)]); % Full screen figure
end
for id_axe=1:1:n_axe
    % Create new axis
    sub_axes=subplot(row_axe,2,id_axe,'Parent',Fig_correlationmap);
    % Relative or absolute values
    if id_axe==1
        cdata=(p.correlation);
    elseif id_axe==2
        cdata=abs(p.correlation);
    elseif id_axe==3
        cdata=(p.correlation_group);
    elseif id_axe==4
        cdata=abs(p.correlation_group);        
    end
    % Plot heat map
    h=heatmap(xvalues,yvalues,cdata);
    % Color limit
    if id_axe==1 || id_axe==3
        h.ColorLimits = [-1 1];
    else
        h.ColorLimits = [0 1];
    end
    % Colormap
    % h.Colormap = jet;
    % Title
    if id_axe==1
        h.Title = 'Kendal''s Tau ranking';
    elseif id_axe==2
        h.Title = 'Kendal''s Tau ranking (absolute)';
    elseif id_axe==3
        h.Title = 'Kendal''s Tau ranking, group average';
    elseif id_axe==4
        h.Title = 'Kendal''s Tau ranking, group average (absolute)';        
    end    
end
if p.save_fig == true % Save figure
    filename= function_remove_emptyandspecialcharacter_string(['Correlation map ' p.figurename name_ext]);
    function_savefig(Fig_correlationmap, p.Save_folder, filename); % Call function
end

end

