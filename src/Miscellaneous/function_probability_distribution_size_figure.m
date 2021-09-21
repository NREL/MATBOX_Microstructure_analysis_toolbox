function [] = function_probability_distribution_size_figure(data_figure,parameters_figure)

font_axe = 12;
font_axistitle = 14;
font_title = 16;

Fig_= figure; % Create figure
Fig_.Name= parameters_figure.figurename;
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', parameters_figure.figureposition)
sub_axes=subplot(2,2,1,'Parent',Fig_); % Create axes
hold(sub_axes,'on'); % Active subplot
t_=title (' ','FontName',['Raw ' parameters_figure.fontname]); % Title
t_.String= parameters_figure.subaxe1_title;
% Curves
[~, n_curve] = size(data_figure.iteration);
for k_curve=1:1:n_curve
    % Curve
    h_=plot(data_figure.iteration(k_curve).psd.cumulative_fct(:,1),data_figure.iteration(k_curve).psd.cumulative_fct(:,2));
    set(h_,'LineWidth',2); 
    str_legend(k_curve).name = ['Voxel size ' num2str(data_figure.iteration(k_curve).voxelsize*1000,'%1.1f') 'nm' ', D_{50}= ' num2str(data_figure.iteration(k_curve).psd.x50,'%1.3f') ' ' data_figure.unit];
end
% Legend
h_legend = legend(sub_axes,str_legend.name,'Location','best');
% x50
max_x50=0; min_x50=0;
for k_curve=1:1:n_curve
    max_x50 = max([max_x50 data_figure.iteration(k_curve).psd.x50]);
    min_x50 = min([min_x50 min(data_figure.iteration(k_curve).psd.cumulative_fct(:,1))]);
end
h_x50_a=plot([min_x50 max_x50],[0.5 0.5]);
set(h_x50_a,'LineStyle','--','Marker','none','Color','k','LineWidth',1);
set(get(get(h_x50_a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
% for k_curve=1:1:n_curve
%     % x50
%     h_x50_a=plot([min(data_figure.iteration(k_curve).psd.cumulative_fct(:,1)) data_figure.iteration(k_curve).psd.x50],[0.5 0.5]);
%     set(h_x50_a,'LineStyle','--','Marker','none','Color','k','LineWidth',1);
%     set(get(get(h_x50_a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
% end
sub_axes.ColorOrderIndex = 1; % Reset color order
for k_curve=1:1:n_curve
    % x50
    h_x50_b=plot([data_figure.iteration(k_curve).psd.x50 data_figure.iteration(k_curve).psd.x50],[0 0.5]);
    set(h_x50_b,'LineStyle','--','Marker','none','LineWidth',1);
    set(get(get(h_x50_b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
end
% Axis label
xlabel(parameters_figure.xlabel);
ylabel('Cumulative function');
set(sub_axes, 'Ylim', [0 1]) % y axis
grid(sub_axes,parameters_figure.grid); % Display grid
set(sub_axes,'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
set(sub_axes,'FontName', parameters_figure.fontname,'FontSize',font_axe); % Fontname and fontsize
t_.FontSize= font_axistitle;
hold(sub_axes,'off'); % Axe is done

sub_axes=subplot(2,2,2,'Parent',Fig_); % Create axes
hold(sub_axes,'on'); % Active subplot
t_=title (' ','FontName',['Raw ' parameters_figure.fontname]); % Title
t_.String= parameters_figure.subaxe2_title;
% Curves
for k_curve=1:1:n_curve
    % Curve
    h_=plot(data_figure.iteration(k_curve).psd.probability_density_fct(:,1),data_figure.iteration(k_curve).psd.probability_density_fct(:,2));
    set(h_,'LineWidth',2); 
    str_legend(k_curve).name = ['Voxel size ' num2str(data_figure.iteration(k_curve).voxelsize*1000,'%1.1f') 'nm' ', D_{50}= ' num2str(data_figure.iteration(k_curve).psd.x50,'%1.3f') ' ' data_figure.unit ', Integral= ' num2str(data_figure.iteration(k_curve).psd.integral_probability_density_fct,'%1.3f')];
end
% Legend
h_legend = legend(sub_axes,str_legend.name,'Location','best');
sub_axes.ColorOrderIndex = 1; % Reset color order
for k_curve=1:1:n_curve
    % x50
    h_x50=plot([data_figure.iteration(k_curve).psd.x50 data_figure.iteration(k_curve).psd.x50],[0 max(data_figure.iteration(k_curve).psd.probability_density_fct(:,2))]);
    set(h_x50,'LineStyle','--','Marker','none','LineWidth',1);
    set(get(get(h_x50,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
end
% - Axis label
xlabel(parameters_figure.xlabel);
ylabel('Distribution function');
grid(sub_axes,parameters_figure.grid); % Display grid
set(sub_axes,'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
set(sub_axes,'FontName',parameters_figure.fontname,'FontSize',font_axe); % - Fontname and fontsize
t_.FontSize= font_axistitle;
hold(sub_axes,'off'); % Axe is done


exist_cumulative_fct = false;
for k_curve=1:1:n_curve
    if ~isempty(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct)
        exist_cumulative_fct = true;
    end
end
if exist_cumulative_fct    
    sub_axes=subplot(2,2,3,'Parent',Fig_); % Create axes
    hold(sub_axes,'on'); % Active subplot
    t_=title (' ','FontName',parameters_figure.fontname); % Title
    t_.String= ['Smoothed ' parameters_figure.subaxe1_title];
    % Curves
    [~, n_curve] = size(data_figure.iteration);
    k_legend=0;
    clear str_legend
    for k_curve=1:1:n_curve
        % Curve
        if ~isempty(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct)
            h_=plot(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct(:,1),data_figure.iteration(k_curve).psd.smoothed_cumulative_fct(:,2));
            set(h_,'LineWidth',2);
            k_legend=k_legend+1;
            str_legend(k_legend).name = ['Voxel size ' num2str(data_figure.iteration(k_curve).voxelsize*1000,'%1.1f') 'nm' ', D_{50}= ' num2str(data_figure.iteration(k_curve).psd.smoothed_x50,'%1.3f') ' ' data_figure.unit];
        end
    end
    % Legend
    h_legend = legend(sub_axes,str_legend.name,'Location','best');
    % x50
    max_x50=0; min_x50=0;
    for k_curve=1:1:n_curve
        if ~isempty(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct)
            max_x50 = max([max_x50 data_figure.iteration(k_curve).psd.smoothed_x50]);
            min_x50 = min([min_x50 min(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct(:,1))]);
        end
    end
    h_x50_a=plot([min_x50 max_x50],[0.5 0.5]);
    set(h_x50_a,'LineStyle','--','Marker','none','Color','k','LineWidth',1);
    set(get(get(h_x50_a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
    % for k_curve=1:1:n_curve
    %     % x50
    %     h_x50_a=plot([min(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct(:,1)) data_figure.iteration(k_curve).psd.smoothed_x50],[0.5 0.5]);
    %     set(h_x50_a,'LineStyle','--','Marker','none','Color','k','LineWidth',1);
    %     set(get(get(h_x50_a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
    % end
    sub_axes.ColorOrderIndex = 1; % Reset color order
    for k_curve=1:1:n_curve
        % x50
        if ~isempty(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct)
            h_x50_b=plot([data_figure.iteration(k_curve).psd.smoothed_x50 data_figure.iteration(k_curve).psd.smoothed_x50],[0 0.5]);
            set(h_x50_b,'LineStyle','--','Marker','none','LineWidth',1);
            set(get(get(h_x50_b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
        end
    end
    % Axis label
    xlabel(parameters_figure.xlabel);
    ylabel('Cumulative function');
    set(sub_axes, 'Ylim', [0 1]) % y axis
    grid(sub_axes,parameters_figure.grid); % Display grid
    set(sub_axes,'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
    set(sub_axes,'FontName', parameters_figure.fontname,'FontSize',font_axe); % Fontname and fontsize
    t_.FontSize= font_axistitle;
    hold(sub_axes,'off'); % Axe is done
    
    sub_axes=subplot(2,2,4,'Parent',Fig_); % Create axes
    hold(sub_axes,'on'); % Active subplot
    t_=title (' ','FontName', parameters_figure.fontname); % Title
    t_.String= ['Smoothed ' parameters_figure.subaxe2_title];
    % Curves
    clear str_legend
    k_legend=0;
    for k_curve=1:1:n_curve
        % Curve
        if ~isempty(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct)
            h_=plot(data_figure.iteration(k_curve).psd.smoothed_probability_density_fct(:,1),data_figure.iteration(k_curve).psd.smoothed_probability_density_fct(:,2));
            set(h_,'LineWidth',2);
            k_legend=k_legend+1;
            str_legend(k_legend).name = ['Voxel size ' num2str(data_figure.iteration(k_curve).voxelsize*1000,'%1.1f') 'nm' ', D_{50}= ' num2str(data_figure.iteration(k_curve).psd.smoothed_x50,'%1.3f') ' ' data_figure.unit ', Integral= ' num2str(data_figure.iteration(k_curve).psd.integral_smoothed_probability_density_fct,'%1.3f')];
        end
    end
    % Legend
    h_legend = legend(sub_axes,str_legend.name,'Location','best');
    sub_axes.ColorOrderIndex = 1; % Reset color order
    for k_curve=1:1:n_curve
        % x50
        if ~isempty(data_figure.iteration(k_curve).psd.smoothed_cumulative_fct)
            h_x50=plot([data_figure.iteration(k_curve).psd.smoothed_x50 data_figure.iteration(k_curve).psd.smoothed_x50],[0 max(data_figure.iteration(k_curve).psd.smoothed_probability_density_fct(:,2))]);
            set(h_x50,'LineStyle','--','Marker','none','LineWidth',1);
            set(get(get(h_x50,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % No legend
        end
    end
    % - Axis label
    xlabel(parameters_figure.xlabel);
    ylabel('Distribution function');
    grid(sub_axes,parameters_figure.grid); % Display grid
    set(sub_axes,'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
    set(sub_axes,'FontName',parameters_figure.fontname,'FontSize',font_axe); % - Fontname and fontsize
    t_.FontSize= font_axistitle;
    hold(sub_axes,'off'); % Axe is done
    
end

if ~isempty(parameters_figure.title)
    sgtitle(Fig_,parameters_figure.title,'FontWeight','bold','FontSize',font_title,'FontName',parameters_figure.fontname);
end
% Save figures
function_savefig(Fig_, parameters_figure.fullpath, parameters_figure.filename);
close(Fig_)

end

