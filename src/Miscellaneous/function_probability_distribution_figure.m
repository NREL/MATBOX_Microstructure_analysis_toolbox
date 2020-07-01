function [Fig_] = function_probability_distribution_figure(data_figure,parameters_figure)

font_axe = 12;
font_axistitle = 14;
font_title = 16;

Fig_= figure; % Create figure
Fig_.Name= parameters_figure.figurename;
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', parameters_figure.figureposition)
sub_axes{1}=subplot(1,2,1,'Parent',Fig_); % Create axes
hold(sub_axes{1},'on'); % Active subplot
t_=title (' ','FontName',parameters_figure.fontname); % Title
t_.String= parameters_figure.subaxe1_title;
% Curves
h_raw=plot(data_figure.cumulative_fct(:,1),data_figure.cumulative_fct(:,2));
if ~isempty(data_figure.smoothed_cumulative_fct)
    h_smooth=plot(data_figure.smoothed_cumulative_fct(:,1),data_figure.smoothed_cumulative_fct(:,2));
end
% Colors
if ~isempty(data_figure.smoothed_cumulative_fct)
   set(h_raw,'LineStyle','-','Marker','none','Color','k','LineWidth',1); 
   set(h_smooth,'LineStyle','-','Marker','none','Color','r','LineWidth',2);
else
   set(h_raw,'LineStyle','-','Marker','none','Color','r','LineWidth',2);
end
% x50
h_raw_x50_a=plot([min(data_figure.cumulative_fct(:,1)) data_figure.x50],[0.5 0.5]);
h_raw_x50_b=plot([data_figure.x50 data_figure.x50],[0 0.5]);
set(h_raw_x50_a,'LineStyle','--','Marker','none','Color','k','LineWidth',1); 
set(h_raw_x50_b,'LineStyle','--','Marker','none','Color','k','LineWidth',1); 
if ~isempty(data_figure.smoothed_cumulative_fct)
    h_smooth_x50_a=plot([min(data_figure.smoothed_cumulative_fct(:,1)) data_figure.smoothed_x50],[0.5 0.5]);
    h_smooth_x50_b=plot([data_figure.smoothed_x50 data_figure.smoothed_x50],[0 0.5]);
    set(h_smooth_x50_a,'LineStyle','--','Marker','none','Color','r','LineWidth',1);
    set(h_smooth_x50_b,'LineStyle','--','Marker','none','Color','r','LineWidth',1);
end
% Axis label
xlabel(parameters_figure.xlabel);
ylabel('Cumulative function');
set(sub_axes{1}, 'Ylim', [0 1]) % y axis
% Legend
if ~isempty(data_figure.smoothed_cumulative_fct)
    legend(sub_axes{1},['Raw cumulative function: x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' data_figure.unit],['Smoothed cumulative function: x_{50}= ' num2str(data_figure.smoothed_x50,'%1.2f') ' ' data_figure.unit],'Location','best');
else
    legend(sub_axes{1},['Cumulative function: x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' data_figure.unit],'Location','best');
end
grid(sub_axes{1},parameters_figure.grid); % Display grid
set(sub_axes{1},'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
set(sub_axes{1},'FontName', parameters_figure.fontname,'FontSize',font_axe); % Fontname and fontsize
t_.FontSize= font_axistitle;
hold(sub_axes{1},'off'); % Axe is done

sub_axes{2}=subplot(1,2,2,'Parent',Fig_); % Create axes
hold(sub_axes{2},'on'); % Active subplot
t_=title (' ','FontName',parameters_figure.fontname); % Title
t_.String= parameters_figure.subaxe2_title;
h_raw=plot(data_figure.probability_density_fct(:,1),data_figure.probability_density_fct(:,2)); % Curves
if ~isempty(data_figure.smoothed_probability_density_fct)
    h_smooth=plot(data_figure.smoothed_probability_density_fct(:,1),data_figure.smoothed_probability_density_fct(:,2));
end
% Colors
if ~isempty(data_figure.smoothed_probability_density_fct)
   set(h_raw,'LineStyle','-','Marker','none','Color','k','LineWidth',1); 
   set(h_smooth,'LineStyle','-','Marker','none','Color','r','LineWidth',2);
else
   set(h_raw,'LineStyle','-','Marker','none','Color','r','LineWidth',2);
end
% x50
h_raw_x50=plot([data_figure.x50 data_figure.x50],[0 max(data_figure.probability_density_fct(:,2))]);
set(h_raw_x50,'LineStyle','--','Marker','none','Color','k','LineWidth',1); 
if ~isempty(data_figure.smoothed_probability_density_fct)
    h_smooth_x50=plot([data_figure.smoothed_x50 data_figure.smoothed_x50],[0 max(data_figure.smoothed_probability_density_fct(:,2))]);
    set(h_smooth_x50,'LineStyle','--','Marker','none','Color','r','LineWidth',1);
end
% - Axis label
xlabel(parameters_figure.xlabel);
ylabel('Distribution function');
% Legend
if ~isempty(data_figure.smoothed_probability_density_fct)
    legend(sub_axes{2},['Distribution function from raw cumulative function' newline 'x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' data_figure.unit ', Integral= ' num2str(data_figure.integral_probability_density_fct,'%1.3f')],['Distribution function from smoothed cumulative function' newline 'x_{50}= ' num2str(data_figure.smoothed_x50,'%1.2f') ' ' data_figure.unit ', Integral= ' num2str(data_figure.integral_smoothed_probability_density_fct,'%1.3f')],'Location','best');
else
    legend(sub_axes{2},['Distribution function' newline 'x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' data_figure.unit ', Integral= ' num2str(data_figure.integral_probability_density_fct,'%1.3f')],'Location','best');
end
grid(sub_axes{2},parameters_figure.grid); % Display grid
set(sub_axes{2},'XMinorGrid',parameters_figure.minorgrid,'YMinorGrid',parameters_figure.minorgrid); % Display grid for minor thicks also
set(sub_axes{2},'FontName',parameters_figure.fontname,'FontSize',font_axe); % - Fontname and fontsize
t_.FontSize= font_axistitle;
hold(sub_axes{2},'off'); % Axe is done

if ~isempty(parameters_figure.title)
    sgtitle(Fig_,parameters_figure.title,'FontWeight','bold','FontSize',font_title,'FontName',parameters_figure.fontname);
end

% Save figures
parameters_figure.save=true;
if parameters_figure.save
    function_savefig(Fig_, parameters_figure.fullpath, parameters_figure.filename);
    close(Fig_)
end

end

