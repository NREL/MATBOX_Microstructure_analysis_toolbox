function [Fig_] = function_probability_distribution_figure(data_figure,p)

Fig_= figure; % Create figure
Fig_.Name= p.figurename;
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', p.figureposition)
sub_axes{1}=subplot(1,2,1,'Parent',Fig_); % Create axes
hold(sub_axes{1},'on'); % Active subplot
t_=title (' ','FontName',p.fontname); % Title
t_.String= p.subaxe1_title;
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
xlabel(p.xlabel);
ylabel('Cumulative function');
set(sub_axes{1}, 'Ylim', [0 1]) % y axis
% Legend
if ~isempty(data_figure.smoothed_cumulative_fct)
    h_lgd = legend(sub_axes{1},['Raw cumulative function: x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' p.unit],['Smoothed cumulative function: x_{50}= ' num2str(data_figure.smoothed_x50,'%1.2f') ' ' p.unit],'Location','best');
else
    h_lgd = legend(sub_axes{1},['Cumulative function: x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' p.unit],'Location','best');
end
grid(sub_axes{1},p.grid); % Display grid
set(sub_axes{1},'XMinorGrid',p.minorgrid,'YMinorGrid',p.minorgrid); % Display grid for minor thicks also
set(sub_axes{1},'FontName', p.fontname,'FontSize',p.axefontsize); % Fontname and fontsize
t_.FontSize= p.titlefontsize;
h_lgd.FontSize = p.legendfontsize;
hold(sub_axes{1},'off'); % Axe is done

sub_axes{2}=subplot(1,2,2,'Parent',Fig_); % Create axes
hold(sub_axes{2},'on'); % Active subplot
t_=title (' ','FontName',p.fontname); % Title
t_.String= p.subaxe2_title;
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
xlabel(p.xlabel);
ylabel('Distribution function');
% Legend
if ~isempty(data_figure.smoothed_probability_density_fct)
    h_lgd = legend(sub_axes{2},['Distribution function from raw cumulative function' newline 'x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' p.unit ', Integral= ' num2str(data_figure.integral_probability_density_fct,'%1.3f')],['Distribution function from smoothed cumulative function' newline 'x_{50}= ' num2str(data_figure.smoothed_x50,'%1.2f') ' ' p.unit ', Integral= ' num2str(data_figure.integral_smoothed_probability_density_fct,'%1.3f')],'Location','best');
else
    h_lgd = legend(sub_axes{2},['Distribution function' newline 'x_{50}= ' num2str(data_figure.x50,'%1.2f') ' ' p.unit ', Integral= ' num2str(data_figure.integral_probability_density_fct,'%1.3f')],'Location','best');
end
grid(sub_axes{2},p.grid); % Display grid
set(sub_axes{2},'XMinorGrid',p.minorgrid,'YMinorGrid',p.minorgrid); % Display grid for minor thicks also
set(sub_axes{2},'FontName',p.fontname,'FontSize',p.axefontsize); % - Fontname and fontsize
t_.FontSize= p.titlefontsize;
h_lgd.FontSize = p.legendfontsize;
hold(sub_axes{2},'off'); % Axe is done

if ~isempty(p.title)
    sgtitle(Fig_,p.title,'FontWeight','bold','FontSize',p.sgtitlefontsize,'FontName',p.fontname);
end

% Save figures
%p.save=true;
if p.save
    function_savefig(Fig_, p.fullpath, p.filename);
end
%p.closefig = true;
if p.closefig
    close(Fig_)
end

end

