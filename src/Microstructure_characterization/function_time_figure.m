function [Results] = function_time_figure(Time_measure,date_start, date_end, Results, Current_folder, filename, propertyname, OPTIONS)

propertyname(1) = upper(propertyname(1));

%% TIME
lasted_time = date_end-date_start;
% Table
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start_date' 'End_date' 'Lasted_time'});
Table_time_measure = array2table(Time_measure,...
    'VariableNames',{'Voxel_number','CPU_time_s' 'Stopwatch_s'});
% Save
Results.date = Table_date;
Results.time_per_calculation = Table_time_measure;
if OPTIONS.save_xls==true
    % Prepare the data
    clear DATA_writetable
    % Time per calculation
    DATA_writetable.sheet(1).name='Time_per_calculation';
    DATA_writetable.sheet(1).table=Table_time_measure;
    % Data : Date
    DATA_writetable.sheet(2).name='Date';
    DATA_writetable.sheet(2).table=Table_date;
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

% Display
if OPTIONS.displaytext==1
    fprintf ('Finished the %s\n\n',date_end);
    fprintf ('Lasted: %s\n\n',lasted_time);
end

scrsz = get(0,'ScreenSize'); % Screen resolution

% Figure
[a_ , ~]= size(Time_measure);
if a_>1
    Fig = figure; % Create figure
    Fig.Name= 'Time measured per calculation'; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)/2]); % Full screen figure
    for id_axe=1:1:2
        if id_axe==1
            sub_axes=subplot(1,2,id_axe,'Parent',Fig);
            h_title=title ('Linear scale'); % Set title font
        else
            sub_axes=subplot(1,2,id_axe,'Parent',Fig,'XScale', 'log', 'YScale', 'log');
            h_title=title ('Logarithmic  scale'); % Set title font
        end
        hold(sub_axes,'on'); % Active subplot
        % Plot graphs
        h_cpu=plot(Time_measure(:,1),Time_measure(:,2));
        h_stopwatch=plot(Time_measure(:,1),Time_measure(:,3));
        set(h_cpu,'Linestyle','none','Linewidth',1.5,'Marker','o','Markersize',12);
        set(h_stopwatch,'Linestyle','none','Linewidth',1.5,'Marker','d','Markersize',12);
        % Axis label
        xlabel('Number of voxels');
        ylabel('Time (s)');
        % Legend
        h_legend = legend(sub_axes,'CPU time','Stopwatch','Location','best');
        % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
        h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
        sgtitle(Fig,[propertyname ' calculation times'],'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
        hold(sub_axes,'off'); % Relase figure
    end
    
    if OPTIONS.save_fig == true % Save figure
        function_savefig(Fig, Current_folder, filename, OPTIONS); % Call function
    end
    if OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end


end

