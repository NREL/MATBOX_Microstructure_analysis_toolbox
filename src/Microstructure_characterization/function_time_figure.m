function [] = function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, filename, propertyname, opt)

propertyname(1) = upper(propertyname(1));
scrsz = get(0,'ScreenSize'); % Screen resolution
[n_volume , ~] = size(timedata_pervolume);
[n_phase, ~] = size(timedata_perphase);
if n_volume>1
    Fig = figure; % Create figure
    Fig.Name= 'Time measured per volume and per phase'; % Figure name
    Fig.Color='white'; % Background colour
    if n_phase==1
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)/2]); % Full screen figure
        n_axe = 2;
    else
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)]); % Full screen figure
        n_axe = 4;
    end
    for id_axe=1:1:n_axe
        if id_axe==1
            sub_axes=subplot(n_axe/2,2,id_axe,'Parent',Fig);
            h_title=title ('Calculation time per volume, linear scale'); % Set title font
            tmp = timedata_pervolume;
        elseif id_axe==2
            sub_axes=subplot(n_axe/2,2,id_axe,'Parent',Fig,'XScale', 'log', 'YScale', 'log');
            h_title=title ('Calculation time per volume, logarithmic  scale'); % Set title font
            tmp = timedata_pervolume;
        elseif id_axe==3
            sub_axes=subplot(n_axe/2,2,id_axe,'Parent',Fig);
            h_title=title ('Calculation time per phase, linear scale'); % Set title font
            tmp = timedata_perphase;
        elseif id_axe==4
            sub_axes=subplot(n_axe/2,2,id_axe,'Parent',Fig,'XScale', 'log', 'YScale', 'log');
            h_title=title ('Calculation time per phase, logarithmic  scale'); % Set title font            
            tmp = timedata_perphase;
        end
        hold(sub_axes,'on'); % Active subplot
        % Plot graphs
        h_cpu=plot(tmp(:,1),tmp(:,2));
        h_stopwatch=plot(tmp(:,1),tmp(:,3));
        set(h_cpu,'Linestyle','none','Linewidth',0.5,'Marker','o','Markersize',8);
        set(h_stopwatch,'Linestyle','none','Linewidth',0.5,'Marker','d','Markersize',8);
        % Axis label
        xlabel('Number of voxels');
        ylabel('Time (s)');
        % Legend
        h_legend = legend(sub_axes,'CPU time','Stopwatch','Location','best');
        % - Grid
        grid(sub_axes,'on'); % Display grid
        set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks
        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
        h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
        sgtitle(Fig,[propertyname ' calculation times'],'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
        hold(sub_axes,'off'); % Relase figure
    end    
    if opt.save.savefig % Save figure
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end
end

end
