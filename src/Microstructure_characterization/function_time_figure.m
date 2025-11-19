function [] = function_time_figure(timedata_pervolume, timedata_perphase, folder, filename, propertyname, input_filename, optformat, optsave)

propertyname(1) = upper(propertyname(1));
[n_volume , ~] = size(timedata_pervolume);
[n_phase, ~] = size(timedata_perphase);

if n_volume>1
    Fig = figure;
    Fig.Name= 'Time measured per volume and per phase';
    Fig.Color='white';

    if n_phase==1
        n_row = 1; n_axe = 2;
    else
        n_row = 2; n_axe = 4;
    end    
    tiledlayout(n_row,2,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);

    h1 = optsave.Height_foroneplot;
    w1 = h1;
    optsave.Height = n_row*h1;
    optsave.Width = 2*w1;

    for id_axe=1:1:n_axe
        ax = nexttile;
        hold(ax,'on');

        if id_axe==1
            h_title=title ('Linear scale (per volume)'); % Set title font
            tmp = timedata_pervolume;
        elseif id_axe==2
            h_title=title ('Logarithmic  scale (per volume)'); % Set title font
            tmp = timedata_pervolume;
        elseif id_axe==3
            h_title=title ('Linear scale (per phase)'); % Set title font
            tmp = timedata_perphase;
        elseif id_axe==4
            h_title=title ('Logarithmic  scale (per phase)'); % Set title font            
            tmp = timedata_perphase;
        end

        h_cpu=plot(tmp(:,1),tmp(:,2));
        h_stopwatch=plot(tmp(:,1),tmp(:,3));
        set(h_cpu,'Linestyle','none','Linewidth',0.5,'Marker','o','Markersize',8,'DisplayName','CPU time');
        set(h_stopwatch,'Linestyle','none','Linewidth',0.5,'Marker','d','Markersize',8,'DisplayName','Stopwatch');

        if id_axe==1 || id_axe==3
            set(ax,'XScale','linear','YScale','linear');
        else
            set(ax,'XScale','log','YScale','log');
        end
        xlabel(ax,'Number of voxels');
        ylabel(ax,'Time (s)');
        ax.XLabel.FontWeight = "bold";
        ax.YLabel.FontWeight = "bold";
        ax.LineWidth = optformat.linewidth;
        ax.FontName = optformat.fontname;
        ax.FontSize = optformat.axefontsize;

        hl = legend;
        hl.FontSize = optformat.legendfontsize;
        hl.Location = 'best';
        hl.IconColumnWidth = 20;
        hl.LineWidth = 1;

        if optformat.grid
            grid(ax,'on');
            set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
        end

        h_title.FontSize = optformat.titlefontsize; % Set title fontsize
        sgtitle(Fig,[propertyname ' calculation times'],'FontWeight','bold','FontSize',optformat.sgtitlefontsize,'FontName',optformat.fontname);

        if optformat.includefilenameintitle
            sgtitle({[propertyname ' calculation times'],input_filename},'Interpreter','none','FontSize',optformat.sgtitlefontsize);
        else
            sgtitle(Fig,[propertyname ' calculation times'],'Interpreter','none','FontSize',optformat.sgtitlefontsize);
        end

        hold(ax,'off'); % Relase figure
    end    

    if optsave.savefig % Save figure
        function_savefig(Fig, folder, filename, optsave); % Call function
    end
    if optsave.autoclosefig
        close(Fig); % Do not keep open figures
    end

    
end

end
