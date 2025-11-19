function [ ] = plot_avg_perphase(phasename,phasecolor,avgs,p,optformat,optsave,folder,filename)

%% FIGURE
fig = figure;
fig.Color = 'w';
fig.Name = p.array_name;

nlabel = length(phasename);
nplot = length(p.plots);
tiledlayout(1,nplot,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);

h1 = optsave.Height_foroneplot;
w1 = h1;
optsave.Width = nplot*w1; % + (nplot-1) * WidthBtw;
optsave.Height = h1;

for k=1:nplot
    ax = nexttile;
    if strcmp(char(p.plots(k)),'bar')
        b=bar(ax,phasename,avgs,1.0);
        b.FaceColor = 'flat';
        b.Labels = round(avgs,3);
        for kphase=1:1:nlabel
            b.CData(kphase,:)=phasecolor(kphase,:);
        end  
        xlabel(ax,"Phase")
        ylabel(ax,p.array_name)
        ysecondarylabel(ax,p.array_unit)
        ax.Box = "off";
        ax.XLabel.FontWeight = "bold";
        ax.YLabel.FontWeight = "bold";
        ax.LineWidth = optformat.linewidth;
        ax.FontName = optformat.fontname;
        ax.FontSize = optformat.axefontsize;
        if optformat.grid
            grid(ax,'on');
            set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
        end
        if optformat.includefilenameintitle
            if nplot==1
                title(p.inputfilename,'Interpreter','none','FontSize',optformat.titlefontsize);
            end
        end

    elseif strcmp(char(p.plots(k)),'pie')
        if p.pie_in_percents
            pc=piechart(round(avgs,3)*100,phasename);
        else
            pc=piechart(round(avgs,3),phasename);
        end
        pc.LabelStyle = "namedata";       
        for kphase=1:1:nlabel
            pc.ColorOrder(kphase,:)=phasecolor(kphase,:);
        end
        pc.FontName = optformat.fontname;
        pc.FontSize = optformat.axefontsize;
        if isempty(p.array_unit)
            title(gca,p.array_name);
        else
            title(gca,[p.array_name ' (' p.array_unit ')']);
        end
        if optformat.includefilenameintitle && nplot==1
            pc.Interpreter = 'none';
            title(gca,{"Volume fractions label-wise",infovol.filename});
        end

    end
end

if optformat.includefilenameintitle
    if nplot>1
        sgtitle(p.inputfilename,'Interpreter','none','FontSize',optformat.sgtitlefontsize);
    end
end

%% SAVE FIGURE
if optsave.savefig % Save figure
    function_savefig(fig, folder, filename, optsave); % Call function
end
if optsave.autoclosefig
    close(fig); % Do not keep open figures
end


end