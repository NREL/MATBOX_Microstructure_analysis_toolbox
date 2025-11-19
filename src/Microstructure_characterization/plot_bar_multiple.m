function [ ] = plot_bar_multiple(phasename,p,optformat,optsave,folder,filename,context)

%% FIGURE
fig = figure;
fig.Color = 'w';
fig.Name = p.name;

nlabel = length(phasename);
nplot = length(p.array_name);
tiledlayout(p.layout(1),p.layout(2),'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);

h1 = optsave.Height_foroneplot;
w1 = h1;
optsave.Width = p.layout(2)*w1;
optsave.Height = p.layout(1)*h1;

for k=1:nplot
    ax = nexttile;
    hold(ax,'on');

    vals = cell2mat(p.vals(k));
    b=bar(ax,phasename,vals,1.0);
    for kb = 1:length(b)
        if length(b)==1
            b(kb).Labels = round(vals(:,kb),3);
        end
        b(kb).FaceColor = 'flat';
    end

    if strcmp(context,'Transport') && k==5 % Bruggeman
        xls = xlim;
        plot([xls(1) xls(2)],[1 1],'Color','k','LineStyle','--','LineWidth',1);
        plot([xls(1) xls(2)],[1.5 1.5],'Color','k','LineStyle','--','LineWidth',1);
        %plot([xls(1) xls(2)],[2.0 2.0],'Color','k','LineStyle','--','LineWidth',1);
        %xt = [xls(1) xls(1) xls(1)];
        %yt = [1 1.5 2.0]+0.1;
        xt = [xls(2) xls(2)];
        yt = [1 1.5];        
        %str = {'Rule of mixture (p=1)','Spheres (p=1.5)','Cylinders (p=2.0)'};
        str = {'Rule of mixture','Spheres'};
        text(xt,yt,str,'FontName',optformat.fontname,'FontSize',optformat.axefontsize);
    end

    xlabel(ax,"Domain")
    ylabel(ax,p.array_name(k))
    ylim(p.yaxis_range(k,:))
    ysecondarylabel(ax,p.array_unit(k))
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

    hl = legend(ax,'String',p.legend);
    hl.FontSize = optformat.legendfontsize;
    hl.IconColumnWidth = 20;
    hl.LineWidth = 1;
    if strcmp(context,'Transport') && k==5
        hl.Location = 'north';
    else
        hl.Location = 'best';
    end

    hold(ax,'off');
end

if optformat.includefilenameintitle
    sgtitle({p.title,p.inputfilename},'Interpreter','none','FontSize',optformat.sgtitlefontsize);
else
    sgtitle(p.title,'FontSize',optformat.sgtitlefontsize)
end

%% SAVE FIGURE
if optsave.savefig % Save figure
    function_savefig(fig, folder, filename, optsave); % Call function
end
if optsave.autoclosefig
    close(fig); % Do not keep open figures
end


end