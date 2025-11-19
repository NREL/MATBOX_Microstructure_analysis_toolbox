function [] = ImageResolution_erroranalysis_figures(p)

number_domain = length(p.domain_name);

optsave = p.opt.save;
optformat = p.opt.format;

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter

% X-axis unit
strunit = p.infovol.unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '\mum';
else
    axisunit = strunit;
end
if strcmp(p.propertynameunit,'um') || strcmp(p.propertynameunit,'micrometer') || strcmp(p.propertynameunit,'Micrometer') || strcmp(p.propertynameunit,'micrometers') || strcmp(p.propertynameunit,'Micrometers')
    p.propertynameunit = '\mum';
end

fig = figure;
fig.Color = 'w';
tiledlayout(1,2,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);
  
h1 = optsave.Height_foroneplot;
w1 = h1;
optsave.Height = h1;
optsave.Width = 2*w1;

fig.Name= [p.propertyname ', voxel size dependence analysis'];

for id_axe=1:1:2
    ax = nexttile;
    hold(ax,'on');

    if id_axe==1
        h_title=title ('Linear scale');
    else
        set(ax,'XScale','log','YScale', 'log')
        h_title=title ('Logarithmic  scale');
    end

    % Plot graphs
    x_=p.property_voxelsizedependence(:,1); % Horizontal axis
    if isempty(p.interpolation_voxelsize)
        for current_domain=1:number_domain % Loop over all phases
            h=plot(x_,p.property_voxelsizedependence(:,current_domain+1));
            set(h, 'Color', p.domain_color(current_domain,:),'MarkerSize',optformat.axefontsize,'Marker','o','LineWidth',optformat.linewidth,'Linestyle','none');
        end
    else
        for current_domain=1:number_domain % Loop over all phases
            h=plot(x_(id_axe:end),p.property_voxelsizedependence(id_axe:end,current_domain+1));
            set(h, 'Color', p.domain_color(current_domain,:),'MarkerSize',optformat.axefontsize,'Marker','o','LineWidth',optformat.linewidth,'Linestyle','none');
        end
        for current_domain=1:1:number_domain % Loop over all phases
            if id_axe==1
                x_fit=linspace(min(x_),max(x_),100);
                y_fit = polyval(p.interpolation_voxelsize.domain(current_domain).pi,x_fit);
            else
                % Linear scale always in log-log
                x=x_(2:end);
                y=p.property_voxelsizedependence(2:end,current_domain+1);
                pp=polyfit(log(x),log(y),1);
                x_fit = linspace(min(x),max(x),100);
                y_fit = x_fit.^(pp(1)).*exp(pp(2));
            end
            h_extrapolated = plot(x_fit,y_fit);
            set(h_extrapolated, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'Linestyle','--');
        end
    end

    yl=ylim;
    h=plot([p.infovol.voxelsize p.infovol.voxelsize],[yl(1) yl(2)]);
    set(h, 'Color','k','LineWidth',optformat.linewidth,'Linestyle','-.');
    ylim([yl(1) yl(2)]); % Reset y axis limit

    xlabel('Voxel size');
    xsecondarylabel(ax,axisunit)
    ylabel(p.ylabel);
    if ~isempty(p.propertynameunit)
        ysecondarylabel(ax,p.propertynameunit)
    end

    % Legend
    idx_initial = find(p.property_voxelsizedependence(:,1) == p.infovol.voxelsize);
    if ~isempty(p.interpolation_voxelsize)
        idx_0 = find(p.property_voxelsizedependence(:,1) == 0);
    end
    for current_domain=1:1:number_domain
        if ~isempty(p.interpolation_voxelsize)
            if isempty(p.propertynameunit)
                str_legend(current_domain).name = [char(p.domain_name(current_domain)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain+1),'%1.3f') ' (' num2str(p.property_voxelsizedependence(idx_0,current_domain+1),'%1.3f') ')'];
            else
                str_legend(current_domain).name = [char(p.domain_name(current_domain)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain+1),'%1.3f') p.propertynameunit ' (' num2str(p.property_voxelsizedependence(idx_0,current_domain+1),'%1.3f') p.propertynameunit  ')'];
            end
        else
            if isempty(p.propertynameunit)
                str_legend(current_domain).name = [char(p.domain_name(current_domain)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain+1),'%1.3f')];
            else
                str_legend(current_domain).name = [char(p.domain_name(current_domain)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain+1),'%1.3f') p.propertynameunit];
            end
        end
    end

    hl = legend(ax,str_legend.name,'Location','best');
    hl.FontSize = optformat.legendfontsize;
    hl.Location = 'best';
    hl.IconColumnWidth = 20;
    hl.LineWidth = 1;

    ax.XLabel.FontWeight = "bold";
    ax.YLabel.FontWeight = "bold";
    ax.LineWidth = optformat.linewidth;
    ax.FontName = optformat.fontname;
    ax.FontSize = optformat.axefontsize;
    if optformat.grid
        grid(ax,'on');
        set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
    end
    hold(ax,'off'); % Relase figure
end

if optformat.includefilenameintitle
    sgtitle(p.inputfilename,'Interpreter','none','FontSize',optformat.sgtitlefontsize);
end

if optsave.savefig % Save figure
    function_savefig(fig, p.Current_folder, p.filename, optsave); % Call function
end
if optsave.autoclosefig
    close(fig); % Do not keep open figures
end

end

