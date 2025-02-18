function [] = ImageResolution_erroranalysis_figures(p)

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter
% Number of phase
number_domain = length(p.todo);

% X-axis unit
strunit = p.infovol.unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
else
    axisunit = ['(' strunit ')'];
end

Fig = figure; % Create figure
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)/2]);
Fig.Color='white'; % Background colour

if isempty(p.method)
    Fig.Name= [p.propertyname ', voxel size dependence analysis']; % Figure name
    sg_title = [p.propertyname ', voxel size dependence analysis'];
else
    Fig.Name= [p.propertyname ', ' p.method ', voxel size dependence analysis']; % Figure name
    sg_title = {[p.propertyname ', ' p.method],'voxel size dependence analysis'};
end

for id_axe=1:1:2
    if id_axe==1
        sub_axes=subplot(1,2,id_axe,'Parent',Fig);
        h_title=title ('Linear scale');
    else
        sub_axes=subplot(1,2,id_axe,'Parent',Fig,'XScale', 'log', 'YScale', 'log');
        h_title=title ('Logarithmic  scale');
    end
    hold(sub_axes,'on');
    % Plot graphs
    x_=p.property_voxelsizedependence(:,1); % Horizontal axis
    if isempty(p.interpolation_voxelsize)
        current_domain_todo = 0;
        for current_domain=1:1:number_domain % Loop over all phases
            if p.todo(current_domain)
                current_domain_todo=current_domain_todo+1;
                h=plot(x_,p.property_voxelsizedependence(:,current_domain_todo+1));
                set(h, 'Color', p.infovol.phasecolor(current_domain,:),'MarkerSize',p.opt.format.axefontsize,'Marker','o','LineWidth',p.opt.format.linewidth,'Linestyle','none');
            end
        end
    else
        current_domain_todo = 0;
        for current_domain=1:1:number_domain % Loop over all phases
            if p.todo(current_domain)
                current_domain_todo=current_domain_todo+1;
                h=plot(x_(id_axe:end),p.property_voxelsizedependence(id_axe:end,current_domain_todo+1));
                set(h, 'Color', p.infovol.phasecolor(current_domain,:),'MarkerSize',p.opt.format.axefontsize,'Marker','o','LineWidth',p.opt.format.linewidth,'Linestyle','none');
            end
        end
        current_domain_todo = 0;
        for current_domain=1:1:number_domain % Loop over all phases
            if p.todo(current_domain)
                current_domain_todo=current_domain_todo+1;
                if id_axe==1
                    x_fit=linspace(min(x_),max(x_),100);
                    y_fit = polyval(p.interpolation_voxelsize.domain(current_domain_todo).pi,x_fit);
                else
                    % Linear scale always in log-log
                    x=x_(2:end);
                    y=p.property_voxelsizedependence(2:end,current_domain_todo+1);
                    pp=polyfit(log(x),log(y),1);
                    x_fit = linspace(min(x),max(x),100);
                    y_fit = x_fit.^(pp(1)).*exp(pp(2));
                end
                h_extrapolated = plot(x_fit,y_fit);
                set(h_extrapolated, 'Color', p.infovol.phasecolor(current_domain,:),'LineWidth',p.opt.format.linewidth,'Linestyle','--');
            end
        end
    end
    yl=ylim;
    h=plot([p.infovol.voxelsize p.infovol.voxelsize],[yl(1) yl(2)]);
    set(h, 'Color','k','LineWidth',p.opt.format.linewidth,'Linestyle','-.');
    ylim([yl(1) yl(2)]); % Reset y axis limit
    % Axis label
    xlabel(['Voxel size ' axisunit]);
    t_ = ylabel(' ');
    t_.String= p.str_ylabel; % Sprintf does not accept greek characters
    % Legend
    idx_initial = find(p.property_voxelsizedependence(:,1) == p.infovol.voxelsize);
    if ~isempty(p.interpolation_voxelsize)
        idx_0 = find(p.property_voxelsizedependence(:,1) == 0);
    end
    current_domain_todo = 0;
    for current_domain=1:1:number_domain
        if p.todo(current_domain)
            current_domain_todo=current_domain_todo+1;
            if ~isempty(p.interpolation_voxelsize)
                if isempty(p.propertynameunit)
                    str_legend(current_domain_todo).name = [char(p.infovol.phasename(current_domain,1)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain_todo+1),'%1.3f') ' (' num2str(p.property_voxelsizedependence(idx_0,current_domain_todo+1),'%1.3f') ')'];
                else
                    str_legend(current_domain_todo).name = [char(p.infovol.phasename(current_domain,1)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain_todo+1),'%1.3f') p.propertynameunit ' (' num2str(p.property_voxelsizedependence(idx_0,current_domain_todo+1),'%1.3f') p.propertynameunit  ')'];
                end
            else
                if isempty(p.propertynameunit)
                    str_legend(current_domain_todo).name = [char(p.infovol.phasename(current_domain,1)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain_todo+1),'%1.3f')];
                else
                    str_legend(current_domain_todo).name = [char(p.infovol.phasename(current_domain,1)) ', ' num2str(p.property_voxelsizedependence(idx_initial,current_domain_todo+1),'%1.3f') p.propertynameunit];
                end
            end
        end
    end

    h_legend = legend(sub_axes,str_legend.name,'Location','best');
    % - Grid
    grid(sub_axes,p.opt.format.grid); % Display grid
    set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
    set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
    h_legend.FontSize = p.opt.format.legendfontsize; % Set title fontsize
    hold(sub_axes,'off'); % Relase figure
end
sgtitle(Fig,sg_title,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);


if p.opt.save.savefig % Save figure
    function_savefig(Fig, p.Current_folder, p.filename, p.opt.save); % Call function
end
if p.opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

end

