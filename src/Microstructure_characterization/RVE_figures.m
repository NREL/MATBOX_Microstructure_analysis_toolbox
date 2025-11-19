function [] = RVE_figures(p)

number_domain = length(p.domain_name);

optsave = p.opt.save;
optformat = p.opt.format;

h1 = optsave.Height_foroneplot;
w1 = h1;

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter
strunit =  p.infovol.unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '\mum'; lgdunit = '\mum';
else
    axisunit = strunit; lgdunit = strunit;
end
if strcmp(p.propertynameunit,'um') || strcmp(p.propertynameunit,'micrometer') || strcmp(p.propertynameunit,'Micrometer') || strcmp(p.propertynameunit,'micrometers') || strcmp(p.propertynameunit,'Micrometers')
    p.propertynameunit = '\mum';
end

idx_size = find(p.Property_subdomains_statistics(1,1:3,1) > 0);
idxcloud_size = find(p.Property_eachsubdomain(1,2:4,1) > 0)+1;

if strcmp(p.RVE.analysis,'Independent subvolumes')
    if length(p.Wholevolume_size)==2
        r_axe=2; n_axe=8;
        n2_axe=2;
    else
        r_axe=1; n_axe=4;
        n2_axe=1;
    end
    optsave.Height = r_axe*h1;
    optsave.Width = 4*w1;

    if strcmp(p.RVE.threshold_subs_choice,'Relative standard deviation')
        idx_threshold = 7;
        str_threshold = 'Relative standard deviation';
    elseif strcmp(p.RVE.threshold_subs_choice,'Standard deviation')
        idx_threshold = 6;
        str_threshold = 'Standard deviation';
    elseif strcmp(p.RVE.threshold_subs_choice,'Maximum - minimum')
        idx_threshold = 10;
        str_threshold = 'Maximum - minimum';
    end

    for current_domain=1:number_domain % Loop over all phases
        Fig = figure; % Create figure
        Fig.Name= ['Representativity analysis: ' p.propertyname ', ' char(p.domain_name(current_domain))];
        Fig.Color='white'; % Background colour
        tiledlayout(r_axe,4,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);

        for id_axe=1:1:n_axe % Iterate over axe
            ax = nexttile;
            hold(ax,'on');

            if id_axe<=4
                str_xlabel=['Subdomain FOV ' char(p.Length_name(1))];
                kx = idx_size(1);
                kxcloud = idxcloud_size(1);
                ksizeRVE = 1;
                kx_whole = 1;
            else
                str_xlabel=['Subdomain FOV ' char(p.Length_name(2))];
                kx = idx_size(2);
                kxcloud = idxcloud_size(2);
                ksizeRVE = 2;
                kx_whole = 2;
            end

            h_title=title(p.propertyname);
            x_=p.Property_subdomains_statistics(:,kx,current_domain);

            if id_axe==2 ||  id_axe==6
                h_subtitle = subtitle('Statistics');
                

                % Mean with error bar (+- standard deviation)
                x_tmp = []; y_tmp = []; e_tmp = [];
                for current_size=1:1:length(x_)
                    x_tmp=[x_tmp p.Property_subdomains_statistics(current_size,kx,current_domain)];
                    y_tmp=[y_tmp p.Property_subdomains_statistics(current_size,5,current_domain)];
                    e_tmp=[e_tmp p.Property_subdomains_statistics(current_size,6,current_domain)];
                end
                h_mean_witherrorbar = errorbar(x_tmp,y_tmp,e_tmp);
                set(h_mean_witherrorbar, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o');

                % 5% interval confidence
                int_min = (p.Property_subdomains_statistics(:,5,current_domain) - p.Property_subdomains_statistics(:,15,current_domain))';
                int_max = (p.Property_subdomains_statistics(:,5,current_domain) + p.Property_subdomains_statistics(:,15,current_domain))';
                x2 = [x_', fliplr(x_')];
                inBetween = [int_min, fliplr(int_max)];
                h_=fill(x2, inBetween, p.domain_color(current_domain,:));
                set(h_,'LineStyle','none','FaceAlpha',0.25);

                % Extremums and Whole volume
                h_min=plot(x_,p.Property_subdomains_statistics(:,8,current_domain));
                xmin = min(x_); xmax=max(x_);
                h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_domain) p.Wholevolume_results(current_domain)]);
                h_max=plot(x_,p.Property_subdomains_statistics(:,9,current_domain));

                % Colors, thickness, markers
                set(h_min, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'LineStyle','--');
                set(h_max, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'LineStyle','--');
                set(h_whole, 'Color', 'k','LineWidth',optformat.linewidth+1,'LineStyle',':');

                % Annotation: number of subdomains
                for current_size=1:1:length(x_)
                    yl = ylim ;
                    y_bottom = p.Property_subdomains_statistics(current_size,5,current_domain)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,5,current_domain)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_domain) ));
                    set(h,'Color','k','FontSize',optformat.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',optformat.fontname);
                end
                % - Legend
                h_legend = legend(ax,'Mean +/- standard deviation','95% interval confidence','Extremums',['Whole volume, ' num2str(p.Wholevolume_size(kx_whole),'%1.1f') ' ' lgdunit],'Location','best');
                ylabel(p.propertyname);
                if ~isempty(p.propertynameunit)
                    ysecondarylabel(ax,p.propertynameunit)
                end

            elseif id_axe==4 || id_axe==8    
                ylabel(str_threshold);
                if strcmp(p.RVE.threshold_subs_choice,'Relative standard deviation')
                    ysecondarylabel('%')
                elseif strcmp(p.RVE.threshold_subs_choice,'Standard deviation')
                    if ~isempty(p.propertynameunit)
                        ysecondarylabel(p.propertynameunit)
                    end
                elseif strcmp(p.RVE.threshold_subs_choice,'Maximum - minimum')
                    if ~isempty(p.propertynameunit)
                        ysecondarylabel(p.propertynameunit)
                    end
                end
                h_subtitle = subtitle(str_threshold);
                
                h=plot(x_,p.Property_subdomains_statistics(:,idx_threshold,current_domain)); % Curves
                set(h, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o'); % Colors
                x_RVE = [x_(1) x_(end)];

                if sum(p.Size_RVE(:,current_domain,2,ksizeRVE)~=0)==1 % Only plot if only one threshold has been selected
                    for k_threshold=1:1:length(p.Criterion)-1
                        if p.Size_RVE(k_threshold,current_domain,2,ksizeRVE)~=0
                            y_RVE = [p.Criterion(k_threshold) p.Criterion(k_threshold)];
                            h_RVE = plot(x_RVE,y_RVE);
                            set(h_RVE, 'Color', 'black','LineStyle','--','LineWidth',optformat.linewidth);
                            h=plot([p.Size_RVE(k_threshold,current_domain,2,ksizeRVE) p.Size_RVE(k_threshold,current_domain,2,ksizeRVE)],[0 p.Criterion(k_threshold)]); % Curves
                            set(h, 'Color', 'k','LineStyle','--','LineWidth',optformat.linewidth); % Colors
                        end
                    end
                end

                % Annotation: number of subdomains
                for current_size=1:1:length(x_)
                    yl = ylim ;
                    y_bottom = p.Property_subdomains_statistics(current_size,idx_threshold,current_domain)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,idx_threshold,current_domain)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_domain) ));
                    set(h,'Color','k','FontSize',optformat.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',optformat.fontname);
                end
                %h_legend.FontSize = optformat.legendfontsize; % Set fontsize
                h_legend = [];

            elseif id_axe==3 || id_axe==7
                h_subtitle = subtitle('95% confidence interval');
                h=plot(x_, p.Property_subdomains_statistics(:,15,current_domain)); % Curves
                set(h, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o'); % Colors
                x_RVE = [x_(1) x_(end)];

                % Annotation: number of subdomains
                for current_size=1:1:length(x_)
                    yl = ylim ;
                    y_bottom = p.Property_subdomains_statistics(current_size,15,current_domain)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,15,current_domain)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_domain) ));
                    set(h,'Color','k','FontSize',optformat.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',optformat.fontname);
                end
                %h_legend.FontSize = optformat.legendfontsize; % Set fontsize
                h_legend = [];
                ylabel([p.propertyname ': +/- 2.5%']);

            elseif id_axe==1 || id_axe==5
                h_subtitle = subtitle('Points cloud');                
                % Mean
                h_mean=plot(x_,p.Property_subdomains_statistics(:,5,current_domain));
                % Point cloud
                x_= p.Property_eachsubdomain(:,kxcloud);
                y_= p.Property_eachsubdomain(:,current_domain+4);
                h_pointcloud = scatter(x_,y_);
                % Whole volume
                xmin = min(x_); xmax=max(x_);
                h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_domain) p.Wholevolume_results(current_domain)]);
                % Colors, thickness, markers
                set(h_mean, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o');
                set(h_pointcloud, 'MarkerEdgeColor', p.domain_color(current_domain,:));
                set(h_whole, 'Color', 'k','LineWidth',optformat.linewidth+1,'LineStyle',':');
                % Annotation: number of subdomains
                x_=p.Property_subdomains_statistics(:,kx,current_domain);
                for current_size=1:1:length(x_)
                    yl = ylim;
                    y_bottom = p.Property_subdomains_statistics(current_size,5,current_domain)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,5,current_domain)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_domain) ));
                    set(h,'Color','k','FontSize',optformat.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',optformat.fontname);
                end
                h_legend = legend(ax,'Mean','Point cloud',['Whole volume, ' num2str(p.Wholevolume_size(kx_whole),'%1.1f') ' ' lgdunit],'Location','best');
                ylabel(p.propertyname);
                if ~isempty(p.propertynameunit)
                    ysecondarylabel(ax,p.propertynameunit)
                end                
            end

            xlabel(str_xlabel);
            xsecondarylabel(ax,axisunit)

            if ~isempty(h_legend)
                h_legend.FontSize = optformat.legendfontsize;
                h_legend.Location = 'best';
                h_legend.IconColumnWidth = 20;
                h_legend.LineWidth = 1;
            end

            ax.XLabel.FontWeight = "bold";
            ax.YLabel.FontWeight = "bold";
            ax.LineWidth = optformat.linewidth;
            ax.FontName = optformat.fontname;
            ax.FontSize = optformat.axefontsize;
            if optformat.grid
                grid(ax,'on');
                set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
            end

            h_title.FontSize = optformat.titlefontsize; % Set title fontsize
            h_subtitle.FontSize = optformat.titlefontsize-2; % Set title fontsize
            
            hold(ax,'off'); % Relase figure
            if r_axe==1
                str_title_1 = [char(p.Analysis_name(1)) ' for ' p.propertyname ' on ' char(p.domain_name(current_domain,1))];
            else
                str_title_1 = ['1st row: ' char(p.Analysis_name(1)) ', 2nd row: ' char(p.Analysis_name(2)) ' for ' p.propertyname ' on ' char(p.domain_name(current_domain,1))];
            end
            str_title_2 = p.RVE.choice_RVE;
            if ~strcmp(p.RVE.Constantdirection,'n/a')
                str_title_2 = [str_title_2 ': ' p.RVE.Constantdirection];
            end
            if ~ischar(p.RVE.Aspectratio_RVE)
                str_title_2 = [str_title_2 ', aspect ratio: ' num2str(p.RVE.Aspectratio_RVE)];
            end
            str_title = {str_title_1,str_title_2};
            sgtitle(Fig,str_title,'FontWeight','bold','FontSize',optformat.sgtitlefontsize,'FontName',optformat.fontname);
        end

        if optsave.savefig % Save figure
            filename= [p.propertyname '_' char(p.domain_name(current_domain,1)) '_' p.RVE.savename];
            function_savefig(Fig, p.savefolder, filename, optsave); % Call function
        end
        if optsave.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end

    % Relative standard deviation, for all phases
    Fig = figure; % Create figure
    Fig.Name= [p.propertyname ' ' str_threshold]; % Figure name
    Fig.Color='white'; % Background colour
    optsave.Height = h1;
    optsave.Width = n2_axe*w1;
    for id_axe=1:1:n2_axe % Iterate over axe
        str_xlabel=['Subdomain FOV ' char(p.Length_name(id_axe))];
        kx = idx_size(id_axe);
        kxRVE = id_axe;
        ax = nexttile;
        hold(ax,'on');
        h_title=title (p.propertyname,'FontName',optformat.fontname,'FontSize',optformat.titlefontsize); % Set title font and string
        % - Plot graphs
        for current_domain=1:1:number_domain
            str_legend(current_domain).name = char(p.domain_name(current_domain));
            x_=p.Property_subdomains_statistics(:,kx,current_domain);
            h=plot(x_,p.Property_subdomains_statistics(:,idx_threshold,current_domain)); % Curves
            set(h, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o'); % Colors
        end
        % RVE standard deviation criterion
        if sum(p.Size_RVE(:,current_domain,2,kxRVE)~=0)==1 % Only plot if only one threshold has been selected
            x_RVE = [x_(1) x_(end)];
            for k_threshold=1:1:length(p.Criterion)-1
                if p.Size_RVE(k_threshold,current_domain,2,kxRVE)~=0
                    y_RVE = [p.Criterion(k_threshold) p.Criterion(k_threshold)];
                    h_RVE = plot(x_RVE,y_RVE);
                    set(h_RVE, 'Color', 'black','LineStyle','--','LineWidth',optformat.linewidth);
                end
            end
            for current_domain=1:1:number_domain % Loop over all phases
                for k_threshold=1:1:length(p.Criterion)-1
                    if p.Size_RVE(k_threshold,current_domain,2,kxRVE)~=0
                        h=plot([p.Size_RVE(k_threshold,current_domain,2,kxRVE) p.Size_RVE(k_threshold,current_domain,2,kxRVE)],[0 p.Criterion(k_threshold)]); % Curves
                        set(h, 'Color', p.domain_color(current_domain,:),'LineStyle','--','LineWidth',optformat.linewidth); % Colors
                    end
                end
            end
        end

        h_legend = legend(ax,str_legend.name);
        h_legend.FontSize = optformat.legendfontsize;
        h_legend.Location = 'best';
        h_legend.IconColumnWidth = 20;
        h_legend.LineWidth = 1;

        % - Axis label
        xlabel(str_xlabel);
        xsecondarylabel(ax,axisunit)
        ylabel(str_threshold);
        if strcmp(p.RVE.threshold_subs_choice,'Relative standard deviation')
            ysecondarylabel('%')
        elseif strcmp(p.RVE.threshold_subs_choice,'Standard deviation')
            if ~isempty(p.propertynameunit)
                ysecondarylabel(p.propertynameunit)
            end
        elseif strcmp(p.RVE.threshold_subs_choice,'Maximum - minimum')
            if ~isempty(p.propertynameunit)
                ysecondarylabel(p.propertynameunit)
            end
        end

        ax.XLabel.FontWeight = "bold";
        ax.YLabel.FontWeight = "bold";
        ax.LineWidth = optformat.linewidth;
        ax.FontName = optformat.fontname;
        ax.FontSize = optformat.axefontsize;
        if optformat.grid
            grid(ax,'on');
            set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
        end
        h_title.FontSize = optformat.titlefontsize; % Set title fontsize
        sgtitle(Fig,str_title_2,'FontWeight','bold','FontSize',optformat.sgtitlefontsize,'FontName',optformat.fontname);
        hold(ax,'off'); % Relase figure
    end
    if optsave.savefig % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, optsave); % Call function
    end
    if optsave.autoclosefig
        close(Fig); % Do not keep open figures
    end

    % RVE as function of threshold
    if sum(sum(p.Size_RVE(:,:,2,1)~=0) >= 2) % At least 2 points
        Fig = figure; % Create figure
        Fig.Name= [p.propertyname ' RVE size as function of ' str_threshold ' threshold']; % Figure name
        Fig.Color='white'; % Background colour
        optsave.Height = h1;
        optsave.Width = n2_axe*w1;
        for id_axe=1:1:n2_axe % Iterate over axe
            kxRVE = id_axe;
            str = char(p.Analysis_name(id_axe));
            idx = find(str == '(');
            short = str(idx+1:end-1);
            str_ylabel = [short ' ' char(p.Length_name(id_axe))];
            ax = nexttile;
            hold(ax,'on');
            h_title=title (p.propertyname,'FontName',optformat.fontname,'FontSize',optformat.titlefontsize); % Set title font and string
            % - Plot graphs
            for current_domain=1:1:number_domain
                str_legend(current_domain).name = char(p.domain_name(current_domain));
                x_= p.Criterion(1:end-1);
                y_ = p.Size_RVE(:,current_domain,2,kxRVE);
                y_(y_==0)=NaN;
                h=plot(x_,y_); % Curves
                set(h, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o'); % Colors
            end

            h_legend = legend(ax,str_legend.name);
            h_legend.FontSize = optformat.legendfontsize;
            h_legend.Location = 'best';
            h_legend.IconColumnWidth = 20;
            h_legend.LineWidth = 1;

            % - Axis label
            xlabel(str_threshold);
            if strcmp(p.RVE.threshold_subs_choice,'Relative standard deviation')
                xsecondarylabel('%')
            elseif strcmp(p.RVE.threshold_subs_choice,'Standard deviation')
                if ~isempty(p.propertynameunit)
                    xsecondarylabel(p.propertynameunit)
                end
            elseif strcmp(p.RVE.threshold_subs_choice,'Maximum - minimum')
                if ~isempty(p.propertynameunit)
                    xsecondarylabel(p.propertynameunit)
                end
            end
            ylabel(str_ylabel);
            ysecondarylabel(axisunit)

            ax.XLabel.FontWeight = "bold";
            ax.YLabel.FontWeight = "bold";
            ax.LineWidth = optformat.linewidth;
            ax.FontName = optformat.fontname;
            ax.FontSize = optformat.axefontsize;
            if optformat.grid
                grid(ax,'on');
                set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
            end
            h_title.FontSize = optformat.titlefontsize; % Set title fontsize
            sgtitle(Fig,str_title_2,'FontWeight','bold','FontSize',optformat.sgtitlefontsize,'FontName',optformat.fontname);
            hold(ax,'off'); % Relase figure
        end
        if optsave.savefig % Save figure
            filename= [p.propertyname '_fctThreshold_' p.RVE.savename];
            function_savefig(Fig, p.savefolder, filename, optsave); % Call function
        end
        if optsave.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end


else

    if length(p.Wholevolume_size)==2
        r_axe=2; n_axe=4;
        n2_axe=2;
    else
        r_axe=1; n_axe=2;
        n2_axe=1;
    end

    optsave.Height = r_axe*h1;
    optsave.Width = 2*w1;

    str_threshold = p.RVE.threshold_onesub_choice;

    % Property evolution for all phases
    Fig = figure; % Create figure
    Fig.Name= ['Convergence analysis:' p.propertyname];
    Fig.Color='white'; % Background colour
    tiledlayout(r_axe,2,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);
    for id_axe=1:1:n_axe % Iterate over axe
        if id_axe<=2
            str_xlabel=['Subdomain FOV ' char(p.Length_name(1))];
            kx = idx_size(1);
            kx_whole=1;
            ksizeRVE = 1;
            kconv = 1;
        else
            str_xlabel=['Subdomain FOV ' char(p.Length_name(2))];
            kx = idx_size(2);
            kx_whole=2;
            ksizeRVE = 2;
            kconv = 2;
        end
        ax = nexttile;
        hold(ax,'on');

        if id_axe==1 || id_axe==3
            h_title=title ([p.propertyname ' evolution'],'FontName',optformat.fontname,'FontSize',optformat.titlefontsize); % Set title font and string
        else
            h_title=title ([p.propertyname ' ' str_threshold],'FontName',optformat.fontname,'FontSize',optformat.titlefontsize); % Set title font and string
        end
        for current_domain=1:1:number_domain
            str_legend(current_domain).name = char(p.domain_name(current_domain));
            if id_axe==1 || id_axe==3
                x_=p.Property_subdomains_statistics(:,kx,current_domain);
                x_=[x_; p.Wholevolume_size(kx_whole)];
                y_ = p.Property_subdomains_statistics(:,5,current_domain);
                y_ = [y_; p.Wholevolume_results(current_domain)];
            elseif id_axe==2
                x_=p.Difference_convergence(:,1);
                y_=p.Difference_convergence(:,current_domain+2);
                x_conv = [x_(1) x_(end)];
            else
                x_=p.Difference_convergence(:,2);
                y_=p.Difference_convergence(:,current_domain+2);
                x_conv = [x_(1) x_(end)];
            end
            h_=plot(x_,y_);
            % Colors, thickness, markers
            set(h_, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o');
        end

        % Other plot with no legend
        if sum(p.Size_convergence(:,current_domain,2,ksizeRVE)~=0)==1 % Only plot if only one threshold has been selected
            if id_axe==2 || id_axe==4
                for k_threshold=1:1:length(p.convergence_criterion)
                    if p.Size_convergence(k_threshold,current_domain,2,kconv)~=0
                        y_conv = [p.convergence_criterion(k_threshold) p.convergence_criterion(k_threshold)];
                        h_conv = plot(x_conv,y_conv);
                        set(h_conv, 'Color', 'black','LineStyle','--','LineWidth',optformat.linewidth);

                        h=plot([p.Size_convergence(k_threshold,current_domain,2,kconv) p.Size_convergence(k_threshold,current_domain,2,kconv)],[0 p.convergence_criterion(k_threshold)]); % Curves
                        set(h, 'Color', 'k','LineStyle','--','LineWidth',optformat.linewidth); % Colors
                    end
                end
            end
        end

        h_legend = legend(ax,str_legend.name);
        h_legend.FontSize = optformat.legendfontsize;
        h_legend.Location = 'best';
        h_legend.IconColumnWidth = 20;
        h_legend.LineWidth = 1;

        % - Axis label
        xlabel(str_xlabel);
        xsecondarylabel(axisunit);
        if id_axe==1 || id_axe==3
            ylabel(p.propertyname)
            if ~isempty(p.propertynameunit)
                ysecondarylabel(p.propertynameunit)
            end
        else
            ylabel([p.propertyname ' ' str_threshold]);
            if strcmp(p.RVE.threshold_onesub_choice,'Relative difference')
                ysecondarylabel('%')
            elseif strcmp(p.RVE.threshold_onesub_choice,'Difference')
                if ~isempty(p.propertynameunit)
                    ysecondarylabel(p.propertynameunit)
                end
            end
        end

        ax.XLabel.FontWeight = "bold";
        ax.YLabel.FontWeight = "bold";
        ax.LineWidth = optformat.linewidth;
        ax.FontName = optformat.fontname;
        ax.FontSize = optformat.axefontsize;
        if optformat.grid
            grid(ax,'on');
            set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
        end
        h_title.FontSize = optformat.titlefontsize; % Set title fontsize
        hold(ax,'off'); % Relase figure
    end

    if r_axe==1
        str_title_1 = [char(p.Analysis_name(1)) ' for ' p.propertyname];
    else
        str_title_1 = ['1st row: ' char(p.Analysis_name(1)) ', 2nd row: ' char(p.Analysis_name(2)) ' for ' p.propertyname];
    end
    str_title_2 = p.RVE.choice_Onsub;
    if ~strcmp(p.RVE.Growthdirection,'n/a')
        str_title_2 = [str_title_2 ': ' p.RVE.Growthdirection];
    end
    if ~ischar(p.RVE.Aspectratio_Onesub)
        str_title_2 = [str_title_2 ', aspect ratio: ' num2str(p.RVE.Aspectratio_Onesub)];
    end
    str_title = {str_title_1,str_title_2};
    sgtitle(Fig,str_title,'FontWeight','bold','FontSize',optformat.sgtitlefontsize,'FontName',optformat.fontname);
    if optsave.savefig % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, optsave); % Call function
    end
    if optsave.autoclosefig
        close(Fig); % Do not keep open figures
    end


    % Convergence length as function of threshold
    if sum(sum(p.Size_convergence(:,:,2,1)~=0) >= 2) % At least 2 points
        Fig = figure; % Create figure
        Fig.Name= [p.propertyname ' Convergence length as function of ' str_threshold ' threshold']; % Figure name
        Fig.Color='white'; % Background colour
        tiledlayout(1,n2_axe,'TileSpacing',optformat.tile_spacing,'Padding',optformat.layout_padding);

        optsave.Height = 1*h1;
        optsave.Width = n2_axe*w1;

        for id_axe=1:1:n2_axe % Iterate over axe
            str_ylabel = ['FOV ' char(p.Length_name(id_axe)) ' to reach convergence'];
            kconv = id_axe;
            ax = nexttile;
            hold(ax,'on');
            str1 = [p.propertyname ' convergence size'];
            str_title = {str1,str_title_2,};
            h_title=title (str_title,'FontName',optformat.fontname,'FontSize',optformat.titlefontsize); % Set title font and string
            % - Plot graphs
            for current_domain=1:1:number_domain
                str_legend(current_domain).name = char(p.domain_name(current_domain));
                x_= p.convergence_criterion;
                y_ = p.Size_convergence(:,current_domain,2,kconv);
                y_(y_==0)=NaN;
                h=plot(x_,y_); % Curves
                set(h, 'Color', p.domain_color(current_domain,:),'LineWidth',optformat.linewidth,'MarkerSize',optformat.axefontsize,'Marker','o'); % Colors
            end

            h_legend = legend(ax,str_legend.name);
            h_legend.FontSize = optformat.legendfontsize;
            h_legend.Location = 'best';
            h_legend.IconColumnWidth = 20;
            h_legend.LineWidth = 1;

            % - Axis label
            xlabel([p.propertyname ' ' str_threshold]);
            if strcmp(p.RVE.threshold_onesub_choice,'Relative difference')
                xsecondarylabel('%')
            elseif strcmp(p.RVE.threshold_onesub_choice,'Difference')
                if ~isempty(p.propertynameunit)
                    xsecondarylabel(p.propertynameunit)
                end
            end
            ylabel(str_ylabel);
            ysecondarylabel(axisunit);

            ax.XLabel.FontWeight = "bold";
            ax.YLabel.FontWeight = "bold";
            ax.LineWidth = optformat.linewidth;
            ax.FontName = optformat.fontname;
            ax.FontSize = optformat.axefontsize;
            if optformat.grid
                grid(ax,'on');
                set(ax,'XMinorGrid',optformat.minorgrid,'YMinorGrid',optformat.minorgrid,'GridLineWidth',1);
            end

            h_title.FontSize = optformat.titlefontsize; % Set title fontsize
            hold(ax,'off'); % Relase figure
        end
        if optsave.savefig % Save figure
            filename= [p.propertyname '_fctThreshold_' p.RVE.savename];
            function_savefig(Fig, p.savefolder, filename, optsave); % Call function
        end
        if optsave.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end


end

end