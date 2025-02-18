function [] = Function_create_figures_RVE2(p)

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter
scrsz = get(0,'ScreenSize'); % Screen resolution

strunit =  p.infovol.unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)'; lgdunit = '\mum';
else
    axisunit = ['(' strunit ')']; lgdunit = strunit;
end 

if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B') || strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'D')

    if strcmp(p.RVE.type,'C') || strcmp(p.RVE.type,'D')
        r_axe=2; n_axe=8; scr=1;
        n2_axe=2; scr2=2/3;
    else
        r_axe=1; n_axe=4; scr=1/2;
        n2_axe=1; scr2=1/3;
    end
    current_phase_todo = 0;
    for current_phase=1:1:p.number_phase % Loop over all phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;

            Fig = figure; % Create figure
            Fig.Name= ['Representative volume element analysis: ' p.propertyname ', ' char(p.infovol.phasename(current_phase,1))];
            Fig.Color='white'; % Background colour
            set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)*scr]); % Full screen figure
            for id_axe=1:1:n_axe % Iterate over axe
                if id_axe<=4
                    if p.dimension == 3
                        str_xlabel=['Subdomain FOV cubic root ' axisunit];
                    else
                        str_xlabel=['Subdomain FOV square root ' axisunit];
                    end
                else
                    if p.dimension == 3
                        if strcmp(p.RVE.type,'C')
                            str_xlabel=['Subdomain FOV square root ' axisunit];
                        elseif strcmp(p.RVE.type,'D')
                            str_xlabel=['Subdomains FOV length ' axisunit];
                        end
                    else
                        if strcmp(p.RVE.type,'C')
                            str_xlabel=['Subdomains FOV length ' axisunit];
                        elseif strcmp(p.RVE.type,'D')
                            str_xlabel=['Subdomains FOV length ' axisunit];
                        end
                    end
                end
                sub_axes=subplot(r_axe,4,id_axe,'Parent',Fig);
                hold(sub_axes,'on'); % Active subplot

                if id_axe==1 ||  id_axe==5
                    if id_axe==1
                        kx = 1; % Use first RVE size
                    else
                        if strcmp(p.RVE.type,'C')
                            kx = 2; % Use second RVE size
                        elseif strcmp(p.RVE.type,'D')
                            kx = 3; % Use third RVE size
                        end
                    end
                    h_title=title ({p.propertyname, 'Statistics'}); % Set title font and string

                    x_=p.Property_subdomains_statistics(:,kx,current_phase_todo);

                    % Mean with error bar (+- standard deviation)
                    x_tmp = []; y_tmp = []; e_tmp = [];
                    for current_size=1:1:length(x_)
                        x_tmp=[x_tmp p.Property_subdomains_statistics(current_size,kx,current_phase_todo)];
                        y_tmp=[y_tmp p.Property_subdomains_statistics(current_size,5,current_phase_todo)];
                        e_tmp=[e_tmp p.Property_subdomains_statistics(current_size,6,current_phase_todo)];
                    end
                    h_mean_witherrorbar = errorbar(x_tmp,y_tmp,e_tmp);
                    set(h_mean_witherrorbar, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o');                    

                    % 5% interval confidence
                    int_min = (p.Property_subdomains_statistics(:,5,current_phase_todo) - p.Property_subdomains_statistics(:,14,current_phase_todo))';
                    int_max = (p.Property_subdomains_statistics(:,5,current_phase_todo) + p.Property_subdomains_statistics(:,14,current_phase_todo))';
                    x2 = [x_', fliplr(x_')];
                    inBetween = [int_min, fliplr(int_max)];
                    h_=fill(x2, inBetween, p.infovol.phasecolor(current_phase,:));
                    set(h_,'LineStyle','none','FaceAlpha',0.25);   

                    % Extremums and Whole volume
                    h_min=plot(x_,p.Property_subdomains_statistics(:,8,current_phase_todo));
                    xmin = min(x_); xmax=max(x_);
                    h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_phase_todo) p.Wholevolume_results(current_phase_todo)]);
                    h_max=plot(x_,p.Property_subdomains_statistics(:,9,current_phase_todo));

                    % Colors, thickness, markers
                    set(h_min, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'LineStyle','--');
                    set(h_max, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'LineStyle','--');
                    set(h_whole, 'Color', 'k','LineWidth',p.opt.format.linewidth+1,'LineStyle',':');

                    % Annotation: number of subdomains
                    for current_size=1:1:length(x_)
                        yl = ylim ;
                        y_bottom = p.Property_subdomains_statistics(current_size,5,current_phase_todo)-0.1*(yl(2)-yl(1));
                        y_top = p.Property_subdomains_statistics(current_size,5,current_phase_todo)+0.1*(yl(2)-yl(1));
                        if y_bottom<yl(1)
                            y_=y_top;
                        else
                            y_=y_bottom;
                        end
                        h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_phase_todo) ));
                        set(h,'Color','k','FontSize',p.opt.format.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',p.opt.format.fontname);
                    end
                    % - Legend
                    h_legend = legend(sub_axes,'Mean +/- standard deviation','95% interval confidence','Extremums',['Whole volume, ' num2str(p.Wholevolume_size(kx),'%1.1f') ' ' lgdunit],'Location','best');
                    h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
                    if isempty(p.propertynameunit)
                        ylabel(p.propertyname);
                    else
                        ylabel([p.propertyname ' (' p.propertynameunit ')']);
                    end

                elseif id_axe==2 || id_axe==6
                    if id_axe==2
                        kx = 1;
                    else
                        if strcmp(p.RVE.type,'C')
                            kx = 2; % Use second RVE size
                        elseif strcmp(p.RVE.type,'D')
                            kx = 2; % Use third RVE size
                        end                        
                    end
                    h_title=title ({p.propertyname,'Relative standard deviation'}); % Set title font and string
                    h=plot(x_,p.Property_subdomains_statistics(:,7,current_phase_todo)); % Curves
                    set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
                    x_RVE = [x_(1) x_(end)];

                    if sum(p.Size_RVE(:,current_phase_todo,2,kx)~=0)==1 % Only plot if only one threshold has been selected
                        for k_threshold=1:1:length(p.Criterion)-1
                            if p.Size_RVE(k_threshold,current_phase_todo,2,kx)~=0
                                y_RVE = [p.Criterion(k_threshold) p.Criterion(k_threshold)];
                                h_RVE = plot(x_RVE,y_RVE);
                                set(h_RVE, 'Color', 'black','LineStyle','--','LineWidth',p.opt.format.linewidth);
                                h=plot([p.Size_RVE(k_threshold,current_phase_todo,2,kx) p.Size_RVE(k_threshold,current_phase_todo,2,kx)],[0 p.Criterion(k_threshold)]); % Curves
                                set(h, 'Color', 'k','LineStyle','--','LineWidth',p.opt.format.linewidth); % Colors
                            end
                        end
                    end

                    % Annotation: number of subdomains
                    for current_size=1:1:length(x_)
                        yl = ylim ;
                        y_bottom = p.Property_subdomains_statistics(current_size,7,current_phase_todo)-0.1*(yl(2)-yl(1));
                        y_top = p.Property_subdomains_statistics(current_size,7,current_phase_todo)+0.1*(yl(2)-yl(1));
                        if y_bottom<yl(1)
                            y_=y_top;
                        else
                            y_=y_bottom;
                        end
                        h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_phase_todo) ));
                        set(h,'Color','k','FontSize',p.opt.format.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',p.opt.format.fontname);
                    end
                    %h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
                    ylabel('Relative standard deviation (%)');


                elseif id_axe==3 || id_axe==7
                    if id_axe==3
                        kx = 1;
                    else
                        if strcmp(p.RVE.type,'C')
                            kx = 2; % Use second RVE size
                        elseif strcmp(p.RVE.type,'D')
                            kx = 3; % Use third RVE size
                        end
                    end
                    h_title=title ({p.propertyname,'95% confidence interval'}); % Set title font and string
                    h=plot(x_, p.Property_subdomains_statistics(:,14,current_phase_todo)); % Curves
                    set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
                    x_RVE = [x_(1) x_(end)];

                    % Annotation: number of subdomains
                    for current_size=1:1:length(x_)
                        yl = ylim ;
                        y_bottom = p.Property_subdomains_statistics(current_size,14,current_phase_todo)-0.1*(yl(2)-yl(1));
                        y_top = p.Property_subdomains_statistics(current_size,14,current_phase_todo)+0.1*(yl(2)-yl(1));
                        if y_bottom<yl(1)
                            y_=y_top;
                        else
                            y_=y_bottom;
                        end
                        h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_phase_todo) ));
                        set(h,'Color','k','FontSize',p.opt.format.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',p.opt.format.fontname);
                    end
                    %h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
                    if isempty(p.propertynameunit)
                        ylabel([p.propertyname ': +/- 2.5%']);
                    else
                        ylabel([p.propertyname ' (' p.propertynameunit '): +/- 2.5%']);
                    end


                elseif id_axe==4 || id_axe==8
                    if id_axe==4
                        kx = 1; kxcloud = 2;
                    else
                        if strcmp(p.RVE.type,'C')
                            kx = 2; % Use second RVE size
                            kxcloud = 3;
                        elseif strcmp(p.RVE.type,'D')
                            kx = 3; % Use third RVE size
                            kxcloud = 4;
                        end
                    end
                    h_title=title ({p.propertyname, 'Points cloud'}); % Set title font and string
                    % Mean
                    h_mean=plot(x_,p.Property_subdomains_statistics(:,5,current_phase_todo));
                    % Point cloud
                    x_= p.Property_eachsubdomain(:,kxcloud);
                    y_= p.Property_eachsubdomain(:,current_phase_todo+4);
                    h_pointcloud = scatter(x_,y_);
                    % Whole volume
                    xmin = min(x_); xmax=max(x_);
                    h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_phase_todo) p.Wholevolume_results(current_phase_todo)]);
                    % Colors, thickness, markers
                    set(h_mean, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o');
                    set(h_pointcloud, 'MarkerEdgeColor', p.infovol.phasecolor(current_phase,:));
                    set(h_whole, 'Color', 'k','LineWidth',p.opt.format.linewidth+1,'LineStyle',':');
                    % Annotation: number of subdomains
                    x_=p.Property_subdomains_statistics(:,kx,current_phase_todo);
                    for current_size=1:1:length(x_)
                        yl = ylim;
                        y_bottom = p.Property_subdomains_statistics(current_size,5,current_phase_todo)-0.1*(yl(2)-yl(1));
                        y_top = p.Property_subdomains_statistics(current_size,5,current_phase_todo)+0.1*(yl(2)-yl(1));
                        if y_bottom<yl(1)
                            y_=y_top;
                        else
                            y_=y_bottom;
                        end
                        h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_phase_todo) ));
                        set(h,'Color','k','FontSize',p.opt.format.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',p.opt.format.fontname);
                    end
                    h_legend = legend(sub_axes,'Mean','Point cloud',['Whole volume, ' num2str(p.Wholevolume_size(kx),'%1.1f') ' ' lgdunit],'Location','best');
                    h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
                    if isempty(p.propertynameunit)
                        ylabel(p.propertyname);
                    else
                        ylabel([p.propertyname ' (' p.propertynameunit ')']);
                    end

                end

                xlabel(str_xlabel);
                grid(sub_axes,p.opt.format.grid); % Display grid
                set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
                set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
                h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
                hold(sub_axes,'off'); % Relase figure
                if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
                    str_title = {['Representative volume element analysis: ' p.propertyname ', ' char(p.infovol.phasename(current_phase,1))],p.RVE.choice_RVE, ['Aspect ratio: ' num2str(p.RVE.Aspectratio_RVE)]};
                else
                    str_title = {['Representative volume element analysis: ' p.propertyname ', ' char(p.infovol.phasename(current_phase,1))],p.RVE.choice_RVE};
                end
                sgtitle(Fig,str_title,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);
            end

            if p.opt.save.savefig % Save figure
                filename= [p.propertyname '_' char(p.infovol.phasename(current_phase,1)) '_' p.RVE.savename];
                function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
            end
            if p.opt.format.autoclosefig
                close(Fig); % Do not keep open figures
            end

        end
    end


%     if strcmp(p.RVE.type,'C')
%         r_axe=2; n_axe=6; scr=1;
%         n2_axe=2; scr2=2/3;
%     else
%         r_axe=1; n_axe=3; scr=1/2;
%         n2_axe=1; scr2=1/3;
%     end
%     current_phase_todo = 0;
%     for current_phase=1:1:p.number_phase % Loop over all phases
%         if p.todo(current_phase)
%             current_phase_todo=current_phase_todo+1;
% 
%             Fig = figure; % Create figure
%             Fig.Name= ['Error analysis: ' p.propertyname ', ' char(p.infovol.phasename(current_phase,1))];
%             Fig.Color='white'; % Background colour
%             set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)*scr]); % Full screen figure
%             for id_axe=1:1:n_axe % Iterate over axe
%                 if id_axe<=3
%                     str_xlabel=['Equivalent length: cubic root of the subdomains ' axisunit];
%                 else
%                     str_xlabel=['Equivalent length: square root of the subdomains FOV ' axisunit];
%                 end
%                 sub_axes=subplot(r_axe,3,id_axe,'Parent',Fig);
%                 hold(sub_axes,'on'); % Active subplot
%                 
%                 if id_axe==1 || id_axe==4
%                     if id_axe==1
%                         kx = 1;
%                     else
%                         kx = 2;
%                     end
%                     h_title=title ({p.propertyname,'95% confidence interval'}); % Set title font and string
%                     x_=p.Property_subdomains_statistics(:,kx,current_phase_todo);
%                     h=plot(x_, p.Property_subdomains_statistics(:,14,current_phase_todo)); % Curves
%                     set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
%                     x_RVE = [x_(1) x_(end)];                    
% 
%                     % Annotation: number of subdomains
%                     for current_size=1:1:length(x_)
%                         yl = ylim ;
%                         y_bottom = 2*p.Property_subdomains_statistics(current_size,14,current_phase_todo)-0.1*(yl(2)-yl(1));
%                         y_top = 2*p.Property_subdomains_statistics(current_size,14,current_phase_todo)+0.1*(yl(2)-yl(1));
%                         if y_bottom<yl(1)
%                             y_=y_top;
%                         else
%                             y_=y_bottom;
%                         end
%                         h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_phase_todo) ));
%                         set(h,'Color','k','FontSize',p.opt.format.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',p.opt.format.fontname);
%                     end
%                     %h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
%                     if isempty(p.propertynameunit)
%                         ylabel([p.propertyname ': +/- 2.5%']);
%                     else
%                         ylabel([p.propertyname ' (' p.propertynameunit '): +/- 2.5%']);
%                     end
% 
%                 elseif id_axe==2 ||  id_axe==5
%                     if id_axe==2
%                         kx = 1; % Use first RVE size
%                     else
%                         kx = 2; % Use second RVE size
%                     end
%                     h_title=title ({p.propertyname, 'Number of subvolumes required to reach a given error on the mean'}); % Set title font and string
% 
%                     x_=p.Property_subdomains_statistics(:,kx,current_phase_todo);
% 
%                     h_1p=plot(x_,p.Property_subdomains_statistics(:,15,current_phase_todo));
%                     h_5p=plot(x_,p.Property_subdomains_statistics(:,16,current_phase_todo));
%                     h_10p=plot(x_,p.Property_subdomains_statistics(:,17,current_phase_todo));
%                     
%                     % Colors, thickness, markers
%                     set(h_1p, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
%                     set(h_5p, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','square'); % Colors
%                     set(h_10p, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','diamond'); % Colors
% 
%                     % Annotation: number of subdomains
%                     for current_size=1:1:length(x_)
%                         yl = ylim ;
%                         y_bottom = p.Property_subdomains_statistics(current_size,16,current_phase_todo)-0.1*(yl(2)-yl(1));
%                         y_top = p.Property_subdomains_statistics(current_size,16,current_phase_todo)+0.1*(yl(2)-yl(1));
%                         if y_bottom<yl(1)
%                             y_=y_top;
%                         else
%                             y_=y_bottom;
%                         end
%                         h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,4,current_phase_todo) ));
%                         set(h,'Color','k','FontSize',p.opt.format.axefontsize,'HorizontalAlignment','center','FontWeight','bold','FontName',p.opt.format.fontname);
%                     end
%                     set(sub_axes, 'YScale', 'log')
%                     % - Legend
%                     h_legend = legend(sub_axes,'1% error','5% error','10% error','Location','best');
%                     h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
%                     ylabel('Number of subvolumes');
% 
%                 elseif id_axe==3 || id_axe==6
%                     if id_axe==3
%                         kx = 1; kxcloud = 2;
%                     else
%                         kx = 2; kxcloud = 3;
%                     end
%                     h_title=title ({p.propertyname, 'Points cloud'}); % Set title font and string
%                     
% 
%                 end
% 
%                 xlabel(str_xlabel);
%                 grid(sub_axes,p.opt.format.grid); % Display grid
%                 set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
%                 set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
%                 h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
%                 hold(sub_axes,'off'); % Relase figure
%                 if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
%                     str_title = {['Error analysis: ' p.propertyname ', ' char(p.infovol.phasename(current_phase,1))],p.RVE.name, ['Aspect ratio: ' num2str(p.RVE.Aspectratio)]};
%                 else
%                     str_title = {['Error analysis: ' p.propertyname ', ' char(p.infovol.phasename(current_phase,1))],p.RVE.name};
%                 end
%                 sgtitle(Fig,str_title,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);
%             end
% 
%             if p.opt.save.savefig % Save figure
%                 filename= [p.propertyname '_' char(p.infovol.phasename(current_phase,1)) '_error_' p.RVE.savename];
%                 function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
%             end
%             if p.opt.format.autoclosefig
%                 close(Fig); % Do not keep open figures
%             end
% 
%         end
%     end


    % Relative standard deviation, for all phases
    Fig = figure; % Create figure
    Fig.Name= [p.propertyname ' relative standard deviation']; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*scr2 scrsz(4)*1/2]); % Full screen figure
    for id_axe=1:1:n2_axe % Iterate over axe
        if id_axe==1
            if p.dimension == 3
                str_xlabel=['Subdomain FOV cubic root ' axisunit];
            else
                str_xlabel=['Subdomain FOV square root ' axisunit];
            end
            kx = 1;
            kxRVE = 1;
        else
            if p.dimension == 3
                if strcmp(p.RVE.type,'C')
                    str_xlabel=['Subdomain FOV square root ' axisunit];
                    kx = 2;
                    kxRVE = 2;
                elseif strcmp(p.RVE.type,'D')
                    str_xlabel=['Subdomains FOV length ' axisunit];
                    kx = 3;
                    kxRVE = 2;
                end
            else
                if strcmp(p.RVE.type,'C')
                    str_xlabel=['Subdomains FOV length ' axisunit];
                    kx = 2;
                    kxRVE = 2;
                elseif strcmp(p.RVE.type,'D')
                    str_xlabel=['Subdomains FOV length ' axisunit];
                    kx = 3;
                    kxRVE = 2;
                end
            end
        end
        sub_axes=subplot(1,n2_axe,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        str1 = [p.propertyname ' relative standard deviation'];
        str2 = p.RVE.choice_RVE;
        if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
            str_title = {str1,str2,['Aspect ratio: ' num2str(p.RVE.Aspectratio_RVE)]};
        else
            str_title = {str1,str2,};
        end
        h_title=title (str_title,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.titlefontsize); % Set title font and string
        % - Plot graphs
        current_phase_todo = 0;
        for current_phase=1:1:p.number_phase
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                str_legend(current_phase_todo).name = char(p.infovol.phasename(current_phase));
                x_=p.Property_subdomains_statistics(:,kx,current_phase_todo);
                h=plot(x_,p.Property_subdomains_statistics(:,7,current_phase_todo)); % Curves
                set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
            end
        end
        % RVE standard deviation criterion
        if sum(p.Size_RVE(:,current_phase_todo,2,kxRVE)~=0)==1 % Only plot if only one threshold has been selected
            x_RVE = [x_(1) x_(end)];
            for k_threshold=1:1:length(p.Criterion)-1
                if p.Size_RVE(k_threshold,current_phase_todo,2,kxRVE)~=0
                    y_RVE = [p.Criterion(k_threshold) p.Criterion(k_threshold)];
                    h_RVE = plot(x_RVE,y_RVE);
                    set(h_RVE, 'Color', 'black','LineStyle','--','LineWidth',p.opt.format.linewidth);
                end
            end
            current_phase_todo = 0;
            for current_phase=1:1:p.number_phase % Loop over all phases
                if p.todo(current_phase)
                    current_phase_todo=current_phase_todo+1;
                    for k_threshold=1:1:length(p.Criterion)-1
                        if p.Size_RVE(k_threshold,current_phase_todo,2,kxRVE)~=0
                            h=plot([p.Size_RVE(k_threshold,current_phase_todo,2,kxRVE) p.Size_RVE(k_threshold,current_phase_todo,2,kxRVE)],[0 p.Criterion(k_threshold)]); % Curves
                            set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineStyle','--','LineWidth',p.opt.format.linewidth); % Colors
                        end
                    end
                end
            end
        end

        h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
        h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
        % - Axis label
        xlabel(str_xlabel);
        ylabel('Relative standard deviation (%)');
        grid(sub_axes,p.opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
        h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    if p.opt.save.savefig % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
    end
    if p.opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end

    % RVE as function of threshold
    if sum(sum(p.Size_RVE(:,:,2,1)~=0) >= 2) % At least 2 points
        Fig = figure; % Create figure
        Fig.Name= [p.propertyname ' RVE size as function of standard deviation threshold']; % Figure name
        Fig.Color='white'; % Background colour
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*scr2 scrsz(4)*1/2]); % Full screen figure
        for id_axe=1:1:n2_axe % Iterate over axe
            if id_axe==1
                if p.dimension == 3
                    str_ylabel=['RVE cubic root ' axisunit];
                else
                    str_ylabel=['RSA square root ' axisunit];                    
                end
                kx = 1;
            else
                if p.dimension == 3
                    if strcmp(p.RVE.type,'C')
                        str_ylabel=['RSA square root ' axisunit];
                        kx = 2;
                    elseif strcmp(p.RVE.type,'D')
                        str_ylabel=['Representative length ' axisunit];
                        kx = 2;
                    end
                else
                    if strcmp(p.RVE.type,'C')
                        str_ylabel=['Representative length ' axisunit];
                        kx = 2;
                    elseif strcmp(p.RVE.type,'D')
                        str_ylabel=['Representative length ' axisunit];
                        kx = 2;
                    end
                end
            end
            sub_axes=subplot(1,n2_axe,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            str1 = [p.propertyname ' RVE size'];
            str2 = p.RVE.choice_RVE;
            if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
                str_title = {str1,str2,['Aspect ratio: ' num2str(p.RVE.Aspectratio_RVE)]};
            else
                str_title = {str1,str2,};
            end
            h_title=title (str_title,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.titlefontsize); % Set title font and string
            % - Plot graphs
            current_phase_todo = 0;
            for current_phase=1:1:p.number_phase
                if p.todo(current_phase)
                    current_phase_todo=current_phase_todo+1;
                    str_legend(current_phase_todo).name = char(p.infovol.phasename(current_phase));
                    x_= p.Criterion(1:end-1);
                    y_ = p.Size_RVE(:,current_phase_todo,2,kx);
                    y_(y_==0)=NaN;
                    h=plot(x_,y_); % Curves
                    set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
                end
            end
            
            h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
            h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
            % - Axis label
            xlabel('Relative standard deviation threshold (%)');
            ylabel(str_ylabel);
            grid(sub_axes,p.opt.format.grid); % Display grid
            set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
        if p.opt.save.savefig % Save figure
            filename= [p.propertyname '_fctThreshold_' p.RVE.savename];
            function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
        end
        if p.opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end


elseif strcmp(p.RVE.type,'E') || strcmp(p.RVE.type,'F') || strcmp(p.RVE.type,'G') || strcmp(p.RVE.type,'H')
   
    if strcmp(p.RVE.type,'G') || strcmp(p.RVE.type,'H')
        r_axe=2; n_axe=4; scr=1;
        n2_axe=2; scr2=2/3;
    else
        r_axe=1; n_axe=2; scr=1/2;
        n2_axe=1; scr2=1/3;
    end

    % Property evolution for all phases
    Fig = figure; % Create figure
    Fig.Name= ['Convergence analysis:' p.propertyname];
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)*scr]); % Full screen figure
    for id_axe=1:1:n_axe % Iterate over axe
        if id_axe==1 || id_axe==2
            if p.dimension == 3
                str_xlabel=['Subdomain FOV cubic root ' axisunit];
            else
                str_xlabel=['Subdomain FOV square root ' axisunit];
            end
            kx = 1; kxx=1;
        else
            if p.dimension == 3
                if strcmp(p.RVE.type,'G')
                    str_xlabel=['Subdomains FOV length ' axisunit];
                    kx = 3;
                elseif strcmp(p.RVE.type,'H')
                    str_xlabel=['Subdomain FOV square root ' axisunit];
                    kx = 2;
                end
            else
                if strcmp(p.RVE.type,'G')
                    str_xlabel=['Subdomains FOV length ' axisunit];
                    kx = 3;
                elseif strcmp(p.RVE.type,'H')
                    str_xlabel=['Subdomains FOV length ' axisunit];
                    kx = 2;
                end
            end
            kxx=2;
        end
        sub_axes=subplot(r_axe,2,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1 || id_axe==3
            h_title=title ([p.propertyname ' evolution'],'FontName',p.opt.format.fontname,'FontSize',p.opt.format.titlefontsize); % Set title font and string
        else
            h_title=title ([p.propertyname ' relative difference'],'FontName',p.opt.format.fontname,'FontSize',p.opt.format.titlefontsize); % Set title font and string
        end
        current_phase_todo = 0;
        for current_phase=1:1:p.number_phase
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                str_legend(current_phase_todo).name = char(p.infovol.phasename(current_phase));
                if id_axe==1 || id_axe==3
                    x_=p.Property_subdomains_statistics(:,kx,current_phase_todo);
                    x_=[x_; p.Wholevolume_size(kxx)];
                    y_ = p.Property_subdomains_statistics(:,5,current_phase_todo);
                    y_ = [y_; p.Wholevolume_results(current_phase_todo)];
                elseif id_axe==2
                    x_=p.relativedifference_convergence(:,1);
                    y_=p.relativedifference_convergence(:,current_phase_todo+2);
                    x_conv = [x_(1) x_(end)];
                else
                    x_=p.relativedifference_convergence(:,2);
                    y_=p.relativedifference_convergence(:,current_phase_todo+2);  
                    x_conv = [x_(1) x_(end)];
                end
                h_=plot(x_,y_);
                % Colors, thickness, markers
                set(h_, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o');
            end
        end

        % Other plot with no legend
        if sum(p.Size_RVE(:,current_phase_todo,2,kxx)~=0)==1 % Only plot if only one threshold has been selected
            if id_axe==2 || id_axe==4
                for k_threshold=1:1:length(p.convergence_criterion)
                    if p.Size_convergence(k_threshold,current_phase_todo,2,kxx)~=0
                        y_conv = [p.convergence_criterion(k_threshold) p.convergence_criterion(k_threshold)];
                        h_conv = plot(x_conv,y_conv);
                        set(h_conv, 'Color', 'black','LineStyle','--','LineWidth',p.opt.format.linewidth);

                        h=plot([p.Size_convergence(k_threshold,current_phase_todo,2,kxx) p.Size_convergence(k_threshold,current_phase_todo,2,kxx)],[0 p.convergence_criterion(k_threshold)]); % Curves
                        set(h, 'Color', 'k','LineStyle','--','LineWidth',p.opt.format.linewidth); % Colors
                    end
                end
            end
        end

        h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
        h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
        % - Axis label
        xlabel(str_xlabel);
        if id_axe==1 || id_axe==3
            if isempty(p.propertynameunit)
                ylabel(p.propertyname);
            else
                ylabel([p.propertyname ' (' p.propertynameunit ')']);
            end
        else
            ylabel('%');
        end
        grid(sub_axes,p.opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
        h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    if strcmp(p.RVE.type,'E') || strcmp(p.RVE.type,'F')
        ti = {['Convergence analysis:' p.propertyname],p.RVE.choice_Onsub,num2str(p.RVE.Aspectratio_Onesub)};
    else
        ti = {['Convergence analysis:' p.propertyname],p.RVE.choice_Onsub,p.RVE.Growthdirection};
    end

    sgtitle(Fig,ti,'FontWeight','bold','FontSize',p.opt.format.sgtitlefontsize,'FontName',p.opt.format.fontname);
    if p.opt.save.savefig % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
    end
    if p.opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end


    % Convergence length as function of threshold
    if sum(sum(p.Size_convergence(:,:,2,1)~=0) >= 2) % At least 2 points
        Fig = figure; % Create figure
        Fig.Name= [p.propertyname ' Convergence length as function of relative difference threshold']; % Figure name
        Fig.Color='white'; % Background colour
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*scr2 scrsz(4)*1/2]); % Full screen figure
        for id_axe=1:1:n2_axe % Iterate over axe
            if id_axe==1
                str_ylabel=['Equivalent length to reach convergence, cubic root ' axisunit];
                kx = 1;
            else
                if strcmp(p.RVE.type,'G')
                    str_ylabel=['Length to reach convergence ' axisunit];
                elseif strcmp(p.RVE.type,'H')
                    str_ylabel=['Square root length to reach convergence ' axisunit];
                end
                kx = 2;
            end

            sub_axes=subplot(1,n2_axe,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            str1 = [p.propertyname ' convergence size'];
            str2 = p.RVE.choice_Onsub;
            str_title = {str1,str2,};
            h_title=title (str_title,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.titlefontsize); % Set title font and string
            % - Plot graphs
            current_phase_todo = 0;
            for current_phase=1:1:p.number_phase
                if p.todo(current_phase)
                    current_phase_todo=current_phase_todo+1;
                    str_legend(current_phase_todo).name = char(p.infovol.phasename(current_phase));
                    x_= p.convergence_criterion;
                    y_ = p.Size_convergence(:,current_phase_todo,2,kx);
                    y_(y_==0)=NaN;
                    h=plot(x_,y_); % Curves
                    set(h, 'Color', p.infovol.phasecolor(current_phase,:),'LineWidth',p.opt.format.linewidth,'MarkerSize',p.opt.format.axefontsize,'Marker','o'); % Colors
                end
            end            
            h_legend = legend(sub_axes,str_legend.name,'Location','best'); % Legend
            h_legend.FontSize = p.opt.format.legendfontsize; % Set fontsize
            % - Axis label
            xlabel('Relative difference threshold (%)');
            ylabel(str_ylabel);
            grid(sub_axes,p.opt.format.grid); % Display grid
            set(sub_axes,'XMinorGrid',p.opt.format.minorgrid,'YMinorGrid',p.opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',p.opt.format.fontname,'FontSize',p.opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = p.opt.format.titlefontsize; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
        if p.opt.save.savefig % Save figure
            filename= [p.propertyname '_fctThreshold_' p.RVE.savename];
            function_savefig(Fig, p.savefolder, filename, p.opt.save); % Call function
        end
        if p.opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end


end

end

