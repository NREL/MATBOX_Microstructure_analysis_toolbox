function [] = Function_create_figures_RVE(p)

if strcmp(p.RVE.type,'C')
    r_axe=2; n_axe=6; scr=1;
    n2_axe=2; scr2=2/3;
else
    r_axe=1; n_axe=3; scr=1/2;
    n2_axe=1; scr2=1/3;
end

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter
scrsz = get(0,'ScreenSize'); % Screen resolution

if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B') || strcmp(p.RVE.type,'C')
        
    for current_phase=1:1:p.number_phase % Loop over all phases
        
        Fig = figure; % Create figure
        Fig.Name= ['Representative volume element analysis:' p.propertyname ', ' p.INFO.phase(current_phase).name];
        Fig.Color='white'; % Background colour
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)*scr]); % Full screen figure
        for id_axe=1:1:n_axe % Iterate over axe
            if id_axe<=3
                str_xlabel='Equivalent cubic length of the subdomains (\mum)';
            else
                str_xlabel='Equivalent section length of the subdomains (\mum)';
            end
            sub_axes=subplot(r_axe,3,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            
            if id_axe==1 ||  id_axe==4
                if id_axe==1
                    kx = 1;
                else
                    kx = 2;
                end
                h_title=title ({p.propertyname, 'Statistics'}); % Set title font and string
                % Mean
                x_=p.Property_subdomains_statistics(:,kx,current_phase);
                h_mean=plot(x_,p.Property_subdomains_statistics(:,4,current_phase));
                % Extremums and Whole volume
                h_min=plot(x_,p.Property_subdomains_statistics(:,7,current_phase));
                
                xmin = min(x_); xmax=max(x_);
                h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_phase) p.Wholevolume_results(current_phase)]);
                h_max=plot(x_,p.Property_subdomains_statistics(:,8,current_phase));
                % Colors, thickness, markers
                set(h_mean, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o');
                set(h_min, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'LineStyle','--');
                set(h_max, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'LineStyle','--');
                set(h_whole, 'Color', 'k','LineWidth',p.OPTIONS.Linewidth+1,'LineStyle',':');
                
                % Mean with error bar (+- standard deviation)
                x_tmp = []; y_tmp = []; e_tmp = [];
                for current_size=1:1:length(x_)
                    %if p.Property_subdomains_statistics(current_size,3,current_phase)>=p.Criterion(2) % Standard deviation is meaningless for a low number of values
                    x_tmp=[x_tmp p.Property_subdomains_statistics(current_size,kx,current_phase)];
                    y_tmp=[y_tmp p.Property_subdomains_statistics(current_size,4,current_phase)];
                    e_tmp=[e_tmp p.Property_subdomains_statistics(current_size,5,current_phase)];
                    %end
                end
                h_mean_witherrorbar = errorbar(x_tmp,y_tmp,e_tmp);
                set(h_mean_witherrorbar, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o');
                h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_phase) p.Wholevolume_results(current_phase)]);
                set(h_whole, 'Color', 'k','LineWidth',p.OPTIONS.Linewidth+1,'LineStyle',':');
                % Annotation: number of subdomains
                for current_size=1:1:length(x_)
                    yl = ylim ;
                    y_bottom = p.Property_subdomains_statistics(current_size,4,current_phase)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,4,current_phase)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,3,current_phase) ));
                    set(h,'Color','k','FontSize',p.OPTIONS.Fontsize_axe,'HorizontalAlignment','center','FontWeight','bold','FontName',p.OPTIONS.fontname);
                end
                % - Legend
                h_legend = legend(sub_axes,'Mean with standard deviation','Extremums',['Whole volume, ' num2str(p.Wholevolume_size(kx),'%1.1f') ' \mum'],'Location','best');
                h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
                if isempty(p.propertynameunit)
                    ylabel(p.propertyname);
                else
                    ylabel([p.propertyname ' (' p.propertynameunit ')']);
                end
                
            elseif id_axe==2 || id_axe==5
                if id_axe==2
                    kx = 1;
                else
                    kx = 2;
                end
                h_title=title ({p.propertyname,'Relative standard deviation'}); % Set title font and string
                x_=p.Property_subdomains_statistics(:,kx,current_phase);
                h=plot(x_,p.Property_subdomains_statistics(:,6,current_phase)); % Curves
                set(h, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o'); % Colors
                x_RVE = [x_(1) x_(end)];
                y_RVE = [p.Criterion(1) p.Criterion(1)];
                h_RVE = plot(x_RVE,y_RVE);
                set(h_RVE, 'Color', 'black','LineStyle','--','LineWidth',p.OPTIONS.Linewidth);
                if p.Size_RVE(2,current_phase,kx)~=0
                    h=plot([p.Size_RVE(2,current_phase,kx) p.Size_RVE(2,current_phase,kx)],[0 p.Criterion(1)]); % Curves
                    set(h, 'Color', 'k','LineStyle','--','LineWidth',p.OPTIONS.Linewidth); % Colors
                end
                % Annotation: number of subdomains
                for current_size=1:1:length(x_)
                    yl = ylim ;
                    y_bottom = p.Property_subdomains_statistics(current_size,6,current_phase)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,6,current_phase)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,3,current_phase) ));
                    set(h,'Color','k','FontSize',p.OPTIONS.Fontsize_axe,'HorizontalAlignment','center','FontWeight','bold','FontName',p.OPTIONS.fontname);
                end
                h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
                ylabel('Relative standard deviation (%)');
                
            elseif id_axe==3 || id_axe==6
                if id_axe==3
                    kx = 1; kxcloud = 2;
                else
                    kx = 2; kxcloud = 3;
                end
                h_title=title ({p.propertyname, 'Points cloud'}); % Set title font and string
                % Mean
                x_=p.Property_subdomains_statistics(:,kx,current_phase);
                h_mean=plot(x_,p.Property_subdomains_statistics(:,4,current_phase));
                % Point cloud
                x_= p.Property_eachsubdomain(:,kxcloud);
                y_= p.Property_eachsubdomain(:,current_phase+3);
                h_pointcloud = scatter(x_,y_);
                % Whole volume
                xmin = min(x_); xmax=max(x_);
                h_whole = plot([xmin+0.9*(xmax-xmin) xmax],[p.Wholevolume_results(current_phase) p.Wholevolume_results(current_phase)]);
                % Colors, thickness, markers
                set(h_mean, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o');
                set(h_pointcloud, 'MarkerEdgeColor', p.INFO.phase(current_phase).color);
                set(h_whole, 'Color', 'k','LineWidth',p.OPTIONS.Linewidth+1,'LineStyle',':');
                % Annotation: number of subdomains
                x_=p.Property_subdomains_statistics(:,kx,current_phase);
                for current_size=1:1:length(x_)
                    yl = ylim;
                    y_bottom = p.Property_subdomains_statistics(current_size,4,current_phase)-0.1*(yl(2)-yl(1));
                    y_top = p.Property_subdomains_statistics(current_size,4,current_phase)+0.1*(yl(2)-yl(1));
                    if y_bottom<yl(1)
                        y_=y_top;
                    else
                        y_=y_bottom;
                    end
                    h=text(x_(current_size),y_,num2str( p.Property_subdomains_statistics(current_size,3,current_phase) ));
                    set(h,'Color','k','FontSize',p.OPTIONS.Fontsize_axe,'HorizontalAlignment','center','FontWeight','bold','FontName',p.OPTIONS.fontname);
                end
                h_legend = legend(sub_axes,'Mean','Point cloud',['Whole volume, ' num2str(p.Wholevolume_size(kx),'%1.1f') ' \mum'],'Location','best');
                h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
                if isempty(p.propertynameunit)
                    ylabel(p.propertyname);
                else
                    ylabel([p.propertyname ' (' p.propertynameunit ')']);
                end                
            end
            
            xlabel(str_xlabel);
            if strcmp(p.OPTIONS.grid,'on') % Grid
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
            end
            
            set(sub_axes,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
            h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
            if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
                str_title = {['Representative volume element analysis:' p.propertyname ', ' p.INFO.phase(current_phase).name],p.RVE.name, ['Aspect ratio: ' p.RVE.Aspectratio_name]};
            else
                str_title = {['Representative volume element analysis:' p.propertyname ', ' p.INFO.phase(current_phase).name],p.RVE.name};
            end
            sgtitle(Fig,str_title,'FontWeight','bold','FontSize',p.OPTIONS.Fontsize_title+2,'FontName',p.OPTIONS.fontname);            
        end
        if p.OPTIONS.save_fig == true % Save figure
            filename= [p.propertyname '_' p.INFO.phase(current_phase).name '_' p.RVE.savename];
            function_savefig(Fig, p.savefolder, filename, p.OPTIONS); % Call function
        end
        if p.OPTIONS.closefigureaftercreation == true
            close(Fig); % Do not keep open figures
        end
        
    end
    
    % Relative standard deviation, for all phases
    Fig = figure; % Create figure
    Fig.Name= [p.propertyname ' relative standard deviation']; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*scr2 scrsz(4)*1/2]); % Full screen figure
    for id_axe=1:1:n2_axe % Iterate over axe
        if id_axe==1
            str_xlabel='Equivalent cubic length of the subdomains (\mum)';
            kx = 1;
        else
            str_xlabel='Equivalent section length of the subdomains (\mum)';
            kx = 2;
        end
        sub_axes=subplot(1,n2_axe,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        str1 = [p.propertyname ' relative standard deviation'];
        str2 = p.RVE.name;
        if strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')
            str_title = {str1,str2,['Aspect ratio: ' p.RVE.Aspectratio_name]};
        else
            str_title = {str1,str2,};
        end
        h_title=title (str_title,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_title); % Set title font and string
        % - Plot graphs
        for current_phase=1:1:p.number_phase
            x_=p.Property_subdomains_statistics(:,kx,current_phase);
            h=plot(x_,p.Property_subdomains_statistics(:,6,current_phase)); % Curves
            set(h, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o'); % Colors
        end
        % RVE standard deviation criterion
        x_RVE = [x_(1) x_(end)];
        y_RVE = [p.Criterion(1) p.Criterion(1)];
        h_RVE = plot(x_RVE,y_RVE);
        set(h_RVE, 'Color', 'black','LineStyle','--','LineWidth',p.OPTIONS.Linewidth);
        % Determination of the Representative Volume Elements
        for current_phase=1:1:p.number_phase
            if p.Size_RVE(2,current_phase,kx)~=0
                h=plot([p.Size_RVE(2,current_phase,kx) p.Size_RVE(2,current_phase,kx)],[0 p.Criterion(1)]); % Curves
                set(h, 'Color', p.INFO.phase(current_phase).color,'LineStyle','--','LineWidth',p.OPTIONS.Linewidth); % Colors
            end
        end
        h_legend = legend(sub_axes,p.INFO.phase.name,'Location','best'); % Legend
        h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
        % - Axis label
        xlabel(str_xlabel);
        ylabel('Relative standard deviation (%)');
        if strcmp(p.OPTIONS.grid,'on') % Grid
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    if p.OPTIONS.save_fig == true % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, p.OPTIONS); % Call function
    end
    if p.OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end

   
elseif strcmp(p.RVE.type,'D')
    % Property evolution for all phases
    Fig = figure; % Create figure
    Fig.Name= ['Representative volume element analysis:' p.propertyname];
    Fig.Color='white'; % Background colour
    str_xlabel='Equivalent cubic length of the subdomains (\mum)';
    kx = 1;
    axe_ = axes('Parent',Fig); % Create axes
    hold(axe_,'on');
    h_title=title ({['Representative volume element analysis:' p.propertyname],p.RVE.name,['Aspect ratio: ' p.RVE.Aspectratio_name]},'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_title); % Set title font and string
    for current_phase=1:1:p.number_phase % Loop over all phases
        % Unique volume
        x_=p.Property_subdomains_statistics(:,kx,current_phase);
        x_=[x_; p.Wholevolume_size(kx)];
        y_ = p.Property_subdomains_statistics(:,4,current_phase);
        y_ = [y_; p.Wholevolume_results(current_phase)];
        h_=plot(x_,y_);
        % Colors, thickness, markers
        set(h_, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o');
    end
    h_legend = legend(axe_,p.INFO.phase.name,'Location','best'); % Legend
    h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
    % - Axis label
    xlabel(str_xlabel);
    if isempty(p.propertynameunit)
        ylabel(p.propertyname);
    else
        ylabel([p.propertyname ' (' p.propertynameunit ')']);
    end
    if strcmp(p.OPTIONS.grid,'on') % Grid
        grid(axe_,'on'); % Display grid
        set(axe_,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
    end
    set(axe_,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
    h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
    hold(axe_,'off'); % Relase figure
    if p.OPTIONS.save_fig == true % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, p.OPTIONS); % Call function
    end
    if p.OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
    
elseif strcmp(p.RVE.type,'E')
    % Property evolution for all phases
    Fig = figure; % Create figure
    Fig.Name= ['Representative volume element analysis:' p.propertyname];
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*2/3 scrsz(4)*1/2]); % Full screen figure
    for id_axe=1:1:2 % Iterate over axe
        if id_axe==1
            str_xlabel='Equivalent cubic length of the subdomains (\mum)';
            kx = 1;
        else
            str_xlabel='Length of the subdomains (\mum)';
            kx = 2;
        end
        sub_axes=subplot(1,2,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        
        h_title=title ([p.propertyname ' evolution'],'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_title); % Set title font and string
        for current_phase=1:1:p.number_phase % Loop over all phases
            % Unique volume
            x_=p.Property_subdomains_statistics(:,kx,current_phase);
            x_=[x_; p.Wholevolume_size(id_axe)];
            y_ = p.Property_subdomains_statistics(:,4,current_phase);
            y_ = [y_; p.Wholevolume_results(current_phase)];
            h_=plot(x_,y_);
            % Colors, thickness, markers
            set(h_, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth,'MarkerSize',p.OPTIONS.Fontsize_axe,'Marker','o');
        end
        h_legend = legend(sub_axes,p.INFO.phase.name,'Location','best'); % Legend
        h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
        % - Axis label
        xlabel(str_xlabel);
        if isempty(p.propertynameunit)
            ylabel(p.propertyname);
        else
            ylabel([p.propertyname ' (' p.propertynameunit ')']);
        end
        if strcmp(p.OPTIONS.grid,'on') % Grid
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,{['Representative volume element analysis:' p.propertyname],p.RVE.name,p.RVE.Growthdirection},...
        'FontWeight','bold','FontSize',p.OPTIONS.Fontsize_title+2,'FontName',p.OPTIONS.fontname);
    if p.OPTIONS.save_fig == true % Save figure
        filename= [p.propertyname '_' p.RVE.savename];
        function_savefig(Fig, p.savefolder, filename, p.OPTIONS); % Call function
    end
    if p.OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end


% % Impact of Aspect ratio
% Fig = figure; % Create figure
% Fig.Name= [p.propertyname ' Aspect ratio']; % Figure name
% Fig.Color='white'; % Background colour
% axe_ = axes('Parent',Fig); % Create axes
% hold(axe_,'on');
% str1 = [p.propertyname ' Aspect ratio x:y:1'];
% str2 = p.RVE.name;
% if strcmp(p.RVE.type,'E')
%     str3 = p.RVE.Growthdirection;
%     h_title=title ({str1,str2,str3},'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_title); % Set title font and string
% elseif strcmp(p.RVE.type,'A') || strcmp(p.RVE.type,'B')  || strcmp(p.RVE.type,'D') 
%     str3 = ['Aspect ratio: ' p.RVE.Aspectratio_name];
%     h_title=title ({str1,str2,str3},'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_title); % Set title font and string
% else
%     h_title=title ({str1,str2},'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_title); % Set title font and string
% end
% % Plot
% for current_phase=1:1:p.number_phase
%     AR1 = [p.Wholevolume_size(3); p.Property_subdomains_statistics(:,9,current_phase)]; % Aspect ratio 1
%     AR2 = [p.Wholevolume_size(4); p.Property_subdomains_statistics(:,10,current_phase)]; % Aspect ratio 2
%     sz  = [p.Wholevolume_size(1);  p.Property_subdomains_statistics(:,1,current_phase)]; % subdomain size
%     c   = [p.Wholevolume_results(current_phase); p.Property_subdomains_statistics(:,4,current_phase)]; % Propriety mean value
%     if current_phase==1
%         scatter(AR1,AR2,5*sz,c,'filled','o');
%     elseif current_phase==2
%         scatter(AR1,AR2,5*sz,c,'filled','p');
%     elseif current_phase==3
%         scatter(AR1,AR2,5*sz,c,'filled','d');
%     elseif current_phase==4
%         scatter(AR1,AR2,5*sz,c,'filled','^');
%     elseif current_phase==5
%         scatter(AR1,AR2,5*sz,c,'filled','v');
%     elseif current_phase==6
%         scatter(AR1,AR2,5*sz,c,'filled','>');
%     elseif current_phase==7
%         scatter(AR1,AR2,5*sz,c,'filled','<');
%     end
% end
% xl=xlim;
% yl=ylim;
% scatter(xl(1)+0.9*(xl(2)-xl(1)) ,yl(1)+0.1*(yl(2)-yl(1)),5*100,'k');
% for k=1:1:p.number_phase
%     tmp(k).name = p.INFO.phase(k).name;
% end
% tmp(k+1).name = 'Size reference 100 \mum';
% h_legend = legend(axe_,tmp.name,'Location','best'); % Legend
% h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set fontsize
% h=colorbar(axe_);
% colormap cool
% if isempty(p.propertynameunit)
%     ylabel(h,p.propertyname);
% else
%     ylabel(h,[p.propertyname ' (' p.propertynameunit ')']);
% end
% set(h,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe);
% % - Axis label
% xlabel('Domain aspect ratio x');
% ylabel('Domain aspect ratio y');
% if strcmp(p.OPTIONS.grid,'on') % Grid
%     grid(axe_,'on'); % Display grid
%     set(axe_,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
% end
% set(axe_,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
% h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
% hold(axe_,'off'); % Relase figure
% if p.OPTIONS.save_fig == true % Save figure
%     filename= [p.propertyname '_Aspectratio_' p.RVE.savename];
%     function_savefig(Fig, p.savefolder, filename, p.OPTIONS); % Call function
% end
% if p.OPTIONS.closefigureaftercreation == true
%     close(Fig); % Do not keep open figures
% end


end

