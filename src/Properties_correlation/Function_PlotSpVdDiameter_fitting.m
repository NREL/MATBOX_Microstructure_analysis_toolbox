function [ ] = Function_PlotSpVdDiameter_fitting(eps, pv, pd, sp, spunit, polyfit_order, names, plotper, volume_choice, group_choice, form, pformat)

% fit Sp=a * vf + b
% fit Sp=a * eps + b;

%% PRE-PROCESS DATA (CORRELATION)
% Group id
if strcmp(plotper,'group')
    nvol = 0;
    groupnames = {}; kg=0;
    group_id = cell2mat(group_choice(:,1));
    volume_group = cell2mat(volume_choice(:,3));
    %number_group = length(unique(volume_group));
    number_group = length(group_id);    
    for kgr = 1:1:number_group
        idg = group_id(kgr);
        idv = find(volume_group==idg);
        if ~isempty(idv)

            epspvpdsp_v = [eps(idv) pv(idv) pd(idv) sp(idv)];
            names_v = names(idv);

            idnan = find(isnan(epspvpdsp_v));
            if ~isempty(idnan)
                [rownan,~]=ind2sub(size(epspvpdsp_v),idnan);
                epspvpdsp_v(rownan,:) = [];
                names_v(rownan) = [];
            end

            if ~isempty(epspvpdsp_v) % We have at least one {porosity,particle volume fraction, particle diameter, specific surface area} quadruplet
                res(kgr).groupname = group_choice(kgr,2);
                res(kgr).volname = names_v;
                res(kgr).avgratio_3overR = mean( 3 ./ (0.5*epspvpdsp_v(:,3)) ); % <3/R>
                if pformat.thlimit
                    epspvpdsp_v = [epspvpdsp_v; 0 1 0 0; 1 0 0 0];
                end                
                res(kgr).vals = epspvpdsp_v;
                kg = kg+1; groupnames(kg) = res(kgr).groupname;

                % Fit
                if strcmp(form,'Porosity')
                    x = epspvpdsp_v(:,1);
                else
                    x = epspvpdsp_v(:,2);
                end
                y = epspvpdsp_v(:,4);
                res(kgr).fit = polyfit(x,y,min([polyfit_order length(x)-1]));
                nvol = max([nvol, length(names_v)]);
            end
        end
    end
    number_group = length(res);
else
    epspvpdsp_v = [eps pv pd sp];
    names_v = names;

    idnan = find(isnan(epspvpdsp_v));
    if ~isempty(idnan)
        [rownan,~]=ind2sub(size(epspvpdsp_v),idnan);
        epspvpdsp_v(rownan,:) = [];
        names_v(rownan) = [];
    end

    if ~isempty(epspvpdsp_v) % We have at least one {porosity,particle volume fraction, particle diameter, specific surface area} quadruplet
        res(1).volname = names_v;
        res(1).avgratio_3overR = mean( 3 ./ (0.5*epspvpdsp_v(:,3)) ); % <3/R>
        if pformat.thlimit
            epspvpdsp_v = [epspvpdsp_v; 0 1 0 0; 1 0 0 0];
        end    
        res(1).vals = epspvpdsp_v;
        % Fit
        if strcmp(form,'Porosity')
            x = epspvpdsp_v(:,1);
        else
            x = epspvpdsp_v(:,2);
        end
        y = epspvpdsp_v(:,4);
        res(1).fit = polyfit(x,y,min([polyfit_order length(x)-1]));
    end

    [number_volume,~] = size(epspvpdsp_v);
    if pformat.thlimit
        number_volume=number_volume-2;
    end    

end

%% PLOT
markerlist = {'o','square','diamond','^','v','<','pentagram','hexagram','+','*','.','x','_','|'};
markerlist = [markerlist markerlist markerlist markerlist markerlist];
markerlist = [markerlist markerlist markerlist markerlist markerlist];
colorlist = [colororder; rand(100,3)];

Fig = figure; % Create figure
Fig.Name= ['Specific surface area - particle diameter and volume fraction, correlation per ' plotper];
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.75*scrsz(3) 0.5*scrsz(4)]);

t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
title(t,'Specific surface area, particle diameter and volume fraction correlation','FontWeight','normal','Fontname',pformat.fontname,'Fontsize',pformat.size.sgtitle_main)
t.Subtitle.String = ['Per ' plotper];
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = pformat.fontname;
t.Subtitle.FontSize = pformat.size.sgtitle_sub;

if strcmp(form,'Porosity')
    c_vf = 1;
    str_xlabel = 'Porosity \epsilon';
    fit_x = '\epsilon';
else
    c_vf = 2;
    str_xlabel = 'Particle volume fraction \epsilon_{s}';
    fit_x = '\epsilon_{s}';    
end

if strcmp(spunit,'um-1') || strcmp(spunit,'micrometer-1') || strcmp(spunit,'Micrometer-1') || strcmp(spunit,'micrometers-1') || strcmp(spunit,'Micrometers-1')
    unit_name = '(\mum^{-1})';
elseif strcmp(spunit,'nm-1') || strcmp(spunit,'nanometer-1') || strcmp(spunit,'Nanometer-1') || strcmp(spunit,'nanometers-1') || strcmp(spunit,'Nanometers-1')
    unit_name = '(nm^{-1})';
else
    unit_name = ['(' num2str(cell2mat(spunit)) ')'];
end

for col=1:1:3
    nexttile

    if col==1
        if strcmp(plotper,'volume') % Plot per volume
            hold on;
            colvol = [];
            for kvol = 1:1:number_volume
                x = res(1).vals(kvol,c_vf);
                y = res(1).vals(kvol,4);
                colvol = [colvol; colorlist(kvol,:)]; % To make sure color is synchronized with second plot
                plot(x,y,"LineStyle","none","Marker",markerlist(kvol),"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",colorlist(kvol,:),'DisplayName',char(names_v(kvol)));
                if pformat.text
                    text(x,y,char(names_v(kvol)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);                    
                end
            end

            if pformat.thlimit
                plot([0 1],[0 0],"LineStyle","none","Marker",'o',"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",'k','DisplayName','Theoritical limit');
            end

            % fit
            min_x = min(res(1).vals(:,c_vf));
            max_x = max(res(1).vals(:,c_vf));
            x_fit = linspace(min_x,max_x,1000);
            y_fit = polyval(res(1).fit,x_fit);
            str_fit = [];
            n_fit = length(res(1).fit)-1;
            for k=1:1:n_fit+1
                str_fit = [str_fit num2str(res(1).fit(k),'%1.3f') fit_x '^{' num2str(n_fit-k+1,'%i') '}'];
                if k<n_fit+1
                    if res(1).fit(k+1)>=0
                        str_fit = [str_fit ' +'];
                    else
                        str_fit = [str_fit ' '];
                    end
                end
            end
            plot(x_fit,y_fit,"LineStyle","--","LineWidth",pformat.linewidth,"Color",'k','DisplayName',str_fit);
      
            % Reference line
            if pformat.referenceline
                if strcmp(form,'Porosity')
                    vs = epspvpdsp_v(:,2); x = epspvpdsp_v(:,1);
                else
                    vs = 0:0.01:1; x = vs;
                end
                Sp_theoritical = res(1).avgratio_3overR .* vs;
                plot(x,Sp_theoritical,"LineStyle",':',"LineWidth",2,"Color",'k','DisplayName','<3/R> \times \epsilon_{s}');
            end 

        elseif strcmp(plotper,'group') % Plot per group
            hold on;
            for kgr = 1:1:number_group
                xy = res(kgr).vals;
                if cell2mat(group_choice(kgr,8))
                    plot(xy(:,c_vf),xy(:,4),"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"MarkerFaceColor",str2num(cell2mat(group_choice(kgr,3))),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',char(res(kgr).groupname));
                else
                    plot(xy(:,c_vf),xy(:,4),"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',char(res(kgr).groupname));
                end
                if pformat.text
                    for kv=1:1:length(res(kgr).volname)
                        x = xy(kv,c_vf);
                        y = xy(kv,4);
                        text(x,y,char(res(kgr).volname(kv)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);
                    end
                end

                % fit
                min_x = min(xy(:,c_vf));
                max_x = max(xy(:,c_vf));
                x_fit = linspace(min_x,max_x,1000);
                y_fit = polyval(res(kgr).fit,x_fit);
                str_fit = [];
                n_fit = length(res(kgr).fit)-1;
                for k=1:1:n_fit+1
                    str_fit = [str_fit num2str(res(kgr).fit(k),'%1.3f') fit_x '^{' num2str(n_fit-k+1,'%i') '}'];
                    if k<n_fit+1
                        if res(kgr).fit(k+1)>=0
                            str_fit = [str_fit ' +'];
                        else
                            str_fit = [str_fit ' '];
                        end
                    end
                end
                plot(x_fit,y_fit,"LineStyle","--","LineWidth",pformat.linewidth,"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',str_fit);

                % Reference line
                if pformat.referenceline
                    if strcmp(form,'Porosity')
                        vs = xy(:,2); x = xy(:,1);
                    else
                        vs = 0:0.01:1; x = vs;
                    end
                    Sp_theoritical = res(kgr).avgratio_3overR .* vs;
                    plot(x,Sp_theoritical,"LineStyle",':',"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName','<3/R> \times \epsilon_{s}');

                end

            end

        end

        xlabel(str_xlabel);
        ylabel(['Specific surface area S_{p} ' unit_name]);
        set(gca,'Fontname',pformat.fontname,'Fontsize',pformat.size.axes)
        grid(gca,pformat.grid); % Display grid
        set(gca,'XMinorGrid',pformat.minorgrid,'YMinorGrid',pformat.minorgrid);
        h_legend = legend('Location','best');
        h_legend.FontSize = pformat.size.legend;

    elseif col==2

        if strcmp(plotper,'volume') % Plot per volume
            hold on;
            colvol = [];
            for kvol = 1:1:number_volume
                x = res(1).vals(kvol,c_vf);
                y = res(1).vals(kvol,4)*res(1).vals(kvol,3)/2; % R*Sp
                colvol = [colvol; colorlist(kvol,:)]; % To make sure color is synchronized with second plot
                plot(x,y,"LineStyle","none","Marker",markerlist(kvol),"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",colorlist(kvol,:),'DisplayName',char(names_v(kvol)));
                if pformat.text
                    text(x,y,char(names_v(kvol)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);
                end
            end

            if pformat.thlimit
                plot([0 1],[0 0],"LineStyle","none","Marker",'o',"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",'k','DisplayName','Theoritical limit');
            end

            % fit
            x = res(1).vals(:,c_vf);
            y = res(1).vals(:,4).*res(1).vals(:,3)/2; % R*Sp            
            p =  polyfit(x, y,min([polyfit_order length(x)-1]));
            x_fit = linspace(min(x),max(x),1000);            
            y_fit = polyval(p,x_fit);
            str_fit = [];
            n_fit = length(p)-1;
            for k=1:1:n_fit+1
                str_fit = [str_fit num2str(p(k),'%1.3f') fit_x '^{' num2str(n_fit-k+1,'%i') '}'];
                if k<n_fit+1
                    if p(k+1)>=0
                        str_fit = [str_fit ' +'];
                    else
                        str_fit = [str_fit ' '];
                    end
                end
            end
            plot(x_fit,y_fit,"LineStyle","--","LineWidth",pformat.linewidth,"Color",'k','DisplayName',str_fit);

            % Reference line
            if pformat.referenceline
                if strcmp(form,'Porosity')
                    vs = epspvpdsp_v(:,2); x = epspvpdsp_v(:,1); y = 3*vs;
                else
                    vs = 0:0.01:1; x = vs; y = 3*x;
                end
                plot(x,y,"LineStyle",':',"LineWidth",2,"Color",'k','DisplayName','3 \times \epsilon_{s}');
            end

        elseif strcmp(plotper,'group') % Plot per group
            hold on;
            for kgr = 1:1:number_group
                xy = res(kgr).vals;
                if cell2mat(group_choice(kgr,8))
                    plot(xy(:,c_vf),xy(:,3).*xy(:,4)/2,"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"MarkerFaceColor",str2num(cell2mat(group_choice(kgr,3))),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',char(res(kgr).groupname));
                else
                    plot(xy(:,c_vf),xy(:,3).*xy(:,4)/2,"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',char(res(kgr).groupname));
                end
                if pformat.text
                    for kv=1:1:length(res(kgr).volname)
                        x = xy(kv,c_vf);
                        y = xy(kv,3).*xy(kv,4)/2;
                        text(x,y,char(res(kgr).volname(kv)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);
                    end
                end

                % fit
                x = xy(:,c_vf);
                y = xy(:,3).*xy(:,4)/2; % R*Sp            
                p =  polyfit(x, y,min([polyfit_order length(x)-1]));
                x_fit = linspace(min(x),max(x),1000);    
                y_fit = polyval(p,x_fit);
                str_fit = [];
                n_fit = length(p)-1;
                for k=1:1:n_fit+1
                    str_fit = [str_fit num2str(p(k),'%1.3f') fit_x '^{' num2str(n_fit-k+1,'%i') '}'];
                    if k<n_fit+1
                        if p(k+1)>=0
                            str_fit = [str_fit ' +'];
                        else
                            str_fit = [str_fit ' '];
                        end
                    end
                end
                plot(x_fit,y_fit,"LineStyle","--","LineWidth",pformat.linewidth,"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',str_fit);

                % Reference line
                if pformat.referenceline
                    if strcmp(form,'Porosity')
                        vs = xy(:,2); x = xy(:,1);  y = 3*vs;
                    else
                        vs = 0:0.01:1; x = vs;  y = 3*x;
                    end
                    plot(x,y,"LineStyle",':',"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName','3 \times \epsilon_{s}');

                end

            end

        end

        xlabel(str_xlabel);
        ylabel('Specific surface area S_{p} \times particle radius R_{s}');
        set(gca,'Fontname',pformat.fontname,'Fontsize',pformat.size.axes)
        grid(gca,pformat.grid); % Display grid
        set(gca,'XMinorGrid',pformat.minorgrid,'YMinorGrid',pformat.minorgrid);
        h_legend = legend('Location','best');
        h_legend.FontSize = pformat.size.legend;


    else

        if strcmp(plotper,'volume') % Plot per volume
            X = categorical(names_v);
            X = reordercats(X,names_v);
            hold on
            for kv=1:1:number_volume
                yy = NaN(number_volume,1); % Creating as many bar object as group so that i can customize color
                yy(kv,1) = res(1).vals(kv,3)*res(1).vals(kv,4)/(2*res(1).vals(kv,2)); % R*Sp/solid vol fraction
                b=bar(X,yy);
                b(1).FaceColor = colvol(kv,:);
            end
            n = number_volume;

        elseif strcmp(plotper,'group') % Plot per group
            X = categorical(groupnames);
            X = reordercats(X,groupnames);
            y = NaN(number_group,nvol);
            for kgr=1:1:number_group
                %y(kgr,:) = res(kgr).vals(1:nvol,3).*res(kgr).vals(1:nvol,4) ./ (2.*res(kgr).vals(1:nvol,2));
                current_nvol = size(res(kgr).vals);
                current_nvol = current_nvol(1)-2;
                for kv=1:1:current_nvol
                    y(kgr,kv) = res(kgr).vals(kv,3).*res(kgr).vals(kv,4) ./ (2.*res(kgr).vals(kv,2));
                end
            end
            hold on
            for kgr=1:1:number_group
                yy = NaN(number_group,nvol); % Creating as many bar object as group so that i can customize color
                yy(kgr,:) = y(kgr,:);
                b=bar(X,yy);
                for kv=1:1:length(res(kgr).volname)
                    b(kv).FaceColor = str2num(cell2mat(group_choice(kgr,3)));
                end
                if pformat.text
                    for kv=1:1:length(res(kgr).volname)
                        xtips1 = b(kv).XEndPoints;
                        ytips1 = b(kv).YEndPoints;
                        %labels1 = string(b(1).YData);
                        labels1 = char(res(kgr).volname(kv));
                        text(xtips1,ytips1,labels1,'HorizontalAlignment','left','VerticalAlignment','middle','Rotation',90,'Fontname',pformat.fontname,'Fontsize',pformat.textsize)
                    end
                end                
            end
            n = number_group;
        end

        % Reference line
        if pformat.referenceline
            plot([1 n], [3 3],"LineStyle",'--',"LineWidth",1,"Color",'k');
            text(n,3,'Spheres','Fontname',pformat.fontname,'Fontsize',pformat.textsize)
        end

        ylabel('R_{s} \times S_{p} / \epsilon_{s}');
        set(gca,'Fontname',pformat.fontname,'Fontsize',pformat.size.axes)
        grid(gca,pformat.grid); % Display grid
        set(gca,'XMinorGrid',pformat.minorgrid,'YMinorGrid',pformat.minorgrid);
    end

end


end