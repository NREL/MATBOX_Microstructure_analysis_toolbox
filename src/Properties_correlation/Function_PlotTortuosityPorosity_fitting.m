function [ ] = Function_PlotTortuosityPorosity_fitting(eps, tau, names, plotper, volume_choice, group_choice, form, pformat)

%% PRE-PROCESS DATA (CORRELATION)
% Group id
if strcmp(plotper,'group')
    nvol = 0;
    groupnames = {}; kg=0;
    group_id = cell2mat(group_choice(:,1));
    volume_group = cell2mat(volume_choice(:,3));
    number_group = length(unique(volume_group));
    for kgr = 1:1:number_group

        idg = group_id(kgr);
        idv = find(volume_group==idg);
        if ~isempty(idv)

            epstau_v = [eps(idv) tau(idv)];
            names_v = names(idv);

            idnan = find(isnan(epstau_v));
            if ~isempty(idnan)
                [rownan,~]=ind2sub(size(epstau_v),idnan);
                epstau_v(rownan,:) = [];
                names_v(rownan) = [];
            end

            if ~isempty(epstau_v) % We have at least one {porosity,tortuosity} doublet
                res(kgr).groupname = group_choice(kgr,2);
                res(kgr).volname = names_v;
                if pformat.thlimit
                    epstau_v = [epstau_v; 1 1];
                end
                res(kgr).vals = epstau_v;
                res(kgr).corr = Function_TortuosityPorosity_fitting(epstau_v(:,1),epstau_v(:,2), form);
                kg = kg+1; groupnames(kg) = res(kgr).groupname;

                % per volume
                nvol = max([nvol, length(names_v)]);
                gammas = zeros(1,length(names_v));
                for kv = 1:1:length(names_v)
                    tmp = Function_TortuosityPorosity_fitting(epstau_v(kv,1),epstau_v(kv,2), 'Bruggeman');
                    gammas(kv) = tmp.gamma;
                end
                res(kgr).gammas = gammas;
            end
        end
    end
    number_group = length(res);
else
    epstau_v = [eps tau];
    names_v = names;

    idnan = find(isnan(epstau_v));
    if ~isempty(idnan)
        [rownan,~]=ind2sub(size(epstau_v),idnan);
        epstau_v(rownan,:) = [];
        names_v(rownan) = [];
    end

    if ~isempty(epstau_v) % We have at least one {porosity,tortuosity} doublet
        res(1).volname = names_v;
        if pformat.thlimit
            epstau_v = [epstau_v; 1 1];
        end
        res(1).vals = epstau_v;
        res(1).corr = Function_TortuosityPorosity_fitting(epstau_v(:,1),epstau_v(:,2), form);
        % per volume
        gammas = zeros(1,length(names_v));
        for kv = 1:1:length(names_v)
            tmp = Function_TortuosityPorosity_fitting(epstau_v(kv,1),epstau_v(kv,2), 'Bruggeman');
            gammas(kv) = tmp.gamma;
        end
        res(1).gammas = gammas;
    end

    [number_volume,~] = size(epstau_v);
    if pformat.thlimit
        number_volume=number_volume-1;
    end

end

%% PLOT
markerlist = {'o','square','diamond','^','v','<','pentagram','hexagram','+','*','.','x','_','|'};
markerlist = [markerlist markerlist markerlist markerlist markerlist];
markerlist = [markerlist markerlist markerlist markerlist markerlist];
colorlist = [colororder; rand(100,3)];

Fig = figure; % Create figure
Fig.Name= ['Tortuosity factor - porosity, correlation per ' plotper];
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.5*scrsz(3) 0.5*scrsz(4)]);

t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
title(t,[form ' tortuosity factor - porosity correlation'],'FontWeight','normal','Fontname',pformat.fontname,'Fontsize',pformat.size.sgtitle_main)
t.Subtitle.String = ['Per ' plotper];
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = pformat.fontname;
t.Subtitle.FontSize = pformat.size.sgtitle_sub;

for col=1:1:2
    nexttile

    if col==1
        if strcmp(plotper,'volume') % Plot per volume
            hold on;
            colvol = [];
            for kvol = 1:1:number_volume
                str_eval = ['strlgd2 = ''\tau = \epsilon^{' num2str(-1*res(1).gammas(kvol),'%1.3f') '}'';'];
                eval(str_eval);
                strlgd1 = [char(names_v(kvol)) ', ' strlgd2];
                x = res(1).vals(kvol,1);
                y = res(1).vals(kvol,2);
                colvol = [colvol; colorlist(kvol,:)]; % To make sure color is synchronized with second plot
                plot(x,y,"LineStyle","none","Marker",markerlist(kvol),"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",colorlist(kvol,:),'DisplayName',strlgd1);
                if pformat.text
                    text(x,y,char(names_v(kvol)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);                    
                end
            end
            if pformat.thlimit
                plot(1,1,"LineStyle","none","Marker",'o',"MarkerSize",pformat.markersize,"LineWidth",pformat.linewidth,"Color",'k','DisplayName','Theoritical limit');
            end

            min_x = min(res(1).vals(:,1));
            max_x = max(res(1).vals(:,1));            
            if strcmp(res(1).corr.equation,'tau = alpha * epsilon^(-gamma)')
                str_eval = ['strlgd3 = ''\tau = ' num2str(res(1).corr.alpha,'%1.3f') ' \times \epsilon^{' num2str(-1*res(1).corr.gamma,'%1.3f') '}'';'];
                eval(str_eval);
            else
                str_eval = ['strlgd3 = ''\tau = \epsilon^{' num2str(-1*res(1).corr.gamma,'%1.3f') '}'';'];
                eval(str_eval);
            end
            if isfield(res(1).corr,'epsilon_fit')
                plot(res(1).corr.epsilon_fit,res(1).corr.tau_fit,"LineStyle",'--',"LineWidth",pformat.linewidth,"Color",'k','DisplayName',strlgd3);
            end

        elseif strcmp(plotper,'group') % Plot per group
            hold on;
            min_x = 1e9; max_x = -1e9;
            for kgr = 1:1:number_group
                xy = res(kgr).vals;
                if strcmp(res(kgr).corr.equation,'tau = alpha * epsilon^(-gamma)')
                    str_eval = ['strlgd2 = ''\tau = ' num2str(res(kgr).corr.alpha,'%1.3f') ' \times \epsilon^{' num2str(-1*res(kgr).corr.gamma,'%1.3f') '}'';'];
                    eval(str_eval);
                else
                    str_eval = ['strlgd2 = ''\tau = \epsilon^{' num2str(-1*res(kgr).corr.gamma,'%1.3f') '}'';'];
                    eval(str_eval);
                end
                strlgd1 = char(res(kgr).groupname);
                if ~isfield(res(kgr).corr,'epsilon_fit')
                    strlgd1 = [strlgd1 ', ' strlgd2];
                end
                if cell2mat(group_choice(kgr,8))
                    plot(xy(:,1),xy(:,2),"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"MarkerFaceColor",str2num(cell2mat(group_choice(kgr,3))),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',strlgd1);
                else
                    plot(xy(:,1),xy(:,2),"LineStyle",char(group_choice(kgr,4)),"Marker",char(group_choice(kgr,6)),"MarkerSize",cell2mat(group_choice(kgr,7)),"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',strlgd1);
                end
                if isfield(res(kgr).corr,'epsilon_fit')
                    plot(res(kgr).corr.epsilon_fit,res(kgr).corr.tau_fit,"LineStyle",'--',"LineWidth",cell2mat(group_choice(kgr,5)),"Color",str2num(cell2mat(group_choice(kgr,3))),'DisplayName',strlgd2);
                end
                min_x = min([min_x min(xy(:,1))]);
                max_x = max([max_x max(xy(:,1))]);
                if pformat.text
                    for kv=1:1:length(res(kgr).volname)
                        x = xy(kv,1);
                        y = xy(kv,2);
                        text(x,y,char(res(kgr).volname(kv)),'Fontname',pformat.fontname,'Fontsize',pformat.textsize);
                    end
                end               
                
            end

        end
        % Reference line
        if pformat.referenceline
            xx = linspace(min_x,max_x,1000);
            tau_ruleofmixture = xx.^(1-1);
            tau_sphere = xx.^(1-1.5);
            tau_cylinders = xx.^(1-2.0);
            if strcmp(plotper,'volume')
                c = [0.5 0.5 0.5];
            else
                c = 'k';
            end
            plot(xx, tau_ruleofmixture,"LineStyle",'-',"LineWidth",2,"Color",c,'DisplayName','Rule of mixture');
            plot(xx, tau_sphere,"LineStyle",'--',"LineWidth",2,"Color",c,'DisplayName','Spheres');
            plot(xx, tau_cylinders,"LineStyle",':',"LineWidth",2,"Color",c,'DisplayName','Cylinders');
        end

        xlabel('Porosity \epsilon');
        ylabel('Tortuosity factor \tau');
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
                yy(kv,1) = 1+res(1).gammas(kv);
                b=bar(X,yy);
                b(1).FaceColor = colvol(kv,:);
            end
            n = number_volume;

        elseif strcmp(plotper,'group') % Plot per group
            X = categorical(groupnames);
            X = reordercats(X,groupnames);
            y = NaN(number_group,nvol);
            for kgr=1:1:number_group
                y(kgr,1:length(res(kgr).gammas)) = 1+res(kgr).gammas;
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
            plot([1 n], [1 1],"LineStyle",'--',"LineWidth",1,"Color",'k');
            plot([1 n], [1.5 1.5],"LineStyle",'--',"LineWidth",1,"Color",'k');
            plot([1 n], [2.0 2.0],"LineStyle",'--',"LineWidth",1,"Color",'k');
            text(n,1,'Rule of mixture','Fontname',pformat.fontname,'Fontsize',pformat.textsize)
            text(n,1.5,'Spheres','Fontname',pformat.fontname,'Fontsize',pformat.textsize)
            text(n,2.0,'Cylinders','Fontname',pformat.fontname,'Fontsize',pformat.textsize)
        end

        ylabel('Bruggeman exponent p, \tau = \epsilon^{1-p}');
        set(gca,'Fontname',pformat.fontname,'Fontsize',pformat.size.axes)
        grid(gca,pformat.grid); % Display grid
        set(gca,'XMinorGrid',pformat.minorgrid,'YMinorGrid',pformat.minorgrid);
    end

end


end