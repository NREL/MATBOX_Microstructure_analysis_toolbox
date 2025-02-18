function [M] = fct_greylevel_analysis(M,p)

sz = size(M);
dimension = length(sz);

if p.excludebackground
    str_sub = [p.sub ' (background excluded)'];
else
    str_sub = p.sub;
end

if strcmp(p.voxelunit,'um')
    p.voxelunit = '\mum';
end

scrsz = get(0,'ScreenSize'); % Screen resolution

Fig = figure;
Fig.Color='white'; % Background colour
Fig.Name= [p.choice ', ' p.filename];

if strcmp(p.choice,'Histogram for the whole volume')
    if p.excludebackground
        if strcmp(p.selectbackground,'from label')
            M(M==p.background_label) = [];
        else
            M(p.background_volume==1) = [];
        end
    end

    ax = axes('Parent',Fig);
    hold(ax,'on');
    if strcmp(p.option,'Bins number: one per grey level')
        nbins = max(max(max(M))) - min(min(min(M))) + 1;
        h=histogram(M,nbins);
    elseif strcmp(p.option,'Bins number: custom')
        h=histogram(M,p.Nbin);
    elseif strcmp(p.option,'Bins number: auto')
        h=histogram(M);
    elseif strcmp(p.option, 'Distribution function (discrete)') || strcmp(p.option, 'Distribution function (continuous)')
        M = reshape(M,[numel(M),1]);
        x = min(M):1:max(M);
        y = sum(M==x);
        if strcmp(p.option, 'Distribution function (discrete)')
            y = y/sum(y);
        else
            y = y/trapz(double(x),double(y));
        end
        if p.movingaverage~=0
            y_smoothed = y;
            for k=1:1:numel(x)
                min_range=max(1,k-p.movingaverage);
                max_range=min(numel(x),k+p.movingaverage);
                y_smoothed(k)= mean(y(min_range:max_range));
            end
            y = y_smoothed;
        end
        h = plot(x,y,'LineWidth',p.width);
    end
    if strcmp(p.option, 'Distribution function (discrete)')
        ylabel('Probability (\Sigma=1)');
    elseif strcmp(p.option, 'Distribution function (continuous)')
        ylabel('Probability (\int=1)');
    else
        set(h,'LineStyle','none','Normalization','probability');
        ylabel('Probability');
    end
    xlabel('Grey level');
    grid(ax,'on');
    if strcmp(p.option, 'Distribution function (discrete)') || strcmp(p.option, 'Distribution function (continuous)')
        axis tight;
        xlim([0 max(x)]);
    end
    set(ax,'FontName',p.fontname,'FontSize',12);
    title(ax,p.filename,'FontName',p.fontname,'FontSize',14);
    subtitle(ax,str_sub,'FontName',p.fontname,'FontSize',12);
    hold(ax,'off');


elseif strcmp(p.choice,'Mean from center (radial evolution)') || strcmp(p.choice,'Radial histogram')
    if dimension==2
        ntile = 1;
    else
        ntile = 3;
        set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*3/5])); % Full screen figure
    end
    t = tiledlayout(1,ntile);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,p.filename,'FontName',p.fontname,'FontSize',14,'FontWeight','bold');
    subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);

    M = double(M);
    if p.excludebackground
        if strcmp(p.selectbackground,'from label')
            idnan = M==p.background_label;
        else
            idnan = p.background_volume==1;
        end
        M(idnan) = NaN;
    else
        idnan=[];
    end

    for ktile=1:1:ntile
        BW=zeros(sz);
        if ktile==1 && dimension==3
            tmp = zeros(sz(2),sz(3));
            tmp(round(sz(2)/2),round(sz(3)/2))=1;
            dmap = bwdist(tmp);
            for x=1:1:sz(1)
                BW(x,:,:) = dmap;
            end
            r_discmax = round(min([sz(2) sz(3)])/2);
            dim = 1;

        elseif ktile==2 && dimension==3
            tmp = zeros(sz(1),sz(3));
            tmp(round(sz(1)/2),round(sz(3)/2))=1;
            dmap = bwdist(tmp);
            for y=1:1:sz(2)
                BW(:,y,:) = dmap;
            end
            r_discmax = round(min([sz(1) sz(3)])/2);
            dim = 2;

        elseif ktile==3 || dimension==2
            tmp = zeros(sz(1),sz(2));
            tmp(round(sz(1)/2),round(sz(2)/2))=1;
            dmap = bwdist(tmp);
            if dimension==2
                BW = dmap;
            else
                for z=1:1:sz(3)
                    BW(:,:,z) = dmap;
                end
            end
            r_discmax = round(min([sz(1) sz(2)])/2);  
            dim = 3;
        end

        BW(idnan) = NaN;
        BW = round(BW);
        unis = unique(BW);
        unis(isnan(unis))=[];
        unis(unis==0)=[];
        r_max = max(unis);
        n = length(unis);

        ax = nexttile(ktile);
        hold(ax,'on');        

        if strcmp(p.choice,'Radial histogram')
            if ~p.allines
                x = round(linspace(1,r_max,p.nlines));
            else
                x = [1:1:r_max];
            end
            x_um = x*p.voxelsize;
            n = length(x);
            eval(['cmap = ' p.colormap '(r_max);']);
            for k=1:1:n
                idx = find(BW==x(k));
                vals = M(idx);
                sl = reshape(vals,[numel(vals),1]);

                if ~isempty(sl)
                    xx = min(sl):1:max(sl);
                    yy = sum(sl==xx);
                    if strcmp(p.option, 'Distribution function (discrete)')
                        yy = yy/sum(yy);
                    else
                        yy = yy/trapz(double(xx),double(yy));
                    end

                    if p.movingaverage~=0
                        yy_smoothed = yy;
                        for kkk=1:1:numel(xx)
                            min_range=max(1,kkk-p.movingaverage);
                            max_range=min(numel(xx),kkk+p.movingaverage);
                            yy_smoothed(kkk)= mean(yy(min_range:max_range));
                        end
                        yy = yy_smoothed;
                    end

                    if strcmp(p.perspective,'2D')
                        h = plot(xx,yy,'Color',cmap(x(k),:),'LineWidth',p.width);
                    else
                        z = ones(1,numel(xx))*x(k);
                        h = plot3(xx,yy,z,'Color',cmap(x(k),:),'LineWidth',p.width);
                    end
                end
            end

            xlabel('Grey level');
            if strcmp(p.option, 'Distribution function (discrete)')
                ylabel('Probability (\Sigma=1)');
            elseif strcmp(p.option, 'Distribution function (continuous)')
                ylabel('Probability (\int=1)');
            end

            if strcmp(p.perspective,'3D')
                zlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
                view(180,-25)
            end

            colormap(cmap); clim([min(x_um) max(x_um)]);

            h=colorbar;
            h.FontName = p.fontname;
            h.FontSize = 12;
            h.Limits = [min(x_um) max(x_um)];
            ylabel(h,['Radial location, normal to axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
            axis tight;

            grid(ax,'on');
            set(ax,'FontName',p.fontname,'FontSize',12);
            hold(ax,'off');

        else

            x = [1:1:r_max]*p.voxelsize;
            y_min = NaN(1,n); y_mean = NaN(1,n); y_max = NaN(1,n); y_std = NaN(1,n);
            for pos=1:1:n
                idx = find(BW==pos);
                vals = M(idx);
                y_min(pos) = min(vals);
                y_mean(pos) = mean(vals);
                y_max(pos) = max(vals);
                if length(unique(vals))>1
                    y_std(pos) = std(double(vals));
                end
            end

            h_max = plot([r_discmax r_discmax],[min(y_min) max(y_max)], 'Color', 'b','LineStyle','--','LineWidth',1,'DisplayName','Max full disc');
            h_mean=plot(x,y_mean, 'Color', 'k','LineWidth',2,'DisplayName','Mean');
            h_min= plot(x,y_min, 'Color', 'b','LineStyle','--','LineWidth',1,'DisplayName','Min');
            h_max= plot(x,y_max, 'Color', 'r','LineStyle','--','LineWidth',1,'DisplayName','Max');

            id_NaN = unique([find(isnan(y_mean)); find(isnan(y_std))]);
            x_tmp = x; y_mean_tmp= y_mean;  y_std_tmp= y_std;
            x_tmp(id_NaN)=[]; y_mean_tmp(id_NaN)=[]; y_std_tmp(id_NaN)=[];

            std_min = (y_mean_tmp - y_std_tmp);
            std_max = (y_mean_tmp + y_std_tmp);
            x2 = [x_tmp, fliplr(x_tmp)];
            inBetween = [std_min, fliplr(std_max)];
            h_=fill(x2, inBetween, [0.5 0.5 0.5],'DisplayName','+/- std');
            set(h_,'LineStyle','none','FaceAlpha',0.25);

            legend(ax)
            xlabel(['Radial location, normal to axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
            ylabel('Grey level');
            axis tight;
            xlim([0 max(x)]);

            grid(ax,'on');
            set(ax,'FontName',p.fontname,'FontSize',12);
            hold(ax,'off');
        end
    end


elseif strcmp(p.choice,'Histogram per slice') || strcmp(p.choice,'Mean along slices')
    set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*3/5])); % Full screen figure

    t = tiledlayout(1,dimension);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,p.filename,'FontName',p.fontname,'FontSize',14,'FontWeight','bold');
    subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);

    for dim=1:1:dimension
        eval(['cmap = ' p.colormap '(sz(dim));']);
   
        ax = nexttile(dim);
        hold(ax,'on');
        if strcmp(p.choice,'Histogram per slice') && ~p.allines
            x = round(linspace(1,sz(dim),p.nlines));
        else
            x = 1:1:sz(dim);
        end
        n = numel(x);
        x_um = x*p.voxelsize;
        y_min = NaN(1,n); y_mean = NaN(1,n); y_max = NaN(1,n); y_std = NaN(1,n);
        for k=1:1:n
            if dim==1
                sl = M(x(k),:,:);
            elseif dim==2
                sl = M(:,x(k),:);
            else
                sl = M(:,:,x(k));
            end
            sl = reshape(sl,[numel(sl),1]);

            if p.excludebackground
                if strcmp(p.selectbackground,'from label')
                    sl(sl==p.background_label) = [];
                else
                    if dim==1
                        sl_background = p.background_volume(x(k),:,:);
                    elseif dim==2
                        sl_background = p.background_volume(:,x(k),:);
                    else
                        sl_background = p.background_volume(:,:,x(k));
                    end
                    sl_background = reshape(sl_background,[numel(sl_background),1]);
                    sl(sl_background==1) = [];
                end
            end

            if strcmp(p.choice,'Histogram per slice')
                if ~isempty(sl)
                    xx = min(sl):1:max(sl);
                    yy = sum(sl==xx);
                    if strcmp(p.option, 'Distribution function (discrete)')
                        yy = yy/sum(yy);
                    else
                        yy = yy/trapz(double(xx),double(yy));
                    end

                    if p.movingaverage~=0
                        yy_smoothed = yy;
                        for kkk=1:1:numel(xx)
                            min_range=max(1,kkk-p.movingaverage);
                            max_range=min(numel(xx),kkk+p.movingaverage);
                            yy_smoothed(kkk)= mean(yy(min_range:max_range));
                        end
                        yy = yy_smoothed;
                    end

                    if strcmp(p.perspective,'2D')
                        h = plot(xx,yy,'Color',cmap(x(k),:),'LineWidth',p.width);
                    else
                        z = ones(1,numel(xx))*x(k);
                        h = plot3(xx,yy,z,'Color',cmap(x(k),:),'LineWidth',p.width);
                    end
                end

            elseif strcmp(p.choice,'Mean along slices')
                if ~isempty(sl)
                    y_min(k) = min(sl);
                    y_mean(k) = mean(sl);
                    y_max(k) = max(sl);
                    if length(unique(sl))>1
                        y_std(k) = std(double(sl));
                    end
                end
            end
        end

        if strcmp(p.choice,'Histogram per slice')
            xlabel('Grey level');
            if strcmp(p.option, 'Distribution function (discrete)')
                ylabel('Probability (\Sigma=1)');
            elseif strcmp(p.option, 'Distribution function (continuous)')
                ylabel('Probability (\int=1)');
            end

            if strcmp(p.perspective,'3D')
                zlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
                view(180,-25)
            end

            colormap(cmap); clim([min(x_um) max(x_um)]);

            h=colorbar;
            h.FontName = p.fontname;
            h.FontSize = 12;
            h.Limits = [min(x_um) max(x_um)];
            ylabel(h,['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
            axis tight;

        elseif strcmp(p.choice,'Mean along slices')
            h_mean=plot(x_um,y_mean, 'Color', 'k','LineWidth',2,'DisplayName','Mean');
            h_min= plot(x_um,y_min, 'Color', 'b','LineStyle','--','LineWidth',1,'DisplayName','Min');
            h_max= plot(x_um,y_max, 'Color', 'r','LineStyle','--','LineWidth',1,'DisplayName','Max');

            id_NaN = unique([find(isnan(y_mean)); find(isnan(y_std))]);
            x_tmp = x_um; y_mean_tmp= y_mean;  y_std_tmp= y_std;
            x_tmp(id_NaN)=[]; y_mean_tmp(id_NaN)=[]; y_std_tmp(id_NaN)=[];

            std_min = (y_mean_tmp - y_std_tmp);
            std_max = (y_mean_tmp + y_std_tmp);
            x2 = [x_tmp, fliplr(x_tmp)];
            inBetween = [std_min, fliplr(std_max)];
            h_=fill(x2, inBetween, [0.5 0.5 0.5],'DisplayName','+/- std');
            set(h_,'LineStyle','none','FaceAlpha',0.25);

            legend(ax)
            xlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
            ylabel('Grey level');
            axis tight;
            xlim([0 max(x_um)]);
        end       
        grid(ax,'on');
        set(ax,'FontName',p.fontname,'FontSize',12);
        hold(ax,'off');
    end

elseif strcmp(p.choice,'2D projection maps (3D volume only)') && dimension==3
    M = double(M);
    if p.excludebackground
        if strcmp(p.selectbackground,'from label')
            M(M==p.background_label) = NaN;
        else
            M(p.background_volume==1) = NaN;
        end
    end

    % Normal to axe 1
    greylevel_normal(1).map =  zeros(sz(2),sz(3));
    for x=1:1:sz(2)
        for y=1:1:sz(3)
            line=reshape(M(:,x,y),[sz(1), 1]);
            greylevel_normal(1).map(x,y)=mean(line,'omitnan');
        end
    end
    
    % Normal to axe 2
    greylevel_normal(2).map =  zeros(sz(1),sz(3));
    for x=1:1:sz(1)
        for y=1:1:sz(3)
            line=reshape(M(x,:,y),[sz(2), 1]);
            greylevel_normal(2).map(x,y)=mean(line,'omitnan');
        end
    end
    
    % Normal to axe 3
    greylevel_normal(3).map =  zeros(sz(1),sz(2));
    for x=1:1:sz(1)
        for y=1:1:sz(2)
            line=reshape(M(x,y,:),[sz(3), 1]);
            greylevel_normal(3).map(x,y)=mean(line,'omitnan');
        end
    end    

    set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*3/5])); % Full screen figure
    t = tiledlayout(1,dimension);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,p.filename,'FontName',p.fontname,'FontSize',14,'FontWeight','bold');
    subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);

    for dim=1:1:3
        ax = nexttile(dim);
        hold(ax,'on');
        %title (sprintf('Averaged grey level, view normal to axe %i',dim),'FontName','Times New Roman','FontSize',16);
        if dim==1
            [Y,X]=meshgrid(1:1:sz(2),1:1:sz(3)); % Create grid
            xlabel('3rd Axis (voxel)'); ylabel('2nd Axis (voxel)');
        elseif dim==2
            [Y,X]=meshgrid(1:1:sz(1),1:1:sz(3));
            xlabel('3rd Axis (voxel)'); ylabel('1st Axis (voxel)');
        elseif dim==3
            [Y,X]=meshgrid(1:1:sz(1),1:1:sz(2));
            xlabel('2nd Axis (voxel)'); ylabel('1st Axis (voxel)');
        end
        surfc(X,Y,greylevel_normal(dim).map','Parent',ax,'LineStyle','none'); % Create surface plot with contour
        colormap(p.colormap);
        h=colorbar('peer',ax); % Create colorbar;
        h.FontName = p.fontname;
        h.FontSize = 12;
        ylabel(h,'Mean grey level');
        set(ax,'FontName',p.fontname,'FontSize',12);
        axis(ax,'equal'); axis(ax,'tight');
        % View
        if strcmp(p.perspective,'2D')
            view(0,90)
        elseif strcmp(p.perspective,'3D')
            view(20,45)
        end
    end




end

end