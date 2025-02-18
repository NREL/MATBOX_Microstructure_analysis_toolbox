function [] = fct_plothistogram_perradialposition(M,p)

Fig = figure;

%% Default options
if nargin==1
    p = [];
end
% Set default option
if ~isfield(p,'linewidth')
    p.linewidth = 1;
end
if ~isfield(p,'fontname')
    p.fontname = 'Times new roman';
end
if ~isfield(p,'movingaverage')
    p.movingaverage = 0;
end
if ~isfield(p,'option')
    p.option = 'Distribution function (continuous)';
end
if ~isfield(p,'excludebackground')
    p.excludebackground = false;
end
if ~isfield(p,'filename')
    Fig.Name= 'Histogram per radial position';
    p.filename = 'Histogram per radial position';
else
    Fig.Name= ['Histogram per radial position, ' p.filename];
end
if ~isfield(p,'sub')
    p.sub = p.option;
end

if ~isfield(p,'voxelsize')
    p.voxelsize = 1;
end
if ~isfield(p,'voxelunit')
    p.voxelunit = 'Voxel';
end
if ~isfield(p,'colormap')
    p.colormap = 'turbo';
end
if ~isfield(p,'perspective')
    p.perspective = '2D';
end
if ~isfield(p,'allines')
    p.allines = true;
end

if ~isfield(p,'xlabel')
    p.xlabel = 'Grey level';
end

%% Plot
sz = size(M);
dimension = length(sz);

if p.excludebackground
    str_sub = [p.sub ' (background excluded)'];
else
    str_sub = p.sub;
end

p.filename = strrep(p.filename,'_',' ');
str_sub = strrep(str_sub,'_',' ');

if strcmp(p.voxelunit,'um')
    p.voxelunit = '\mum';
end

scrsz = get(0,'ScreenSize'); % Screen resolution

Fig.Color='white'; % Background colour

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
    if strcmp(p.selectbackground,'From label')
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

    ax = nexttile(ktile);
    hold(ax,'on');

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
                h = plot(xx,yy,'Color',cmap(x(k),:),'LineWidth',p.linewidth);
            else
                z = ones(1,numel(xx))*x(k);
                h = plot3(xx,yy,z,'Color',cmap(x(k),:),'LineWidth',p.linewidth);
            end
        end
    end

    xlabel(p.xlabel);
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

end

if isfield(p,'fullpath')
    [filepath,~,~] = fileparts(p.fullpath);
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    savefig(Fig,p.fullpath)
    saveas(Fig,p.fullpath,'png');
end

end