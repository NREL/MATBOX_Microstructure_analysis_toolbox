function [] = fct_plothistogram_perslice(M,p)

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
    Fig.Name= 'Histogram per slice';
    p.filename = 'Histogram per slice';
else
    Fig.Name= ['Histogram per slice, ' p.filename];
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
    if p.allines
        x = 1:1:sz(dim);
    else
        x = round(linspace(1,sz(dim),p.nlines));
    end

    n = numel(x);
    x_um = x*p.voxelsize;
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
            if strcmp(p.selectbackground,'From label')
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
                plot(xx,yy,'Color',cmap(x(k),:),'LineWidth',p.linewidth);
            else
                z = ones(1,numel(xx))*x(k);
                plot3(xx,yy,z,'Color',cmap(x(k),:),'LineWidth',p.linewidth);
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
    ylabel(h,['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);

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