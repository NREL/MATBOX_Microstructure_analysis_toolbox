function [] = fct_plotgreylevel_perradialposition(M,p)

Fig = figure;

%% Default options
if nargin==1
    p = [];
end
% Set default option
if ~isfield(p,'fontname')
    p.fontname = 'Times new roman';
end
if ~isfield(p,'excludebackground')
    p.excludebackground = false;
end
if ~isfield(p,'filename')
    Fig.Name= 'Grey level per radial position';
    p.filename = 'Grey level per radial position';
else
    Fig.Name= ['Grey level per radial position, ' p.filename];
end
if ~isfield(p,'sub')
    p.sub = [];
end

if ~isfield(p,'voxelsize')
    p.voxelsize = 1;
end
if ~isfield(p,'voxelunit')
    p.voxelunit = 'Voxel';
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
if ~isempty(str_sub)
    subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);
end

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
    n = length(unis);

    ax = nexttile(ktile);
    hold(ax,'on');

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

    h_maxdisc = plot([r_discmax r_discmax]*p.voxelsize,[min(y_min) max(y_max)], 'Color', 'k','LineStyle','--','LineWidth',1,'DisplayName','Max full disc');
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
    ylabel(p.xlabel);
    axis tight;
    xlim([0 max(x)]);

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