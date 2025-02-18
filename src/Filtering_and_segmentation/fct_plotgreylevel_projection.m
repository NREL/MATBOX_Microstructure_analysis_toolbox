function [] = fct_plotgreylevel_projection(M,p)

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
    Fig.Name= '2D projection maps';
    p.filename = '2D projection maps';
else
    Fig.Name= ['2D projection maps, ' p.filename];
end
if ~isfield(p,'sub')
    p.sub = [];
end
if ~isfield(p,'colormap')
    p.colormap = 'turbo';
end
if ~isfield(p,'perspective')
    p.perspective = '2D';
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

scrsz = get(0,'ScreenSize'); % Screen resolution

Fig.Color='white'; % Background colour

if dimension==3
    M = double(M);
    if p.excludebackground
        if strcmp(p.selectbackground,'From label')
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
    if ~isempty(str_sub)
        subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);
    end

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
        ylabel(h,['Mean ' p.xlabel]);
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

if isfield(p,'fullpath')
    [filepath,~,~] = fileparts(p.fullpath);
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    savefig(Fig,p.fullpath)
    saveas(Fig,p.fullpath,'png');
end

end