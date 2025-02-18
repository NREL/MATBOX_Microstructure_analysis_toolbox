function [] = fct_plotallseparability_wholevolume(M,p)

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
    Fig.Name= 'Separability for the whole volume';
    p.filename = 'Separability of the whole array';
else
    Fig.Name= ['Separability for the whole volume' ', ' p.filename];
end
if ~isfield(p,'xlabel')
    p.xlabel = 'Grey level';
end

%% Plot
if p.excludebackground
    str_sub = '(background excluded)';
else
    str_sub = [];
end

p.filename = strrep(p.filename,'_',' ');

Fig.Color='white'; % Background colour

if p.excludebackground
    if strcmp(p.selectbackground,'From label')
        M(M==p.background_label) = [];
    else
        M(p.background_volume==1) = [];
    end
end

% Calculation
m_min = min(min(min( M)));
m_max = max(max(max( M)));
[histogram_M] = function_calculate_histogram(M);

% Otsu's algorithm applied to the whole volume
Otsu_result_wholevolume=Function_otsu_algorithm(histogram_M,p.n_phase);
% Best threshold
all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:p.n_phase);
all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
best_ratio_wholevolume = max(all_ratio_wholevolume);
idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
best_threshold_wholevolume = all_threshold_wholevolume(idx(1),:);
thresholds = [m_min; best_threshold_wholevolume'; m_max];


ax = axes('Parent',Fig);
hold(ax,'on');
if p.n_phase==2

    plot(all_threshold_wholevolume,all_ratio_wholevolume,'LineWidth',2);
    xlabel([p.xlabel ' threshold']);
    ylabel('Separability criterion \eta');
    grid(ax,'on');
    set(ax,'FontName',p.fontname,'FontSize',12);
    title(ax,p.filename,'FontName',p.fontname,'FontSize',14);
    if ~isempty(str_sub)
        subtitle(ax,str_sub,'FontName',p.fontname,'FontSize',12);
    end
    hold(ax,'off');
elseif p.n_phase==3
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*3/5 scrsz(4)*3/5])); % Full screen figure
    x = min(min(all_threshold_wholevolume)):1:max(max(all_threshold_wholevolume));
    x_min = min(x);
    [X,Y] = meshgrid(x,x);
    Z = NaN(size(X));
    for k=1:1:length(all_ratio_wholevolume)
        xs = all_threshold_wholevolume(k,:)-x_min+1;
        Z(xs(2),xs(1)) = all_ratio_wholevolume(k);
    end
    levels = linspace(min(all_ratio_wholevolume),max(all_ratio_wholevolume),20);
    levels(end) = levels(end)-0.001;
    max(all_ratio_wholevolume)
    contour(X,Y,Z,levels,'ShowText','on','LabelFormat','%0.2f')
    xlabel([p.xlabel ' 1st threshold']);
    ylabel([p.xlabel ' 2nd threshold']);
    colormap turbo;
    h = colorbar;
    h.FontName = p.fontname;
    h.FontSize = 12;
    ylabel(h,'Separability criterion \eta');
end

grid(ax,'on');
set(ax,'FontName',p.fontname,'FontSize',12);
title(ax,p.filename,'FontName',p.fontname,'FontSize',14);
if ~isempty(str_sub)
    subtitle(ax,str_sub,'FontName',p.fontname,'FontSize',12);
end
hold(ax,'off');

if isfield(p,'fullpath')
    [filepath,~,~] = fileparts(p.fullpath);
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    savefig(Fig,p.fullpath)
    saveas(Fig,p.fullpath,'png');
end


end