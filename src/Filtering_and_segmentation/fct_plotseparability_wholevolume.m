function [] = fct_plotseparability_wholevolume(M,p)

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

n_phases = p.n_min:1:p.n_max;
max_separability = zeros(length(n_phases),1);

for k = 1:1:length(n_phases)
    n_phase = n_phases(k);
   
    % Otsu's algorithm applied to the whole volume
    Otsu_result_wholevolume=Function_otsu_algorithm(histogram_M,n_phase);
    % Best threshold
    all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:n_phase);
    all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
    best_ratio_wholevolume = max(all_ratio_wholevolume);
    idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
    best_threshold_wholevolume = all_threshold_wholevolume(idx(1),:);
    thresholds(k).vals = [m_min; best_threshold_wholevolume'; m_max];
    max_separability(k) = best_ratio_wholevolume; 
end

ax = axes('Parent',Fig);
hold(ax,'on');
yyaxis left
plot(n_phases,max_separability,'Color','k','LineStyle','-','LineWidth',2,'Marker','o','MarkerSize',10);
ylabel('Maximum separability criterion \eta');
ax.YColor = 'k';

yyaxis right
for k = 1:1:length(n_phases)
    x = ones(1,length(thresholds(k).vals))*n_phases(k);
    plot(x,thresholds(k).vals,'Color','b','LineStyle','none','LineWidth',2,'Marker','_','MarkerSize',10);
end
ylabel([p.xlabel ' thresholds']);
ax.YColor = 'b';

xlabel('Number of phases');
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