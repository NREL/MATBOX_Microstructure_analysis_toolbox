function [] = fct_plothistogram_wholevolume(M,p)

Fig = figure;

%% Default options
if nargin==1
    p = [];
end
% Set default option
if ~isfield(p,'linewidth')
    p.linewidth = 2;
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
    Fig.Name= 'Histogram for the whole volume';
    p.filename = 'Histogram of the whole array';
else
    Fig.Name= ['Histogram for the whole volume' ', ' p.filename];
end
if ~isfield(p,'sub')
    p.sub = p.option;
end

if ~isfield(p,'xlabel')
    p.xlabel = 'Grey level';
end


%% Plot
if p.excludebackground
    str_sub = [p.sub ' (background excluded)'];
else
    str_sub = p.sub;
end

p.filename = strrep(p.filename,'_',' ');
str_sub = strrep(str_sub,'_',' ');

Fig.Color='white'; % Background colour

if p.excludebackground
    if strcmp(p.selectbackground,'From label')
        M(M==p.background_label) = [];
    else
        M(p.background_volume==1) = [];
    end
end

ax = axes('Parent',Fig);
hold(ax,'on');
if strcmp(p.option,'Bins number: one per grey level')
    nbins = max(max(max(M))) - min(min(min(M))) + 1;
    if nbins~=round(nbins)
        nbins = length(unique(M));
    end
    h=histogram(M,nbins);
elseif strcmp(p.option,'Bins number: custom')
    h=histogram(M,p.Nbin);
elseif strcmp(p.option,'Bins number: auto')
    h=histogram(M);
elseif strcmp(p.option, 'Distribution function (discrete)') || strcmp(p.option, 'Distribution function (continuous)')
    M = reshape(M,[numel(M),1]);
    x = min(M):1:max(M);
    
    unis = unique(M);
    if sum(unis==round(unis)) ~= length(unis)
        M = round(M,1);
        x=min(M):0.1:max(M);
    end
    try
        y = sum(M==x); % Parallelized, RAM heavy
    catch ME
        if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
            for k=1:1:numel(x)
                y(k) = sum(M==x(k));
            end
        end
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
    if strcmp(p.option, 'Distribution function (discrete)')
        y = y/sum(y);
    else
        y = y/trapz(double(x),double(y));
    end    
    h = plot(x,y,'LineWidth',p.linewidth);
end
if strcmp(p.option, 'Distribution function (discrete)')
    ylabel('Probability (\Sigma=1)');
elseif strcmp(p.option, 'Distribution function (continuous)')
    ylabel('Probability (\int=1)');
else
    set(h,'LineStyle','none','Normalization','probability');
    ylabel('Probability');
end
xlabel(p.xlabel);
grid(ax,'on');
if strcmp(p.option, 'Distribution function (discrete)') || strcmp(p.option, 'Distribution function (continuous)')
    axis tight;
    xlim([0 max(x)]);
end
set(ax,'FontName',p.fontname,'FontSize',12);
title(ax,p.filename,'FontName',p.fontname,'FontSize',14);
subtitle(ax,str_sub,'FontName',p.fontname,'FontSize',12);
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