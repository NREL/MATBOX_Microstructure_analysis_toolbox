function [] = fct_plotseparability_perslice(M,p)

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
    Fig.Name= 'Separability per slice';
    p.filename = 'Separability per slice';
else
    Fig.Name= ['Separability per slice, ' p.filename];
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
    str_sub = '(background excluded)';
else
    str_sub = [];
end

p.filename = strrep(p.filename,'_',' ');

if strcmp(p.voxelunit,'um')
    p.voxelunit = '\mum';
end

n_phases = p.n_min:1:p.n_max;

%min_ = double(min(min(min(M))));
%max_ = double(max(max(max(M))));


for dim=1:1:dimension
    for kp = 1:1:length(n_phases)
        n_phase = n_phases(kp);
        thresholds(dim).phases(kp).vals = zeros(sz(dim),n_phase-1);
    end
end

for dim=1:1:dimension

    x = 1:1:sz(dim);

    best_ratio(dim).vals = zeros(sz(dim),length(n_phases));

    for k=1:1:sz(dim)

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

        [histogram_M] = function_calculate_histogram(sl);
        for kp = 1:1:length(n_phases)
            n_phase = n_phases(kp);
            if ~isempty(sl) && length(unique(sl))>=n_phase
                % Otsu's algorithm applied to the whole volume
                Otsu_result_wholevolume=Function_otsu_algorithm(histogram_M,n_phase);
                % Best threshold
                all_threshold=Otsu_result_wholevolume.allpermuation(:,2:n_phase);
                all_ratio=Otsu_result_wholevolume.allratio;
                best_ratio(dim).vals(k,kp) = max(all_ratio);
                idx = find(all_ratio==best_ratio(dim).vals(k,kp));
                best_threshold = all_threshold(idx(1),:);
                thresholds(dim).phases(kp).vals(k,:) = best_threshold';
            else
                best_ratio(dim).vals(k,kp) = NaN;
                %sl = linspace(min_,max_,n_phase+1);
                %thresholds(dim).phases(kp).vals(k,:) = sl(2:end-1);
                thresholds(dim).phases(kp).vals(k,:) = NaN;
            end
        end

    end

end


scrsz = get(0,'ScreenSize'); % Screen resolution
Fig.Color='white'; % Background colour
set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*3/5])); % Full screen figure
t = tiledlayout(2,dimension);
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t,p.filename,'FontName',p.fontname,'FontSize',14,'FontWeight','bold');
if ~isempty(str_sub)
    subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);
end

col = colororder;

for dim=1:1:dimension
    x = 1:1:sz(dim);
    x_um = x*p.voxelsize;

    ax = nexttile(dim);
    hold(ax,'on');    
    for kp = 1:1:length(n_phases)
        plot(x_um, best_ratio(dim).vals(:,kp),'Color',col(kp,:),'LineWidth',2,'DisplayName',[num2str(n_phases(kp)) ' phases']);
    end
    legend(ax,'Location','best');

    xlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
    ylabel('Maximum separability criterion \eta');
    axis tight;
    xlim([0 max(x_um)]);

    grid(ax,'on');
    set(ax,'FontName',p.fontname,'FontSize',12);
    hold(ax,'off');
end

for dim=1:1:dimension
    x = 1:1:sz(dim);
    x_um = x*p.voxelsize;

    ax = nexttile(dimension+dim);
    hold(ax,'on');
    for kp = 1:1:length(n_phases)
        n_phase = n_phases(kp);
        for k=1:1:n_phase-1
            plot(x_um, thresholds(dim).phases(kp).vals(:,k),'Color',col(kp,:),'LineWidth',2,'DisplayName',[num2str(n_phases(kp)) ' phases, threshold ' num2str(k)]);
        end
    end
    legend(ax,'Location','best');

    xlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
    ylabel([p.xlabel ' thresholds']);
    axis tight;
    xlim([0 max(x_um)]);

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