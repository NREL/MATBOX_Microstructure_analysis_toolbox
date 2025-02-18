function [] = fct_volumefractions(M,p)

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
    Fig.Name= 'Volume fractions per slice';
    p.filename = 'Volume fractions per slice';
else
    Fig.Name= ['Volume fractions per slice, ' p.filename];
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


%% GLOBAL
labels = unique(M);
n_labels = length(labels);

if p.excludebackground
    vf = cell(n_labels,3);
else
    vf = cell(n_labels,2);
end

ntot = numel(M);

if p.excludebackground 
    if strcmp(p.selectbackground,'From label')
        n_background = sum(sum(sum( M==p.background_label )));
    else
        idbackground = find(p.background_volume==1);
        Mwoback = M;
        Mwoback(idbackground)=max(labels+1);
        n_background = length(idbackground);
    end
end

for k=1:1:n_labels
    if p.excludebackground && strcmp(p.selectbackground,'From label') && p.background_label == labels(k)
        vf(k,1) = {[num2str(labels(k),'%i') '(background)']};
    else
        vf(k,1) = {labels(k)};
    end

    tmp = sum(sum(sum( M==labels(k) )));

    if p.excludebackground
        if strcmp(p.selectbackground,'From label')
            if p.background_label == labels(k)
                vf(k,3) = {'n/a'};
            else
                vf(k,3) = {tmp / (ntot-n_background)};
            end
        else
            vf(k,3) = {sum(sum(sum( Mwoback==labels(k) ))) / (ntot-n_background)};
        end
    end

    vf(k,2) = {tmp / ntot};
end

if p.excludebackground 
    Table_vf = cell2table(vf,"VariableNames",["Labels" "Volume fractions" "Volume fractions (w/o background)"]);
else
    Table_vf = cell2table(vf,"VariableNames",["Labels" "Volume fractions"]);
end

disp([p.filename ', ' p.sub]);
disp(Table_vf);

if isfield(p,'fullpath')
    [filepath,filename,~] = fileparts(p.fullpath);
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    filename = [filename ' volume fractions']; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Volume fractions';
    DATA_writetable.sheet(1).table=Table_vf;
    % Save function
    Function_Writetable(filepath,filename,DATA_writetable)
end

%% LOCAL PER SLICE

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

col = colororder;
col = [col; rand(n_labels,3)];

sz = size(M);
dimension = length(sz);

scrsz = get(0,'ScreenSize'); % Screen resolution
Fig.Color='white'; % Background colour
set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*3/5])); % Full screen figure

t = tiledlayout(1,dimension);
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t,p.filename,'FontName',p.fontname,'FontSize',14,'FontWeight','bold');
if ~isempty(str_sub)
    subtitle(t,str_sub,'FontName',p.fontname,'FontSize',12);
end

for dim=1:1:dimension
    ax = nexttile(dim);
    hold(ax,'on');

    x = 1:1:sz(dim);
    n = numel(x);
    x_um = x*p.voxelsize;
    vf_loc = NaN(n,n_labels,2);

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
            sl_woback = sl;
            if strcmp(p.selectbackground,'From label')
                sl_woback(sl==p.background_label) = [];
            else
                if dim==1
                    sl_background = p.background_volume(x(k),:,:);
                elseif dim==2
                    sl_background = p.background_volume(:,x(k),:);
                else
                    sl_background = p.background_volume(:,:,x(k));
                end
                sl_background = reshape(sl_background,[numel(sl_background),1]);
                sl_woback(sl_background==1) = [];
            end
        end

        for kl = 1:1:n_labels
            vf_loc(k,kl,1) = sum(sum( sl == labels(kl) )) / numel(sl);
            if p.excludebackground && ~isempty(sl_woback)
                vf_loc(k,kl,2) = sum(sum( sl_woback == labels(kl) )) / numel(sl_woback);
            end
        end
    end

    for kl = 1:1:n_labels
        if p.excludebackground && strcmp(p.selectbackground,'From label') && p.background_label == labels(kl)
            plot(x_um,vf_loc(:,kl,1),'Color', col(kl,:),'LineWidth',2,'DisplayName',['Label ' num2str(labels(kl),'%i') ' (background), <\epsilon>=' num2str(cell2mat(vf(kl,2)),'%1.3f')]);
        else
            plot(x_um,vf_loc(:,kl,1),'Color', col(kl,:),'LineWidth',2,'DisplayName',['Label ' num2str(labels(kl),'%i') ', <\epsilon>=' num2str(cell2mat(vf(kl,2)),'%1.3f')]);
        end
        if p.excludebackground
            if strcmp(p.selectbackground,'From label') && p.background_label == labels(kl)
                foo=1;
            else
                plot(x_um,vf_loc(:,kl,2),'Color', col(kl,:),'LineWidth',2,'LineStyle','--','DisplayName',['Label ' num2str(labels(kl),'%i') ' (w/o background), <\epsilon>=' num2str(cell2mat(vf(kl,3)),'%1.3f')]);
            end
        end
    end

    legend(ax)

    xlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
    ylabel('Volume fractions');
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