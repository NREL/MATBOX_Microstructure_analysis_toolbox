function [] = fct_noisecalculation(M,p)

Fig = figure;

%% Default options
if nargin==1
    p = [];
end
% Set default option
if ~isfield(p,'fontname')
    p.fontname = 'Times new roman';
end

if ~isfield(p,'filename')
    Fig.Name= 'Blind noise estimation per slice';
    p.filename = 'Blind noise estimation per slice';
else
    Fig.Name= ['Blind noise estimation per slice, ' p.filename];
end

if ~isfield(p,'voxelsize')
    p.voxelsize = 1;
end
if ~isfield(p,'voxelunit')
    p.voxelunit = 'Voxel';
end

if strcmp(p.voxelunit,'um')
    p.voxelunit = '\mum';
end

p.filename = strrep(p.filename,'_',' ');

%% Calculation and plot

% modified_img = img + randn(size(img)) * level;
% true noise = level;
% [calculated noise,~] = NoiseLevel(modified_img);
sz=size(M);
dimension = length(sz);
if dimension==3
    array = im2double(M); % Noise function not supported for integer
    for direction=1:1:3
        Allnoise(direction).nlevel=zeros(1,sz(direction));
        for k=1:1:sz(direction)
            if direction==1
                slice_ = squeeze(array(k,:,:));
            elseif direction==2
                slice_ = squeeze(array(:,k,:));
            elseif direction==3
                slice_ = array(:,:,k);
            end
            [nlevel, th] = NoiseLevel(slice_);
            Allnoise(direction).nlevel(1,k) = nlevel / mean(mean(slice_)); % Normalized
        end
    end

    Fig.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*3/5]));

    t = tiledlayout(1,3);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    title(t,p.filename,'FontName',p.fontname,'FontSize',14,'FontWeight','bold');

    for dim=1:1:3 % Iterate over axe
        ax = nexttile(dim);
        hold(ax,'on');

        x = 1:1:sz(dim);
        n = numel(x);
        x_um = x*p.voxelsize;

        grid(ax,'on'); % Display grid
        set(ax,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
        y = Allnoise(dim).nlevel(1,:);
        plot(x,y,'LineWidth',2);
        legend(ax,['Mean noise: ' num2str(sum(y)/length(x))],'Location','best');

        xlabel(['Position along axe ' num2str(dim,'%i') ' (' p.voxelunit  ')' ]);
        ylabel('Normalized noise per slice');
        axis(ax,'tight');
        xlim([0 max(x_um)]);
        set(ax,'FontName',p.fontname,'FontSize',12);
        hold(ax,'off');
    end
    sgtitle(Fig,Fig.Name,'FontWeight','bold','FontSize',16,'FontName',p.fontname);

    if isfield(p,'fullpath')
        [filepath,~,~] = fileparts(p.fullpath);
        if ~exist(filepath,'dir')
            mkdir(filepath);
        end
        savefig(Fig,p.fullpath)
        saveas(Fig,p.fullpath,'png');
    end

else
    array = double(M); % Noise function not supported for integer
    [nlevel, th] = NoiseLevel(array);
    noise = nlevel / mean(mean(array)); % Normalized
    fprintf('Image normalized noise is %e',noise);

    if isfield(p,'fullpath')
        [filepath,~,~] = fileparts(p.fullpath);
        if ~exist(filepath,'dir')
            mkdir(filepath);
        end
        writematrix(noise,p.fullpath,'Delimiter',';')
    end

end