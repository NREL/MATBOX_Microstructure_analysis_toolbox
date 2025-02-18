function [] = fct_particleconnectivity(M,p)

sz = size(M);
dimension = length(sz);
if dimension == 2

    %% Default options
    if nargin==1
        p = [];
    end
    % Set default option
    if ~isfield(p,'fontname')
        p.fontname = 'Times new roman';
    end

    if dimension == 2
        if ~isfield(p,'filename')
            p.filename = 'Particle connectivity';
            Fig_strName = 'Particle connectivity';
        else
            Fig_strName = ['Particle connectivity, ' p.filename];
        end
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

    %% CONNECTIVITY MATRIX
    [Connectivity_particle,~] = Function_particleconnectivity(M);

    %% WATERSHED LINES, CENTROID, GRAPH
    if dimension == 2
        labels = unique(M);
        n_labels = length(labels);

        scrsz = get(0,'ScreenSize'); % Screen resolution
        Fig = figure; % Create figure
        Fig.Color='white'; % Background colour
        Fig.Name = Fig_strName;
        set(Fig,'position',scrsz); % Full screen figure
        ax_ = axes('Parent',Fig);
        hold(ax_,'on');

        lake_id_RGB_color = 0.1 + (0.9-0.1).*rand(1e6,3);
        lake_id_grey_color = sum(lake_id_RGB_color,2)/3;
        slice_color = zeros(sz(1),sz(2),3); % RGB color map
        slice_r = zeros(sz(1),sz(2)); % Red color map
        slice_g = zeros(sz(1),sz(2)); % Green color map
        slice_b = zeros(sz(1),sz(2)); % Blue color map
        slice_grey = zeros(sz(1),sz(2)); % Blue color map
        % Lake, centroid
        lake_centroid = zeros( n_labels-1,2);
        lake_radius = zeros( n_labels-1,1);
        centroid_loc = zeros( n_labels-1,2);
        lake_id = zeros( n_labels-1,1);
        kk=0;
        for k=1:1:n_labels
            idx=find(M==labels(k));
            if labels(k)==0
                slice_r(idx) = 1; slice_g(idx) = 1; slice_b(idx) = 1; % Background
                slice_grey(idx) = 1;
            else
                % Equivalent radius
                area_ = length(idx);
                equivalent_radius = (area_/pi)^(1/2);
                slice_r(idx) = lake_id_RGB_color(k,1); slice_g(idx) = lake_id_RGB_color(k,2); slice_b(idx) = lake_id_RGB_color(k,3); % Background
                slice_grey(idx) = lake_id_grey_color(k);
                % Centroid
                kk=kk+1;
                [IX,IJ,IK] = ind2sub(sz,idx);
                centroid_loc(kk,1) = mean(IX-0.5); centroid_loc(kk,2) = mean(IJ-0.5); centroid_loc(kk,3) = mean(IK-0.5);
                lake_centroid(kk,:) = [centroid_loc(kk,2),centroid_loc(kk,1)];
                lake_radius(kk,1) = equivalent_radius;
                lake_id(kk,1) = M(idx(1));
            end
        end
        % Edges
        BW = zeros(sz);
        BW(M~=0)=1;
        [index_border_phase,~,~,~] = Function_identify_edges(BW);
        slice_r(index_border_phase) = 0; slice_g(index_border_phase) = 0; slice_b(index_border_phase) = 0;
        slice_grey(index_border_phase) = 0;
        % Watershed lines
        % Edge detection
        background = 0;
        edgewithbackground = false;
        [index_border_label,~,~,~] = Function_identify_labelsedges(M, background, edgewithbackground);
        tmp = slice_grey; tmp(index_border_label) = 1;
        slice_color(:,:,1)=tmp; % Attribute RGB color
        tmp = slice_grey; tmp(index_border_label) = 0;
        slice_color(:,:,2)=tmp;
        tmp = slice_grey; tmp(index_border_label) = 0;
        slice_color(:,:,3)=tmp;
        slice_image = image(slice_color,'parent',ax_); % Display the slice
        % Equivalent diameter
        viscircles(lake_centroid, lake_radius,'Color',[0.4660    0.6740    0.1880]','LineWidth',1.5); % Equivalent diamter
        % Particle connection
        for r=2:1:length(labels)
            for c=2:1:r-1
                if Connectivity_particle(r,c)>0
                    xx = [lake_centroid(r-1,1) lake_centroid(c-1,1)];
                    yy = [lake_centroid(r-1,2) lake_centroid(c-1,2)];
                    plot(xx,yy,'Color','b','LineWidth',2);
                end
            end
        end
        % Centroid
        plot(centroid_loc(:,2), centroid_loc(:,1),'+','MarkerSize',8,'Color','r','LineWidth',1.5);
        set(ax_,'YDir','normal')
        title('Particle boundaries, centroids, and equivalent diameters','FontName','Times New Roman','FontSize',14,'Parent',ax_);

        xlabel(ax_,'Voxels'); % Label
        ylabel(ax_,'Voxels');
        axis(ax_,'tight'); % Fit the axes box
        axis(ax_,'equal'); % Aspect ratio is 1:1
        box(ax_,'on')
        set(ax_,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        xlim(ax_,[0, sz(2)]);
        ylim(ax_,[0, sz(1)]);

        hold(ax_,'off'); % Active subplot
    end

    if isfield(p,'fullpath')
        [filepath,~,~] = fileparts(p.fullpath);
        if ~exist(filepath,'dir')
            mkdir(filepath);
        end
        savefig(Fig,p.fullpath)
        saveas(Fig,p.fullpath,'png');

        [filepath,filename,~] = fileparts(p.fullpath);
        if ~exist(filepath,'dir')
            mkdir(filepath);
        end
        filename = [filename ' connectivitymatrix']; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Connectivity';
        DATA_writetable.sheet(1).table=array2table(Connectivity_particle);
        % Save function
        Function_Writetable(filepath,filename,DATA_writetable)
    end

end

end