function [] = Function_visualize_particle_identification(binary_phase, Label_lake, D_PSD_particle_size, parameters)

scrsz = get(0,'ScreenSize'); % Screen resolution
if ~isfield(parameters, 'save_path')
    desktop=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
    parameters.save_path = desktop;
end

Domain_size=size(binary_phase); % Domain size of the microstructure
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end

Fig = figure; % Create figure
Fig.Name= 'Particle identification'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',scrsz); % Full screen figure
hsg = sgtitle(Fig,'Particle identification','FontName','Times New Roman','FontSize',18);

%% LABEL ID
subaxe{1} = subplot(1,3,1,'Parent',Fig); % Create axes
hold(subaxe{1},'on'); % Active subplot
rnd_label_lake = Label_lake;
currentvalues = unique(Label_lake); tmp = currentvalues; tmp(tmp==0)=[];
for k=1:1:length(currentvalues)
    if currentvalues(k)~=0
        idx = randi(length(tmp));
        rnd_label_lake( Label_lake==currentvalues(k) ) = idx;
        tmp(idx)=[];
    end
end
imagesc(subaxe{1},rnd_label_lake);
custom_colormap = parula; custom_colormap(1,:)=[1 1 1];
colormap(subaxe{1},custom_colormap);    
title(['Particle label (', num2str(length(currentvalues)-1) ')'],'FontName','Times New Roman','FontSize',14,'Parent',subaxe{1});
hold(subaxe{1},'off'); % Active subplot

%% PARTICLE SIZE
subaxe{2} = subplot(1,3,2,'Parent',Fig); % Create axes
hold(subaxe{2},'on'); % Active subplot
imagesc(subaxe{2},D_PSD_particle_size);
custom_colormap = jet; custom_colormap(1,:)=[1 1 1];
colormap(subaxe{2},custom_colormap);
h=colorbar;
ylabel(h, 'Equivalent diameter in voxel length');
set(h,'FontSize',14,'FontName','Times New Roman');
title('Particle size','FontName','Times New Roman','FontSize',14,'Parent',subaxe{2});
hold(subaxe{2},'off'); % Active subplot

%% CALCULATE PARTICLE CONNECTION
[Connectivity_particle,~] = Function_particleconnectivity(Label_lake);

%% WATERSHED LINES, CENTROID, GRAPH
subaxe{3} = subplot(1,3,3,'Parent',Fig); % Create axes
hold(subaxe{3},'on'); % Active subplot
lake_id_RGB_color = 0.1 + (0.9-0.1).*rand(1e6,3);
lake_id_grey_color = sum(lake_id_RGB_color,2)/3;
slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
slice_grey = zeros(Domain_size(1),Domain_size(2)); % Blue color map
% Lake, centroid
lake_centroid = zeros( length(currentvalues)-1,2);
lake_radius = zeros( length(currentvalues)-1,1);
centroid_loc = zeros( length(currentvalues)-1,2);
lake_id = zeros( length(currentvalues)-1,1);
kk=0;
for k=1:1:length(currentvalues)
    idx=find(Label_lake==currentvalues(k));
    if currentvalues(k)==0
        slice_r(idx) = 1; slice_g(idx) = 1; slice_b(idx) = 1; % Background
        slice_grey(idx) = 1;
    else
        slice_r(idx) = lake_id_RGB_color(k,1); slice_g(idx) = lake_id_RGB_color(k,2); slice_b(idx) = lake_id_RGB_color(k,3); % Background
        slice_grey(idx) = lake_id_grey_color(k);
        % Centroid
        kk=kk+1;
        [IX,IJ,IK] = ind2sub(Domain_size,idx);
        centroid_loc(kk,1) = mean(IX-0.5); centroid_loc(kk,2) = mean(IJ-0.5); centroid_loc(kk,3) = mean(IK-0.5);
        lake_centroid(kk,:) = [centroid_loc(kk,2),centroid_loc(kk,1)];
        lake_radius(kk,1) = D_PSD_particle_size(idx(1))/2;
        lake_id(kk,1) = Label_lake(idx(1));
    end
end
% Edges
[index_border_phase,~,~,~] = Function_identify_edges(binary_phase);
slice_r(index_border_phase) = 0; slice_g(index_border_phase) = 0; slice_b(index_border_phase) = 0;
slice_grey(index_border_phase) = 0;
% Watershed lines
% Edge detection
background = 0;
edgewithbackground = false;
[index_border_label,~,~,~] = Function_identify_labelsedges(Label_lake, background, edgewithbackground);
tmp = slice_grey; tmp(index_border_label) = 1;
slice_color(:,:,1)=tmp; % Attribute RGB color
tmp = slice_grey; tmp(index_border_label) = 0;
slice_color(:,:,2)=tmp;
tmp = slice_grey; tmp(index_border_label) = 0;
slice_color(:,:,3)=tmp;
slice_image = image(slice_color,'parent',subaxe{3}); % Display the slice
% Equivalent diameter
viscircles(lake_centroid, lake_radius,'Color',[0.4660    0.6740    0.1880]','LineWidth',1.5); % Equivalent diamter
% Particle connection
for r=2:1:length(currentvalues)
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
set(subaxe{3},'YDir','normal')
title('Particle boundaries, centroids, and equivalent diameters','FontName','Times New Roman','FontSize',14,'Parent',subaxe{3});
hold(subaxe{3},'off'); % Active subplot

%% Both figure
for k=1:1:3
    xlabel(subaxe{k},'Voxels'); % Label
    ylabel(subaxe{k},'Voxels');
    axis(subaxe{k},'tight'); % Fit the axes box
    axis(subaxe{k},'equal'); % Aspect ratio is 1:1
    box(subaxe{k},'on')
    set(subaxe{k},'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
end
xlim(subaxe{3},[0, Domain_size(2)]);
ylim(subaxe{3},[0, Domain_size(1)]);

function_savefig(Fig, [desktop '\'], 'Particle identification'); % Save figure

