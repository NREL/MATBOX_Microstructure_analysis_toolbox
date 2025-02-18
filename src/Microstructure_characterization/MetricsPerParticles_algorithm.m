clearvars
close all
clc

tmp = function_load_tif('C:\Users\fussegli\Desktop\Articles\tLDRD seg\Dualfilter_Multilayer_segmentation\Paul SEM\Cropped\WF13_Semantic final\S2_Manual reassignment.tif');
Semantic = zeros(size(tmp));
Semantic(tmp==0)=0;
Semantic(tmp==1)=3;
Semantic(tmp==2)=1;
Semantic(tmp==3)=2;

% Label requirements: label 0  must be pore, label 1 must be cracks within particle shells, label 2 must be particles, and subsequent labels can be any other domains (e.g., additives)';
Instance = function_load_tif('C:\Users\fussegli\Desktop\Articles\tLDRD seg\Dualfilter_Multilayer_segmentation\Paul SEM\Cropped\WF6_Particle shells (bis)\S1_Combination (multiplication).tif');

figure; imagesc(Semantic); axis equal; axis tight; colormap gray;
figure; imagesc(Instance); axis equal; axis tight; colormap turbo;

%% PARAMETERS
p.background = 0;
p.voxelsize = 0.5;

p.poreId = 0;
p.cracksId = 1;
p.ParticlesId = 2;

p.adjacency_dist = 5;


%% ALGORITHM
sz = size(Instance);
dimension = length(sz);

Instancelabels = unique(Instance);
Instancelabels(Instancelabels==p.background)=[];
n_instance = length(Instancelabels);

% Initialize
nvoxels = zeros(n_instance,1);
instance_volume = zeros(n_instance,1);
eq_diameter = zeros(n_instance,1);
crack_volume = zeros(n_instance,1);
crack_ratio = zeros(n_instance,1);
n_adjacent = zeros(n_instance,1);
connectivity = [];
centroids = zeros(n_instance,dimension);

adjacency_distum = round(p.adjacency_dist/p.voxelsize);

% Characterize
for ki=1:1:n_instance
    bool = Instance == Instancelabels(ki);

    % Area or volume
    nvoxels(ki) = sum(sum(sum( bool )));
    instance_volume(ki) = nvoxels(ki) * p.voxelsize^dimension;

    % Equivalent diameter
    if dimension == 2
        eq_diameter(ki) = ((3*instance_volume(ki))/(4*pi))^(1/3);
    else
        eq_diameter(ki) = sqrt(instance_volume(ki)/pi);
    end

    % Crack volume and ratio
    idx = find(bool);
    labelsInInstance = Semantic(idx);
    n_cracks = sum(sum(sum(labelsInInstance==p.cracksId)));
    n_particle = sum(sum(sum(labelsInInstance==p.ParticlesId)));
    crack_volume(ki) = n_cracks * p.voxelsize^dimension;
    crack_ratio(ki) = n_cracks/(n_cracks+n_particle);

    % Adjacent particles
    if dimension == 2
        [IX,IY] = ind2sub(sz,idx);
        z_min = 1;
        z_max = 1;
    else
        [IX,IY,IZ] = ind2sub(sz,idx);
        z_min = max([1 min(IZ)-adjacency_distum]);
        z_max = min([sz(3) max(IZ)+adjacency_distum]);
    end

    % Centroids
    centroids(ki,1) = mean(IX);
    centroids(ki,2) = mean(IY);
    if dimension==3
        centroids(ki,3) = mean(IZ);
    end

    x_min = max([1 min(IX)-adjacency_distum]);
    x_max = min([sz(1) max(IX)+adjacency_distum]);
    y_min = max([1 min(IY)-adjacency_distum]);
    y_max = min([sz(2) max(IY)+adjacency_distum]);
    sub = Instance(x_min:x_max,y_min:y_max,z_min:z_max);
    % Find nearby particles
    BW = zeros(size(sub));
    BW(sub==Instancelabels(ki))=1;
    dmap = bwdist(BW);
    sub_adjacent = sub;
    sub_adjacent(dmap>adjacency_distum) = 0;
    sub_adjacent(sub_adjacent==Instancelabels(ki))= 0;
    tmp = unique(sub_adjacent);
    tmp(tmp==0)=[];
    n_adjacent(ki) = length(tmp) - 1;
    connectivity(ki).adjacentinstances = tmp;

    % figure; imagesc(sub); axis equal; axis tight; colormap turbo;
    % 
    %     figure; imagesc(dmap); axis equal; axis tight; colormap turbo;
    %     figure; imagesc(sub_adjacent); axis equal; axis tight; colormap turbo;
    % 
    % ddd


end

largest_adjacent_diameter = zeros(n_instance,1);
avg_adjacent_diameter = zeros(n_instance,1);
for ki=1:1:n_instance
    labels = connectivity(ki).adjacentinstances;
    if ~isempty(labels)
        for kl=1:1:length(labels)
            id = find(Instancelabels == labels(kl));
            largest_adjacent_diameter(ki) = max([eq_diameter(id) largest_adjacent_diameter(ki)]);
            avg_adjacent_diameter(ki) = avg_adjacent_diameter(ki) + eq_diameter(id);
        end
        avg_adjacent_diameter(ki) = avg_adjacent_diameter(ki)/length(labels);
    else
        largest_adjacent_diameter(ki)=0;
        avg_adjacent_diameter(ki)=0;
    end
end

hold on
for ki=1:1:n_instance
    plot(centroids(ki,2),centroids(ki,1),'Marker','+','MarkerSize',8,'Color','k')
    if ~isempty(connectivity(ki).adjacentinstances)
        x0 = centroids(ki,2);
        y0 = centroids(ki,1);
        for k=1:1:length(connectivity(ki).adjacentinstances)
            x1 = centroids(connectivity(ki).adjacentinstances(k),2);
            y1 = centroids(connectivity(ki).adjacentinstances(k),1);
            plot([x0 x1],[y0 y1],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2);
        end
    end
end

idcrack = 1;
Charact_CrackParticle_orientation(Semantic, idcrack, Instance, centroids, connectivity)

figure;
plot(largest_adjacent_diameter,crack_ratio,'LineStyle','none','Marker','o');
figure;
plot(avg_adjacent_diameter,crack_volume,'LineStyle','none','Marker','o');
dddd





% Connectivity
[Connectivity_matrix, Interface_label] = Function_connectivitymatrix(Instance, p.background, 1, 2);
figure; imagesc(Interface_label); axis equal; axis tight; colormap gray;



figure;
plot(eq_diameter,crack_ratio,'LineStyle','none','Marker','o');
figure;
plot(eq_diameter,crack_volume,'LineStyle','none','Marker','o');

figure;
plot(n_adjacent,crack_ratio,'LineStyle','none','Marker','o');


figure
histogram2(eq_diameter,crack_volume)

nbins = 10;
figure
histogram2(eq_diameter,crack_volume,nbins,'DisplayStyle','tile','ShowEmptyBins','off',...
    'XBinLimits',[0 20],'YBinLimits',[0 1]);
colorbar
xlabel('Red Values')
ylabel('Blue Values')
title('Blue vs. Red Pixel Components')
ax = gca;
ax.CLim = [0 500];



%% CRACK ALIGNMENT and PARTICLE-TO-PARTICLE ALIGNMENT




