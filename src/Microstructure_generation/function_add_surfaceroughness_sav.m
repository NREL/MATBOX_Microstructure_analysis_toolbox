
% BW=function_load_tif('C:\Users\fussegli\Desktop\Articles\Specific_surface_area\Test_geometry\Sphere.tif');
% Microstructure_basic_visualization_interface(BW)
% BWsav = BW;

%%

BW = zeros(50,50,50);
BW(25,25,25)=1;
dmap=bwdist(BW);
BW(dmap<=75/4)=1;
BWsav = BW;

ratio_outside_sr = 0.5;
ratio_inside_sr = 0.5;


%%
close all
clc

niter=4;

for iter=1:1:niter

    if ratio_outside_sr>0
        dmap_0to1=bwdist(BW);
        cond1=dmap_0to1< (2-0.001);
        cond2=BW==0;
        id_0t01 = find(cond1.*cond2);
        n=length(id_0t01);
        p = randperm(n,round(ratio_outside_sr*n));
        BW(id_0t01(p))=1;
    end

    if ratio_inside_sr>0
        dmap_1to0=bwdist(~BW);
        cond1=dmap_0to1< (2-0.001);
        cond2=BW==0;
        id_1t00 = find(cond1.*cond2);
        n=length(id_1t00);
        p = randperm(n,round(ratio_outside_sr*n));
        BW(id_1t00(p))=0;
    end

    % Set parameters
    parameters_scaling.scaling_factor = 1/2;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = 0;
    % Scale
    BWscaled = function_scaling(BW,parameters_scaling);

    BW=BWscaled;
end

C=bwlabeln(~BWscaled,6);
[GC,GR] = groupcounts(reshape(C,[numel(C), 1]));
tmp = [GC,GR];
id0 = find(GR==0);
tmp(id0,:)=[];
tmp=sortrows(tmp,1,"descend");
n_cluster = length(tmp(:,1));
if n_cluster>1
    BWscaled = ones(size(BWscaled));
    BWscaled(C==tmp(1,2))=0;
%     for k=2:1:n_cluster
%         BWscaled(C==GR(k+1))=1;
%     end
end
% 
C=bwlabeln(BWscaled,6);
[GC,GR] = groupcounts(reshape(C,[numel(C), 1]));
tmp = [GC,GR];
id0 = find(GR==0);
tmp(id0,:)=[];
tmp=sortrows(tmp,1,"descend");
n_cluster = length(tmp(:,1));
if n_cluster>1
    BWscaled = zeros(size(BWscaled));
    BWscaled(C==tmp(1,2))=1;
    %for k=2:1:n_cluster
    %    BWscaled(C==GR(k+1))=0;
    %end
end

% Uc=unique(C);
% if length(Uc)>2
%     for k=3:1:length(Uc)
%         BWscaled(C==Uc(k))=1;
%     end
% end
% C=bwlabeln(BWscaled,6);
% Uc=unique(C);
% if length(Uc)>2
%     for k=3:1:length(Uc)
%         BWscaled(C==Uc(k))=0;
%     end
% end

Microstructure_comparison_visualization_interface(BWsav,BW,BWscaled);

sz=size(BWscaled);
viewer = viewer3d();
Vol = volshow(BWscaled( round(sz(1)/2):end, round(sz(2)/2):end , round(sz(3)/2):end),Parent=viewer);


%% MESH

opts.radbound = 1;
opts.distbound =  1;
opts.method_surfacemesh = 'lowpass';
opts.iteration_smoothing = 2;
opts.useralpha = 0.5;
opts.keepratio = 1;
opts.maxvol = 10;
[node_, face_, elem_, subdomain_] = function_iso2mesh_from_array(app.Volume+1, opts);


current_phase = 2;
idx = find(subdomain_==current_phase);
vertices_idx_phase = unique(elem_(idx,:));
vertices_phase = node_(vertices_idx_phase,:);
check_face = ismember(face_,vertices_idx_phase);
faces_idx_phase = find( sum(check_face,2)==3 );

%[Lia,Locb] = ismember(round(node_,5), round(vertices_group,5),'rows'); % Precision with 5 digits (tolerance)
[Lia,Locb] = ismembertol(node_, vertices_phase,1e-6,'ByRows',true);
old_id = find(Lia==1);
new_id = Locb(old_id);

meshphase(1).node_ = vertices_phase;
meshphase(1).face_ = my_changem(face_(faces_idx_phase,:), new_id, old_id);
meshphase(1).elem_ = my_changem(elem_(idx,:), new_id, old_id);
meshphase(1).subdomain_ = subdomain_(idx);
meshphase(1).name = 'Surface roughness';

options.show = 'Volume';
options.lighting = 'flat';
options.edges = 'Not visible';
options.coloris = 'z-axis';
options.colormap = 'Default';
options.transparency = 0;
options.grid = 0;
options.axislabel = 0;
options.savefig = 0;
options.savepng = 1;
options.indexstart_zero = 1;
options.folder = 'C:\Users\fussegli\Desktop\T\mesh\';

function_showmesh(meshphase(1).node_, meshphase(1).elem_, meshphase(1).face_, meshphase(1).subdomain_, meshphase(1).name, options);




%%



%function_save_tif(uint8(BWscaled),'C:\Users\fussegli\Desktop\Articles\Specific_surface_area\Test_geometry\Sphere_surfaceroughness_800.tif')


% 
% 
% 
% % Set parameters
% parameters_scaling.scaling_factor = 1/2;
% parameters_scaling.label_or_greylevel = 'Label';
% parameters_scaling.background = 0;
% % Scale
% BWscaled = function_scaling(BW,parameters_scaling);
% 
% 
% parameters_scaling.scaling_factor = 1/3;
% BWscaled2 = function_scaling(BW,parameters_scaling);
% 
% 
% Microstructure_comparison_visualization_interface(BWsav,BW,BWscaled,BWscaled2);
% 




