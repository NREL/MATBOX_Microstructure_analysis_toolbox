clearvars
close all
clc

% Author: Francois Usseglio-Viretta, NREL
% File used to create geometry tests for journal article:
% F. L. E. Usseglio-Viretta et al., "Quantitative Relationships Between Pore Tortuosity, Pore Topology, and Solid Particle Morphology Using a Novel Discrete Particle Size Algorithm"
% Journal of The Electrochemical Society, 167 100513, 2020
% https://doi.org/10.1149/1945-7111/ab913b


%%
%% TABLE III: Investigating impact of symmetry on oversegmantion
%%

%% #1: IDEAL SPHERE
Domain_size=[200 200 1];
radius = 60;
binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

figure
imagesc(binary_phase)
axis equal

%% #2: SPHERE COMPLEMENT

figure
imagesc(~binary_phase)
axis equal

%% #3: OVERLAPPING SPHERES

Domain_size=[200 300 1]; % Inter res
radius = 60;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(3.1*Domain_size(2)/10),round(Domain_size(3)/2))=1;
binary_phase(round(Domain_size(1)/2),round(6.9*Domain_size(2)/10),round(Domain_size(3)/2))=1;

distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

figure
imagesc(binary_phase)
axis equal

%% #4: OVERLAPPING SHPERES COMPLEMENT

Domain_size=[150 235 1]; % Inter res
radius = 60;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));

binary_phase(1,1,1)=1;
binary_phase(1,round(Domain_size(2)/2),1)=1;
binary_phase(1,Domain_size(2),1)=1;

binary_phase(Domain_size(1),1,1)=1;
binary_phase(Domain_size(1),round(Domain_size(2)/2),1)=1;
binary_phase(Domain_size(1),Domain_size(2),1)=1;

distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

figure
imagesc(binary_phase)
axis equal

%% #5: SPHERES AND CHANNELS

Domain_size=[100 500 1];
radius = 40;
binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(2*Domain_size(2)/10),round(Domain_size(3)/2))=1;
binary_phase(round(Domain_size(1)/2),round(8*Domain_size(2)/10),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;
binary_phase(round(9*Domain_size(1)/20):round(11*Domain_size(1)/20),:,:)=1;

figure
imagesc(binary_phase)
axis equal


%% #6: SPEHRE WITH SPHERICAL INCLUSION

Domain_size=[200 200 1];
radius = 60;
radius_inclusion=10;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

binary_phase2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase2(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase2,'euclidean');
binary_phase2(distance_map<=radius_inclusion)=1;
binary_phase(binary_phase2==1)=0;

figure
imagesc(binary_phase)
axis equal

%% #7: SPHERE WITH ORIENTED INCLUSION

Domain_size=[200 200 1];
radius = 60;
radius_inclusion=20;
length_line = 20;
thickness_line = 8;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

binary_phase(round(Domain_size(1)/2)-radius_inclusion-length_line:round(Domain_size(1)/2)+radius_inclusion+length_line,round(Domain_size(2)/2-thickness_line/2):round(Domain_size(2)/2+thickness_line/2),round(Domain_size(3)/2))=0;
binary_phase(round(Domain_size(1)/2-thickness_line/2):round(Domain_size(1)/2+thickness_line/2),round(Domain_size(2)/2)-radius_inclusion-length_line:round(Domain_size(2)/2)+radius_inclusion+length_line,round(Domain_size(3)/2))=0;

binary_phase2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase2(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase2,'euclidean');
binary_phase2(distance_map<=radius_inclusion)=1;
binary_phase(binary_phase2==1)=0;

figure
imagesc(binary_phase)
axis equal


%% #8 SPHERE WITH OFFSET INCLUSION

Domain_size=[200 200 1];
radius = 60;
radius_inclusion=20;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

binary_phase2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase2(round(Domain_size(1)/2),2*round(Domain_size(2)/3),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase2,'euclidean');
binary_phase2(distance_map<=radius_inclusion)=1;
binary_phase(binary_phase2==1)=0;

figure
imagesc(binary_phase)
axis equal


%% #9: SPEHRE WITH OPEN CRACK

Domain_size=[200 200 1];
radius = 60;
binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

crack_radius = 5;
tmp=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
tmp(round(Domain_size(1)/2),1:round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(tmp,'euclidean');
binary_phase(distance_map<=crack_radius)=0;

figure
imagesc(binary_phase)
axis equal


%% #10: SPEHRE WTIH CLOSED CRACK

Domain_size=[200 200 1];
radius = 60;
binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

crack_radius = 5;
tmp=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
tmp(round(Domain_size(1)/2),100-45:round(Domain_size(2)/2),round(Domain_size(3)/2))=1;
distance_map=bwdist(tmp,'euclidean');
binary_phase(distance_map<=crack_radius)=0;

figure
imagesc(binary_phase)
axis equal


%%
%% FIGURE 10: Imvestigating impact of surface rougness on oversegmentation
%%

%% FIGURE 10a w/o surface roughness

Domain_size=[100 100 1];
width = 10;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% binary_phase(1:round(Domain_size(1)/2),round(Domain_size(2)/2),:)=1;
% binary_phase(round(Domain_size(1)/2),1:round(Domain_size(2)/2),:)=1;

binary_phase(1:round(3*Domain_size(1)/4),round(3*Domain_size(2)/4),:)=1;
binary_phase(round(3*Domain_size(1)/4),1:round(3*Domain_size(2)/4),:)=1;

[D,~] = bwdist(binary_phase,'Chessboard');
binary_phase(D<=width)=1;

figure
imagesc(binary_phase)
axis equal

%% FIGURE 10a w/ surface roughness

Domain_size=[100 100 1];
width = 10;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(1:round(3*Domain_size(1)/4),round(3*Domain_size(2)/4),:)=1;
binary_phase(round(3*Domain_size(1)/4),1:round(3*Domain_size(2)/4),:)=1;

[D,~] = bwdist(binary_phase,'Chessboard');
binary_phase(D<=width)=1;

[D,~] = bwdist(binary_phase,'Chessboard');
binary_phase(D==1)=2;
idx = find(D==1);
for k=1:1:length(idx)
    binary_phase(idx(k))=randi(2)-1;
end

figure
imagesc(binary_phase)
axis equal


%% FIGURE 10 b with smooth surface

Domain_size=[200 300 1];
radius = 60;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(3.1*Domain_size(2)/10),round(Domain_size(3)/2))=1;
binary_phase(round(Domain_size(1)/2),round(6.9*Domain_size(2)/10),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

figure
imagesc(binary_phase)
axis equal


%% FIGURE 10 b with few surface roughness

Domain_size=[200 300 1];
radius = 60;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(3.1*Domain_size(2)/10),round(Domain_size(3)/2))=1;
binary_phase(round(Domain_size(1)/2),round(6.9*Domain_size(2)/10),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

[D,~] = bwdist(binary_phase,'Chessboard');
binary_phase(D==1)=2;
idx = find(D==1);
for k=1:1:length(idx)
    binary_phase(idx(k))=randi(2)-1;
end

figure
imagesc(binary_phase)
axis equal



%% FIGURE 10 c with large surface roughness

Domain_size=[20 30 1];
radius = 6;

binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
binary_phase(round(Domain_size(1)/2),round(3.1*Domain_size(2)/10),round(Domain_size(3)/2))=1;
binary_phase(round(Domain_size(1)/2),round(6.9*Domain_size(2)/10),round(Domain_size(3)/2))=1;
distance_map=bwdist(binary_phase,'euclidean');
binary_phase(distance_map<=radius)=1;

[D,~] = bwdist(binary_phase,'Chessboard');
binary_phase(D==1)=2;
idx = find(D==1);
for k=1:1:length(idx)
    binary_phase(idx(k))=randi(2)-1;
end

binary_phase=imresize(binary_phase,4,'bicubic');
binary_phase(binary_phase>=0.5)=1;
binary_phase(binary_phase~=1)=0;

[L,n]=bwlabel(binary_phase,4);
cluster_=zeros(n,1);
for k=1:1:n
    cluster_(k,1)=sum(sum(sum(L==k)));
end
max_idx = find(cluster_==max(cluster_));
binary_phase(L~=max_idx)=0;

figure
imagesc(binary_phase)
axis equal



