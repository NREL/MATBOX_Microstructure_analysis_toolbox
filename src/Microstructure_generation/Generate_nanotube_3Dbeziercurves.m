clearvars
close all
clc

%% THIRD-PARTY LICENSE

% I am using a slighly modified function "hobbysplines" to create smooth 3D bezier curves from:

% Will Robertson (2023). Smooth 3D bezier curves with implicit control points
% https://www.mathworks.com/matlabcentral/fileexchange/42302-smooth-3d-bezier-curves-with-implicit-control-points
% MATLAB Central File Exchange. Retrieved March 29, 2023.
% 
% Copyright (c) 2013, Will Robertson
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% PARAMETERS

domain_size = [400 400 400];
%domain_size = [100 100 100];

domain_size = [200 200 200];


% Tube angular deviation
max_deviation_angle = 50; % deg (0 for straight lines)
sigma_deviation_angle = 5;
distance_between_angle_update = 25;

max_deviation_angle = 15; % deg (0 for straight lines)
sigma_deviation_angle = 5;
distance_between_angle_update = 25;

% Tube dimension
mean_radius = 7;
max_deviation_radius = 3;
sigma_deviation_radius = 1.0;

% Number of tube
n_tube = 60;
n_tube = 30;
n_tube = 10;


% Min length
min_point = round(domain_size(1)/distance_between_angle_update);

% Remove smaller tube
n_remove = 0;

%% TODO
% probe up to x points
%       generate tube, if possible, keep probing for the next x points
%                      if not possible, re-try
% Get all possible control points up to x points, choose among them, of fail remove this possibility and keep trying

% Add mask for existing material: same approach than with parrticles (AT LEAST, NOT MORE, WITHIN)

%% TUBE ANGULAR DEVIATION AND RADIUS
% Univariate Gaussian Distribution
% Probability density function of a gaussian distribution
pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

x_deviation=[-max_deviation_angle:0.1:max_deviation_angle];
pdf_deviation=pdf_GD(x_deviation,0,sigma_deviation_angle);
cdf_deviation = pdf2cdf(x_deviation,pdf_deviation);

x_radius=[mean_radius-max_deviation_radius:0.1:mean_radius+max_deviation_radius];
pdf_radius=pdf_GD(x_radius,mean_radius,sigma_deviation_radius);
cdf_radius = pdf2cdf(x_radius,pdf_radius);

figure
plot(x_deviation,pdf_deviation)
figure
plot(x_deviation,cdf_deviation)

figure
plot(x_radius,pdf_radius)
figure
plot(x_radius,cdf_radius)

%% GENERATION
ext_domain_size = domain_size+50;
ext_domain_size = domain_size+0;

Fig=figure;
hold on;

mask = zeros(ext_domain_size);
mask_id = zeros(ext_domain_size);
edges = mask;
edges(:,:,1)=1;
edges(:,:,end)=1;
edges(:,1,:)=1;
edges(:,end,:)=1;
edges(1,:,:)=1;
edges(end,:,:)=1;
idseeds = find(edges);
all_radius=[];
k_tube = 0;
%new_tube = true;
while k_tube < n_tube
    new_tube = false;
    % Select radius within a distribution
    r=rand;
    radius = interp1(cdf_radius,x_radius,r);
    % Select angle randomly
    azimuth = rand*(2*pi)-pi; % -pi +pi
    elevation = rand*pi-pi/2; % -pi/2 pi/2

    %azimuth = azimuth/24; % -15 +15
    %elevation = elevation/15; % -15 +15

    % Select seed
    idseed = max([1 round(numel(idseeds)*rand)]);
    [x_seed, y_seed, z_seed] = ind2sub(ext_domain_size,idseeds(idseed));

    % Initialize
    points = cell([1,1e6]);
    points(1) = {[x_seed y_seed z_seed]};

    % Get first next point, 0-ref
    [xnew,ynew,znew] = sph2cart(azimuth,elevation,distance_between_angle_update);
    xnew =  round(x_seed+xnew);
    ynew =  round(y_seed+ynew);
    znew =  round(z_seed+znew);
    if xnew>=1 && ynew>=1 && znew>=1 && xnew<=ext_domain_size(1) && ynew<=ext_domain_size(2) && znew<=ext_domain_size(3)
        points(2) = {[xnew ynew znew]};
        npoint = 2;
        keep_advancing = true;
        while keep_advancing
            azimuth = azimuth + deg2rad(interp1(cdf_deviation,x_deviation,rand));
            elevation = elevation + deg2rad(interp1(cdf_deviation,x_deviation,rand));
            [xnew2,ynew2,znew2] = sph2cart(azimuth,elevation,distance_between_angle_update);
            xnew =  round(xnew2+xnew);
            ynew =  round(ynew2+ynew);
            znew =  round(znew2+znew);
            if xnew>=1 && ynew>=1 && znew>=1 && xnew<=ext_domain_size(1) && ynew<=ext_domain_size(2) && znew<=ext_domain_size(3)
                npoint = npoint +1;
                points(npoint) = {[xnew ynew znew]};
            else
                keep_advancing = false;
            end
        end
        if npoint<min_point
            continue % Re-try with a new seed
        end
        points(npoint+1:end)=[];
        Q = hobbysplines(points,'debug',false,'cycle',false,'tension',0.5,'linestyle','-');
        Q = round(Q);
        if min(Q(:,1)) >= 1 && min(Q(:,2)) >= 1 && min(Q(:,3)) >= 1 && max(Q(:,1)) <= ext_domain_size(1) && max(Q(:,2)) <= ext_domain_size(3) && max(Q(:,3)) <= ext_domain_size(3)
            idQ = sub2ind(ext_domain_size,Q(:,1),Q(:,2),Q(:,3));
            idQ(isnan(idQ))=[];
            if k_tube==0
                mask(idQ)=1;
                dmap = bwdist(mask);
                mask(dmap<=radius)=1;
                mask_id(mask==1)=1;
                new_tube = true;
            else
                tmpmask = zeros(ext_domain_size);
                tmpmask(idQ)=1;
                dmap = bwdist(tmpmask);
                tmpmask(dmap<=radius)=1;
                overlapping = mask+tmpmask;
                if max(max(max(overlapping)))==1
                    mask = overlapping;
                    mask_id(tmpmask==1)=k_tube+1;
                    new_tube = true;
                else
                    continue % Re-try with a new seed
                end
            end
        else
            continue % Re-try with a new seed
        end
    else
        continue % Re-try with a new seed
    end

    if new_tube
        k_tube = k_tube +1;
        tubes(k_tube).radius = radius;
        tubes(k_tube).points = points;
        %all_radius(k_tube)=radius;
        k_tube
    end
end

%% SAVE INFO
pars.cycle = false;
pars.tension = 0.5;
save('C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\Test_upscale\TestCNT.mat','tubes','pars','ext_domain_size');



%% REMOVE TINY TUBES
size_tube = zeros(n_tube,2);
for k=1:1:n_tube
    size_tube(k,1) = k;
    size_tube(k,2) = sum(sum(sum(mask_id==k)));
end
size_tube = sortrows(size_tube,2);

mask_tmp = mask;
mask_id_tmp = mask_id;
if n_remove>0
    for k=1:1:n_remove
        mask_tmp( mask_id==size_tube(k,1) )=0;
        mask_id_tmp( mask_id==size_tube(k,1) )=0;
    end
end
mask_tmp = uint8(mask_tmp);
mask_id_tmp = uint8(mask_id_tmp);
vol = volshow(mask_tmp);
vol = volshow(mask_id_tmp);


%% VISUALIZATION
mask = uint8(mask);
mask_id = uint8(mask_id);

Microstructure_basic_visualization_interface(mask);
vol = volshow(mask);     
vol = volshow(mask_id);     

%% EXPORT
%function_save_tif(mask_tmp,'C:\Users\fussegli\Desktop\T\TestGen\Mask_25r.tif');
%function_save_tif(mask_id_tmp,'C:\Users\fussegli\Desktop\T\TestGen\mask_id_25r.tif');

%function_save_tif(mask_tmp,'C:\Users\fussegli\Desktop\T\TestGen\Mask_Large_400_n60.tif');
%function_save_tif(mask_id_tmp,'C:\Users\fussegli\Desktop\T\TestGen\mask_id_Large_400_n60.tif');

% function_save_tif(mask_tmp,'C:\Users\fussegli\Desktop\T\TestGen\CNT_400_n30.tif');
% function_save_tif(mask_id_tmp,'C:\Users\fussegli\Desktop\T\TestGen\CNT_Id_400_n30.tif');

function_save_tif(mask_tmp,'C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\Test_upscale\CNT_n10.tif');
function_save_tif(mask_id_tmp,'C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\Test_upscale\CNT_Id_n10.tif');

%% FUNCTION
function c = pdf2cdf(x,pdf)
    n=length(x);
    c=zeros(1,n);
    for k=2:1:n
        c(k) = trapz(x(1:k),pdf(1:k));
    end
    if c(end)<1
        c(end)=1;
    end
end


