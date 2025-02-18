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
Bezier_pars.cycle = false;
Bezier_pars.tension = 0.5;

n_attempt_max = 10;

domain_size = [400 400 400];
%domain_size = [100 100 100];

domain_size = [200 200 200];
domain_size = [50 50 50];


% Tube angular deviation
max_deviation_angle = 50; % deg (0 for straight lines)
sigma_deviation_angle = 5;
distance_between_angle_update = 25;

max_deviation_angle = 15; % deg (0 for straight lines)
sigma_deviation_angle = 5;
distance_between_angle_update = 25;

max_deviation_angle = 50; % deg (0 for straight lines)
sigma_deviation_angle = 25;
distance_between_angle_update = 25;

max_deviation_angle = 50; % deg (0 for straight lines)
sigma_deviation_angle = 25;
distance_between_angle_update = 5;


% Tube dimension
mean_radius = 7;
max_deviation_radius = 3;
sigma_deviation_radius = 1.0;

mean_radius = 2;
max_deviation_radius = 1;
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


x_azimuth=[mean_azimuth-max_deviation_azimuth:0.1:mean_radius+max_deviation_azimuth];
pdf_azimuth=pdf_GD(x_azimuth,mean_radius,sigma_deviation_radius);
cdf_radius = pdf2cdf(x_azimuth,pdf_radius);



% figure
% plot(x_deviation,pdf_deviation)
% figure
% plot(x_deviation,cdf_deviation)
% 
figure
plot(x_radius,pdf_radius)
figure
plot(x_radius,cdf_radius)

%% GENERATION
ext_domain_size = domain_size+50;
ext_domain_size = domain_size+0;

Fig=figure;
hold on;

tube_phase = zeros(ext_domain_size);
tube_id = zeros(ext_domain_size);
tube_skeleton = zeros(ext_domain_size);
mask = zeros(ext_domain_size);
%mask(100,100,100)=1;
%dmap = bwdist(mask);
%mask(dmap<50)=1;

% Locate seeds at domain's edge
edges = tube_phase;
edges(:,:,1)=1;
edges(:,:,end)=1;
edges(:,1,:)=1;
edges(:,end,:)=1;
edges(1,:,:)=1;
edges(end,:,:)=1;
idseeds = find(edges);

k_tube = 0;
while k_tube < n_tube
    k_tube

    % Select radius within a distribution
    r=rand;
    radius = interp1(cdf_radius,x_radius,r);

    % Select initial angle randomly
    azimuth = rand*(2*pi)-pi; % -pi +pi
    elevation = rand*pi-pi/2; % -pi/2 pi/2
    %azimuth = azimuth/24; % -15 +15
    %elevation = elevation/15; % -15 +15

    % Select seed
    idseed = max([1 round(numel(idseeds)*rand)]);
    [x_previous, y_previous, z_previous] = ind2sub(ext_domain_size,idseeds(idseed));

    % Initialize
    npoint = 1;
    points = cell([1,1e6]);
    points(npoint) = {[x_previous y_previous z_previous]};

    n_attempt = 0; % reset
    while true % Exit only through break statment
        tmp_points = points;

        if npoint>1 % New orientation
            azimuth = azimuth + deg2rad(interp1(cdf_deviation,x_deviation,rand));
            elevation = elevation + deg2rad(interp1(cdf_deviation,x_deviation,rand));
        end
        [xstep,ystep,zstep] = sph2cart(azimuth,elevation,distance_between_angle_update);
        xnew =  round(x_previous+xstep);
        ynew =  round(y_previous+ystep);
        znew =  round(z_previous+zstep);

        % Check if out of bounds
        if xnew>=1 && ynew>=1 && znew>=1 && xnew<=ext_domain_size(1) && ynew<=ext_domain_size(2) && znew<=ext_domain_size(3)
            tmp_points(npoint+1) = {[xnew ynew znew]};

            tmp_points(npoint+2:end)=[];
            Q = hobbysplines(tmp_points,'debug',false,'cycle',Bezier_pars.cycle,'tension',Bezier_pars.tension,'linestyle','-');
            Q = round(Q);
            
            if min(Q(:,1)) >= 1 && min(Q(:,2)) >= 1 && min(Q(:,3)) >= 1 && max(Q(:,1)) <= ext_domain_size(1) && max(Q(:,2)) <= ext_domain_size(3) && max(Q(:,3)) <= ext_domain_size(3)
                idQ = sub2ind(ext_domain_size,Q(:,1),Q(:,2),Q(:,3));
                idQ(isnan(idQ))=[];
                % Create tube
                tube_tmp = zeros(ext_domain_size);
                tube_tmp(idQ) = 1;
                dmap = bwdist(tube_tmp);
                tube_tmp(dmap<=radius)=1;
                % Check overlapping 
                overlapping = tube_tmp + tube_phase + mask;

                if max(max(max(overlapping)))==1 % No overlapping with other tube and with mask
                    % Add point 
                    npoint = npoint+1;
                    points(npoint) = tmp_points(npoint);     
                    x_previous = xnew;
                    y_previous = ynew;
                    z_previous = znew;                    
                    n_attempt = 0; % reset
                else
                    n_attempt = n_attempt+1;
                    if n_attempt < n_attempt_max
                        continue % Try with another orientation
                    else
                        break % Try with another seed
                    end
                end

            else
                if npoint==1 % We need another seed
                    break
                else
                    n_attempt = n_attempt+1;
                    if n_attempt < n_attempt_max
                        continue % Try with another orientation
                    else
                        break % Try with another seed
                    end
                end
            end

        else % Out of bounds
            if npoint==1 % We need another seed
                break
            else
                if npoint<min_point % Discard, not long enough. Try with another seed
                    break
                else % Keep tube
                    k_tube = k_tube + 1;
                    tube_id = tube_id + tube_tmp*k_tube;
                    tube_skeleton(idQ) = k_tube;
                    tube_phase = tube_phase + tube_tmp;
                    tubes(k_tube).radius = radius;
                    points(npoint+1:end)=[];
                    tubes(k_tube).points = points;
                    break
                end
            end
        end

    end

end

Microstructure_comparison_visualization_interface(tube_phase,tube_id)

tube_phase = uint8(tube_phase);
tube_id = uint8(tube_id);
tube_skeleton = uint8(tube_skeleton);
vol = volshow(tube_phase);     
vol = volshow(tube_id);     
vol = volshow(tube_skeleton);     

                   %ttt = tube_phase+uint8(mask);
                   % vol = volshow(ttt);   




%% SAVE INFO
save('C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\Test_upscale\TestCNT.mat','tubes','Bezier_pars','ext_domain_size');
dddddddd


%% REMOVE TINY TUBES
size_tube = zeros(n_tube,2);
for k=1:1:n_tube
    size_tube(k,1) = k;
    size_tube(k,2) = sum(sum(sum(tube_id==k)));
end
size_tube = sortrows(size_tube,2);

mask_tmp = tube_phase;
mask_id_tmp = tube_id;
if n_remove>0
    for k=1:1:n_remove
        mask_tmp( tube_id==size_tube(k,1) )=0;
        mask_id_tmp( tube_id==size_tube(k,1) )=0;
    end
end
mask_tmp = uint8(mask_tmp);
mask_id_tmp = uint8(mask_id_tmp);
vol = volshow(mask_tmp);
vol = volshow(mask_id_tmp);


%% VISUALIZATION
tube_phase = uint8(tube_phase);
tube_id = uint8(tube_id);

Microstructure_basic_visualization_interface(tube_phase);
vol = volshow(tube_phase);     
vol = volshow(tube_id);     

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


