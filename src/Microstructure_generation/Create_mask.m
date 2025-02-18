clearvars
close all
clc

domain_size = [100 100 100];

%% Vertical line
% mask = zeros(domain_size);
% mask( round(domain_size(1)/2), round(domain_size(2)/2), :)=1;
% dmapmask = bwdist(mask);
% mask(dmapmask<25)=1;

%% Spiral
% mask = zeros(domain_size);
% x_center = round(domain_size(1)/2);
% y_center = round(domain_size(2)/2);
% z = 1;
% mask( x_center, y_center, z)=1;
% theta = 0;
% r_spiral = 0;
% for z=2:1:domain_size(3)
%     theta = theta+5.0;
%     if theta>=360
%         theta=0;
%     end
%     r_spiral = r_spiral + 0.5;
%     x = round(x_center + r_spiral*cos(deg2rad(theta)));
%     y = round(y_center + r_spiral*sin(deg2rad(theta)));
%     mask( x, y, z)=1;
% end
% dmapmask = bwdist(mask);
% mask(dmapmask<15)=1;

%% Complex
points = {...
 [0 0 0],[1 1 1],[3 0 0],[4 -1 1]...
};
Q = hobbysplines(points,'debug',false,'cycle',false,'tension',0.5,'linestyle','.');
Fig = figure;
plot3(Q(:,1),Q(:,2),Q(:,3),'--');

% Parameters
ntube = 2;
r = 5;
e = 20;

max_angle = 25;
%step = 

% Univariate Gaussian Distribution
% Probability density function of a gaussian distribution
pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

x=[-4:0.1:4];
mu=0;
sigma=1
y=pdf_GD(x,mu,sigma);
figure
plot(x,y)

cdf = pdf2cdf(x,y);
figure
plot(x,cdf)

% mask = zeros(domain_size);
% ktube = 0;
% tube_generated = true;
% while ktube <ntube
%     if tube_generated % reset
%         % select new seed
%         x = randi([e,domain_size(1)-e]);
%         y = randi([e,domain_size(2)-e]);
%         z = randi([e,domain_size(3)-e]);
%         if ktube > 1
%             ddd
%         else
%             forbidden_id = [];
%         end
% 
% 
% 
% 
% 
% 
% 
% 
%     end
% 
% 
% 
% 
%     if tube_generated
%         tube_generated = false;
%         ktube=ktube+1;
%     end
% end

%% Visualization
% mask = uint8(mask);
% Microstructure_basic_visualization_interface(mask);

%% Save
%function_save_tif(mask,'C:\Users\fussegli\Desktop\T\TestGen\Mask.tif');


function c = pdf2cdf(x,pdf)
    n=length(x);
    c=zeros(1,n);
    for k=2:1:n
        c(k) = trapz(x(1:k),pdf(1:k));
    end
end