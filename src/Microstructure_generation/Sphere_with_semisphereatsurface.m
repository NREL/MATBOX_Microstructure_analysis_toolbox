clear all
close all
clc

%% PARAMETERS
res = 200;
r1=75; % Main sphere radius
n=100; % Number of semisphere

res = 70;
r1=25; % Main sphere radius
n=100; % Number of semisphere

%% ALGORITHM
% Create main sphere
M=zeros(res,res,res);
M(round(res/2),round(res/2),round(res/2))=1;
dmap=bwdist(M);
M(dmap<=r1)=1;

% Create equidistant semispheres
% More complicated than it looks like...
% https://medium.com/@oscarsc/four-ways-to-create-a-mesh-for-a-sphere-d7956b825db4
 
% % Not equidistant
% M2=zeros(200,200,200);
% thetas = linspace(0,pi,n);
% phis = linspace(0,2*pi,2*n-1);
% for ktheta = 1:1:length(thetas)
%     theta=thetas(ktheta);
%     for kphi = 1:1:length(phis)
%         phi=phis(kphi);
%         x = round(r1*cos(phi)*sin(theta) + 100);
%         y = round(r1*sin(phi)*sin(theta) + 100);
%         z = round(r1*cos(theta) + 100);
%         M2(x,y,z)=1;
%     end
% end
% r2=r1*tan(thetas(2)); % Semisphere radius
% dmap=bwdist(M2);
% M(dmap<=r2)=1;

% Equidistant points at the surface of a unit sphere -1/+1
M2=zeros(res,res,res);
[X,Y,Z,N_new] = mySphere(n);
% Convert in indexes
a = ((round(res/2)+r1) - (round(res/2)-r1))/(1+1);
b = (round(res/2)+r1) - a*1;
XX = round(a*X+b);
YY = round(a*Y+b);
ZZ = round(a*Z+b);
for k=1:1:length(XX)
    M2(XX(k),YY(k),ZZ(k))=1;
end

% Distance
all_dist = pdist2([XX' YY' ZZ'],[XX' YY' ZZ']);
all_dist(all_dist==0)=1e9;
r2 = mean(min(all_dist)); % Semisphere radius 

dmap=bwdist(M2);
M(dmap<=r2)=1;

Microstructure_basic_visualization_interface(M)
ddd
Fig_3D=figure;
col_ = bone;
volshow(M,'Parent', Fig_3D,'BackgroundColor','w','Colormap',col_,'Renderer','VolumeRendering');

function_save_tif(uint8(M),'C:\Users\fussegli\Desktop\Articles\Specific_surface_area\Test_geometry\Sphere_semispheres.tif');


function [X,Y,Z,N_new] = mySphere(N)
% Generate Node xyz positions
% Used 2004 paper by Markus Deserno, Max-Planck-Institut:
% "How to generate equidistributed points on the surface of a sphere"
% Enforces constant intervales d_theta ~ d_phi
% Assumes unit radius
% Does not replace MATLAB "sphere" function

% https://www.mathworks.com/matlabcentral/fileexchange/57877-mysphere-n

% Create Sphere 3D Geometry Centered at (x,y,z) = (0,0,0)
%
% N: target number of nodes
% N_new: final number of nodes
% X,Y,Z: column vectors of length N_new containing node coordinates

r_unit = 1;

Area = 4*pi*r_unit^2/N;
Distance = sqrt(Area);
M_theta = round(pi/Distance);
d_theta = pi/M_theta;
d_phi = Area/d_theta;

N_new = 0;
for m = 0:M_theta-1
    
    Theta = pi*(m+0.5)/M_theta;
    M_phi = round(2*pi*sin(Theta)/d_phi); % not exact
    
    for n = 0:M_phi-1        
        Phi = 2*pi*n/M_phi;    
        
        N_new = N_new + 1;
        
        X(N_new) = sin(Theta)*cos(Phi);
        Y(N_new) = sin(Theta)*sin(Phi);
        Z(N_new) = cos(Theta);
        
    end
end

end