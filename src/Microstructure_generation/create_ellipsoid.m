function [binary_ellipsoid] = create_ellipsoid(dx,dy,dz)
if ~mod(dx,2) && ~mod(dy,2) && ~mod(dz,2) && dx==dy && dx==dz
    % Sphere, even diameter
    binary_ellipsoid = zeros(dx,dx,dx);
    center1 = dx/2;
    center2 = center1+1;
    binary_ellipsoid(center1,center1,center1) = 1;
    binary_ellipsoid(center2,center1,center1) = 1;
    binary_ellipsoid(center1,center2,center1) = 1;
    binary_ellipsoid(center2,center2,center1) = 1;
    binary_ellipsoid(center1,center1,center2) = 1;
    binary_ellipsoid(center2,center1,center2) = 1;
    binary_ellipsoid(center1,center2,center2) = 1;
    binary_ellipsoid(center2,center2,center2) = 1;    
    dmap_sphere = bwdist(binary_ellipsoid);
    radius = (dx-2)/2;
    binary_ellipsoid(dmap_sphere<=radius) = 1;

else
    % Exact for sphere and ellipsoid with odd diameters. Approximation if one diameter is even
    % Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix
    Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
    [X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
    R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
    binary_ellipsoid = zeros(size(X));
    %binary_ellipsoid(R <= 1 ) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
    idr = find(R<=1);
    if isempty(idr)
        idr = find(R==min(min(min(R))));
    end
    binary_ellipsoid(idr) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
end

% % Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix
% Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
% [X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
% R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
% binary_ellipsoid = zeros(size(X));
% 
% idr = find(R<=1);
% if isempty(idr)
%     idr = find(R==min(min(min(R))));
% end
% binary_ellipsoid(idr) = 1; % Assign 1 for ellipsoid, 0 for complementary volume

end

