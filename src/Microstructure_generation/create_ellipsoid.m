function [binary_ellipsoid] = create_ellipsoid(dx,dy,dz)
% Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix
Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
[X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
binary_ellipsoid = zeros(size(X));

idr = find(R<=1);
if isempty(idr)
    idr = find(R==min(min(min(R))));
end
binary_ellipsoid(idr) = 1; % Assign 1 for ellipsoid, 0 for complementary volume

end

