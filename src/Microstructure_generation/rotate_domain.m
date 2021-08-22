function [binary_ellipsoid] = rotate_domain(binary_ellipsoid,angle_x_deg, angle_y_deg, angle_z_deg)
% Angle rotation in radians
angle_x_rad=deg2rad(angle_x_deg);
angle_y_rad=deg2rad(angle_y_deg);
angle_z_rad=deg2rad(angle_z_deg);

% Rotation matrix
% https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
% Transformation matrix that rotates the image around the x-axis
tx = [1 0 0 0
    0 cos(angle_x_rad)  sin(angle_x_rad) 0
    0 -sin(angle_x_rad) cos(angle_x_rad) 0
    0             0              0       1];
% Transformation matrix that rotates the image around the y-axis
ty = [cos(angle_y_rad)  0      -sin(angle_y_rad)   0
    0             1              0     0
    sin(angle_y_rad)    0       cos(angle_y_rad)   0
    0             0              0     1];
% Transformation matrix that rotates the image around the z-axis
tz = [cos(angle_z_rad) sin(angle_z_rad) 0 0
    -sin(angle_z_rad)  cos(angle_z_rad) 0 0
    0 0 1 0
    0 0 0 1];
% Then pass the matrix to the affine3d object constructor.
tx_form = affine3d(tx);
ty_form = affine3d(ty);
tz_form = affine3d(tz);

% Apply the transformation to the image
if angle_x_deg~=0
    binary_ellipsoid = imwarp(binary_ellipsoid,tx_form,'interp','nearest');
end
if angle_y_deg~=0
    binary_ellipsoid = imwarp(binary_ellipsoid,ty_form,'interp','nearest');
end
if angle_z_deg~=0
    binary_ellipsoid = imwarp(binary_ellipsoid,tz_form,'interp','nearest');
end

% Crop
idx = find(binary_ellipsoid~=0);
[I,J,K]=ind2sub(size(binary_ellipsoid),idx);
binary_ellipsoid = binary_ellipsoid(min(I):max(I),min(J):max(J),min(K):max(K));
end

