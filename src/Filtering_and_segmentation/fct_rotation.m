function [M,newtype,p] = fct_rotation(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

if p.addition
    M = M+1;
end

% https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
%Transformation matrix that rotates the image around the x-axis
angle_rad = deg2rad(p.angle);
if dimension==3
    if strcmp(p.axis,'Axis 1')
        t = [cos(angle_rad)  0      -sin(angle_rad)   0
            0             1              0     0
            sin(angle_rad)    0       cos(angle_rad)   0
            0             0              0     1];
    elseif strcmp(p.axis,'Axis 2')
        t = [1 0 0 0
            0 cos(-angle_rad)  sin(-angle_rad) 0
            0 -sin(-angle_rad) cos(-angle_rad) 0
            0             0              0       1];
    elseif strcmp(p.axis,'Axis 3')
        %Transformation matrix that rotates the image around the z-axis
        t = [cos(angle_rad) sin(angle_rad) 0 0
            -sin(angle_rad)  cos(angle_rad) 0 0
            0 0 1 0
            0 0 0 1];
    end
else
    %Transformation matrix that rotates the image around the z-axis
    t = [cos(angle_rad) sin(angle_rad) 0
        -sin(angle_rad)  cos(angle_rad) 0
        0 0 1];
end
if ~isempty(t)
    % Then pass the matrix to the affine3d object constructor.
    if dimension==3
        t_form = affine3d(t);
    else
        t_form = affine2d(t);
    end

    % Interpolation option
    if strcmp(p.datatype,'label') || strcmp(p.datatype,'Segmented (phase)') || strcmp(p.datatype,'Segmented (instance)')
        M = imwarp(M,t_form,'nearest');
    else
        M = imwarp(M,t_form,'linear'); % Grey level
    end
end

[M] = fct_intconvert(M);

end