function [x,y,z] = sph2cart2(az,elev,r)
%SPH2CART Transform spherical to Cartesian coordinates.
%   [X,Y,Z] = SPH2CART(TH,PHI,R) transforms corresponding elements of data
%   stored in spherical coordinates (azimuth TH, elevation PHI, radius R)
%   to Cartesian coordinates X,Y,Z. The arrays TH, PHI, and R must have
%   compatible sizes. In the simplest cases, they can be the same size or
%   one can be a scalar. Two inputs have compatible sizes if, for every
%   dimension, the dimension sizes of the inputs are either the same or one
%   of them is 1. TH and PHI must be in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
%   Class support for inputs TH,PHI,R:
%      float: double, single
%
%   See also CART2SPH, CART2POL, POL2CART.

%   Copyright 1984-2021 The MathWorks, Inc. 

z = r .* sin(elev);
rcoselev = r .* cos(elev);
x = rcoselev .* cos(az);
y = rcoselev .* sin(az);

% rr=rand;
% if rr>0.5
%     rrr=rand;
%     if rrr>0.5
%         tmp = x;
%         x=z;
%         z=tmp;
%     else
%         tmp = y;
%         y=z;
%         z=tmp;
%     end
% end