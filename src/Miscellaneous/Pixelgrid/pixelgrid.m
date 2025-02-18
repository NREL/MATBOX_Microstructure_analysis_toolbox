function hout = pixelgrid(h)
%pixelgrid Superimpose a grid of pixel edges on an image
%   pixelgrid superimposes a grid of pixel edges on the image in the
%   current axes.
%
%   pixelgrid(h_ax) superimposes the grid on the image in the specified
%   axes.
%
%   pixelgrid(h_im) superimposes the grid on the specified image.
%
%   group = pixelgrid(___) returns an hggroup object that contains the two
%   lines that are used to draw the grid. One line is thick and darker, and
%   the other is thin and lighter. Using two contrasting lines in this way
%   guarantees that the grid will be visible regardless of pixel colors.
%
%   EXAMPLES
%
%   Superimpose pixel grid on color image. Zoom in.
%
%       rgb = imread('peppers.png');
%       imshow(rgb)
%       pixelgrid
%       axis([440 455 240 250])
%
%   Change the colors and line widths of the pixel grid.
%
%       rgb = imread('peppers.png');
%       imshow(rgb)
%       h = pixelgrid;
%       axis([440 455 240 250])
%       h.Children(1).Color = [178 223 138]/255;
%       h.Children(1).LineWidth = 2;
%       h.Children(2).Color = [31 120 180]/255;
%       h.Children(2).LineWidth = 4;
%
%   LIMITATIONS
%
%   This function is intended for use when looking at a zoomed-in image
%   region with a relatively small number of rows and columns. If you use
%   this function on a typical, full-size image without zooming in, the
%   image will not be visible under the grid lines.

%   Steve Eddins
%   Copyright 2017-2019 The MathWorks, Inc.

if nargin < 1
    him = findobj(gca,'type','image');
elseif strcmp(h.Type,'axes')
    him = findobj(h,'type','image');
elseif strcmp(h.Type,'image')
    him = h;
else
    error('Invalid graphics object.')
end

if isempty(him)
    error('Image not found.');
end

hax = ancestor(him,'axes');

xdata = him.XData;
ydata = him.YData;
[M,N,~] = size(him.CData);

if M > 1
    pixel_height = diff(ydata) / (M-1);
else
    % Special case. Assume unit height.
    pixel_height = 1;
end

if N > 1
    pixel_width = diff(xdata) / (N-1);
else
    % Special case. Assume unit width.
    pixel_width = 1;
end

y_top = ydata(1) - (pixel_height/2);
y_bottom = ydata(2) + (pixel_height/2);
y = linspace(y_top, y_bottom, M+1);

x_left = xdata(1) - (pixel_width/2);
x_right = xdata(2) + (pixel_width/2);
x = linspace(x_left, x_right, N+1);

% Construct xv1 and yv1 to draw all the vertical line segments. Separate
% the line segments by NaN to avoid drawing diagonal line segments from the
% bottom of one line to the top of the next line over.
xv1 = NaN(1,3*numel(x));
xv1(1:3:end) = x;
xv1(2:3:end) = x;
yv1 = repmat([y(1) y(end) NaN], 1, numel(x));

% Construct xv2 and yv2 to draw all the horizontal line segments.
yv2 = NaN(1,3*numel(y));
yv2(1:3:end) = y;
yv2(2:3:end) = y;
xv2 = repmat([x(1) x(end) NaN], 1, numel(y));

% Put all the vertices together so that they can be drawn with a single
% call to line().
xv = [xv1(:) ; xv2(:)];
yv = [yv1(:) ; yv2(:)];

hh = hggroup(hax);
dark_gray = [0.3 0.3 0.3];
light_gray = [0.8 0.8 0.8];

bottom_line_width = 2;
top_line_width = 1;

% When creating the lines, use AlignVertexCenters to avoid antialias
% effects that would cause some lines in the grid to appear brighter than
% others.
line(...
    'Parent',hh,...
    'XData',xv,...
    'YData',yv,...
    'LineWidth',bottom_line_width,...
    'Color',dark_gray,...
    'LineStyle','-',...
    'AlignVertexCenters','on');
line(...
    'Parent',hh,...
    'XData',xv,...
    'YData',yv,...
    'LineWidth',top_line_width,...
    'Color',light_gray,...
    'LineStyle','-',...
    'AlignVertexCenters','on');

% Only return an output if requested.
if nargout > 0
    hout = hh;
end

