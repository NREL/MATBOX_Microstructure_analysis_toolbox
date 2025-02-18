function View3D(M,colmap,alphamap,p)

% Size and quality
hFig = uifigure;
scrsz = get(0,'ScreenSize'); % Screen resolution
set(hFig,'position',[scrsz(1) scrsz(2) scrsz(3)*1 scrsz(4)*1]); % Full screen figuredrawnow;
viewer = viewer3d('Parent',hFig);
viewer.RenderingQuality = 'High';

% Box and scale bar
viewer.Box = p.box;
viewer.ScaleBar = p.scalebar;
viewer.ScaleBarUnits = 'Voxels';
viewer.OrientationAxes = p.orientationaxes;

% Background, gradient and lighting
viewer.BackgroundColor = p.backgroundcolor;
viewer.GradientColor = p.gradientcolor;
viewer.BackgroundGradient = p.gradient;
viewer.Lighting="on";

% Volume
Vol = volshow(M,Parent=viewer);

Vol.Colormap = gray;
Vol.Colormap = colmap;
Vol.Alphamap = zeros(1,256);
Vol.RenderingStyle="GradientOpacity";
Vol.OverlayData=M;
Vol.OverlayRenderingStyle="LabelOverlay";
Vol.OverlayThreshold = 0.001;
Vol.LightScatteringQuality = 0.5;

% Phase color and opacity
Vol.OverlayColormap = colmap;
Vol.OverlayAlphamap = alphamap;

% New with 2024b
Vol.OverlayDisplayRange=[min(min(min(M))) max(max(max(M)))];

% Camera position
if p.rotation==0
    p.rotation = 0.01;
end
if p.elevation==0
    p.elevation = 0.01;
end
sz = size(M);
center = sz/2 + 0.5;
viewer.CameraTarget = center;
dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
vec = pi/(p.rotation); 
myPosition = center + ([cos(vec) sin(vec) 1/p.elevation]*dist);
viewer.CameraPosition = myPosition(1,:);

end