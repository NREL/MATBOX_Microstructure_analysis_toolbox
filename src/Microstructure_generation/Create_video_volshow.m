clearvars
close all
clc

%% IMPORT
%M = function_load_tif('C:\Users\fussegli\Desktop\T\TestGen_agglomerates\Phaselabel_run_1_bridge_01.tif');

M = function_load_tif('C:\Users\fussegli\Desktop\T\NanoTube_illustration\Agglomerate\Phaselabel_run_1_bridge_01.tif');
M_tubeagglomerate = function_load_tif('C:\Users\fussegli\Desktop\T\NanoTube_illustration\EvenLargerSphere2_OnTopOfTube\Phaselabel_run_1.tif');

M =uint8(M);
M_tubeagglomerate =uint8(M_tubeagglomerate);

M_tubeagglomerate(M==2)=2; % Add tube

sz = size(M);

M_tube = zeros(sz,"uint8");
M_tubesphere = zeros(sz,"uint8");

M_tube(M==2)=2;

M_tubesphere(M==1)=1;
M_tubesphere(M==2)=2;

% Sceneconfig and Object config
load('C:\Users\fussegli\Desktop\T\NanoTube_illustration\Volviewer_config.mat');

%% INITIALIZE VISUALIZATION

% Viewer3D fields
viewer = viewer3d();
viewer.BackgroundColor = sceneConfig.BackgroundColor;
viewer.GradientColor = sceneConfig.GradientColor;
viewer.BackgroundGradient = sceneConfig.BackgroundGradient;
viewer.Lighting="on";

% volshow fields
Vol = volshow(M_tube,Parent=viewer);
Vol.Colormap = objectConfig.Colormap;
alphamap = zeros(1,256);
%Vol.Alphamap = objectConfig.Alphamap;
Vol.Alphamap = alphamap;
Vol.RenderingStyle="GradientOpacity";
Vol.OverlayData=M_tube;
Vol.OverlayRenderingStyle="LabelOverlay";
Vol.OverlayThreshold = 0.001;
Vol.OverlayColormap = [0 0 0; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556];
Vol.OverlayAlphamap = [0; 1; 1; 0.15];

%% VIDEO
scrsz = get(0,'ScreenSize'); % Screen resolution
hFig = viewer.Parent;
set(hFig,'position',[scrsz(1) scrsz(2) scrsz(3)*5/6 scrsz(4)*5/6]); % Full screen figuredrawnow;

%%
save_folder = 'C:\Users\fussegli\Desktop\T\TestGen_agglomerates\';
video_filename = 'Video2';
video_format = 'mpeg-4';
numberOfFrames = 360;
fps = 25;

numberOfFrames = 80;
fps = 5;
complete_turn = 4;

video_handle = VideoWriter([save_folder video_filename],video_format);
set(video_handle,'Quality',100); % Set video quality
set(video_handle,'FrameRate',fps);
open(video_handle)
frame_number = 0;

%vec = linspace(0,2*pi,numberOfFrames)';
vec = linspace(pi/4,complete_turn*2*pi+pi/4,numberOfFrames)'; 
sz = size(M);
center = sz/2 + 0.5;
viewer.CameraTarget = center;
dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
myPosition = center + ([cos(vec) sin(vec) ones(size(vec))]*dist);

frame_per_turn = numberOfFrames/complete_turn;

scenes = round(linspace(1,numberOfFrames,complete_turn+1));
for idx=1:1:numberOfFrames
    if idx<scenes(2)

    elseif idx<scenes(3)


    elseif idx<scenes(4)

    else

    end


 
 


    viewer.CameraPosition = myPosition(idx,:);
    frame_number = frame_number+1;
    stored_frame(frame_number) = getframe(hFig);
    writeVideo(video_handle,stored_frame(frame_number))
    pause(0.01);
end
close(video_handle)


