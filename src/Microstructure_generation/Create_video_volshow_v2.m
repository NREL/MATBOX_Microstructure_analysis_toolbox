clearvars
close all
clc

%% IMPORT
% %M = function_load_tif('C:\Users\fussegli\Desktop\T\TestGen_agglomerates\Phaselabel_run_1_bridge_01.tif');

M = function_load_tif('C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\NanoTube_illustration\Agglomerate\Phaselabel_run_1_bridge_02.tif');
M_tubeagglomerate = function_load_tif('C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\NanoTube_illustration\EvenLargerSphere2_OnTopOfTube\Phaselabel_run_1.tif');

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
load('C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\NanoTube_illustration\Volviewer_config.mat');

save_folder = 'C:\Users\fussegli\Desktop\NanoTube_illustration\';
video_filename = 'Video_360_bis';


%% IMPORT
%M = function_load_tif('C:\Users\fussegli\Desktop\T\TestGen_agglomerates\Phaselabel_run_1_bridge_01.tif');

% M = function_load_tif('C:\Users\fussegli\Desktop\Ankit_SSE_proposal_CNT\CNT_particles_additives.tif');
% M_tubeagglomerate = function_load_tif('C:\Users\fussegli\Desktop\Ankit_SSE_proposal_CNT\CNT_agglomerates.tif');
% 
% M =uint8(M);
% M_tubeagglomerate =uint8(M_tubeagglomerate);
% 
% id1_M = find(M==1);
% id2_M = find(M==2);
% M(id1_M)=2;
% M(id2_M)=1;
% id1_M = find(M_tubeagglomerate==1);
% id2_M = find(M_tubeagglomerate==2);
% M_tubeagglomerate(id1_M)=2;
% M_tubeagglomerate(id2_M)=1;
% 
% 
% M_tubeagglomerate(M==2)=2; % Add tube
% 
% sz = size(M);
% 
% M_tube = zeros(sz,"uint8");
% M_tubesphere = zeros(sz,"uint8");
% 
% M_tube(M==2)=2;
% 
% M_tubesphere(M==1)=1;
% M_tubesphere(M==2)=2;
% 
% % Sceneconfig and Object config
% load('C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\NanoTube_illustration\Volviewer_config.mat');
% 
% save_folder = 'C:\Users\fussegli\Desktop\Ankit_SSE_proposal_CNT\';
% video_filename = 'Aligned_CNT';


%% INITIALIZE VISUALIZATION

% Viewer3D fields
viewer = viewer3d();
viewer.BackgroundColor = sceneConfig.BackgroundColor;
viewer.GradientColor = sceneConfig.GradientColor;
viewer.BackgroundGradient = sceneConfig.BackgroundGradient;
viewer.Lighting="on";

% A = M;
A = zeros(size(M));
A(M==2)=2;

% volshow fields
Vol = volshow(A,Parent=viewer);
Vol.Colormap = objectConfig.Colormap;
alphamap = zeros(1,256);
%Vol.Alphamap = objectConfig.Alphamap;
Vol.Alphamap = alphamap;
Vol.RenderingStyle="GradientOpacity";
Vol.OverlayData=A;
Vol.OverlayRenderingStyle="LabelOverlay";
Vol.OverlayThreshold = 0.001;
Vol.OverlayColormap = [0 0 0; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556];
Vol.OverlayAlphamap = [0; 1; 1; 0.15];

dddd


%% VIDEO
scrsz = get(0,'ScreenSize'); % Screen resolution
hFig = viewer.Parent;
set(hFig,'position',[scrsz(1) scrsz(2) scrsz(3)*1 scrsz(4)*1]); % Full screen figuredrawnow;

%%
video_format = 'mpeg-4';
numberOfFrames = 360;
fps = 25;

numberOfFrames = 80;
fps = 5;
complete_turn = 4;

numberOfFrames = 2*complete_turn*360;
fps = 60;

video_handle = VideoWriter([save_folder video_filename],video_format);
set(video_handle,'Quality',100); % Set video quality
set(video_handle,'FrameRate',fps);
open(video_handle)
frame_number = 0;

%vec = linspace(0,2*pi,numberOfFrames)';
vec = linspace(pi/4,complete_turn*2*pi+pi/4,numberOfFrames)'; 

vec = linspace(pi/8,complete_turn*2*pi+pi/4,numberOfFrames)'; 


sz = size(M);
center = sz/2 + 0.5;
viewer.CameraTarget = center;
dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
myPosition = center + ([cos(vec) sin(vec) ones(size(vec))]*dist);

frame_per_turn = numberOfFrames/complete_turn;

scenes = round(linspace(1,numberOfFrames,complete_turn+1));
for idx=1:1:numberOfFrames
    if idx<scenes(2)
        current_M = M_tube;
    elseif idx<scenes(3)
        ida = scenes(2);
        idb = scenes(3);
        idm = round( (ida+idb)/2 );
        if idx >= idm
            current_M = M_tubeagglomerate;
        else
            a = (sz(3)-1) / (idm-ida);
            b = sz(3)-a*idm;
            z = round(a*idx + b);
            z = max([z 1]);
            z = min([z sz(3)]);
            if z==1
                current_M = M_tube;
            elseif z == sz(3)
                current_M = M_tubeagglomerate;
            else
                current_M = zeros(sz,"uint8");
                current_M(:,:,1:z) = M_tubeagglomerate(:,:,1:z);
                current_M(:,:,z+1:end) = M_tube(:,:,z+1:end);
            end
        end

    elseif idx<scenes(4)
        ida = scenes(3);
        idb = scenes(4);
        idm = round( (ida+idb)/2 );
        if idx >= idm
            current_M = M_tubesphere;
        else
            a = (1-sz(3)) / (idm-ida);
            b = 1-a*idm;
            z = round(a*idx + b);
            z = max([z 1]);
            z = min([z sz(3)]);
            if z==1
                current_M = M_tubesphere;
            elseif z == sz(3)
                current_M = M_tubeagglomerate;
            else
                current_M = zeros(sz,"uint8");
                current_M(:,:,1:z) = M_tubeagglomerate(:,:,1:z);
                current_M(:,:,z+1:end) = M_tubesphere(:,:,z+1:end);
            end
        end

    else
        ida = scenes(4);
        idb = scenes(5);
        idm = round( (ida+idb)/2 );
        if idx >= idm
            current_M = M;
        else
            a = (sz(3)-1) / (idm-ida);
            b = sz(3)-a*idm;
            z = round(a*idx + b);
            z = max([z 1]);
            z = min([z sz(3)]);
            if z==1
                current_M = M_tubesphere;
            elseif z == sz(3)
                current_M = M;
            else
                current_M = zeros(sz,"uint8");
                current_M(:,:,1:z) = M(:,:,1:z);
                current_M(:,:,z+1:end) = M_tubesphere(:,:,z+1:end);
            end
        end
    end

    % Add color within a tube (thus invisible) so that we can share colormap. 
    current_M(80,310,2)=1;
    current_M(80,311,2)=3;
    

    Vol.Data = current_M;
    Vol.OverlayData=current_M;



  
  
    viewer.CameraPosition = myPosition(idx,:);
    frame_number = frame_number+1;
    stored_frame(frame_number) = getframe(hFig);
    writeVideo(video_handle,stored_frame(frame_number))
    pause(0.01);

    dddd

end
close(video_handle)


