function microstructure_visualization_GUI

clear all
%close all
clc


%% FIGURE RENDERING

% Software or hardware (OpenGL) figure rendering
%opengl 'software'
opengl 'hardware'

if ispc
    folder_separation = '\';
else
    folder_separation = '/';
end

%% GLOBAL VARIABLES

global greyvolume segmentedvolume domain_size_grey domain_size_segmented overlay_ Save_folder FileName_grey FileName_segmented

OPTIONS.savefig_informat = {'png'};
OPTIONS.savefig_infig = false;

%%
%% CREATE GUI: OBJECT
%%

%% MAIN FIGURE
main_figure = figure; % Create figure
main_figure.Name= 'Microstructure visualization'; % Set name
main_figure.NumberTitle='off'; % Remove number from title name
main_figure.Color='white'; % Background colour
main_figure.MenuBar='none'; % Remove menubar and toolbar
main_figure.ToolBar='none';
main_figure.Units='normalized'; % Set unit
main_figure.Position=[0.2 0.15 0.6 0.75]; % Set Position


%%
%% PRESELECTED OPTIONS AND VARIABLES
%%

Domain_size = [100 100 100];
volume = [];
is_loaded = false;
voxel_size=1; % nm
field_results={};
current_result_folder=[];
results_filenames={};
property_string=[];
volume_name=[];
result_id_mat=[];
plot_3Dslices=1;
choice_colormap={'gray','jet','parula','hsv','cool','copper','bone'};
sav_volume=volume;
NaN_background = 0;
volume_wo_background =[];
global_max = 0;
global_min = 0;

% % FOR: tab_segmented

is_loaded_greyvolume = 0;
is_loaded_segmentedvolume = 0;
is_correct_voxelsize = 1;
is_domainsize_coherent = 0;
is_savefolder_selected = 0;
figure_initial_filename = 'Segmentation_check';
video_quality = 100;
video_slicepersecond = 10;
back_ground=[];

initial_color_phase=[
    0 114 189;
    217 83 25
    237 177 32];
tmp=randi(255,1000,3);
initial_color_phase=[initial_color_phase;tmp];
color_phase=initial_color_phase/255; % Normalized for Matlab
overlay_ = 0.2;

% Display array normal to direction
direction=3;
Position_slice=[1 1 1];

% Background ?
back_ground_color=0;
remove_background=false;


%%
%% GUI DISPLAY OPTIONS
%%

% These position varialbe are used in various GUI element.
% They are centralized here to assure visual coherence.

% % Font
font_name_GUI ='Times New Roman';

% % Tab description
% Position
description_tab_fromleft = 0.025;
description_tab_frombottom = 0.925;
description_tab_xlenght = 1-2*description_tab_fromleft;
description_tab_ylenght = 0.045;
% Font size
font_size_large_GUI =14;
% Background and text color
background_description_tab = [0 0.5 0];
ForegroundColor_description_tab = [1 1 1];
background_tab = [1 1 1];

error_message_red = [1 0 0];
error_message_orange = [1 0.64 0];

% % Common position
pos_x_start = description_tab_fromleft;
pos_delta_x = 0.1;
pos_y_start = 0.85;
pos_delta_y = 0.075;
thickness_y = 0.05;

% % Normal text
% Font size
font_size_small_GUI =10;
font_size_large_GUI =14;
font_size_large2_GUI=16;
font_size_large3_GUI=18;
font_size_axe = font_size_small_GUI;


%% USER INTERFACE WITH TABBED PANELS

% Create tabgroup
table_group_1 = uitabgroup('Parent', main_figure);
% Set panels location
table_group_1.TabLocation='Left';
% Create tabbed panels
tab_volumeview = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Volume');
tab_segmented = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Segmented volume');

% Preselected tab
table_group_1.SelectedTab = tab_volumeview;


%%
%% TAB VOLUMEVIEW
%%

color_map = gray;
default_colormap = color_map;
shading_3dslice = 'flat';
video_slicepersecond = 25;
is_volume_savefolder_selected = false;
videoquality = 100;
use_adapthisteq=false;
distribution_adapthisteq='uniform';
White_background = 0;
Dark_background =0;
plot_colorbar=0;

colorbar_minmax_global =0;
colorbar_minmax_local=0;
colorbar_minmax_custom=0;
colorbar_min_custom = -9e9;
colorbar_max_custom = 9e9;


%% GUI OBJECTS

% Description and instructions
Text_tab_volumeview = uicontrol('Parent', tab_volumeview, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','3D volume visualization');

% Load file
Text_instructions = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Choose between a .tif file or a volume folder (from the characterization module).',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start pos_y_start 8*pos_delta_x thickness_y]);

Button_import_volume = uicontrol('Parent', tab_volumeview, 'Style', 'pushbutton', 'String', 'Choice 1: import a .tif image',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start pos_y_start-0.6*pos_delta_y 3*pos_delta_x thickness_y],...
    'Callback',{@button_importvolume_Callback});
Text_volume = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','No volume loaded',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+3.25*pos_delta_x pos_y_start-0.8*pos_delta_y 5*pos_delta_x thickness_y]);

Button_import_result = uicontrol('Parent', tab_volumeview, 'Style', 'pushbutton', 'String', 'Choice 2: import a volume result',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start pos_y_start-1.6*pos_delta_y 3*pos_delta_x thickness_y],...
    'Callback',{@button_importresult_Callback});
Text_result = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','No result loaded',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+3.25*pos_delta_x pos_y_start-1.8*pos_delta_y 5*pos_delta_x thickness_y]);


Text_instructions_results = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Select the property to plot, then choose the phase.',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.6 0.8 thickness_y],'Visible','off');
Popup_selectresult = uicontrol('Parent', tab_volumeview,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', 'None','Units','normalized','Position', [pos_x_start 0.57 0.2 thickness_y],'Callback',{@popup_selectresult_Callback},'enable','off','Visible','off');
Popup_selectphase = uicontrol('Parent', tab_volumeview,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', 'None','Units','normalized','Position', [pos_x_start 0.53 0.2 thickness_y],'Callback',{@popup_selectphase_Callback},'enable','off','Visible','off');

Text_instructions_volume = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Enter a voxel size, unit name, and then create the figure.',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.18 0.8 thickness_y]);
Text_volumevoxelsize = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Voxel size in nanometers:',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.135 1.7*pos_delta_x thickness_y]);
Edit_volumevoxelsize = uicontrol('Parent', tab_volumeview,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',voxel_size,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+1.7*pos_delta_x 0.15 0.5*pos_delta_x thickness_y],'Callback', @edit_volumevoxelsize_Callback);

Text_arrayunit = uicontrol('Parent', tab_volumeview,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Unit:',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+2.4*pos_delta_x 0.135 0.05 thickness_y]);
Edit_arrayunit = uicontrol('Parent', tab_volumeview,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+2.4*pos_delta_x+0.05 0.15 0.2 thickness_y]);

checkbox_plot_3Dslices = uicontrol('Parent', tab_volumeview, 'Style', 'checkbox','Units','normalized','Position',[pos_x_start+2.4*pos_delta_x+0.28 0.15 0.4 thickness_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','Plot 3D slices (CPU and Memory expensive)','Value',plot_3Dslices);

Button_create_volumefigure = uicontrol('Parent', tab_volumeview, 'Style', 'pushbutton', 'String', 'Create figure','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start 0.08 3*pos_delta_x thickness_y],...
    'Callback',{@button_createvolume_Callback});

Text_error = uicontrol('Parent', tab_volumeview, 'Style', 'text','Units','normalized','Position',[0.01 0.01 0.98 0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Error message','Visible','off','Enable','off');

%% FIGURE WITH GUI OBJECTS

volume_figure = figure; % Create figure
volume_figure.Name= 'Volume visualization'; % Set name
volume_figure.NumberTitle='off'; % Remove number from title name
volume_figure.Color='white'; % Background colour
volume_figure.MenuBar='none'; % Remove menubar and toolbar
volume_figure.ToolBar='none';
volume_figure.Units='normalized'; % Set unit
scrsz = get(0,'ScreenSize'); % Screen resolution
volume_figure.Position=[0.05 0.05 0.8 0.9];
volume_figure.Visible ='off';

% Figure for volume view
dx  = 0.475;
x0 = 0.025;
dxx = (1-2*dx-2*x0);

dy  = 0.35;
y0 = 0.1;
dyy = (1-2*dy-2*y0);

dsy=0.025;

Text_volumename = uicontrol('Parent', volume_figure, 'Style', 'text','Units','normalized','Position',[0.1 0.925 0.8 0.05],...
    'FontName',font_name_GUI,'FontSize',font_size_large3_GUI,'BackgroundColor','w','String','None');

axes_view_normal_1 = axes('Parent', volume_figure,'FontName',font_name_GUI,'Units','normalized','Position', [x0 y0+dy+dyy dx dy]);
set(axes_view_normal_1,'xtick',[],'ytick',[]); % Remove tick and label
set(axes_view_normal_1,'xticklabel',[],'yticklabel',[]);
box(axes_view_normal_1,'on'); % Box on
axis(axes_view_normal_1,'tight'); % Fit the axes box
axis(axes_view_normal_1,'equal'); % Aspect ratio is 1:1
xlabel(axes_view_normal_1,'Second axis (\mum)');
ylabel(axes_view_normal_1,'Third axis (\mum)');
set(axes_view_normal_1,'FontName',font_name_GUI,'FontSize',font_size_axe); % - Fontname and fontsize
pos1=get(axes_view_normal_1,'Position');
Text_axes_view_normal_1 = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','View normal to first axis',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos1(1) pos1(2)+pos1(4) pos1(3) dsy]);
Slider_axes_view_normal_1 = uicontrol('Parent', volume_figure,'Style', 'slider','Min',1,'Max',Domain_size(1),'Value',1,'SliderStep', [1/(Domain_size(1)-1), 0.1],'Units','normalized','Position', [pos1(1)+0.15*pos1(3) pos1(2)-2.25*dsy 0.35*pos1(3) dsy],'Callback', @slider_axes_view_normal_1_Callback);
str_slider_axes_view_normal_1 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers %s',Position_slice(1), Domain_size(1), Position_slice(1)*voxel_size/1000, Domain_size(1)*voxel_size/1000);
Text_slider_axes_view_normal_1 = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_small_GUI+1,'FontName',font_name_GUI,'String',str_slider_axes_view_normal_1,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [Slider_axes_view_normal_1.Position(1)+Slider_axes_view_normal_1.Position(3)+0.025*pos1(3) pos1(2)-2.5*dsy 0.4*pos1(3) dsy]);

axes_view_normal_2 = axes('Parent', volume_figure,'FontName',font_name_GUI,'Units','normalized','Position', [x0+dx+dxx y0+dy+dyy dx dy]);
set(axes_view_normal_2,'xtick',[],'ytick',[]); % Remove tick and label
set(axes_view_normal_2,'xticklabel',[],'yticklabel',[]);
box(axes_view_normal_2,'on'); % Box on
axis(axes_view_normal_2,'tight'); % Fit the axes box
axis(axes_view_normal_2,'equal'); % Aspect ratio is 1:1
xlabel(axes_view_normal_2,'First axis (\mum)');
ylabel(axes_view_normal_2,'Third axis (\mum)');
set(axes_view_normal_2,'FontName',font_name_GUI,'FontSize',font_size_axe); % - Fontname and fontsize
pos2=get(axes_view_normal_2,'Position');
Text_axes_view_normal_2 = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','View normal to second axis',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos2(1) pos2(2)+pos2(4) pos2(3) dsy]);
Slider_axes_view_normal_2 = uicontrol('Parent', volume_figure,'Style', 'slider','Min',1,'Max',Domain_size(2),'Value',1,'SliderStep', [1/(Domain_size(2)-1), 0.1],'Units','normalized','Position', [pos2(1)+0.15*pos2(3) pos2(2)-2.25*dsy 0.35*pos2(3) dsy],'Callback', @slider_axes_view_normal_2_Callback);
str_slider_axes_view_normal_2 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers %s',Position_slice(2), Domain_size(2), Position_slice(2)*voxel_size/1000, Domain_size(2)*voxel_size/1000);
Text_slider_axes_view_normal_2 = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_small_GUI+1,'FontName',font_name_GUI,'String',str_slider_axes_view_normal_2,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [Slider_axes_view_normal_2.Position(1)+Slider_axes_view_normal_2.Position(3)+0.025*pos2(3) pos2(2)-2.5*dsy 0.4*pos2(3) dsy]);

axes_view_normal_3 = axes('Parent', volume_figure,'FontName',font_name_GUI,'Units','normalized','Position', [x0 y0 dx dy]);
set(axes_view_normal_3,'xtick',[],'ytick',[]); % Remove tick and label
set(axes_view_normal_3,'xticklabel',[],'yticklabel',[]);
box(axes_view_normal_3,'on'); % Box on
axis(axes_view_normal_3,'tight'); % Fit the axes box
axis(axes_view_normal_3,'equal'); % Aspect ratio is 1:1
xlabel(axes_view_normal_3,'First axis (\mum)');
ylabel(axes_view_normal_3,'Second axis (\mum)');
set(axes_view_normal_3,'FontName',font_name_GUI,'FontSize',font_size_axe); % - Fontname and fontsize
pos3=get(axes_view_normal_3,'Position');
Text_axes_view_normal_3 = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','View normal to third axis',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos3(1) pos3(2)+pos3(4) pos3(3) dsy]);
Slider_axes_view_normal_3 = uicontrol('Parent', volume_figure,'Style', 'slider','Min',1,'Max',Domain_size(3),'Value',1,'SliderStep', [1/(Domain_size(3)-1), 0.1],'Units','normalized','Position', [pos3(1)+0.15*pos3(3) pos3(2)-2.25*dsy 0.35*pos3(3) dsy],'Callback', @slider_axes_view_normal_3_Callback);
str_slider_axes_view_normal_3 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers %s',Position_slice(3), Domain_size(3), Position_slice(3)*voxel_size/1000, Domain_size(3)*voxel_size/1000);
Text_slider_axes_view_normal_3 = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_small_GUI+1,'FontName',font_name_GUI,'String',str_slider_axes_view_normal_3,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [Slider_axes_view_normal_3.Position(1)+Slider_axes_view_normal_3.Position(3)+0.025*pos3(3) pos3(2)-2.5*dsy 0.4*pos3(3) dsy]);

axes_3dslices = axes('Parent', volume_figure,'FontName',font_name_GUI,'Units','normalized','Position', [x0+dx+dxx y0 dx dy]);
set(axes_3dslices,'xtick',[],'ytick',[]); % Remove tick and label
set(axes_3dslices,'xticklabel',[],'yticklabel',[]);
box(axes_3dslices,'on'); % Box on
axis(axes_3dslices,'tight'); % Fit the axes box
axis(axes_3dslices,'equal'); % Aspect ratio is 1:1
pos4=get(axes_3dslices,'Position');
xlabel(axes_3dslices,'Second axis (\mum)');
ylabel(axes_3dslices,'First axis (\mum)');
zlabel(axes_3dslices,'Third axis (\mum)');
grid(axes_3dslices,'on'); % Display grid
shading(shading_3dslice)
colormap(color_map)
set(axes_3dslices,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
set(axes_3dslices,'FontName',font_name_GUI,'FontSize',font_size_axe); % - Fontname and fontsize
Text_axes_view_3dslices = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Slices 3D view',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos4(1) pos4(2)+pos4(4) pos4(3) dsy]);

% Colormap
Text_colormap = uicontrol('Parent', volume_figure, 'Style', 'text','Units','normalized','Position',[0.9 0.85 0.1 0.025],...
    'FontName',font_name_GUI,'HorizontalAlignment','center','FontSize',font_size_large_GUI,'BackgroundColor','w','String','Colormap');
Popup_colormap = uicontrol('Parent', volume_figure,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', choice_colormap,'Units','normalized','Position', [0.9 0.825 0.1 0.025],'Callback',{@popup_selectcolormap_Callback});
checkbox_whitebackground = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.8 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','1st color is whiter','Value',White_background,'Callback',{@checkbox_whitebackground_Callback});
checkbox_darkbackground = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.775 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','1st color is darker','Value',Dark_background,'Callback',{@checkbox_darkbackground_Callback});
checkbox_NaNbackground = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.75 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','NaN background','Value',NaN_background,'Callback',{@checkbox_NaNbackground_Callback});

% Colorbar
Text_colorbar = uicontrol('Parent', volume_figure, 'Style', 'text','Units','normalized','Position',[0.9 0.7 0.1 0.025],...
    'FontName',font_name_GUI,'HorizontalAlignment','center','FontSize',font_size_large_GUI,'BackgroundColor','w','String','Colorbar');
checkbox_colorbar = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.675 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','Colorbar','Value',White_background,'Callback',{@checkbox_colorbar_Callback});

checkbox_colorbar_minmax_global = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.65 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','Scale: global','Value',colorbar_minmax_global,'Callback',{@checkbox_colorbar_minmax_global_Callback});
checkbox_colorbar_minmax_local = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.625 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','Scale: local','Value',colorbar_minmax_local,'Callback',{@checkbox_colorbar_minmax_local_Callback});
checkbox_colorbar_minmax_custom = uicontrol('Parent', volume_figure, 'Style', 'checkbox','Units','normalized','Position',[0.9 0.6 0.1 0.025],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor','w',...
    'HorizontalAlignment','left','String','Scale: custom','Value',colorbar_minmax_custom,'Callback',{@checkbox_colorbar_minmax_custom_Callback});
Edit_colorbar_min_custom = uicontrol('Parent', volume_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',colorbar_min_custom,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.9 0.575 0.05 0.025]);
Edit_colorbar_max_custom = uicontrol('Parent', volume_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',colorbar_max_custom,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.95 0.575 0.05 0.025]);

% Save folder
Button_volume_savefolder = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Click to select save folder','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
    'Position', [0.5 0.035 0.25 0.025],'Callback',{@button_volume_savefolder_Callback});
Text_volume_filename = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Filename w/o extension:',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.5 0.005 0.125 0.025]);
Edit_volume_filename = uicontrol('Parent', volume_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','None',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.625 0.005 0.125 0.025]);

% Save figure
Button_volume_save_figure = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Save figure','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
    'Position', [0.76 0.005 0.075 0.055],'enable','off',...
    'Callback',{@button_volume_savefigure_Callback});
Button_volume_save_video = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Save video','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
    'Position', [0.76+0.076 0.035 0.15 0.025],'enable','off',...
    'Callback',{@button_volume_savevideo_Callback});
Text_volume_videoframe = uicontrol('Parent', volume_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Slice per second:',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.76+0.076 0.005 0.075 0.025]);
Edit_volume_videoframe = uicontrol('Parent', volume_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',video_slicepersecond,...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.76+0.076+0.075 0.005 0.075 0.025]);

% Close figure
Button_volume_close_figure = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Close figure','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
    'BackgroundColor',error_message_orange,'Position', [0.9 0.95 0.1 0.05],'Callback',{@button_volume_closefigure_Callback});


%% CALLBACK FUNCTIONS

    function checkbox_colorbar_minmax_global_Callback(~,~)
        colorbar_minmax_global = logical(checkbox_colorbar_minmax_global.Value);
        if colorbar_minmax_global
            colorbar_minmax_local=0; checkbox_colorbar_minmax_local.Value=colorbar_minmax_local;
            colorbar_minmax_custom=0; checkbox_colorbar_minmax_custom.Value=colorbar_minmax_custom;
        end
        
        button_createvolume_Callback
    end

    function checkbox_colorbar_minmax_local_Callback(~,~)
        colorbar_minmax_local = logical(checkbox_colorbar_minmax_local.Value);
        if colorbar_minmax_local
            colorbar_minmax_global=0; checkbox_colorbar_minmax_global.Value=colorbar_minmax_global;
            colorbar_minmax_custom=0; checkbox_colorbar_minmax_custom.Value=colorbar_minmax_custom;
        end
        button_createvolume_Callback
    end

    function checkbox_colorbar_minmax_custom_Callback(~,~)
        colorbar_minmax_custom = logical(checkbox_colorbar_minmax_custom.Value);
        if colorbar_minmax_custom
            colorbar_minmax_local=0; checkbox_colorbar_minmax_local.Value=colorbar_minmax_local;
            colorbar_minmax_global=0; checkbox_colorbar_minmax_global.Value=colorbar_minmax_global;
        end
        button_createvolume_Callback
    end

    function checkbox_colorbar_Callback(~,~)
        plot_colorbar = logical(checkbox_colorbar.Value);
        button_createvolume_Callback
    end

    function checkbox_whitebackground_Callback(~,~)
        White_background = logical(checkbox_whitebackground.Value);
        Dark_background=0;
        checkbox_darkbackground.Value=Dark_background;
        if Dark_background
            color_map(1,:) = 0; % background color is dark
        elseif White_background
            color_map(1,:) = 1; % background color is white
        else
            color_map=default_colormap;
        end
        if NaN_background
            volume=volume_wo_background;
        else
            volume=sav_volume;
        end
        button_createvolume_Callback
    end

    function checkbox_darkbackground_Callback(~,~)
        Dark_background = logical(checkbox_darkbackground.Value);
        White_background=0;
        checkbox_whitebackground.Value=White_background;
        if Dark_background
            color_map(1,:) = 0; % background color is dark
        elseif White_background
            color_map(1,:) = 1; % background color is white
        else
            color_map=default_colormap;
        end
        if NaN_background
            volume=volume_wo_background;
        else
            volume=sav_volume;
        end
        button_createvolume_Callback
    end

    function checkbox_NaNbackground_Callback(~,~)
        NaN_background = logical(checkbox_NaNbackground.Value);
        if Dark_background
            color_map(1,:) = 0; % background color is dark
        elseif White_background
            color_map(1,:) = 1; % background color is white
        else
            color_map=default_colormap;
        end
        if NaN_background
            volume=volume_wo_background;
            
        else
            volume=sav_volume;
        end
        button_createvolume_Callback
    end

    function button_importresult_Callback(~,~)
        str_dialogbox = 'Select the main folder of the volume you want to investigate'; % Set string of the dialog box
        volume_main_folder = uigetdir(pwd,str_dialogbox); % Open dialog box to choose file path
        volume=[]; is_loaded=false; set(Button_create_volumefigure,'enable','off'); % Deactivate related ui
        if volume_main_folder==0
            % User clicked cancel button or closed the dialog box
        else
            current_result_folder = [volume_main_folder folder_separation 'Visualization' folder_separation];
            tmp_ = find(volume_main_folder==folder_separation);
            volume_name = volume_main_folder(tmp_(end)+1:length(volume_main_folder));
            % Find all mat files that start with 'Correlation'
            MyFolderInfo=dir(current_result_folder);
            n_files=length(MyFolderInfo);
            detected_files=0;
            results_filenames={};
            for k_file=1:1:n_files
                [~,filename,ext] = fileparts(MyFolderInfo(k_file).name);
                if strcmp(ext,'.mat')
                    sub_name = filename(1:13);
                    if strcmp(sub_name,'Visualization')
                        detected_files=detected_files+1;
                        results_filenames(detected_files)={[filename '.mat']};
                    end
                end
            end
            if detected_files
                set(Text_error,'String','Error message','BackgroundColor',error_message_red,'Visible','off','Enable','off');
                set(Text_result,'String',['Folder selected: ' volume_main_folder]); % Set name
                set(Text_volume,'String','No volume loaded'); % Set name
                
                number_results = length(results_filenames); % Number of results to summarize
                fields=[];
                for current_matfile=1:1:number_results % Loop over all saved results
                    current_result=char(results_filenames(current_matfile)); % Load result
                    pathname = [current_result_folder current_result];
                    datamat = load(pathname);
                    if isfield(datamat,'results_visualization')
                        data = datamat.results_visualization;
                        fields = [fields; fieldnames(data)];
                    end
                end
                fields = unique(fields);
                field_results={};
                kk=0;
                for k=1:1:length(fields)
                    if ~strcmp(char(fields(k)),'name')
                        kk=kk+1;
                        field_results(kk,1)=fields(k);
                    end
                end
                set(Text_instructions_results,'Visible','on');
                set(Popup_selectresult,'String', field_results,'enable','on','Visible','on');
                
            else
                set(Text_error,'String','Wrong folder! ...\Main_folder_to_select\Visualization\Visualization*.mat','BackgroundColor',error_message_orange,'Visible','on','Enable','on');
            end
        end
    end

    function [] = popup_selectcolormap_Callback(source,~)
        color_map = char(Popup_colormap.String(source.Value));
        color_map = eval(color_map);
        default_colormap=color_map;
        if Dark_background
            color_map(1,:) = 0; % background color is dark
        elseif White_background
            color_map(1,:) = 1; % background color is white
        else
            color_map=default_colormap;
        end
        if NaN_background
            volume=volume_wo_background;
        else
            volume=sav_volume;
        end
        button_createvolume_Callback
    end

    function [] = popup_selectresult_Callback(source,~)
        property_string = char(Popup_selectresult.String(source.Value));
        number_results = length(results_filenames); % Number of results to summarize
        result_id_mat=[];
        for current_matfile=1:1:number_results % Loop over all saved results
            current_result=char(results_filenames(current_matfile)); % Load result
            pathname = [current_result_folder current_result];
            datamat = load(pathname);
            if isfield(datamat,'results_visualization')
                data = datamat.results_visualization;
                fields = fieldnames(data);
                for k=1:1:length(fields)
                    if strcmp( char(fields(k)),property_string)
                        result_id_mat=current_matfile;
                        [~, number_phase] = size(data);
                        phase_results={};
                        for kk=1:1:number_phase
                            phase_results(kk,1)={data(kk).name};
                        end
                        set(Popup_selectphase,'String', phase_results,'enable','on','Visible','on');
                    end
                end
            end
        end
    end

    function [] = popup_selectphase_Callback(source,~)
        phase_string = char(Popup_selectphase.String(source.Value));
        name = [volume_name ', ' property_string ', ' phase_string];
        current_result=char(results_filenames(result_id_mat)); % Load result
        pathname = [current_result_folder current_result];
        datamat = load(pathname);
        data = datamat.results_visualization;
        % Select phase
        [~, number_phase] = size(data);
        for k=1:1:number_phase
            if strcmp(data(k).name, phase_string)
                k_phase=k;
            end
        end
        % Select array
        volume = data(k_phase).(property_string);
        
        sav_volume = volume;
        volume_wo_background = volume;
        volume_wo_background(volume==0)=NaN;
        global_max = nanmax(nanmax(nanmax(volume_wo_background)));
        global_min = nanmin(nanmin(nanmin(volume_wo_background)));
        Edit_colorbar_min_custom.String = num2str(global_min);
        Edit_colorbar_max_custom.String = num2str(global_max);
        
        is_loaded = true;
        Domain_size = size(volume);
        Position_slice=round(Domain_size/2);
        set(Slider_axes_view_normal_1,'Min',1,'Max',Domain_size(1))
        set(Slider_axes_view_normal_2,'Min',1,'Max',Domain_size(2))
        set(Slider_axes_view_normal_3,'Min',1,'Max',Domain_size(3))
        set(Slider_axes_view_normal_1,'SliderStep', [1/(Domain_size(1)-1), 0.1]);
        set(Slider_axes_view_normal_2,'SliderStep', [1/(Domain_size(2)-1), 0.1]);
        set(Slider_axes_view_normal_3,'SliderStep', [1/(Domain_size(3)-1), 0.1]);
        set(Edit_volume_filename,'ForegroundColor','k','String',name);
        set(Text_volumename,'ForegroundColor','k','String',name);
    end

    function [] = GUI_volume_figure(statut_GUI)
        set(Slider_axes_view_normal_1,'Visible',statut_GUI)
        set(Slider_axes_view_normal_2,'Visible',statut_GUI)
        set(Slider_axes_view_normal_3,'Visible',statut_GUI)
        set(Button_volume_savefolder,'Visible',statut_GUI)
        set(Text_volume_filename,'Visible',statut_GUI)
        set(Edit_volume_filename,'Visible',statut_GUI)
        set(Button_volume_save_figure,'Visible',statut_GUI)
        set(Button_volume_save_video,'Visible',statut_GUI)
        set(Text_volume_videoframe,'Visible',statut_GUI)
        set(Edit_volume_videoframe,'Visible',statut_GUI)
        set(Button_volume_close_figure,'Visible',statut_GUI)
        
        set(Text_colormap,'Visible',statut_GUI)
        set(Popup_colormap,'Visible',statut_GUI)
        set(checkbox_whitebackground,'Visible',statut_GUI)
        set(checkbox_darkbackground,'Visible',statut_GUI)
        set(checkbox_NaNbackground,'Visible',statut_GUI)
        set(Text_colorbar,'Visible',statut_GUI)
        set(checkbox_colorbar,'Visible',statut_GUI)
        set(checkbox_colorbar_minmax_global,'Visible',statut_GUI)
        set(checkbox_colorbar_minmax_local,'Visible',statut_GUI)
        set(checkbox_colorbar_minmax_custom,'Visible',statut_GUI)
        set(Edit_colorbar_min_custom,'Visible',statut_GUI)
        set(Edit_colorbar_max_custom,'Visible',statut_GUI)        
        
        if strcmp(statut_GUI,'on')
            set(Text_slider_axes_view_normal_1,'Position',[Slider_axes_view_normal_1.Position(1)+Slider_axes_view_normal_1.Position(3)+0.025*pos1(3) pos1(2)-2.5*dsy 0.4*pos1(3) dsy],'HorizontalAlignment','left')
            set(Text_slider_axes_view_normal_2,'Position',[Slider_axes_view_normal_2.Position(1)+Slider_axes_view_normal_2.Position(3)+0.025*pos2(3) pos2(2)-2.5*dsy 0.4*pos2(3) dsy],'HorizontalAlignment','left')
            set(Text_slider_axes_view_normal_3,'Position',[Slider_axes_view_normal_3.Position(1)+Slider_axes_view_normal_3.Position(3)+0.025*pos3(3) pos3(2)-2.5*dsy 0.4*pos3(3) dsy],'HorizontalAlignment','left')
        else
            set(Text_slider_axes_view_normal_1,'Position',[pos1(1) pos1(2)-2.5*dsy pos1(3) dsy],'HorizontalAlignment','center')
            set(Text_slider_axes_view_normal_2,'Position',[pos2(1) pos2(2)-2.5*dsy pos2(3) dsy],'HorizontalAlignment','center')
            set(Text_slider_axes_view_normal_3,'Position',[pos3(1) pos3(2)-2.5*dsy pos3(3) dsy],'HorizontalAlignment','center')
        end
    end

    function button_volume_savefigure_Callback(~,~)
        GUI_volume_figure('off')
        figure_filename=Edit_volume_filename.String;
        if ispc
            function_savefig(volume_figure, [Save_folder '\'], figure_filename, OPTIONS);
        else
            function_savefig(volume_figure, [Save_folder '/'], figure_filename, OPTIONS);
        end
        GUI_volume_figure('on')
    end

    function button_volume_closefigure_Callback(~,~)
        set(volume_figure,'Visible','off')
    end

    function button_volume_savevideo_Callback(~,~)
        GUI_volume_figure('off')
        figure_filename=Edit_volume_filename.String;
        if ispc
            video_handle = VideoWriter([[Save_folder '\'] figure_filename],'mpeg-4');
        else
            video_handle = VideoWriter([[Save_folder '/'] figure_filename],'mpeg-4');
        end
        set(video_handle,'Quality',videoquality); % Set video quality
        set(video_handle,'FrameRate',str2double(Edit_volume_videoframe.String)); % Set video framerate
        % Open video
        open(video_handle)
        
        total_frame_number=0;
        for direction_video=1:1:3
            Position_slice=round(Domain_size/2);
            for frame_number=1:1:Domain_size(direction_video)
                Position_slice(direction_video)=frame_number;
                % Update figure
                str_slider_axes_view_normal_1 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(1), Domain_size(1), Position_slice(1)*voxel_size/1000, Domain_size(1)*voxel_size/1000);
                str_slider_axes_view_normal_2 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(2), Domain_size(2), Position_slice(2)*voxel_size/1000, Domain_size(2)*voxel_size/1000);
                str_slider_axes_view_normal_3 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(3), Domain_size(3), Position_slice(3)*voxel_size/1000, Domain_size(3)*voxel_size/1000);
                set(Text_slider_axes_view_normal_1,'String',str_slider_axes_view_normal_1);
                set(Text_slider_axes_view_normal_2,'String',str_slider_axes_view_normal_2);
                set(Text_slider_axes_view_normal_3,'String',str_slider_axes_view_normal_3);
                str_xaxis = ['Third axis ' num2str(Domain_size(3)*voxel_size/1000) '\mum'];
                str_yaxis = ['Second axis ' num2str(Domain_size(2)*voxel_size/1000) '\mum'];
                update_axes_view_normal(volume,axes_view_normal_1,1,str_xaxis,str_yaxis)
                str_xaxis = ['Third axis ' num2str(Domain_size(3)*voxel_size/1000) '\mum'];
                str_yaxis = ['First axis ' num2str(Domain_size(1)*voxel_size/1000) '\mum'];
                update_axes_view_normal(volume,axes_view_normal_2,2,str_xaxis,str_yaxis)
                str_xaxis = ['Second axis ' num2str(Domain_size(2)*voxel_size/1000) '\mum'];
                str_yaxis = ['First axis ' num2str(Domain_size(1)*voxel_size/1000) '\mum'];
                update_axes_view_normal(volume,axes_view_normal_3,3,str_xaxis,str_yaxis)
                update_axes_view_3dslices
                % Save frame
                total_frame_number=total_frame_number+1;
                stored_frame(total_frame_number) = getframe(volume_figure);
                writeVideo(video_handle,stored_frame(total_frame_number))
            end
        end
        % Close video
        close(video_handle)
        GUI_volume_figure('on')
    end

    function button_volume_savefolder_Callback(~,~)
        str_dialogbox = 'Select folder where your figure and video will be saved'; % Set string of the dialog box
        Foldername = uigetdir(matlabroot,str_dialogbox); % Open dialog box to choose file path
        if Foldername==0 % User clicked cancel button or closed the dialog box
            is_volume_savefolder_selected=false;
        else
            is_volume_savefolder_selected=true; % Save path
            Save_folder = Foldername;
        end
        if is_volume_savefolder_selected
            set(Button_volume_save_figure,'enable','on');
            set(Button_volume_save_video,'enable','on');
        else
            set(Button_volume_save_figure,'enable','off');
            set(Button_volume_save_video,'enable','off');
        end
    end

    function button_importvolume_Callback(~,~)
        str_dialogbox = 'Select volume to visualize'; % Set string of the dialog box
        [FileName,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)'},str_dialogbox); % Open dialog box to choose file path
        volume=[]; is_loaded=false; set(Button_create_volumefigure,'enable','off'); % Deactivate related ui
        set(Text_instructions_results,'Visible','off');
        set(Popup_selectresult,'String', 'None','enable','off','Visible','off');
        set(Popup_selectphase,'String', 'None','enable','off','Visible','off');
        if FileName==0
            % User clicked cancel button or closed the dialog box
            set(Text_volume,'ForegroundColor','k','String','No volume loaded');
            is_loaded = false;
        else
            full_path_greylevelvolume = [PathName FileName]; % Full path
            [volume, outcome] = function_loadvolume(full_path_greylevelvolume, 'uint16', 'none' ); % Load new volume
            if outcome.success % Success to import
                is_loaded = true;
                Domain_size = size(volume);
                
                sav_volume = volume;
                volume_wo_background = volume;
                volume_wo_background(volume==0)=NaN;
                global_max = nanmax(nanmax(nanmax(volume_wo_background)));
                global_min = nanmin(nanmin(nanmin(volume_wo_background)));
                Edit_colorbar_min_custom.String = num2str(global_min);
                Edit_colorbar_max_custom.String = num2str(global_max);
                
                Position_slice=round(Domain_size/2);
                set(Slider_axes_view_normal_1,'Min',1,'Max',Domain_size(1))
                set(Slider_axes_view_normal_2,'Min',1,'Max',Domain_size(2))
                set(Slider_axes_view_normal_3,'Min',1,'Max',Domain_size(3))
                set(Slider_axes_view_normal_1,'SliderStep', [1/(Domain_size(1)-1), 0.1]);
                set(Slider_axes_view_normal_2,'SliderStep', [1/(Domain_size(2)-1), 0.1]);
                set(Slider_axes_view_normal_3,'SliderStep', [1/(Domain_size(3)-1), 0.1]);                
                 set(Text_volume,'ForegroundColor','k','String',FileName);
                set(Text_result,'String','No results loaded'); % Set name
                [~,name,~] = fileparts(FileName);
                set(Edit_volume_filename,'ForegroundColor','k','String',name);
                set(Text_volumename,'ForegroundColor','k','String',name);
            else
                set(Text_volume,'ForegroundColor','r','String','Failed to load the volume!');
                is_loaded = false;
            end
        end
        if is_loaded
            set(Button_create_volumefigure,'enable','on'); % Activate related ui
        else
            set(Button_create_volumefigure,'enable','off'); % Deactivate related ui
        end
    end

    function edit_volumevoxelsize_Callback(~,~)
        voxel_size=str2double(Edit_volumevoxelsize.String); % Get voxel size
        if isnan(voxel_size) || voxel_size<=0 % Check entered value are correct
            set(Text_volumevoxelsize,'ForegroundColor','r','String','You entered an incorrect value: voxel size must be a real positive number');
            is_correct_voxelsize=0;
        else
            set(Text_volumevoxelsize,'ForegroundColor','k','String','Voxel size in nanometers');
            is_correct_voxelsize=1;
        end
        if is_loaded==1 && is_correct_voxelsize==1
            set(Button_create_volumefigure,'enable','on');% Activate related ui
        else
            set(Button_create_volumefigure,'enable','off');% Deactivate related ui
        end
    end

    function button_createvolume_Callback(~,~)
        str_slider_axes_view_normal_1 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(1), Domain_size(1), Position_slice(1)*voxel_size/1000, Domain_size(1)*voxel_size/1000);
        str_slider_axes_view_normal_2 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(2), Domain_size(2), Position_slice(2)*voxel_size/1000, Domain_size(2)*voxel_size/1000);
        str_slider_axes_view_normal_3 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(3), Domain_size(3), Position_slice(3)*voxel_size/1000, Domain_size(3)*voxel_size/1000);
        set(Text_slider_axes_view_normal_1,'String',str_slider_axes_view_normal_1);
        set(Text_slider_axes_view_normal_2,'String',str_slider_axes_view_normal_2);
        set(Text_slider_axes_view_normal_3,'String',str_slider_axes_view_normal_3);
        set(Slider_axes_view_normal_1,'Value',Position_slice(1));
        set(Slider_axes_view_normal_2,'Value',Position_slice(2));
        set(Slider_axes_view_normal_3,'Value',Position_slice(3));
        update_axes_view_normal(volume,axes_view_normal_1,1,['Third axis ' num2str(Domain_size(3)*voxel_size/1000) '\mum'],['Second axis ' num2str(Domain_size(2)*voxel_size/1000) '\mum'])
        update_axes_view_normal(volume,axes_view_normal_2,2,['Third axis ' num2str(Domain_size(3)*voxel_size/1000) '\mum'],['First axis ' num2str(Domain_size(1)*voxel_size/1000) '\mum'])
        update_axes_view_normal(volume,axes_view_normal_3,3,['First axis ' num2str(Domain_size(1)*voxel_size/1000) '\mum'],['Second axis ' num2str(Domain_size(2)*voxel_size/1000) '\mum'])
        update_axes_view_3dslices
        set(volume_figure,'Visible','on');
    end

    function slider_axes_view_normal_1_Callback(source,~)
        pos_=round(source.Value); % Get position value
        Position_slice(1)=pos_;% Update position array
        % Update text
        str_slider_axes_view_normal_1 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(1), Domain_size(1), Position_slice(1)*voxel_size/1000, Domain_size(1)*voxel_size/1000);
        set(Text_slider_axes_view_normal_1,'String',str_slider_axes_view_normal_1);
        set(Slider_axes_view_normal_1,'SliderStep', [1/(Domain_size(1)-1), 0.1]);
        % Update figure
        str_xaxis = ['Third axis ' num2str(Domain_size(3)*voxel_size/1000) '\mum'];
        str_yaxis = ['Second axis ' num2str(Domain_size(2)*voxel_size/1000) '\mum'];
        update_axes_view_normal(volume,axes_view_normal_1,1,str_xaxis,str_yaxis)
        update_axes_view_3dslices
    end

    function slider_axes_view_normal_2_Callback(source,~)
        pos_=round(source.Value); % Get position value
        Position_slice(2)=pos_;% Update position array
        % Update text
        str_slider_axes_view_normal_2 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(2), Domain_size(2), Position_slice(2)*voxel_size/1000, Domain_size(2)*voxel_size/1000);
        set(Text_slider_axes_view_normal_2,'String',str_slider_axes_view_normal_2);
        set(Slider_axes_view_normal_2,'SliderStep', [1/(Domain_size(2)-1), 0.1]);
        % Update figure
        str_xaxis = ['Third axis ' num2str(Domain_size(3)*voxel_size/1000) '\mum'];
        str_yaxis = ['First axis ' num2str(Domain_size(1)*voxel_size/1000) '\mum'];
        update_axes_view_normal(volume,axes_view_normal_2,2,str_xaxis,str_yaxis)
        update_axes_view_3dslices
    end

    function slider_axes_view_normal_3_Callback(source,~)
        pos_=round(source.Value); % Get position value
        Position_slice(3)=pos_;% Update position array
        % Update text
        str_slider_axes_view_normal_3 = sprintf('Slice: %i/%i, %1.1f/%1.1f micrometers',Position_slice(3), Domain_size(3), Position_slice(3)*voxel_size/1000, Domain_size(3)*voxel_size/1000);
        set(Text_slider_axes_view_normal_3,'String',str_slider_axes_view_normal_3);
        set(Slider_axes_view_normal_3,'SliderStep', [1/(Domain_size(3)-1), 0.1]);
        % Update figure
        str_xaxis = ['Second axis ' num2str(Domain_size(2)*voxel_size/1000) '\mum'];
        str_yaxis = ['First axis ' num2str(Domain_size(1)*voxel_size/1000) '\mum'];
        update_axes_view_normal(volume,axes_view_normal_3,3,str_xaxis,str_yaxis)
        update_axes_view_3dslices
    end

    function [] = update_axes_view_normal(volume,axe_,direction,str_xaxis,str_yaxis)
        if direction==1
            slice_ = squeeze(volume(Position_slice(1),:,:));
        elseif direction==2
            slice_ = squeeze(volume(:,Position_slice(2),:));
        elseif direction==3
            slice_ = volume(:,:,Position_slice(3));
        end
        if use_adapthisteq
            imagesc( adapthisteq(slice_(:,:,1),'Distribution',distribution_adapthisteq) ,'CDataMapping','scaled','Parent', axe_);
        else
            h=pcolor(slice_(:,:,1),'Parent', axe_); % Pcolor: NaN value are transparent
            set(h,'EdgeColor','none');
        end
        colormap(axe_,color_map)
        if plot_colorbar
            h=colorbar(axe_);
            ylabel(h, Edit_arrayunit.String);
            set(h,'FontName',font_name_GUI,'FontSize',font_size_axe);
        end
        if colorbar_minmax_global
            caxis(axe_,[global_min global_max])
        elseif colorbar_minmax_local
            local_min = nanmin(nanmin(nanmin(slice_)));
            local_max = nanmax(nanmax(nanmax(slice_)));
            caxis(axe_,[local_min local_max])
        elseif colorbar_minmax_custom
            custom_min = str2num(Edit_colorbar_min_custom.String);
            custom_max = str2num(Edit_colorbar_max_custom.String);
            caxis(axe_,[custom_min custom_max])
        else
            caxis(axe_,'auto');
        end
        xlabel(axe_,str_xaxis);
        ylabel(axe_,str_yaxis);
        set(axe_,'FontName',font_name_GUI,'FontSize',font_size_axe); % - Fontname and fontsize
        axe_position =axe_.Position; % Get axis position
        set(axe_ ,'Position', axe_position);
        set(axe_,'xtick',[],'ytick',[]); % Remove tick and label
        set(axe_,'xticklabel',[],'yticklabel',[]);
        axis(axe_,'tight'); % Fit the axes box
        axis(axe_,'equal'); % Aspect ratio is 1:1
    end

    function [] = update_axes_view_3dslices
        cla(axes_3dslices)
        hold(axes_3dslices,'on');
        view(axes_3dslices,3)
        %box(axes_3dslices,'off');
        if checkbox_plot_3Dslices.Value
            slice(double(volume),Position_slice(2),Position_slice(1),Position_slice(3),'Parent', axes_3dslices);
            if min(min(default_colormap == gray)) || min(min(default_colormap == bone))
                c = [get(0, 'DefaultAxesColorOrder')];
                line_color1=c(1,:);
                line_color2=c(2,:);
                line_color3=c(3,:);
            else
                if White_background || NaN_background
                    line_color1='k';
                    line_color2='k';
                    line_color3='k';
                else
                    line_color1='w';
                    line_color2='w';
                    line_color3='w';
                end
            end
        else
            X = [0 Domain_size(2) Domain_size(2) 0];
            Y = [Position_slice(1) Position_slice(1) Position_slice(1) Position_slice(1)];
            Z = [0 0 Domain_size(3) Domain_size(3)];
            fill3(X,Y,Z,[0.7 0.7 0.7 0.7],'FaceAlpha',.5,'Parent', axes_3dslices)
            
            Y = [0 Domain_size(1) Domain_size(1) 0];
            X = [Position_slice(2) Position_slice(2) Position_slice(2) Position_slice(2)];
            Z = [0 0 Domain_size(3) Domain_size(3)];
            fill3(X,Y,Z,[0.6 0.6 0.6 0.6],'FaceAlpha',.5,'Parent', axes_3dslices)
            
            X = [0 Domain_size(2) Domain_size(2) 0];
            Y = [0 0 Domain_size(1) Domain_size(1)];
            Z = [Position_slice(3) Position_slice(3) Position_slice(3) Position_slice(3)];
            fill3(X,Y,Z,[0.5 0.5 0.5 0.5],'FaceAlpha',0.5,'Parent', axes_3dslices)
            
            line_color1='k';
            line_color2='k';
            line_color3='k';
        end
        xlabel(axes_3dslices,'Second axis (\mum)');
        ylabel(axes_3dslices,'First axis (\mum)');
        zlabel(axes_3dslices,'Third axis (\mum)');
        X=[Position_slice(2) Position_slice(2)];
        Y=[0 Domain_size(1)];
        Z=[Position_slice(3) Position_slice(3)];
        plot3(X,Y,Z,'--','Color',line_color1,'LineWidth',2,'Parent', axes_3dslices)
        X=[0 Domain_size(2)];
        Y=[Position_slice(1) Position_slice(1)];
        Z=[Position_slice(3) Position_slice(3)];
        plot3(X,Y,Z,'--','Color',line_color2','LineWidth',2,'Parent', axes_3dslices)
        X=[Position_slice(2) Position_slice(2)];
        Y=[Position_slice(1) Position_slice(1)];
        Z=[0 Domain_size(3)];
        plot3(X,Y,Z,'--','Color',line_color3,'LineWidth',2,'Parent', axes_3dslices)
        grid(axes_3dslices,'on'); % Display grid
        set(axes_3dslices,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
        set(axes_3dslices,'FontName',font_name_GUI,'FontSize',font_size_axe); % - Fontname and fontsize
        axis(axes_3dslices,'tight'); % Fit the axes box
        axis(axes_3dslices,'equal'); % Aspect ratio is 1:1
        if checkbox_plot_3Dslices.Value
            colormap(color_map)
            if colorbar_minmax_global
                caxis(axes_3dslices,[global_min global_max])
            elseif colorbar_minmax_local
                caxis(axes_3dslices,[global_min global_max]) % No local min-max (3 slices)
            elseif colorbar_minmax_custom
                custom_min = str2num(Edit_colorbar_min_custom.String);
                custom_max = str2num(Edit_colorbar_max_custom.String);
                caxis(axes_3dslices,[custom_min custom_max])
            else
                caxis(axes_3dslices,'auto');
            end
            
        else
            colormap(axes_3dslices,gray)
        end
        %light(axes_3dslices,'Position',[-Domain_size(1) -Domain_size(2) 2*Domain_size(3)],'Style','local')
        %lighting(axes_3dslices,'flat');
        shading(axes_3dslices,shading_3dslice)
        set(axes_3dslices,'XTickMode','auto')
        set(axes_3dslices,'YTickMode','auto')
        set(axes_3dslices,'ZTickMode','auto')        
        x_value = get(axes_3dslices,'XTick');
        set(axes_3dslices,'XtickLabel',round(x_value*voxel_size/1000,1));
        y_value = get(axes_3dslices,'YTick');
        set(axes_3dslices,'YtickLabel',round(y_value*voxel_size/1000,1));
        z_value = get(axes_3dslices,'ZTick');
        set(axes_3dslices,'ZtickLabel',round(z_value*voxel_size/1000,1));
        hold(axes_3dslices,'off');
    end


%%
%% TAB SEGMENTED
%%

%% GUI

% Description and instructions
Text_tab_segmented = uicontrol('Parent', tab_segmented, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Select grey level and segmented volume files to compare.');

% Load file
Text_instructions = uicontrol('Parent', tab_segmented,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Select both files, enter a voxel size and then create figure.',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start pos_y_start 4*pos_delta_x thickness_y]);

Button_import_greylevelimage = uicontrol('Parent', tab_segmented, 'Style', 'pushbutton', 'String', 'Click to import grey-level image',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start pos_y_start-0.6*pos_delta_y 3*pos_delta_x thickness_y],'UserData',0,...
    'Callback',{@import_greylevelimage_Callback});
Text_greylevelvolume = uicontrol('Parent', tab_segmented,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','No grey-level volume loaded',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+3.25*pos_delta_x pos_y_start-0.8*pos_delta_y 5*pos_delta_x thickness_y]);

Button_import_segmentedimage = uicontrol('Parent', tab_segmented, 'Style', 'pushbutton', 'String', 'Click to import segmented image',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start pos_y_start-1.6*pos_delta_y 3*pos_delta_x thickness_y],'UserData',0,...
    'Callback',{@import_segmentedimage_Callback});
Text_segmentedvolume = uicontrol('Parent', tab_segmented,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','No segmented volume loaded',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+3.25*pos_delta_x pos_y_start-1.8*pos_delta_y 5*pos_delta_x thickness_y]);

Edit_voxelsize = uicontrol('Parent', tab_segmented,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',voxel_size,...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start pos_y_start-2.6*pos_delta_y 3*pos_delta_x thickness_y],'Callback', @edit_voxelsize_Callback);
Text_voxel_size = uicontrol('Parent', tab_segmented,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Voxel size in nanometers',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [pos_x_start+3.25*pos_delta_x pos_y_start-2.8*pos_delta_y 5*pos_delta_x thickness_y]);

Button_create_figure = uicontrol('Parent', tab_segmented, 'Style', 'pushbutton', 'String', 'Create figure','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start pos_y_start-3.5*pos_delta_y 3*pos_delta_x thickness_y],'UserData',0,...
    'Callback',{@create_segmentationfigure_Callback});


%% CALLBACKS

    function import_greylevelimage_Callback(~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select grey-level volume';
        % Open dialog box to choose file path
        [FileName_grey,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)'},str_dialogbox);
        if FileName_grey==0
            % User clicked cancel button or closed the dialog box
            set(Text_greylevelvolume,'ForegroundColor','k','String','No grey-level volume loaded');
            is_loaded_greyvolume = 0;
        else
            full_path_greylevelvolume = [PathName FileName_grey]; % Full path
            [greyvolume, outcome] = function_loadvolume(full_path_greylevelvolume, 'uint16', 'none' ); % Load new volume
            if outcome.success % Success to import
                % Update GUI text field and user data, get remove button enable, update number of volume
                is_loaded_greyvolume = 1;
                domain_size_grey = size(greyvolume);
                if is_loaded_greyvolume == 1 && is_loaded_segmentedvolume==1
                    if domain_size_grey == domain_size_segmented
                        is_domainsize_coherent=1;
                        set(Text_greylevelvolume,'ForegroundColor','k','String',FileName_grey);
                    else
                        is_domainsize_coherent=0;
                        set(Text_greylevelvolume,'ForegroundColor','r','String','Volumes have different fields of view!');
                    end
                else
                    set(Text_greylevelvolume,'ForegroundColor','k','String',FileName_grey);
                end
            else
                set(Text_greylevelvolume,'ForegroundColor','r','String','Failed to load the volume!');
                is_loaded_greyvolume = 0;
            end
        end
        if is_loaded_greyvolume==1 && is_loaded_segmentedvolume==1 && is_correct_voxelsize==1 && is_domainsize_coherent==1
            % Activate related ui
            set(Button_create_figure,'enable','on');
        else
            % Deactivate related ui
            set(Button_create_figure,'enable','off');
        end
    end

    function import_segmentedimage_Callback(~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select segmented volume';
        % Open dialog box to choose file path
        [FileName_segmented,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)'},str_dialogbox);
        if FileName_segmented==0
            % User clicked cancel button or closed the dialog box
            set(Text_segmentedvolume,'ForegroundColor','k','String','No segmented volume loaded');
            is_loaded_segmentedvolume = 0;
        else
            full_path_segmentedvolume = [PathName FileName_segmented]; % Full path
            [segmentedvolume, outcome] = function_loadvolume(full_path_segmentedvolume, 'uint16', 'none' ); % Load new volume
            if outcome.success % Success to import
                % Update GUI text field and user data, get remove button enable, update number of volume
                is_loaded_segmentedvolume = 1;
                domain_size_segmented = size(segmentedvolume);
                if remove_background==true
                    for k=1:1:domain_size_segmented(3)
                        back_ground(k).idx = find(segmentedvolume(:,:,k)==back_ground_color);
                    end
                end
                if is_loaded_greyvolume == 1 && is_loaded_segmentedvolume==1
                    if domain_size_grey == domain_size_segmented
                        is_domainsize_coherent=1;
                        set(Text_segmentedvolume,'ForegroundColor','k','String',FileName_segmented);
                    else
                        is_domainsize_coherent=0;
                        set(Text_segmentedvolume,'ForegroundColor','r','String','Volumes have different fields of view!');
                    end
                else
                    set(Text_segmentedvolume,'ForegroundColor','k','String',FileName_segmented);
                end
            else
                set(Text_segmentedvolume,'ForegroundColor','r','String','Failed to load the volume!');
                is_loaded_segmentedvolume = 0;
            end
        end
        if is_loaded_greyvolume==1 && is_loaded_segmentedvolume==1 && is_correct_voxelsize==1 && is_domainsize_coherent==1
            % Activate related ui
            set(Button_create_figure,'enable','on');
        else
            % Deactivate related ui
            set(Button_create_figure,'enable','off');
        end
    end

    function edit_voxelsize_Callback(~,~)
        % Get voxel size
        voxel_size=str2double(Edit_voxelsize.String);
        % Check entered value are correct
        if isnan(voxel_size) || voxel_size<=0
            % Wrong entry
            set(Text_voxel_size,'ForegroundColor','r','String','You entered an incorrect value: voxel size must be a real positive number');
            is_correct_voxelsize=0;
        else
            set(Text_voxel_size,'ForegroundColor','k','String','Voxel size in nanometers');
            is_correct_voxelsize=1;
        end
        if is_loaded_greyvolume==1 && is_loaded_segmentedvolume==1 && is_correct_voxelsize==1
            % Activate related ui
            set(Button_create_figure,'enable','on');
        else
            % Deactivate related ui
            set(Button_create_figure,'enable','off');
        end
    end

    function create_segmentationfigure_Callback(~,~)
        % Rest initial values
        direction=3;
        Position_slice=[1 1 1];
        color_phase=initial_color_phase/255; % Normalized for Matlab
        overlay_ = 0.2;
        
        segmented_grey_figure = figure; % Create figure
        segmented_grey_figure.Name= 'Visual comparison between grey-level and segmented image'; % Set name
        segmented_grey_figure.NumberTitle='off'; % Remove number from title name
        segmented_grey_figure.Color='white'; % Background colour
        % Remove menubar and toolbar
        segmented_grey_figure.MenuBar='none';
        segmented_grey_figure.ToolBar='none';
        % Set unit
        segmented_grey_figure.Units='normalized';
        % Set Position
        scrsz = get(0,'ScreenSize'); % Screen resolution
        aspect_ratio = scrsz(3)/scrsz(4);
        if aspect_ratio==1.6
            segmented_grey_figure.Position=[0.05 0.05 0.9/1.4 0.9];
        else
            % Default
            segmented_grey_figure.Position=[0.05 0.05 0.7 0.9];
        end
        
        
        % Max grey level
        max_grey = double(max(max(max(greyvolume))));
        % Domain size
        Domain_size=size(greyvolume);
        % Phase information
        Phase_code = unique(segmentedvolume);
        number_phase  = length(Phase_code);
        % Get volume fraction
        total_number_voxel=prod(Domain_size);
        volumefraction=zeros(number_phase,1);
        for n=1:1:number_phase
            volumefraction(n,1)=sum(sum(sum(segmentedvolume==Phase_code(n))))/total_number_voxel;
        end
        % Set default name
        for n=1:1:number_phase
            Phase(n).name=['Phase ' num2str(Phase_code(n))];
        end
        % Set default colors
        for current_phase=1:1:number_phase
            RGB_red.index(current_phase)=color_phase(current_phase,1);
            RGB_green.index(current_phase)=color_phase(current_phase,2);
            RGB_blue.index(current_phase)=color_phase(current_phase,3);
            RGB_phase.index(current_phase).rgb = [color_phase(current_phase,1) color_phase(current_phase,2) color_phase(current_phase,3)];
        end
        
        % Figure for volume view
        axes_greylevel = axes('Parent', segmented_grey_figure,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.5+0.025 0.5-0.05 0.5-0.05]);
        set(axes_greylevel,'xtick',[],'ytick',[]); % Remove tick and label
        set(axes_greylevel,'xticklabel',[],'yticklabel',[]);
        box(axes_greylevel,'on'); % Box on
        axis(axes_greylevel,'tight'); % Fit the axes box
        axis(axes_greylevel,'equal'); % Aspect ratio is 1:1
        axes_greylevel.YDir='normal';

        axes_segmented = axes('Parent', segmented_grey_figure,'FontName',font_name_GUI,'Units','normalized','Position', [0.525 0.5+0.025 0.5-0.05 0.5-0.05]);
        set(axes_segmented,'xtick',[],'ytick',[]); % Remove tick and label
        set(axes_segmented,'xticklabel',[],'yticklabel',[]);
        box(axes_segmented,'on'); % Box on
        axis(axes_segmented,'tight'); % Fit the axes box
        axis(axes_segmented,'equal'); % Aspect ratio is 1:1
        axes_segmented.YDir='normal';
        
        axes_overlay_ = axes('Parent', segmented_grey_figure,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.025 0.5-0.05 0.5-0.05]);
        set(axes_overlay_,'xtick',[],'ytick',[]); % Remove tick and label
        set(axes_overlay_,'xticklabel',[],'yticklabel',[]);
        box(axes_overlay_,'on'); % Box on
        axis(axes_overlay_,'tight'); % Fit the axes box
        axis(axes_overlay_,'equal'); % Aspect ratio is 1:1
        axes_overlay_.YDir='normal';
       
        
        % Text
        Text_directionslider = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Select direction normal to the 2D view',...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.5 0.42 0.5 0.05]);
        % Popup
        Popup_slider = uicontrol('Parent', segmented_grey_figure,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
            'String', {'Direction 1','Direction 2','Direction 3'},'Value',3,'Units','normalized','Position', [0.65 0.4 0.2 0.03],'Callback', @popup_axes2dview_Callback);
        % Slider
        Slider_axes_2dview = uicontrol('Parent', segmented_grey_figure,'Style', 'slider','Min',1,'Max',Domain_size(direction),'Value',1,'Units','normalized','Position', [0.65 0.37 0.2 0.03],'Callback', @slider_axes2dview_Callback);
        str_slider = sprintf('Slice: %i/%i, %1.3f/%1.3f micrometers %s',Position_slice(direction), Domain_size(direction), Position_slice(direction)*voxel_size/1000, Domain_size(direction)*voxel_size/1000);
        Text_slider = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',str_slider,...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.65 0.33 0.2 0.03]);
        
        % Phase table
        Text_phasetable = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Choose phases name and color (RGB)',...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.5 0.28 0.5 0.05]);
        table_phase = uitable('Parent', segmented_grey_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.55 0.195 0.4 0.1],'CellEditCallback',@cellsection_tablePhase_Callback);
        table_phase.ColumnName = {'Phase','Volume fraction','Red','Green','Blue','Name'}; % Column name
        table_phase.ColumnEditable = [false false true true true true]; % Select column editable
        table_phase.ColumnWidth = {'auto', 'auto', 50, 50, 50, 100}; % Auto width
        table_phase.RowName = []; % Remove row name
        phase_data = [num2cell(Phase_code) num2cell(volumefraction) num2cell(RGB_red.index') num2cell(RGB_green.index') num2cell(RGB_blue.index') {Phase.name}'];
        set(table_phase,'Data',phase_data);
        
        % overlay_
        Text_overlay = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Overlay opacity:',...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.75-0.15 0.12 0.15 0.05]);
        Edit_overlay = uicontrol('Parent', segmented_grey_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',overlay_,...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.75 0.145 0.11 0.03],'Callback', @edit_overlay_Callback);
        
        % Save folder
        Button_savefolder = uicontrol('Parent', segmented_grey_figure, 'Style', 'pushbutton', 'String', 'Click to choose save folder','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
            'Position', [0.75-0.2 0.06 0.15 0.04],'UserData',0,...
            'Callback',{@button_savefolder_Callback});
        % Save figure
        Button_save_figure = uicontrol('Parent', segmented_grey_figure, 'Style', 'pushbutton', 'String', 'Save figure','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
            'Position', [0.75-0.2 0.02 0.075 0.04],'enable','off',...
            'Callback',{@button_savefigure_Callback});
        Button_save_video = uicontrol('Parent', segmented_grey_figure, 'Style', 'pushbutton', 'String', 'Save video','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized',...
            'Position', [0.75-0.125 0.02 0.075 0.04],'enable','off',...
            'Callback',{@button_savevideo_Callback});
        
        Text_filename = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Filename w/o extension:',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.71 0.075 0.12 0.02]);
        Text_videoquality = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Video quality (0-100):',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.71 0.045 0.12 0.02]);
        Text_videoframe = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Video slice per second:',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.71 0.015 0.12 0.02]);
        
        Edit_filename = uicontrol('Parent', segmented_grey_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',figure_initial_filename,...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.83 0.075 0.16 0.03]);
        Edit_videoquality = uicontrol('Parent', segmented_grey_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',video_quality,...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.83 0.045 0.16 0.03]);
        Edit_videoframe = uicontrol('Parent', segmented_grey_figure,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',video_slicepersecond,...
            'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.83 0.015 0.16 0.03]);
        
        
        Text_figure_files = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large2_GUI,'FontWeight','bold','FontName',font_name_GUI,'String','Files',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.55 0.45 0.4 0.025],'visible','off');
        line1 = annotation(segmented_grey_figure,'line',[0.55 0.95],[0.445 0.445],'visible','off');
        str_grey = sprintf('Grey-level image  : %s',FileName_grey);
        Text_figure_file_grey = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_grey,...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.41 0.4 0.03],'visible','off');
        str_segmented = sprintf('Segmented image : %s',FileName_segmented);
        Text_figure_file_segmented = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_segmented,...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.38 0.4 0.03],'visible','off');
        
        Text_figure_domain = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large2_GUI,'FontWeight','bold','FontName',font_name_GUI,'String','Domain',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.55 0.33 0.4 0.03],'visible','off');
        line2 = annotation(segmented_grey_figure,'line',[0.55 0.95],[0.325 0.325],'visible','off');
        str_domain = sprintf('Field of view: %1.2f x %1.2f x %1.2f um3',Domain_size(1)*voxel_size/1000,Domain_size(2)*voxel_size/1000,Domain_size(3)*voxel_size/1000);
        Text_figure_fov = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_domain,...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.29 0.4 0.03],'visible','off');
        str_voxel = sprintf('Voxel size %1.2f nm',voxel_size);
        Text_figure_voxel = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_voxel,...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.26 0.4 0.03],'visible','off');
        
        Text_figure_view = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large2_GUI,'FontWeight','bold','FontName',font_name_GUI,'String','View',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.55 0.21 0.4 0.03],'visible','off');
        line3 = annotation(segmented_grey_figure,'line',[0.55 0.95],[0.205 0.205],'visible','off');
        Text_figure_slider = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_slider,...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.17 0.4 0.03],'visible','off');
        str_viewnormal = sprintf('View normal to direction %i',direction);
        Text_figure_normal = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_viewnormal,...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.14 0.4 0.03],'visible','off');
        
        if number_phase<=6
            y_=0.12;
            for current_phase=1:1:number_phase
                y_=y_-0.019;
                ui_phase(current_phase).rectangle = annotation('rectangle',[0.55 y_ .1 .015],'FaceColor',[RGB_phase.index(current_phase).rgb(1) RGB_phase.index(current_phase).rgb(2) RGB_phase.index(current_phase).rgb(3)],'visible','off');
                str_ = sprintf('%s, volume fraction %1.3f',Phase(current_phase).name, volumefraction(current_phase,1));
                ui_phase(current_phase).text = uicontrol('Parent', segmented_grey_figure,'Style', 'text','FontSize',font_size_large_GUI-2,'FontName',font_name_GUI,'String',str_,...
                    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.66 y_-0.01 0.34 0.026],'visible','off');
            end
        end
        
        function button_savefolder_Callback(~,~)
            % Set string of the dialog box
            str_dialogbox = 'Select folder where your figure and video will be saved';
            % Open dialog box to choose file path
            Foldername = uigetdir(matlabroot,str_dialogbox);
            if Foldername==0
                % User clicked cancel button or closed the dialog box
                is_savefolder_selected=0;
            else
                % Save path
                is_savefolder_selected=1;
                Save_folder = Foldername;
            end
            if is_savefolder_selected==1
                set(Button_save_figure,'enable','on');
                set(Button_save_video,'enable','on');
            else
                set(Button_save_figure,'enable','off');
                set(Button_save_video,'enable','off');
            end
        end
        
        function GUI_off
            % Visible off GUI element
            set(Text_directionslider,'enable','off','visible','off');
            set(Popup_slider,'enable','off','visible','off');
            set(Slider_axes_2dview,'enable','off','visible','off');
            set(Text_slider,'enable','off','visible','off');
            set(Text_phasetable,'enable','off','visible','off');
            set(table_phase,'enable','off','visible','off');
            set(Text_overlay,'enable','off','visible','off');
            set(Edit_overlay,'enable','off','visible','off');
            set(Button_savefolder,'enable','off','visible','off');
            set(Button_save_figure,'enable','off','visible','off');
            set(Button_save_video,'enable','off','visible','off');
            set(Text_filename,'enable','off','visible','off');
            set(Text_videoquality,'enable','off','visible','off');
            set(Text_videoframe,'enable','off','visible','off');
            set(Edit_filename,'enable','off','visible','off');
            set(Edit_videoquality,'enable','off','visible','off');
            set(Edit_videoframe,'enable','off','visible','off');
            % Visible on information only for figure/video
            set(Text_figure_files,'visible','on');
            set(line1,'visible','on');
            set(Text_figure_file_grey,'visible','on');
            set(Text_figure_file_segmented,'visible','on');
            set(Text_figure_domain,'visible','on');
            set(line2,'visible','on');
            set(Text_figure_fov,'visible','on');
            set(Text_figure_voxel,'visible','on');
            set(Text_figure_view,'visible','on');
            set(line3,'visible','on');
            set(Text_figure_slider,'visible','on');
            set(Text_figure_normal,'visible','on');
            if number_phase<=6
                for current_phase=1:1:number_phase
                    set(ui_phase(current_phase).rectangle,'visible','on');
                    set(ui_phase(current_phase).text,'visible','on');
                end
            end
        end
        
        function GUI_on
            % Visible on GUI element
            set(Text_directionslider,'enable','on','visible','on');
            set(Popup_slider,'enable','on','visible','on');
            set(Slider_axes_2dview,'enable','on','visible','on');
            set(Text_slider,'enable','on','visible','on');
            set(Text_phasetable,'enable','on','visible','on');
            set(table_phase,'enable','on','visible','on');
            set(Text_overlay,'enable','on','visible','on');
            set(Edit_overlay,'enable','on','visible','on');
            set(Button_savefolder,'enable','on','visible','on');
            set(Button_save_figure,'enable','on','visible','on');
            set(Button_save_video,'enable','on','visible','on');
            set(Text_filename,'enable','on','visible','on');
            set(Text_videoquality,'enable','on','visible','on');
            set(Text_videoframe,'enable','on','visible','on');
            set(Edit_filename,'enable','on','visible','on');
            set(Edit_videoquality,'enable','on','visible','on');
            set(Edit_videoframe,'enable','on','visible','on');
            % Visible off information only for figure/video
            set(Text_figure_files,'visible','off');
            set(line1,'visible','off');
            set(Text_figure_file_grey,'visible','off');
            set(Text_figure_file_segmented,'visible','off');
            set(Text_figure_domain,'visible','off');
            set(line2,'visible','off');
            set(Text_figure_fov,'visible','off');
            set(Text_figure_voxel,'visible','off');
            set(Text_figure_view,'visible','off');
            set(line3,'visible','off');
            set(Text_figure_slider,'visible','off');
            set(Text_figure_normal,'visible','off');
            for current_phase=1:1:number_phase
                set(ui_phase(current_phase).rectangle,'visible','off');
                set(ui_phase(current_phase).text,'visible','off');
            end
        end
        
        function button_savefigure_Callback(~,~)
            GUI_off
            figure_filename=Edit_filename.String;
            function_savefig(segmented_grey_figure, [Save_folder '\'], figure_filename, OPTIONS);
            GUI_on
        end
        
        function button_savevideo_Callback(~,~)
            GUI_off
            figure_filename=Edit_filename.String;
            
            video_handle = VideoWriter([[Save_folder '\'] figure_filename],'mpeg-4');
            % Set video quality
            set(video_handle,'Quality',str2double(Edit_videoquality.String));
            % Set video framerate
            set(video_handle,'FrameRate',str2double(Edit_videoframe.String));
            % Open video
            open(video_handle)
            for frame_number=1:1:Domain_size(direction)
                % Update figure
                Position_slice(direction) = frame_number;
                str_slider = sprintf('Slice: %i/%i, %1.3f/%1.3f micrometers %s',Position_slice(direction), Domain_size(direction), Position_slice(direction)*voxel_size/1000, Domain_size(direction)*voxel_size/1000);
                set(Text_figure_slider,'String',str_slider);
                update_allfigures
                % Save frame
                stored_frame(frame_number) = getframe(segmented_grey_figure);
                writeVideo(video_handle,stored_frame(frame_number))
            end
            % Close video
            close(video_handle)
            GUI_on
        end
        
        % Create figure
        update_allfigures
        
        function update_allfigures
            rgb_array = update_RGB_figure(segmentedvolume,axes_segmented);
            grey_array = update_greylevel_figure(greyvolume,axes_greylevel);
            update_overlay_figure(rgb_array, grey_array, axes_overlay_);
        end
        
        function update_overlay_figure(rgb_array, grey_array, axe_)
            %C = imfuse(rgb_array,grey_array,'blend','Scaling','joint'); % Alternative
            % Convert grey image in RGB
            grey_in_RGB = cat(3, grey_array, grey_array, grey_array);
            grey_in_RGB=double(grey_in_RGB)/max_grey; % Scale
            % Linear combination
            C = overlay_*rgb_array + (1-overlay_)*grey_in_RGB;
            
            if remove_background==true
                channelR = C(:,:,1);
                channelG = C(:,:,2);
                channelB = C(:,:,3);
                channelR(back_ground(Position_slice(3)).idx)=double(grey_array(back_ground(Position_slice(3)).idx))/max_grey;
                channelG(back_ground(Position_slice(3)).idx)=double(grey_array(back_ground(Position_slice(3)).idx))/max_grey;
                channelB(back_ground(Position_slice(3)).idx)=double(grey_array(back_ground(Position_slice(3)).idx))/max_grey;
                C(:,:,1) = channelR;
                C(:,:,2) = channelG;
                C(:,:,3) = channelB;
            end
            
            % Display
            image(C,'Parent', axe_);
            % Get axis position
            axe_position =axe_.Position;
            set(axe_ ,'Position', axe_position);
            % Remove tick and label
            set(axe_,'xtick',[],'ytick',[]);
            set(axe_,'xticklabel',[],'yticklabel',[]);
            % Fit the axes box
            axis(axe_,'tight');
            % Aspect ratio is 1:1
            axis(axe_,'equal');
            axes_greylevel.YDir='normal';
            axes_segmented.YDir='normal';
            axes_overlay_.YDir='normal';            
        end
        
        function [slice_] = update_greylevel_figure(volume,axe_)
            if direction==1
                slice_ = squeeze(volume(Position_slice(1),:,:));
            elseif direction==2
                slice_ = squeeze(volume(:,Position_slice(2),:));
            elseif direction==3
                slice_ = volume(:,:,Position_slice(3));
            end
            image(slice_(:,:,1),'CDataMapping','scaled','Parent', axe_);
            colormap(gray)
            % Get axis position
            axe_position =axe_.Position;
            set(axe_ ,'Position', axe_position);
            % Remove tick and label
            set(axe_,'xtick',[],'ytick',[]);
            set(axe_,'xticklabel',[],'yticklabel',[]);
            % Fit the axes box
            axis(axe_,'tight');
            % Aspect ratio is 1:1
            axis(axe_,'equal');
            axes_greylevel.YDir='normal';
            axes_segmented.YDir='normal';
            axes_overlay_.YDir='normal';
        end
        
        function [slice_color] = update_RGB_figure(volume,axe_)
            if direction==1
                % Initializaion
                slice_color = zeros(Domain_size(2),Domain_size(3),3); % RGB color map
                slice_r = zeros(Domain_size(2),Domain_size(3)); % Red color map
                slice_g = zeros(Domain_size(2),Domain_size(3)); % Green color map
                slice_b = zeros(Domain_size(2),Domain_size(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:number_phase
                    code_tmp =Phase_code(current_phase); % Current phase code
                    slice_r(volume(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(volume(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(volume(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            elseif direction==2
                % Initializaion
                slice_color = zeros(Domain_size(1),Domain_size(3),3); % RGB color map
                slice_r = zeros(Domain_size(1),Domain_size(3)); % Red color map
                slice_g = zeros(Domain_size(1),Domain_size(3)); % Green color map
                slice_b = zeros(Domain_size(1),Domain_size(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:number_phase
                    code_tmp =Phase_code(current_phase); % Current phase code
                    slice_r(volume(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(volume(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(volume(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            elseif direction==3
                % Initializaion
                slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
                slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
                slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
                slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:number_phase
                    code_tmp =Phase_code(current_phase); % Current phase code
                    slice_r(volume(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(volume(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(volume(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            end
            slice_color(:,:,1)=slice_r; % Attribute RGB color
            slice_color(:,:,2)=slice_g;
            slice_color(:,:,3)=slice_b;
            % Display the slice
            image(slice_color,'parent',axe_);
            % Get axis position
            axe_position =axe_.Position;
            set(axe_ ,'Position', axe_position);
            % Remove tick and label
            set(axe_,'xtick',[],'ytick',[]);
            set(axe_,'xticklabel',[],'yticklabel',[]);
            % Fit the axes box
            axis(axe_,'tight');
            % Aspect ratio is 1:1
            axis(axe_,'equal');
            axes_greylevel.YDir='normal';
            axes_segmented.YDir='normal';
            axes_overlay_.YDir='normal';            
        end
        
        function edit_overlay_Callback(~,~)
            % Get voxel size
            overlay_=str2double(Edit_overlay.String);
            % Check entered value are correct
            if isnan(overlay_) || overlay_<0 || overlay_ >1
                overlay_=0.5;
                set(Edit_overlay,'String',overlay_);
            end
            update_allfigures
        end
        
        % Slider
        function slider_axes2dview_Callback(source,~)
            % Get position value
            pos_=round(source.Value);
            % Update position array
            Position_slice(direction)=pos_;
            % Update text
            str_slider = sprintf('Slice: %i/%i, %1.3f/%1.3f micrometers %s',Position_slice(direction), Domain_size(direction), Position_slice(direction)*voxel_size/1000, Domain_size(direction)*voxel_size/1000);
            set(Text_slider,'String',str_slider);
            set(Text_figure_slider,'String',str_slider);
            % Update figure
            update_allfigures
        end
        
        % Select direction
        function popup_axes2dview_Callback(source,~)
            % Get direction
            direction=source.Value;
            % Set slider min, max
            minor_step = 1/(Domain_size(direction)-1);
            major_step = 0.1;
            set(Slider_axes_2dview,'Min',1,'Max',Domain_size(direction),'SliderStep', [minor_step, major_step],'Value',1);
            % Update text
            str_slider = sprintf('Slice: %i/%i, %1.3f/%1.3f micrometers %s',Position_slice(direction), Domain_size(direction), Position_slice(direction)*voxel_size/1000, Domain_size(direction)*voxel_size/1000);
            set(Text_slider,'String',str_slider);
            set(Text_figure_slider,'String',str_slider);
            str_viewnormal = sprintf('View normal to direction %i',direction);
            set(Text_figure_normal,'String',str_viewnormal);
            % Update figure
            update_allfigures
        end
        
        function cellsection_tablePhase_Callback(~,~)
            % Get back data of the table
            table_data = table_phase.Data;
            % Check color
            Color_=table_data(:,3:5);
            Color_ = cell2mat(Color_);
            % Only number
            detect_nan = isnan(Color_);
            tmp=sum(sum(detect_nan));
            if tmp==0
                Color_hasnot_nan=1;
            else
                Color_hasnot_nan=0;
            end
            % Only between 0 and 1
            if Color_hasnot_nan==1
                tmp1=Color_; tmp2=Color_;
                tmp1(Color_>=0)=1; tmp2(Color_<=1)=1;
                tmp=tmp1.*tmp2;
                if (length(unique(tmp))==1 && unique(tmp)==1)
                    correct_color=1;
                else
                    correct_color=0;
                end
            else
                correct_color=0;
            end
            if correct_color==1
                % Update color
                color_phase=Color_;
                for current_phase=1:1:number_phase
                    RGB_red.index(current_phase)=color_phase(current_phase,1);
                    RGB_green.index(current_phase)=color_phase(current_phase,2);
                    RGB_blue.index(current_phase)=color_phase(current_phase,3);
                    RGB_phase.index(current_phase).rgb = [color_phase(current_phase,1) color_phase(current_phase,2) color_phase(current_phase,3)];
                end
                % Update figure
                update_allfigures
            end
            
            phase_name = table_data(:,6);
            for current_phase=1:1:number_phase
                str_ = sprintf('%s, volume fraction %1.3f',char(phase_name(current_phase)), volumefraction(current_phase,1));
                set(ui_phase(current_phase).text,'string',str_);
                set(ui_phase(current_phase).rectangle,'FaceColor',[RGB_phase.index(current_phase).rgb(1) RGB_phase.index(current_phase).rgb(2) RGB_phase.index(current_phase).rgb(3)]);
            end
        end
    end
end

