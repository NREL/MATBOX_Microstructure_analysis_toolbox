function microstructure_characterization_GUI
%Graphic user interface of the microstructure characterisation toolbox module

%%
%% GUI DISPLAY OPTIONS
%%
scrsz = get(0,'ScreenSize');

% These position varialbe are used in various GUI element.
% They are centralized here to assure visual coherence.

% Tab position
description_tab_fromleft = 0.025;
description_tab_frombottom = 0.925;
description_tab_xlenght = 1-2*description_tab_fromleft;
description_tab_ylenght = 0.045;
% Tab background and text color
background_description_tab = [0 0.5 0];
ForegroundColor_description_tab = [1 1 1]; 
background_tab = [1 1 1];

% Common position
pos_x_start = description_tab_fromleft;
pos_delta_x = 0.1;
pos_y_start = 0.85;
pos_delta_y = 0.1;

% Font
font_name_GUI ='Times New Roman';
% Font size
font_size_small_GUI =10;
font_size_large_GUI =14;


%%
%% PRESELECTED OPTIONS AND VARIABLES
%%

color_phase = [get(0, 'DefaultAxesColorOrder'); rand(100,3);]; % Color attributed to each phase. One row per phase.
number_phase=0;

number_volume_to_analyse = 0; % Number of volume to analyse
new_volume=1;
inputvolume=[];

Direction(1).name = 'Direction 1';
Direction(2).name = 'Direction 2';
Direction(3).name = 'Direction 3';
direction=1; % Microstructure view normal to direction
Position_slice=[1 1 1]; % Microstructure view slice position

initial_voxel_size=1; % Voxel size initialized
new_voxel_size=1; % Voxel upscale or downscale size initialized

ROI=0; % Region Of Interest (ROI) initialized
correct_ROI=true; % 
correct_assignedcode=true;
correct_ratio_voxel_size=true;
correct_folder=false;
correct_color=true;
correct_voxel_size=true;
resultfolder_location=false;

% RVE
RVE_choice = {'Independant subvolumes of same size + keep initial aspect ratio (A)',...
    'Independant subvolumes of same size + user-defined aspect ratio (B)',...
    'Independant subvolumes of same size + constant length (C)',...
    'One subvolume + growing from volume center (D)',...
    'One subvolume + growing from volume edge (E)'};

Slice_choice = {'Slice parameter is number of thick slice','Slice parameter is slice thickness (micrometers)'};


RVEparameters=[];
default_threshold_std = 5;
defaut_threshold_numbersubvolumes = 4;

% Default parameters for volume fractions 
PROPERTY=[];
PROPERTY.volumefractions.todo = true;
PROPERTY.volumefractions.voxel_size_dependence.todo = true;
PROPERTY.volumefractions.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.volumefractions.number_RVE = 1;
PROPERTY.volumefractions.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.volumefractions.RVE(1).savename = 'SubInitialAR';
PROPERTY.volumefractions.RVE(1).type = 'A';
PROPERTY.volumefractions.RVE(1).divisions = [2 3 4 5];
PROPERTY.volumefractions.RVE(1).subs2 = true;
PROPERTY.volumefractions.RVE(1).subs4 = true;
PROPERTY.volumefractions.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.volumefractions.RVE(1).Constantdirection = 'n/a';
PROPERTY.volumefractions.RVE(1).Growthdirection = 'n/a';
PROPERTY.volumefractions.RVE(1).Growthperstep = 'n/a';
PROPERTY.volumefractions.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.volumefractions.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.volumefractions.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.volumefractions.RVE(1).threshold_std = default_threshold_std;
PROPERTY.volumefractions.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

% Default parameters for specific surface area (direct method)
PROPERTY.specificsurfacearea_directmethod.todo = true;
PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.todo = true;
PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.specificsurfacearea_directmethod.number_RVE = 1;
PROPERTY.specificsurfacearea_directmethod.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.specificsurfacearea_directmethod.RVE(1).savename = 'SubInitialAR';
PROPERTY.specificsurfacearea_directmethod.RVE(1).type = 'A';
PROPERTY.specificsurfacearea_directmethod.RVE(1).divisions = [2 3 4 5];
PROPERTY.specificsurfacearea_directmethod.RVE(1).subs2 = true;
PROPERTY.specificsurfacearea_directmethod.RVE(1).subs4 = true;
PROPERTY.specificsurfacearea_directmethod.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.specificsurfacearea_directmethod.RVE(1).Constantdirection = 'n/a';
PROPERTY.specificsurfacearea_directmethod.RVE(1).Growthdirection = 'n/a';
PROPERTY.specificsurfacearea_directmethod.RVE(1).Growthperstep = 'n/a';
PROPERTY.specificsurfacearea_directmethod.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.specificsurfacearea_directmethod.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.specificsurfacearea_directmethod.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.specificsurfacearea_directmethod.RVE(1).threshold_std = default_threshold_std;
PROPERTY.specificsurfacearea_directmethod.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

specificsurfacearea_choice = {'Phase surface divided by domain volume','Phase surface divided by phase volume'};
default_specificsurfacearea_definition = 'Phase surface divided by domain volume';
default_correctivefactor_spdirect = 2/3;
PROPERTY.specificsurfacearea_directmethod.definition = default_specificsurfacearea_definition;
PROPERTY.specificsurfacearea_directmethod.correctivefactor = default_correctivefactor_spdirect;

% Default parameters for specific interface area (direct method)
PROPERTY.specificinterfacearea_directmethod.todo = true;
PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.todo = true;
PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.specificinterfacearea_directmethod.number_RVE = 1;
PROPERTY.specificinterfacearea_directmethod.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.specificinterfacearea_directmethod.RVE(1).savename = 'SubInitialAR';
PROPERTY.specificinterfacearea_directmethod.RVE(1).type = 'A';
PROPERTY.specificinterfacearea_directmethod.RVE(1).divisions = [2 3 4 5];
PROPERTY.specificinterfacearea_directmethod.RVE(1).subs2 = true;
PROPERTY.specificinterfacearea_directmethod.RVE(1).subs4 = true;
PROPERTY.specificinterfacearea_directmethod.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.specificinterfacearea_directmethod.RVE(1).Constantdirection = 'n/a';
PROPERTY.specificinterfacearea_directmethod.RVE(1).Growthdirection = 'n/a';
PROPERTY.specificinterfacearea_directmethod.RVE(1).Growthperstep = 'n/a';
PROPERTY.specificinterfacearea_directmethod.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.specificinterfacearea_directmethod.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.specificinterfacearea_directmethod.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.specificinterfacearea_directmethod.RVE(1).threshold_std = default_threshold_std;
PROPERTY.specificinterfacearea_directmethod.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

default_correctivefactor_spidirect = 2/3;
PROPERTY.specificinterfacearea_directmethod.correctivefactor = default_correctivefactor_spidirect;

% Defaut parameters for particle size (C-PSD)
PROPERTY.particlesize_cpsd.todo = true;
PROPERTY.particlesize_cpsd.voxel_size_dependence.todo = true;
PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.particlesize_cpsd.number_RVE = 1;
PROPERTY.particlesize_cpsd.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.particlesize_cpsd.RVE(1).savename = 'SubInitialAR';
PROPERTY.particlesize_cpsd.RVE(1).type = 'A';
PROPERTY.particlesize_cpsd.RVE(1).divisions = [2 3 4 5];
PROPERTY.particlesize_cpsd.RVE(1).subs2 = true;
PROPERTY.particlesize_cpsd.RVE(1).subs4 = true;
PROPERTY.particlesize_cpsd.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.particlesize_cpsd.RVE(1).Constantdirection = 'n/a';
PROPERTY.particlesize_cpsd.RVE(1).Growthdirection = 'n/a';
PROPERTY.particlesize_cpsd.RVE(1).Growthperstep = 'n/a';
PROPERTY.particlesize_cpsd.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.particlesize_cpsd.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.particlesize_cpsd.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.particlesize_cpsd.RVE(1).threshold_std = default_threshold_std;
PROPERTY.particlesize_cpsd.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

% Defaut parameters for particle size (distance map)
PROPERTY.particlesize_dmap.todo = true;
PROPERTY.particlesize_dmap.voxel_size_dependence.todo = true;
PROPERTY.particlesize_dmap.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.particlesize_dmap.number_RVE = 1;
PROPERTY.particlesize_dmap.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.particlesize_dmap.RVE(1).savename = 'SubInitialAR';
PROPERTY.particlesize_dmap.RVE(1).type = 'A';
PROPERTY.particlesize_dmap.RVE(1).divisions = [2 3 4 5];
PROPERTY.particlesize_dmap.RVE(1).subs2 = true;
PROPERTY.particlesize_dmap.RVE(1).subs4 = true;
PROPERTY.particlesize_dmap.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.particlesize_dmap.RVE(1).Constantdirection = 'n/a';
PROPERTY.particlesize_dmap.RVE(1).Growthdirection = 'n/a';
PROPERTY.particlesize_dmap.RVE(1).Growthperstep = 'n/a';
PROPERTY.particlesize_dmap.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.particlesize_dmap.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.particlesize_dmap.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.particlesize_dmap.RVE(1).threshold_std = default_threshold_std;
PROPERTY.particlesize_dmap.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;


% Defaut parameters for particle size (watershed)
PROPERTY.particlesize_watershed.todo = true;
PROPERTY.particlesize_watershed.voxel_size_dependence.todo = false;
PROPERTY.particlesize_watershed.voxel_size_dependence.voxel = [2 3 4 5];

% PROPERTY.particlesize_watershed.number_RVE = 1;
% PROPERTY.particlesize_watershed.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
% PROPERTY.particlesize_watershed.RVE(1).savename = 'SubInitialAR';
% PROPERTY.particlesize_watershed.RVE(1).type = 'A';
% PROPERTY.particlesize_watershed.RVE(1).divisions = [2 3 4 5];
% PROPERTY.particlesize_watershed.RVE(1).subs2 = true;
% PROPERTY.particlesize_watershed.RVE(1).subs4 = true;
% PROPERTY.particlesize_watershed.RVE(1).Aspectratio = 'n/a'; 
% PROPERTY.particlesize_watershed.RVE(1).Constantdirection = 'n/a';
% PROPERTY.particlesize_watershed.RVE(1).Growthdirection = 'n/a';
% PROPERTY.particlesize_watershed.RVE(1).Growthperstep = 'n/a';
% PROPERTY.particlesize_watershed.RVE(1).firstuniquevolume_size='n/a';
% PROPERTY.particlesize_watershed.RVE(1).firstuniquevolume_unit='n/a';
% PROPERTY.particlesize_watershed.RVE(1).Growthrelativeto = 'n/a';
% PROPERTY.particlesize_watershed.RVE(1).threshold_std = default_threshold_std;
% PROPERTY.particlesize_watershed.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

PROPERTY.particlesize_watershed.number_RVE = [];
PROPERTY.particlesize_watershed.RVE(1).name = []';
PROPERTY.particlesize_watershed.RVE(1).savename = [];
PROPERTY.particlesize_watershed.RVE(1).type = [];
PROPERTY.particlesize_watershed.RVE(1).divisions = [];
PROPERTY.particlesize_watershed.RVE(1).subs2 = [];
PROPERTY.particlesize_watershed.RVE(1).subs4 = [];
PROPERTY.particlesize_watershed.RVE(1).Aspectratio = []; 
PROPERTY.particlesize_watershed.RVE(1).Constantdirection = [];
PROPERTY.particlesize_watershed.RVE(1).Growthdirection = [];
PROPERTY.particlesize_watershed.RVE(1).Growthperstep = [];
PROPERTY.particlesize_watershed.RVE(1).firstuniquevolume_size=[];
PROPERTY.particlesize_watershed.RVE(1).firstuniquevolume_unit=[];
PROPERTY.particlesize_watershed.RVE(1).Growthrelativeto = [];
PROPERTY.particlesize_watershed.RVE(1).threshold_std = [];
PROPERTY.particlesize_watershed.RVE(1).threshold_numbersubvolumes = [];

PROPERTY.particlesize_watershed.customfunction = true;
PROPERTY.particlesize_watershed.cpsd_refining = true;
PROPERTY.particlesize_watershed.details_convergence = false;



% Defaut parameters for particle size (PCRF)
PROPERTY.particlesize_PCRF.todo = true;
PROPERTY.particlesize_PCRF.voxel_size_dependence.todo = false;
PROPERTY.particlesize_PCRF.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.particlesize_PCRF.number_RVE = [];
PROPERTY.particlesize_PCRF.RVE(1).name = []';
PROPERTY.particlesize_PCRF.RVE(1).savename = [];
PROPERTY.particlesize_PCRF.RVE(1).type = [];
PROPERTY.particlesize_PCRF.RVE(1).divisions = [];
PROPERTY.particlesize_PCRF.RVE(1).subs2 = [];
PROPERTY.particlesize_PCRF.RVE(1).subs4 = [];
PROPERTY.particlesize_PCRF.RVE(1).Aspectratio = []; 
PROPERTY.particlesize_PCRF.RVE(1).Constantdirection = [];
PROPERTY.particlesize_PCRF.RVE(1).Growthdirection = [];
PROPERTY.particlesize_PCRF.RVE(1).Growthperstep = [];
PROPERTY.particlesize_PCRF.RVE(1).firstuniquevolume_size=[];
PROPERTY.particlesize_PCRF.RVE(1).firstuniquevolume_unit=[];
PROPERTY.particlesize_PCRF.RVE(1).Growthrelativeto = [];
PROPERTY.particlesize_PCRF.RVE(1).threshold_std = [];
PROPERTY.particlesize_PCRF.RVE(1).threshold_numbersubvolumes = [];

PROPERTY.particlesize_PCRF.cpsd_refining = true;
PROPERTY.particlesize_PCRF.details_convergence = false;




% Defaut parameters for tortuosity factor (tau factor)
PROPERTY.tortuosity_taufactor.todo = true;
PROPERTY.tortuosity_taufactor.voxel_size_dependence.todo = true;
PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.tortuosity_taufactor.number_RVE = 1;
PROPERTY.tortuosity_taufactor.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.tortuosity_taufactor.RVE(1).savename = 'SubInitialAR';
PROPERTY.tortuosity_taufactor.RVE(1).type = 'A';
PROPERTY.tortuosity_taufactor.RVE(1).divisions = [2 3 4 5];
PROPERTY.tortuosity_taufactor.RVE(1).subs2 = true;
PROPERTY.tortuosity_taufactor.RVE(1).subs4 = true;
PROPERTY.tortuosity_taufactor.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.tortuosity_taufactor.RVE(1).Constantdirection = 'n/a';
PROPERTY.tortuosity_taufactor.RVE(1).Growthdirection = 'n/a';
PROPERTY.tortuosity_taufactor.RVE(1).Growthperstep = 'n/a';
PROPERTY.tortuosity_taufactor.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.tortuosity_taufactor.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.tortuosity_taufactor.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.tortuosity_taufactor.RVE(1).threshold_std = default_threshold_std;
PROPERTY.tortuosity_taufactor.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

PROPERTY.tortuosity_taufactor.todo_slice = [ {true}; {true}; {true}];
PROPERTY.tortuosity_taufactor.sliceparameter = [{3};{3};{3}];
PROPERTY.tortuosity_taufactor.slicechoice=1;

% Defaut parameters for connectivity
PROPERTY.connectivity.todo = true;
PROPERTY.connectivity.voxel_size_dependence.todo = true;
PROPERTY.connectivity.voxel_size_dependence.voxel = [2 3 4 5];

PROPERTY.connectivity.number_RVE = 1;
PROPERTY.connectivity.RVE(1).name = 'Independant subvolumes of same size + keep initial aspect ratio (A)';
PROPERTY.connectivity.RVE(1).savename = 'SubInitialAR';
PROPERTY.connectivity.RVE(1).type = 'A';
PROPERTY.connectivity.RVE(1).divisions = [2 3 4 5];
PROPERTY.connectivity.RVE(1).subs2 = true;
PROPERTY.connectivity.RVE(1).subs4 = true;
PROPERTY.connectivity.RVE(1).Aspectratio = 'n/a'; 
PROPERTY.connectivity.RVE(1).Constantdirection = 'n/a';
PROPERTY.connectivity.RVE(1).Growthdirection = 'n/a';
PROPERTY.connectivity.RVE(1).Growthperstep = 'n/a';
PROPERTY.connectivity.RVE(1).firstuniquevolume_size='n/a';
PROPERTY.connectivity.RVE(1).firstuniquevolume_unit='n/a';
PROPERTY.connectivity.RVE(1).Growthrelativeto = 'n/a';
PROPERTY.connectivity.RVE(1).threshold_std = default_threshold_std;
PROPERTY.connectivity.RVE(1).threshold_numbersubvolumes = defaut_threshold_numbersubvolumes;

PROPERTY.connectivity.todo_slice = [ {true}; {true}; {true}];
PROPERTY.connectivity.sliceparameter = [{3};{3};{3}];
PROPERTY.connectivity.slicechoice=1;


%%
%% CREATE GUI: OBJECT
%%

%% MAIN FIGURE
main_figure = figure; % Create figure
main_figure.Name= 'Microstructure characterization'; % Set name
main_figure.NumberTitle='off'; % Remove number from title name
main_figure.Color='white'; % Background colour
main_figure.MenuBar='none'; % Remove menubar and toolbar
main_figure.ToolBar='none';
main_figure.Units='normalized'; % Set unit
main_figure.Position=[0.1 0.1 0.75 0.8]; % Set Position

%% USER INTERFACE WITH TABBED PANELS
table_group_1 = uitabgroup('Parent', main_figure); % Create tabgroup
table_group_1.TabLocation='Left'; % Set panels location
% Create tabbed panels
tab_import = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Import segmented volume');
tab_imported = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Imported volume(s)');
tab_volume_fraction = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Volume fraction options');
tab_connectivity = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Connectivity options');
tab_tortuosity_taufactor = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Tortuosity factor (tau-factor) options');
tab_specificsurfacearea = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Specific surface area (direct) options');
tab_specificinterfacearea = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Specific interface area (direct) options');
%tab_covariance = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Surface and particle size (covariance) options');
tab_particlesize_CPSD = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Particle size (C-PSD) options');
tab_particlesize_dmap = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Particle size (EDMF) options');
tab_particlesize_watershed = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Particle size (watershed) options');
tab_particlesize_PCRF = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Particle size (PCRF) options');
tab_savedisplayoptions = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Display/save options');
tab_run = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Run ALL calculations');
%tab_help = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Help');

% Tab colors
tab_volume_fraction.ForegroundColor  ='b';
tab_connectivity.ForegroundColor  ='b';
tab_tortuosity_taufactor.ForegroundColor  ='b';
tab_specificsurfacearea.ForegroundColor  =[0 0.3922 0];
tab_specificinterfacearea.ForegroundColor  =[0 0.3922 0];
%tab_covariance.ForegroundColor  =[0 0.3922 0];
tab_particlesize_CPSD.ForegroundColor  =[1 0.549 0];
tab_particlesize_dmap.ForegroundColor  =[1 0.549 0];
tab_particlesize_watershed.ForegroundColor  =[1 0.549 0];
tab_particlesize_PCRF.ForegroundColor  =[1 0.549 0];
tab_run.ForegroundColor  ='r';

% Preselected tab
table_group_1.SelectedTab = tab_import;

%%
%% TAB IMPORT
%%

%% GUI

% Description and instructions
Text_tab_import = uicontrol('Parent', tab_import, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Select segmented volume file to import. You can select additional volumes after.');

% % SECTION FOR (ADD) NEW VOLUME
% Add volume push button
Button_add_volume = uicontrol('Parent', tab_import, 'Style', 'pushbutton', 'String', '>>>  Click to add volume  <<<',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [pos_x_start pos_y_start 3*pos_delta_x 0.5*pos_delta_y],'UserData',0,...
    'Callback',{@add_volume_Callback});
Text_new_volume = uicontrol('Parent', tab_import, 'Style', 'text','enable','off','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','No volume selected','UserData',0,...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos_x_start 0.79 3*pos_delta_x 0.5*pos_delta_y]);
Text_volumenumber = uicontrol('Parent', tab_import, 'Style', 'text','Visible','off','enable','off','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Volume numero: ','UserData',0,...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos_x_start 0.75 3*pos_delta_x 0.5*pos_delta_y]); 
 
% Figure for volume view
axes_2dview = axes('Parent', tab_import,'Visible','off','FontName',font_name_GUI,'Units','normalized','Position', [0 0.381 3.5*pos_delta_x 3.7*pos_delta_x]);
% Remove tick and label
set(axes_2dview,'xtick',[],'ytick',[]);
set(axes_2dview,'xticklabel',[],'yticklabel',[]);
% Box on
box(axes_2dview,'on');
% Fit the axes box
axis(axes_2dview,'tight');
% Aspect ratio is 1:1
axis(axes_2dview,'equal');
% Align
align([Button_add_volume axes_2dview Text_new_volume],'HorizontalAlignment','Center');
% Slider
Slider_axes_2dview = uicontrol('Parent', tab_import,'Style', 'slider','Min',1,'Max',100,'Value',1,'Units','normalized','Position', [pos_x_start 0.325 0.3 0.04],'Callback', @slider_axes2dview_Callback,...
                               'Visible','off','enable','off'); 
% Text (slider position)
Text_slider = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Position: -/-','Visible','off','enable','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos_x_start+0.15 0.245 1.1*pos_delta_x 0.5*pos_delta_y]);
% Popup menu (slider direction)
Popup_slider = uicontrol('Parent', tab_import,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', {'Direction 1','Direction 2','Direction 3'},'Units','normalized','Position', [pos_x_start 0.255 0.15 0.05],'enable','off','Visible','off','Callback', @popup_axes2dview_Callback);    
       
% Voxel size
% text
Text_voxelsize = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','1) Enter voxel size (in nanometers). New voxel size can be lower (downscaling) or higher (upscaling) compared with initial voxel size.','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.385 0.8 0.5365 0.1]);
Text_initialvoxelsize = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Initial voxel size (nm):','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.37 0.795 0.2 0.05]); 
Edit_initialvoxelsize = uicontrol('Parent', tab_import,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',initial_voxel_size,'Visible','off',...
     'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.55 0.815 0.075 0.04],'Callback', @edit_voxelsize_Callback); 
Text_newvoxelsize = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','New voxel size (nm):','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.64 0.795 0.15 0.05]); 
Edit_newvoxelsize = uicontrol('Parent', tab_import,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',new_voxel_size,'Visible','off',...
     'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.79 0.815 0.075 0.04],'Callback', @edit_voxelsize_Callback);  
 
% Axe/direction
% text
Text_directiontable = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','2) Direction properties: choose Region Of Interest (ROI) and name later used for figures','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.385 0.75 0.5365 0.045]);

 % Table ROI
table_input_ROI = uitable('Parent', tab_import,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.385 0.675 0.565 0.095],'CellEditCallback',@cellsection_tableROI_Callback, 'CellSelectionCallback',@cellsection_tableROI_Callback );
table_input_ROI.ColumnName = {'Direction','Number of voxel','ROI start','ROI end','Name'}; % Column name
table_input_ROI.ColumnEditable = [false false true true true]; % Select column editable
table_input_ROI.ColumnWidth = {'auto', 'auto', 'auto', 'auto', 241}; % Auto width
table_input_ROI.RowName = []; % Remove row name
input_ROI_data = [num2cell([1;2;3]) num2cell([1;2;3]) num2cell([1;2;3]) num2cell([1;2;3]) {Direction(1).name; Direction(2).name; Direction(3).name}]; 
table_input_ROI.Data=input_ROI_data;
Text_InitialROI = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Initial volume and voxel size','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.59 0.62 0.25 0.04]);
Text_Numberofvoxel = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Number of voxels','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.39 0.60 0.25 0.04]);
Text_Domainsize = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Domain size','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.39 0.58 0.4 0.04]);
Text_Equivalentdomainsize = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Equivalent cubic domain size','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.39 0.56 0.2 0.04]);
Text_Numberofvoxel_bis = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',':','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.60 0.25 0.04]);
Text_Domainsize_bis = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',':','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.58 0.4 0.04]);
Text_Equivalentdomainsize_bis = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',':','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.56 0.56 0.2 0.04]);
Text_Numberofvoxel_values = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.59 0.60 0.25 0.04]);
Text_Domainsize_values = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.59 0.58 0.4 0.04]);
Text_Equivalentdomainsize_values = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.59 0.56 0.2 0.04]); 
Text_userROI = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','User-choice (ROI)','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.785 0.62 0.25 0.04]); 
Text_Numberofvoxel_uservalues = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.785 0.60 0.25 0.04]);
Text_Domainsize_uservalues = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.785 0.58 0.4 0.04]);
Text_Equivalentdomainsize_uservalues = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.785 0.56 0.2 0.04]); 

% Phase
Text_phasetable = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','3) Phase properties: choose name and color (RGB) for future figures','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.385 0.525 0.5365 0.04]);
table_input_phase = uitable('Parent', tab_import,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.385 0.39 0.565 0.15],'CellEditCallback',@cellsection_tablePhase_Callback);
table_input_phase.ColumnName = {'Phase','Volume fraction','Red','Green','Blue','Name','Assign code'}; % Column name
table_input_phase.ColumnEditable = [false false true true true true true]; % Select column editable
table_input_phase.ColumnWidth = {'auto', 'auto', 50, 50, 50, 168, 'auto'}; % Auto width
table_input_phase.RowName = []; % Remove row name

Text_ErrorColor= uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Error: RGB colors must be >=0 and <=1','Visible','off',...
     'BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.077 0.5 0.2 0.1]);
align([Text_ErrorColor axes_2dview],'HorizontalAlignment','Center');

Text_volumeinfo = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','4) Write volume information','Visible','off',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.385 0.315 0.5365 0.05]);
table_input_info = uitable('Parent', tab_import,'enable','on','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.385 0.11 0.565 0.23]);
table_input_info.ColumnName = {'Information','Text'}; % Column name
table_input_info.ColumnEditable = [false true]; % Select column editable
table_input_info.RowName = []; % Remove row name
inputinfo(1).name = 'Material source'; inputinfo(1).text = 'Unknown';
inputinfo(2).name = 'Material name'; inputinfo(2).text = 'Unknown';
inputinfo(3).name = 'Material group name'; inputinfo(3).text = 'Unknown';
inputinfo(4).name = 'Volume description'; inputinfo(4).text = 'No volume description';
inputinfo(5).name = 'Raw data source'; inputinfo(5).text = 'Unknown';
inputinfo(6).name = 'Raw data technique'; inputinfo(6).text = 'Unknown';
inputinfo(7).name = 'Segmentation source'; inputinfo(7).text = 'Unknown';
inputinfo(8).name = 'Segmentation technique'; inputinfo(8).text = 'Unknown';
inputinfo(9).name = 'Analysis decription'; inputinfo(9).text = 'No analysis description';
table_input_info.ColumnWidth = {140, 415}; % Auto width
input_info_data = [{inputinfo(1).name; inputinfo(2).name; inputinfo(3).name; inputinfo(4).name; inputinfo(5).name; inputinfo(6).name; inputinfo(7).name; inputinfo(8).name; inputinfo(9).name},...
                   {inputinfo(1).text; inputinfo(2).text; inputinfo(3).text; inputinfo(4).text; inputinfo(5).text; inputinfo(6).text; inputinfo(7).text; inputinfo(8).text; inputinfo(9).text}];
set(table_input_info,'Data',input_info_data);

Text_save_volume = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','5) Select location where the result folder will be saved, and write the result folder name','Visible','off',...
      'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.385 0.035 0.5365 0.05]);
Button_select_savefolder = uicontrol('Parent', tab_import, 'Style', 'pushbutton', 'String', 'Select result folder location','Visible','off','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
     'enable','off','Units','normalized','Position', [0.375 0.02 0.2 0.04],'Callback',{@select_savefolder_Callback}); 
Edit_resultfoldername = uicontrol('Parent', tab_import,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','-','Visible','off',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.58 0.0225 0.37 0.04]);
Text_Voxelsize_statut = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Voxel size: CORRECT','Visible','off',...
    'BackgroundColor','w','HorizontalAlignment','center','ForegroundColor',[34/255 177/255 76/255],'Units','normalized','Position', [0.05 0.2 0.25 0.04]);
Text_ROI_statut = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Region of interest (ROI): CORRECT','Visible','off',...
    'BackgroundColor','w','HorizontalAlignment','center','ForegroundColor',[34/255 177/255 76/255],'Units','normalized','Position', [0.05 0.175 0.25 0.04]);
Text_Color_statut = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Colors: CORRECT','Visible','off',...
    'BackgroundColor','w','HorizontalAlignment','center','ForegroundColor',[34/255 177/255 76/255],'Units','normalized','Position', [0.05 0.15 0.25 0.04]);
Text_Assigncode_statut = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Assigned code: CORRECT','Visible','off',...
    'BackgroundColor','w','HorizontalAlignment','center','ForegroundColor',[34/255 177/255 76/255],'Units','normalized','Position', [0.05 0.125 0.25 0.04]);
Text_Folder_statut = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Folder location: NOT DEFINED','Visible','off',...
    'BackgroundColor','w','HorizontalAlignment','center','ForegroundColor','r','Units','normalized','Position', [0.05 0.1 0.25 0.04]);
Button_Oncefinisehd = uicontrol('Parent', tab_import,'Style', 'pushbutton','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','6) Once steps 1-5 done, click to save your choices','Visible','off','enable','off',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.07 0.3 0.04],'Callback',{@Savevolume_and_options_Callback});
Button_Reset = uicontrol('Parent', tab_import,'Style', 'pushbutton','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Reset selected options and current volume','Visible','off','enable','off',...
    'BackgroundColor',[255 178 102]/255,'HorizontalAlignment','center','Units','normalized','Position', [0.025 0.02 0.3 0.04],'Callback',{@Reset_Callback});
  
A{1,1} = 'Volume information has been saved';                                                           
A{1,2} = ' You can add new volume or move to other tabs:';
A{1,3} = ' - ''Imported volume(s)'' tab: see volume(s) that will be investigate';
A{1,4} = ' - ''* options'' tabs: choose what property to calculate, and how';
A{1,5} = ' - ''Display/save options'' tab: choose your display and save options';
A{1,6} = ' - ''Run ALL calculations'' tab: once you have selected all your options, run ALL calculations for ALL volumes';
A{1,7} = ' - ''Help'' tab: link to module documentation';
multiple_line_string = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',A{1,1},A{1,2},A{1,3},A{1,4},A{1,5},A{1,6},A{1,7}); 
Text_saveconfirmation = uicontrol('Parent', tab_import,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String',multiple_line_string,'Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.15 0.1 0.75 0.6]); 
 
% Align
align([Text_Voxelsize_statut Text_ROI_statut Text_Color_statut Text_Folder_statut],'HorizontalAlignment','Center');  
 % Align
align([Text_voxelsize Text_directiontable Text_phasetable Text_save_volume],'HorizontalAlignment','Center');

%% CALLBACKS

% Add volume
    function add_volume_Callback(~,~)
        % Clear memory
        clear current_volume Phase
        % Set string of the dialog box
        str_dialogbox = ['Select volume numero ' num2str(number_volume_to_analyse+1)]; 
        % Open dialog box to choose file path
        [FileName,PathName,~] = uigetfile({'*.tif;*.tiff','Tif image (*.tif, *.tiff)';'*.m;*.mat','MATLAB File (*.m,*.mat)'},str_dialogbox);
        if FileName==0
            % User clicked cancel button or closed the dialog box
        else
            full_path_new_volume = [PathName FileName]; % Full path
            [new_volume, outcome] = function_loadvolume(full_path_new_volume, 'uint8', 'none' ); % Load new volume
            if outcome.success % Success to import
                % Update GUI text field and user data, get remove button enable, update number of volume
                Text_new_volume.String=FileName;
                % Activate related ui
                set(Text_new_volume,'enable','on');
                set(Text_slider,'enable','on','Visible','on');
                set(Popup_slider,'enable','on','Visible','on');
                set(Slider_axes_2dview,'enable','on','Visible','on');
                set(Text_volumenumber,'enable','on','Visible','on');
                set(Text_voxelsize,'Visible','on');
                set(Text_initialvoxelsize,'Visible','on');
                set(Text_newvoxelsize,'Visible','on');
                set(Edit_initialvoxelsize,'Visible','on','enable','on');
                set(Edit_newvoxelsize,'Visible','on','enable','on');
                set(Button_add_volume,'enable','off');
                set(axes_2dview,'Visible','on');
                set(Text_volumeinfo,'visible','on');
                set(table_input_info,'visible','on');
                
                % Update volume number
                Text_volumenumber.String=['Volume numero: ' num2str(number_volume_to_analyse+1)];
                % Remove save confirmation
                set(Text_saveconfirmation,'visible','off');
                % Get domain size
                Domain_size=size(new_volume);
                % Get phase code
                Phase_code=unique(new_volume);
                % Get number of phase
                number_phase=length(Phase_code);
                
                % Save in a structure information on the volume
                % Save them in the userdata field of this ui object (thus shared memory between all other ui objects)
                current_volume.new_volume=new_volume;
                current_volume.Domain_size=Domain_size;
                current_volume.Phase_code=Phase_code;
                current_volume.number_phase=number_phase;
                current_volume.full_path_new_volume=full_path_new_volume;
                set(Button_add_volume,'UserData',current_volume);

                % Get volume fraction
                total_number_voxel=prod(Domain_size);
                volumefraction=zeros(number_phase,1);
                for n=1:1:number_phase
                    volumefraction(n,1)=sum(sum(sum(new_volume==Phase_code(n))))/total_number_voxel;
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

                % Set slider min, max
                minor_step = 1/(Domain_size(direction)-1);
                major_step = 0.1;
                set(Slider_axes_2dview,'Min',1,'Max',Domain_size(direction),'SliderStep', [minor_step, major_step],'Value',1);
                % Update text
                set(Text_slider,'String',['Position: ' num2str(1) '/' num2str(Domain_size(direction))]);

                % Update direction table
                set(Text_directiontable,'visible','on');
                set(table_input_ROI,'visible','on');set(table_input_ROI,'enable','on');
                input_ROI_data = [num2cell([1;2;3]) num2cell([Domain_size(1);Domain_size(2);Domain_size(3)]) num2cell([1;1;1]) num2cell([Domain_size(1);Domain_size(2);Domain_size(3)]) {Direction(1).name; Direction(2).name; Direction(3).name}]; 
                set(table_input_ROI,'Data',input_ROI_data);
                % Update phase table
                set(Text_phasetable,'visible','on');
                set(table_input_phase,'visible','on');set(table_input_phase,'enable','on');
                input_phase_data = [num2cell(Phase_code) num2cell(volumefraction) num2cell(RGB_red.index') num2cell(RGB_green.index') num2cell(RGB_blue.index') {Phase.name}' num2cell(Phase_code)]; 
                set(table_input_phase,'Data',input_phase_data);
                % Align
                align([Text_directiontable Text_phasetable],'HorizontalAlignment','Center');

                % Update domain size text
                size_direction_1 = Domain_size(1)*initial_voxel_size/1000;
                size_direction_2 = Domain_size(2)*initial_voxel_size/1000;
                size_direction_3 = Domain_size(3)*initial_voxel_size/1000;
                size_cubeequivalent = (size_direction_1*size_direction_2*size_direction_3)^(1/3);
                str_1 = sprintf('%0.1f',size_direction_1);
                str_2 = sprintf('%0.1f',size_direction_2);
                str_3 = sprintf('%0.1f',size_direction_3);
                str_4 = sprintf('%0.1f',size_cubeequivalent);
                str_sav = sprintf('%0.0f',size_cubeequivalent);

                % Update ROI
                ROI=[1 Domain_size(1);1 Domain_size(2);1 Domain_size(3)];

                set(Text_Numberofvoxel_values,'String',num2str(total_number_voxel));
                set(Text_Domainsize_values,'String',[str_1 ' x ' str_2 ' x ' str_3 ' um3']);
                set(Text_Equivalentdomainsize_values,'String',[str_4 ' x ' str_4 ' x ' str_4 ' um3']);
                align([Text_InitialROI Text_Numberofvoxel_values Text_Domainsize_values Text_Equivalentdomainsize_values],'HorizontalAlignment','Center');
                set(Text_InitialROI,'Visible','on'); 
                set(Text_Numberofvoxel_values,'Visible','on'); set(Text_Domainsize_values,'Visible','on');  set(Text_Equivalentdomainsize_values,'Visible','on');
                set(Text_Numberofvoxel,'Visible','on'); set(Text_Numberofvoxel_bis,'Visible','on');
                set(Text_Domainsize,'Visible','on'); set(Text_Domainsize_bis,'Visible','on');
                set(Text_Equivalentdomainsize,'Visible','on'); set(Text_Equivalentdomainsize_bis,'Visible','on');

                set(Text_Numberofvoxel_uservalues,'String',num2str(total_number_voxel));
                set(Text_Domainsize_uservalues,'String',[str_1 ' x ' str_2 ' x ' str_3 ' um3']);
                set(Text_Equivalentdomainsize_uservalues,'String',[str_4 ' x ' str_4 ' x ' str_4 ' um3']);
                align([Text_userROI Text_Numberofvoxel_uservalues Text_Domainsize_uservalues Text_Equivalentdomainsize_uservalues],'HorizontalAlignment','Center');
                set(Text_userROI,'Visible','on'); 
                set(Text_Numberofvoxel_uservalues,'Visible','on');
                set(Text_Domainsize_uservalues,'Visible','on');
                set(Text_Equivalentdomainsize_uservalues,'Visible','on');            

                % Update save
                % Select filename without extension
                [~, name_, ~] = fileparts(Text_new_volume.String);
                savefolder_name=[name_ '_Vol_' str_sav 'um_Voxel_1nm'];
                set(Edit_resultfoldername,'String',savefolder_name);
                set(Text_save_volume,'Visible','on');
                set(Button_select_savefolder,'Visible','on','enable','on');
                set(Edit_resultfoldername,'Visible','on','enable','on');

                % Update statut
                set(Text_Voxelsize_statut,'Visible','on');
                set(Text_ROI_statut,'Visible','on');
                set(Text_Color_statut,'Visible','on');
                set(Text_Assigncode_statut,'Visible','on');
                set(Text_Folder_statut,'Visible','on');
                set(Button_Oncefinisehd,'Visible','on');
                set(Button_Reset,'Visible','on','enable','on');

                % Update figure
                update_figure
            end
        end
    end

% Phase table
    function cellsection_tablePhase_Callback(~,~)
        % Get back data of the table
        table_data = table_input_phase.Data;

        % % Check assigned code
        Assign_code=table_data(:,7);
        Assign_code = cell2mat(Assign_code);
        unique_assign = unique(Assign_code);
        % Only number
        detect_nan = isnan(unique_assign); 
        tmp=sum(sum(detect_nan));
        if tmp==0
            Codeassing_hasnot_nan=1;
        else
            Codeassing_hasnot_nan=0;
        end
        % Superior or equal to 0
        if Codeassing_hasnot_nan==1
            if min(unique_assign)>=0
                correct_code_sup0=1;
            else
                correct_code_sup0=0;
            end
        else
            correct_code_sup0=0;
        end
        % Only integer
        if correct_code_sup0==1
            check_integer = sum(mod(unique_assign,1));
            if check_integer==0
                correct_code_int=1;
            else
                correct_code_int=0;
            end
        else
            correct_code_int=0;
        end        
        correct_assignedcode = Codeassing_hasnot_nan*correct_code_sup0*correct_code_int;
     
        if correct_assignedcode==1
            % Update display
            set(Text_Assigncode_statut,'ForegroundColor',[34/255 177/255 76/255],'String','Assigned code: CORRECT');
        else
            % Update display
            set(Text_Assigncode_statut,'ForegroundColor','r','String','Assigned code: ERROR (>=0 integers)');
        end        

        % % Check color
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
            % Update display
            set(Text_Color_statut,'ForegroundColor',[34/255 177/255 76/255],'String','Colors: CORRECT');
            set(Text_ErrorColor,'Visible','off');
            set(axes_2dview,'Visible','on');
            set(axes_2dview.Children,'Visible','on');
            % Update figure color
            color_phase=Color_;
            update_figure
        else
            % Update display
            set(Text_Color_statut,'ForegroundColor','r','String','Colors: ERROR');
            set(Text_ErrorColor,'Visible','on');
            set(axes_2dview,'Visible','off');
            set(axes_2dview.Children,'Visible','off');
        end
        
        % Check save is possible
        save_possible =  correct_color*correct_assignedcode*correct_ROI*correct_ratio_voxel_size*correct_folder;
        if save_possible==1
            set(Button_Oncefinisehd,'enable','on');
        else
            set(Button_Oncefinisehd,'enable','off');
        end
    end

% ROI table
    function cellsection_tableROI_Callback(~,~)
        % Get back data of the table
        table_data = table_input_ROI.Data;
        ROI_=table_data(:,3:4);
        ROI_ = cell2mat(ROI_);
        max_roiend=table_data(:,2);
        max_roiend=cell2mat(max_roiend);
        
        % % Check validity of user-defined ROI
        % Only number
        detect_nan = isnan(ROI_);
        tmp=sum(sum(detect_nan));
        if tmp==0
            ROI_hasnot_nan=1;
        else
            ROI_hasnot_nan=0;
        end
        % Only integer
        integer_ROI = round(ROI_);
        tmp = integer_ROI-ROI_;
        tmp2=nansum(nansum(tmp));
        if tmp2==0
            ROI_has_integer=1;
        else
            ROI_has_integer=0;
        end
        % Only positive
        if ROI_hasnot_nan==1
            tmp=ROI_;
            tmp(ROI_>0)=1;
            if (length(unique(tmp))==1 && unique(tmp)==1)
                ROI_has_onlypositive=1;
            else
                ROI_has_onlypositive=0;
            end
        else
            ROI_has_onlypositive=0;
        end
        % Correctly ordered
        if ROI_hasnot_nan==1
            tmp=ROI_(:,2)-ROI_(:,1);
            tmp2=ROI_(:,1)*0;
            tmp2(tmp>=0)=1;
            if (length(unique(tmp2))==1 && unique(tmp2)==1)
                ROI_has_correctorder=1;
            else
                ROI_has_correctorder=0;
            end
        else
            ROI_has_correctorder=0;
        end
        % Within field of view
        if ROI_hasnot_nan==1
            tmp=max_roiend-ROI_(:,2);
            tmp2=max_roiend*0;
            tmp2(tmp>=0)=1;
            if (length(unique(tmp2))==1 && unique(tmp2)==1)
                ROI_withinfieldview=1;
            else
                ROI_withinfieldview=0;
            end
        else
            ROI_withinfieldview=0;
        end
        % Thus, is the userdefined ROI is correct?
        correct_ROI = ROI_hasnot_nan*ROI_has_integer*ROI_has_onlypositive*ROI_has_correctorder*ROI_withinfieldview;
        
        if correct_ROI==1
            set(Text_ROI_statut,'ForegroundColor',[34/255 177/255 76/255],'String','Region of interest (ROI): CORRECT');
            % Update domain size
            ROI=ROI_;
            if correct_ratio_voxel_size==1
                % Domain size
                new_domain_size=ROI(:,2)-ROI(:,1)+1;
                % Number of voxel
                ROI_number_of_voxel=prod(new_domain_size);
                ROIvoxelsize_number_of_voxel = floor(ROI_number_of_voxel*((initial_voxel_size/new_voxel_size)^3));
                % Update domain size text
                size_direction_1 = new_domain_size(1)*initial_voxel_size/1000;
                size_direction_2 = new_domain_size(2)*initial_voxel_size/1000;
                size_direction_3 = new_domain_size(3)*initial_voxel_size/1000;
                size_cubeequivalent = (size_direction_1*size_direction_2*size_direction_3)^(1/3);
                str_1 = sprintf('%0.1f',size_direction_1);
                str_2 = sprintf('%0.1f',size_direction_2);
                str_3 = sprintf('%0.1f',size_direction_3);
                str_4 = sprintf('%0.1f',size_cubeequivalent);
                % Update display
                set(Text_userROI,'ForegroundColor','k');
                set(Text_Numberofvoxel_uservalues,'ForegroundColor','k','String',['~ ' num2str(ROIvoxelsize_number_of_voxel)]);
                set(Text_Domainsize_uservalues,'ForegroundColor','k','String',['~ ' str_1 ' x ' str_2 ' x ' str_3 ' um3']);
                set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','k','String',['~ ' str_4 ' x ' str_4 ' x ' str_4 ' um3']);
                
                % Update save
                % Select filename without extension
                [~, name_, ~] = fileparts(Text_new_volume.String);
                % Equivalent volume
                str_sav1 = sprintf('%0.0f',size_cubeequivalent);
                % New voxel size
                str_sav2 = sprintf('%0.0f',new_voxel_size);
                % Edit
                savefolder_name=[name_ '_Vol_' str_sav1 'um_Voxel_' str_sav2 'nm'];
                set(Edit_resultfoldername,'String',savefolder_name);
                               
            else
                % Update display
                set(Text_userROI,'ForegroundColor','r');
                set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','Voxel size error');
                set(Text_Domainsize_uservalues,'ForegroundColor','r','String','Voxel size error');
                set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','Voxel size error');
            end
        else
            % Update display
            set(Text_userROI,'ForegroundColor','r');
            set(Text_ROI_statut,'ForegroundColor','r','String','Region of interest (ROI): ERROR');
            if correct_ratio_voxel_size==1
                set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','ROI error');
                set(Text_Domainsize_uservalues,'ForegroundColor','r','String','ROI error');
                set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','ROI error');
            else
                set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                set(Text_Domainsize_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
            end
            % Update domain size
            ROI=0;
        end
        
        % Check save is possible
        save_possible =  correct_color*correct_ROI*correct_ratio_voxel_size*correct_folder;
        if save_possible==1
            set(Button_Oncefinisehd,'enable','on');
        else
            set(Button_Oncefinisehd,'enable','off');
        end        
        
    end

% Voxel size
    function edit_voxelsize_Callback(~,~)
        % Get initial voxel size
        initial_voxel_size=str2double(Edit_initialvoxelsize.String);
        % Get new voxel size
        new_voxel_size=str2double(Edit_newvoxelsize.String);
        % Get initial domain size
        Domain_size=Button_add_volume.UserData.Domain_size;        
        
        % Check entered value are correct
        if isnan(initial_voxel_size) || isnan(new_voxel_size) || initial_voxel_size<=0 || new_voxel_size<=0
            % Wrong entry
            correct_voxel_size=1;
            set(Text_Voxelsize_statut,'ForegroundColor','r','String','Voxel size: ERROR');
            if correct_ROI==1
                set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','Voxel size error');
                set(Text_Domainsize_uservalues,'ForegroundColor','r','String','Voxel size error');
                set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','Voxel size error');
            else
                set(Text_userROI,'ForegroundColor','r');
                set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                set(Text_Domainsize_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
            end
        else
            % Correct entry
            correct_voxel_size=1;
            % Update
            if correct_voxel_size==1
                set(Text_Voxelsize_statut,'ForegroundColor',[34/255 177/255 76/255],'String','Voxel size: CORRECT');
                if correct_ROI==1
                    % Domain size
                    new_domain_size=ROI(:,2)-ROI(:,1)+1;
                    % Number of voxel
                    ROI_number_of_voxel=prod(new_domain_size);
                    ROIvoxelsize_number_of_voxel = floor(ROI_number_of_voxel*((initial_voxel_size/new_voxel_size)^3));
                     % Update domain size text
                    size_direction_1 = new_domain_size(1)*initial_voxel_size/1000;
                    size_direction_2 = new_domain_size(2)*initial_voxel_size/1000;
                    size_direction_3 = new_domain_size(3)*initial_voxel_size/1000;
                    size_cubeequivalent = (size_direction_1*size_direction_2*size_direction_3)^(1/3);
                    str_1 = sprintf('%0.1f',size_direction_1);
                    str_2 = sprintf('%0.1f',size_direction_2);
                    str_3 = sprintf('%0.1f',size_direction_3);
                    str_4 = sprintf('%0.1f',size_cubeequivalent);
                    set(Text_userROI,'ForegroundColor','k');
                    set(Text_Numberofvoxel_uservalues,'ForegroundColor','k','String',['~ ' num2str(ROIvoxelsize_number_of_voxel)]);
                    set(Text_Domainsize_uservalues,'ForegroundColor','k','String',['~ ' str_1 ' x ' str_2 ' x ' str_3 ' um3']);
                    set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','k','String',['~ ' str_4 ' x ' str_4 ' x ' str_4 ' um3']);
                    
                    % Update save
                    % Select filename without extension
                    [~, name_, ~] = fileparts(Text_new_volume.String);
                    % Equivalent volume
                    str_sav1 = sprintf('%0.0f',size_cubeequivalent);
                    % New voxel size
                    str_sav2 = sprintf('%0.0f',new_voxel_size);
                    % Edit
                    savefolder_name=[name_ '_Vol_' str_sav1 'um_Voxel_' str_sav2 'nm'];
                    set(Edit_resultfoldername,'String',savefolder_name);
                else
                    set(Text_userROI,'ForegroundColor','r');
                    set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','ROI error');
                    set(Text_Domainsize_uservalues,'ForegroundColor','r','String','ROI error');
                    set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','ROI error');
                end
            else
                set(Text_Voxelsize_statut,'ForegroundColor','r','String','Voxel size: ERROR');
                if correct_ROI==1
                    set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','Voxel size error');
                    set(Text_Domainsize_uservalues,'ForegroundColor','r','String','Voxel size error');
                    set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','Voxel size error');
                else
                    set(Text_userROI,'ForegroundColor','r');
                    set(Text_Numberofvoxel_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                    set(Text_Domainsize_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                    set(Text_Equivalentdomainsize_uservalues,'ForegroundColor','r','String','ROI and voxel size error');
                end
            end
        end
        % Check entered value are correct
        if ~isnan(initial_voxel_size) && initial_voxel_size>0
            % Update initial voxel size display
            size_direction_1 = Domain_size(1)*initial_voxel_size/1000;
            size_direction_2 = Domain_size(2)*initial_voxel_size/1000;
            size_direction_3 = Domain_size(3)*initial_voxel_size/1000;
            size_cubeequivalent = (size_direction_1*size_direction_2*size_direction_3)^(1/3);
            str_1 = sprintf('%0.1f',size_direction_1);
            str_2 = sprintf('%0.1f',size_direction_2);
            str_3 = sprintf('%0.1f',size_direction_3);
            str_4 = sprintf('%0.1f',size_cubeequivalent);
            set(Text_Domainsize_values,'ForegroundColor','k','String',[str_1 ' x ' str_2 ' x ' str_3 ' um3']);
            set(Text_Equivalentdomainsize_values,'ForegroundColor','k','String',[str_4 ' x ' str_4 ' x ' str_4 ' um3']);
        else
            set(Text_Voxelsize_statut,'ForegroundColor','r','String','Voxel size: ERROR');
            set(Text_Domainsize_values,'ForegroundColor','r','String','Voxel size error');
            set(Text_Equivalentdomainsize_values,'ForegroundColor','r','String','Voxel size error');
        end
        
        % Check save is possible
        save_possible =  correct_color*correct_ROI*correct_voxel_size*correct_folder;
        if save_possible==1
            set(Button_Oncefinisehd,'enable','on');
        else
            set(Button_Oncefinisehd,'enable','off');
        end        
        
    end

    % Slider
    function slider_axes2dview_Callback(source,~)
        % Get position value
        pos_=round(source.Value);
        % Get domain size
        Domain_size=Button_add_volume.UserData.Domain_size;            
        % Update text
        set(Text_slider,'String',['Position: ' num2str(pos_) '/' num2str(Domain_size(direction))]);
        % Update position array
        Position_slice(direction)=pos_;
        % Update figure
        update_figure
    end

    % Select direction
    function popup_axes2dview_Callback(source,~)
        % Get direction
        direction=source.Value;
        % Get domain size
        Domain_size=Button_add_volume.UserData.Domain_size;        
        % Set slider min, max
        minor_step = 1/(Domain_size(direction)-1);
        major_step = 0.1;
        set(Slider_axes_2dview,'Min',1,'Max',Domain_size(direction),'SliderStep', [minor_step, major_step],'Value',1);
        % Update text
        set(Text_slider,'String',['Position: ' num2str(1) '/' num2str(Domain_size(direction))]);        
        % Update figure
        update_figure
    end

% Select savefolder
    function select_savefolder_Callback (~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select location where the result folder will be created';
        % Open dialog box to choose folder path
        resultfolder_location = uigetdir(matlabroot,str_dialogbox);
        if resultfolder_location==0
            % User clicked cancel button or closed the dialog box
            correct_folder=0;
            set(Text_Folder_statut,'String','Folder location: NOT DEFINED','ForegroundColor','r');        
        else
            % Update savefolder
            correct_folder=1;
            if ispc
                resultfolder_location = [resultfolder_location '\'];
            else
                resultfolder_location = [resultfolder_location '/'];
            end
            set(Text_Folder_statut,'String','Folder location: DEFINED','ForegroundColor',[34/255 177/255 76/255]);        
        end
        % Check save is possible
        save_possible =  correct_color*correct_ROI*correct_voxel_size*correct_folder;
        if save_possible==1
            set(Button_Oncefinisehd,'enable','on');
        else
            set(Button_Oncefinisehd,'enable','off');
        end        
    end

% Save and reset
    function Savevolume_and_options_Callback(~,~)
        % Increment number
        number_volume_to_analyse=number_volume_to_analyse+1;

        % Save:
        % Order of calculation
        inputvolume.volume(number_volume_to_analyse).orderofcalculation=number_volume_to_analyse;
        % loading path
        inputvolume.volume(number_volume_to_analyse).loadingpath=Button_add_volume.UserData.full_path_new_volume;
        % Input file name
        [~, name_, ~] = fileparts(Text_new_volume.String);
        inputvolume.volume(number_volume_to_analyse).filename_input=name_;
        % ROI
        inputvolume.volume(number_volume_to_analyse).ROI=ROI;
        % ROI and equivalent cubic domain size in micrometers
        new_domain_size=ROI(:,2)-ROI(:,1)+1;
        % Number of voxel
        ROI_number_of_voxel=prod(new_domain_size);
        ROIvoxelsize_number_of_voxel = floor(ROI_number_of_voxel*((initial_voxel_size/new_voxel_size)^3));
        inputvolume.volume(number_volume_to_analyse).ROI_number_of_voxel = ROI_number_of_voxel;
        inputvolume.volume(number_volume_to_analyse).ROI_theoritical_number_of_voxel_voxelresize = ROIvoxelsize_number_of_voxel;
        
        % Update domain size text
        size_direction_1 = new_domain_size(1)*initial_voxel_size/1000;
        size_direction_2 = new_domain_size(2)*initial_voxel_size/1000;
        size_direction_3 = new_domain_size(3)*initial_voxel_size/1000;
        size_cubeequivalent = (size_direction_1*size_direction_2*size_direction_3)^(1/3);
        inputvolume.volume(number_volume_to_analyse).size_direction_1=size_direction_1;   
        inputvolume.volume(number_volume_to_analyse).size_direction_2=size_direction_2;        
        inputvolume.volume(number_volume_to_analyse).size_direction_3=size_direction_3;        
        inputvolume.volume(number_volume_to_analyse).equivalentcubicsize=size_cubeequivalent;        
        % Initial voxel size
        inputvolume.volume(number_volume_to_analyse).initialvoxelsize=initial_voxel_size;
        % New voxel size
        inputvolume.volume(number_volume_to_analyse).newvoxelsize=new_voxel_size;
        % Result folder location
        inputvolume.volume(number_volume_to_analyse).resultfolderlocation=resultfolder_location;
        % Result folder
        inputvolume.volume(number_volume_to_analyse).resultfoldername=Edit_resultfoldername.String;
        % Direction name
        table_data = table_input_ROI.Data;
        direction_=table_data(:,5);
        inputvolume.volume(number_volume_to_analyse).directioname=direction_;
        % Phase name
        table_data = table_input_phase.Data;
        phase_=table_data(:,6);
        inputvolume.volume(number_volume_to_analyse).phasename=phase_;        
        % Phase color
        table_data = table_input_phase.Data;
        color_=table_data(:,3:5);
        color_ = cell2mat(color_);
        inputvolume.volume(number_volume_to_analyse).phasecolor=color_;  
        % Initial code phase
        code_=table_data(:,1);
        code_ = cell2mat(code_);        
        inputvolume.volume(number_volume_to_analyse).phasecode=code_;        
        % Assigned phase
        assignedcode_=table_data(:,7);
        assignedcode_ = cell2mat(assignedcode_);        
        inputvolume.volume(number_volume_to_analyse).assignedcode=assignedcode_;
                
        % Number of phase
        inputvolume.volume(number_volume_to_analyse).initial_number_phase=number_phase;  
        inputvolume.volume(number_volume_to_analyse).assignedcode_number_phase=length(unique(assignedcode_));  
       
        % Volume information
        table_data = table_input_info.Data;
        inputvolume.volume(number_volume_to_analyse).volumeinformation = table_data;

        % Display save confirmation
        set(Text_saveconfirmation,'Visible','on');
       
        % Update other tab!
        set(Text_nodata_tabloaded,'Visible','off');
        set(Text_nodata_runcalculation,'Visible','off');
        set(Text_tab_imported,'Visible','on');
        set(Text_instruction_tabloaded,'Visible','on');
        set(table_organize_volume,'Visible','on','enable','on');
        set(Button_changeorganise,'Visible','on','enable','on');
        set(Button_runcalculation,'Visible','on','enable','on');
        
        % Update table_organize_volume
        new_volume=1;
        Update_tab_importedvolume_Callback
        % Reset
        Reset_Callback
    end

% Update axes_2dview
    function update_figure
        Domain_size=Button_add_volume.UserData.Domain_size; % Get domain size
        number_phase=Button_add_volume.UserData.number_phase; % Get number of phase
        Phase_code=Button_add_volume.UserData.Phase_code; % Get phase code        
        new_volume=Button_add_volume.UserData.new_volume; % Get volume         
        % Colors
        for current_phase=1:1:number_phase
            RGB_phase.index(current_phase).rgb = [color_phase(current_phase,1) color_phase(current_phase,2) color_phase(current_phase,3)];
        end
        if direction==1
            % Initializaion
            slice_color = zeros(Domain_size(2),Domain_size(3),3); % RGB color map
            slice_r = zeros(Domain_size(2),Domain_size(3)); % Red color map
            slice_g = zeros(Domain_size(2),Domain_size(3)); % Green color map
            slice_b = zeros(Domain_size(2),Domain_size(3)); % Blue color map
            % Attribute RGB colors for each voxel
            for current_phase=1:1:number_phase
                code_tmp =Phase_code(current_phase); % Current phase code
                slice_r(new_volume(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                slice_g(new_volume(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                slice_b(new_volume(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
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
                slice_r(new_volume(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                slice_g(new_volume(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                slice_b(new_volume(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
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
                slice_r(new_volume(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                slice_g(new_volume(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                slice_b(new_volume(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
            end
        end
        slice_color(:,:,1)=slice_r; % Attribute RGB color
        slice_color(:,:,2)=slice_g;
        slice_color(:,:,3)=slice_b;
        % Display the slice
        slice_image = image(slice_color,'parent',axes_2dview);
        % Get axis position
        axe_position =axes_2dview.Position;
        set(axes_2dview ,'Position', axe_position);
        % Remove tick and label
        set(axes_2dview,'xtick',[],'ytick',[]);
        set(axes_2dview,'xticklabel',[],'yticklabel',[]);
        % Fit the axes box
        axis(axes_2dview,'tight');
        % Aspect ratio is 1:1
        axis(axes_2dview,'equal');
        % Align
        align([Button_add_volume axes_2dview Text_new_volume],'HorizontalAlignment','Center');
    end

% Reset
    function Reset_Callback(~,~)
        % Reset variables
        Position_slice=[1 1 1];
        initial_voxel_size=1;
        new_voxel_size=1;
        ROI=0;
        correct_ROI=1;
        correct_voxel_size=1;
        correct_folder=0;
        correct_color=1;
        resultfolder_location=0;
        
        % Reactivate add volume button
        set(Button_add_volume,'enable','on');
        
        % Reset ui
        set(Text_new_volume,'enable','off','String','No volume selected');
        set(Text_slider,'enable','off','Visible','off');
        set(Popup_slider,'enable','off','Visible','off');
        set(Slider_axes_2dview,'enable','off','Visible','off');
        set(Text_volumenumber,'enable','off','Visible','off');
        set(Text_voxelsize,'Visible','off');
        set(Text_initialvoxelsize,'Visible','off');
        set(Text_newvoxelsize,'Visible','off');
        set(Edit_initialvoxelsize,'Visible','off','enable','off');
        set(Edit_newvoxelsize,'Visible','off','enable','off');
        
        set(Text_InitialROI,'Visible','off');
        set(Text_Numberofvoxel_values,'Visible','off'); set(Text_Domainsize_values,'Visible','off');  set(Text_Equivalentdomainsize_values,'Visible','off');
        set(Text_Numberofvoxel,'Visible','off'); set(Text_Numberofvoxel_bis,'Visible','off');
        set(Text_Domainsize,'Visible','off'); set(Text_Domainsize_bis,'Visible','off');
        set(Text_Equivalentdomainsize,'Visible','off'); set(Text_Equivalentdomainsize_bis,'Visible','off');
        
        set(Text_userROI,'Visible','off');
        set(Text_Numberofvoxel_uservalues,'Visible','off');
        set(Text_Domainsize_uservalues,'Visible','off');
        set(Text_Equivalentdomainsize_uservalues,'Visible','off');
        
        set(Text_save_volume,'Visible','off');
        set(Button_select_savefolder,'Visible','off','enable','off');
        set(Edit_resultfoldername,'Visible','off','enable','off');
        
        set(Text_Voxelsize_statut,'Visible','off','String','Voxel size: CORRECT','ForegroundColor',[34/255 177/255 76/255]);
        set(Text_ROI_statut,'Visible','off','String','Region of interest (ROI): CORRECT','ForegroundColor',[34/255 177/255 76/255]);
        set(Text_Color_statut,'Visible','off','String','Colors: CORRECT','ForegroundColor',[34/255 177/255 76/255]);
        set(Text_Assigncode_statut,'Visible','off','String','Assigned code: CORRECT','ForegroundColor',[34/255 177/255 76/255]);
        set(Text_Folder_statut,'Visible','off','String','Folder location: NOT DEFINED','ForegroundColor','r');
        set(Button_Oncefinisehd,'Visible','off');
        set(Button_Reset,'Visible','off','enable','off');
        
        set(Text_directiontable,'visible','off');
        set(table_input_ROI,'visible','off');set(table_input_ROI,'enable','off');
        set(Text_phasetable,'visible','off');
        set(table_input_phase,'visible','off');set(table_input_phase,'enable','off');
        
        set(Text_volumeinfo,'visible','off');
        set(table_input_info,'visible','off');
        
        set(axes_2dview,'Visible','off');
        set(axes_2dview.Children,'Visible','off');
        
    end

%%
%% TAB IMPORTED VOLUME(S)
%%

Text_nodata_tabloaded = uicontrol('Parent', tab_imported,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','No data. Go to ''Import segmented volume'' tab','Visible','on',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0 0.4 1 0.2]);

% Description and instructions
Text_tab_imported = uicontrol('Parent', tab_imported, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],'Visible','off',...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Review volume(s) to analyze. Remove volume(s), change order of calculations if wished.');

Text_instruction_tabloaded = uicontrol('Parent', tab_imported,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','1) Check information is correct. You can only change the order of calculation (starting from 1) in this tab. If you want to remove a volume, enter a 0 value in the ''order of calculation'' associated cell.','Visible','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [description_tab_fromleft 0.8 description_tab_xlenght 0.1]);
Button_changeorganise = uicontrol('Parent', tab_imported,'Style', 'pushbutton','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','2) Once you have entered all your modifications, click to update the table.','Visible','off','enable','off',...
     'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.2 0.025 0.6 0.04],'Callback',{@Update_tab_importedvolume_Callback}); 

% Organise
table_organize_volume = uitable('Parent', tab_imported,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.05 0.1 0.9 0.73],'CellEditCallback',@celledit_organisevolume_Callback);
table_organize_volume.ColumnName = {'File name input','Equivalent cubic ROI size (um)','Chosen voxel size (nm)','Order of calculation'}; % Column name
table_organize_volume.ColumnEditable = [false false false true]; % Select column editable
table_organize_volume.ColumnWidth = {500, 'auto', 'auto', 'auto'}; % Auto width
table_organize_volume.RowName = []; % Remove row name
    
A{1,1} = 'Incorrect order of calculation';                                                           
A{1,2} = ' - Values must be integer >=0';
A{1,3} = ' - No duplicate values (except 0)';
A{1,4} = ' - Values must be continous with no gaps (i.e. 1,2,3,...)';
multiple_line_string = sprintf('%s\n%s\n%s\n%s\n%s\n%s',A{1,1},A{1,2},A{1,3},A{1,4}); 
Text_wrongorderofcalculation = uicontrol('Parent', tab_imported,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String',multiple_line_string,'Visible','off',...
     'ForegroundColor','r','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.15 0.1 0.75 0.6]);  
 
%%
%% TAB VOLUME FRACTIONS
%%

%% GUI
Text_tab_volumefraction = uicontrol('Parent', tab_volume_fraction, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the volume fractions.'); 

checkbox_volumefraction_todo = uicontrol('Parent', tab_volume_fraction, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 3.5*pos_delta_x 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate volume fractions.','Value',PROPERTY.volumefractions.todo,'Callback',{@checkbox_volumefraction_todo_Callback}); 
annotation(tab_volume_fraction,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_volumefraction_todo_voxelsize = uicontrol('Parent', tab_volume_fraction, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.volumefractions.voxel_size_dependence.todo,'Callback',{@checkbox_volumefraction_todo_voxelsize_Callback}); 

table_options_volumefractions_voxelsize = uitable('Parent', tab_volume_fraction,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tablevolumefractionvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tablevolumefractionvoxelsize_Callback);
table_options_volumefractions_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_volumefractions_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_volumefractions_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_volumefraction_todo_voxelsize.Position(3)}; % Auto width
table_options_volumefractions_voxelsize.RowName = []; % Remove row name
table_options_volumefractions_voxelsize.Data=PROPERTY.volumefractions.voxel_size_dependence.voxel';

Button_options_volumefractions_voxelsize_add = uicontrol('Parent', tab_volume_fraction, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_volumefractions_voxelsize_add_Callback});
Button_options_volumefractions_voxelsize_remove = uicontrol('Parent', tab_volume_fraction, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_volumefractions_voxelsize_remove_Callback});
Button_options_volumefractions_voxelsize_update = uicontrol('Parent', tab_volume_fraction, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_volumefractions_voxelsize_update_Callback});

text_volumefraction_todo_RVE = uicontrol('Parent', tab_volume_fraction, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_volumefraction = uicontrol('Parent', tab_volume_fraction,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_volumefractions_RVE_new = uicontrol('Parent', tab_volume_fraction, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_volumefractions_RVE_add_Callback});
Button_options_volumefractions_RVE_reset = uicontrol('Parent', tab_volume_fraction, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_volumefractions_RVE_reset_Callback});
       
table_options_volumefractions_RVE_subvolumes = uitable('Parent', tab_volume_fraction,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_volumefractions_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_volumefractions_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_volumefractions_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.volumefractions.RVE(1).type};
column2={num2str(PROPERTY.volumefractions.RVE(1).divisions)};
column3={PROPERTY.volumefractions.RVE(1).subs2};
column4={PROPERTY.volumefractions.RVE(1).subs4};
column5={PROPERTY.volumefractions.RVE(1).Aspectratio};
column6={PROPERTY.volumefractions.RVE(1).Constantdirection};
column7={PROPERTY.volumefractions.RVE(1).Growthdirection};
column8={[PROPERTY.volumefractions.RVE(1).Growthperstep ' ' PROPERTY.volumefractions.RVE(1).Growthrelativeto]};
column9={[PROPERTY.volumefractions.RVE(1).firstuniquevolume_size ' ' PROPERTY.volumefractions.RVE(1).firstuniquevolume_unit]};
table_options_volumefractions_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];

%% CALLBACKS

    function pushbutton_options_volumefractions_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_volumefraction.String(Popup_RVE_volumefraction.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'volumefractions', table_options_volumefractions_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'volumefractions', table_options_volumefractions_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'volumefractions', table_options_volumefractions_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'volumefractions', table_options_volumefractions_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'volumefractions', table_options_volumefractions_RVE_subvolumes);
        end
    end

    function pushbutton_options_volumefractions_RVE_reset_Callback(~,~)
        PROPERTY.volumefractions.RVE=[];
        PROPERTY.volumefractions.number_RVE = 0; % no RVE to perform
        table_options_volumefractions_RVE_subvolumes.Data=[];
    end

    function checkbox_volumefraction_todo_Callback(~,~)
        PROPERTY.volumefractions.todo = checkbox_volumefraction_todo.Value;
        if PROPERTY.volumefractions.todo
            set(checkbox_volumefraction_todo_voxelsize,'enable','on');
            checkbox_volumefraction_todo_voxelsize_Callback
            set(text_volumefraction_todo_RVE,'enable','on');            
            set(Popup_RVE_volumefraction,'enable','on');            
            set(Button_options_volumefractions_RVE_new,'enable','on');            
            set(Button_options_volumefractions_RVE_reset,'enable','on');            
            set(table_options_volumefractions_RVE_subvolumes,'enable','on');            
        else
            set(checkbox_volumefraction_todo_voxelsize,'enable','off');
            set(table_options_volumefractions_voxelsize,'enable','off');
            set(Button_options_volumefractions_voxelsize_add,'enable','off');
            set(Button_options_volumefractions_voxelsize_remove,'enable','off');            
            set(Button_options_volumefractions_voxelsize_update,'enable','off');            
            set(text_volumefraction_todo_RVE,'enable','off');            
            set(Popup_RVE_volumefraction,'enable','off');            
            set(Button_options_volumefractions_RVE_new,'enable','off');            
            set(Button_options_volumefractions_RVE_reset,'enable','off');            
            set(table_options_volumefractions_RVE_subvolumes,'enable','off');                 
        end
    end

    function checkbox_volumefraction_todo_voxelsize_Callback(~,~)
        PROPERTY.volumefractions.voxel_size_dependence.todo = checkbox_volumefraction_todo_voxelsize.Value;
        if PROPERTY.volumefractions.voxel_size_dependence.todo
            set(table_options_volumefractions_voxelsize,'enable','on');
            set(Button_options_volumefractions_voxelsize_add,'enable','on');
            set(Button_options_volumefractions_voxelsize_remove,'enable','on');
            set(Button_options_volumefractions_voxelsize_update,'enable','on');
        else
            set(table_options_volumefractions_voxelsize,'enable','off');
            set(Button_options_volumefractions_voxelsize_add,'enable','off');
            set(Button_options_volumefractions_voxelsize_remove,'enable','off');
            set(Button_options_volumefractions_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_volumefractions_voxelsize_add_Callback(~,~)
        val_ = table_options_volumefractions_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_volumefractions_voxelsize.Data = val_;
        PROPERTY.volumefractions.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_volumefractions_voxelsize_remove_Callback(~,~)
        val_ = table_options_volumefractions_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_volumefractions_voxelsize.Data = val_;
        PROPERTY.volumefractions.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_volumefractions_voxelsize_update_Callback(~,~)
        val_ = table_options_volumefractions_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_volumefractions_voxelsize.Data = val_;
        PROPERTY.volumefractions.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tablevolumefractionvoxelsize_Callback(~,~)
        val_ = table_options_volumefractions_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_volumefraction_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.volumefractions.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_volumefraction_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end

%%
%% TAB SPECIFIC SURFACE AREA (DIRECT)
%%

%% GUI
Text_tab_specificsurfacearea = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the specific surface area (direct counting).'); 

checkbox_specificsurfacearea_direct_todo = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate specific surface area with the direct counting method.','Value',PROPERTY.specificsurfacearea_directmethod.todo,'Callback',{@checkbox_specificsurfacearea_direct_todo_Callback}); 
annotation(tab_specificsurfacearea,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_specificsurfacearea_direct_todo_voxelsize = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.todo,'Callback',{@checkbox_specificsurfacearea_direct_todo_voxelsize_Callback}); 

table_options_specificsurfacearea_direct_voxelsize = uitable('Parent', tab_specificsurfacearea,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tablespecificsurfacearea_directvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tablespecificsurfacearea_directvoxelsize_Callback);
table_options_specificsurfacearea_direct_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_specificsurfacearea_direct_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_specificsurfacearea_direct_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_specificsurfacearea_direct_todo_voxelsize.Position(3)}; % Auto width
table_options_specificsurfacearea_direct_voxelsize.RowName = []; % Remove row name
table_options_specificsurfacearea_direct_voxelsize.Data=PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel';

Button_options_specificsurfacearea_direct_voxelsize_add = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_sp_direct_voxelsize_add_Callback});
Button_options_specificsurfacearea_direct_voxelsize_remove = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_sp_direct_voxelsize_remove_Callback});
Button_options_specificsurfacearea_direct_voxelsize_update = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_sp_direct_voxelsize_update_Callback});

text_specificsurfacearea_direct_todo_RVE = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_specificsurfacearea_direct = uicontrol('Parent', tab_specificsurfacearea,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_specificsurfacearea_direct_RVE_new = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_sp_direct_RVE_add_Callback});
Button_options_specificsurfacearea_direct_RVE_reset = uicontrol('Parent', tab_specificsurfacearea, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_sp_direct_RVE_reset_Callback});
       
table_options_specificsurfacearea_direct_RVE_subvolumes = uitable('Parent', tab_specificsurfacearea,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_specificsurfacearea_direct_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_specificsurfacearea_direct_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_specificsurfacearea_direct_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.specificsurfacearea_directmethod.RVE(1).type};
column2={num2str(PROPERTY.specificsurfacearea_directmethod.RVE(1).divisions)};
column3={PROPERTY.specificsurfacearea_directmethod.RVE(1).subs2};
column4={PROPERTY.specificsurfacearea_directmethod.RVE(1).subs4};
column5={PROPERTY.specificsurfacearea_directmethod.RVE(1).Aspectratio};
column6={PROPERTY.specificsurfacearea_directmethod.RVE(1).Constantdirection};
column7={PROPERTY.specificsurfacearea_directmethod.RVE(1).Growthdirection};
column8={[PROPERTY.specificsurfacearea_directmethod.RVE(1).Growthperstep ' ' PROPERTY.specificsurfacearea_directmethod.RVE(1).Growthrelativeto]};
column9={[PROPERTY.specificsurfacearea_directmethod.RVE(1).firstuniquevolume_size ' ' PROPERTY.specificsurfacearea_directmethod.RVE(1).firstuniquevolume_unit]};
table_options_specificsurfacearea_direct_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];


% Options
annotation(tab_specificsurfacearea,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.5 0.5],'visible','on');

Text_sp_definition = uicontrol('Parent', tab_specificsurfacearea,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Definition of specific surface area',...
     'HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.425 2.5*pos_delta_x 0.3*pos_delta_y]);
 
popup_sp_definition = uicontrol('Parent', tab_specificsurfacearea,'Style', 'popup','String',specificsurfacearea_choice,'Value',1,'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Visible','on','enable','on',...
    'Units','normalized','Position', [pos_x_start 0.425-0.5*pos_delta_y 2.5*pos_delta_x 0.5*pos_delta_y],'Callback', @popup_sp_method_Callback);  
    function popup_sp_method_Callback(~,~)
        PROPERTY.specificsurfacearea_directmethod.definition=char(specificsurfacearea_choice(popup_sp_definition.Value));
    end

Text_spdirect_correctivefactor = uicontrol('Parent', tab_specificsurfacearea,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Corrective factor (cubic representation overestimation)',...
     'HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.35 3.5*pos_delta_x 0.3*pos_delta_y]);
Edit_spdirect_correctivefactor = uicontrol('Parent', tab_specificsurfacearea,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',default_correctivefactor_spdirect,...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos_x_start 0.35-0.3*pos_delta_y 3.5*pos_delta_x 0.3*pos_delta_y],'Callback', @edit_spdirect_correctivefactor_Callback); 
    function edit_spdirect_correctivefactor_Callback(~,~)
        PROPERTY.specificsurfacearea_directmethod.correctivefactor=str2double(Edit_spdirect_correctivefactor.String);
    end



%% CALLBACKS

    function pushbutton_options_sp_direct_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_specificsurfacearea_direct.String(Popup_RVE_specificsurfacearea_direct.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'specificsurfacearea_directmethod', table_options_specificsurfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'specificsurfacearea_directmethod', table_options_specificsurfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'specificsurfacearea_directmethod', table_options_specificsurfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'specificsurfacearea_directmethod', table_options_specificsurfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'specificsurfacearea_directmethod', table_options_specificsurfacearea_direct_RVE_subvolumes);
        end
    end

    function pushbutton_options_sp_direct_RVE_reset_Callback(~,~)
        PROPERTY.specificsurfacearea_directmethod.RVE=[];
        PROPERTY.specificsurfacearea_directmethod.number_RVE = 0; % no RVE to perform
        table_options_specificsurfacearea_direct_RVE_subvolumes.Data=[];
    end

    function checkbox_specificsurfacearea_direct_todo_Callback(~,~)
        PROPERTY.specificsurfacearea_directmethod.todo = checkbox_specificsurfacearea_direct_todo.Value;
        if PROPERTY.specificsurfacearea_directmethod.todo
            set(checkbox_specificsurfacearea_direct_todo_voxelsize,'enable','on');
            checkbox_specificsurfacearea_direct_todo_voxelsize_Callback
            set(text_specificsurfacearea_direct_todo_RVE,'enable','on');            
            set(Popup_RVE_specificsurfacearea_direct,'enable','on');            
            set(Button_options_specificsurfacearea_direct_RVE_new,'enable','on');            
            set(Button_options_specificsurfacearea_direct_RVE_reset,'enable','on');            
            set(table_options_specificsurfacearea_direct_RVE_subvolumes,'enable','on'); 
            set(Text_sp_definition,'enable','on'); 
            set(popup_sp_definition,'enable','on'); 
            set(Text_spdirect_correctivefactor,'enable','on'); 
            set(Edit_spdirect_correctivefactor,'enable','on'); 
            
        else
            set(checkbox_specificsurfacearea_direct_todo_voxelsize,'enable','off');
            set(table_options_specificsurfacearea_direct_voxelsize,'enable','off');
            set(Button_options_specificsurfacearea_direct_voxelsize_add,'enable','off');
            set(Button_options_specificsurfacearea_direct_voxelsize_remove,'enable','off');            
            set(Button_options_specificsurfacearea_direct_voxelsize_update,'enable','off');            
            set(text_specificsurfacearea_direct_todo_RVE,'enable','off');            
            set(Popup_RVE_specificsurfacearea_direct,'enable','off');            
            set(Button_options_specificsurfacearea_direct_RVE_new,'enable','off');            
            set(Button_options_specificsurfacearea_direct_RVE_reset,'enable','off');            
            set(table_options_specificsurfacearea_direct_RVE_subvolumes,'enable','off');    
            set(Text_sp_definition,'enable','off'); 
            set(popup_sp_definition,'enable','off'); 
            set(Text_spdirect_correctivefactor,'enable','off'); 
            set(Edit_spdirect_correctivefactor,'enable','off');             
        end
    end

    function checkbox_specificsurfacearea_direct_todo_voxelsize_Callback(~,~)
        PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.todo = checkbox_specificsurfacearea_direct_todo_voxelsize.Value;
        if PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.todo
            set(table_options_specificsurfacearea_direct_voxelsize,'enable','on');
            set(Button_options_specificsurfacearea_direct_voxelsize_add,'enable','on');
            set(Button_options_specificsurfacearea_direct_voxelsize_remove,'enable','on');
            set(Button_options_specificsurfacearea_direct_voxelsize_update,'enable','on');
        else
            set(table_options_specificsurfacearea_direct_voxelsize,'enable','off');
            set(Button_options_specificsurfacearea_direct_voxelsize_add,'enable','off');
            set(Button_options_specificsurfacearea_direct_voxelsize_remove,'enable','off');
            set(Button_options_specificsurfacearea_direct_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_sp_direct_voxelsize_add_Callback(~,~)
        val_ = table_options_specificsurfacearea_direct_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_specificsurfacearea_direct_voxelsize.Data = val_;
        PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_sp_direct_voxelsize_remove_Callback(~,~)
        val_ = table_options_specificsurfacearea_direct_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_specificsurfacearea_direct_voxelsize.Data = val_;
        PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_sp_direct_voxelsize_update_Callback(~,~)
        val_ = table_options_specificsurfacearea_direct_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_specificsurfacearea_direct_voxelsize.Data = val_;
        PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tablespecificsurfacearea_directvoxelsize_Callback(~,~)
        val_ = table_options_specificsurfacearea_direct_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_specificsurfacearea_direct_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.specificsurfacearea_directmethod.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_specificsurfacearea_direct_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end

%%
%% TAB SPECIFIC INTERFACE AREA (DIRECT)
%%

%% GUI
Text_tab_specificinterfacearea = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the specific interface area (direct counting).'); 

checkbox_specificinterfacearea_direct_todo = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate specific surface area with the direct counting method.','Value',PROPERTY.specificinterfacearea_directmethod.todo,'Callback',{@checkbox_specificinterfacearea_direct_todo_Callback}); 
annotation(tab_specificinterfacearea,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_specificinterfacearea_direct_todo_voxelsize = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.todo,'Callback',{@checkbox_specificinterfacearea_direct_todo_voxelsize_Callback}); 

table_options_specificinterfacearea_direct_voxelsize = uitable('Parent', tab_specificinterfacearea,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tablespecificinterfacearea_directvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tablespecificinterfacearea_directvoxelsize_Callback);
table_options_specificinterfacearea_direct_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_specificinterfacearea_direct_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_specificinterfacearea_direct_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_specificinterfacearea_direct_todo_voxelsize.Position(3)}; % Auto width
table_options_specificinterfacearea_direct_voxelsize.RowName = []; % Remove row name
table_options_specificinterfacearea_direct_voxelsize.Data=PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.voxel';

Button_options_specificinterfacearea_direct_voxelsize_add = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_spi_direct_voxelsize_add_Callback});
Button_options_specificinterfacearea_direct_voxelsize_remove = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_spi_direct_voxelsize_remove_Callback});
Button_options_specificinterfacearea_direct_voxelsize_update = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_spi_direct_voxelsize_update_Callback});

text_specificinterfacearea_direct_todo_RVE = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_specificinterfacearea_direct = uicontrol('Parent', tab_specificinterfacearea,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_specificinterfacearea_direct_RVE_new = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_spi_direct_RVE_add_Callback});
Button_options_specificinterfacearea_direct_RVE_reset = uicontrol('Parent', tab_specificinterfacearea, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_spi_direct_RVE_reset_Callback});
       
table_options_specificinterfacearea_direct_RVE_subvolumes = uitable('Parent', tab_specificinterfacearea,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_specificinterfacearea_direct_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_specificinterfacearea_direct_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_specificinterfacearea_direct_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.specificinterfacearea_directmethod.RVE(1).type};
column2={num2str(PROPERTY.specificinterfacearea_directmethod.RVE(1).divisions)};
column3={PROPERTY.specificinterfacearea_directmethod.RVE(1).subs2};
column4={PROPERTY.specificinterfacearea_directmethod.RVE(1).subs4};
column5={PROPERTY.specificinterfacearea_directmethod.RVE(1).Aspectratio};
column6={PROPERTY.specificinterfacearea_directmethod.RVE(1).Constantdirection};
column7={PROPERTY.specificinterfacearea_directmethod.RVE(1).Growthdirection};
column8={[PROPERTY.specificinterfacearea_directmethod.RVE(1).Growthperstep ' ' PROPERTY.specificinterfacearea_directmethod.RVE(1).Growthrelativeto]};
column9={[PROPERTY.specificinterfacearea_directmethod.RVE(1).firstuniquevolume_size ' ' PROPERTY.specificinterfacearea_directmethod.RVE(1).firstuniquevolume_unit]};
table_options_specificinterfacearea_direct_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];


% Options
annotation(tab_specificinterfacearea,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.5 0.5],'visible','on');

Text_spidirect_correctivefactor = uicontrol('Parent', tab_specificinterfacearea,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Corrective factor (cubic representation overestimation)',...
     'HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.425 3.5*pos_delta_x 0.3*pos_delta_y]);
Edit_spidirect_correctivefactor = uicontrol('Parent', tab_specificinterfacearea,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',default_correctivefactor_spdirect,...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [pos_x_start 0.425-0.3*pos_delta_y 3.5*pos_delta_x 0.3*pos_delta_y],'Callback', @edit_spdirect_correctivefactor_Callback); 
    function edit_spidirect_correctivefactor_Callback(~,~)
        PROPERTY.specificinterfacearea_directmethod.correctivefactor=str2double(Edit_spidirect_correctivefactor.String);
        PROPERTY.specificinterfacearea_directmethod.correctivefactor
    end



%% CALLBACKS

    function pushbutton_options_spi_direct_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_specificinterfacearea_direct.String(Popup_RVE_specificinterfacearea_direct.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'specificinterfacearea_directmethod', table_options_specificinterfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'specificinterfacearea_directmethod', table_options_specificinterfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'specificinterfacearea_directmethod', table_options_specificinterfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'specificinterfacearea_directmethod', table_options_specificinterfacearea_direct_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'specificinterfacearea_directmethod', table_options_specificinterfacearea_direct_RVE_subvolumes);
        end
    end

    function pushbutton_options_spi_direct_RVE_reset_Callback(~,~)
        PROPERTY.specificinterfacearea_directmethod.RVE=[];
        PROPERTY.specificinterfacearea_directmethod.number_RVE = 0; % no RVE to perform
        table_options_specificinterfacearea_direct_RVE_subvolumes.Data=[];
    end

    function checkbox_specificinterfacearea_direct_todo_Callback(~,~)
        PROPERTY.specificinterfacearea_directmethod.todo = checkbox_specificinterfacearea_direct_todo.Value;
        if PROPERTY.specificinterfacearea_directmethod.todo
            set(checkbox_specificinterfacearea_direct_todo_voxelsize,'enable','on');
            checkbox_specificinterfacearea_direct_todo_voxelsize_Callback
            set(text_specificinterfacearea_direct_todo_RVE,'enable','on');            
            set(Popup_RVE_specificinterfacearea_direct,'enable','on');            
            set(Button_options_specificinterfacearea_direct_RVE_new,'enable','on');            
            set(Button_options_specificinterfacearea_direct_RVE_reset,'enable','on');            
            set(table_options_specificinterfacearea_direct_RVE_subvolumes,'enable','on'); 
            set(Text_spidirect_correctivefactor,'enable','on'); 
            set(Edit_spidirect_correctivefactor,'enable','on'); 
            
        else
            set(checkbox_specificinterfacearea_direct_todo_voxelsize,'enable','off');
            set(table_options_specificinterfacearea_direct_voxelsize,'enable','off');
            set(Button_options_specificinterfacearea_direct_voxelsize_add,'enable','off');
            set(Button_options_specificinterfacearea_direct_voxelsize_remove,'enable','off');            
            set(Button_options_specificinterfacearea_direct_voxelsize_update,'enable','off');            
            set(text_specificinterfacearea_direct_todo_RVE,'enable','off');            
            set(Popup_RVE_specificinterfacearea_direct,'enable','off');            
            set(Button_options_specificinterfacearea_direct_RVE_new,'enable','off');            
            set(Button_options_specificinterfacearea_direct_RVE_reset,'enable','off');            
            set(table_options_specificinterfacearea_direct_RVE_subvolumes,'enable','off');    
            set(Text_spidirect_correctivefactor,'enable','off'); 
            set(Edit_spidirect_correctivefactor,'enable','off');             
        end
    end

    function checkbox_specificinterfacearea_direct_todo_voxelsize_Callback(~,~)
        PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.todo = checkbox_specificinterfacearea_direct_todo_voxelsize.Value;
        if PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.todo
            set(table_options_specificinterfacearea_direct_voxelsize,'enable','on');
            set(Button_options_specificinterfacearea_direct_voxelsize_add,'enable','on');
            set(Button_options_specificinterfacearea_direct_voxelsize_remove,'enable','on');
            set(Button_options_specificinterfacearea_direct_voxelsize_update,'enable','on');
        else
            set(table_options_specificinterfacearea_direct_voxelsize,'enable','off');
            set(Button_options_specificinterfacearea_direct_voxelsize_add,'enable','off');
            set(Button_options_specificinterfacearea_direct_voxelsize_remove,'enable','off');
            set(Button_options_specificinterfacearea_direct_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_spi_direct_voxelsize_add_Callback(~,~)
        val_ = table_options_specificinterfacearea_direct_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_specificinterfacearea_direct_voxelsize.Data = val_;
        PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_spi_direct_voxelsize_remove_Callback(~,~)
        val_ = table_options_specificinterfacearea_direct_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_specificinterfacearea_direct_voxelsize.Data = val_;
        PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_spi_direct_voxelsize_update_Callback(~,~)
        val_ = table_options_specificinterfacearea_direct_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_specificinterfacearea_direct_voxelsize.Data = val_;
        PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tablespecificinterfacearea_directvoxelsize_Callback(~,~)
        val_ = table_options_specificinterfacearea_direct_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_specificinterfacearea_direct_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.specificinterfacearea_directmethod.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_specificinterfacearea_direct_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end


%%
%% TAB PARTICLE SIZE (C-PSD)
%%


%% GUI
Text_tab_particlesize_cpsd = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the particle size (Continuum Particle-Size Distribution, C-PSD).'); 

checkbox_particlesizeCpsd_todo = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate particle size with C-PSD method.','Value',PROPERTY.particlesize_cpsd.todo,'Callback',{@checkbox_particlesizeCpsd_todo_Callback}); 
annotation(tab_particlesize_CPSD,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_particlesizeCpsd_todo_voxelsize = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.particlesize_cpsd.voxel_size_dependence.todo,'Callback',{@checkbox_particlesizeCpsd_todo_voxelsize_Callback}); 

table_options_particlesizeCpsd_voxelsize = uitable('Parent', tab_particlesize_CPSD,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tableparticlesizeCpsdvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tableparticlesizeCpsdvoxelsize_Callback);
table_options_particlesizeCpsd_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_particlesizeCpsd_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_particlesizeCpsd_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_particlesizeCpsd_todo_voxelsize.Position(3)}; % Auto width
table_options_particlesizeCpsd_voxelsize.RowName = []; % Remove row name
table_options_particlesizeCpsd_voxelsize.Data=PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel';

Button_options_particlesizeCpsd_voxelsize_add = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_particlesizeCpsd_voxelsize_add_Callback});
Button_options_particlesizeCpsd_voxelsize_remove = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_particlesizeCpsd_voxelsize_remove_Callback});
Button_options_particlesizeCpsd_voxelsize_update = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_particlesizeCpsd_voxelsize_update_Callback});

text_particlesizeCpsd_todo_RVE = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_particlesizeCpsd = uicontrol('Parent', tab_particlesize_CPSD,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_particlesizeCpsd_RVE_new = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_particlesizeCpsd_RVE_add_Callback});
Button_options_particlesizeCpsd_RVE_reset = uicontrol('Parent', tab_particlesize_CPSD, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_particlesizeCpsd_RVE_reset_Callback});
       
table_options_particlesizeCpsd_RVE_subvolumes = uitable('Parent', tab_particlesize_CPSD,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_particlesizeCpsd_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_particlesizeCpsd_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_particlesizeCpsd_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.particlesize_cpsd.RVE(1).type};
column2={num2str(PROPERTY.particlesize_cpsd.RVE(1).divisions)};
column3={PROPERTY.particlesize_cpsd.RVE(1).subs2};
column4={PROPERTY.particlesize_cpsd.RVE(1).subs4};
column5={PROPERTY.particlesize_cpsd.RVE(1).Aspectratio};
column6={PROPERTY.particlesize_cpsd.RVE(1).Constantdirection};
column7={PROPERTY.particlesize_cpsd.RVE(1).Growthdirection};
column8={[PROPERTY.particlesize_cpsd.RVE(1).Growthperstep ' ' PROPERTY.particlesize_cpsd.RVE(1).Growthrelativeto]};
column9={[PROPERTY.particlesize_cpsd.RVE(1).firstuniquevolume_size ' ' PROPERTY.particlesize_cpsd.RVE(1).firstuniquevolume_unit]};
table_options_particlesizeCpsd_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];

%% CALLBACKS

    function pushbutton_options_particlesizeCpsd_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_particlesizeCpsd.String(Popup_RVE_particlesizeCpsd.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'particlesize_cpsd', table_options_particlesizeCpsd_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'particlesize_cpsd', table_options_particlesizeCpsd_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'particlesize_cpsd', table_options_particlesizeCpsd_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'particlesize_cpsd', table_options_particlesizeCpsd_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'particlesize_cpsd', table_options_particlesizeCpsd_RVE_subvolumes);
        end
    end

    function pushbutton_options_particlesizeCpsd_RVE_reset_Callback(~,~)
        PROPERTY.particlesize_cpsd.RVE=[];
        PROPERTY.particlesize_cpsd.number_RVE = 0; % no RVE to perform
        table_options_particlesizeCpsd_RVE_subvolumes.Data=[];
    end

    function checkbox_particlesizeCpsd_todo_Callback(~,~)
        PROPERTY.particlesize_cpsd.todo = checkbox_particlesizeCpsd_todo.Value;
        if PROPERTY.particlesize_cpsd.todo
            set(checkbox_particlesizeCpsd_todo_voxelsize,'enable','on');
            checkbox_particlesizeCpsd_todo_voxelsize_Callback
            set(text_particlesizeCpsd_todo_RVE,'enable','on');            
            set(Popup_RVE_particlesizeCpsd,'enable','on');            
            set(Button_options_particlesizeCpsd_RVE_new,'enable','on');            
            set(Button_options_particlesizeCpsd_RVE_reset,'enable','on');            
            set(table_options_particlesizeCpsd_RVE_subvolumes,'enable','on');            
        else
            set(checkbox_particlesizeCpsd_todo_voxelsize,'enable','off');
            set(table_options_particlesizeCpsd_voxelsize,'enable','off');
            set(Button_options_particlesizeCpsd_voxelsize_add,'enable','off');
            set(Button_options_particlesizeCpsd_voxelsize_remove,'enable','off');            
            set(Button_options_particlesizeCpsd_voxelsize_update,'enable','off');            
            set(text_particlesizeCpsd_todo_RVE,'enable','off');            
            set(Popup_RVE_particlesizeCpsd,'enable','off');            
            set(Button_options_particlesizeCpsd_RVE_new,'enable','off');            
            set(Button_options_particlesizeCpsd_RVE_reset,'enable','off');            
            set(table_options_particlesizeCpsd_RVE_subvolumes,'enable','off');                 
        end
    end

    function checkbox_particlesizeCpsd_todo_voxelsize_Callback(~,~)
        PROPERTY.particlesize_cpsd.voxel_size_dependence.todo = checkbox_particlesizeCpsd_todo_voxelsize.Value;
        if PROPERTY.particlesize_cpsd.voxel_size_dependence.todo
            set(table_options_particlesizeCpsd_voxelsize,'enable','on');
            set(Button_options_particlesizeCpsd_voxelsize_add,'enable','on');
            set(Button_options_particlesizeCpsd_voxelsize_remove,'enable','on');
            set(Button_options_particlesizeCpsd_voxelsize_update,'enable','on');
        else
            set(table_options_particlesizeCpsd_voxelsize,'enable','off');
            set(Button_options_particlesizeCpsd_voxelsize_add,'enable','off');
            set(Button_options_particlesizeCpsd_voxelsize_remove,'enable','off');
            set(Button_options_particlesizeCpsd_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_particlesizeCpsd_voxelsize_add_Callback(~,~)
        val_ = table_options_particlesizeCpsd_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_particlesizeCpsd_voxelsize.Data = val_;
        PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_particlesizeCpsd_voxelsize_remove_Callback(~,~)
        val_ = table_options_particlesizeCpsd_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_particlesizeCpsd_voxelsize.Data = val_;
        PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_particlesizeCpsd_voxelsize_update_Callback(~,~)
        val_ = table_options_particlesizeCpsd_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_particlesizeCpsd_voxelsize.Data = val_;
        PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tableparticlesizeCpsdvoxelsize_Callback(~,~)
        val_ = table_options_particlesizeCpsd_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_particlesizeCpsd_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.particlesize_cpsd.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_particlesizeCpsd_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end



%%
%% TAB PARTICLE SIZE (DISTANCE MAP)
%%


%% GUI
Text_tab_particlesize_dmap = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the particle size (Euclidean Distance Map Fitting EDMF).'); 

checkbox_particlesizeDmap_todo = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate particle size with distance map cumulative function fitting method.','Value',PROPERTY.particlesize_dmap.todo,'Callback',{@checkbox_particlesizeDmap_todo_Callback}); 
annotation(tab_particlesize_dmap,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_particlesizeDmap_todo_voxelsize = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.particlesize_dmap.voxel_size_dependence.todo,'Callback',{@checkbox_particlesizeDmap_todo_voxelsize_Callback}); 

table_options_particlesizeDmap_voxelsize = uitable('Parent', tab_particlesize_dmap,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tableparticlesizeDmapvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tableparticlesizeDmapvoxelsize_Callback);
table_options_particlesizeDmap_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_particlesizeDmap_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_particlesizeDmap_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_particlesizeDmap_todo_voxelsize.Position(3)}; % Auto width
table_options_particlesizeDmap_voxelsize.RowName = []; % Remove row name
table_options_particlesizeDmap_voxelsize.Data=PROPERTY.particlesize_dmap.voxel_size_dependence.voxel';

Button_options_particlesizeDmap_voxelsize_add = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_particlesizeDmap_voxelsize_add_Callback});
Button_options_particlesizeDmap_voxelsize_remove = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_particlesizeDmap_voxelsize_remove_Callback});
Button_options_particlesizeDmap_voxelsize_update = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_particlesizeDmap_voxelsize_update_Callback});

text_particlesizeDmap_todo_RVE = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_particlesizeDmap = uicontrol('Parent', tab_particlesize_dmap,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_particlesizeDmap_RVE_new = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_particlesizeDmap_RVE_add_Callback});
Button_options_particlesizeDmap_RVE_reset = uicontrol('Parent', tab_particlesize_dmap, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_particlesizeDmap_RVE_reset_Callback});
       
table_options_particlesizeDmap_RVE_subvolumes = uitable('Parent', tab_particlesize_dmap,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_particlesizeDmap_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_particlesizeDmap_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_particlesizeDmap_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.particlesize_dmap.RVE(1).type};
column2={num2str(PROPERTY.particlesize_dmap.RVE(1).divisions)};
column3={PROPERTY.particlesize_dmap.RVE(1).subs2};
column4={PROPERTY.particlesize_dmap.RVE(1).subs4};
column5={PROPERTY.particlesize_dmap.RVE(1).Aspectratio};
column6={PROPERTY.particlesize_dmap.RVE(1).Constantdirection};
column7={PROPERTY.particlesize_dmap.RVE(1).Growthdirection};
column8={[PROPERTY.particlesize_dmap.RVE(1).Growthperstep ' ' PROPERTY.particlesize_dmap.RVE(1).Growthrelativeto]};
column9={[PROPERTY.particlesize_dmap.RVE(1).firstuniquevolume_size ' ' PROPERTY.particlesize_dmap.RVE(1).firstuniquevolume_unit]};
table_options_particlesizeDmap_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];

%% CALLBACKS

    function pushbutton_options_particlesizeDmap_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_particlesizeDmap.String(Popup_RVE_particlesizeDmap.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'particlesize_dmap', table_options_particlesizeDmap_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'particlesize_dmap', table_options_particlesizeDmap_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'particlesize_dmap', table_options_particlesizeDmap_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'particlesize_dmap', table_options_particlesizeDmap_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'particlesize_dmap', table_options_particlesizeDmap_RVE_subvolumes);
        end
    end

    function pushbutton_options_particlesizeDmap_RVE_reset_Callback(~,~)
        PROPERTY.particlesize_dmap.RVE=[];
        PROPERTY.particlesize_dmap.number_RVE = 0; % no RVE to perform
        table_options_particlesizeDmap_RVE_subvolumes.Data=[];
    end

    function checkbox_particlesizeDmap_todo_Callback(~,~)
        PROPERTY.particlesize_dmap.todo = checkbox_particlesizeDmap_todo.Value;
        if PROPERTY.particlesize_dmap.todo
            set(checkbox_particlesizeDmap_todo_voxelsize,'enable','on');
            checkbox_particlesizeDmap_todo_voxelsize_Callback
            set(text_particlesizeDmap_todo_RVE,'enable','on');            
            set(Popup_RVE_particlesizeDmap,'enable','on');            
            set(Button_options_particlesizeDmap_RVE_new,'enable','on');            
            set(Button_options_particlesizeDmap_RVE_reset,'enable','on');            
            set(table_options_particlesizeDmap_RVE_subvolumes,'enable','on');            
        else
            set(checkbox_particlesizeDmap_todo_voxelsize,'enable','off');
            set(table_options_particlesizeDmap_voxelsize,'enable','off');
            set(Button_options_particlesizeDmap_voxelsize_add,'enable','off');
            set(Button_options_particlesizeDmap_voxelsize_remove,'enable','off');            
            set(Button_options_particlesizeDmap_voxelsize_update,'enable','off');            
            set(text_particlesizeDmap_todo_RVE,'enable','off');            
            set(Popup_RVE_particlesizeDmap,'enable','off');            
            set(Button_options_particlesizeDmap_RVE_new,'enable','off');            
            set(Button_options_particlesizeDmap_RVE_reset,'enable','off');            
            set(table_options_particlesizeDmap_RVE_subvolumes,'enable','off');                 
        end
    end

    function checkbox_particlesizeDmap_todo_voxelsize_Callback(~,~)
        PROPERTY.particlesize_dmap.voxel_size_dependence.todo = checkbox_particlesizeDmap_todo_voxelsize.Value;
        if PROPERTY.particlesize_dmap.voxel_size_dependence.todo
            set(table_options_particlesizeDmap_voxelsize,'enable','on');
            set(Button_options_particlesizeDmap_voxelsize_add,'enable','on');
            set(Button_options_particlesizeDmap_voxelsize_remove,'enable','on');
            set(Button_options_particlesizeDmap_voxelsize_update,'enable','on');
        else
            set(table_options_particlesizeDmap_voxelsize,'enable','off');
            set(Button_options_particlesizeDmap_voxelsize_add,'enable','off');
            set(Button_options_particlesizeDmap_voxelsize_remove,'enable','off');
            set(Button_options_particlesizeDmap_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_particlesizeDmap_voxelsize_add_Callback(~,~)
        val_ = table_options_particlesizeDmap_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_particlesizeDmap_voxelsize.Data = val_;
        PROPERTY.particlesize_dmap.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_particlesizeDmap_voxelsize_remove_Callback(~,~)
        val_ = table_options_particlesizeDmap_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_particlesizeDmap_voxelsize.Data = val_;
        PROPERTY.particlesize_dmap.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_particlesizeDmap_voxelsize_update_Callback(~,~)
        val_ = table_options_particlesizeDmap_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_particlesizeDmap_voxelsize.Data = val_;
        PROPERTY.particlesize_dmap.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tableparticlesizeDmapvoxelsize_Callback(~,~)
        val_ = table_options_particlesizeDmap_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_particlesizeDmap_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.particlesize_dmap.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_particlesizeDmap_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end


%%
%% TAB D-PSD (WATERSHED)
%%

%% GUI
Text_tab_particlesize_watershed = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the particle size (watershed immersion). (BETA: slow)'); 

checkbox_particlesize_watershed_todo = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate particle size with the watershed immersion method.','Value',PROPERTY.particlesize_watershed.todo,'Callback',{@checkbox_particlesize_watershed_todo_Callback}); 
annotation(tab_particlesize_watershed,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_particlesize_watershed_todo_voxelsize = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','enable','off',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.particlesize_watershed.voxel_size_dependence.todo,'Callback',{@checkbox_particlesize_watershed_todo_voxelsize_Callback}); 

table_options_particlesize_watershed_voxelsize = uitable('Parent', tab_particlesize_watershed,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tableparticlesize_watershedvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tableparticlesize_watershedvoxelsize_Callback,'enable','off');
table_options_particlesize_watershed_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_particlesize_watershed_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_particlesize_watershed_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_particlesize_watershed_todo_voxelsize.Position(3)}; % Auto width
table_options_particlesize_watershed_voxelsize.RowName = []; % Remove row name
table_options_particlesize_watershed_voxelsize.Data=PROPERTY.particlesize_watershed.voxel_size_dependence.voxel';

Button_options_particlesize_watershed_voxelsize_add = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'pushbutton', 'String', '+','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_watershed_voxelsize_add_Callback});
Button_options_particlesize_watershed_voxelsize_remove = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'pushbutton', 'String', '-','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_watershed_voxelsize_remove_Callback});
Button_options_particlesize_watershed_voxelsize_update = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'pushbutton', 'String', 'update','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_watershed_voxelsize_update_Callback});

text_particlesize_watershed_todo_RVE = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','enable','off',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_particlesize_watershed = uicontrol('Parent', tab_particlesize_watershed,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'enable','off',...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_particlesize_watershed_RVE_new = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'pushbutton', 'String', 'Add one RVE','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_watershed_RVE_add_Callback});
Button_options_particlesize_watershed_RVE_reset = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'pushbutton', 'String', 'Remove all RVE','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_watershed_RVE_reset_Callback});
       
table_options_particlesize_watershed_RVE_subvolumes = uitable('Parent', tab_particlesize_watershed,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'enable','off','Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_particlesize_watershed_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_particlesize_watershed_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_particlesize_watershed_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.particlesize_watershed.RVE(1).type};
column2={num2str(PROPERTY.particlesize_watershed.RVE(1).divisions)};
column3={PROPERTY.particlesize_watershed.RVE(1).subs2};
column4={PROPERTY.particlesize_watershed.RVE(1).subs4};
column5={PROPERTY.particlesize_watershed.RVE(1).Aspectratio};
column6={PROPERTY.particlesize_watershed.RVE(1).Constantdirection};
column7={PROPERTY.particlesize_watershed.RVE(1).Growthdirection};
column8={[PROPERTY.particlesize_watershed.RVE(1).Growthperstep ' ' PROPERTY.particlesize_watershed.RVE(1).Growthrelativeto]};
column9={[PROPERTY.particlesize_watershed.RVE(1).firstuniquevolume_size ' ' PROPERTY.particlesize_watershed.RVE(1).firstuniquevolume_unit]};
table_options_particlesize_watershed_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];


% Options
annotation(tab_particlesize_watershed,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.5 0.5],'visible','on');

checkbox_particlesize_watershed_custom = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.425 5*pos_delta_x 0.3*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w','enable','off',...
    'String','Use custom watershed function (it not, use MATLAB built-in watershed function).','Value',PROPERTY.particlesize_watershed.customfunction,'Callback',{@checkbox_particlesize_watershed_customfunction_Callback}); 
    function checkbox_particlesize_watershed_customfunction_Callback(~,~)
        if checkbox_particlesize_watershed_custom.Value
            set(checkbox_particlesize_watershed_customdetails,'enable','on');
        else
            set(checkbox_particlesize_watershed_customdetails,'enable','off');
        end
        PROPERTY.particlesize_watershed.customfunction = checkbox_particlesize_watershed_custom.Value;
    end

checkbox_particlesize_watershed_customdetails = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft+5*pos_delta_x 0.425 3.5*pos_delta_x 0.3*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Detail convergence','Value',PROPERTY.particlesize_watershed.details_convergence,'Callback',{@checkbox_particlesize_watershed_customdetail_Callback}); 
    function checkbox_particlesize_watershed_customdetail_Callback(~,~)
        PROPERTY.particlesize_watershed.details_convergence = checkbox_particlesize_watershed_customdetails.Value;
    end

checkbox_particlesize_watershed_cpsdrefining = uicontrol('Parent', tab_particlesize_watershed, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.425-0.3*pos_delta_y 5*pos_delta_x 0.3*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Correct oversegmentation using spherical particle size as lower bound','Value',PROPERTY.particlesize_watershed.cpsd_refining,'Callback',{@checkbox_particlesize_watershed_cpsdrefining_Callback}); 
    function checkbox_particlesize_watershed_cpsdrefining_Callback(~,~)
        PROPERTY.particlesize_watershed.cpsd_refining = checkbox_particlesize_watershed_cpsdrefining.Value;
    end


%% CALLBACKS

    function pushbutton_options_watershed_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_particlesize_watershed.String(Popup_RVE_particlesize_watershed.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'particlesize_watershed', table_options_particlesize_watershed_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'particlesize_watershed', table_options_particlesize_watershed_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'particlesize_watershed', table_options_particlesize_watershed_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'particlesize_watershed', table_options_particlesize_watershed_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'particlesize_watershed', table_options_particlesize_watershed_RVE_subvolumes);
        end
    end

    function pushbutton_options_watershed_RVE_reset_Callback(~,~)
        PROPERTY.particlesize_watershed.RVE=[];
        PROPERTY.particlesize_watershed.number_RVE = 0; % no RVE to perform
        table_options_particlesize_watershed_RVE_subvolumes.Data=[];
    end

    function checkbox_particlesize_watershed_todo_Callback(~,~)
        PROPERTY.particlesize_watershed.todo = checkbox_particlesize_watershed_todo.Value;
        if PROPERTY.particlesize_watershed.todo
            %set(checkbox_particlesize_watershed_todo_voxelsize,'enable','on');
            %checkbox_particlesize_watershed_todo_voxelsize_Callback
            %set(text_particlesize_watershed_todo_RVE,'enable','on');            
            %set(Popup_RVE_particlesize_watershed,'enable','on');            
            %set(Button_options_particlesize_watershed_RVE_new,'enable','on');            
            %set(Button_options_particlesize_watershed_RVE_reset,'enable','on');            
            %set(table_options_particlesize_watershed_RVE_subvolumes,'enable','on'); 
            %set(checkbox_particlesize_watershed_custom,'enable','on');
            set(checkbox_particlesize_watershed_customdetails,'enable','on');
            set(checkbox_particlesize_watershed_cpsdrefining,'enable','on');
        else
            %set(checkbox_particlesize_watershed_todo_voxelsize,'enable','off');
            %set(table_options_particlesize_watershed_voxelsize,'enable','off');
            %set(Button_options_particlesize_watershed_voxelsize_add,'enable','off');
            %set(Button_options_particlesize_watershed_voxelsize_remove,'enable','off');            
            %set(Button_options_particlesize_watershed_voxelsize_update,'enable','off');            
            %set(text_particlesize_watershed_todo_RVE,'enable','off');            
            %set(Popup_RVE_particlesize_watershed,'enable','off');            
            %set(Button_options_particlesize_watershed_RVE_new,'enable','off');            
            %set(Button_options_particlesize_watershed_RVE_reset,'enable','off');            
            %set(table_options_particlesize_watershed_RVE_subvolumes,'enable','off');    
            %set(checkbox_particlesize_watershed_custom,'enable','off');
            set(checkbox_particlesize_watershed_customdetails,'enable','off');
            set(checkbox_particlesize_watershed_cpsdrefining,'enable','off');            
        end
    end

    function checkbox_particlesize_watershed_todo_voxelsize_Callback(~,~)
        PROPERTY.particlesize_watershed.voxel_size_dependence.todo = checkbox_particlesize_watershed_todo_voxelsize.Value;
        if PROPERTY.particlesize_watershed.voxel_size_dependence.todo
            set(table_options_particlesize_watershed_voxelsize,'enable','on');
            set(Button_options_particlesize_watershed_voxelsize_add,'enable','on');
            set(Button_options_particlesize_watershed_voxelsize_remove,'enable','on');
            set(Button_options_particlesize_watershed_voxelsize_update,'enable','on');
        else
            set(table_options_particlesize_watershed_voxelsize,'enable','off');
            set(Button_options_particlesize_watershed_voxelsize_add,'enable','off');
            set(Button_options_particlesize_watershed_voxelsize_remove,'enable','off');
            set(Button_options_particlesize_watershed_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_watershed_voxelsize_add_Callback(~,~)
        val_ = table_options_particlesize_watershed_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_particlesize_watershed_voxelsize.Data = val_;
        PROPERTY.particlesize_watershed.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_watershed_voxelsize_remove_Callback(~,~)
        val_ = table_options_particlesize_watershed_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_particlesize_watershed_voxelsize.Data = val_;
        PROPERTY.particlesize_watershed.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_watershed_voxelsize_update_Callback(~,~)
        val_ = table_options_particlesize_watershed_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_particlesize_watershed_voxelsize.Data = val_;
        PROPERTY.particlesize_watershed.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tableparticlesize_watershedvoxelsize_Callback(~,~)
        val_ = table_options_particlesize_watershed_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_particlesize_watershed_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.particlesize_watershed.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_particlesize_watershed_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end



%%
%% TAB D-PSD (PCRF)
%%

%% GUI
Text_tab_particlesize_PCRF = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the particle size (Pseudo Coulomb Repulsive Field, PCRF). (BETA: slow)'); 

checkbox_particlesize_PCRF_todo = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate particle size with the PCRF method.','Value',PROPERTY.particlesize_PCRF.todo,'Callback',{@checkbox_particlesize_PCRF_todo_Callback}); 
annotation(tab_particlesize_PCRF,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_particlesize_PCRF_todo_voxelsize = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','enable','off',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.particlesize_PCRF.voxel_size_dependence.todo,'Callback',{@checkbox_particlesize_PCRF_todo_voxelsize_Callback}); 

table_options_particlesize_PCRF_voxelsize = uitable('Parent', tab_particlesize_PCRF,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tableparticlesize_PCRFvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tableparticlesize_PCRFvoxelsize_Callback,'enable','off');
table_options_particlesize_PCRF_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_particlesize_PCRF_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_particlesize_PCRF_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_particlesize_PCRF_todo_voxelsize.Position(3)}; % Auto width
table_options_particlesize_PCRF_voxelsize.RowName = []; % Remove row name
table_options_particlesize_PCRF_voxelsize.Data=PROPERTY.particlesize_PCRF.voxel_size_dependence.voxel';

Button_options_particlesize_PCRF_voxelsize_add = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'pushbutton', 'String', '+','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_PCRF_voxelsize_add_Callback});
Button_options_particlesize_PCRF_voxelsize_remove = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'pushbutton', 'String', '-','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_PCRF_voxelsize_remove_Callback});
Button_options_particlesize_PCRF_voxelsize_update = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'pushbutton', 'String', 'update','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_PCRF_voxelsize_update_Callback});

text_particlesize_PCRF_todo_RVE = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','enable','off',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_particlesize_PCRF = uicontrol('Parent', tab_particlesize_PCRF,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'enable','off',...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_particlesize_PCRF_RVE_new = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'pushbutton', 'String', 'Add one RVE','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_PCRF_RVE_add_Callback});
Button_options_particlesize_PCRF_RVE_reset = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'pushbutton', 'String', 'Remove all RVE','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_PCRF_RVE_reset_Callback});
       
table_options_particlesize_PCRF_RVE_subvolumes = uitable('Parent', tab_particlesize_PCRF,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'enable','off','Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_particlesize_PCRF_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_particlesize_PCRF_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_particlesize_PCRF_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.particlesize_PCRF.RVE(1).type};
column2={num2str(PROPERTY.particlesize_PCRF.RVE(1).divisions)};
column3={PROPERTY.particlesize_PCRF.RVE(1).subs2};
column4={PROPERTY.particlesize_PCRF.RVE(1).subs4};
column5={PROPERTY.particlesize_PCRF.RVE(1).Aspectratio};
column6={PROPERTY.particlesize_PCRF.RVE(1).Constantdirection};
column7={PROPERTY.particlesize_PCRF.RVE(1).Growthdirection};
column8={[PROPERTY.particlesize_PCRF.RVE(1).Growthperstep ' ' PROPERTY.particlesize_PCRF.RVE(1).Growthrelativeto]};
column9={[PROPERTY.particlesize_PCRF.RVE(1).firstuniquevolume_size ' ' PROPERTY.particlesize_PCRF.RVE(1).firstuniquevolume_unit]};
table_options_particlesize_PCRF_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];


% Options
annotation(tab_particlesize_PCRF,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.5 0.5],'visible','on');

checkbox_particlesize_PCRF_customdetails = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft+5*pos_delta_x 0.425 3.5*pos_delta_x 0.3*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Detail convergence','Value',PROPERTY.particlesize_PCRF.details_convergence,'Callback',{@checkbox_particlesize_PCRF_customdetail_Callback}); 
    function checkbox_particlesize_PCRF_customdetail_Callback(~,~)
        PROPERTY.particlesize_PCRF.details_convergence = checkbox_particlesize_PCRF_customdetails.Value;
    end

checkbox_particlesize_PCRF_cpsdrefining = uicontrol('Parent', tab_particlesize_PCRF, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.425-0.3*pos_delta_y 5*pos_delta_x 0.3*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Correct oversegmentation using spherical particle size as lower bound','Value',PROPERTY.particlesize_PCRF.cpsd_refining,'Callback',{@checkbox_particlesize_PCRF_cpsdrefining_Callback}); 
    function checkbox_particlesize_PCRF_cpsdrefining_Callback(~,~)
        PROPERTY.particlesize_PCRF.cpsd_refining = checkbox_particlesize_PCRF_cpsdrefining.Value;
    end


%% CALLBACKS

    function pushbutton_options_PCRF_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_particlesize_PCRF.String(Popup_RVE_particlesize_PCRF.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'particlesize_PCRF', table_options_particlesize_PCRF_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'particlesize_PCRF', table_options_particlesize_PCRF_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'particlesize_PCRF', table_options_particlesize_PCRF_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'particlesize_PCRF', table_options_particlesize_PCRF_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'particlesize_PCRF', table_options_particlesize_PCRF_RVE_subvolumes);
        end
    end

    function pushbutton_options_PCRF_RVE_reset_Callback(~,~)
        PROPERTY.particlesize_PCRF.RVE=[];
        PROPERTY.particlesize_PCRF.number_RVE = 0; % no RVE to perform
        table_options_particlesize_PCRF_RVE_subvolumes.Data=[];
    end

    function checkbox_particlesize_PCRF_todo_Callback(~,~)
        PROPERTY.particlesize_PCRF.todo = checkbox_particlesize_PCRF_todo.Value;
        if PROPERTY.particlesize_PCRF.todo
            %set(checkbox_particlesize_PCRF_todo_voxelsize,'enable','on');
            %checkbox_particlesize_PCRF_todo_voxelsize_Callback
            %set(text_particlesize_PCRF_todo_RVE,'enable','on');            
            %set(Popup_RVE_particlesize_PCRF,'enable','on');            
            %set(Button_options_particlesize_PCRF_RVE_new,'enable','on');            
            %set(Button_options_particlesize_PCRF_RVE_reset,'enable','on');            
            %set(table_options_particlesize_PCRF_RVE_subvolumes,'enable','on'); 
            set(checkbox_particlesize_PCRF_customdetails,'enable','on');
            set(checkbox_particlesize_PCRF_cpsdrefining,'enable','on');
        else
            %set(checkbox_particlesize_PCRF_todo_voxelsize,'enable','off');
            %set(table_options_particlesize_PCRF_voxelsize,'enable','off');
            %set(Button_options_particlesize_PCRF_voxelsize_add,'enable','off');
            %set(Button_options_particlesize_PCRF_voxelsize_remove,'enable','off');            
            %set(Button_options_particlesize_PCRF_voxelsize_update,'enable','off');            
            %set(text_particlesize_PCRF_todo_RVE,'enable','off');            
            %set(Popup_RVE_particlesize_PCRF,'enable','off');            
            %set(Button_options_particlesize_PCRF_RVE_new,'enable','off');            
            %set(Button_options_particlesize_PCRF_RVE_reset,'enable','off');            
            %set(table_options_particlesize_PCRF_RVE_subvolumes,'enable','off');    
            set(checkbox_particlesize_PCRF_customdetails,'enable','off');
            set(checkbox_particlesize_PCRF_cpsdrefining,'enable','off');            
        end
    end

    function checkbox_particlesize_PCRF_todo_voxelsize_Callback(~,~)
        PROPERTY.particlesize_PCRF.voxel_size_dependence.todo = checkbox_particlesize_PCRF_todo_voxelsize.Value;
        if PROPERTY.particlesize_PCRF.voxel_size_dependence.todo
            set(table_options_particlesize_PCRF_voxelsize,'enable','on');
            set(Button_options_particlesize_PCRF_voxelsize_add,'enable','on');
            set(Button_options_particlesize_PCRF_voxelsize_remove,'enable','on');
            set(Button_options_particlesize_PCRF_voxelsize_update,'enable','on');
        else
            set(table_options_particlesize_PCRF_voxelsize,'enable','off');
            set(Button_options_particlesize_PCRF_voxelsize_add,'enable','off');
            set(Button_options_particlesize_PCRF_voxelsize_remove,'enable','off');
            set(Button_options_particlesize_PCRF_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_PCRF_voxelsize_add_Callback(~,~)
        val_ = table_options_particlesize_PCRF_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_particlesize_PCRF_voxelsize.Data = val_;
        PROPERTY.particlesize_PCRF.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_PCRF_voxelsize_remove_Callback(~,~)
        val_ = table_options_particlesize_PCRF_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_particlesize_PCRF_voxelsize.Data = val_;
        PROPERTY.particlesize_PCRF.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_PCRF_voxelsize_update_Callback(~,~)
        val_ = table_options_particlesize_PCRF_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_particlesize_PCRF_voxelsize.Data = val_;
        PROPERTY.particlesize_PCRF.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tableparticlesize_PCRFvoxelsize_Callback(~,~)
        val_ = table_options_particlesize_PCRF_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_particlesize_PCRF_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.particlesize_PCRF.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_particlesize_PCRF_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end


%%
%% TAB TORTUOSITY FACTOR (TAU FACTOR)
%%

%% GUI
Text_tab_tortuosity_taufactor = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the tortuosity factor (Tau factor, by Dr. Samuel Cooper).'); 

checkbox_tortuosity_taufactor_todo = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate tortuosity factor with tau factor.','Value',PROPERTY.tortuosity_taufactor.todo,'Callback',{@checkbox_tortuosity_taufactor_todo_Callback}); 
annotation(tab_tortuosity_taufactor,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_tortuosity_taufactor_todo_voxelsize = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.tortuosity_taufactor.voxel_size_dependence.todo,'Callback',{@checkbox_tortuosity_taufactor_todo_voxelsize_Callback}); 

table_options_tortuosity_taufactor_voxelsize = uitable('Parent', tab_tortuosity_taufactor,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tabletortuosity_taufactorvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tabletortuosity_taufactorvoxelsize_Callback);
table_options_tortuosity_taufactor_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_tortuosity_taufactor_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_tortuosity_taufactor_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_tortuosity_taufactor_todo_voxelsize.Position(3)}; % Auto width
table_options_tortuosity_taufactor_voxelsize.RowName = []; % Remove row name
table_options_tortuosity_taufactor_voxelsize.Data=PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel';

Button_options_tortuosity_taufactor_voxelsize_add = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_taufactor_voxelsize_add_Callback});
Button_options_tortuosity_taufactor_voxelsize_remove = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_taufactor_voxelsize_remove_Callback});
Button_options_tortuosity_taufactor_voxelsize_update = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_taufactor_voxelsize_update_Callback});

text_tortuosity_taufactor_todo_RVE = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_tortuosity_taufactor = uicontrol('Parent', tab_tortuosity_taufactor,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_tortuosity_taufactor_RVE_new = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_taufactor_RVE_add_Callback});
Button_options_tortuosity_taufactor_RVE_reset = uicontrol('Parent', tab_tortuosity_taufactor, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_taufactor_RVE_reset_Callback});
       
table_options_tortuosity_taufactor_RVE_subvolumes = uitable('Parent', tab_tortuosity_taufactor,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_tortuosity_taufactor_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_tortuosity_taufactor_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_tortuosity_taufactor_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.tortuosity_taufactor.RVE(1).type};
column2={num2str(PROPERTY.tortuosity_taufactor.RVE(1).divisions)};
column3={PROPERTY.tortuosity_taufactor.RVE(1).subs2};
column4={PROPERTY.tortuosity_taufactor.RVE(1).subs4};
column5={PROPERTY.tortuosity_taufactor.RVE(1).Aspectratio};
column6={PROPERTY.tortuosity_taufactor.RVE(1).Constantdirection};
column7={PROPERTY.tortuosity_taufactor.RVE(1).Growthdirection};
column8={[PROPERTY.tortuosity_taufactor.RVE(1).Growthperstep ' ' PROPERTY.tortuosity_taufactor.RVE(1).Growthrelativeto]};
column9={[PROPERTY.tortuosity_taufactor.RVE(1).firstuniquevolume_size ' ' PROPERTY.tortuosity_taufactor.RVE(1).firstuniquevolume_unit]};
table_options_tortuosity_taufactor_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];


% Options
annotation(tab_tortuosity_taufactor,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.5 0.5],'visible','on');

Text_taufactor_perslice = uicontrol('Parent', tab_tortuosity_taufactor,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Tortuosity factor calculated thick slice per thick slice',...
     'HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.425 3.5*pos_delta_x 0.3*pos_delta_y]);
table_options_tortuosity_taufactor_slice = uitable('Parent', tab_tortuosity_taufactor,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[pos_x_start 0.425-0.125 3.5*pos_delta_x 0.125],'CellEditCallback',@cellsection_tabletortuosity_taufactorperslice_Callback, 'CellSelectionCallback',@cellsection_tabletortuosity_taufactorperslice_Callback);
table_options_tortuosity_taufactor_slice.ColumnName = {'Direction','Calculate','Slice parameter'}; % Column name
table_options_tortuosity_taufactor_slice.ColumnEditable = [false true true]; % Select column editable
tmp = 0.95*scrsz(4)*table_options_tortuosity_taufactor_slice.Position(3);
table_options_tortuosity_taufactor_slice.ColumnWidth = {tmp/4 tmp/4 tmp/2}; % Auto width
table_options_tortuosity_taufactor_slice.RowName = []; % Remove row name
column_1 = [{1};{2};{3}];
column_2 =PROPERTY.tortuosity_taufactor.todo_slice;
column_3 = PROPERTY.tortuosity_taufactor.sliceparameter;
Popup_tortuosity_taufactor_slice = uicontrol('Parent', tab_tortuosity_taufactor,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', Slice_choice,'Units','normalized','Position', [pos_x_start 0.3-0.025 3.5*pos_delta_x 0.025],'Callback', @popup_taufactor_slice_Callback);    
table_options_tortuosity_taufactor_slice.Data=[column_1 column_2 column_3];

% Text_taufactor_choicephases = uicontrol('Parent', tab_tortuosity_taufactor,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Calculate for the phases',...
%      'HorizontalAlignment','left','Units','normalized','Position', [0.7 0.425 0.95-0.7 0.3*pos_delta_y]);
% table_options_tortuosity_taufactor_phases = uitable('Parent', tab_tortuosity_taufactor,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.7 0.425-0.125 0.95-0.7 0.125],'CellEditCallback',@cellsection_tabletortuosity_phases_Callback, 'CellSelectionCallback',@cellsection_tabletortuosity_phases_Callback);

    function cellsection_tabletortuosity_taufactorperslice_Callback(~,~)
        val_ = table_options_tortuosity_taufactor_slice.Data;
        PROPERTY.tortuosity_taufactor.todo_slice = val_(:,2);
        PROPERTY.tortuosity_taufactor.sliceparameter = val_(:,3);
    end
    function popup_taufactor_slice_Callback(~,~)
        PROPERTY.tortuosity_taufactor.slicechoice = Popup_tortuosity_taufactor_slice.Value;
    end
%     function cellsection_tabletortuosity_phases_Callback(~,~)
%         foo=1;
%     end

%% CALLBACKS

    function pushbutton_options_taufactor_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_tortuosity_taufactor.String(Popup_RVE_tortuosity_taufactor.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'tortuosity_taufactor', table_options_tortuosity_taufactor_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'tortuosity_taufactor', table_options_tortuosity_taufactor_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'tortuosity_taufactor', table_options_tortuosity_taufactor_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'tortuosity_taufactor', table_options_tortuosity_taufactor_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'tortuosity_taufactor', table_options_tortuosity_taufactor_RVE_subvolumes);
        end
    end

    function pushbutton_options_taufactor_RVE_reset_Callback(~,~)
        PROPERTY.tortuosity_taufactor.RVE=[];
        PROPERTY.tortuosity_taufactor.number_RVE = 0; % no RVE to perform
        table_options_tortuosity_taufactor_RVE_subvolumes.Data=[];
    end

    function checkbox_tortuosity_taufactor_todo_Callback(~,~)
        PROPERTY.tortuosity_taufactor.todo = checkbox_tortuosity_taufactor_todo.Value;
        if PROPERTY.tortuosity_taufactor.todo
            set(checkbox_tortuosity_taufactor_todo_voxelsize,'enable','on');
            checkbox_tortuosity_taufactor_todo_voxelsize_Callback
            set(text_tortuosity_taufactor_todo_RVE,'enable','on');            
            set(Popup_RVE_tortuosity_taufactor,'enable','on');            
            set(Button_options_tortuosity_taufactor_RVE_new,'enable','on');            
            set(Button_options_tortuosity_taufactor_RVE_reset,'enable','on');            
            set(table_options_tortuosity_taufactor_RVE_subvolumes,'enable','on'); 
            set(Text_taufactor_perslice,'enable','on');             
            set(table_options_tortuosity_taufactor_slice,'enable','on');             
            set(Popup_tortuosity_taufactor_slice,'enable','on');             
            %set(Text_taufactor_choicephases,'enable','on');             
            %set(table_options_tortuosity_taufactor_phases,'enable','on');             
        else
            set(checkbox_tortuosity_taufactor_todo_voxelsize,'enable','off');
            set(table_options_tortuosity_taufactor_voxelsize,'enable','off');
            set(Button_options_tortuosity_taufactor_voxelsize_add,'enable','off');
            set(Button_options_tortuosity_taufactor_voxelsize_remove,'enable','off');            
            set(Button_options_tortuosity_taufactor_voxelsize_update,'enable','off');            
            set(text_tortuosity_taufactor_todo_RVE,'enable','off');            
            set(Popup_RVE_tortuosity_taufactor,'enable','off');            
            set(Button_options_tortuosity_taufactor_RVE_new,'enable','off');            
            set(Button_options_tortuosity_taufactor_RVE_reset,'enable','off');            
            set(table_options_tortuosity_taufactor_RVE_subvolumes,'enable','off');    
            set(Text_taufactor_perslice,'enable','off');             
            set(table_options_tortuosity_taufactor_slice,'enable','off');             
            set(Popup_tortuosity_taufactor_slice,'enable','off');             
            %set(Text_taufactor_choicephases,'enable','off');             
            %set(table_options_tortuosity_taufactor_phases,'enable','off');                
        end
    end

    function checkbox_tortuosity_taufactor_todo_voxelsize_Callback(~,~)
        PROPERTY.tortuosity_taufactor.voxel_size_dependence.todo = checkbox_tortuosity_taufactor_todo_voxelsize.Value;
        if PROPERTY.tortuosity_taufactor.voxel_size_dependence.todo
            set(table_options_tortuosity_taufactor_voxelsize,'enable','on');
            set(Button_options_tortuosity_taufactor_voxelsize_add,'enable','on');
            set(Button_options_tortuosity_taufactor_voxelsize_remove,'enable','on');
            set(Button_options_tortuosity_taufactor_voxelsize_update,'enable','on');
        else
            set(table_options_tortuosity_taufactor_voxelsize,'enable','off');
            set(Button_options_tortuosity_taufactor_voxelsize_add,'enable','off');
            set(Button_options_tortuosity_taufactor_voxelsize_remove,'enable','off');
            set(Button_options_tortuosity_taufactor_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_taufactor_voxelsize_add_Callback(~,~)
        val_ = table_options_tortuosity_taufactor_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_tortuosity_taufactor_voxelsize.Data = val_;
        PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_taufactor_voxelsize_remove_Callback(~,~)
        val_ = table_options_tortuosity_taufactor_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_tortuosity_taufactor_voxelsize.Data = val_;
        PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_taufactor_voxelsize_update_Callback(~,~)
        val_ = table_options_tortuosity_taufactor_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_tortuosity_taufactor_voxelsize.Data = val_;
        PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tabletortuosity_taufactorvoxelsize_Callback(~,~)
        val_ = table_options_tortuosity_taufactor_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_tortuosity_taufactor_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_tortuosity_taufactor_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end


%%
%% TAB CONNECTIVITY
%%

%% GUI
Text_tab_connectivity = uicontrol('Parent', tab_connectivity, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Setup options for the connectivity (percolation).'); 

checkbox_connectivity_todo = uicontrol('Parent', tab_connectivity, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft 0.87 description_tab_xlenght 0.4*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Calculate connectivity.','Value',PROPERTY.connectivity.todo,'Callback',{@checkbox_connectivity_todo_Callback}); 
annotation(tab_connectivity,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.87 0.87],'visible','on');

stepy = 0.05;
checkbox_connectivity_todo_voxelsize = uicontrol('Parent', tab_connectivity, 'Style', 'checkbox','Units','normalized','Position',[description_tab_fromleft pos_y_start-stepy 2.5*pos_delta_x 0.5*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Perform voxel size dependence analysis.','Value',PROPERTY.connectivity.voxel_size_dependence.todo,'Callback',{@checkbox_connectivity_todo_voxelsize_Callback}); 

table_options_connectivity_voxelsize = uitable('Parent', tab_connectivity,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.64-stepy 2.5*pos_delta_x 0.2],'CellEditCallback',@cellsection_tableconnectivityvoxelsize_Callback, 'CellSelectionCallback',@cellsection_tableconnectivityvoxelsize_Callback);
table_options_connectivity_voxelsize.ColumnName = {'Voxel size multiplier'}; % Column name
table_options_connectivity_voxelsize.ColumnEditable = [true]; % Select column editable
table_options_connectivity_voxelsize.ColumnWidth = {0.95*scrsz(4)*checkbox_connectivity_todo_voxelsize.Position(3)}; % Auto width
table_options_connectivity_voxelsize.RowName = []; % Remove row name
table_options_connectivity_voxelsize.Data=PROPERTY.connectivity.voxel_size_dependence.voxel';

Button_options_connectivity_voxelsize_add = uicontrol('Parent', tab_connectivity, 'Style', 'pushbutton', 'String', '+',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_connectivity_voxelsize_add_Callback});
Button_options_connectivity_voxelsize_remove = uicontrol('Parent', tab_connectivity, 'Style', 'pushbutton', 'String', '-',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(0.8+0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_connectivity_voxelsize_remove_Callback});
Button_options_connectivity_voxelsize_update = uicontrol('Parent', tab_connectivity, 'Style', 'pushbutton', 'String', 'update',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [description_tab_fromleft+(pos_delta_x*(2*0.8+2*0.05)) 0.61-stepy pos_delta_x*0.8 0.03],'Callback',{@pushbutton_options_connectivity_voxelsize_update_Callback});

text_connectivity_todo_RVE = uicontrol('Parent', tab_connectivity, 'Style', 'text','Units','normalized','Position',[0.3 pos_y_start+0.3*pos_delta_y-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2*pos_delta_y],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Representative Volume Element analysis (RVA). Choose below the RVA you want to perform and add it.','HorizontalAlignment','left'); 

Popup_RVE_connectivity = uicontrol('Parent', tab_connectivity,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', RVE_choice,'Units','normalized','Position', [0.3 pos_y_start-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.3*pos_delta_y]);    
    
Button_options_connectivity_RVE_new = uicontrol('Parent', tab_connectivity, 'Style', 'pushbutton', 'String', 'Add one RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_connectivity_RVE_add_Callback});
Button_options_connectivity_RVE_reset = uicontrol('Parent', tab_connectivity, 'Style', 'pushbutton', 'String', 'Remove all RVE',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.4 0.81-stepy pos_delta_x*0.9 0.03],'Callback',{@pushbutton_options_connectivity_RVE_reset_Callback});
       
table_options_connectivity_RVE_subvolumes = uitable('Parent', tab_connectivity,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.3 0.6-stepy description_tab_xlenght-(0.3-description_tab_fromleft) 0.2]);
table_options_connectivity_RVE_subvolumes.ColumnName = {'RVE','Divisions','2 Subs','4 Subs','Aspect ratio','Constant direction length','Growth direction','Growth per step','First volume size'}; % Column name
table_options_connectivity_RVE_subvolumes.ColumnEditable = [false false false false false false false false false]; % Select column editable
table_options_connectivity_RVE_subvolumes.RowName = []; % Remove row name

column1={PROPERTY.connectivity.RVE(1).type};
column2={num2str(PROPERTY.connectivity.RVE(1).divisions)};
column3={PROPERTY.connectivity.RVE(1).subs2};
column4={PROPERTY.connectivity.RVE(1).subs4};
column5={PROPERTY.connectivity.RVE(1).Aspectratio};
column6={PROPERTY.connectivity.RVE(1).Constantdirection};
column7={PROPERTY.connectivity.RVE(1).Growthdirection};
column8={[PROPERTY.connectivity.RVE(1).Growthperstep ' ' PROPERTY.connectivity.RVE(1).Growthrelativeto]};
column9={[PROPERTY.connectivity.RVE(1).firstuniquevolume_size ' ' PROPERTY.connectivity.RVE(1).firstuniquevolume_unit]};
table_options_connectivity_RVE_subvolumes.Data=[column1 column2 column3 column4 column5 column6 column7 column8 column9];


% Options
annotation(tab_connectivity,'line',[description_tab_fromleft description_tab_fromleft+description_tab_xlenght],[0.5 0.5],'visible','on');

Text_connectivity_perslice = uicontrol('Parent', tab_connectivity,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Connectivity calculated thick slice per thick slice',...
     'HorizontalAlignment','left','Units','normalized','Position', [pos_x_start 0.425 3.5*pos_delta_x 0.3*pos_delta_y]);
table_options_connectivity_slice = uitable('Parent', tab_connectivity,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[pos_x_start 0.425-0.125 3.5*pos_delta_x 0.125],'CellEditCallback',@cellsection_tableconnectivityperslice_Callback, 'CellSelectionCallback',@cellsection_tableconnectivityperslice_Callback);
table_options_connectivity_slice.ColumnName = {'Direction','Calculate','Slice parameter'}; % Column name
table_options_connectivity_slice.ColumnEditable = [false true true]; % Select column editable
tmp = 0.95*scrsz(4)*table_options_connectivity_slice.Position(3);
table_options_connectivity_slice.ColumnWidth = {tmp/4 tmp/4 tmp/2}; % Auto width
table_options_connectivity_slice.RowName = []; % Remove row name
column_1 = [{1};{2};{3}];
column_2 =PROPERTY.connectivity.todo_slice;
column_3 = PROPERTY.connectivity.sliceparameter;
Popup_connectivity_slice = uicontrol('Parent', tab_connectivity,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
           'String', Slice_choice,'Units','normalized','Position', [pos_x_start 0.3-0.025 3.5*pos_delta_x 0.025],'Callback', @popup_connectivity_slice_Callback);    
table_options_connectivity_slice.Data=[column_1 column_2 column_3];

% Text_connectivity_choicephases = uicontrol('Parent', tab_connectivity,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Calculate for the phases',...
%      'HorizontalAlignment','left','Units','normalized','Position', [0.7 0.425 0.95-0.7 0.3*pos_delta_y]);
% table_options_connectivity_phases = uitable('Parent', tab_connectivity,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.7 0.425-0.125 0.95-0.7 0.125],'CellEditCallback',@cellsection_tabletortuosity_phases_Callback, 'CellSelectionCallback',@cellsection_tabletortuosity_phases_Callback);

    function cellsection_tableconnectivityperslice_Callback(~,~)
        val_ = table_options_connectivity_slice.Data;
        PROPERTY.connectivity.todo_slice = val_(:,2);
        PROPERTY.connectivity.sliceparameter = val_(:,3);
    end
    function popup_connectivity_slice_Callback(~,~)
        PROPERTY.connectivity.slicechoice = Popup_connectivity_slice.Value;
    end
%     function cellsection_tabletortuosity_phases_Callback(~,~)
%         foo=1;
%     end

%% CALLBACKS

    function pushbutton_options_connectivity_RVE_add_Callback(~,~)
        str_ = char(Popup_RVE_connectivity.String(Popup_RVE_connectivity.Value));
        if strcmp(str_,'Independant subvolumes of same size + keep initial aspect ratio (A)')
            GUI_options_RVE('A', 'connectivity', table_options_connectivity_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + user-defined aspect ratio (B)')
            GUI_options_RVE('B', 'connectivity', table_options_connectivity_RVE_subvolumes);
        elseif strcmp(str_,'Independant subvolumes of same size + constant length (C)')
            GUI_options_RVE('C', 'connectivity', table_options_connectivity_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume center (D)')
            GUI_options_RVE('D', 'connectivity', table_options_connectivity_RVE_subvolumes);
        elseif strcmp(str_,'One subvolume + growing from volume edge (E)')
            GUI_options_RVE('E', 'connectivity', table_options_connectivity_RVE_subvolumes);
        end
    end

    function pushbutton_options_connectivity_RVE_reset_Callback(~,~)
        PROPERTY.connectivity.RVE=[];
        PROPERTY.connectivity.number_RVE = 0; % no RVE to perform
        table_options_connectivity_RVE_subvolumes.Data=[];
    end

    function checkbox_connectivity_todo_Callback(~,~)
        PROPERTY.connectivity.todo = checkbox_connectivity_todo.Value;
        if PROPERTY.connectivity.todo
            set(checkbox_connectivity_todo_voxelsize,'enable','on');
            checkbox_connectivity_todo_voxelsize_Callback
            set(text_connectivity_todo_RVE,'enable','on');            
            set(Popup_RVE_connectivity,'enable','on');            
            set(Button_options_connectivity_RVE_new,'enable','on');            
            set(Button_options_connectivity_RVE_reset,'enable','on');            
            set(table_options_connectivity_RVE_subvolumes,'enable','on'); 
            set(Text_connectivity_perslice,'enable','on');             
            set(table_options_connectivity_slice,'enable','on');             
            set(Popup_connectivity_slice,'enable','on');             
            %set(Text_connectivity_choicephases,'enable','on');             
            %set(table_options_connectivity_phases,'enable','on');             
        else
            set(checkbox_connectivity_todo_voxelsize,'enable','off');
            set(table_options_connectivity_voxelsize,'enable','off');
            set(Button_options_connectivity_voxelsize_add,'enable','off');
            set(Button_options_connectivity_voxelsize_remove,'enable','off');            
            set(Button_options_connectivity_voxelsize_update,'enable','off');            
            set(text_connectivity_todo_RVE,'enable','off');            
            set(Popup_RVE_connectivity,'enable','off');            
            set(Button_options_connectivity_RVE_new,'enable','off');            
            set(Button_options_connectivity_RVE_reset,'enable','off');            
            set(table_options_connectivity_RVE_subvolumes,'enable','off');    
            set(Text_connectivity_perslice,'enable','off');             
            set(table_options_connectivity_slice,'enable','off');             
            set(Popup_connectivity_slice,'enable','off');             
            %set(Text_connectivity_choicephases,'enable','off');             
            %set(table_options_connectivity_phases,'enable','off');                
        end
    end

    function checkbox_connectivity_todo_voxelsize_Callback(~,~)
        PROPERTY.connectivity.voxel_size_dependence.todo = checkbox_connectivity_todo_voxelsize.Value;
        if PROPERTY.connectivity.voxel_size_dependence.todo
            set(table_options_connectivity_voxelsize,'enable','on');
            set(Button_options_connectivity_voxelsize_add,'enable','on');
            set(Button_options_connectivity_voxelsize_remove,'enable','on');
            set(Button_options_connectivity_voxelsize_update,'enable','on');
        else
            set(table_options_connectivity_voxelsize,'enable','off');
            set(Button_options_connectivity_voxelsize_add,'enable','off');
            set(Button_options_connectivity_voxelsize_remove,'enable','off');
            set(Button_options_connectivity_voxelsize_update,'enable','off');
        end
    end

    function pushbutton_options_connectivity_voxelsize_add_Callback(~,~)
        val_ = table_options_connectivity_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_= [val_; val_(end)+1];
        table_options_connectivity_voxelsize.Data = val_;
        PROPERTY.connectivity.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_connectivity_voxelsize_remove_Callback(~,~)
        val_ = table_options_connectivity_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        val_(end)=[];
        if isempty(val_)
            val_=2;
        end
        table_options_connectivity_voxelsize.Data = val_;
        PROPERTY.connectivity.voxel_size_dependence.voxel = val_;
    end
    function pushbutton_options_connectivity_voxelsize_update_Callback(~,~)
        val_ = table_options_connectivity_voxelsize.Data;
        val_(val_==1)=[]; % Un-necessary
        val_ = unique(val_);
        val_=sort(val_);
        table_options_connectivity_voxelsize.Data = val_;
        PROPERTY.connectivity.voxel_size_dependence.voxel = val_;
    end

    function cellsection_tableconnectivityvoxelsize_Callback(~,~)
        val_ = table_options_connectivity_voxelsize.Data;
        no_nan = ~isnan(val_);
        is_positive=1;
        required_update=false;
        if no_nan
            positive = val_>0;
            if min(positive)==0
                is_positive=0;
            end
        end
        if no_nan*is_positive
            set(checkbox_connectivity_todo_voxelsize,'String','Perform voxel size dependence analysis.','BackgroundColor',[0.9400 0.9400 0.9400]);
            PROPERTY.connectivity.voxel_size_dependence.voxel = val_;
        else
            set(checkbox_connectivity_todo_voxelsize,'String','Error! Must be >0','BackgroundColor','r');
        end
    end


%%
%% TAL ALL LOADED VOLUME
%%

% Update table
    function Update_tab_importedvolume_Callback(~,~)

        % Check if some data exist
        if ~isempty(table_organize_volume.Data) && new_volume==0
          
            % Get back order of calculation
            table_data = table_organize_volume.Data;
            Orderofcalculation=table_data(:,4);
            Orderofcalculation = cell2mat(Orderofcalculation);
            % Copy of inputvolume
            copy_inputvolume=inputvolume;
            inputvolume=[];
            
            % Loop over order of calculation
            number_volume_to_analyse=0;
            for k=1:1:length(Orderofcalculation)
                % Different from 0
                if Orderofcalculation(k)~=0
                    number_volume_to_analyse=number_volume_to_analyse+1;
                    inputvolume.volume(number_volume_to_analyse)=copy_inputvolume.volume(k);
                    inputvolume.volume(number_volume_to_analyse).orderofcalculation=Orderofcalculation(k);
                end
            end
        end
        
        % Update GUI
        if number_volume_to_analyse==0
            % no data to analyse
            set(Text_nodata_tabloaded,'Visible','on');
            set(Text_instruction_tabloaded,'Visible','off');
            set(Text_nodata_runcalculation,'Visible','on');
            set(Text_tab_imported,'Visible','off');
            set(table_organize_volume,'Visible','off','enable','off');
            set(Button_changeorganise,'Visible','off','enable','off');
            set(Button_runcalculation,'Visible','off','enable','off');
            table_organize_volume.Data=[]; % Empty
        else
            % Extract data from inputvolume to be displayed on this table
            organisevolume_data = [{inputvolume.volume.filename_input}' {inputvolume.volume.equivalentcubicsize}' {inputvolume.volume.newvoxelsize}' {inputvolume.volume.orderofcalculation}'];
            % Update table
            table_organize_volume.Data=organisevolume_data;
            
            % Update table used to select phases for which parameters are calculated
            
            
            
        end
        % Reset
        new_volume=0;
    end

% Edit table
  function celledit_organisevolume_Callback(~,~)
        % Get back data of the table
        table_data = table_organize_volume.Data;
        Orderofcalculation=table_data(:,4);
        Orderofcalculation = cell2mat(Orderofcalculation);
        
        % Check validity of order of calculation
        % Only number
        detect_nan = isnan(Orderofcalculation);
        tmp=sum(sum(detect_nan));
        if tmp==0
            Orderofcalculation_hasnot_nan=1;
        else
            Orderofcalculation_hasnot_nan=0;
        end
 
        % Check integer
        if Orderofcalculation_hasnot_nan==1
            integer_orderofcalculation = round(Orderofcalculation);
            tmp = integer_orderofcalculation-Orderofcalculation;
            tmp2=nansum(nansum(tmp));
            if tmp2==0
                Orderofcalculation_has_integer=1;
            else
                Orderofcalculation_has_integer=1;
            end
        else
            Orderofcalculation_has_integer=0;
        end
        
        % Check positive or zero
        if  Orderofcalculation_hasnot_nan==1
            tmp=Orderofcalculation;
            tmp(Orderofcalculation>=0)=1;
            if (length(unique(tmp))==1 && unique(tmp)==1)
                Orderofcalculation_has_onlypositive=1;
            else
                Orderofcalculation_has_onlypositive=0;
            end
        else
            Orderofcalculation_has_onlypositive=0;
        end
        
        % Check different order of calculation
        if  Orderofcalculation_hasnot_nan==1
            % Select only non zero values
            nonzerovalues=Orderofcalculation(Orderofcalculation~=0);           
            unique_nonzerovalues=unique(nonzerovalues);
            if length(nonzerovalues)==length(unique_nonzerovalues)
                Orderofcalculation_hasuniquevalues=1;
            else
                Orderofcalculation_hasuniquevalues=0;
            end
        else
            Orderofcalculation_hasuniquevalues=0;
        end
        
        % CHeck order of calculation start at 1 and is continous (for non zero values)
        if Orderofcalculation_hasnot_nan==1 && Orderofcalculation_hasuniquevalues==1
            % Select only non zero values
            nonzerovalues=Orderofcalculation(Orderofcalculation~=0);
            if isempty(nonzerovalues)
                Orderofcalculation_hasorder=1;
            else
                if (min(nonzerovalues)==1 && max(Orderofcalculation)==length(nonzerovalues))
                    Orderofcalculation_hasorder=1;
                else
                    Orderofcalculation_hasorder=0;
                end
            end
        else
            Orderofcalculation_hasorder=0;
        end
        
        % Thus, is order of calculaion is correct?
        correct_orderofcalculation = Orderofcalculation_hasnot_nan*Orderofcalculation_has_integer*Orderofcalculation_has_onlypositive*Orderofcalculation_hasuniquevalues*Orderofcalculation_hasorder;
        
        % Display
        if correct_orderofcalculation==1
            % Re-enable update button
            set(Button_changeorganise,'enable','on');
            set(Text_wrongorderofcalculation,'Visible','off');
        else
            % Disabled update button
            set(Button_changeorganise,'enable','off');
            set(Text_wrongorderofcalculation,'Visible','on');
        end
        
  end

%%
%% TAB DISPLAY/SAVE OPTIONS
%%

Text_tab_import = uicontrol('Parent', tab_savedisplayoptions, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Select options for your results and figures.');


%% REPORT GENERATOR OPTIONS (BETA)
% Report options
import mlreportgen.dom.*;
OPTIONS.report.todo = true;
OPTIONS.report.rpt_type = 'docx';
OPTIONS.report.style.title = {Color('black'), FontFamily('Times'), FontSize('28pt'), Bold(true),  HAlign('Center'), LineSpacing(2)};
OPTIONS.report.style.title2 = {Color('blue'), FontFamily('Times'), FontSize('22pt'), Bold(false), HAlign('Center'), LineSpacing(1.5)};
OPTIONS.report.style.standardtext = {Color('black'), FontFamily('Times'), FontSize('12pt')};

%% HARDCODED
OPTIONS.arrayformat = 'uint8';
OPTIONS.structurefieldname = 'none';

%% GUI
table_displayoptions_results = uitable('Parent', tab_savedisplayoptions,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.75 0.3 0.15]);
table_displayoptions_results.ColumnName = {'Parameters','Choice'}; % Column name
table_displayoptions_results.ColumnEditable = [false, true]; % Select column editable
table_displayoptions_results.ColumnWidth = {220, 'auto'}; % Auto width
table_displayoptions_results.RowName = []; % Remove row name
table_displayoptions_results.Data=[[{'Save results in .mat'}; {'Save results in .xls'}; {'Save figures in .fig'}; {'Save figures in .png'}; {'Print results in the command window'}] [{true}; {true}; {true}; {true}; {true}]];

table_displayoptions_figure = uitable('Parent', tab_savedisplayoptions,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[description_tab_fromleft 0.5 0.3 0.2]);
table_displayoptions_figure.ColumnName = {'Parameters','Choice'}; % Column name
table_displayoptions_figure.ColumnEditable = [false, true]; % Select column editable
table_displayoptions_figure.ColumnWidth = {220, 'auto'}; % Auto width
table_displayoptions_figure.RowName = []; % Remove row name
table_displayoptions_figure.Data=[[{'Grid'}; {'Minor grid'}; {'Linewidth'}; {'Fontsize axe'}; {'Fontsize legend'}; {'Fontsize title'}; {'Close figure once displayed and/or saved'}] [{true}; {false}; {2}; {12}; {12}; {14}; {true}]];

text_displayoption_fontname = uicontrol('Parent', tab_savedisplayoptions, 'Style', 'text','Units','normalized','Position',[0.375 0.87 0.2 0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Select fontname','HorizontalAlignment','left');
Popup_displayoptions_fontname = uicontrol('Parent', tab_savedisplayoptions,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', listfonts,'Value',203,'Units','normalized','Position', [0.375 0.77 0.2 0.1]);
           

%%
%% TAB RUN CALCULATIONS
%%

%% GUI
Text_nodata_runcalculation = uicontrol('Parent', tab_run,'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','No data. Go to ''Import segmented volume'' tab','Visible','on',...
    'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0 0.4 1 0.2]);

% Description and instructions
Text_tab_run = uicontrol('Parent', tab_run, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],'Visible','off',...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Click to run calculations');
 
% Run calculation
Button_runcalculation = uicontrol('Parent', tab_run, 'Style', 'pushbutton', 'String', '>>>  Click to run ALL calculation on ALL volumes, according to their order of calculation  <<<',...
    'Visible','off','enable','off','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.2 0.5 0.6 0.1],'Callback',{@pushbutton_runallcalculations_Callback});

%% CALLBACK
    function pushbutton_runallcalculations_Callback(~,~)
        
        % Set OPTIONS
        val_ = table_displayoptions_results.Data(:,2);
        % Save options
        OPTIONS.save_resultsmat = logical(cell2mat(val_(1)));
        OPTIONS.save_xls = logical(cell2mat(val_(2)));
        OPTIONS.save_fig = logical(cell2mat(val_(3)));
        if logical(cell2mat(val_(4)))
            OPTIONS.savefig_informat = {'png'};
        else
            OPTIONS.savefig_informat = [];
        end
        OPTIONS.displaytext = logical(cell2mat(val_(5)));
        OPTIONS.savefig_infig = logical(logical(cell2mat(val_(3))) + logical(cell2mat(val_(4))));
        
        val_ = table_displayoptions_figure.Data(:,2);
        if logical(cell2mat(val_(1)))
            OPTIONS.grid = 'on';
        else
            OPTIONS.grid = 'off';
        end
        if logical(cell2mat(val_(2)))
            OPTIONS.minorgrid = 'on';
        else
            OPTIONS.minorgrid = 'off';
        end
        OPTIONS.Linewidth =  cell2mat(val_(3));
        OPTIONS.Fontsize_axe =  cell2mat(val_(4));
        OPTIONS.Fontsize_legend =  cell2mat(val_(5));
        OPTIONS.Fontsize_title =  cell2mat(val_(6));
        OPTIONS.closefigureaftercreation =  logical(cell2mat(val_(7)));
        OPTIONS.fontname = char(Popup_displayoptions_fontname.String(Popup_displayoptions_fontname.Value));
        
        for k=1:1:number_volume_to_analyse % Loop over the different file
            clear INFO
            
            % Input/Output
            OPTIONS.loadingpath = inputvolume.volume(k).loadingpath;
            OPTIONS.filename = inputvolume.volume(k).filename_input;
            OPTIONS.mainsavefolder = [inputvolume.volume(k).resultfolderlocation inputvolume.volume(k).resultfoldername];            
            
            % Prepare data required for the calculation
            % % INFO: volume related information
            INFO.RegionOfInterest = inputvolume.volume(k).ROI;
            INFO.initial_voxelsize = inputvolume.volume(k).initialvoxelsize;
            INFO.asked_voxelsize = inputvolume.volume(k).newvoxelsize;
            INFO.volumeinformation = inputvolume.volume(k).volumeinformation;
            
            % Phase name and code
            INFO.initial_phasecode = inputvolume.volume(k).phasecode;
            INFO.initial_phasename = inputvolume.volume(k).phasename;
            INFO.assigned_phasecode = inputvolume.volume(k).assignedcode;
            unique_assigncode = unique(inputvolume.volume(k).assignedcode);
            number_assignedcode = length(unique_assigncode);
            name_length=zeros(number_assignedcode,1);
            for current_phase = 1:1:number_assignedcode
                code_ = unique_assigncode(current_phase);
                id_ = find(inputvolume.volume(k).assignedcode == code_);
                INFO.phase(current_phase).code = code_;    
                INFO.phase(current_phase).color = inputvolume.volume(k).phasecolor(id_(1),:); 
                if length(id_)==1
                    INFO.phase(current_phase).name = char(inputvolume.volume(k).phasename(id_));
                else
                    str_ = char(inputvolume.volume(k).phasename(id_(1))); % initialization
                    for k_id=2:1:length(id_)
                        str_ = [str_ ' and ' char(inputvolume.volume(k).phasename(id_(k_id)))];
                    end
                    INFO.phase(current_phase).name = str_;
                end
                name_length(current_phase,1)=length(INFO.phase(current_phase).name);
                INFO.phasename(current_phase,1)={INFO.phase(current_phase).name};
                % Without space
                INFO.phase(current_phase).filename = function_remove_emptyandspecialcharacter_string( INFO.phase(current_phase).name );                
            end
            max_length = max(name_length);
            for current_phase=1:1:number_assignedcode
                % String with same length
                INFO.phase(current_phase).consistentname = function_enforcesamelength_string( INFO.phase(current_phase).name ,max_length);
            end
            % Store phase code
            n_ = length(INFO.phase);
            tmp=zeros(n_,2);
            for l=1:1:n_
                tmp(l,1)=INFO.phase(l).code;
            end
            INFO.phaseinfo=tmp; clear tmp;
            
            % Direction name
            number_direction = length(inputvolume.volume(k).directioname);
            name_length=zeros(number_direction,1);
            for current_direction=1:1:number_direction
                % For figure
                INFO.direction(current_direction).name = char(inputvolume.volume(k).directioname(current_direction));
                name_length(current_direction,1)=length(INFO.direction(current_direction).name);
                % Without space
                INFO.direction(current_direction).filename = function_remove_emptyandspecialcharacter_string( INFO.direction(current_direction).name );
            end
            max_length = max(name_length);
            for current_direction=1:1:number_direction
                % String with same length
                INFO.direction(current_direction).consistentname = function_enforcesamelength_string( INFO.direction(current_direction).name ,max_length);
            end                        
            % Call function
            function_microstructure_characterization_segmentedvolume(OPTIONS, INFO, PROPERTY)
        end
    end


%%
%% RVE GUI
%%

    function GUI_options_RVE(typeRVE, p, h)
        
        RVEparameters.type=typeRVE;
        RVEparameters.divisions='n/a';
        RVEparameters.subs2='n/a';
        RVEparameters.subs4='n/a';
        RVEparameters.Aspectratio='n/a';
        RVEparameters.Constantdirection='n/a';
        RVEparameters.Growthdirection='n/a';
        RVEparameters.Growthperstep='n/a';
        RVEparameters.Growthrelativeto='n/a';
        RVEparameters.firstuniquevolume_size='n/a';
        RVEparameters.firstuniquevolume_unit='n/a';
        RVEparameters.threshold_std='n/a';
        RVEparameters.threshold_numbersubvolumes='n/a';
        
        fig_RVE = figure;
        fig_RVE.Name= 'RVE options'; % Set name
        fig_RVE.NumberTitle='off'; % Remove number from title name
        fig_RVE.Color='white'; % Background colour
        fig_RVE.MenuBar='none'; % Remove menubar and toolbar
        fig_RVE.ToolBar='none';
        fig_RVE.Units='normalized'; % Set unit
        fig_RVE.Position=[0.3 0.3 0.4 0.4]; % Set Position
        
        if strcmp(typeRVE,'A') || strcmp(typeRVE,'B') || strcmp(typeRVE,'C')
            table_RVE_divisions = uitable('Parent', fig_RVE,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.05 0.55 0.5 0.4],'CellEditCallback',@cellsection_RVEdivision_Callback, 'CellSelectionCallback',@cellsection_RVEdivision_Callback);
            table_RVE_divisions.ColumnName = {'Volume division (n*n*n) or (n*n*1)'}; % Column name
            table_RVE_divisions.ColumnEditable = [true]; % Select column editable
            table_RVE_divisions.ColumnWidth = {0.99*0.4*scrsz(3)*table_RVE_divisions.Position(3)}; % Auto width
            table_RVE_divisions.RowName = []; % Remove row name
            table_RVE_divisions.Data=[2];
            
            Button_RVE_add = uicontrol('Parent', fig_RVE, 'Style', 'pushbutton', 'String', '+',...
                'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.05 0.45 0.15 0.08],'Callback',{@pushbutton_RVE_add_Callback});
            Button_RVE_remove = uicontrol('Parent', fig_RVE, 'Style', 'pushbutton', 'String', '-',...
                'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.075+0.15 0.45 0.15 0.08],'Callback',{@pushbutton_RVE_remove_Callback});
            Button_RVE_update = uicontrol('Parent', fig_RVE, 'Style', 'pushbutton', 'String', 'update',...
                'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1+0.3 0.45 0.15 0.08],'Callback',{@pushbutton_RVE_update_Callback});
            
            if strcmp(typeRVE,'A') || strcmp(typeRVE,'C')
                checkbox_RVE_2subs = uicontrol('Parent', fig_RVE, 'Style', 'checkbox','Units','normalized','Position',[0.05 0.35 0.5 0.08],...
                    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,...
                    'String','Calculate properties with 2 subvolumes','Value',0,'Callback',{@checkbox_RVE_sub2_Callback});
                RVEparameters.subs2=false;
            end
            
            if strcmp(typeRVE,'A')
                RVEparameters.subs4=false;  
                checkbox_RVE_4subs = uicontrol('Parent', fig_RVE, 'Style', 'checkbox','Units','normalized','Position',[0.05 0.25 0.5 0.08],...
                    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,...
                    'String','Calculate properties with 4 subvolumes','Value',0,'Callback',{@checkbox_RVE_sub4_Callback});
            end
            text_RVE_criterion = uicontrol('Parent', fig_RVE, 'Style', 'text','Units','normalized','Position',[0.6 0.7 0.35 0.05],...
                'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
                'String','RVE criterion','HorizontalAlignment','left');
            table_RVE_criterion = uitable('Parent', fig_RVE,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.6 0.6 0.35 0.1],'CellEditCallback',@cellsection_RVE_criterion_Callback, 'CellSelectionCallback',@cellsection_RVE_criterion_Callback);
            table_RVE_criterion.ColumnName = {'Std < ?%','Min. number of subvolumes'}; % Column name
            table_RVE_criterion.ColumnEditable = [true true]; % Select column editable
            table_RVE_criterion.RowName = []; % Remove row name
            table_RVE_criterion.Data=[default_threshold_std defaut_threshold_numbersubvolumes];
            RVEparameters.threshold_std=default_threshold_std;
            RVEparameters.threshold_numbersubvolumes=defaut_threshold_numbersubvolumes;
        end
        
        if strcmp(typeRVE,'C')
            text_RVE_constantlength = uicontrol('Parent', fig_RVE, 'Style', 'text','Units','normalized','Position',[0.6 0.9 0.35 0.05],...
                'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
                'String','Constant direction length','HorizontalAlignment','left');
            Popup_RVEconstantlength = uicontrol('Parent', fig_RVE,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
                'String', {' ', 'Direction 1','Direction 2', 'Direction 3'},'Units','normalized','Position', [0.6 0.8 0.35 0.1],'Callback', @popup_RVEconstantlength_Callback);
        end
        
        if strcmp(typeRVE,'D') || strcmp(typeRVE,'E')
            text_RVEdirection = uicontrol('Parent', fig_RVE, 'Style', 'text','Units','normalized','Position',[0.05 0.9 0.5 0.05],...
                'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
                'String','Choose below the direction of growth.','HorizontalAlignment','left');
            if strcmp(typeRVE,'D')
                Popup_RVEdirection = uicontrol('Parent', fig_RVE,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'enable','off',...
                    'String', {'From center'},'Units','normalized','Position', [0.05 0.8 0.5 0.1]);
            elseif strcmp(typeRVE,'E')
                Popup_RVEdirection = uicontrol('Parent', fig_RVE,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
                    'String', {' ', 'Direction 1, from x min to x max','Direction 1, from x max to x min','Direction 2, from x min to x max','Direction 2, from x max to x min','Direction 3, from x min to x max','Direction 3, from x max to x min'},'Units','normalized','Position', [0.05 0.8 0.5 0.1],'Callback', @popup_RVEdirection_Callback);
            end
            text_growthRVE = uicontrol('Parent', fig_RVE, 'Style', 'text','Units','normalized','Position',[0.05 0.7 0.9 0.05],...
                'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
                'String','Volume growth, expressed in percents of the total volume, or of the subvolume','HorizontalAlignment','left');            
            Edit_growperstep = uicontrol('Parent', fig_RVE,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
                'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.05 0.6 0.1 0.08],'Callback', @edit_RVEgrowthstep_Callback); 
            Popup_RVErelativegrowth = uicontrol('Parent', fig_RVE,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
                'String', {' ', '% of total volume','% of current subvolume','micrometers'},'Units','normalized','Position', [0.175 0.6 0.3 0.08],'Callback', @popup_RVErelativegrowth_Callback);

            text_fistvolumeRVE = uicontrol('Parent', fig_RVE, 'Style', 'text','Units','normalized','Position',[0.05 0.5 0.9 0.05],...
                'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
                'String','Size of the first subvolume, expressed in percents of the total volume, or in equivalent cubic length (um)','HorizontalAlignment','left');               
            Edit_RVEfirstvolumesize = uicontrol('Parent', fig_RVE,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
                'BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.05 0.4 0.1 0.08],'Callback', @edit_fistvolumeRVE_Callback);
            Popup_RVEfirstvolumesizeunit = uicontrol('Parent', fig_RVE,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
                'String', {' ', '% of total volume','micrometers'},'Units','normalized','Position', [0.175 0.4 0.3 0.08],'Callback', @popup_RVEfirstvolumesizeunit_Callback);
        end
        
        if strcmp(typeRVE,'B') || strcmp(typeRVE,'D')
            text_RVE_AR = uicontrol('Parent', fig_RVE, 'Style', 'text','Units','normalized','Position',[0.6 0.9 0.35 0.05],...
                'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
                'String','Subvolume Aspect ratio','HorizontalAlignment','left');
            table_RVE_AR = uitable('Parent', fig_RVE,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.6 0.8 0.35 0.1],'CellEditCallback',@cellsection_RVE_AR_Callback, 'CellSelectionCallback',@cellsection_RVE_AR_Callback);
            table_RVE_AR.ColumnName = {'Direction 1','Direction 2','Direction 3'}; % Column name
            table_RVE_AR.ColumnEditable = [true true true]; % Select column editable
            table_RVE_AR.RowName = []; % Remove row name
            table_RVE_AR.Data=[1 1 1];
            RVEparameters.Aspectratio=[1 1 1];
        end

        Button_RVE_save = uicontrol('Parent', fig_RVE, 'Style', 'pushbutton', 'String', 'Save and close',...
            'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.05 0.05 0.4 0.05],'Callback',{@pushbutton_RVE_saveclose_Callback});
        Button_RVE_cancel = uicontrol('Parent', fig_RVE, 'Style', 'pushbutton', 'String', 'Cancel and close',...
            'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.55 0.05 0.4 0.05],'Callback',{@pushbutton_RVE_cancelclose_Callback});

        function edit_fistvolumeRVE_Callback(~,~)
            val_ = str2num(Edit_RVEfirstvolumesize.String);
            no_nan = ~isnan(val_);
            is_positive=1;
            if no_nan
                positive = val_>0;
                if min(positive)==0
                    is_positive=0;
                end
            end
            if no_nan*is_positive
                set(text_fistvolumeRVE,'String','Size of the first subvolume, expressed in percents of the total volume, or in equivalent cubic length (um)','BackgroundColor',[0.9400 0.9400 0.9400]);
                RVEparameters.firstuniquevolume_size = val_;
            else
                set(text_fistvolumeRVE,'String','Error! Must be >0','BackgroundColor','r');
            end
        end           
        
        function popup_RVEfirstvolumesizeunit_Callback(~,~)
            RVEparameters.firstuniquevolume_unit = char(Popup_RVEfirstvolumesizeunit.String(Popup_RVEfirstvolumesizeunit.Value));
        end        
        
        function popup_RVEconstantlength_Callback(~,~)
            RVEparameters.Constantdirection = char(Popup_RVEconstantlength.String(Popup_RVEconstantlength.Value));
        end                
               
        function popup_RVErelativegrowth_Callback(~,~)
            RVEparameters.Growthrelativeto = char(Popup_RVErelativegrowth.String(Popup_RVErelativegrowth.Value));
        end
        
        function edit_RVEgrowthstep_Callback(~,~)
            val_ = str2num(Edit_growperstep.String);
            no_nan = ~isnan(val_);
            is_positive=1;
            if no_nan
                positive = val_>0;
                if min(positive)==0
                    is_positive=0;
                end
            end
            if no_nan*is_positive
                set(text_growthRVE,'String','Volume growth, expressed in percents of the total volume, or of the subvolume','BackgroundColor',[0.9400 0.9400 0.9400]);
                RVEparameters.Growthperstep = val_;
            else
                set(text_growthRVE,'String','Error! Must be >0','BackgroundColor','r');
            end
        end        
        
        function checkbox_RVE_sub2_Callback(~,~)
            RVEparameters.subs2 = logical(checkbox_RVE_2subs.Value);
        end
        function checkbox_RVE_sub4_Callback(~,~)
            RVEparameters.subs4 = logical(checkbox_RVE_4subs.Value);
        end
        
        function cellsection_RVE_criterion_Callback(~,~)
            val_ = table_RVE_criterion.Data;
            no_nan = ~isnan(val_);
            is_positive=1;
            if no_nan
                positive = val_>0;
                if min(positive)==0
                    is_positive=0;
                end
            end
            if no_nan*is_positive
                set(text_RVE_criterion,'String','RVE criterion','BackgroundColor',[0.9400 0.9400 0.9400]);
                RVEparameters.threshold_std = val_(1,1);
                RVEparameters.threshold_numbersubvolumes = val_(1,2);
            else
                set(text_RVE_criterion,'String','Error! Must be >0','BackgroundColor','r');
            end
        end
        
        function cellsection_RVE_AR_Callback(~,~)
            val_ = table_RVE_AR.Data;
            no_nan = ~isnan(val_);
            is_positive=1;
            if no_nan
                positive = val_>0;
                if min(positive)==0
                    is_positive=0;
                end
            end
            if no_nan*is_positive
                set(text_RVE_AR,'String','Subvolume Aspect ratio','BackgroundColor',[0.9400 0.9400 0.9400]);
                val_=val_/val_(3); % Normalize
                RVEparameters.Aspectratio = val_;
            else
                set(text_RVE_AR,'String','Error! Must be >0','BackgroundColor','r');
            end
        end
        
        function popup_RVEdirection_Callback(~,~)
            RVEparameters.Growthdirection = char(Popup_RVEdirection.String(Popup_RVEdirection.Value));
        end
        
        function cellsection_RVEdivision_Callback(~,~)
            val_ = table_RVE_divisions.Data;
            no_nan = ~isnan(val_);
            is_positive=1;
            if no_nan
                positive = val_>0;
                if min(positive)==0
                    is_positive=0;
                end
            end
            if no_nan*is_positive
                RVEparameters.divisions = val_';
            end
        end
        
        function pushbutton_RVE_add_Callback(~,~)
            val_ = table_RVE_divisions.Data;
            val_(val_==1)=[]; % Un-necessary
            val_ = unique(val_);
            val_=sort(val_);
            val_= [val_; val_(end)+1];
            table_RVE_divisions.Data = val_;
            RVEparameters.divisions = val_';
        end
        function pushbutton_RVE_remove_Callback(~,~)
            val_ = table_RVE_divisions.Data;
            val_(val_==1)=[]; % Un-necessary
            val_ = unique(val_);
            val_=sort(val_);
            val_(end)=[];
            if isempty(val_)
                val_=2;
            end
            table_RVE_divisions.Data = val_;
            RVEparameters.divisions = val_';
        end
        function pushbutton_RVE_update_Callback(~,~)
            val_ = table_RVE_divisions.Data;
            val_(val_==1)=[]; % Un-necessary
            val_ = unique(val_);
            val_=sort(val_);
            table_RVE_divisions.Data = val_;
            RVEparameters.divisions = val_';
        end

        function pushbutton_RVE_saveclose_Callback(~,~) % Update
            if strcmp(typeRVE,'A')
                RVE_name='Independant subvolumes of same size + keep initial aspect ratio (A)';
                RVE_savename = 'SubInitialAR';
            elseif strcmp(typeRVE,'B')
                RVE_name='Independant subvolumes of same size + user-defined aspect ratio (B)';
                %RVE_savename = ['SubCustomAR_' num2str(RVEparameters.Aspectratio,'%1.3f\t')];
                RVE_savename = ['SubCustomAR_' num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
            elseif strcmp(typeRVE,'C')
                RVE_name='Independant subvolumes of same size + constant length (C)';
                tmp = num2str(RVEparameters.Constantdirection);
                RVE_savename = ['SubConstL_' tmp(end)];
            elseif strcmp(typeRVE,'D')
                RVE_name='One subvolume + growing from volume center (D)';
                %RVE_savename = ['GrowthCenter_' num2str(RVEparameters.Aspectratio,'%1.3f\t')];
                RVE_savename = ['GrowthCenter_' num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
            elseif strcmp(typeRVE,'E')
                RVE_name='One subvolume + growing from volume edge (E)';
                if strcmp(RVEparameters.Growthdirection, 'Direction 1, from x min to x max')
                    RVE_savename = 'GrowthEdge_D1_minmax';
                elseif strcmp(RVEparameters.Growthdirection,'Direction 1, from x max to x min')
                    RVE_savename = 'GrowthEdge_D1_maxmin';
                elseif strcmp(RVEparameters.Growthdirection,'Direction 2, from x min to x max')
                    RVE_savename = 'GrowthEdge_D2_minmax';                    
                elseif strcmp(RVEparameters.Growthdirection,'Direction 2, from x max to x min')
                    RVE_savename = 'GrowthEdge_D2_maxmin';
                elseif strcmp(RVEparameters.Growthdirection,'Direction 3, from x min to x max')
                    RVE_savename = 'GrowthEdge_D3_minmax';
                elseif strcmp(RVEparameters.Growthdirection,'Direction 3, from x max to x min')
                    RVE_savename = 'GrowthEdge_D3_maxmin';                    
                end
            end
            [ RVE_savename ] = function_remove_emptyandspecialcharacter_string(RVE_savename);
            PROPERTY.(p).number_RVE = PROPERTY.(p).number_RVE+1;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).name = RVE_name;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).savename = RVE_savename;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).type = RVEparameters.type;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).divisions = RVEparameters.divisions;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).subs2 = RVEparameters.subs2;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).subs4 = RVEparameters.subs4;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).Aspectratio = RVEparameters.Aspectratio; 
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).Constantdirection = RVEparameters.Constantdirection;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).Growthdirection =  RVEparameters.Growthdirection;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).Growthperstep = RVEparameters.Growthperstep;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).Growthrelativeto = RVEparameters.Growthrelativeto;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).threshold_std = RVEparameters.threshold_std;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).threshold_numbersubvolumes = RVEparameters.threshold_numbersubvolumes;     
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).firstuniquevolume_size = RVEparameters.firstuniquevolume_size;
            PROPERTY.(p).RVE(PROPERTY.(p).number_RVE).firstuniquevolume_unit = RVEparameters.firstuniquevolume_unit;
            
            column_1={RVEparameters.type;};
            column_2={num2str(RVEparameters.divisions)};
            column_3={RVEparameters.subs2};
            column_4={RVEparameters.subs4};
            column_5={num2str(RVEparameters.Aspectratio)};
            column_6={RVEparameters.Constantdirection};
            column_7={RVEparameters.Growthdirection};
            column_8={[num2str(RVEparameters.Growthperstep) ' ' RVEparameters.Growthrelativeto]};
            column_9={[num2str(RVEparameters.firstuniquevolume_size) ' ' RVEparameters.firstuniquevolume_unit]};
            
            h.Data= [h.Data ;[column_1 column_2 column_3 column_4 column_5 column_6 column_7 column_8 column_9]];

            close(fig_RVE);
        end
        function pushbutton_RVE_cancelclose_Callback(~,~)
            close(fig_RVE);
        end
    end




end