function microstructure_generation_ellipsoids_GUI
% Graphic user interface of the microstructure generation toolbox module 
% GUI by Prehit Patel
% Generation algorithm by Francois Usseglio-Viretta (NREL)
%    see: function_generate_ellipsoid_microstructure.m

close all
clear all
clc

%% GUI DISPLAY OPTIONS
font_name_GUI ='Arial'; % Font
font_type_plot = 'Arial'; % plot font
% Font size
font_size_plot_label_large = 14;
font_size_plot_label_small = 8;
font_size_plot_title = 16;
font_size_small_GUI =8;
font_size_medium_GUI =10;
font_size_large_GUI =14;
% Background and text color
background_description_tab = [0 0.5 0];
ForegroundColor_description_tab = [1 1 1];
error_message_red = [1 0 0]; % Error that prevent user to gon
error_message_orange = [1 0.64 0]; % Warning message

%% ERROR STATUTS

% Tab 1 error
error_tab1_voxel_scaling = false; error_tab1_num_voxel = false; error_tab1_num_phase = false; error_tab1_code_check = false; error_tab1_num_slice = false; error_tab1_table3_column1 = false;
error_tab1_table3_porosity = false;

% Tab 2 error
error_tab2_para_num_slices_num_dia = false; error_tab2_table1_a_sum = false; error_tab2_table1_b_column1 = false;
error_tab2_table2_a_sum = false; error_tab2_table2_b_column1 = false; error_tab2_table3_a_sum = false; error_tab2_table3_b_column1 = false;

% Tab 3 error 
error_tab3_para_num_slices_num_dia = false; error_tab2_table4_a_sum = false; error_tab2_table4_b_column1 = false; error_tab2_table4_c_minmax = false; error_tab2_table4_d_minmax_value_match = false;
error_tab2_table5_a_sum = false; error_tab2_table5_b_column1 = false; error_tab2_table5_c_minmax = false; error_tab2_table5_d_minmax_value_match = false;
error_tab2_table6_a_sum = false; error_tab2_table6_b_column1 = false; error_tab2_table6_c_minmax = false; error_tab2_table6_d_minmax_value_match = false;

%Tab 3 error
error_tab3_num_pass = false; error_tab3_table3_fillratio = false;

%Tab 5 error
error_tab5_num_of_run = false;
%% PRESELECTED OPTIONS AND VARIABLES
scrsz = get(0,'ScreenSize');

% % NOTE: to prepare save path in both Windows and Unix system

if ispc
    folder_separation = '\';
else
    folder_separation = '/';
end

% % NOTE: put here all preallocated values so we do not have to change

default_voxelsize = 50;
default_scalingfactor = 2;
Number_of_Voxel = [100;100;150];
current_phase_sizeorientation_tab2 =[];
statuts_size_orientation_tab2 =[];
current_phase_sizeorientation_tab3 =[];
statuts_size_orientation_tab3 =[];
domain_size=[];
phase =[];
Maximum_interpenetration=[];
Minimum_particle_volume_conservated=[];
voxel_size_nm = [];
Save_folder = 'None';
OPTIONS.fontname = 'Arial';
OPTIONS.axe_fontsize=10;
OPTIONS.legend_fontsize=8;
OPTIONS.title_fontsize=14;
OPTIONS.figure_Linewidth=2;
OPTIONS.grid = true;
OPTIONS.minorgrid = true;
OPTIONS.save_fig = true;
OPTIONS.savefig_infig = true;
OPTIONS.savefig_informat = {'png'};
OPTIONS.closefigureaftercreation = false;
OPTIONS.save_video=true; % Creating video takes times.
OPTIONS.video_quality=100;
OPTIONS.video_framerate=15;
OPTIONS.video_colormap_phase=copper;

%% MAIN FIGURE
main_figure = figure; % Create figure
main_figure.Name= 'Microstructure generation'; % Set name
main_figure.NumberTitle='off'; % Remove number from title name
main_figure.Color='white'; % Background colour
main_figure.MenuBar='none'; % Remove menubar and toolbar
main_figure.ToolBar='none';
main_figure.Units='normalized'; % Set unit
lfigure=[0.2 0.2 0.7 0.6];
main_figure.Position=lfigure; % Set Position

%Tab group
table_group_1 = uitabgroup('Parent', main_figure,'TabLocation','Left');
tab_phase_volumefractions = uitab('Parent', table_group_1,'Title','Phase and volume Fractions','Units','normalized');
tab_particlesize = uitab('parent',table_group_1,'Title','Particle size','Units','normalized');
tab_orientation = uitab('parent',table_group_1,'Title','Particle orientation','Units','normalized');
tab_overlapping = uitab('parent',table_group_1,'Title','Particle overlapping','Units','normalized');
tab_pp = uitab('parent',table_group_1,'Title','Post-processing and Save options','Units','normalized');
tab_run = uitab('parent',table_group_1,'Title','Run code','Units','normalized');
tab_Generated_micro = uitab('parent',table_group_1,'Title','Generated microstructures','Units','normalized');

% Tab colors
tab_phase_volumefractions.ForegroundColor  ='b';
tab_particlesize.ForegroundColor  ='b';
tab_orientation.ForegroundColor = 'b';
tab_overlapping.ForegroundColor  ='b';
tab_pp.ForegroundColor  ='b';
tab_run.ForegroundColor  ='r';
tab_Generated_micro.ForegroundColor = 'b';

% Preselected tab
table_group_1.SelectedTab = tab_phase_volumefractions;

% Tabgroup width estimation
tabgroup_width=0.16; % % NOTE: i couldnt find a command that provides this value...

%% GUI OBJECTS
%% Phase and Volume Fraction Calculation tab 1

% Title Phase and Volume Fraction tab
Title_tab1 = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','Phase and Volume Fraction Settings','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]); % Heading

% Sub heading Voxel Size and Scaling Factor
Text_voxelsize_1 = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','Voxel Size (nm)','Units','normalized','Position',[0,0.92,0.1,0.03],'FontName',font_name_GUI,'FontSize',font_size_small_GUI);
Scaling_Factor_1 = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','Scaling Factor','Units','normalized','Position',[0,0.87,0.1,0.03],'FontName',font_name_GUI,'FontSize',font_size_small_GUI);
Text_voxelsize_2 = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','This is the voxel size(nm) after rescaling','Units','normalized','Position',[0.15,0.92,0.2,0.028],'FontName',font_name_GUI,'FontSize',font_size_small_GUI);
Scaling_Factor_2 = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','Microstructure generation is performed on a coarse domain to save CPU time. Then, the obtained microstructure is upscaled.','Units','normalized','Position',[0.16,0.85,0.35,0.05],'HorizontalAlignment','left','FontName',font_name_GUI,'FontSize',font_size_small_GUI);

% Input Value box (Voxel Size and Scaling Factor)
Edit_voxelsize = uicontrol('Parent',tab_phase_volumefractions,'Style','edit','String',num2str(default_voxelsize),'Units','normalized','Position',[0.1,0.92,0.05,0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'Callback',{@edit_voxelsize_scalingfactor_Callback});
Edit_Scalingfactor = uicontrol('Parent',tab_phase_volumefractions,'Style','edit','String',num2str(default_scalingfactor),'Units','normalized','Position',[0.1,0.87,0.05,0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'Callback',{@edit_voxelsize_scalingfactor_Callback});

% Table Input Direction
ltable_direction = [0.01 0.7 0.5 0.14];
table_direction = uitable('Parent', tab_phase_volumefractions,'Units','normalized','Position',ltable_direction,'CellEditCallback',@cellsection_table_direction_Callback,'FontName',font_name_GUI,'FontSize',font_size_small_GUI);
table_direction.ColumnName = {'Direction','Number of voxel',' After re-scaling','length (um)'};
Direction(1).name='Direction 1';
Direction(2).name='Direction 2';
Direction(3).name='Direction 3';
input_data = [{Direction(1).name; Direction(2).name; Direction(3).name} num2cell(Number_of_Voxel) num2cell([200;200;300]) num2cell([10;10;15])];
table_direction.Data=input_data;
table_direction.ColumnEditable = [false true false false];
table_direction.RowName = []; % Remove row name
ltable=scrsz(3)*lfigure(3)*(1-tabgroup_width)*ltable_direction(3);
table_direction.ColumnWidth = {ltable*0.3, ltable*0.25, ltable*0.25, ltable*0.2}; % Fit table length

Text_Direction3 = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','The convention used in this generation algorithm is that direction 3 corresponds to the volume ‘thickness’, for which different volume fractions, particle size and particle rotation can be set to generate a graded microstructure if desired.','Units','normalized','Position',[0.01,0.63,0.5,0.05],'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'HorizontalAlignment','left');

%Sub heading Number of Phase
Num_phas = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','Number of Solid Phase','Units','normalized','Position',[0,0.55,0.12,0.024],'FontName',font_name_GUI,'FontSize',font_size_small_GUI);

% Input Value box (Number of Phase)
Edit_Num_phas = uicontrol('Parent',tab_phase_volumefractions,'Style','edit','String',[],'Units','normalized','Position',[0.12,0.55,0.05,0.028],...
    'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'Callback',@edit_Num_phas);

% Table 2 created based on input given in Number of phase
table_phas_num = uitable('Parent', tab_phase_volumefractions,'Units','normalized','Position',[0.01 0.38 0.25 0.14],'Visible','off','Enable','off','CellEditCallback',@cellsection_table_phas_num_Callback,'FontName',font_name_GUI,'FontSize',font_size_small_GUI);
table_phas_num.ColumnName = {'Phase Number','Code','Name'};
table_phas_num.RowName = [];
table_phas_num.ColumnWidth = {ltable*0.2, ltable*0.15, ltable*0.15};
table_phas_num.ColumnEditable = [false true true];

A{1,1} = '“Phase number” is a hard coded value used for numbering solid phases.';
A{1,2} = '“Code” is the user-defined positive integer value that will be assigned to each voxel that belong to the corresponding phase. The value 0 is reserved for the background, i.e. the pore phase.';
A{1,3} = '';
A{1,4} = '“Name” is the string used to label graphs. It is purely an esthetic parameter.';
multiple_line_string_tab1 = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',A{1,1},A{1,2},A{1,3},A{1,4});
Text_saveconfirmation_tab1 = uicontrol('Parent',tab_phase_volumefractions,'Style', 'text','String',multiple_line_string_tab1,'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Visible','on'...
    ,'HorizontalAlignment','left','Units','normalized','Position',[0.01 0.2 0.5 0.15],'Visible','off');

%Split screen
tab1 = annotation(tab_phase_volumefractions,'line',[0.52 0.52],[0.06 0.95]);

% Sub Heading Number of Slices
Num_Slices = uicontrol('Parent',tab_phase_volumefractions,'Style','text','String','In the table below, you can specify volume fractions at different normalized positions along the microstructure thickness (i.e., along direction 3). Volumes fractions must be between 0 and 1, and set so that their sum is below 1 as the remaining volume is reserved for the background (i.e., the pore) phase. You do not need to specify values for every positions as volume fractions for intermediate positions will be linearly interpolated. At least the first slice (normalized position = 0) and the last slice (normalized position = 1) must have defined volume fractions.'...
    ,'Units','normalized','HorizontalAlignment','left','Position',[0.55,0.77,0.38,0.15],'Visible','off','Enable','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI);

% Input value box (Number of Slices)
Edit_Num_slices = uicontrol('Parent',tab_phase_volumefractions,'Style','edit','String',[],'Units','normalized','Position',[0.935,0.87,0.05,0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'Visible','off','Enable','off','Callback',{@edit_Num_slice});

% Table 3 created based on input given in Number of slice
table_position = uitable('Parent', tab_phase_volumefractions,'Units','normalized','Position',[0.55 0.63 0.43 0.13],'Visible','off','Enable','off','CellEditCallBack',@cellsection_table_position_Callback,'FontName',font_name_GUI,'FontSize',font_size_small_GUI);
table_position.ColumnWidth = {'auto'}; % Auto width
table_position.RowName = []; % Remove row name

Volumefraction_Plot_View_button = uicontrol( 'Parent', tab_phase_volumefractions,'Style','pushbutton','String','Click to plot volume fraction','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.55,0.6,0.43,0.03],'Visible','off','Callback',{@Callback_plot});
str = sprintf('Click to plot volume fraction in direction 3…');
Volumefraction_Plot_View_button.TooltipString = str;
volume_fraction_plot = axes('Parent', tab_phase_volumefractions,'Visible','off','Units','normalized','Position', [0.58 0.14 0.4 0.41]);
tab1_plot = axes('Visible','off');
tab1_plot_legend = legend('Visible','off');

Save_Input_button_tab1 = uicontrol('Parent', tab_phase_volumefractions,'Style','pushbutton','String','Click to save all the input parametes','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.55,0.02,0.43,0.035],'Visible','off','Callback',@Callback_save_input_data,'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);

Test_error_tab1 = uicontrol('Parent', tab_phase_volumefractions, 'Style', 'text','Units','normalized','Position',[0.01 0.01 0.98 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Error message','Visible','off');

%% Particle size and orienation tab 2

% Title Particle size and orienation tab
Title_tab2 = uicontrol('Parent',tab_particlesize,'Style','text','String','Particle Size Settings','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]);

% Popup menu
Select_Phase_tab2 = uicontrol('Parent',tab_particlesize,'Style','text','String','Select Phase','Units','normalized','Position',[0.01,0.91,0.08,0.028]);

Configured_tab2 = uicontrol('Parent',tab_particlesize,'Style','text','String',' ','Units','normalized','Position',[0.22,0.915,0.2,0.025],'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'BackgroundColor','w','ForegroundColor','r');

A{1,1} = 'Step 1 - Select the phase to setup';
A{1,2} = 'Step 2 - Edit Number of Slices, Diameter, Elongation, and Orientation';
A{1,3} = 'Step 3 - Edit all the 3 tables below (diameter in number of voxel)';
A{1,4} = 'Step 4 - Click the Plot Area buttons to verify';
A{1,5} = 'Step 5 - Click to save current phase parameters';
A{1,6} = 'Step 6 - Repeat untill all the phases have been configured, then move to the next step';
multiple_line_string = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',A{1,1},A{1,2},A{1,3},A{1,4},A{1,5},A{1,6});
Text_saveconfirmation = uicontrol('Parent',tab_particlesize,'Style', 'text','String',multiple_line_string,'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Visible','on',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position',[0.56 0.74 0.45 0.13]);

popupmenu_phase_selection_tab2 = uicontrol('Parent',tab_particlesize,'Style','popupmenu','Units','normalized','Position',[0.1,0.92,0.1,0.028],'Visible','off','Enable','off', 'CallBack', @Callback_popupmenu);

parameter_table_position_tab2 = [0.01 0.74 0.54 0.13];
parameter_table_tab2 = uitable('Parent',tab_particlesize,'Units','normalized','Position',parameter_table_position_tab2,'Visible','off','Enable','off','CellEditCallBack',@parameters_callback);

Diameter_Dz_along_3rd_axis = uicontrol('Parent',tab_particlesize,'Style','text','String','Diameter Dx along 3rd axis','Units','normalized','Position',[0.01 0.55 0.25 0.15],'Visible','off');
Diameter_Dz_3rdaxis_table = uitable('Parent',tab_particlesize,'Units','normalized','Position',[0.01 0.53 0.25 0.15],'Visible','off','Enable','off','CellEditCallBack',@Dia_Dz_3rdaxis_table);
Diameter_Dz_3rdaxis_table.RowName = [];

Diameter_elongation_Dx_Dy_3rd_axis = uicontrol('Parent',tab_particlesize,'Style','text','String','Diameter Elongation Dx/Dy along 3rd axis','Units','normalized','Position',[0.01 0.35 0.25 0.15],'Visible','off');
Diameter_Elongation_Dx_Dy_3rd_axis_table = uitable('Parent',tab_particlesize,'Units','normalized','Position',[0.01 0.33 0.25 0.15],'Visible','off','Enable','off','CellEditCallBack',@Dia_elongation_Dx_Dy_axis_table);
Diameter_Elongation_Dx_Dy_3rd_axis_table.RowName = [];

Diameter_elongation_Dx_Dz_3rd_axis = uicontrol('Parent',tab_particlesize,'Style','text','String','Diameter Elongation Dx/Dz along 3rd axis','Units','normalized','Position',[0.01 0.15 0.25 0.15],'Visible','off');
Diameter_Elongation_Dx_Dz_3rd_axis_table = uitable('Parent',tab_particlesize,'Units','normalized','Position',[0.01 0.13 0.25 0.15],'Visible','off','Enable','off','CellEditCallBack',@Dia_elongation_Dx_Dz_axis_table);
Diameter_Elongation_Dx_Dz_3rd_axis_table.RowName = [];

plot_area_view_button_diameter = uicontrol('Parent',tab_particlesize,'Style','pushbutton','String','Plot Area','Units','normalized','Position',[0.32,0.681,0.06,0.028],'Visible','off','Enable','off','CallBack',@Plot_Diameter_area_callback);
str = sprintf('Click to plot area in direction 3….');
plot_area_view_button_diameter.TooltipString = str;
plot_area_view_button_diameter_plot_1 = axes('Parent', tab_particlesize,'Visible','off','Units','normalized','Position', [0.29 0.55 0.12 0.12]);
plot_area_view_button_diameter_plot_2 = axes('Parent', tab_particlesize,'Visible','off','Units','normalized','Position', [0.29 0.35 0.12 0.12]);
plot_area_view_button_diameter_plot_3 = axes('Parent', tab_particlesize,'Visible','off','Units','normalized','Position', [0.29 0.14 0.12 0.12]);

A1 = axes('Visible','off');
A2 = axes('Visible','off');
A3 = axes('Visible','off');
a1= legend('Visible','off');
a2= legend('Visible','off');
a3= legend('Visible','off');

visual_check_tab2 = uicontrol('Parent',tab_particlesize,'Style','text','String','Particle Diameter Visualization Check','Units','normalized','Position',[0.55 0.55 0.25 0.15],'Visible','off');
visual_check_table_tab2 = uitable('Parent',tab_particlesize,'Units','normalized','Position',[0.55 0.55 0.242 0.125],'Visible','off','Enable','off','CellEditCallBack',@visual_plot_table_tab2);
visual_check_table_tab2.RowName = [];
visual_check_table_tab2.ColumnName = {'Direction', 'Diameter'}; 
visual_plot_particle_button_tab2 = uicontrol('Parent',tab_particlesize,'Style','pushbutton','String','Generate particle','Units','normalized','Position',[0.85,0.6,0.08,0.035],'Visible','off','Enable','off','CallBack',@Visual_Plot_button_callback_tab2);
axes1_tab2 = axes('Parent',tab_particlesize,'Units','normalized','Position', [0.58 0.08 0.4 0.45],'Visible','off');

Save_Input_button_tab2 = uicontrol( 'Parent', tab_particlesize,'Style','pushbutton','String','Save Current Phase','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.79,0.91,0.2,0.035],'Visible','off','Callback',@Callback_save_input_data_tab2);
str = sprintf('Click to save all the input parametes...');
Save_Input_button_tab2.TooltipString = str;

Cancel_button_tab2 = uicontrol( 'Parent', tab_particlesize,'Style','pushbutton','String','Cancel Current Phase','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.58,0.91,0.2,0.035],'Visible','off','Callback',@Callback_cancel_tab2);
str = sprintf('Click to cancel current phase...');
Cancel_button_tab2.TooltipString = str;

Test_error_tab2 = uicontrol('Parent',tab_particlesize, 'Style', 'text','Units','normalized','Position',[0.01 0.01 0.98 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Error message','Visible','off');

%% Tab 3

% Title Particle size and orienation tab
Title_tab2 = uicontrol('Parent',tab_orientation,'Style','text','String','Particle Orientation Settings','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]);

% Popup menu
Select_Phase = uicontrol('Parent',tab_orientation,'Style','text','String','Select Phase','Units','normalized','Position',[0.01,0.91,0.08,0.028]);

Configured_tab3 = uicontrol('Parent',tab_orientation,'Style','text','String',' ','Units','normalized','Position',[0.22,0.915,0.2,0.025],'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'BackgroundColor','w','ForegroundColor','r');

A{1,1} = 'Step 1 - Select the phase to setup';
A{1,2} = 'Step 2 - Edit Number of Slices, Diameter, Elongation, and Orientation';
A{1,3} = 'Step 3 - Edit all the 6 tables below (diameter in number of voxel and orientation in degrees)';
A{1,4} = 'Step 4 - Click to both the Plot Area buttons to verify';
A{1,5} = 'Step 5 - Click to save current phase parameters';
A{1,6} = 'Step 6 - Repeat untill all the phases have been configured, then move to the next step';
multiple_line_string_tab3 = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',A{1,1},A{1,2},A{1,3},A{1,4},A{1,5},A{1,6});
Text_saveconfirmation_tab3 = uicontrol('Parent',tab_orientation,'Style', 'text','String',multiple_line_string_tab3,'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Visible','on',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position',[0.56 0.74 0.45 0.13]);

popupmenu_phase_selection_tab3 = uicontrol('Parent',tab_orientation,'Style','popupmenu','Units','normalized','Position',[0.1,0.92,0.1,0.028],'Visible','off','Enable','off', 'CallBack', @Callback_popupmenu_tab3);

parameter_table_position_tab3 = [0.01 0.74 0.54 0.13];
parameter_table_tab3 = uitable('Parent',tab_orientation,'Units','normalized','Position',parameter_table_position_tab3,'Visible','off','Enable','off','CellEditCallBack',@parameters_callback_tab3);

Orientation_normal_Dx_3rd_axis = uicontrol('Parent',tab_orientation,'Style','text','String','Orientation Normal to Ox Axis','Units','normalized','Position',[0.08 0.55 0.25 0.15],'Visible','off');
Orientation_normal_Dx_3rd_axis_table = uitable('Parent',tab_orientation,'Units','normalized','Position',[0.01 0.53 0.38 0.15],'Visible','off','Enable','off','CellEditCallBack',@Orientation_Dx_axis_callback);
Orientation_normal_Dx_3rd_axis_table.RowName = [];

Orientation_normal_Dy_3rd_axis = uicontrol('Parent',tab_orientation,'Style','text','String','Orientation Normal to Oy Axis','Units','normalized','Position',[0.08 0.35 0.25 0.15],'Visible','off');
Orientation_normal_Dy_3rd_axis_table = uitable('Parent',tab_orientation,'Units','normalized','Position',[0.01 0.33 0.38 0.15],'Visible','off','Enable','off','CellEditCallBack',@Orientation_Dy_axis_callback);
Orientation_normal_Dy_3rd_axis_table.RowName = [];

Orientation_normal_Dz_3rd_axis = uicontrol('Parent',tab_orientation,'Style','text','String','Orientation Normal to Oz Axis','Units','normalized','Position',[0.08 0.15 0.25 0.15],'Visible','off');
Orientation_normal_Dz_3rd_axis_table = uitable('Parent',tab_orientation,'Units','normalized','Position',[0.01 0.13 0.38 0.15],'Visible','off','Enable','off','CellEditCallBack',@Orientation_Dz_axis_callback);
Orientation_normal_Dz_3rd_axis_table.RowName = [];

plot_area_view_button_orientation = uicontrol('Parent',tab_orientation,'Style','pushbutton','String','Plot Area','Units','normalized','Position',[0.45,0.681,0.06,0.028],'Visible','off','Enable','off','CallBack',@Plot_orientation_area_callback);
str = sprintf('Click to plot area in direction 3….');
plot_area_view_button_orientation.TooltipString = str;
plot_area_view_button_orientation_plot_1 = axes('Parent', tab_orientation,'Visible','off','Units','normalized','Position', [0.42 0.55 0.12 0.12]);
plot_area_view_button_orientation_plot_2 = axes('Parent', tab_orientation,'Visible','off','Units','normalized','Position', [0.42 0.35 0.12 0.12]);
plot_area_view_button_orientation_plot_3 = axes('Parent', tab_orientation,'Visible','off','Units','normalized','Position', [0.42 0.14 0.12 0.12]);

A4 = axes('Visible','off');
A5 = axes('Visible','off');
A6 = axes('Visible','off');
a4= legend('Visible','off');
a5= legend('Visible','off');
a6= legend('Visible','off');

visual_check_tab3 = uicontrol('Parent',tab_orientation,'Style','text','String','Particle Orientation Visualization Check','Units','normalized','Position',[0.55 0.55 0.25 0.15],'Visible','off');
visual_check_table_tab3 = uitable('Parent',tab_orientation,'Units','normalized','Position',[0.56 0.55 0.242 0.125],'Visible','off','Enable','off','CellEditCallBack',@visual_plot_table_tab3);
visual_check_table_tab3.RowName = [];
visual_check_table_tab3.ColumnName = {'Direction', 'Degree'}; 
visual_plot_particle_button_tab3 = uicontrol('Parent',tab_orientation,'Style','pushbutton','String','Generate particle','Units','normalized','Position',[0.85,0.6,0.08,0.035],'Visible','off','Enable','off','CallBack',@Visual_Plot_button_callback_tab3);
axes2_tab3 = axes('Parent',tab_orientation,'Units','normalized','Position', [0.58 0.08 0.4 0.45],'Visible','off');


Save_Input_button_tab3 = uicontrol( 'Parent', tab_orientation,'Style','pushbutton','String','Save Current Phase','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.79,0.91,0.2,0.035],'Visible','off','Callback',@Callback_save_input_data_tab3);
str = sprintf('Click to save all the input parametes...');
Save_Input_button_tab3.TooltipString = str;

Cancel_button_tab3 = uicontrol( 'Parent', tab_orientation,'Style','pushbutton','String','Cancel Current Phase','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.58,0.91,0.2,0.035],'Visible','off','Callback',@Callback_cancel_tab3);
str = sprintf('Click to cancel current phase...');
Cancel_button_tab3.TooltipString = str;

Test_error_tab3 = uicontrol('Parent',tab_orientation, 'Style', 'text','Units','normalized','Position',[0.01 0.01 0.98 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Error message','Visible','off');

%% Particle Overlapping tab 4
Title_tab3 = uicontrol('Parent',tab_overlapping,'Style','text','String','Particle Overlapping and Fill Ratio Settings','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]);

Phase_overlapping_table = uicontrol('Parent',tab_overlapping,'Style','text','String','Phase Overlapping Data','Units','normalized','Position',[0.01 0.9 0.3 0.03],'Visible','off');
Maximum_overlapping_table_position = [0.01 0.75 0.3 0.15];
Table_maximum_overlapping = uitable('Parent',tab_overlapping,'Units','normalized','Position',Maximum_overlapping_table_position,'Visible','off','Enable','off','CellEditCallBack',@Table_maximum_overlapping_callback);

Table_min_par_vol = uicontrol('Parent',tab_overlapping,'Style','text','String','Normalized Minimum Particle Volume','Units','normalized','Position',[0.01 0.68 0.3 0.03],'Visible','off');
Min_par_volume_table_position = [0.01 0.53 0.3 0.15];
Table_minimum_particle_volume_conservated = uitable('Parent',tab_overlapping,'Units','normalized','Position',Min_par_volume_table_position,'Visible','off','Enable','off','CellEditCallBack',@Table_maximum_overlapping_callback);
Table_minimum_particle_volume_conservated.ColumnEditable = true;

% Fill ratio table
Number_of_generation_pass = uicontrol('Parent',tab_overlapping,'Style','text','String','Number of generation pass','Units','normalized','Position',[0.01,0.45,0.15,0.028],'Visible','off');
Edit_fillratio = uicontrol('Parent',tab_overlapping,'Style','edit','String',[],'Units','normalized','Position',[0.15,0.45,0.05,0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'Callback',@edit_fillratio,'Visible','off','Enable','off');

fill_ratio = uicontrol('Parent',tab_overlapping,'Style','text','String','Fill Ratio Data','Units','normalized','Position',[0.01 0.38 0.3 0.03],'Visible','off');
fillratio_table_position = [0.01 0.23 0.3 0.15];
fill_ratio_table = uitable('Parent',tab_overlapping,'Units','normalized','Position',fillratio_table_position,'Visible','off','Enable','off','CellEditCallBack',@fillratio);

Save_Input_button_tab4 = uicontrol( 'Parent', tab_overlapping,'Style','pushbutton','String','Click to save all the input parametes','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.01,0.15,0.3,0.03],'Visible','off','Callback',@Callback_save_tab2);

% overlapping_imgData = imread('overlapping.jpg');
% overlapping_image = uicontrol('Parent',tab_overlapping,'CData',overlapping_imgData,'Units','normalized','Position',[0.6 0.75 0.3 0.185],'Visible','off');

% min_particle_volume_imgData = imread('overlapping.jpg');
% min_particle_volume_image = uicontrol('Parent',tab_overlapping,'CData',min_particle_volume_imgData,'Units','normalized','Position',[0.6 0.51 0.3 0.185],'Visible','off');


%fill_ration_imgData = imread('fill_ratio.jpg');
% fill_image = uicontrol('Parent',tab_overlapping,'CData',fill_ration_imgData,'Units','normalized','Position',[0.6 0.20 0.35 0.22],'Visible','off');

%% PostProcessing and Save Option Tab 5

Title_tab4 = uicontrol('Parent',tab_pp,'Style','text','String','Post-Processing and Save Options','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]);

Display_Save_Option_text = uicontrol('Parent',tab_pp,'Style','text','String','Display and save options','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'FontWeight','bold',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.01 0.9 0.15 0.04]);

Button_save_folder = uicontrol('Parent',tab_pp,'Style', 'pushbutton','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'String','Select save folder',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.01 0.85 0.15 0.04],'Callback',{@button_select_savefolder_Callback});

Text_savefolder = uicontrol('Parent',tab_pp,'Style', 'text','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'String',Save_folder,...
    'HorizontalAlignment','left','BackgroundColor','w','Units','normalized','Position', [0.16 0.85 0.345 0.04]);

Select_fontname = uicontrol('Parent',tab_pp,'Style','text','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'String','Select Fontname',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.01 0.8 0.1 0.02]);

Popup_displayoptions_fontname = uicontrol('Parent', tab_pp,'Style', 'popup','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,...
    'String', listfonts,'Value',3,'Units','normalized','Position', [0.11 0.8 0.1 0.025],'CallBack',@Font_name_CallBack);

Font_Grid_options = uitable('Parent',tab_pp,'Units','normalized','Position',[0.01 0.68 0.5 0.07],'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'CellEditCallBack',@Font_Grid_Callback);
Font_Grid_options.ColumnName = {'Axe font size','Legend font size','Title font size','Linewidth','Grid','Minor Grid'};
Font_Grid_options.ColumnWidth = {100,100,100,90,90,90};
Font_Grid_options.RowName = [];
Font_Grid_options.Data = [{OPTIONS.axe_fontsize} {OPTIONS.legend_fontsize} {OPTIONS.title_fontsize} {OPTIONS.figure_Linewidth} {OPTIONS.grid} {OPTIONS.minorgrid}];
Font_Grid_options_format =cell(1,6);
Font_Grid_options_format(1,:) ={'numeric'};
Font_Grid_options_format(1,5:6) ={'logical'};
Font_Grid_options.ColumnFormat = Font_Grid_options_format;
Font_Grid_options.ColumnEditable = [true true true true true true];

% Save_tif = uicontrol('Parent',tab_pp,'Style','text','String','Save .tif',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.7 0.85 0.15 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% Save_tif_checkbox = uicontrol('Parent',tab_pp,'Style','checkbox',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.85 0.85 0.05 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% 
% Save_tif_rescaling = uicontrol('Parent',tab_pp,'Style','text','String','Save .tif (re-scale)',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.7 0.8 0.15 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% Save_tif_rescaling_checkbox = uicontrol('Parent',tab_pp,'Style','checkbox',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.85 0.8 0.05 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% 
% Save_para_mat = uicontrol('Parent',tab_pp,'Style','text','String','Save Parameters (.mat)',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.7 0.75 0.15 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% Save_para_mat_checkbox = uicontrol('Parent',tab_pp,'Style','checkbox',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.85 0.75 0.05 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% 
% Save_para = uicontrol('Parent',tab_pp,'Style','text','String','Save Parameters ',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.7 0.7 0.8 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
% Save_para_checkbox = uicontrol('Parent',tab_pp,'Style','checkbox',...
%     'HorizontalAlignment','left','Units','normalized','Position', [0.85 0.7 0.05 0.025],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);

%Split screen
tab4 = annotation(tab_pp,'line',[0 1],[0.6 0.6]);

% Post Processing
Post_Processing_text = uicontrol('Parent',tab_pp,'Style','text','String','Post - Processing','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'FontWeight','bold',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.01 0.5 0.15 0.04]);

Post_processing_table = uitable('Parent',tab_pp,'Units','normalized','Position',[0.01 0.25 0.25 0.20],'FontName',font_name_GUI,'FontSize',font_size_medium_GUI);
Post_processing_table.ColumnName = {'Options','To Do'};
Post_processing_table.ColumnWidth = {190,93};
Post_processing_table.RowName = [];
Post_processing_input_data = num2cell(zeros(1,2));
Post_processing_input_data(:,2) = {false};
Post_processing_input_data(1,1) = {'Calculate Tortuosity'};

%Post_processing_input_data(1,1) = {'Checking Volume Fraction'};
%Post_processing_input_data(2,1) = {'Checking Particle Size'};
%Post_processing_input_data(3,1) = {'Calculate Connectivity'};
%Post_processing_input_data(4,1) = {'Calculate Tortuosity'};
%Post_processing_input_data(5,1) = {'Calculate Particle Size(C-PSD)'};

Post_processing_table_format =cell(1,2);
Post_processing_table_format(1,1) ={'numeric'};
Post_processing_table_format(1,2) ={'logical'};
Post_processing_table.ColumnFormat = Post_processing_table_format;
Post_processing_table.Data = Post_processing_input_data;
Post_processing_table.ColumnEditable = [false true];

%% Run Code tab 6
Title_tab6 = uicontrol('Parent',tab_run,'Style','text','String','Numerical generation','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]); % Heading

Display_run_text = uicontrol('Parent',tab_run,'Style','text','String','How many times do you want to generate the microstructure with the same parameters?','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'FontWeight','bold',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.01 0.9 0.9 0.04]);

Num_of_run_text = uicontrol('Parent',tab_run,'Style','text','String','Number of Runs','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.01,0.8,0.1,0.025],'Visible','on','Enable','on');
Num_of_run_box = uicontrol('Parent',tab_run,'Style','edit','String',[],'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.11,0.8,0.05,0.03],'Visible','on','Enable','on','CallBack',@callback_num_of_runs);

%Popupmenu_save_option = uicontrol('Parent',tab_run,'Style','text','String','Save microstructure','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.01,0.73,0.11,0.025],'Visible','off');
%popupmenu_tab6 = uicontrol('Parent',tab_run,'Style','popupmenu','String',{'Save all the microstructures';'Keep only the best microstructure based on input'},'Units','normalized','Position',[0.12,0.73,0.2,0.03],'Visible','off','FontSize',font_size_medium_GUI,'FontName',font_name_GUI);%'CallBack', @Callback_popupmenu);

Microstructure_generation_button = uicontrol( 'Parent',tab_run,'Style','pushbutton','String','Generate Microstructure','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'Units','normalized','Position',[0.4,0.55,0.2,0.05],'Visible','off','Enable','on','Callback',@Generate_Micro_Struct);
str = sprintf('Click to generate microstructure....');
Microstructure_generation_button.TooltipString = str;


Test_error_tab5 = uicontrol('Parent',tab_run, 'Style', 'text','Units','normalized','Position',[0.01 0.01 0.98 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Error message','Visible','off');

%% Generate Microstructure tab 7
Title_tab7 = uicontrol('Parent',tab_Generated_micro,'Style','text','String','Results','FontName',font_name_GUI,'FontSize',font_size_large_GUI,...
    'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,'Units','normalized','Position',[0,0.96,1,0.04]); % Heading

Display_results_text = uicontrol('Parent',tab_Generated_micro,'Style','text','String','Table will be updated for each finished run','FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'FontWeight','bold',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.01 0.9 0.9 0.04]);

% Table 1 
ltable_direction = [0.1 0.25 0.8 0.65];
table_output = uitable('Parent',tab_Generated_micro,'Units','normalized','Position',ltable_direction,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'CellEditCallBack',@edit_table_output,'Visible','off','Enable','off');
table_output.ColumnName = {'Attempt','Volume fraction','Tortuosity factor','Elapsed time'};
table_output.RowName = [];
table_output.ColumnWidth = {ltable*0.4, ltable*0.4, ltable*0.4, ltable*0.4};
%% CALLBACK FUNCTIONS for Phase and Volume Fraction Calculation tab1

% Callback function for voxel size and scaling factor
    function edit_voxelsize_scalingfactor_Callback(~,~)
        voxel_size_nm=str2double(Edit_voxelsize.String); % Get voxel size
        scaling_factor=str2double(Edit_Scalingfactor.String); % Get scaling factor
        tmp = table_direction.Data;
        number_voxels = cell2mat(tmp(:,2)); % Get dimension in voxels
        after_rescaling = number_voxels*scaling_factor;
        length_um = after_rescaling*voxel_size_nm/1000;
        column_1 = tmp(:,1);
        column_2 = tmp(:,2);
        column_3 = num2cell(after_rescaling);
        column_4 = num2cell(length_um);
        table_direction.Data=[column_1 column_2 column_3 column_4];
        if max(isnan(scaling_factor)) || max(isnan(voxel_size_nm))|| max(scaling_factor<0)|| max(voxel_size_nm<0)
            error_tab1_voxel_scaling = true;
            set(Test_error_tab1,'String','Error: Scaling factor and/or voxel size should be positive value','Visible','on');
        else
            error_tab1_voxel_scaling = false;
            check_all_error_tab1
        end
    end

% Callback function for direction table
    function cellsection_table_direction_Callback(~,~)
        edit_voxelsize_scalingfactor_Callback
        tmp = table_direction.Data;
        number_voxels = cell2mat(tmp(:,2));
        if max(isnan(number_voxels)) || max(number_voxels<0)
            error_tab1_num_voxel = true;
            set(Test_error_tab1,'String','Error: Number of voxels should be positive integer','Visible','on');
        else
            error_tab1_num_voxel = false;
            check_all_error_tab1
        end
    end

% Callback function for number of phase
    function edit_Num_phas (~,~)
        num_phas = str2double(Edit_Num_phas.String);
        min=1;
        if  rem(num_phas,1)==0 && num_phas >= min
            for k = 1:1:num_phas
                column_1(k,1)= {num2str(k)};
                column_2(k,1)= {num2str(k)};
                column_3(k,1)= {[ 'Phase ' num2str(k)]};
            end
            table_phas_num.Data= [column_1 column_2 column_3];
            set(table_phas_num,'Visible','on','Enable','on');
            set(Num_Slices,'Visible','on','Enable','on');
            set(Edit_Num_slices,'Visible','on','Enable','on');
            tmp1 = table_phas_num.Data;
            Phase_Name = tmp1(:,3)';
            Col_name = {'Position along the direction 3'};
            Col_edit = true;
            for k = 1:1:num_phas
                Col_name(1,k+1) = Phase_Name(k);
                Col_edit(1,k+1) = true;
            end
            Col_name(1,num_phas+2)={'Porosity'};
            Col_edit(1,k+2) = false;
            table_position.ColumnName = Col_name;
            table_position.ColumnEditable = Col_edit;
            set(table_position,'Visible','off','Enable','off')
            statuts_size_orientation_tab2 = zeros(num_phas,1);
            statuts_size_orientation_tab3 = zeros(num_phas,1);
            error_tab1_num_phase = false;
            check_all_error_tab1
            set(Text_saveconfirmation_tab1,'Visible','on')
        else
            set(Text_saveconfirmation_tab1,'Visible','off')
            error_tab1_num_phase = true;
            set(Test_error_tab1,'String','Error: Number of Phase cannot be less than 1 and has to be whole number','Visible','on');
        end
    end

% Callback function for phase number table
    function cellsection_table_phas_num_Callback (~,~)
        num_phas = str2double(Edit_Num_phas.String);
        tmp1 = table_phas_num.Data;
        Phase_Name = tmp1(:,3)';
        Code_check = str2double(tmp1(:,2));
        Col_name = {'Position along the direction'};
        for k = 1:1:num_phas
            Col_name(1,k+1) = Phase_Name(k);
            phase_popupmenu(k) = Phase_Name(k);
        end
        Col_name(1,num_phas+2)={'Porosity'};
        table_position.ColumnName = Col_name;
        set(popupmenu_phase_selection_tab2,'Style','popupmenu','String',phase_popupmenu,'Visible','off','Enable','off')
        set(popupmenu_phase_selection_tab3,'Style','popupmenu','String',phase_popupmenu,'Visible','off','Enable','off')
        set(table_position,'Visible','on','Enable','on')
        clear legend_x_volumefrac
        for k = 1:1:num_phas
            legend_x_volumefrac(k) = Phase_Name(k);
        end
        grid 'on'
        title('Volume Fraction Along Direction 3')
        % legend(legend_x_volumefrac,'Location','best')
        xlabel('Distance')
        ylabel('Volume Fraction')
        fig.color = 'w';
        set(volume_fraction_plot,'Visible','on')
        if max(isnan(Code_check))|| max(Code_check<0)
            error_tab1_code_check = true;
            set(Test_error_tab1,'String','Error: Code value should be positive integer','Visible','on');
        else
            error_tab1_code_check = false;
            check_all_error_tab1
        end
    end

% Callback function for number of slices
    function edit_Num_slice (~,~)
        num_slices = str2double(Edit_Num_slices.String);
        num_phas = str2double(Edit_Num_phas.String);
        min=2;
        set(Volumefraction_Plot_View_button,'Visible','on')
        if rem(num_slices,1)==0 && num_slices >= min
            table_position_input = zeros(num_slices, (num_phas + 2));
            table_position_input(:,1) = linspace(0,1,num_slices);
            table_position_input(:,2:(end-1)) = 1/(num_phas+1);
            table_position_input(:,end) = round(1 - sum(table_position_input(:,2:end-1),2),2);
            set(table_position,'Visible','on','Enable','on','Data',table_position_input);
            set(Volumefraction_Plot_View_button,'Enable','on');
            set(volume_fraction_plot, 'Visible','off')
            set(tab1_plot,'Visible','off');
            set(tab1_plot_legend,'Visible','off');
            error_tab1_num_slice = false;
            check_all_error_tab1
        else
            error_tab1_num_slice = true;
            set(Test_error_tab1,'String','Error: Number of Slices cannot be less than 2 and has to be whole number','Visible','on');
            set(Volumefraction_Plot_View_button,'Enable','off');
            set(volume_fraction_plot, 'Visible','off')
            set(tab1_plot,'Visible','off');
            set(tab1_plot_legend,'Visible','off');
        end
    end

% Callback function for position table
    function cellsection_table_position_Callback(~,~)
        
        tmp2 = table_position.Data;
        table_position_input1 = tmp2(:,2:end-1);
        tmp2(:,end) = round(1 - sum(table_position_input1,2),2);
        table_position.Data=tmp2;
        error_position = tmp2(:,1);
        %tmp2(:,end);
        if min(min(tmp2(:,2:end) >= 0))==0 || min(min(tmp2(:,2:end) <= 1))==0
            error_tab1_table3_porosity = true;
            set(Test_error_tab1,'String','Error: Porosity and/or Phase value cannot be less than 0 or more than 1' ,'Visible','on')
            set(Volumefraction_Plot_View_button, 'Enable','off')
        else
            error_tab1_table3_porosity = false;
            check_all_error_tab1
            set(Volumefraction_Plot_View_button, 'Enable','on')
            
            h =error_tab1_table3_column1;
            string = ('Error: Position along the direction 3 column values should be in ascending order');
            Enable = Volumefraction_Plot_View_button;
            current_tab = Test_error_tab1;
            all_error_check = @check_all_error_tab1;
            [h] = function_check_position(error_position,h,string,current_tab,all_error_check,Enable);
        end
    end

% Callback function for the volumn fraction plot
    function Callback_plot(~,~)
        set(Save_Input_button_tab1,'Visible','on')
        tmp3 = table_position.Data;
        num_phas = str2double(Edit_Num_phas.String);
        Marker_Counter = 1;
        Markers={'+','o','*','x','v','d','^','s','>','<'};
        x = tmp3(:,1);
        y = tmp3(:,2:end-1);
        y1 = tmp3(:,end);
        tab1_plot = plot(x,y,strcat('-',Markers{Marker_Counter}),x,y1,'--','LineWidth',2,'parent',volume_fraction_plot);
        Marker_Counter =Marker_Counter+1;
        ylim(volume_fraction_plot,[0 1])
        tmp1 = table_phas_num.Data;
        Phase_Name = tmp1(:,3)';
        set(volume_fraction_plot,'Visible','on')
        clear legend_x_volumefrac
        for k = 1:1:num_phas
            legend_x_volumefrac(k) = Phase_Name(k);
        end
        legend_x_volumefrac(k+1) = {'Porosity'};
        grid (volume_fraction_plot,'on');
        grid (volume_fraction_plot,'minor');
        title(volume_fraction_plot,'Volume Fraction Along Direction 3','FontSize',font_size_plot_title,'FontName',font_type_plot)
        tab1_plot_legend = legend(volume_fraction_plot,legend_x_volumefrac,'Location','best','FontName',font_type_plot);
        xlabel(volume_fraction_plot,'Distance','FontSize',font_size_plot_label_large,'FontName',font_type_plot)
        ylabel(volume_fraction_plot,'Volume Fraction','FontSize',font_size_plot_label_large,'FontName',font_type_plot)
        
    end

% Callback function to save all the input tab 1
    function Callback_save_input_data (~,~)
        % Input Parameters saved
        tmp = table_direction.Data;
        tmp1 = table_phas_num.Data;
        table_position.Data;
        tmp2 = table_position.Data;
        domain_size = cell2mat(tmp(:,2)');% get domain size
        voxel_size_nm=str2double(Edit_voxelsize.String); % Get voxel size
        scaling_factor=str2double(Edit_Scalingfactor.String); % Get Scaling factor
        num_phas = str2double(Edit_Num_phas.String);
        Phase_Name = tmp1(:,3);
        for k = 1:1:num_phas
            phase(k).name = cell2mat(tmp1(k,3));
            phase(k).code = str2double(cell2mat(tmp1(k,2)));
            phase(k).volumefraction.along_3rd_axis = [tmp2(:,1)'; tmp2(:,k+1)'];
            phase_popupmenu(k) = Phase_Name(k);
            col_name(1,k) = Phase_Name(k);
            Col_edit(1,k) = true;
        end
        phase.name; % Phase Name
        phase.code; % Phase Code
        phase(1).volumefraction.along_3rd_axis; % Volume Fraction
        
        % Tab 2 Parameters
        set(popupmenu_phase_selection_tab2,'Visible','on','Enable','on')
        set(popupmenu_phase_selection_tab2,'Style','popupmenu','String',phase_popupmenu)
        
        %Tab 3 Parameters
        set(popupmenu_phase_selection_tab3,'Visible','on','Enable','on')
        set(popupmenu_phase_selection_tab3,'Style','popupmenu','String',phase_popupmenu)
        % Tab 3 Parameters
        Table_maximum_overlapping.ColumnName = col_name; % Table 2 max overlapping
        Table_maximum_overlapping.RowName = col_name;
        Table_maximum_overlapping.ColumnEditable = Col_edit;
        Table_maximum_overlapping_input = zeros(num_phas,num_phas);
        Table_maximum_overlapping_input(:,:)= 0.5;
        Table_maximum_overlapping.Data = Table_maximum_overlapping_input;
        
        Table_minimum_particle_volume_conservated.RowName = col_name; % Table 1 min particle volume
        Table_minimum_particle_volume_conservated.ColumnName = {'Normalized Min Particle Volume'};
        Table_minimum_particle_volume_conservated_input = zeros(num_phas,1);
        Table_minimum_particle_volume_conservated_input(:,:)= 0.5;
        Table_minimum_particle_volume_conservated.Data = Table_minimum_particle_volume_conservated_input;
        
        fill_ratio_table.ColumnName = col_name;
        fill_ratio_table.ColumnEditable = Col_edit;
        
        set(Phase_overlapping_table,'Visible','on')
        set(Table_min_par_vol,'Visible','on')
        set(Number_of_generation_pass,'Visible','on')
        set(Edit_fillratio,'Visible','on','Enable','on')
        set(Table_maximum_overlapping,'Visible','on','Enable','on')
        set(Table_minimum_particle_volume_conservated,'Visible','on','Enable','on')
        set (fill_ratio_table, 'Visible','off','Enable','on')
        %set(overlapping_image,'Visible','on')
        %set(min_particle_volume_image,'Visible','on')
        
        statuts_size_orientation_tab2(current_phase_sizeorientation_tab2,1)=true;
        str = ['Configured_tab2 phase ' num2str(sum(statuts_size_orientation_tab2)) ' / ' num2str(length(statuts_size_orientation_tab2))];
        set(Configured_tab2,'String',str);
        
        statuts_size_orientation_tab3(current_phase_sizeorientation_tab3,1)=true;
        str = ['Configured_tab2 phase ' num2str(sum(statuts_size_orientation_tab3)) ' / ' num2str(length(statuts_size_orientation_tab3))];
        set(Configured_tab3,'String',str);
    end

%% CALLBACK FUNCTIONS for particle size and orientation tab2

%Callback function for the popupmenu
    function Callback_popupmenu(~,~)
        parameter_table_tab2.ColumnName = {'Parameters', 'Number of Slices', 'Number of Diameter/Elongation/Orientation'};
        Row(1).name='Direction Dx along 3rd Axis';
        Row(2).name='Direction Elongation Dx/Dy along 3rd Axis';
        Row(3).name='Direction Elongation Dx/Dz along 3rd Axis';
        
        NumberofSlices= [2;2;2];
        Numberofdiameter= [1;1;1];
        diameter= [10;1;1];
        parameter_table_tab2_input = [ {Row(1).name; Row(2).name; Row(3).name} num2cell(NumberofSlices) num2cell(Numberofdiameter)];
        parameter_table_tab2.Data = parameter_table_tab2_input;
        parameter_table_tab2.ColumnEditable = [false true true];
        parameter_table_tab2.RowName = [];
        ltable_1=scrsz(3)*lfigure(3)*(1-tabgroup_width)* parameter_table_position_tab2(3);
        parameter_table_tab2.ColumnWidth = {ltable_1*0.35, ltable_1*0.2, ltable_1*0.45};
        set(parameter_table_tab2,'Visible','on','Enable','on','Units','normalized');
        visual_check_table_input= [{Row(1).name; Row(2).name; Row(3).name} num2cell(diameter)];
        visual_check_table_tab2.Data = visual_check_table_input;
        visual_check_table_tab2.ColumnEditable = [false true];
        visual_check_table_tab2.ColumnWidth = {ltable_1*0.35, ltable_1*0.1};
        set(Diameter_Dz_along_3rd_axis,'Visible','off')
        set(Diameter_elongation_Dx_Dy_3rd_axis,'Visible','off')
        set(Diameter_elongation_Dx_Dz_3rd_axis,'Visible','off')
        set(Diameter_Dz_3rdaxis_table,'Visible','off','Enable','off')
        set(Diameter_Elongation_Dx_Dy_3rd_axis_table,'Visible','off','Enable','off')
        set(Diameter_Elongation_Dx_Dz_3rd_axis_table,'Visible','off','Enable','off')
        set(plot_area_view_button_diameter,'Visible','off','Enable','off')
        set(plot_area_view_button_diameter_plot_1,'Visible','off')
        set(plot_area_view_button_diameter_plot_2,'Visible','off')
        set(plot_area_view_button_diameter_plot_3,'Visible','off')
        set(A1,'Visible','off');
        set(A2,'Visible','off');
        set(A3,'Visible','off');
        set(a1,'Visible','off');
        set(a2,'Visible','off');
        set(a3,'Visible','off');
        set(Save_Input_button_tab2,'Visible','off');
        set(Test_error_tab2,'Visible','off')
        set(Cancel_button_tab2,'Visible','on')
        current_phase_sizeorientation_tab2 = get(popupmenu_phase_selection_tab2,'Value');
        
    end

%Callback function for the parameters table
    function parameters_callback(~,~)
        
        set(Diameter_Dz_along_3rd_axis,'Visible','on')
        set(Diameter_elongation_Dx_Dy_3rd_axis,'Visible','on')
        set(Diameter_elongation_Dx_Dz_3rd_axis,'Visible','on')
        set(visual_check_tab2,'Visible','on')
        set(visual_check_table_tab2,'Visible','on','Enable','on')
        set(visual_plot_particle_button_tab2,'Visible','on','Enable','on')
        % Direction along 3rd axis table
        min_1 =2;
        min_2 =1;
        tmp_1= parameter_table_tab2.Data;
        Numberofslices= cell2mat(tmp_1(:,2));
        Numberofdiameter= cell2mat(tmp_1(:,3));
        Num_of_slices_dia_3rdaxis = cell2mat(tmp_1(1,2));
        Dir_3rdaxis_diameter = cell2mat(tmp_1(1,3));
        % Column Name for the table 1-1 (number of slices and Diameter along
        % 3rd axis)
        if   min(Numberofslices)>= min_1 && max(rem(Numberofslices,1))== 0 && min(Numberofdiameter)>= min_2 && max(rem(Numberofdiameter,1))== 0
            Col_name = {'Position along 3rd axis'};
            Col_edit = true;
            for k = 1:1:Dir_3rdaxis_diameter
                Col_name(1,k+1) = {['Dx ' num2str(k)]};
                Col_edit(1,k+1) = true;
            end
            Col_name(1,k+2) = {'Total'};
            Col_edit(1,k+2) = false;
            Diameter_Dz_3rdaxis_table.ColumnName = Col_name;
            Diameter_Dz_3rdaxis_table.ColumnEditable = Col_edit;
            
            % Data for the table 1-1 (number of slices and Diameter along
            % 3rd axis)
            Diameter_3rdaxis_table_input = zeros(Num_of_slices_dia_3rdaxis+1,Dir_3rdaxis_diameter+2);
            Diameter_3rdaxis_table_input(2:end,1) = linspace(0,1,Num_of_slices_dia_3rdaxis);
            Diameter_3rdaxis_table_input(2:end,2:(end-1)) = 100/(Dir_3rdaxis_diameter);
            Diameter_3rdaxis_table_input(2:end,end) = round( sum(Diameter_3rdaxis_table_input(2:end,2:end-1),2),2);
            Diameter_3rdaxis_table_input(1,1)=NaN;
            Diameter_3rdaxis_table_input(1,end)=NaN;
            Diameter_3rdaxis_table_input(1,2:Dir_3rdaxis_diameter)=11;
            Diameter_3rdaxis_table_input(1,2:Dir_3rdaxis_diameter+1)=17;
            set(Diameter_Dz_3rdaxis_table,'Visible','on','Enable','on','Data',Diameter_3rdaxis_table_input);
            error_tab2_para_num_slices_num_dia = false;
            check_all_error_tab2
            set(plot_area_view_button_diameter,'Visible','on','Enable','on')
        else
            error_tab2_para_num_slices_num_dia = true;
            set(Diameter_Dz_3rdaxis_table,'Visible','off','Enable','off');
            set(Test_error_tab2,'String','Error: Number of Slices should be positive integer and greater than 2 and diameter should be greater than 1','Visible','on');
            set(plot_area_view_button_diameter,'Visible','on','Enable','off')
        end
        
        
        % Direction elongation 1st axis/ 2nd axis table
        Num_of_slices_dia_1st_2nd_axis = cell2mat(tmp_1(2,2));
        Dir_1st_2nd_diameter = cell2mat(tmp_1(2,3));
        % Column Name for the table 1-2 (number of slices and Diameter along
        % 1st and 2nd axis)
        if   min(Numberofslices)>= min_1 && max(rem(Numberofslices,1))== 0 && min(Numberofdiameter)>= min_2 && max(rem(Numberofdiameter,1))== 0
            Col_name = {'Position along 3rd axis'};
            Col_edit = true;
            for k = 1:1:Dir_1st_2nd_diameter
                Col_name(1,k+1) = {['Dx/Dy ' num2str(k)]};
                Col_edit(1,k+1) = true;
            end
            Col_name(1,k+2) = {'Total'};
            Col_edit(1,k+2) = false;
            Diameter_Elongation_Dx_Dy_3rd_axis_table.ColumnName = Col_name;
            Diameter_Elongation_Dx_Dy_3rd_axis_table.ColumnEditable = Col_edit;
            
            % Data for the table 1-2 (number of slices and Diameter along
            % 3rd axis)
            Diameter_Elongation_1st_2nd_axis_input = zeros(Num_of_slices_dia_1st_2nd_axis+1,Dir_1st_2nd_diameter+2);
            Diameter_Elongation_1st_2nd_axis_input(2:end,1) = linspace(0,1,Num_of_slices_dia_1st_2nd_axis);
            Diameter_Elongation_1st_2nd_axis_input(2:end,2:(end-1)) = 100/(Dir_1st_2nd_diameter);
            Diameter_Elongation_1st_2nd_axis_input(2:end,end) = round( sum(Diameter_Elongation_1st_2nd_axis_input(2:end,2:end-1),2),2);
            Diameter_Elongation_1st_2nd_axis_input(1,1)=NaN;
            Diameter_Elongation_1st_2nd_axis_input(1,end)=NaN;
            Diameter_Elongation_1st_2nd_axis_input(1,2:Dir_1st_2nd_diameter+1)=1;
            set(Diameter_Elongation_Dx_Dy_3rd_axis_table,'Visible','on','Enable','on','Data',Diameter_Elongation_1st_2nd_axis_input);
            error_tab2_para_num_slices_num_dia = false;
            check_all_error_tab2
        else
            error_tab2_para_num_slices_num_dia = true;
            set(Diameter_Elongation_Dx_Dy_3rd_axis_table,'Visible','off','Enable','off');
            set(Test_error_tab2,'String','Error: Number of Slices should be positive integer and greater than 2 and diameter should be greater than 1','Visible','on');
        end
        
        
        % Direction elongation 1st axis/ 3rd axis table
        Num_of_slices_dia_1st_3rd_axis = cell2mat(tmp_1(3,2));
        Dir_1st_3rd_diameter = cell2mat(tmp_1(3,3));
        % Column Name for the table 1-3 (number of slices and Diameter along
        % 1st and 2nd axis)
        if   min(Numberofslices)>= min_1 && max(rem(Numberofslices,1))== 0 && min(Numberofdiameter)>= min_2 && max(rem(Numberofdiameter,1))== 0
            Col_name = {'Position along 3rd axis'};
            Col_edit = true;
            for k = 1:1:Dir_1st_3rd_diameter
                Col_name(1,k+1) = {['Dx/Dz ' num2str(k)]};
                Col_edit(1,k+1) = true;
            end
            Col_name(1,k+2) = {'Total'};
            Col_edit(1,k+2) = false;
            Diameter_Elongation_Dx_Dz_3rd_axis_table.ColumnName = Col_name;
            Diameter_Elongation_Dx_Dz_3rd_axis_table.ColumnEditable = Col_edit;
            
            % Data for the table 1-3 (number of slices and Diameter along
            % 3rd axis)
            Diameter_Elongation_1st_3rd_axis_input = zeros(Num_of_slices_dia_1st_3rd_axis+1,Dir_1st_3rd_diameter+2);
            Diameter_Elongation_1st_3rd_axis_input(2:end,1) = linspace(0,1,Num_of_slices_dia_1st_3rd_axis);
            Diameter_Elongation_1st_3rd_axis_input(2:end,2:(end-1)) = 100/(Dir_1st_3rd_diameter);
            Diameter_Elongation_1st_3rd_axis_input(2:end,end) = round( sum(Diameter_Elongation_1st_3rd_axis_input(2:end,2:end-1),2),2);
            Diameter_Elongation_1st_3rd_axis_input(1,1)=NaN;
            Diameter_Elongation_1st_3rd_axis_input(1,end)=NaN;
            Diameter_Elongation_1st_3rd_axis_input(1,2:Dir_1st_3rd_diameter+1)=0.5;
            set(Diameter_Elongation_Dx_Dz_3rd_axis_table,'Visible','on','Enable','on','Data',Diameter_Elongation_1st_3rd_axis_input);
            error_tab2_para_num_slices_num_dia = false;
            check_all_error_tab2
        else
            error_tab2_para_num_slices_num_dia = true;
            set(Diameter_Elongation_Dx_Dz_3rd_axis_table,'Visible','off','Enable','off');
            set(Test_error_tab2,'String','Error: Number of Slices should be positive integer and greater than 2 and diameter should be greater than 1','Visible','on');
        end
    end

%Callback function fot the Direction along the 3rd axis table
    function Dia_Dz_3rdaxis_table(~,~)
        tmp1 = Diameter_Dz_3rdaxis_table.Data;
        tmp_1= parameter_table_tab2.Data;
        Num_of_slices_dia_3rdaxis= cell2mat(tmp_1(1,2));
        tmp1(1,1)=NaN;
        tmp1(1,end)=NaN;
        tmp2 = tmp1(2:end,2:end-1);
        tmp1(2:end,end)= round(sum(tmp2,2),3);
        Diameter_Dz_3rdaxis_table.Data = tmp1;
        if round(tmp1(2:end,end),3) ==100
            error_tab2_table1_a_sum = false;
            check_all_error_tab2
            set(plot_area_view_button_diameter,'Enable','on')
            for k = 1:1:Num_of_slices_dia_3rdaxis
                if max(tmp1(k,1)>tmp1(k+1,1))==1
                    error_tab2_table1_b_column1 = true;
                    set(Test_error_tab2,'String','Error: Position along the direction 3 column values should be in ascending order. Check Diameter along Dz axis table' ,'Visible','on')
                    set(plot_area_view_button_diameter, 'Enable','off')
                else
                    error_tab2_table1_b_column1 = false;
                    check_all_error_tab2
                    set(plot_area_view_button_diameter, 'Enable','on')
                    
                end
            end
        else
            error_tab2_table1_a_sum = true;
            set(Test_error_tab2,'String','Error: Sum of all the Diameter has to be equal to 100. Check Diameter along 3rd axis table','Visible','on');
            set(plot_area_view_button_diameter,'Enable','off')
        end
        
        
    end

%Callback function fot the Direction elongation 1st/2nd axis table
    function Dia_elongation_Dx_Dy_axis_table(~,~)
        tmp1 = Diameter_Elongation_Dx_Dy_3rd_axis_table.Data;
        tmp_1= parameter_table_tab2.Data;
        Dir_1st_2nd_diameter= cell2mat(tmp_1(2,2));
        tmp1(1,1)=NaN;
        tmp1(1,end)=NaN;
        tmp2 = tmp1(2:end,2:end-1);
        tmp1(2:end,end)= round(sum(tmp2,2),3);
        Diameter_Elongation_Dx_Dy_3rd_axis_table.Data = tmp1;
        if round(tmp1(2:end,end),3) ==100
            error_tab2_table2_a_sum = false;
            check_all_error_tab2
            set(plot_area_view_button_diameter,'Enable','on');
            for k = 1:1:Dir_1st_2nd_diameter
                if max(tmp1(k,1)>tmp1(k+1,1))==1
                    error_tab2_table2_b_column1 = true;
                    set(Test_error_tab2,'String','Error: Position along the direction 3 column values should be in ascending order. Check Diameter Elongation Dx/Dy axis table' ,'Visible','on')
                    set(plot_area_view_button_diameter, 'Enable','off')
                else
                    error_tab2_table2_b_column1 = false;
                    check_all_error_tab2
                    set(plot_area_view_button_diameter, 'Enable','on')
                end
            end
        else
            error_tab2_table2_a_sum = true;
            set(Test_error_tab2,'String','Error: Sum of all the diameter has to be equal to 100. Check Diameter Elongation 1st/2nd axis table','Visible','on');
            set(plot_area_view_button_diameter,'Enable','off');
        end
    end

%Callback function fot the Direction elongation 1st/3rd axis table
    function Dia_elongation_Dx_Dz_axis_table(~,~)
        tmp1 = Diameter_Elongation_Dx_Dz_3rd_axis_table.Data;
        tmp_1= parameter_table_tab2.Data;
        Dir_1st_3rd_diameter= cell2mat(tmp_1(3,2));
        tmp1(1,1)=NaN;
        tmp1(1,end)=NaN;
        tmp2 = tmp1(2:end,2:end-1);
        tmp1(2:end,end)= round(sum(tmp2,2),3);
        Diameter_Elongation_Dx_Dz_3rd_axis_table.Data = tmp1;
        if round(tmp1(2:end,end),3) ==100
            error_tab2_table3_a_sum = false;
            check_all_error_tab2
            set(plot_area_view_button_diameter,'Enable','on');
            for k = 1:1:Dir_1st_3rd_diameter
                if max(tmp1(k,1)>tmp1(k+1,1))==1
                    error_tab2_table3_b_column1 = true;
                    set(plot_area_view_button_diameter, 'Enable','off')
                else
                    error_tab2_table3_b_column1 = false;
                    check_all_error_tab2
                    set(plot_area_view_button_diameter, 'Enable','on')
                end
            end
        else
            error_tab2_table3_a_sum = true;
            set(Test_error_tab2,'String','Error: Sum of all the diameter has to be equal to 100. Check Diameter Elongation 1st/3rd axis table','Visible','on');
            set(plot_area_view_button_diameter,'Enable','off');
        end
    end


%Callback function for visual plot table tab2

    function visual_plot_table_tab2(~,~)
        tmp1 = visual_check_table_tab2.Data;
        tmp = cell2mat(tmp1(:,2));
    end

%Callback function for particle direction plot button tab2 
    function Visual_Plot_button_callback_tab2(~,~)
        tmp1 = visual_check_table_tab2.Data;
        tmp = cell2mat(tmp1(:,2));
        dx = tmp(1,:);
        dy = dx/tmp(2,:); 
        dz = dx/tmp(3,:); 
        mincheck = [dx,dy,dx];
        linelength = (1/3)*min(mincheck);       
        [binary_ellipsoid] = create_ellipsoid(dx,dy,dz);
        [f,v] = function_patch_3Darray(binary_ellipsoid);
        cla(axes1_tab2)
        %axes1_tab2 = axes('Parent',tab_particlesize,'LineWidth',linewidth,'Units','normalized','Position', [0.58 0.08 0.4 0.45]);
        hold(axes1_tab2,'on');
        plot3([0 linelength],[0 0],[0 0],'DisplayName','Axe 1','Parent',axes1_tab2);
        plot3([0 0],[0 linelength],[0 0],'DisplayName','Axe 2','Parent',axes1_tab2);
        plot3([0 0],[0 0],[0 linelength],'DisplayName','Axe 3','Parent',axes1_tab2);
        col=v(:,3);
        patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','interp','DisplayName','Particle','Parent',axes1_tab2);
        h=colorbar(axes1_tab2);
        ylabel(h, 'z position');
        xlabel('x - 1st axis')
        ylabel('y - 2nd axis')
        zlabel('z - 3rd axis (thickness)')
        axis(axes1_tab2,'equal');
        axis(axes1_tab2,'equal');
        legend_axes2_tab2 = {'x - 1st axis','y - 2nd axis','z - 3rd axis (thickness)','Particle'}';
        legend(axes1_tab2,legend_axes2_tab2,'Location','best');
        view(axes1_tab2,3)
        hold(axes1_tab2,'off');
    end 

%Callback function for area plot
    function Plot_Diameter_area_callback(~,~)
        
        % Area plot for Diameter along the 3rd axis
        tmp1 = Diameter_Dz_3rdaxis_table.Data;
        x1 = tmp1(2:end,1);
        y1 = tmp1(2:end,2:end-1);
        A1 = area(x1,y1,'parent',plot_area_view_button_diameter_plot_1);
        
        tmp_1= parameter_table_tab2.Data;
        Dir_3rdaxis_diameter = cell2mat(tmp_1(1,3));
        clear legend_D
        for k = 1:1: Dir_3rdaxis_diameter
            legend_D(k) = {['Dx' num2str(k)]};
        end
        a1 = legend(plot_area_view_button_diameter_plot_1,legend_D,'Location','best','FontName',font_type_plot);
        xlabel(plot_area_view_button_diameter_plot_1,'Distance along 3rd axis','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        ylabel(plot_area_view_button_diameter_plot_1,'%','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        set(plot_area_view_button_diameter_plot_1,'Visible','on');
        
        %Area plot for Diameter elongation 1st/2nd axis
        tmp2 = Diameter_Elongation_Dx_Dy_3rd_axis_table.Data;
        x2 = tmp2(2:end,1);
        y2 = tmp2(2:end,2:end-1);
        A2 = area(x2,y2,'parent',plot_area_view_button_diameter_plot_2);
        
        Dir_1st_2nd_diameter = cell2mat(tmp_1(2,3));
        clear legend_D
        for k = 1:1: Dir_1st_2nd_diameter
            legend_D(k) = {['Dx/Dy' num2str(k)]};
        end
        a2 = legend(plot_area_view_button_diameter_plot_2,legend_D,'Location','best','FontName',font_type_plot);
        xlabel(plot_area_view_button_diameter_plot_2,'Distance along 3rd axis','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        ylabel(plot_area_view_button_diameter_plot_2,'%','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        set(plot_area_view_button_diameter_plot_2,'Visible','on');
        
        %Area plot for Diameter elongation 1st/3rd axis
        tmp3 = Diameter_Elongation_Dx_Dz_3rd_axis_table.Data;
        x3 = tmp3(2:end,1);
        y3 = tmp3(2:end,2:end-1);
        A3 = area(x3,y3,'parent',plot_area_view_button_diameter_plot_3);
        
        Dir_1st_3rd_diameter = cell2mat(tmp_1(3,3));
        clear legend_D
        for k = 1:1: Dir_1st_3rd_diameter
            legend_D(k) = {['Dx/Dz' num2str(k)]};
        end
        a3 = legend(plot_area_view_button_diameter_plot_3,legend_D,'Location','best','FontName',font_type_plot);
        xlabel(plot_area_view_button_diameter_plot_3,'Distance along 3rd axis','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        ylabel(plot_area_view_button_diameter_plot_3,'%','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        set(plot_area_view_button_diameter_plot_3,'Visible','on');
        
        set(Save_Input_button_tab2,'Visible','on')
    end

%Callback function for cancel button
    function Callback_cancel_tab2(~,~)
        Callback_popupmenu
        statuts_size_orientation_tab2(current_phase_sizeorientation_tab2,1)=false;
        str = ['Configured phase ' num2str(sum(statuts_size_orientation_tab2)) ' / ' num2str(length(statuts_size_orientation_tab2))];
        set(Configured_tab2,'String',str);
    end

%Callback function to save all the input tab 2
    function Callback_save_input_data_tab2 (~,~)
        % Input Parameters saved
        tmp_1_1 = Diameter_Dz_3rdaxis_table.Data;
        tmp_1_2 = Diameter_Elongation_Dx_Dy_3rd_axis_table.Data;
        tmp_1_3 = Diameter_Elongation_Dx_Dz_3rd_axis_table.Data;
        
        phase(current_phase_sizeorientation_tab2).size_histogram.along_3rd_axis = (tmp_1_1(:,1:end-1));
        phase(current_phase_sizeorientation_tab2).elongation_histogram_dx_over_dy.along_3rd_axis = tmp_1_2(:,1:end-1);
        phase(current_phase_sizeorientation_tab2).elongation_histogram_dx_over_dz.along_3rd_axis = tmp_1_3(:,1:end-1);
        
        statuts_size_orientation_tab2(current_phase_sizeorientation_tab2,1)=true;
        str = ['Configured phase ' num2str(sum(statuts_size_orientation_tab2)) ' / ' num2str(length(statuts_size_orientation_tab2))];
        if sum(statuts_size_orientation_tab2) == length(statuts_size_orientation_tab2)
            set (Configured_tab2,'ForegroundColor',[1 1 1],'BackgroundColor',[0 0.5 0])
        else
            set (Configured_tab2,'ForegroundColor','r')
        end
        set(Configured_tab2,'String',str);
        
    end


%% CALLBACK FUNCTIONS for particle size and orientation tab3

%Callback function for the popupmenu
    function Callback_popupmenu_tab3(~,~)
        parameter_table_tab3.ColumnName = {'Parameters', 'Number of Slices', 'Number of Diameter/Elongation/Orientation'};
        Row(1).name='Orientation Normal to Ox Axis';
        Row(2).name='Orientation Normal to Oy Axis';
        Row(3).name='Orientation Normal to Oz Axis';
        NumberofSlices= [2;2;2];
        Numberofdiameter= [1;1;1];
        degree = [180;180;180];
        parameter_table_tab3_input = [ {Row(1).name; Row(2).name; Row(3).name} num2cell(NumberofSlices) num2cell(Numberofdiameter)];
        parameter_table_tab3.Data = parameter_table_tab3_input;
        parameter_table_tab3.ColumnEditable = [false true true];
        parameter_table_tab3.RowName = [];
        ltable_1=scrsz(3)*lfigure(3)*(1-tabgroup_width)* parameter_table_position_tab3(3);
        parameter_table_tab3.ColumnWidth = {ltable_1*0.35, ltable_1*0.2, ltable_1*0.45};
        set(parameter_table_tab3,'Visible','on','Enable','on','Units','normalized');
        visual_check_table_input_tab3= [{Row(1).name; Row(2).name; Row(3).name} num2cell(degree)];
        visual_check_table_tab3.Data = visual_check_table_input_tab3;
        visual_check_table_tab3.ColumnEditable = [false true];
        visual_check_table_tab3.ColumnWidth = {ltable_1*0.35, ltable_1*0.1};
        set(Orientation_normal_Dx_3rd_axis,'Visible','off')
        set(Orientation_normal_Dy_3rd_axis,'Visible','off')
        set(Orientation_normal_Dz_3rd_axis,'Visible','off')
        set(Orientation_normal_Dx_3rd_axis_table,'Visible','off','Enable','off')
        set(Orientation_normal_Dy_3rd_axis_table,'Visible','off','Enable','off')
        set(Orientation_normal_Dz_3rd_axis_table,'Visible','off','Enable','off')
        set(plot_area_view_button_orientation,'Visible','off','Enable','off')
        set(plot_area_view_button_orientation_plot_1,'Visible','off')
        set(plot_area_view_button_orientation_plot_2,'Visible','off')
        set(plot_area_view_button_orientation_plot_3,'Visible','off')
        set(A4,'Visible','off');
        set(A5,'Visible','off');
        set(A6,'Visible','off');
        set(a4,'Visible','off');
        set(a5,'Visible','off');
        set(a6,'Visible','off');
        set(Save_Input_button_tab3,'Visible','off');
        set(Test_error_tab3,'Visible','off')
        set(Cancel_button_tab3,'Visible','on')
        current_phase_sizeorientation_tab3 = get(popupmenu_phase_selection_tab3,'Value');
        
    end

%Callback function for the parameters table
    function parameters_callback_tab3(~,~)
        
        set(Orientation_normal_Dx_3rd_axis,'Visible','on')
        set(Orientation_normal_Dy_3rd_axis,'Visible','on')
        set(Orientation_normal_Dz_3rd_axis,'Visible','on')
        set(visual_check_tab3,'Visible','on')
        set(visual_check_table_tab3,'Visible','on','Enable','on')
        set(visual_plot_particle_button_tab3,'Visible','on','Enable','on')
        %Direction along 3rd axis table
        min_1 =2;
        min_2 =1;
        tmp_1= parameter_table_tab3.Data;
        Numberofslices= cell2mat(tmp_1(:,2));
        Numberofdiameter= cell2mat(tmp_1(:,3));
        
        %%Orientation (normal to 1st axis)
        Num_of_slices_normal_1st_axis = cell2mat(tmp_1(1,2));
        Numberorientation_normal_1st_axis = cell2mat(tmp_1(1,3));
        % Column Name for the table 2-1
        if   min(Numberofslices)>= min_1 && max(rem(Numberofslices,1))== 0 && min(Numberofdiameter)>= min_2 && max(rem(Numberofdiameter,1))== 0
            Col_name = {'Position along 3rd axis'};
            Col_edit = true;
            
            for k = 1:1:Numberorientation_normal_1st_axis
                Col_name(1,k+1) = {['O' num2str(k)]};
                Col_edit(1,k+1) = true;
            end
            Col_name(1,k+2) = {'Total'};
            Col_edit(1,k+2) = false;
            Col_name(1,k+3) = {'Randomize Orientation'};
            Col_edit(1,k+3) = true;
            Col_name(1,k+4) = {'Min'};
            Col_edit(1,k+4) = true;
            Col_name(1,k+5) = {'Max'};
            Col_edit(1,k+5) = true;
            
            Orientation_normal_Dx_3rd_axis_table.ColumnName = Col_name;
            Orientation_normal_Dx_3rd_axis_table.ColumnEditable = Col_edit;
            
            %Table 2-1 Input
            Degrees_1st_axis = linspace(0,180,Numberorientation_normal_1st_axis+1);
            Orientation_normal_1st_axis_table_input = zeros(Num_of_slices_normal_1st_axis+1,Numberorientation_normal_1st_axis+5);
            Orientation_normal_1st_axis_table_input(2:end,1) = linspace(0,1,Num_of_slices_normal_1st_axis);
            Orientation_normal_1st_axis_table_input(2:end,2:(end-4)) = 100/(Numberorientation_normal_1st_axis);
            Orientation_normal_1st_axis_table_input(2:end,end-3) = round( sum(Orientation_normal_1st_axis_table_input(2:end,2:end-4),2),2);
            Orientation_normal_1st_axis_table_input(1,2:end-4)= Degrees_1st_axis(2:end);
            Orientation_normal_1st_axis_table_input(1,1)=NaN;
            Orientation_normal_1st_axis_table_input(1,end-3)=NaN;
            Orientation_normal_1st_axis_table_input(2:end,end-1)=Degrees_1st_axis(2);
            Orientation_normal_1st_axis_table_input(2:end,end)=Degrees_1st_axis(end);
            Orientation_normal_1st_axis_table_input(1,end-1)=NaN;
            Orientation_normal_1st_axis_table_input(1,end)=NaN;
            Orientation_normal_1st_axis_table_input = num2cell(Orientation_normal_1st_axis_table_input);
            Orientation_normal_1st_axis_table_input(2:end,end-2)= {false};
            Orientation_normal_1st_axis_table_format =cell(1,Numberorientation_normal_1st_axis+5);
            Orientation_normal_1st_axis_table_format(:) ={'numeric'};
            Orientation_normal_1st_axis_table_format(end-2) ={'logical'};
            
            set(Orientation_normal_Dx_3rd_axis_table,'Visible','on','Enable','on','Data',Orientation_normal_1st_axis_table_input,'ColumnFormat',Orientation_normal_1st_axis_table_format)
            error_tab3_para_num_slices_num_dia = false;
            check_all_error_tab3
            set(plot_area_view_button_orientation,'Visible','on','Enable','on')
        else
            error_tab3_para_num_slices_num_dia = true;
            set(Orientation_normal_Dx_3rd_axis_table,'Visible','off','Enable','off')
            set(Test_error_tab3,'String','Error: Number of Slices should be positive integer and greater than 2 and diameter should be greater than 1','Visible','on');
            set(plot_area_view_button_orientation,'Visible','on','Enable','off')
        end
        
        
        %Orientation (normal to 2nd axis)
        Num_of_slices_normal_2nd_axis = cell2mat(tmp_1(2,2));
        Numberorientation_normal_2nd_axis = cell2mat(tmp_1(2,3));
        % Column Name for the table 2-2
        if   min(Numberofslices)>= min_1 && max(rem(Numberofslices,1))== 0 && min(Numberofdiameter)>= min_2 && max(rem(Numberofdiameter,1))== 0
            Col_name = {'Position along 3rd axis'};
            Col_edit = true;
            for k = 1:1:Numberorientation_normal_2nd_axis
                Col_name(1,k+1) = {['O' num2str(k)]};
                Col_edit(1,k+1) = true;
            end
            Col_name(1,k+2) = {'Total'};
            Col_edit(1,k+2) = false;
            Col_name(1,k+3) = {'Randomize Orientation'};
            Col_edit(1,k+3) = true;
            Col_name(1,k+4) = {'Min'};
            Col_edit(1,k+4) = true;
            Col_name(1,k+5) = {'Max'};
            Col_edit(1,k+5) = true;
            Orientation_normal_Dy_3rd_axis_table.ColumnName = Col_name;
            Orientation_normal_Dy_3rd_axis_table.ColumnEditable = Col_edit;
            
            %Table 2-2 Input
            Degrees_2nd_axis = linspace(0,180,Numberorientation_normal_2nd_axis+1);
            Orientation_normal_2nd_axis_table_input = zeros(Num_of_slices_normal_2nd_axis+1,Numberorientation_normal_2nd_axis+5);
            Orientation_normal_2nd_axis_table_input(2:end,1) = linspace(0,1,Num_of_slices_normal_2nd_axis);
            Orientation_normal_2nd_axis_table_input(2:end,2:(end-4)) = 100/(Numberorientation_normal_2nd_axis);
            Orientation_normal_2nd_axis_table_input(2:end,end-3) = round( sum(Orientation_normal_2nd_axis_table_input(2:end,2:end-4),2),2);
            Orientation_normal_2nd_axis_table_input(1,2:end-4)= Degrees_2nd_axis(2:end);
            Orientation_normal_2nd_axis_table_input(1,1)=NaN;
            Orientation_normal_2nd_axis_table_input(1,end-3)=NaN;
            Orientation_normal_2nd_axis_table_input(2:end,end-1)=Degrees_2nd_axis(2);
            Orientation_normal_2nd_axis_table_input(2:end,end)=Degrees_2nd_axis(end);
            Orientation_normal_2nd_axis_table_input(1,end-1)=NaN;
            Orientation_normal_2nd_axis_table_input(1,end)=NaN;
            Orientation_normal_2nd_axis_table_input = num2cell(Orientation_normal_2nd_axis_table_input);
            Orientation_normal_2nd_axis_table_input(2:end,end-2)= {false};
            Orientation_normal_2nd_axis_table_format =cell(1,Numberorientation_normal_2nd_axis+5);
            Orientation_normal_2nd_axis_table_format(:) ={'numeric'};
            Orientation_normal_2nd_axis_table_format(end-2) ={'logical'};
            
            set(Orientation_normal_Dy_3rd_axis_table,'Visible','on','Enable','on','Data',Orientation_normal_2nd_axis_table_input,'ColumnFormat',Orientation_normal_2nd_axis_table_format)
            error_tab3_para_num_slices_num_dia = false;
            check_all_error_tab3
            set(plot_area_view_button_orientation,'Visible','on','Enable','on')
        else
            error_tab3_para_num_slices_num_dia = true;
            set(Orientation_normal_Dy_3rd_axis_table,'Visible','off','Enable','off')
            set(Test_error_tab3,'String','Error: Number of Slices should be positive integer and greater than 2 and diameter should be greater than 1','Visible','on');
            set(plot_area_view_button_orientation,'Visible','on','Enable','off')
        end
        
        
        %Orientation (normal to 2nd axis)
        Num_of_slices_normal_3rd_axis = cell2mat(tmp_1(3,2));
        Numberorientation_normal_3rd_axis = cell2mat(tmp_1(3,3));
        % Column Name for the table 2-3
        if   min(Numberofslices)>= min_1 && max(rem(Numberofslices,1))== 0 && min(Numberofdiameter)>= min_2 && max(rem(Numberofdiameter,1))== 0
            Col_name = {'Position along 3rd axis'};
            Col_edit = true;
            for k = 1:1:Numberorientation_normal_3rd_axis
                Col_name(1,k+1) = {['O' num2str(k)]};
                Col_edit(1,k+1) = true;
            end
            Col_name(1,k+2) = {'Total'};
            Col_edit(1,k+2) = false;
            Col_name(1,k+3) = {'Randomize Orientation'};
            Col_edit(1,k+3) = true;
            Col_name(1,k+4) = {'Min'};
            Col_edit(1,k+4) = true;
            Col_name(1,k+5) = {'Max'};
            Col_edit(1,k+5) = true;
            Orientation_normal_Dz_3rd_axis_table.ColumnName = Col_name;
            Orientation_normal_Dz_3rd_axis_table.ColumnEditable = Col_edit;
            
            % Table 2-3 Input
            Degrees_3rd_axis = linspace(0,180,Numberorientation_normal_3rd_axis+1);
            Orientation_normal_3rd_axis_table_input = zeros(Num_of_slices_normal_3rd_axis+1,Numberorientation_normal_3rd_axis+5);
            Orientation_normal_3rd_axis_table_input(2:end,1) = linspace(0,1,Num_of_slices_normal_3rd_axis);
            Orientation_normal_3rd_axis_table_input(2:end,2:(end-4)) = 100/(Numberorientation_normal_3rd_axis);
            Orientation_normal_3rd_axis_table_input(2:end,end-3) = round( sum(Orientation_normal_3rd_axis_table_input(2:end,2:end-4),2),2);
            Orientation_normal_3rd_axis_table_input(1,2:end-4)= Degrees_3rd_axis(2:end);
            Orientation_normal_3rd_axis_table_input(1,1)=NaN;
            Orientation_normal_3rd_axis_table_input(1,end-3)=NaN;
            Orientation_normal_3rd_axis_table_input(2:end,end-1)=Degrees_3rd_axis(2);
            Orientation_normal_3rd_axis_table_input(2:end,end)=Degrees_3rd_axis(end);
            Orientation_normal_3rd_axis_table_input(1,end-1)=NaN;
            Orientation_normal_3rd_axis_table_input(1,end)=NaN;
            Orientation_normal_3rd_axis_table_input = num2cell(Orientation_normal_3rd_axis_table_input);
            Orientation_normal_3rd_axis_table_input(2:end,end-2)= {false};
            Orientation_normal_3rd_axis_table_format =cell(1,Numberorientation_normal_3rd_axis+5);
            Orientation_normal_3rd_axis_table_format(:) ={'numeric'};
            Orientation_normal_3rd_axis_table_format(end-2) ={'logical'};
            
            set(Orientation_normal_Dz_3rd_axis_table,'Visible','on','Enable','on','Data',Orientation_normal_3rd_axis_table_input,'ColumnFormat',Orientation_normal_3rd_axis_table_format)
            error_tab3_para_num_slices_num_dia = false;
            check_all_error_tab3
            set(plot_area_view_button_orientation,'Visible','on','Enable','on')
        else
            error_tab3_para_num_slices_num_dia = true;
            set(Orientation_normal_Dz_3rd_axis_table,'Visible','off','Enable','off')
            set(Test_error_tab3,'String','Error: Number of Slices should be positive integer and greater than 2 and diameter should be greater than 1','Visible','on');
            set(plot_area_view_button_orientation,'Visible','on','Enable','off')
        end
        
    end


%Callback function for orientation 1st axis table 1
    function Orientation_Dx_axis_callback (~,~)
        tmp= parameter_table_tab3.Data;
        tmp_1 = Orientation_normal_Dx_3rd_axis_table.Data;
        tmp1 = cell2mat(tmp_1(:,1:end-3));
        Numberorientation_normal_1st_axis = cell2mat(tmp(1,3));
        Num_of_slices_normal_1st_axis = cell2mat(tmp(1,2));
%         Degrees_1st_axis = linspace(0,180,Numberorientation_normal_1st_axis+1);
%         tmp1(1,2:end-1)= Degrees_1st_axis(2:end);
        tmp1(1,1) = NaN;
        tmp1(1,end) = NaN;
        tmp2 = tmp1(2:end,2:end-1);
        tmp1(2:end, end) = round(sum(tmp2,2),3);
        tmp_1(:,1:end-3)=num2cell(tmp1);
        Orientation_normal_Dx_3rd_axis_table.Data = tmp_1;
        error_position =tmp1(2:end,1);
        tmp3 = cell2mat(tmp_1(2:end,end-1:end));
        max_min_values = tmp1(1,2:end-1);
        check_box = cell2mat(tmp_1(2:end,end-2));
        
        if round(tmp1(2:end,end),3) ==100
            error_tab2_table4_a_sum = false;
            check_all_error_tab3
            set(Test_error_tab3,'String','Error: Sum of all the Diameter has to be equal to 100. Check Orientation (Normal to 1st axis) table','Visible','off');
            set(plot_area_view_button_orientation,'Enable','on')
            
            if max( tmp3(:,1) > tmp3(:,2) )
                error_tab2_table4_c_minmax = true;
                set(Test_error_tab3,'String','Error: The Min value cannot be greater than Max value. Check Orientation Normal to Ox axis table','Visible','on');
                set(plot_area_view_button_orientation,'Enable','off')
            else
                error_tab2_table4_c_minmax = false;
                check_all_error_tab3
                set(Test_error_tab3,'String','Error: The Min value cannot be greater than Max value','Visible','off')
                set(plot_area_view_button_orientation,'Enable','on');
                if min(min(ismember(tmp3, max_min_values)))==0
                    error_tab2_table4_d_minmax_value_match = true;
                    set(Test_error_tab3,'String','Error: The Min and Max values do not correspond to the orientation degree values. Check Orientation Normal to Ox axis table','Visible','on')
                    set(plot_area_view_button_orientation,'Enable','off');
                else
                    error_tab2_table4_d_minmax_value_match = false;
                    check_all_error_tab3
                    set(plot_area_view_button_orientation,'Enable','on');
                    
                    h =error_tab2_table4_b_column1;
                    string = 'Error: Position along the direction 3 column values should be in ascending order. Check Orientation Normal to Ox axis table';
                    Enable = plot_area_view_button_orientation;
                    current_tab = Test_error_tab3;
                    all_error_check = @check_all_error_tab3;
                    [h] = function_check_position(error_position,h,string,current_tab,all_error_check,Enable);
                    
                    for k = 1:1:Num_of_slices_normal_1st_axis
                        if check_box(k,1)== 1
                            a = min(max_min_values >= tmp3(k,1),max_min_values <= tmp3(k,2));
                            b = 100/sum(a);
                            c = find(a==1);
                            d = find(a==0);
                            tmp2(k,c)=b;
                            tmp2(k,d)=0;
                            tmp1(2:end,2:end-1) = tmp2;
                            tmp1(2:end, end) = round(sum(tmp2,2),3);
                            tmp_1(:,1:end-3)=num2cell(tmp1);
                            Orientation_normal_Dx_3rd_axis_table.Data = tmp_1;
                        end
                    end
                end
            end
        else
            error_tab2_table4_a_sum = true;
            set(Test_error_tab3,'String','Error: Sum of all the Diameter has to be equal to 100. Check Orientation Normal to Ox axis table','Visible','on');
            set(plot_area_view_button_orientation,'Enable','off')
        end
    end

%Callback function for orientation 2nd axis table 2
    function Orientation_Dy_axis_callback (~,~)
        tmp= parameter_table_tab3.Data;
        tmp_1 = Orientation_normal_Dy_3rd_axis_table.Data;
        tmp1 = cell2mat(tmp_1(:,1:end-3));
        Num_of_slices_normal_1st_axis = cell2mat(tmp(2,2));
        Numberorientation_normal_1st_axis = cell2mat(tmp(2,3));
%         Degrees_1st_axis = linspace(0,180,Numberorientation_normal_1st_axis+1);
%         tmp1(1,2:end-1)= Degrees_1st_axis(2:end);
        tmp1(1,1) = NaN;
        tmp1(1,end) = NaN;
        tmp2 = tmp1(2:end,2:end-1);
        tmp1(2:end, end) = round(sum(tmp2,2),3);
        tmp_1(:,1:end-3)=num2cell(tmp1);
        Orientation_normal_Dy_3rd_axis_table.Data = tmp_1;
        error_position =tmp1(2:end,1);
        tmp3 = cell2mat(tmp_1(2:end,end-1:end));
        max_min_values = tmp1(1,2:end-1);
        check_box = cell2mat(tmp_1(2:end,end-2)) ;
        
        if round(tmp1(2:end,end),3) ==100
            error_tab2_table5_a_sum = false;
            check_all_error_tab3
            set(plot_area_view_button_orientation,'Enable','on')
            
            if max( tmp3(:,1) > tmp3(:,2) )
                error_tab2_table5_c_minmax = true;
                set(Test_error_tab3,'String','Error: The Min value cannot be greater than Max value. Check Orientation Normal to Oy axis table','Visible','on');
                set(plot_area_view_button_orientation,'Enable','off')
            else
                error_tab2_table5_c_minmax = false;
                check_all_error_tab3
                set(plot_area_view_button_orientation,'Enable','on');
                if min(min(ismember(tmp3, max_min_values)))==0
                    error_tab2_table5_d_minmax_value_match = true;
                    set(Test_error_tab3,'String','Error: The Min and Max values do not correspond to the orientation degree values. Check Orientation Normal to Oy axis table','Visible','on')
                    set(plot_area_view_button_orientation,'Enable','off');
                else
                    error_tab2_table5_d_minmax_value_match = false;
                    check_all_error_tab3
                    set(plot_area_view_button_orientation,'Enable','on');
                    
                    h =error_tab2_table5_b_column1;
                    string = 'Error: Position along the direction 3 column values should be in ascending order. Check Orientation Normal to Oy axis table';
                    Enable = plot_area_view_button_orientation;
                    current_tab = Test_error_tab3;
                    all_error_check = @check_all_error_tab3;
                    [h] = function_check_position(error_position,h,string,current_tab,all_error_check,Enable);
                    
                    for k = 1:1:Num_of_slices_normal_1st_axis
                        if check_box(k,1)== 1
                            a = min(max_min_values >= tmp3(k,1),max_min_values <= tmp3(k,2));
                            b = 100/sum(a);
                            c = find(a==1);
                            d = find(a==0);
                            tmp2(k,c)=b;
                            tmp2(k,d)=0;
                            tmp1(2:end,2:end-1) = tmp2;
                            tmp1(2:end, end) = round(sum(tmp2,2),3);
                            tmp_1(:,1:end-3)=num2cell(tmp1);
                            Orientation_normal_Dy_3rd_axis_table.Data = tmp_1;
                        end
                    end
                    
                end
            end
            
        else
            error_tab2_table5_a_sum = true;
            set(Test_error_tab3,'String','Error: Sum of all the Diameter has to be equal to 100. Check Orientation Normal to Oy axis table','Visible','on');
            set(plot_area_view_button_orientation,'Enable','off')
        end
        
    end

%Callback function for orientation 3rd axis table 3
    function Orientation_Dz_axis_callback (~,~)
        tmp= parameter_table_tab3.Data;
        tmp_1 = Orientation_normal_Dz_3rd_axis_table.Data;
        tmp1 = cell2mat(tmp_1(:,1:end-3));
        Num_of_slices_normal_1st_axis = cell2mat(tmp(3,2));
        Numberorientation_normal_1st_axis = cell2mat(tmp(3,3));
%         Degrees_1st_axis = linspace(0,180,Numberorientation_normal_1st_axis+1);
%         tmp1(1,2:end-1)= Degrees_1st_axis(2:end);
        tmp1(1,1) = NaN;
        tmp1(1,end) = NaN;
        tmp2 = tmp1(2:end,2:end-1);
        tmp1(2:end, end) = round(sum(tmp2,2),3);
        tmp_1(:,1:end-3)=num2cell(tmp1);
        Orientation_normal_Dz_3rd_axis_table.Data = tmp_1;
        error_position =tmp1(2:end,1);
        tmp3 = cell2mat(tmp_1(2:end,end-1:end));
        max_min_values = tmp1(1,2:end-1);
        check_box = cell2mat(tmp_1(2:end,end-2));
        
        if round(tmp1(2:end,end),3) ==100
            error_tab2_table6_a_sum = false;
            check_all_error_tab3
            set(plot_area_view_button_orientation,'Enable','on')
            
            if max( tmp3(:,1) > tmp3(:,2) )
                error_tab2_table6_c_minmax = true;
                set(Test_error_tab3,'String','Error: The Min value cannot be greater than Max value. Check Orientation Normal to Oz axis table','Visible','on');
                set(plot_area_view_button_orientation,'Enable','off')
            else
                error_tab2_table6_c_minmax = false;
                check_all_error_tab3
                set(plot_area_view_button_orientation,'Enable','on');
                
                if min(min(ismember(tmp3, max_min_values)))==0
                    error_tab2_table6_d_minmax_value_match = true;
                    set(Test_error_tab3,'String','Error: The Min and Max values do not correspond to the orientation degree values. Check Orientation Normal to Oz axis table','Visible','on')
                    set(plot_area_view_button_orientation,'Enable','off');
                else
                    error_tab2_table6_d_minmax_value_match = false;
                    check_all_error_tab3
                    set(plot_area_view_button_orientation,'Enable','on');
                    
                    h =error_tab2_table6_b_column1;
                    string = 'Error: Position along the direction 3 column values should be in ascending order. Check Orientation Normal to Oz axis table';
                    Enable = plot_area_view_button_orientation;
                    current_tab = Test_error_tab3;
                    all_error_check = @check_all_error_tab3;
                    [h] = function_check_position(error_position,h,string,current_tab,all_error_check,Enable);
                    
                    
                    for k = 1:1:Num_of_slices_normal_1st_axis
                        if check_box(k,1)== 1
                            a = min(max_min_values >= tmp3(k,1),max_min_values <= tmp3(k,2));
                            b = 100/sum(a);
                            c = find(a==1);
                            d = find(a==0);
                            tmp2(k,c)=b;
                            tmp2(k,d)=0;
                            tmp1(2:end,2:end-1) = tmp2;
                            tmp1(2:end, end) = round(sum(tmp2,2),3);
                            tmp_1(:,1:end-3)=num2cell(tmp1);
                            Orientation_normal_Dz_3rd_axis_table.Data = tmp_1;
                        end
                    end
                end
                
            end
            
        else
            error_tab2_table6_a_sum = true;
            set(Test_error_tab3,'String','Error: Sum of all the Diameter has to be equal to 100. Check Orientation Normal to Oz axis table','Visible','on');
            set(plot_area_view_button_orientation,'Enable','off')
        end
    end

%Callback function for visual plot table tab3 
    function visual_plot_table_tab3(~,~)
        tmp2 = visual_check_table_tab3.Data;
        tmp3 = cell2mat(tmp2(:,2));
    end 

%Callback function for particle orentation plot button tab3
    function Visual_Plot_button_callback_tab3(~,~)
        tmp1 = visual_check_table_tab2.Data;
        tmp = cell2mat(tmp1(:,2));
        tmp2 = visual_check_table_tab3.Data;
        tmp3 = cell2mat(tmp2(:,2));        
        dx = tmp(1,:);
        dy = dx/tmp(2,:); 
        dz = dx/tmp(3,:);         
        mincheck = [dx,dy,dx];
        linelength = (1/3)*min(mincheck);         
        angle_x_deg = tmp3(1,:);
        angle_y_deg = tmp3(2,:);
        angle_z_deg = tmp3(3,:);
        [binary_ellipsoid] = create_ellipsoid(dx,dy,dz);
        [binary_ellipsoid] = rotate_domain(binary_ellipsoid,angle_x_deg, angle_y_deg, angle_z_deg);
        [f,v] = function_patch_3Darray(binary_ellipsoid);
        %set(axes2_tab3,'Visible','on');
        cla(axes2_tab3)
        %axes2_tab3 = axes('Parent',tab_orientation,'Units','normalized','Position', [0.58 0.08 0.4 0.45]);
        hold(axes2_tab3,'on');
        plot3([0 linelength],[0 0],[0 0],'DisplayName','Axe 1','Parent',axes2_tab3);
        plot3([0 0],[0 linelength],[0 0],'DisplayName','Axe 2','Parent',axes2_tab3);
        plot3([0 0],[0 0],[0 linelength],'DisplayName','Axe 3','Parent',axes2_tab3);
        col=v(:,3);
        patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','interp','DisplayName','Particle','Parent',axes2_tab3);
        h=colorbar(axes2_tab3);
        ylabel(h, 'z position');
        xlabel('x - 1st axis')
        ylabel('y - 2nd axis')
        zlabel('z - 3rd axis (thickness)')
        axis(axes2_tab3,'equal');
        legend_axes2_tab3 = {'x - 1st axis','y - 2nd axis','z - 3rd axis (thickness)','Particle'}';
        legend(axes2_tab3,legend_axes2_tab3,'Location','best');
        view(axes2_tab3,3)
        hold(axes2_tab3,'off');
    end


%Callback function for area plot 
    function Plot_orientation_area_callback(~,~)
        
        % Area plot for Diameter along the 3rd axis
        tmp1 = Orientation_normal_Dx_3rd_axis_table.Data;
        x1 = cell2mat(tmp1(2:end,1));
        y1 = cell2mat(tmp1(2:end,2:end-4));
        A4 = area(x1,y1,'parent',plot_area_view_button_orientation_plot_1);
        
        tmp_1= parameter_table_tab3.Data;
        Dir_3rdaxis_diameter = cell2mat(tmp_1(1,3));
        clear legend_D
        for k = 1:1: Dir_3rdaxis_diameter
            legend_D(k) = {['Ox' num2str(k)]};
        end
        a4 = legend(plot_area_view_button_orientation_plot_1,legend_D,'Location','best','FontName',font_type_plot);
        xlabel(plot_area_view_button_orientation_plot_1,'Distance along 3rd axis','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        ylabel(plot_area_view_button_orientation_plot_1,'%','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        set(plot_area_view_button_orientation_plot_1,'Visible','on');
        
        %Area plot for Diameter elongation 1st/2nd axis
        tmp2 = Orientation_normal_Dy_3rd_axis_table.Data;
        x2 = cell2mat(tmp2(2:end,1));
        y2 = cell2mat(tmp2(2:end,2:end-4));
        A5 = area(x2,y2,'parent',plot_area_view_button_orientation_plot_2);
        
        Dir_1st_2nd_diameter = cell2mat(tmp_1(2,3));
        clear legend_D
        for k = 1:1: Dir_1st_2nd_diameter
            legend_D(k) = {['Oy' num2str(k)]};
        end
        a5 = legend(plot_area_view_button_orientation_plot_2,legend_D,'Location','best','FontName',font_type_plot);
        xlabel(plot_area_view_button_orientation_plot_2,'Distance along 3rd axis','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        ylabel(plot_area_view_button_orientation_plot_2,'%','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        set(plot_area_view_button_orientation_plot_2,'Visible','on');
        
        %Area plot for Diameter elongation 1st/3rd axis
        tmp3 = Orientation_normal_Dz_3rd_axis_table.Data;
        x3 = cell2mat(tmp3(2:end,1));
        y3 = cell2mat(tmp3(2:end,2:end-4));
        A6 = area(x3,y3,'parent',plot_area_view_button_orientation_plot_3);
        
        Dir_1st_3rd_diameter = cell2mat(tmp_1(3,3));
        clear legend_D
        for k = 1:1: Dir_1st_3rd_diameter
            legend_D(k) = {['Oz' num2str(k)]};
        end
        a6 = legend(plot_area_view_button_orientation_plot_3,legend_D,'Location','best','FontName',font_type_plot);
        xlabel(plot_area_view_button_orientation_plot_3,'Distance along 3rd axis','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        ylabel(plot_area_view_button_orientation_plot_3,'%','FontSize',font_size_plot_label_small,'FontName',font_type_plot)
        set(plot_area_view_button_orientation_plot_3,'Visible','on');
        
        set(Save_Input_button_tab3,'Visible','on')
    end

%Callback function for cancel button
    function Callback_cancel_tab3(~,~)
        Callback_popupmenu_tab3
        statuts_size_orientation_tab3(current_phase_sizeorientation_tab3,1)=false;
        str = ['Configured phase ' num2str(sum(statuts_size_orientation_tab3)) ' / ' num2str(length(statuts_size_orientation_tab3))];
        set(Configured_tab3,'String',str);
    end

%Callback function to save all the input tab 1
    function Callback_save_input_data_tab3 (~,~)
        % Input Parameters saved
        tmp_2_1 = Orientation_normal_Dx_3rd_axis_table.Data;
        tmp_2_2 = Orientation_normal_Dy_3rd_axis_table.Data;
        tmp_2_3 = Orientation_normal_Dz_3rd_axis_table.Data;
        
        phase(current_phase_sizeorientation_tab3).orientation_histogram_angledeg_x.along_3rd_axis = cell2mat(tmp_2_1(:,1:end-4));
        phase(current_phase_sizeorientation_tab3).orientation_histogram_angledeg_y.along_3rd_axis = cell2mat(tmp_2_2(:,1:end-4));
        phase(current_phase_sizeorientation_tab3).orientation_histogram_angledeg_z.along_3rd_axis = cell2mat(tmp_2_3 (:,1:end-4));
        
        statuts_size_orientation_tab3(current_phase_sizeorientation_tab3,1)=true;
        str = ['Configured phase ' num2str(sum(statuts_size_orientation_tab3)) ' / ' num2str(length(statuts_size_orientation_tab3))];
        if sum(statuts_size_orientation_tab3) == length(statuts_size_orientation_tab3)
            set (Configured_tab3,'ForegroundColor',[1 1 1],'BackgroundColor',[0 0.5 0])
        else
            set (Configured_tab3,'ForegroundColor','r')
        end
        set(Configured_tab3,'String',str);
        
    end
%% CALLBACK FUNCTION for Particle Overlapping and fill ratio tab 4

%Callback function for Table maximum overlapping table
    function Table_maximum_overlapping_callback(~,~)
        
        tmp=Table_maximum_overlapping.Data;
        num_phas = str2double(Edit_Num_phas.String);
        for i = 1:1:num_phas
            for j = 1:1:num_phas
                tmp(i,j) = tmp(j,i);
                tmp(j,i) = tmp(i,j);
            end
        end
        Table_maximum_overlapping.Data=tmp;
    end

%Callback function for Number of generation
    function  edit_fillratio (~,~)
        pass = str2double(Edit_fillratio.String);
        num_phas = str2double(Edit_Num_phas.String);
        for k = 1:1:pass
            rowname(k,1) = {['Pass ' num2str(k)]} ;
        end
        fill_ratio_table.RowName = rowname;
        fill_ratio_table_data = zeros(pass,num_phas);
        lin = linspace(0.5,1,pass);
        for k = 1:1:num_phas
            fill_ratio_table_data(:,k) = lin;
        end
        fill_ratio_table.Data = fill_ratio_table_data;
        
        set(fill_ratio,'Visible','on')
        set(fill_ratio_table, 'Visible','on','Enable','on')
        set(Save_Input_button_tab4,'Visible','on')
    end

%Callback function for fill ratio table
    function fillratio (~,~)
        
    end

%Callback function for save button tab3
    function Callback_save_tab2(~,~)
        
        num_phas = str2double(Edit_Num_phas.String);
        tmp = Table_maximum_overlapping.Data;
        tmp2 = fill_ratio_table.Data;
        
        Maximum_interpenetration=tmp;
        Minimum_particle_volume_conservated = Table_minimum_particle_volume_conservated.Data';
        
        for k = 1:1:num_phas
            phase(k).fillratio = tmp2(:,k);
        end
        phase.fillratio;
    end

%% CALLBACK FUNCTION for Post-processing and save options tab 4

% Callback function for save option
    function button_select_savefolder_Callback(~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select the save folder where figure and results will be saved';
        % Open dialog box to choose file path
        Save_folder_tmp = uigetdir(str_dialogbox);
        if Save_folder_tmp==0
            Save_folder = 'None';
        else
            Save_folder = [Save_folder_tmp folder_separation];
        end
        set(Text_savefolder,'String',Save_folder);
    end

% Callback function for Fontstyle Fontsize and Grid
    function Font_Grid_Callback (~,~)
        Font_Grid_options.Data;
        tmp1 = cell2mat(Font_Grid_options.Data(1,1:end-2));
        tmp2 = cell2mat(Font_Grid_options.Data(1,end-1:end));
        OPTIONS.axe_fontsize = tmp1(1,1);
        OPTIONS.legend_fontsize = tmp1(1,2);
        OPTIONS.title_fontsize = tmp1(1,3);
        OPTIONS.figure_Linewidth = tmp1(1,4);
        OPTIONS.grid = tmp2(1,1);
        OPTIONS.minorgrid = tmp2(1,2);
        
    end

% Callback function for font name
    function Font_name_CallBack (~,~)
        OPTIONS.fontname = Popup_displayoptions_fontname.FontName;
        
    end


%% CALLBACK FUNCTION for Run Code tab 6

% callback num of runs
    function callback_num_of_runs (~,~)
        num_of_run = str2double(Num_of_run_box.String);
        min = 1;
        if  rem(num_of_run,1)==0 && num_of_run >= min
            error_tab5_num_of_run = false;
            check_all_error_tab5
            set(Microstructure_generation_button,'Visible','on')
            %set(Popupmenu_save_option,'Visible','on')
            %set(popupmenu_tab6,'Visible','on')
            set(table_output,'Visible','on','Enable','on')
                          

        else
            error_tab5_num_of_run = true;
            set(Test_error_tab5,'String','Error: Number of Runs cannot be less than 1 and has to be whole number','Visible','on');
            set(Microstructure_generation_button,'Visible','off')
            %set(Popupmenu_save_option,'Visible','off')
            %set(popupmenu_tab6,'Visible','off')
        end
        
    end

   
    function Generate_Micro_Struct(~,~)
        num_of_run = str2double(Num_of_run_box.String);
        set(table_output,'Visible','on','Enable','on')
        for k = 1:1:num_of_run
            % Generatinig Microstructure
            tic;

            [microstructure3D, phase] = function_generate_ellipsoid_microstructure(domain_size,phase,Maximum_interpenetration,Minimum_particle_volume_conservated);
            t = toc; % Wall clock time
            
            % Volume Fraction for Generated Microstructure table
            
            domain_size = size(microstructure3D.phase);
            number_phase = (length(phase)) ;
            voxel_number = prod(domain_size);
            
            choices = cell2mat(Post_processing_table.Data(:,2));
            [Volume_fraction,Tortuosity_factor] = simpleevaluation_microstructure(number_phase,voxel_number,microstructure3D,choices);
            column_1(k,1)= {num2str(k)};
            column_2(k,1)= {num2str(Volume_fraction)};
            column_3(k,1)= {num2str(Tortuosity_factor)};
            column_4(k,1)= {num2str(t)};
            table_output.Data = [column_1 column_2 column_3 column_4];
            
            % function_check_generatedmicrostructure(microstructure3D,phase,voxel_size_nm,OPTIONS);
            
            str_filename = [Text_savefolder.String 'Volume_' num2str(k)];
            function_save_tif( uint8(microstructure3D.phase),[str_filename '_phase.tif']);
            function_save_tif( uint16(microstructure3D.particle_id),[str_filename '_particle.tif']);
            
            % Rescale
            parameters_scaling.scaling_factor = 1/(str2double(Edit_Scalingfactor.String)); % Get scaling factor
            if parameters_scaling.scaling_factor~=1
                parameters_scaling.label_or_greylevel = 'Label';
                parameters_scaling.background = 0;
                [Microstructure_resized] = uint8(function_scaling(uint8(microstructure3D.phase),parameters_scaling));
                str_filename = [Text_savefolder.String 'Volume_' num2str(k) '_phase_rescaled.tif'];
                function_save_tif(Microstructure_resized,str_filename);
            end
        end
        
        
        
    end

%% Error check functions

% Error to check ascending order in column
    function [h] = function_check_position(error_position,h,string,current_tab,all_error_check,Enable)
        if sum(isnan(error_position))==0
            if error_position(1)== 0 && error_position(end) ==1
                if sum(error_position == sort(error_position))==length(error_position)
                    h = false;
                    all_error_check;
                    set(Enable, 'Enable','on')
                else
                    h = true;
                    set(current_tab,'String',string,'Visible','on')
                    set(Enable, 'Enable','off')
                end
            else
                h = true;
                set(current_tab,'String',string,'Visible','on')
                set(Enable, 'Enable','off')
            end
        else
            h = true;
            set(current_tab,'String',string,'Visible','on')
            set(Enable, 'Enable','off')
        end
        
    end

% Error function for tab 1
    function check_all_error_tab1
        if (error_tab1_voxel_scaling + error_tab1_num_voxel + error_tab1_num_phase + error_tab1_num_slice + error_tab1_table3_column1 + error_tab1_table3_porosity) == 0
            set (Test_error_tab1,'Visible','off');
        end
    end

% Error function for tab 2
    function check_all_error_tab2
        if (error_tab2_para_num_slices_num_dia + error_tab2_table1_a_sum + error_tab2_table1_b_column1 + error_tab2_table2_a_sum + error_tab2_table2_b_column1 + error_tab2_table3_a_sum + error_tab2_table3_b_column1 + error_tab2_table4_a_sum + error_tab2_table4_b_column1 + error_tab2_table4_c_minmax + error_tab2_table4_d_minmax_value_match + error_tab2_table5_a_sum + error_tab2_table5_b_column1 + error_tab2_table5_c_minmax + error_tab2_table5_d_minmax_value_match + error_tab2_table6_a_sum + error_tab2_table6_b_column1 + error_tab2_table6_c_minmax + error_tab2_table6_d_minmax_value_match) ==0
            set (Test_error_tab2,'Visible','off')
        end
    end

% Error function for tab 3
    function check_all_error_tab3
        if (error_tab3_para_num_slices_num_dia + error_tab2_table4_a_sum + error_tab2_table4_b_column1 + error_tab2_table4_c_minmax + error_tab2_table4_d_minmax_value_match + error_tab2_table5_a_sum +  error_tab2_table5_b_column1 +  error_tab2_table5_c_minmax +  error_tab2_table5_d_minmax_value_match + error_tab2_table6_a_sum +  error_tab2_table6_b_column1 +  error_tab2_table6_c_minmax +  error_tab2_table6_d_minmax_value_match) == 0
            set (Test_error_tab3,'Visible','off')
        end
    end

% Error function for tab 5
    function check_all_error_tab5
        if (error_tab5_num_of_run) ==0
            set (Test_error_tab5,'Visible','off')
        end
    end

%% Elliposid function 

% Direction Function 
function [binary_ellipsoid] = create_ellipsoid(dx,dy,dz)
% Source: https://www.mathworks.com/matlabcentral/answers/58885-creating-a-spherical-matrix

Rx = (dx)/2; Ry = (dy)/2; Rz = (dz)/2; % Three radius
[X,Y,Z] = ndgrid(linspace(-Rx,Rx,dx),linspace(-Ry,Ry,dy),linspace(-Rz,Rz,dz));
R = sqrt((X/Rx).^2 + (Y/Ry).^2 + (Z/Rz).^2);
binary_ellipsoid = zeros(size(X));
binary_ellipsoid(R <= 1 ) = 1; % Assign 1 for ellipsoid, 0 for complementary volume
% Remove one voxel tip, if any?

% Visualization
%slice(Y,X,Z,binary_ellipsoid,0,0,0);
%axis equal vis3d; colorbar;
%imagesc(binary_ellipsoid(:,:,1))
end

% Angle Function 
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

%% Generated Microstructure tab 7

    function [Volume_fraction,Tortuosity_factor] = simpleevaluation_microstructure(number_phase,voxel_number,microstructure3D,choices)
        for current_phase = 1:1:number_phase
            % Volume Fraction
            Volume_fraction(current_phase+1) = sum(sum(sum( microstructure3D.phase == phase(current_phase).code))) / voxel_number;
        end
        Volume_fraction(1) = sum(sum(sum( microstructure3D.phase ==0))) / voxel_number; 
        %Tau factor
        
        if choices(1)
            binary_phase = zeros(size(microstructure3D.phase));
            binary_phase(microstructure3D.phase==0)=1;
            for current_direction = 1:3
                if current_direction==1
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                    if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
                        Tau_factor_result.Tau_W1.Tau = NaN;
                    end
                    Tortuosity_factor(current_direction)=Tau_factor_result.Tau_W1.Tau;
                elseif current_direction==2
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                    if Tau_factor_result.Tau_W2.Tau == Inf
                        Tau_factor_result.Tau_W2.Tau = NaN;
                    end
                    Tortuosity_factor(current_direction)=Tau_factor_result.Tau_W2.Tau;
                elseif current_direction==3
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                    if Tau_factor_result.Tau_W3.Tau == Inf
                        Tau_factor_result.Tau_W3.Tau = NaN;
                    end
                    Tortuosity_factor(current_direction)=Tau_factor_result.Tau_W3.Tau;
                end
            end
        else
            Tortuosity_factor(1)=0; Tortuosity_factor(2)=0; Tortuosity_factor(3)=0;
        end
    end

end

