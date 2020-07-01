function Microstructure_analysis_toolbox_mainmenu
%Main menu of the microstructure analysis toolbox

%% MEMORY
clear global % Clear global variables
clearvars % Clear variables
%close all % Close figure
clc % clean command window

%% FIGURE RENDERING
% See https://www.mathworks.com/help/matlab/ref/opengl.html for more information
%opengl 'software' % Software figure rendering. More stable, use it if MATLAB experiences figures-related crash
opengl 'hardware' % Hardware (OpenGL) figure rendering.

%%
%% MAIN MENU INTERFACE
%%

%% TEXT OPTIONS
font_name_GUI ='Times New Roman'; % Font
font_size_large_GUI=20; % Font size
font_size_medium_GUI=16;
font_size_small_GUI=14;

% color
figure_color=[0.94 0.94 0.94];
background_main_title = [238	201	0]/255;
ForegroundColor_main_title = [0 0 0]; 
ForegroundColor_introductiontext = [0 0 0]; 
background_pushbutton = [0 139 139]/255;
ForegroundColor_pushbutton = [1 1 1]; 
background_smallpushbutton = [1 1 1];
ForegroundColor_smallpushbutton = [0 0 0]; 


%% CREATE FIGURE
main_menu = figure; % Create figure
main_menu.Name= 'MATBOX Microstructure analytical toolbox, Main menu'; % Set name
main_menu.NumberTitle='off'; % Remove number from title name
main_menu.Color=figure_color; % Background colour
main_menu.MenuBar='none'; % Remove menubar and toolbar
main_menu.ToolBar='none';
main_menu.Units='normalized'; % Set unit
main_menu.Position=[0.25 0.25 0.3 0.65]; % Set Position

%% CREATE FIGURE OBJECT

% Title
Text_main_menu_titlebackground = uicontrol('Parent', main_menu, 'Style', 'text','Units','normalized','Position',[0 0.9 1 0.1],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_main_title,'ForegroundColor',background_main_title,...
    'String','');
Text_main_menu_title = uicontrol('Parent', main_menu, 'Style', 'text','Units','normalized','Position',[0 0.925 1 0.05],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_main_title,'ForegroundColor',ForegroundColor_main_title,'FontWeight','bold',...
    'String','NREL Microstructure analysis toolbox');

% Introduction text
Text_main_menu_introductiontext = uicontrol('Parent', main_menu, 'Style', 'text','Units','normalized','Position',[0 0.825 1 0.05],...
    'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'ForegroundColor',ForegroundColor_introductiontext,...
    'String','Select your activity');

% Create line
annotation(main_menu,'line',[0.1 0.9],[0.82 0.82]);

% Push button: microstructure generation
Button_microstructure_generation = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', '>>> Microstructure generation <<<',...
    'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.72 0.8 0.075],...
    'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_generation});
str = sprintf('Numerical generation of microstructures\nbased on probability (stochastic) method.');
Button_microstructure_generation.TooltipString = str;

% Push button: Microstructure filtering and segmentation
Button_microstructure_filteringsegmentation = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', '>>> Filtering and segmentation <<<',...
    'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.61 0.8 0.075],...
    'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_filteringsegmentation});
str = sprintf('Calculate image quality and apply filter to reduce noise to help further segmentation.\nApply basic segmentation (manual and Otsu) and calculate microstructure parameters sensitivity with the segmentation threshold.');
Button_microstructure_filteringsegmentation.TooltipString = str;

% Push button: Microstructure characterisation
Button_microstructure_characterisation = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', '>>> Microstructure characterization <<<',...
    'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.5 0.8 0.075],...
    'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_characterization});
str = sprintf('Calculate classic properties for battery and fuel cells macromodels: volume fraction, connectivity, specific surface area, triple phase boundary, particle diameter and tortuosity factor.\nChoose among a set of different numerical methods to finely evaluate properties value.\nPerform Representative Volume Element (RVE) analysis and image resolution sensitivity analysis to ckeck the calculation and estimate the error.\nCalculate morphology parameters: particle shpericity and particle elongation.\nCalculate geometric tortuosity with the graph representation of the microstructure.\nGenerate python files to be used with FEniCS Finite Element Solver for homogeneisation calculation.');
Button_microstructure_characterisation.TooltipString = str;

% Push button: microstructue visualisation
Button_microstructure_visualisation = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', '>>> Microstructure and results visualisation <<<',...
    'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.39 0.8 0.075],...
    'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_visualisation});
str = sprintf('Use to visualize volumes and results from the microstructure characterisation');
Button_microstructure_visualisation.TooltipString = str;

% Push button: correlation
Button_microstructure_correlation = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', '>>> Properties correlation <<<',...
    'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.28 0.8 0.07],...
    'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_correlation});
str = sprintf('Plot microstructure parameters obtain from different volumes as a function of each other to estimate their correlation');
Button_microstructure_correlation.TooltipString = str;

% Push button: Meshing
Button_microstructure_mesh = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', '>>> Create mesh for FEM <<<',...
    'FontSize',font_size_medium_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.17 0.8 0.075],...
    'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_meshing});
str = sprintf('Build tetrahedron-based mesh for one phase or a full cell.\nApply smoothing algorithm and control the mesh density.');
Button_microstructure_mesh.TooltipString = str;

% Create line
annotation(main_menu,'line',[0.1 0.9],[0.12 0.12]);

% Push button: About
Button_about = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', 'About',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.025 0.15 0.065],...
    'BackgroundColor',background_smallpushbutton,'ForegroundColor',ForegroundColor_smallpushbutton,'Callback',{@Callback_about});

% Push button: Repository
Button_repo = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', 'Open repository',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.3 0.025 0.325 0.065],...
    'BackgroundColor',background_smallpushbutton,'ForegroundColor',ForegroundColor_smallpushbutton,'Callback',{@Callback_openrepo});

% Push button: Documentation
Button_documentation = uicontrol('Parent', main_menu, 'Style', 'pushbutton', 'String', 'Documentation',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.675 0.025 0.225 0.065],...
    'BackgroundColor',background_smallpushbutton,'ForegroundColor',ForegroundColor_smallpushbutton,'Callback',{@Callback_documentation});


%% CALLBACK FUNCTIONS

function Callback_microstructure_generation(~,~)
    % Sub menu
    Generation_menu = figure; % Create figure
    Generation_menu.Name= 'Microstructure genration menu'; % Set name
    Generation_menu.NumberTitle='off'; % Remove number from title name
    Generation_menu.Color=figure_color; % Background colour
    Generation_menu.MenuBar='none'; % Remove menubar and toolbar
    Generation_menu.ToolBar='none';
    Generation_menu.Units='normalized'; % Set unit
    Generation_menu.Position=[0.55 0.65 0.3 0.25]; % Set Position
    
    % Title
    Text_generation_menu_titlebackground = uicontrol('Parent', Generation_menu, 'Style', 'text','Units','normalized','Position',[0 0.8 1 0.2],...
        'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'BackgroundColor',background_main_title,'ForegroundColor',background_main_title,...
        'String','');
    Text_generation_menu_title = uicontrol('Parent', Generation_menu, 'Style', 'text','Units','normalized','Position',[0 0.84 1 0.1],...
        'FontName',font_name_GUI,'FontSize',font_size_medium_GUI,'BackgroundColor',background_main_title,'ForegroundColor',ForegroundColor_main_title,'FontWeight','bold',...
        'String','Microstructure generation menu');
    % Introduction text
    Text_generation_menu_introductiontext = uicontrol('Parent', Generation_menu, 'Style', 'text','Units','normalized','Position',[0 0.65 1 0.1],...
        'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'ForegroundColor',ForegroundColor_introductiontext,...
        'String','Select generation task');
    
    % Create line
    annotation(Generation_menu,'line',[0.1 0.9],[0.625 0.625]);
    
    % Push button: ellipsoid based generation
    Button_microstructure_generation_sphere = uicontrol('Parent', Generation_menu, 'Style', 'pushbutton', 'String', 'Ellipsoids-based',...
        'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.1 0.4 0.35 0.15],...
        'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_generation_ellipsoids});
    str = sprintf('Numerical generation of microstructures\nbased on probability (stochastic) method, using ellipsoids particle shape.');
    Button_microstructure_generation_sphere.TooltipString = str;
    
    % Push button: additives generation
    Button_microstructure_generation_additives = uicontrol('Parent', Generation_menu, 'Style', 'pushbutton', 'String', 'Additives',...
        'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.55 0.4 0.35 0.15],...
        'BackgroundColor',background_pushbutton,'ForegroundColor',ForegroundColor_pushbutton,'Callback',{@Callback_microstructure_generation_additives});
    str = sprintf('Numerical generation of additive phase, such as carbon-binder for battery electodes\non top of an existing pore and solid domain.');
    Button_microstructure_generation_additives.TooltipString = str;    
    
    
end

function Callback_microstructure_generation_ellipsoids(~,~)
    % Run program
    microstructure_generation_ellipsoids_GUI
end

function Callback_microstructure_generation_additives(~,~)
    % Run program
    Microstructure_generation_additives
end

function Callback_microstructure_filteringsegmentation(~,~)
    % Run program
    Segmentation
end

function Callback_microstructure_characterization(~,~)
    % Run program
    microstructure_characterization_GUI
end

function Callback_microstructure_visualisation(~,~)
    % Run program
    microstructure_visualization_GUI
end

function Callback_microstructure_correlation(~,~)
    % Run program
    microstructure_correlation_GUI
end

function Callback_microstructure_meshing(~,~)
    % Run program
    microstructure_meshing_GUI
end


function Callback_about(~,~)
    web('https://www.nrel.gov/transportation/data-tools.html');
end
function Callback_openrepo(~,~)
    web('https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox');
end
function Callback_documentation(~,~)
    path = matlab.desktop.editor.getActiveFilename; % Path of active file
    higherlevelfolder = extractBetween(path,path(1:5),'MATBOX_Microstructure_analysis_toolbox\','Boundaries','inclusive');
    if ispc
        documentation_path = [char(higherlevelfolder) 'Documentation\MATBOX_Microstructure_analysis_toolbox_documentation.pdf'];
    else
        documentation_path = [char(higherlevelfolder) 'Documentation/MATBOX_Microstructure_analysis_toolbox_documentation.pdf'];
    end
    if exist(documentation_path,'file')
        open(documentation_path);
    else
        disp 'MATLAB did not find the file NREL_Microstructure_analysis_toolbox_documentation.pdf'. 
        disp 'Default location is \Microstructure_analysis_toolbox\Documentation\';
    end
end


end