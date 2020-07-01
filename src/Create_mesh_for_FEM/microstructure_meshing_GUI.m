function microstructure_meshing_GUI
%microstructure_meshing_GUI Graphic user interface for 3D microstructure meshing for a half or full-cell battery.
%Author: Francois L.E. Usseglio-Viretta, National Renewable Energy Laboratory.
%Program contains third-pary code (Iso2mesh). See http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home
%   You need to download it and add it to your Matlab path.
%LICENCSES are included within the code, in the section LICENSES.
%To use it: type microstructure_meshing_GUI or run it from the main GUI Microstructure_analysis_tool.m
%Minimal instructions to use microstructure_meshing_GUI:
%Setup each tab from top to bottom:
%- "Folder and save options" tab: Select your save folder ("Click to select main save folder" button). Figures and data will be saved there.
%                                 You can keep the default saving options.                             
%- "Microstructures" tab: Select the microstructure to setup (Left electrode, Separator, or Right electrode using the popup menu).
%                         Click on "New/reset electrode".
%                         Choose between imported heterogeneous medium (you will need a tif image) and homogeneous medium.
%                         Follow instructions.
%                         Once done, click on "Save electrode".
%                         Repeat the procedure for the other two microstructures. Note you can skip the separator if you intend to generate it in the next tab.
% "Generate separator" tab: Do this tab only if you did not setup the separator in the previous tab.
%                           Select the separator geometry to generate using the popup menu.
%                           Click on "New/reset separator".
%                           Fill parameters.
%                           Click on "Click to generate separator".
%                           Once done, click on "Click to save separator".
% "Cell options" tab: Enter voxel size in micrometers.
%                     If texts are colored in red in the section "1) Check dimensions are matching", check the box "crop in-plane dimensions of larger domains".
%                     Erosion + dilatation is optional and should be used only if the mesh creation fails.
%                     Convert each phase in a unique cluster and clean it: highly recommended.
%                     Tables below "3) Choose phase id that will be used in the cell for each microstructure": enter values for the coloumn "Cell id" (different integer values, >0, not =0).
%                     This is basically a re-assignment. If you enter correct values, you will be able to assin phase in the section "4) Assign phase".
%                     For each entry of the "Cell id" popup menu, select the correct entry of the "Cell assignment" popup menu. Then click on "Validate assignment".
%                     Once all phases have been assigned, click on "Assemble cell", then "Save cell and go on".
% "Iso2mesh options" tab: Select methods in the popup menu (recommended: lowpass)
%                         Enter parameters (you can keep the default values).
%                         See http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/FunctionList for details.
% "Create mesh" tab: Click on "Click to generate mesh".
%                    In case of a Iso2mesh error message: "two surfaces...are found intersecting each other".
%                       See http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/FAQ#I_am_getting_a_Two_subfaces_are_found_intersecting_each_other_error_what_should_I_do
%                       Possible solution: Iso2mesh tab: use "cgalmesh", reduce "radbound", reduce "Iteration number".
%                                          Cell options tab: Use "Erosion + dilatation" and "Convert each phase in a unique cluster and clean it".
%                    You will need a lot of RAM memory: a workstation is recommended. If windows fails due to a memory error, try with unix (it does better).
%                    Click on Inspect regions. It is likely one region will be quite incorrect (usually the region 1). Note that there is some RNG in action here, result will differ at each new try.
%                    Read the instructions below "Inspect regions" to decide how to correct the regions.
%                    For each entry of the "Region number" popup menu, select the correct entry of the "Cell assignment" popup menu. Then click on "Validate assignment".
%                    Click on "Assign solid and electrolyte".
%                    Optional: Click on "Analyse mesh (optional)"
%                    Click on "Save meshes".
%Note: if buttons are grey out (i.e., enable,'off'), it means you need to enter other parameters elsewhere before moving forward.
%Once done, to re-create the volumetric mesh of * you need to use the two files:
%   /Mesh_data/Nodes_*.mat: a n*3 array. x,y,z coordinates of the n verteces of *. Coordinates are adimensional (1 = voxel length)
%   /Mesh_data/Tetrahedron_*.mat: a m*4 array. Verteces index of the m tetrahedron cell.
%   Example with FEniCS version 2017.2.0:
% 	def create_mesh_inserial(node_,tetrahedroncell_): # Function defintion (python)
% 		number_vertices=len(node_) # Number of vertices
% 		number_cells=len(tetrahedroncell_) # Number of cells
% 		mesh = Mesh() # Create an empty mesh object
% 		editor = MeshEditor() # Select the editor
% 		editor.open(mesh, 'tetrahedron', 3, 3, 1) # Open the mesh (mesh to edit, cell type, topological dimension, geometrical dimension, polynomial degree)
% 		editor.init_vertices(number_vertices) # Set number of vertices
% 		editor.init_cells(number_cells) # Set number of cells
% 		# Set vertex location (vertex id,x,y)
% 		for vertex_id in range(number_vertices):
% 			x_ = node_[vertex_id][0] # Vertex coordinate 0
% 			y_ = node_[vertex_id][1] # Vertex coordinate 1
% 			z_ = node_[vertex_id][2] # Vertex coordinate 2
% 			vertex_id=int(vertex_id) # The vertex id must be an integer
% 			editor.add_vertex(vertex_id, x_, y_, z_) # Add the vertex
% 		# Set cell (cell id, vertices id), for tetrahedrons
% 		for cell_id in range(number_cells):
% 			v0 = tetrahedroncell_[cell_id][0] # Vertex indice 0
% 			v1 = tetrahedroncell_[cell_id][1] # Vertex indice 1
% 			v2 = tetrahedroncell_[cell_id][2] # Vertex indice 2
% 			v3 = tetrahedroncell_[cell_id][3] # Vertex indice 3
% 			cell_id=int(cell_id); v0=int(v0); v1=int(v1);	v2=int(v2); v3=int(v3) # The ids must be integers
% 			editor.add_cell(cell_id, v0, v1, v2, v3) # Add the cell
% 		editor.close() # Close the editor
% 		return mesh
%   Node_Electrolyte = spio.loadmat(Nodes_Electrolyte.mat)['node_region'] # node_region is the variable name saved in Matlab
%   Tetrahedron_Electrolyte = spio.loadmat('Tetrahedron_Electrolyte.mat')['elem_'] # elem_ is the variable name saved in Matlab
%	mesh_electrolyte =	create_mesh_inserial(Node_Electrolyte,Tetrahedron_Electrolyte) # Call the function

%% TO-DO LIST
% WIP: code and licences not up to date.
%      code not yet ready to be published - currently only for NREL internal use (i.e. me)

% % Must-do before publishing
% tab About and Licence
% Remove these hard-coded values (used for morphology opening): Imported heterogeneous domains must be 0 for pore, 1 for solid. 
%         celldata(1).codesolid = 1;
%         celldata(2).codesolid = 1;
%         celldata(3).codesolid = 1;
%         celldata(1).codepore = 0;
%         celldata(2).codepore = 0;
%         celldata(3).codepore = 0;
% Morphology opening may have its own tab.
% Check mesh with isolated clusters.
% Fix sub_markers
% Fix hard-coded value of section, with a better integration with the main code %% FROM JEFFEREY ALLEN
% Fix color_phase
% Clean code (re factoring)

% % Can be added after publishing
% Save volume fraction, especially before/after each step of morphology opening
% Clarify saving options
% Visualize function: remove redundant functions (keep only one)
% Improve Iso2mesh region correction (automatic?)
% Cell visualization: assign color, create video, take screenshot.
% Allow saving the cell setup, and re-load it
% Mesh subvolumes to perform subsequent RVE analysis
% Reduce number of global variables (or better remove all of them)

%% MICROSTRUCTURE_MESHING_GUI
% Author: Francois L.E. Usseglio-Viretta
%         National Renewable Energy Laboratory
% ...
%         Jefferey Allen: contribution (section %% FROM JEFFEREY ALLEN)
%         National Renewable Energy Laboratory
% ...

%% THIRD-PARTY CODE: ISO2MESH
% Mesh generation performed in this code is performed with the third-party toolbox Iso2mesh.
% Functions from Iso2mesh used in this code: 
% v2s, meshconn, smoothsurf, surf2mesh, plotmesh, plottetra
% If you are interested in the mesh creation, look at the function mesh_generation_withIso2mesh.m

% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/README
% "
% Author: Qianqian Fang <q.fang at neu.edu>
%   Department of Bioengineering
%   Northeastern University
%   360 Huntington Ave, Boston, MA 02115
% Version: 1.7.9 (Deviled Egg - beta)
% License: GPL v2 or later (see COPYING)
%   (this license does not cover the binaries under the bin/
%   directory, see Section III for more details)
% URL: http://iso2mesh.sf.net
% "
%
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download#Important_Note_on_Licenses
% "
% The iso2mesh toolbox is licensed under GNU Public License (GPL). However, it included a number of external tools (under bin/ directory) to perform certain functionality by means of pipelines. These external tools are subjected to their upstream licenses and are not necessarily GPL (or GPL compatible). The complete list of these external commands and their author info and licenses can be found from the README page.
% Particularly, the 3D meshing tool, tetgen, is licensed under a non-free license: it can be freely used, modified, redistributed only for research and academic purposes, any commercial utility of tetgen requires a permission from its original author. iso2mesh calls tetgen in the background to produce 3D mesh, that means if anyone needs to uses the 3D mesh produced by iso2mesh in a commercial product, you MUST contact the author of tetgen to get permission. Processing binary images and produce surfaces are not subject to this limitation.
% In additional to the licenses, if you use this tool in your research, we are greatly appreciated if you can add iso2mesh to your references:
% Qianqian Fang, iso2mesh: a matlab-based 3D tetrahedral mesh generator, URL: http://iso2mesh.sourceforge.net, 2008
%  or
% Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary and gray-scale images," 
% Proceedings of IEEE International Symposium on Biomedical Imaging 2009, pp. 1142-1145, 2009
% If you used cgalsurf or cgalmesh options in your mesh preparation, you should also acknowledge CGAL publications. If you generated your mesh with the tetgen module in iso2mesh instead of CGAL mesher, you may need to acknowledge tetgen in your publications, you can find more references from this link.
% "
%
% Iso2mesh interacts with a number external meshing tools to perform the essential functionalities. These tools are listed below:
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/README#Acknowledgement

%%
%% GUI DISPLAY OPTIONS
%%

% Tab description
description_tab_fromleft = 0.025; % Position
description_tab_frombottom = 0.925;
description_tab_xlenght = 1-2*description_tab_fromleft;
description_tab_ylenght = 0.045;
background_description_tab = [0 0.5 0]; % Background and text color
ForegroundColor_description_tab = [1 1 1]; 
background_tab = [1 1 1];
% Text
font_name_GUI ='Times New Roman'; % Font
font_size_small_GUI =10; % Font size
font_size_large_GUI =14; % Font size

color_phase=[ % Color
    0 114 189;
    217 83 25
    237 177 32];
tmp=randi(255,100,3);
color_phase=[color_phase;tmp];
color_phase=color_phase/255; % Normalized for Matlab

%% GLOBAL VARIABLES
global Domain_size array_initial array_ROI array_orientation array_resize...
       Phase_code number_phase  min_inplane cell_array options node elem face str_region...
       node_leftelectrode elem_leftelectrode face_leftelectrode...
       node_electrolyte elem_electrolyte face_electrolyte...
       node_rightelectrode elem_rightelectrode face_rightelectrode...
       Fig_visualizemicrostructure Fig_visualizecell Fig_visualizeseparator...
       voxel_size_um separator_array region_
   
%% HARD-CODED OPTIONS
% Folders
if ispc
    options.save.viewalldomainfolder = 'View_alldomains\';
    options.save.vieweachdomainfolder = 'View_eachdomains\';
    options.save.viewregion = 'View_eachregions\';
    options.save.meshdatafolder = 'Mesh_data\';
    options.save.stlfolder = 'Stl\';
    options.save.analysemeshfolder = 'Analyse_meshes\';
    options.domainfolder_left = 'Left_electrode\';
    options.domainfolder_electrolyte = 'Electrolyte\';
    options.domainfolder_right = 'Right_electrode\';
else
    options.save.viewalldomainfolder = 'View_alldomains/';
    options.save.vieweachdomainfolder = 'View_eachdomains/';
    options.save.viewregion = 'View_eachregions/';
    options.save.meshdatafolder = 'Mesh_data/';
    options.save.stlfolder = 'Stl/';
    options.save.analysemeshfolder = 'Analyse_meshes/';
    options.domainfolder_left = 'Left_electrode/';
    options.domainfolder_electrolyte = 'Electrolyte/';
    options.domainfolder_right = 'Right_electrode/';
end
% Analyse meshes
options.analysemesh.left = true;
options.analysemesh.electrolyte = true;
options.analysemesh.right = true;
% Find common verteces at the interface
options.tolerance_interface=1e-6;
% Font
options.fontname = 'Times New Roman';
% Smoothing of cumulative and distribution functions
options.smooth_cumulative_fct = true;
options.minimum_array_length = 5;
options.number_point = 250;
options.moving_rangeratio = 0.05;
options.moving_range = 0;
options.enforce_samelength = false;
options.origin = 'symmetrical';
options.boundary_behavior = 'keep origin';
% Figure 
options.grid = 'on';
options.minorgrid = 'on';
  
%%
%% CREATE MAIN GUI
%%

%% MAIN FIGURE
main_figure = figure; % Create figure
main_figure.Name= 'Microstructure meshing'; % Set name
main_figure.NumberTitle='off'; % Remove number from title name
main_figure.Color='white'; % Background colour
main_figure.MenuBar='none'; % Remove menubar and toolbar
main_figure.ToolBar='none';
main_figure.Units='normalized'; % Set unit
main_figure.Position=[0.2 0.15 0.6 0.75]; % Set Position

%% USER INTERFACE WITH TABBED PANELS
table_group_1 = uitabgroup('Parent', main_figure); % Create tabgroup
table_group_1.TabLocation='Left'; % Set panels location
% Create tabbed panels
tab_save = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Folder and save options');
tab_microstructures = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Microstructures');
tab_separator = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Generate separator');
tab_cell = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Cell options');
tab_Iso2mesh = uitab('Parent', table_group_1,'BackgroundColor',background_tab,'Title', 'Iso2mesh options');
tab_generate = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'Create mesh');
tab_about = uitab('Parent', table_group_1,'BackgroundColor',background_tab, 'Title', 'About');
% Tab colors
tab_save.ForegroundColor  ='b';
tab_microstructures.ForegroundColor  ='b';
tab_separator.ForegroundColor  ='b';
tab_cell.ForegroundColor  ='b';
tab_Iso2mesh.ForegroundColor  ='b';
tab_generate.ForegroundColor  ='r';
tab_about.ForegroundColor  ='k';
% Preselected tab
table_group_1.SelectedTab = tab_save;

%%
%% "FOLDER AND SAVE OPTIONS" TAB
%%

%% PRE-ASSINGNED PARAMETERS
% % Default values
% Save 3D array
options.save.microstructure_tif = 1;
options.save.microstructure_mat = 1;
options.save.cell_tif = 1;
options.save.cell_mat = 1;
% Save mesh data
options.save.microstructuremesh = 1;
options.save.microstructurestl = 1;
options.save.microstructurebinarystl = 0;
options.save.cellmesh = 1;
% Interface data
options.save.nodeinterface_fig = 1;
options.save.nodeinterface_png = 1;
options.save.nodeinterface_mat = 1;
% 3D view region
options.save.regions_fig = 1;
options.save.regions_png = 1;
% 3D view all domain
options.save.initialsurfacemesh = 0;
options.save.smoothedsurfacemesh = 0;
options.save.volumetricmesh = 1;
options.save.volumetricmesh_slices = 1;
options.save.alldomain_fig = 0;
options.save.alldomain_png = 1;
options.save.alldomain_slice_fig = 1;
options.save.alldomain_slice_png = 1;
% 3D view each domain
options.save.dommain_left = 1;
options.save.dommain_electrolyte = 1;
options.save.dommain_right = 1;
options.save.dommain_transparent = 1;
options.save.dommain_plain = 1;
options.save.dommain_light = 0;
options.save.dommain_lightnomesh = 1;
options.save.dommain_lightnomesh_jet = 0;
options.save.dommain_slice = 0;
options.save.dommain_fig = 0;
options.save.domain_png = 1;

%% GUI
% Title
Text_tab_save = uicontrol('Parent', tab_save, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Folder and I/O options');

% Select folder
Button_select_mainfolder = uicontrol('Parent', tab_save, 'Style', 'pushbutton', 'String', '>>>  Click to select main save folder  <<<',...
    'FontSize',font_size_large_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.87 0.4 0.04],...
    'visible','on','enable','on','Callback',{@select_savefolder_Callback});
Text_Folder_statut = uicontrol('Parent', tab_save, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Save folder location: NOT DEFINED','ForegroundColor','r',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.81 0.5 0.05]);
 
% 3D array 
Text_save_3Darray = uicontrol('Parent', tab_save, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Save arrays data (hold Ctrl to select multiple choices)',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.775 0.5 0.04]);
listbox_save_3Darray = uicontrol('Parent', tab_save, 'Style', 'listbox', 'Max',4, 'String', {'Each microstructure (.tif)','Each microstructure (.mat)','Cell (.tif)','Cell (.mat)'},...
    'Value',[1 2 3 4],'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.7 0.3 0.075],'visible','on','enable','on','Callback',{@listbox_save_3Darray_Callback});

% Mesh data
Text_save_meshdata = uicontrol('Parent', tab_save, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Save mesh data (hold Ctrl to select multiple choices)',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.625 0.5 0.04]);
listbox_save_meshdata = uicontrol('Parent', tab_save, 'Style', 'listbox', 'Max',4, 'String', {'Each microstructure (.mat)','Each microstructure (.stl)','Each microstructure (.binarystl)','Cell (.mat)'},...
    'Value',[1 2 4],'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.55 0.3 0.075],'visible','on','enable','on','Callback',{@listbox_save_meshdata_Callback});

% Interface data
Text_save_interfacedata = uicontrol('Parent', tab_save, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Save interface data (hold Ctrl to select multiple choices)',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.475 0.5 0.04]);
listbox_save_interfacedata = uicontrol('Parent', tab_save, 'Style', 'listbox', 'Max',4, 'Min',0, 'String', {'Vertece coordinates at the interface solid-electrolyte'},...
    'Value',[1],'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.4 0.4 0.075],'visible','on','enable','on','Callback',{@listbox_save_interfacedata_Callback});

% View
Text_save_view = uicontrol('Parent', tab_save, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Save 3D view (hold Ctrl to select multiple choices)',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.325 0.5 0.04]);
listbox_save_view = uicontrol('Parent', tab_save, 'Style', 'listbox', 'Max',99, 'Min',0, 'String', {'Iso2mesh regions (.fig)','Iso2mesh regions (.png)','Cell: surface mesh','Cell: smoothed surface mesh','Cell: volumetric mesh (w/ transparency)','Cell (.fig)','Cell (.png)','Cell: slices','Cell slices (.fig)','Cell slices (.png)',...
    'Domain: left electrode (solid)','Domain: Electrolyte','Domain: Right electrode (solid)','Domain: plain image','Domain: transparent image','Domain: image w/ lighting','Domain: image w/o edges','Domain: image w/ jet colormap','Domain: slices','Domain: .fig','Domain: .png'},...
    'Value',[1 2 5 7 8 9 10 11 12 13 14 15 17 21],'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.025 0.4 0.3],'visible','on','enable','on','Callback',{@listbox_save_view_Callback});
 
%% CALLBACKS
     function select_savefolder_Callback (~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select location where the save folder will be created';
        % Open dialog box to choose folder path
        resultfolder_location = uigetdir(matlabroot,str_dialogbox);
        if resultfolder_location==0
            % User clicked cancel button or closed the dialog box
            set(Text_Folder_statut,'String','Save folder location: NOT DEFINED','ForegroundColor','r');
        else
            % Update savefolder
            if ispc
                options.save.mainfolder = [resultfolder_location '\'];
            else
                options.save.mainfolder = [resultfolder_location '/'];
            end
            set(Text_Folder_statut,'String',['Save folder location: ' options.save.mainfolder],'ForegroundColor','k');
        end
    end

    function listbox_save_3Darray_Callback(~,~)
        choices = listbox_save_3Darray.Value;
        if sum(choices==1)==1
            options.save.microstructure_tif=1;
        else
            options.save.microstructure_tif=0;
        end        
        if sum(choices==2)==1
            options.save.microstructure_mat=1;
        else
            options.save.microstructure_mat=0;
        end
        if sum(choices==3)==1
            options.save.cell_tif=1;
        else
            options.save.cell_tif=0;
        end
        if sum(choices==4)==1
            options.save.cell_mat=1;
        else
            options.save.cell_mat=0;
        end        
    end

    function listbox_save_meshdata_Callback(~,~)
        choices = listbox_save_meshdata.Value;
        if sum(choices==1)==1
            options.save.microstructuremesh=1;
        else
            options.save.microstructuremesh=0;
        end
        if sum(choices==2)==1
            options.save.microstructurestl=1;
        else
            options.save.microstructurestl=0;
        end
        if sum(choices==3)==1
            options.save.microstructurebinarystl=1;
        else
            options.save.microstructurebinarystl=0;
        end
        if sum(choices==4)==1
            options.save.cellmesh=1;
        else
            options.save.cellmesh=0;
        end
    end

    function listbox_save_interfacedata_Callback(~,~)
        choices = listbox_save_interfacedata.Value;
        if sum(choices==1)==1
            options.save.nodeinterface_fig = 1;
            options.save.nodeinterface_png = 1;
            options.save.nodeinterface_mat = 1;
        else
            options.save.nodeinterface_fig = 0;
            options.save.nodeinterface_png = 0;
            options.save.nodeinterface_mat = 0;
        end
    end

    function listbox_save_view_Callback(~,~)
        choices = listbox_save_view.Value;
        if sum(choices==1)==1; options.save.regions_fig = 1; else; options.save.regions_fig = 0; end
        if sum(choices==2)==1; options.save.regions_png = 1; else; options.save.regions_png = 0; end
        if sum(choices==3)==1; options.save.initialsurfacemesh = 1; else; options.save.initialsurfacemesh = 0; end
        if sum(choices==4)==1; options.save.smoothedsurfacemesh = 1; else; options.save.smoothedsurfacemesh = 0; end
        if sum(choices==5)==1; options.save.volumetricmesh = 1; else; options.save.volumetricmesh = 0; end
        if sum(choices==6)==1; options.save.alldomain_fig = 1; else; options.save.alldomain_fig = 0; end
        if sum(choices==7)==1; options.save.alldomain_png = 1; else; options.save.alldomain_png = 0; end
        if sum(choices==8)==1; options.save.volumetricmesh_slices = 1; else; options.save.volumetricmesh_slices = 0; end
        if sum(choices==9)==1; options.save.alldomain_slice_fig = 1; else; options.save.alldomain_slice_fig = 0; end
        if sum(choices==10)==1; options.save.alldomain_slice_png = 1; else; options.save.alldomain_slice_png = 0; end
        if sum(choices==11)==1; options.save.dommain_left = 1; else; options.save.dommain_left = 0; end
        if sum(choices==12)==1; options.save.dommain_electrolyte = 1; else; options.save.dommain_electrolyte = 0; end
        if sum(choices==13)==1; options.save.dommain_right = 1; else; options.save.dommain_right = 0; end
        if sum(choices==14)==1; options.save.dommain_plain = 1; else; options.save.dommain_plain = 0; end
        if sum(choices==15)==1; options.save.dommain_transparent = 1; else; options.save.dommain_transparent = 0; end
        if sum(choices==16)==1; options.save.dommain_light = 1; else; options.save.dommain_light = 0; end
        if sum(choices==17)==1; options.save.dommain_lightnomesh = 1; else; options.save.dommain_lightnomesh = 0; end
        if sum(choices==18)==1; options.save.dommain_lightnomesh_jet = 1; else; options.save.dommain_lightnomesh_jet = 0; end
        if sum(choices==19)==1; options.save.dommain_slice = 1; else; options.save.dommain_slice = 0; end
        if sum(choices==20)==1; options.save.dommain_fig = 1; else; options.save.dommain_fig = 0; end
        if sum(choices==21)==1; options.save.domain_png = 1; else; options.save.domain_png = 0; end
    end

%%
%% "MICROSTRUCTURES" TAB
%%

%% PRE-ASSIGNED PARAMETERS
Fig_visualizemicrostructure=0;
choice_cell = ["Left electrode" "Separator" "Right electrode"];
choice_representation = [" " "Imported heterogenous medium" "Homogeneous medium"];
str_geometry = ' ';
initial_length = 0;
% 1,2,3 left electrode, electrolyte, right electrode
celldata(1).matlabname = 'Left electrode'; % Used ?
celldata(2).matlabname = 'Separator'; % Used ?
celldata(3).matlabname = 'Right electrode'; % Used ?
celldata(1).exist = false;
celldata(2).exist = false;
celldata(3).exist = false;
% 4 full or half-cell
celldata(4).exist = false;
celldata(4).array = -1;
celldata(4).iso2mesh_options_set = false;
% WARNING: Hard-coded
celldata(1).codesolid = 1;
celldata(2).codesolid = 1;
celldata(3).codesolid = 1;
celldata(1).codepore = 0;
celldata(2).codepore = 0;
celldata(3).codepore = 0;

%% GUI
% Title
Text_tab_microstructures = uicontrol('Parent', tab_microstructures, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Microstructures');

% New/Save geometry
pushebutton_microstructures_new = uicontrol('Parent', tab_microstructures,'Style', 'pushbutton','Value',0,'String','New/reset electrode','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.87 0.2 0.04],'Callback', @newelectrode_Callback);
pushebutton_microstructures_save = uicontrol('Parent', tab_microstructures,'Style', 'pushbutton','Value',0,'String','Save electrode','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.23 0.87 0.2 0.04],'Callback', @saveelectrode_Callback);

% Choose electrode to mesh
Popupmenu_chooseelectrode = uicontrol('Parent', tab_microstructures, 'Style', 'popupmenu', 'String', choice_cell,...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.82 0.4050 0.04]);
Text_ = uicontrol('Parent', tab_microstructures, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Choose the microstructure to mesh. Once done, save it and choose a new one. If you wish to generate a separator, do only the left and right electrode on this tab.',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.745 0.405 0.08]);

% Choose representation
Text_chooserepresentation = uicontrol('Parent', tab_microstructures, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Choose the electrode representation. Then choose your parameters.',...
    'visible','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.5 0.86 0.475 0.04]);
Popupmenu_choice_representation = uicontrol('Parent', tab_microstructures, 'Style', 'popupmenu', 'String', choice_representation,...
    'visible','off','enable','off','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.5 0.83 0.475 0.04],...
    'Callback',{@choice_representation_Callback});

Text_geometry = uicontrol('Parent', tab_microstructures, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String',str_geometry,...
    'visible','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.7 0.5 0.04]);

% % Homogeneous
Text_length_homogeneous = uicontrol('Parent', tab_microstructures, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Length in voxel:',...
    'visible','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.65 0.2 0.04]);
Edit_length_homogenous = uicontrol('Parent', tab_microstructures,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',initial_length,'Visible','on',...
    'visible','off','enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.135 0.66 0.05 0.04],'Callback', @edit_length_homogenous_Callback);

% % Heterogeneous
% Import
Button_import_microstructure = uicontrol('Parent', tab_microstructures, 'Style', 'pushbutton', 'String', '>>>  Click to import volume  <<<',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.65 0.3 0.05],...
    'visible','off','enable','off','Callback',{@import_microstructure_Callback});
Text_volume = uicontrol('Parent', tab_microstructures, 'Style', 'text','visible','off','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','No volume selected',...
     'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.6 0.2 0.04]);

% Volume fraction 
Text_volumefractions = uicontrol('Parent', tab_microstructures, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Volume fractions',...
    'visible','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.57 0.2 0.04]);
% Table volume fractions
table_volumefraction = uitable('Parent', tab_microstructures,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.025 0.47 0.3 0.12]);
table_volumefraction.ColumnName = {'Phase','Initial','ROI','Voxel resize'}; % Column name
table_volumefraction.ColumnEditable = [false false false false]; % Select column editable
table_volumefraction.ColumnWidth = {'auto', 'auto', 'auto', 'auto'}; % Auto width
table_volumefraction.RowName = []; % Remove row name

% ROI Instructions
Text_ROI_instructions= uicontrol('Parent', tab_microstructures, 'Visible','off','Style', 'text','FontWeight','bold','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','1) Select the Region Of Interest (ROI)',...
    'BackgroundColor','w','ForegroundColor','k','HorizontalAlignment','left','Units','normalized','Position', [0.35 0.675 0.4 0.02]);
% Table ROI
table_ROI = uitable('Parent', tab_microstructures,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.35 0.546 0.5 0.12],'CellEditCallback',@cell_table_ROI_Callback);
table_ROI.ColumnName = {'Array axis','Number of voxels','Start','End','Number of voxels (ROI)'}; % Column name
table_ROI.ColumnEditable = [false false true true false]; % Select column editable
table_ROI.ColumnWidth = {'auto', 'auto', 'auto', 'auto', 'auto'}; % Auto width
table_ROI.RowName = []; % Remove row name
% ROI error message
Text_ROI_error= uicontrol('Parent', tab_microstructures, 'Visible','off','Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Region of interest (ROI) is correct',...
    'BackgroundColor','w','ForegroundColor','k','HorizontalAlignment','center','Units','normalized','Position', [0.875 0.55 0.1 0.1]);

% Orientation instructions
Text_orientation_instructions= uicontrol('Parent', tab_microstructures, 'Visible','off','Style','text','FontWeight','bold','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','2) Set orientation',...
     'BackgroundColor','w','ForegroundColor','k','HorizontalAlignment','left','Units','normalized','Position', [0.35 0.51 0.4 0.02]);
% Table orientations
table_orientations = uitable('Parent', tab_microstructures,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.35 0.38 0.3070 0.12]);
table_orientations.ColumnName = {'Cell direction','Array axis','Number of voxels'}; % Column name
table_orientations.ColumnEditable = [false false false]; % Select column editable
table_orientations.ColumnWidth = {'auto', 'auto', 'auto'}; % Auto width
table_orientations.RowName = []; % Remove row name
% Swap direction buttons  
togglebutton_shiftdirection_1_with_2 = uicontrol('Parent', tab_microstructures,'Style', 'togglebutton','Value',0,'Visible','off','String','1 - 2 (Un-swaped)',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.675 0.465 0.15 0.03],'Callback', @Shiftdirection_1_with_2_Callback); 
togglebutton_shiftdirection_1_with_3 = uicontrol('Parent', tab_microstructures,'Style', 'togglebutton','Value',0,'Visible','off','String','1 - 3 (Un-swaped)',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.675 0.425 0.15 0.03],'Callback', @Shiftdirection_1_with_3_Callback);
togglebutton_shiftdirection_2_with_3 = uicontrol('Parent', tab_microstructures,'Style', 'togglebutton','Value',0,'Visible','off','String','2 - 3 (Un-swaped)',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.675 0.385 0.15 0.03],'Callback', @Shiftdirection_2_with_3_Callback);
% Reverse direction buttons
togglebutton_inversedirection_1 = uicontrol('Parent', tab_microstructures,'Style', 'togglebutton','Value',0,'Visible','off','String','1 (Un-changed)',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.825 0.465 0.15 0.03],'Callback', @Inversedirection_1_Callback); 
togglebutton_inversedirection_2 = uicontrol('Parent', tab_microstructures,'Style', 'togglebutton','Value',0,'Visible','off','String','1 (Un-changed)',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.825 0.425 0.15 0.03],'Callback', @Inversedirection_2_Callback); 
togglebutton_inversedirection_3 = uicontrol('Parent', tab_microstructures,'Style', 'togglebutton','Value',0,'Visible','off','String','1 (Un-changed)',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.825 0.385 0.15 0.03],'Callback', @Inversedirection_3_Callback); 

% Voxel size instructions
Text_voxelsize_instructions= uicontrol('Parent', tab_microstructures, 'Visible','off','Style','text','FontWeight','bold','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','3) Set voxel size scaling factor',...
     'BackgroundColor','w','ForegroundColor','k','HorizontalAlignment','left','Units','normalized','Position', [0.35 0.35 0.4 0.02]);
% Table Voxel size
table_voxelresize = uitable('Parent', tab_microstructures,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.35 0.22 0.34 0.12]);
table_voxelresize.ColumnName = {'Cell direction','Number of voxels, before','After'}; % Column name
table_voxelresize.ColumnEditable = [false false false]; % Select column editable
table_voxelresize.ColumnWidth = {'auto', 'auto', 'auto'}; % Auto width
table_voxelresize.RowName = []; % Remove row name
% Scaling factor
Edit_scalingfactor = uicontrol('Parent', tab_microstructures,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','1',...
    'visible','off','enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.7 0.3 0.1375 0.04]);
Button_scalingfactor = uicontrol('Parent', tab_microstructures, 'Style', 'pushbutton', 'String', 'Apply re-size',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.8375 0.3 0.1375 0.04],...
    'visible','off','enable','off','Callback',{@Scalingfactor_Callback});
Text_scaling_factor= uicontrol('Parent', tab_microstructures, 'Visible','off','Style','text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',{'Choose a real number >0 such as:','>1: upscaling','<1: downscaling'},...
     'BackgroundColor','w','ForegroundColor','k','HorizontalAlignment','left','Units','normalized','Position', [0.7 0.19 0.275 0.1]);
Text_heterogeneous_instructions= uicontrol('Parent', tab_microstructures, 'Visible','off','Style','text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',{'Any modifcations done in 1) reset 2) and 3). Any modifications done in 2) reset 3).','Once done, you can visualize the microstructure by clicking on the button below.','If you are satisfied with your microstructure, click on save electrode.','Then click on New/reset electrode to setup another microstructure.','Once all the geometries have been setup and saved, go to generate separator (if you have not setup the separator) then go to cell options.'}',...
     'BackgroundColor','w','ForegroundColor','k','HorizontalAlignment','left','Units','normalized','Position', [0.35 0.09 0.625 0.12]);
Button_visualize_microstructure = uicontrol('Parent', tab_microstructures, 'Style', 'pushbutton', 'String', 'Click to visualize microstructure',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.35 0.05 0.625 0.04],...
    'visible','off','enable','off','Callback',{@Visualize_microstructure_Callback});
 
%% CALLBACKS
    function newelectrode_Callback(~,~)
        electrode_numero = Popupmenu_chooseelectrode.Value; % Electrode selected
        [celldata] = reset_electrode (celldata,electrode_numero); % Erase electrode information
        Text_volume.String='No volume selected';
        turn_GUI_choicegeometry('on')
        turn_GUI_parametergeometry_heterogenous('off')
        turn_GUI_parametergeometry_homogenous('off')
        set(pushebutton_microstructures_save,'enable','off');
        if ismember(findall(0,'type','figure'),Fig_visualizemicrostructure)
            close Figure_visualizemicrostructure;
        end
    end

    function choice_representation_Callback(~,~)
        choice_representation = Popupmenu_choice_representation.Value;
        electrode_numero = Popupmenu_chooseelectrode.Value; % Electrode selected
        if choice_representation==2
            str_geometry = 'Imported heterogenous medium';
            turn_GUI_parametergeometry_homogenous('off')
            turn_GUI_parametergeometry_heterogenous('on')
            celldata(electrode_numero).heterogeneous = true;
            Text_volume.String='No volume selected';
        elseif choice_representation==3
            str_geometry = 'Homogeneous medium';
            turn_GUI_parametergeometry_heterogenous('off')
            turn_GUI_parametergeometry_homogenous('on')
            celldata(electrode_numero).heterogeneous = false;
            Text_volume.String='No volume selected';
        else
            turn_GUI_parametergeometry_heterogenous('off')
            turn_GUI_parametergeometry_homogenous('off')            
        end
        set(Text_geometry,'String',str_geometry);
        set(pushebutton_microstructures_save,'enable','on');
        if ismember(findall(0,'type','figure'),Fig_visualizemicrostructure)
            close Figure_visualizemicrostructure;
        end
    end

    function edit_length_homogenous_Callback(~,~)
        choice_representation = Popupmenu_choice_representation.Value;
        electrode_numero = Popupmenu_chooseelectrode.Value; % Electrode selected
        celldata(electrode_numero).length_voxel = str2double(Edit_length_homogenous.String);
    end

    function import_microstructure_Callback(~,~)
        % Open dialog box to choose file path
        [FileName,PathName,~] = uigetfile({'*.tif','Tif image (*.tif)';'*.m;*.mat','MATLAB File (*.m,*.mat)'},'Select volume');
        if FileName==0
            % User clicked cancel button or closed the dialog box
        else
            full_path_new_volume = [PathName FileName]; % Full path
            [array_initial, outcome] = function_loadvolume(full_path_new_volume, 'uint8', 'none' ); % Load new volume
            if outcome.success % Success to import
                % Update GUI
                Text_volume.String=FileName;
                set(Text_ROI_error,'Visible','on','String','Region of interest (ROI) is correct');
                set(Text_ROI_instructions,'Visible','on');
                % Update table
                Domain_size=size(array_initial); % Get domain size      
                ROI_data = [num2cell([1;2;3]) num2cell(Domain_size') num2cell([1;1;1]) num2cell(Domain_size') num2cell(Domain_size')];
                table_ROI.Data=ROI_data;
                set(table_ROI,'Visible','on','enable','on');
                
                array_initial(array_initial~=0)=1;
                
                % Phase
                Phase_code=unique(array_initial); % Get phase code
                number_phase=length(Phase_code); % Get number of phase
                % volume fraction
                total_number_voxel=prod(Domain_size);
                volumefraction=zeros(number_phase,1);
                for n=1:1:number_phase
                    volumefraction(n,1)=sum(sum(sum(array_initial==Phase_code(n))))/total_number_voxel;
                end
                Volumefraction_data = [num2cell(Phase_code) num2cell(volumefraction) num2cell(volumefraction) num2cell(volumefraction)];
                table_volumefraction.Data=Volumefraction_data;
                set(Text_volumefractions,'Visible','on');
                set(table_volumefraction,'Visible','on','enable','on');
     
                % Orientation table
                Orientation_data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell([1;2;3]) num2cell(Domain_size')];
                table_orientations.Data=Orientation_data;
                set(Text_orientation_instructions,'Visible','on');
                set(table_orientations,'Visible','on','enable','on');                
                set(togglebutton_shiftdirection_1_with_2,'Visible','on','enable','on');
                set(togglebutton_shiftdirection_1_with_3,'Visible','on','enable','on');
                set(togglebutton_shiftdirection_2_with_3,'Visible','on','enable','on');
                
                % Reverse direction
                set(togglebutton_inversedirection_1,'Visible','on','enable','on');   
                set(togglebutton_inversedirection_2,'Visible','on','enable','on');   
                set(togglebutton_inversedirection_3,'Visible','on','enable','on');   
                                
                % Resize
                Voxelresize_data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(Domain_size') num2cell(Domain_size')];
                table_voxelresize.Data=Voxelresize_data;                
                set(Edit_scalingfactor,'Visible','on','enable','on');
                set(Text_voxelsize_instructions,'Visible','on');
                set(table_voxelresize,'Visible','on','enable','on');
                set(Edit_scalingfactor,'Visible','on','enable','on','String','1');
                set(Button_scalingfactor,'Visible','on','enable','on');
                set(Text_scaling_factor,'Visible','on');
                
                % Visualize
                set(Button_visualize_microstructure,'Visible','on','enable','on');
                set(Text_heterogeneous_instructions,'Visible','on');

                % Initialize other array
                array_ROI = array_initial;
                array_orientation = array_initial;
                array_resize = array_initial;
            end
        end
    end

    function cell_table_ROI_Callback(~,~)
        % Get back data of the table
        table_data = table_ROI.Data;
        ROI_=table_data(:,3:4);
        ROI_ = cell2mat(ROI_);
        max_roiend=table_data(:,2);
        max_roiend=cell2mat(max_roiend);
        % Check validity of user-defined ROI
        [correct_ROI] = check_validity_ROI (ROI_, max_roiend);
        if correct_ROI==1
            set(Text_ROI_error,'ForegroundColor','k','String','Region of interest (ROI) is correct'); % Update error message
            Voxel_ROI = ROI_(:,2) - ROI_(:,1) + 1; % Update ROI value
            table_data(:,5) = num2cell(Voxel_ROI);
            table_ROI.Data=table_data; % Update ui table
            array_ROI = array_initial( ROI_(1,1):ROI_(1,2) , ROI_(2,1):ROI_(2,2) , ROI_(3,1):ROI_(3,2)); % Crop domain
            array_orientation = array_ROI; % Rest other array
            array_resize = array_ROI;
            Domain_size=size(array_ROI); % Domain size      
            % Reset orientation
            Orientation_data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell([1;2;3]) num2cell(Domain_size')];
            table_orientations.Data=Orientation_data;
            togglebutton_shiftdirection_1_with_2.Value=0;
            togglebutton_shiftdirection_1_with_3.Value=0;
            togglebutton_shiftdirection_2_with_3.Value=0;
            togglebutton_shiftdirection_1_with_2.String='1 - 2 (Un-swaped)';
            togglebutton_shiftdirection_1_with_2.BackgroundColor='w';
            togglebutton_shiftdirection_1_with_3.String='1 - 3 (Un-swaped)';
            togglebutton_shiftdirection_1_with_3.BackgroundColor='w';
            togglebutton_shiftdirection_2_with_3.String='2 - 3 (Un-swaped)';
            togglebutton_shiftdirection_2_with_3.BackgroundColor='w';
            togglebutton_inversedirection_1.Value=0;
            togglebutton_inversedirection_2.Value=0;
            togglebutton_inversedirection_3.Value=0;
            togglebutton_inversedirection_1.String='1 (Un-changed)';
            togglebutton_inversedirection_1.BackgroundColor='w';
            togglebutton_inversedirection_2.String='1 (Un-changed)';
            togglebutton_inversedirection_2.BackgroundColor='w';
            togglebutton_inversedirection_3.String='1 (Un-changed)';
            togglebutton_inversedirection_3.BackgroundColor='w';
            % Rest voxel resize
            Voxelresize_data = [{'Through-plane direction';'In-plane direction 1';'In-plane direction 2'} num2cell(Domain_size') num2cell(Domain_size')];
            table_voxelresize.Data=Voxelresize_data;
            set(Edit_scalingfactor,'String','1');
                        
            % Update volume fractions
            total_number_voxel=prod(Domain_size);
            volumefraction=zeros(number_phase,1);
            for n=1:1:number_phase
                volumefraction(n,1)=sum(sum(sum(array_ROI==Phase_code(n))))/total_number_voxel;
            end
            Volumefraction_data = table_volumefraction.Data;
            Volumefraction_data(:,3) = num2cell(volumefraction);
            Volumefraction_data(:,4) = Volumefraction_data(:,3);
            table_volumefraction.Data=Volumefraction_data;
        else
            set(Text_ROI_error,'ForegroundColor','r','String','Region of interest (ROI) is incorrect'); % Update error message
        end
    end


    function Shiftdirection_1_with_2_Callback(~,~)
        button_pushed=togglebutton_shiftdirection_1_with_2.Value;
        if button_pushed==0
            togglebutton_shiftdirection_1_with_2.String='1 - 2 (Un-swaped)';
            togglebutton_shiftdirection_1_with_2.BackgroundColor='w';
        else
             togglebutton_shiftdirection_1_with_2.String='1 - 2 (Swaped)';
             togglebutton_shiftdirection_1_with_2.BackgroundColor=[0.75 0.75 0.75];
        end
        % Switch axe 1 with axe 2
        Domain_size = size(array_orientation);
        New_Domain_size=zeros(1,3);
        New_Domain_size(1)=Domain_size(2); New_Domain_size(2)=Domain_size(1); New_Domain_size(3)=Domain_size(3);
        New_array=zeros(New_Domain_size);
        for k=1:1:Domain_size(2)
            slice=array_orientation(:,k,:);
            New_array(k,:,:)=slice;
        end
        array_orientation=New_array;
        Domain_size=size(array_orientation);
        % Update table
        Orientation_data = table_orientations.Data;
        axis_orientation = cell2mat(Orientation_data(:,2)); axis_orientation_tmp = axis_orientation;
        axis_orientation(2) = axis_orientation_tmp(1);
        axis_orientation(1) = axis_orientation_tmp(2);
        Orientation_data(:,2) = num2cell(axis_orientation);
        Orientation_data(:,3) = num2cell(Domain_size');
        table_orientations.Data=Orientation_data;   
        % Update resize
        array_resize = array_orientation;
        Voxelresize_data = table_voxelresize.Data;
        Voxelresize_data(:,2) = Orientation_data(:,3);
        Voxelresize_data(:,3) = Voxelresize_data(:,2);
        table_voxelresize.Data=Voxelresize_data;
        set(Edit_scalingfactor,'String','1');
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = Volumefraction_data(:,3);
        table_volumefraction.Data=Volumefraction_data;
        
    end 
    function Shiftdirection_1_with_3_Callback(~,~)
        button_pushed=togglebutton_shiftdirection_1_with_3.Value;
        if button_pushed==0
            togglebutton_shiftdirection_1_with_3.String='1 - 3 (Un-swaped)';
            togglebutton_shiftdirection_1_with_3.BackgroundColor='w';
        else
            togglebutton_shiftdirection_1_with_3.String='1 - 3 (Swaped)';
            togglebutton_shiftdirection_1_with_3.BackgroundColor=[0.75 0.75 0.75];
        end
        % Switch axe 3 with axe 1
        Domain_size = size(array_orientation);
        New_Domain_size=zeros(1,3);
        New_Domain_size(1)=Domain_size(3); New_Domain_size(2)=Domain_size(2); New_Domain_size(3)=Domain_size(1);
        New_array=zeros(New_Domain_size);
        for k=1:1:Domain_size(3)
            slice=array_orientation(:,:,k)';
            New_array(k,:,:)=slice;
        end
        array_orientation=New_array;
        Domain_size=size(array_orientation);
        % Update table
        Orientation_data = table_orientations.Data;
        axis_orientation = cell2mat(Orientation_data(:,2)); axis_orientation_tmp = axis_orientation;
        axis_orientation(3) = axis_orientation_tmp(1);
        axis_orientation(1) = axis_orientation_tmp(3);
        Orientation_data(:,2) = num2cell(axis_orientation);        
        Orientation_data(:,3) = num2cell(Domain_size');
        table_orientations.Data=Orientation_data;   
        % Update resize
        array_resize = array_orientation;
        Voxelresize_data = table_voxelresize.Data;
        Voxelresize_data(:,2) = Orientation_data(:,3);
        Voxelresize_data(:,3) = Voxelresize_data(:,2);
        table_voxelresize.Data=Voxelresize_data;
        set(Edit_scalingfactor,'String','1');        
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = Volumefraction_data(:,3);
        table_volumefraction.Data=Volumefraction_data;        
    end
    function Shiftdirection_2_with_3_Callback(~,~)
        button_pushed=togglebutton_shiftdirection_2_with_3.Value;
        if button_pushed==0
            togglebutton_shiftdirection_2_with_3.String='2 - 3 (Un-swaped)';
            togglebutton_shiftdirection_2_with_3.BackgroundColor='w';
        else
            togglebutton_shiftdirection_2_with_3.String='2 - 3 (Swaped)';
            togglebutton_shiftdirection_2_with_3.BackgroundColor=[0.75 0.75 0.75];
        end
        % Switch axe 3 with axe 2
        Domain_size = size(array_orientation);
        New_Domain_size=zeros(1,3);
        New_Domain_size(1)=Domain_size(1); New_Domain_size(2)=Domain_size(3); New_Domain_size(3)=Domain_size(2);
        New_array=zeros(New_Domain_size);
        for k=1:1:Domain_size(3)
            slice=array_orientation(:,:,k);
            New_array(:,k,:)=slice;
        end
        array_orientation=New_array;
        Domain_size=size(array_orientation);
        % Update table
        Orientation_data = table_orientations.Data;
        Orientation_data(:,3) = num2cell(Domain_size');
        axis_orientation = cell2mat(Orientation_data(:,2)); axis_orientation_tmp = axis_orientation;
        axis_orientation(3) = axis_orientation_tmp(2);
        axis_orientation(2) = axis_orientation_tmp(3);
        Orientation_data(:,2) = num2cell(axis_orientation);        
        Orientation_data(:,3) = num2cell(Domain_size');        
        table_orientations.Data=Orientation_data;   
        % Update resize
        array_resize = array_orientation;
        Voxelresize_data = table_voxelresize.Data;
        Voxelresize_data(:,2) = Orientation_data(:,3);
        Voxelresize_data(:,3) = Voxelresize_data(:,2);
        table_voxelresize.Data=Voxelresize_data;
        set(Edit_scalingfactor,'String','1');      
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = Volumefraction_data(:,3);
        table_volumefraction.Data=Volumefraction_data;        
    end

   function Inversedirection_1_Callback(~,~)
        button_pushed=togglebutton_inversedirection_1.Value;
        if button_pushed==0
            togglebutton_inversedirection_1.String='1 (Un-changed)';
            togglebutton_inversedirection_1.BackgroundColor='w';
        else
             togglebutton_inversedirection_1.String='1 (reverse)';
             togglebutton_inversedirection_1.BackgroundColor=[0.75 0.75 0.75];
        end
        Domain_size = size(array_orientation);
        tmp=zeros(Domain_size);
        for k=1:1:Domain_size(1)
            slice=array_orientation(k,:,:);
            tmp(Domain_size(1)-k+1,:,:)=slice;
        end
        array_orientation=tmp;
        % Update resize
        array_resize = array_orientation;
        Voxelresize_data = table_voxelresize.Data;
        Orientation_data = table_orientations.Data;
        Orientation_data(:,3) = num2cell(Domain_size');        
        Voxelresize_data(:,2) = Orientation_data(:,3);
        Voxelresize_data(:,3) = Voxelresize_data(:,2);
        table_voxelresize.Data=Voxelresize_data;
        set(Edit_scalingfactor,'String','1');     
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = Volumefraction_data(:,3);
        table_volumefraction.Data=Volumefraction_data;        
    end 
    function Inversedirection_2_Callback(~,~)
        button_pushed=togglebutton_inversedirection_2.Value;
        if button_pushed==0
            togglebutton_inversedirection_2.String='2 (Un-changed)';
            togglebutton_inversedirection_2.BackgroundColor='w';
        else
            togglebutton_inversedirection_2.String='2 (reverse)';
            togglebutton_inversedirection_2.BackgroundColor=[0.75 0.75 0.75];
        end
        Domain_size = size(array_orientation);
        tmp=zeros(Domain_size);
        for k=1:1:Domain_size(2)
            slice=array_orientation(:,k,:);
            tmp(:,Domain_size(2)-k+1,:)=slice;
        end
        array_orientation=tmp;
        % Update resize
        array_resize = array_orientation;
        Voxelresize_data = table_voxelresize.Data;
        Orientation_data = table_orientations.Data;
        Orientation_data(:,3) = num2cell(Domain_size');          
        Voxelresize_data(:,2) = Orientation_data(:,3);
        Voxelresize_data(:,3) = Voxelresize_data(:,2);
        table_voxelresize.Data=Voxelresize_data;
        set(Edit_scalingfactor,'String','1');      
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = Volumefraction_data(:,3);
        table_volumefraction.Data=Volumefraction_data;        
    end
    function Inversedirection_3_Callback(~,~)
        button_pushed=togglebutton_inversedirection_3.Value;
        if button_pushed==0
            togglebutton_inversedirection_3.String='3 (Un-changed)';
            togglebutton_inversedirection_3.BackgroundColor='w';
        else
            togglebutton_inversedirection_3.String='3 (reverse)';
            togglebutton_inversedirection_3.BackgroundColor=[0.75 0.75 0.75];
        end
        Domain_size = size(array_orientation);
        tmp=zeros(Domain_size);
        for k=1:1:Domain_size(3)
            slice=array_orientation(:,:,k);
            tmp(:,:,Domain_size(3)-k+1)=slice;
        end
        array_orientation=tmp;
        % Update resize
        array_resize = array_orientation;
        Voxelresize_data = table_voxelresize.Data;
        Orientation_data = table_orientations.Data;
        Orientation_data(:,3) = num2cell(Domain_size');          
        Voxelresize_data(:,2) = Orientation_data(:,3);
        Voxelresize_data(:,3) = Voxelresize_data(:,2);
        table_voxelresize.Data=Voxelresize_data;
        set(Edit_scalingfactor,'String','1');  
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = Volumefraction_data(:,3);
        table_volumefraction.Data=Volumefraction_data;        
    end

    function saveelectrode_Callback(~,~)
        electrode_numero = Popupmenu_chooseelectrode.Value; % Electrode selected
        celldata(electrode_numero).exist = true;
        choice_representation = Popupmenu_choice_representation.Value;
        if choice_representation==2
            celldata(electrode_numero).array = array_resize;
            celldata(electrode_numero).heterogeneous = true;
        else
            celldata(electrode_numero).heterogeneous = false;
            celldata(electrode_numero).length_voxel = str2double(Edit_length_homogenous.String);         
        end
        turn_GUI_choicegeometry('off')
        turn_GUI_parametergeometry_heterogenous('off')
        turn_GUI_parametergeometry_homogenous('off')
        set(pushebutton_microstructures_save,'enable','off');
        Text_volume.String='No volume selected';
        if ismember(findall(0,'type','figure'),Fig_visualizemicrostructure)
            close Figure_visualizemicrostructure;
        end
        if celldata(1).exist && celldata(2).exist && celldata(3).exist % Update GUI of cell tab
            set(Popup_cell_id,'value',1,'String', {'Cell id '});
            set(Popup_cell_assign,'value',1,'String', {'Cell assignment','Left electrode: solid','Left electrode: electrolyte','Separator: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'});
            Checkdimension_cell
            phaseid_cell
        end
        if celldata(1).exist && celldata(3).exist % Separator can be generated if both electrodes have been setup
            set(pushbutton_separator_new,'enable','on');
        else
            set(pushbutton_separator_new,'enable','off');
        end
    end

    function Scalingfactor_Callback(~,~)
        scaling_factor = str2double(Edit_scalingfactor.String);
        array_resize=imresize3(array_orientation,scaling_factor,'nearest');
        Domain_size = size(array_resize);
        % Update volume fractions
        total_number_voxel=prod(Domain_size);
        volumefraction=zeros(number_phase,1);
        for n=1:1:number_phase
            volumefraction(n,1)=sum(sum(sum(array_resize==Phase_code(n))))/total_number_voxel;
        end
        Volumefraction_data = table_volumefraction.Data;
        Volumefraction_data(:,4) = num2cell(volumefraction);
        table_volumefraction.Data=Volumefraction_data;       
        % Update resize
        Voxelresize_data = table_voxelresize.Data;
        Voxelresize_data(:,3) = num2cell(Domain_size');
        table_voxelresize.Data=Voxelresize_data;
    end

    function Visualize_microstructure_Callback(~,~)
        if ismember(findall(0,'type','figure'),Fig_visualizemicrostructure)
            close Figure_visualizemicrostructure;
        end
        Fig_visualizemicrostructure = figure; % Create figure
        Fig_visualizemicrostructure.Name= 'Visualize microstructure';
        Fig_visualizemicrostructure.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig_visualizemicrostructure,'position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]); % Full screen figure
        % Figure for volume view
        axes_2dview = axes('Parent', Fig_visualizemicrostructure,'Visible','off','FontName',font_name_GUI,'Units','normalized','Position', [0 0.1 1 0.85]);
        % Remove tick and label
        set(axes_2dview,'xtick',[],'ytick',[]);
        set(axes_2dview,'xticklabel',[],'yticklabel',[]);
        % Box on
        box(axes_2dview,'on');
        % Fit the axes box
        axis(axes_2dview,'tight');
        % Aspect ratio is 1:1
        axis(axes_2dview,'equal');
        
        % Slider
        Slider_axes_2dview = uicontrol('Parent', Fig_visualizemicrostructure,'Style', 'slider','Min',1,'Max',100,'Value',1,'Units','normalized','Position', [0.3 0.025 0.4 0.04],'Callback', @slider_axes2dview_Callback,...
            'Visible','on','enable','on');
        % Text (slider position)
        Text_slider = uicontrol('Parent', Fig_visualizemicrostructure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Position: -/-','Visible','on','enable','on',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.725 0.022 0.15 0.04]);
        % Popup menu (slider direction)
        Popup_slider = uicontrol('Parent', Fig_visualizemicrostructure,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
            'String', {'Through-plane direction','In-plane direction 1','In-plane direction 2'},'Value',1,'Units','normalized','Position', [0.05 0.025 0.215 0.04],'enable','on','Visible','on','Callback', @popup_axes2dview_Callback);
        
        direction = 1; pos_=1;   
        set(Text_slider,'String',['Position: ' num2str(pos_) '/' num2str(Domain_size(direction))]);
        minor_step = 1/(Domain_size(direction)-1);
        major_step = 0.1;
        set(Slider_axes_2dview,'Min',1,'Max',Domain_size(direction),'SliderStep', [minor_step, major_step],'Value',1);
        
        % Slider
        function slider_axes2dview_Callback(source,~)
            % Get position value
            pos_=round(source.Value);
            direction = Popup_slider.Value;
            % Update text
            set(Text_slider,'String',['Position: ' num2str(pos_) '/' num2str(Domain_size(direction))]);
            % Update position array
            % Position_slice(direction)=pos_;
            % Update figure
            update_figure
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
            set(Text_slider,'String',['Position: ' num2str(1) '/' num2str(Domain_size(direction))]);
            % Update figure
            update_figure
        end
        
        function update_figure
            % Colors
            for current_phase=1:1:number_phase
                RGB_phase.index(current_phase).rgb = [color_phase(current_phase,1) color_phase(current_phase,2) color_phase(current_phase,3)];
            end
            direction = Popup_slider.Value;
            pos_ = round(Slider_axes_2dview.Value);
            Position_slice(direction)=pos_;
            if direction==1
                % Initializaion
                slice_color = zeros(Domain_size(2),Domain_size(3),3); % RGB color map
                slice_r = zeros(Domain_size(2),Domain_size(3)); % Red color map
                slice_g = zeros(Domain_size(2),Domain_size(3)); % Green color map
                slice_b = zeros(Domain_size(2),Domain_size(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:number_phase
                    code_tmp =Phase_code(current_phase); % Current phase code
                    slice_r(array_resize(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(array_resize(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(array_resize(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
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
                    slice_r(array_resize(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(array_resize(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(array_resize(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
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
                    slice_r(array_resize(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(array_resize(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(array_resize(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            end
            slice_color(:,:,1)=slice_r; % Attribute RGB color
            slice_color(:,:,2)=slice_g;
            slice_color(:,:,3)=slice_b;
            % Display the slice
            image(slice_color,'parent',axes_2dview);
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
        end
               
        % Initialize
        update_figure
        
    end

%% FUNCTIONS
    function turn_GUI_choicegeometry(value)
        set(Text_chooserepresentation,'Visible',value);
        set(Popupmenu_choice_representation,'Visible',value,'enable',value);
    end

    function turn_GUI_parametergeometry_heterogenous(value)
        set(Text_geometry,'Visible',value);
        set(Text_volume,'Visible',value);
        set(Button_import_microstructure,'Visible',value,'enable',value);       
        %value
        if  strcmp(value,'off')
            set(table_ROI,'Visible',value,'enable',value);
            set(Text_ROI_error,'Visible',value);
            set(Text_ROI_instructions,'Visible',value);
            set(Text_volumefractions,'Visible',value);
            set(table_volumefraction,'Visible',value,'enable',value);
            set(table_orientations,'Visible',value,'enable',value);
            set(Text_orientation_instructions,'Visible',value);
            set(togglebutton_shiftdirection_1_with_2,'Visible',value,'enable',value);
            set(togglebutton_shiftdirection_1_with_3,'Visible',value,'enable',value);
            set(togglebutton_shiftdirection_2_with_3,'Visible',value,'enable',value);
            set(togglebutton_inversedirection_1,'Visible',value,'enable',value);
            set(togglebutton_inversedirection_2,'Visible',value,'enable',value);
            set(togglebutton_inversedirection_3,'Visible',value,'enable',value);
            set(Text_voxelsize_instructions,'Visible',value);
            set(table_voxelresize,'Visible',value,'enable',value);
            set(Edit_scalingfactor,'Visible',value,'enable',value);
            set(Button_scalingfactor,'Visible',value,'enable',value);
            set(Text_scaling_factor,'Visible',value);
            set(Button_visualize_microstructure,'Visible',value,'enable',value);
            set(Text_heterogeneous_instructions,'Visible',value);
        end
    end
    function turn_GUI_parametergeometry_homogenous(value)
        set(Text_geometry,'Visible',value);
        set(Text_length_homogeneous,'Visible',value);
        set(Edit_length_homogenous,'Visible',value,'enable',value);

    end

    function [celldata] = reset_electrode (celldata,electrode)
        celldata(electrode).name = '';
        celldata(electrode).array = -1;
        celldata(electrode).exist = false;
        celldata(electrode).heterogeneous = -1;
        celldata(electrode).length_voxel = -1;
    end

    function [correct_ROI] = check_validity_ROI (ROI_, max_ROI)
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
            tmp=max_ROI-ROI_(:,2);
            tmp2=max_ROI*0;
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
    end

%%
%% "GENERATE SEPARATOR" TAB
%%

%% PRE-ASSIGNED PARAMETERS
Fig_visualizeseparator=0;
choice_separator = {' ', 'Connected plates','Strips'};
% Default connected plates parameters
separator_length = 50;
parameters_separator_connectedplates.distancefromedge = 0;
parameters_separator_connectedplates.centered = 1;
parameters_separator_connectedplates.inplanelength = 9;
parameters_separator_connectedplates.inplanedistance = 3;
parameters_separator_connectedplates.thicknesslayers135 = 2;
parameters_separator_connectedplates.thicknesslayers246 = 2;
parameters_separator_connectedplates.thicknessnolayer = 1;
parameters_separator_connectedplates.only_alllayers = 1;
parameters_separator_connectedplates.centered_alllayers = 1;
% Default strips parameters
parameters_separator_strips.alignement = 2;
parameters_separator_strips.distancefromedge = 0;
parameters_separator_strips.centered = 1;
parameters_separator_strips.stripwidth = 15;
parameters_separator_strips.striplength = 0.9;
parameters_separator_strips.inplanedistance = round(parameters_separator_strips.stripwidth/3);
parameters_separator_strips.thicknesslayers1357 = 3;
parameters_separator_strips.thicknesslayers2468 = 2;
parameters_separator_strips.thicknessnolayer = 3;
parameters_separator_strips.only_alllayers = 0;
parameters_separator_strips.centered_alllayers = 1;

%% GUI
% Title
Text_tab_separator = uicontrol('Parent', tab_separator, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Generate separator');

% New/Save separator
pushbutton_separator_new = uicontrol('Parent', tab_separator,'Style', 'pushbutton','Value',0,'String','New/reset separator','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.87 0.2 0.04],'Callback', @newseparator_Callback);
pushbutton_separator_save = uicontrol('Parent', tab_separator,'Style', 'pushbutton','Value',0,'String','Save separator','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.23 0.87 0.2 0.04],'Callback', @saveseparator_Callback);

% Choose separator type
Popupmenu_chooseseparator = uicontrol('Parent', tab_separator, 'Style', 'popupmenu', 'String', choice_separator,...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.82 0.4050 0.04]);
Text_ = uicontrol('Parent', tab_separator, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Choose the separator geometry to mesh. You must set up the two electrodes before, in the tab "Microstructures". Once done, save it.',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.745 0.405 0.08]);

% Enforced parameters
Text_separator_dimension = uicontrol('Parent', tab_separator, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',' ',...
    'visible','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.5 0.87 0.49 0.04]);
Text_separator_direction = uicontrol('Parent', tab_separator, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Imposed parameter: Direction = 1',...
    'visible','off','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.5 0.84 0.49 0.04]);

% Separator parameters 
table_separator_parameters = uitable('Parent', tab_separator,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.025 0.4 0.95 0.325]);
table_separator_parameters.ColumnName = {'Parameter                                                                                                          ','Value'}; % Column name
table_separator_parameters.ColumnEditable = [false true]; % Select column editable
table_separator_parameters.ColumnWidth = {'auto', 'auto'}; % Auto width
table_separator_parameters.RowName = []; % Remove row name
% Connected plates 
data_separator_connectedplates = [{'Length in voxels';'In-plane distance from edge w/o plates in voxels';'Center plate in-plane pattern. 1 (true) or 0 (false)';'In-plane plate length in voxels';'In-plane distance between plate in voxels';'Plate thickness in voxels';'Distance between plates in voxels';'Distance w/o plates from the border';'Only mesh a complete pattern. 1 (true) or 0 (false)';'Center plate pattern along thickness. 1 (true) or 0 (false)'},...
         {separator_length; parameters_separator_connectedplates.distancefromedge;parameters_separator_connectedplates.centered;parameters_separator_connectedplates.inplanelength;parameters_separator_connectedplates.inplanedistance;parameters_separator_connectedplates.thicknesslayers135;parameters_separator_connectedplates.thicknesslayers246;parameters_separator_connectedplates.thicknessnolayer;parameters_separator_connectedplates.only_alllayers;parameters_separator_connectedplates.centered_alllayers}];
% Strips     
data_separator_strips = [{'Length in voxels';'In-plane alignment (2 or 3)';'In-plane distance from edge w/o strips in voxels';'Center plate in-plane pattern. 1 (true) or 0 (false)';'In-plane strip width in voxels'; 'In-plane strip length, expressed in ratio of the domain''s length (>0, <1)';'In-plane distance between strips in voxels';'Strips thickness in voxels';'Distance between strips in voxels';'Distance w/o strips from the border';'Only mesh a complete pattern. 1 (true) or 0 (false)';'Center strips pattern along thickness. 1 (true) or 0 (false)'},...
         {separator_length; parameters_separator_strips.alignement; parameters_separator_strips.distancefromedge; parameters_separator_strips.centered; parameters_separator_strips.stripwidth; parameters_separator_strips.striplength; parameters_separator_strips.inplanedistance; parameters_separator_strips.thicknesslayers1357; parameters_separator_strips.thicknesslayers2468; parameters_separator_strips.thicknessnolayer; parameters_separator_strips.only_alllayers; parameters_separator_strips.centered_alllayers}];
table_separator_parameters.Data = [];

% Generate separator
Button_generate_separator = uicontrol('Parent', tab_separator, 'Style', 'pushbutton', 'String', 'Click to generate separator',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.05 0.25 0.04],...
    'visible','on','enable','off','Callback',{@Generate_separator_Callback});

% Visualize separator
Button_visualize_separator = uicontrol('Parent', tab_separator, 'Style', 'pushbutton', 'String', 'Click to visualize separator',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.375 0.05 0.25 0.04],...
    'visible','on','enable','off','Callback',{@Visualize_separator_Callback});

% Save separator
Button_save_separator = uicontrol('Parent', tab_separator, 'Style', 'pushbutton', 'String', 'Click to save separator',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.725 0.05 0.25 0.04],...
    'visible','on','enable','off','Callback',{@Save_separator_Callback});

%% CALLBACKS
    function newseparator_Callback(~,~)
        Checkdimension_forseparator
        set(Text_separator_dimension,'visible','on','String',['Imposed parameter: In-plane dimension = ' num2str(min_inplane(1),'%1i') ' and ' num2str(min_inplane(2),'%1i')]);
        set(Text_separator_direction,'visible','on');
        if Popupmenu_chooseseparator.Value==1
            set(table_separator_parameters,'Visible','off','enable','off');
            set(Button_generate_separator,'enable','off')
        else
            if Popupmenu_chooseseparator.Value==2
                table_separator_parameters.Data = data_separator_connectedplates;
            elseif Popupmenu_chooseseparator.Value==3
                table_separator_parameters.Data = data_separator_strips;
            end
            set(table_separator_parameters,'Visible','on','enable','on');
            set(Button_generate_separator,'enable','on')
        end
    end

    function Generate_separator_Callback(~,~)
        if  Popupmenu_chooseseparator.Value==2
            % Set dimension and direction parameter
            tmp = cell2mat(table_separator_parameters.Data(:,2));
            parameters_separator_connectedplates.dimension =[tmp(1) min_inplane(1) min_inplane(2)];
            parameters_separator_connectedplates.direction = 1;
            [separator_array] = Function_generate_connectedplates(parameters_separator_connectedplates);
        elseif Popupmenu_chooseseparator.Value==3
            % Set dimension and direction parameter
            tmp = cell2mat(table_separator_parameters.Data(:,2));
            parameters_separator_strips.dimension =[tmp(1) min_inplane(1) min_inplane(2)];
            parameters_separator_strips.direction = 1;
            [separator_array] = Function_generate_strips(parameters_separator_strips);           
        end
        set(Button_visualize_separator,'enable','on')
        set(Button_save_separator,'enable','on')
    end

    function Save_separator_Callback(~,~)
        celldata(2).heterogeneous = true;
        celldata(2).exist = true;
        celldata(2).array = separator_array;
        if ismember(findall(0,'type','figure'),Fig_visualizeseparator)
            close Fig_visualizeseparator;
        end
        if celldata(1).exist && celldata(2).exist && celldata(3).exist % Update GUI of cell tab
            set(Popup_cell_id,'value',1,'String', {'Cell id '});
            set(Popup_cell_assign,'value',1,'String', {'Cell assignment','Left electrode: solid','Left electrode: electrolyte','Separator: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'});
            Checkdimension_cell
            phaseid_cell
        end
    end

    function Visualize_separator_Callback(~,~)
        if ismember(findall(0,'type','figure'),Fig_visualizeseparator)
            close Fig_visualizeseparator;
        end
        Fig_visualizeseparator = figure; % Create figure
        Fig_visualizeseparator.Name= 'Visualize microstructure';
        Fig_visualizeseparator.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig_visualizeseparator,'position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]); % Full screen figure
        % Figure for volume view
        axes_2dview = axes('Parent', Fig_visualizeseparator,'Visible','off','FontName',font_name_GUI,'Units','normalized','Position', [0 0.1 1 0.85]);
        % Remove tick and label
        set(axes_2dview,'xtick',[],'ytick',[]);
        set(axes_2dview,'xticklabel',[],'yticklabel',[]);
        % Box on
        box(axes_2dview,'on');
        % Fit the axes box
        axis(axes_2dview,'tight');
        % Aspect ratio is 1:1
        axis(axes_2dview,'equal');
        
        % Slider
        Slider_axes_separator_2dview = uicontrol('Parent', Fig_visualizeseparator,'Style', 'slider','Min',1,'Max',100,'Value',1,'Units','normalized','Position', [0.3 0.025 0.4 0.04],'Callback', @slider_axes2dview_separator_Callback,...
            'Visible','on','enable','on');
        % Text (slider position)
        Text_separator_slider = uicontrol('Parent', Fig_visualizeseparator,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Position: -/-','Visible','on','enable','on',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.725 0.022 0.15 0.04]);
        % Popup menu (slider direction)
        Popup_separator_slider = uicontrol('Parent', Fig_visualizeseparator,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
            'String', {'Through-plane direction','In-plane direction 1','In-plane direction 2'},'Value',1,'Units','normalized','Position', [0.05 0.025 0.215 0.04],'enable','on','Visible','on','Callback', @popup_axes2dview_separator_Callback);
        
        dim_ = size(separator_array);
        direction = 1; pos_=1;
        set(Text_separator_slider,'String',['Position: ' num2str(pos_) '/' num2str(dim_(direction))]);
        minor_step = 1/(dim_(direction)-1);
        major_step = 0.1;
        set(Slider_axes_separator_2dview,'Min',1,'Max',dim_(direction),'SliderStep', [minor_step, major_step],'Value',1);
        
        % Slider
        function slider_axes2dview_separator_Callback(source,~)
            % Get position value
            pos_=round(source.Value);
            direction = Popup_separator_slider.Value;
            % Update text
            set(Text_separator_slider,'String',['Position: ' num2str(pos_) '/' num2str(dim_(direction))]);
            % Update position array
            % Position_slice(direction)=pos_;
            % Update figure
            update_separtor_figure
        end
        
        % Select direction
        function popup_axes2dview_separator_Callback(source,~)
            dim_ = size(separator_array);
            % Get direction
            direction=source.Value;
            % Set slider min, max
            minor_step = 1/(dim_(direction)-1);
            major_step = 0.1;
            set(Slider_axes_separator_2dview,'Min',1,'Max',dim_(direction),'SliderStep', [minor_step, major_step],'Value',1);
            % Update text
            set(Text_separator_slider,'String',['Position: ' num2str(1) '/' num2str(dim_(direction))]);
            % Update figure
            update_separtor_figure
        end
        
        function update_separtor_figure
            % Colors
            for current_phase=1:1:1
                RGB_phase.index(current_phase).rgb = [color_phase(current_phase,1) color_phase(current_phase,2) color_phase(current_phase,3)];
            end
            dim_ = size(separator_array);
            direction = Popup_separator_slider.Value;
            pos_ = round(Slider_axes_separator_2dview.Value);
            Position_slice(direction)=pos_;
            if direction==1
                % Initializaion
                slice_color = zeros(dim_(2),dim_(3),3); % RGB color map
                slice_r = zeros(dim_(2),dim_(3)); % Red color map
                slice_g = zeros(dim_(2),dim_(3)); % Green color map
                slice_b = zeros(dim_(2),dim_(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:1
                    code_tmp =current_phase; % Current phase code
                    slice_r(separator_array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(separator_array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(separator_array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            elseif direction==2
                % Initializaion
                slice_color = zeros(dim_(1),dim_(3),3); % RGB color map
                slice_r = zeros(dim_(1),dim_(3)); % Red color map
                slice_g = zeros(dim_(1),dim_(3)); % Green color map
                slice_b = zeros(dim_(1),dim_(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:1
                    code_tmp =current_phase; % Current phase code
                    slice_r(separator_array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(separator_array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(separator_array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            elseif direction==3
                % Initializaion
                slice_color = zeros(dim_(1),dim_(2),3); % RGB color map
                slice_r = zeros(dim_(1),dim_(2)); % Red color map
                slice_g = zeros(dim_(1),dim_(2)); % Green color map
                slice_b = zeros(dim_(1),dim_(2)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:1
                    code_tmp =current_phase; % Current phase code
                    slice_r(separator_array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(separator_array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(separator_array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            end
            slice_color(:,:,1)=slice_r; % Attribute RGB color
            slice_color(:,:,2)=slice_g;
            slice_color(:,:,3)=slice_b;
            % Display the slice
            image(slice_color,'parent',axes_2dview);
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
        end
        
        % Initialize
        update_separtor_figure
    end

%% FUNCTIONS
    function Checkdimension_forseparator
        inplane1 =[];
        inplane2 =[];
        if celldata(1).heterogeneous
            dim_ = size(celldata(1).array);
            inplane1=[inplane1 dim_(2)];
            inplane2=[inplane2 dim_(3)];
        end
        if celldata(3).heterogeneous
            dim_ = size(celldata(3).array);
            inplane1=[inplane1 dim_(2)];
            inplane2=[inplane2 dim_(3)];            
        end        
        % Check dimension identical
        unique_inplane1 = unique(inplane1);
        unique_inplane2 = unique(inplane2);
        min_inplane = [min(unique_inplane1) min(unique_inplane2)];
    end

%%
%% "CELL OPTIONS" TAB
%%

%% PRE-ASSIGNED PARAMETERS
Fig_visualizecell=0;
voxel_size_um = 1;

%% GUI
% Title
Text_tab_cell = uicontrol('Parent', tab_cell, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Cell options');

% Voxel size
Text_voxelsize = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Voxel size in um:',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.88 0.175 0.03]);
Edit_voxelsize = uicontrol('Parent', tab_cell,'Style', 'edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',num2str(voxel_size_um),'Visible','on',...
    'visible','on','enable','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.86 0.175 0.03],'Callback', @edit_voxelsize_Callback);

% Microstructures
Text_leftelectrode = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Left electrode',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.2 0.88 0.2667 0.03]);
Text_separator = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Separator',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.4667 0.88 0.2667 0.03]);
Text_rightelectrode = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','Right electrode',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.7333 0.88 0.2667 0.03]);

% Dimensions
Text_dimension = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','1) Check dimensions are matching',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.82 0.8 0.03]);
line_dimension = annotation(tab_cell,'line',[description_tab_fromleft 0.975],[0.82 0.82],'visible','on');
Text_dimension_TP = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Through-plane direction',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.78 0.175 0.03]);
Text_dimension_IP1 = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','In-plane direction 1',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.75 0.175 0.03]);
Text_dimension_IP2 = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','In-plane direction 2',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.72 0.175 0.03]);

Text_dimension_TP_left = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.2 0.78 0.2667 0.03]);
Text_dimension_IP1_left = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.2 0.75 0.2667 0.03]);
Text_dimension_IP2_left = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.2 0.72 0.2667 0.03]);
Text_dimension_TP_separator = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.4667 0.78 0.2667 0.03]);
Text_dimension_IP1_separator= uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.4667 0.75 0.2667 0.03]);
Text_dimension_IP2_separator = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.4667 0.72 0.2667 0.03]);
Text_dimension_TP_right = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.7333 0.78 0.2667 0.03]);
Text_dimension_IP1_right= uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.7333 0.75 0.2667 0.03]);
Text_dimension_IP2_right = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','not defined',...
    'visible','on','BackgroundColor','w','ForegroundColor','r','HorizontalAlignment','center','Units','normalized','Position', [0.7333 0.72 0.2667 0.03]);

Text_crop = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Crop in-plane dimensions of larger domains:',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.67 0.3 0.03]);
Checkbox_crop = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.325 0.675 0.1 0.03],'Callback', @checkbox_crop_Callback); 

% Morphology opening
Text_morphologyopening = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','2) Apply morphology opening',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.63 0.8 0.03]);
line_morphologyopening = annotation(tab_cell,'line',[description_tab_fromleft 0.975],[0.63 0.63],'visible','on');
Text_erosiondilatation = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Erosion + dilatation',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.59 0.175 0.03]);
Text_cluster = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Convert each phase in a unique cluster and clean it',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.55 0.25 0.04]);

Checkbox_erosiondilatation_left = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.325 0.595 0.1 0.03]); 
Checkbox_erosiondilatation_separator = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.6 0.595 0.1 0.03]); 
Checkbox_erosiondilatation_right = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.8666 0.595 0.1 0.03]); 
Checkbox_uniquecluster_left = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.325 0.565 0.1 0.03]); 
Checkbox_uniquecluster_separator = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.6 0.565 0.1 0.03]); 
Checkbox_uniquecluster_right = uicontrol('Parent', tab_cell,'Style', 'checkbox','Value',0,'Visible','on',...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.8666 0.565 0.1 0.03]); 

pushbutton_applymorphologyopening = uicontrol('Parent', tab_cell,'Style', 'pushbutton','Value',0,'String','Apply','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.5 0.15 0.04],'Callback', @applymorphologyopening_Callback);
Text_morphologyopening_left = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',' ',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.2 0.5 0.2667 0.03]);
Text_morphologyopening_separator= uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',' ',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.4667 0.5 0.2667 0.03]);
Text_morphologyopening_right = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',' ',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.7333 0.5 0.2667 0.03]);

% Phase
Text_phase = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','3) Choose phase id that will be used in the cell for each microstructure',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.455 0.8 0.03]);
line_phasecode = annotation(tab_cell,'line',[description_tab_fromleft 0.975],[0.455 0.455],'visible','on');
Text_phaseinstructions = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Choose integer >= 1 for cell id (do not choose 0). Use different value for each phase, even for the electrolyte of the anode and cathode.',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.395 0.95 0.04]);

table_phaseleft = uitable('Parent', tab_cell,'enable','on','Visible','on','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.195 0.28 0.25 0.12],'CellEditCallback', @table_cellid_Callback);
table_phaseleft.ColumnName = {'Phase id','Volume fraction','Cell id'}; % Column name
table_phaseleft.ColumnEditable = [false false true]; % Select column editable
table_phaseleft.ColumnWidth = {'auto', 'auto', 'auto'}; % Auto width
table_phaseleft.RowName = []; % Remove row name
table_phaseseparator = uitable('Parent', tab_cell,'enable','on','Visible','on','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.46 0.28 0.25 0.12],'CellEditCallback', @table_cellid_Callback);
table_phaseseparator.ColumnName = {'Phase id','Volume fraction','Cell id'}; % Column name
table_phaseseparator.ColumnEditable = [false false true]; % Select column editable
table_phaseseparator.ColumnWidth = {'auto', 'auto', 'auto'}; % Auto width
table_phaseseparator.RowName = []; % Remove row name
table_phaseright = uitable('Parent', tab_cell,'enable','on','Visible','on','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.725 0.28 0.25 0.12],'CellEditCallback', @table_cellid_Callback);
table_phaseright.ColumnName = {'Phase id','Volume fraction','Cell id'}; % Column name
table_phaseright.ColumnEditable = [false false true]; % Select column editable
table_phaseright.ColumnWidth = {'auto', 'auto', 'auto'}; % Auto width
table_phaseright.RowName = []; % Remove row name

% Assign cell
Text_assigncell = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','4) Assign phase',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.225 0.8 0.03]);
line_assigncell = annotation(tab_cell,'line',[description_tab_fromleft 0.975],[0.225 0.225],'visible','on');

Popup_cell_id = uicontrol('Parent', tab_cell,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', {'Cell id '},'Units','normalized','Position', [0.025 0.17 0.1 0.04],'enable','on','Visible','on');
Popup_cell_assign = uicontrol('Parent', tab_cell,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', {'Cell assignment','Left electrode: solid','Left electrode: electrolyte','Separator: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'},'Units','normalized','Position', [0.135 0.17 0.215 0.04],'enable','on','Visible','on');
pushbutton_validate_cellassign = uicontrol('Parent', tab_cell,'Style', 'pushbutton','Value',0,'String','Validate assignment','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.36 0.18 0.2 0.03],'Callback', @validatecellassign_Callback);
table_cell_assign = uitable('Parent', tab_cell,'enable','on','Visible','on','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.57 0.09 0.405 0.12]);
table_cell_assign.ColumnName = {'Cell id','Cell assignment'}; % Column name
table_cell_assign.ColumnEditable = [false false]; % Select column editable
table_cell_assign.ColumnWidth = {'auto', 'auto', 'auto'}; % Auto width
table_cell_assign.RowName = []; % Remove row name

% Assemble cell
Text_assemblecell = uicontrol('Parent', tab_cell, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','5) Assemble all microstructures to create cell',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.05 0.8 0.03]);
line_assemblecell = annotation(tab_cell,'line',[description_tab_fromleft 0.975],[0.05 0.05],'visible','on');
pushbutton_assemblecell = uicontrol('Parent', tab_cell,'Style', 'pushbutton','Value',0,'String','Assemble cell','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.005 0.25 0.04],'Callback', @assemblecell_Callback);
pushbutton_visualizecell = uicontrol('Parent', tab_cell,'Style', 'pushbutton','Value',0,'String','Visualize cell','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.375 0.005 0.25 0.04],'Callback', @visualizecell_Callback);
pushbutton_savecell = uicontrol('Parent', tab_cell,'Style', 'pushbutton','Value',0,'String','Save cell and go on','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.725 0.005 0.25 0.04],'Callback', @savecell_Callback);

%% CALLBACKS
    function edit_voxelsize_Callback(~,~)
        voxel_size_um = str2num(Edit_voxelsize.String);
    end

    function table_cellid_Callback(~,~)
        Popup_cell_assign.Value=1;
        enable_validate = 'on';
        left_ = cell2mat(table_phaseleft.Data(:,3));
        separator = cell2mat(table_phaseseparator.Data(:,3));
        right_ = cell2mat(table_phaseright.Data(:,3));
        cell_id = unique([left_; separator; right_]);
        if length(unique(cell_id))==4
            if length(unique(left_))==1 && length(unique(right_))==2
                str_region = {'Left electrode: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'};
            elseif length(unique(left_))==2 && length(unique(right_))==1
                str_region = {'Left electrode: solid','Left electrode: electrolyte','Separator: electrolyte','Right electrode: solid'};
            end
        elseif length(cell_id)==5
            if length(unique(left_))==2 && length(unique(right_))==2
                str_region = {'Left electrode: solid','Left electrode: electrolyte','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'};
            elseif length(unique(left_))==1 && length(unique(separator))==2 && length(unique(right_))==2
                str_region = {'Left electrode: solid','Separator: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'};
            end
        elseif length(cell_id)==6
             str_region = {'Left electrode: solid','Left electrode: electrolyte','Separator: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'};
        else
            str_region = {'Cell id values are incorrect'};
            enable_validate = 'off';
        end
        Popup_cell_assign.String = ['Cell assignment' str_region];
        Popup_cell_id.String = [{'Cell id '} num2cell(cell_id')]; % Update GUI string
        set(pushbutton_validate_cellassign,'enable',enable_validate);
        table_cell_assign.Data=[]; % Reset table
        
    end

    function validatecellassign_Callback(~,~)
        cell_id_string = cell2mat(Popup_cell_id.String(2:end));
        cell_id_value = Popup_cell_id.Value-1;
        cell_assignment_string = Popup_cell_assign.String(2:end);
        cell_assignment_value = Popup_cell_assign.Value-1;
        if cell_id_value>=1 && cell_assignment_value>=1
            val_ = table_cell_assign.Data;
            if isempty(val_)
                val_ = [num2cell(cell_id_string(cell_id_value)) cell_assignment_string(cell_assignment_value)];
            else
                if sum(unique(str2num(cell2mat(val_(:,1)))) == str2num(cell_id_string(cell_id_value)))== 0 % Not yet assigned
                    if ~ismember(cell_assignment_string(cell_assignment_value),unique(val_(:,2))) % Not yet assigned
                        new_line = [num2cell(cell_id_string(cell_id_value)) cell_assignment_string(cell_assignment_value)];
                        val_=[val_ ; new_line];
                    end
                end
            end
            table_cell_assign.Data= val_; % Update table
            if length(table_cell_assign.Data) == length(cell_id_string)
                set(pushbutton_assemblecell,'enable','on');
            else
                set(pushbutton_assemblecell,'enable','off');
                set(pushbutton_visualizecell,'enable','off');
                set(pushbutton_savecell,'enable','off');
            end
        end
    end

    function savecell_Callback(~,~)
        % Save folder
        if ~exist(options.save.mainfolder,'dir')
            mkdir(options.save.mainfolder);
        end
        if options.save.microstructure_tif
            if celldata(1).heterogeneous
                function_save_tif(celldata(1).array, [options.save.mainfolder 'Array_leftelectrode.tif'])
            end
            if celldata(2).heterogeneous
                function_save_tif(celldata(2).array, [options.save.mainfolder 'Array_separator.tif'])
            end
            if celldata(3).heterogeneous
                function_save_tif(celldata(3).array, [options.save.mainfolder 'Array_rightelectrode.tif'])
            end
        end
        if options.save.microstructure_mat
            if celldata(1).heterogeneous
                tmp = celldata(1).array;
                save([options.save.mainfolder 'Array_leftelectrode.mat'],'tmp')
            end
            if celldata(2).heterogeneous
                tmp = celldata(2).array;
                save([options.save.mainfolder 'Array_separator.mat'],'tmp')
            end
            if celldata(3).heterogeneous
                tmp = celldata(3).array;
                save([options.save.mainfolder 'Array_rightelectrode.mat'],'tmp')
            end
            clear tmp;
        end
                
        if options.save.cell_tif
            function_save_tif(cell_array, [options.save.mainfolder 'Array_cell.tif'])
        end
        if options.save.cell_mat
            save([options.save.mainfolder 'Array_cell.mat'],'cell_array')
        end
        celldata(4).exist = true;
        celldata(4).array = cell_array;
        if celldata(4).iso2mesh_options_set
            set(pushbutton_generatemesh,'enable','on');
        end
        
        % Information below help determining boundary conditions location in FEniCS
        if ~exist([options.save.mainfolder options.save.meshdatafolder],'dir')
            mkdir([options.save.mainfolder options.save.meshdatafolder]);
        end
        tmp = celldata(4).Microstructure_dimension;
        save([options.save.mainfolder options.save.meshdatafolder 'Microstructure_dimension.txt'],'tmp','-ascii');
        tmp = celldata(4).Separator_bounds;
        save([options.save.mainfolder options.save.meshdatafolder 'Separator_bounds.txt'],'tmp','-ascii');
        tmp = celldata(4).leftelectrode_min_max;
        save([options.save.mainfolder options.save.meshdatafolder 'leftelectrode_coor_min_max.txt'],'tmp','-ascii');
        tmp = celldata(4).electrolyte_min_max;
        save([options.save.mainfolder options.save.meshdatafolder 'electrolyte_coor_min_max.txt'],'tmp','-ascii');
        tmp = celldata(4).rightelectrode_min_max ;
        save([options.save.mainfolder options.save.meshdatafolder 'rightelectrode_coor_min_max.txt'],'tmp','-ascii');
    end

    function visualizecell_Callback(~,~)
        if ismember(findall(0,'type','figure'),Fig_visualizecell)
            close Figure_visualizemicrostructure;
        end
        Fig_visualizecell = figure; % Create figure
        Fig_visualizecell.Name= 'Visualize cell';
        Fig_visualizecell.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig_visualizecell,'position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]); % Full screen figure
        % Figure for volume view
        axes_cell_2dview = axes('Parent', Fig_visualizecell,'Visible','off','FontName',font_name_GUI,'Units','normalized','Position', [0 0.1 1 0.85]);
        % Remove tick and label
        set(axes_cell_2dview,'xtick',[],'ytick',[]);
        set(axes_cell_2dview,'xticklabel',[],'yticklabel',[]);
        % Box on
        box(axes_cell_2dview,'on');
        % Fit the axes box
        axis(axes_cell_2dview,'tight');
        % Aspect ratio is 1:1
        axis(axes_cell_2dview,'equal');
        
        % Slider
        Slider_axes_cell_2dview = uicontrol('Parent', Fig_visualizecell,'Style', 'slider','Min',1,'Max',100,'Value',1,'Units','normalized','Position', [0.3 0.025 0.4 0.04],'Callback', @slider_axes2dview_cell_Callback,...
            'Visible','on','enable','on');
        % Text (slider position)
        Text_cell_slider = uicontrol('Parent', Fig_visualizecell,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Position: -/-','Visible','on','enable','on',...
            'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.725 0.022 0.15 0.04]);
        % Popup menu (slider direction)
        Popup_cell_slider = uicontrol('Parent', Fig_visualizecell,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
            'String', {'Through-plane direction','In-plane direction 1','In-plane direction 2'},'Value',1,'Units','normalized','Position', [0.05 0.025 0.215 0.04],'enable','on','Visible','on','Callback', @popup_axes2dview_cell_Callback);
        
        dim_ = size(cell_array);
        direction = 1; pos_=1;
        set(Text_cell_slider,'String',['Position: ' num2str(pos_) '/' num2str(dim_(direction))]);
        minor_step = 1/(dim_(direction)-1);
        major_step = 0.1;
        set(Slider_axes_cell_2dview,'Min',1,'Max',dim_(direction),'SliderStep', [minor_step, major_step],'Value',1);
                        
        % Slider
        function slider_axes2dview_cell_Callback(source,~)
            dim_ = size(cell_array);
            % Get position value
            pos_=round(source.Value);
            direction = Popup_cell_slider.Value;
            set(Text_cell_slider,'String',['Position: ' num2str(pos_) '/' num2str(dim_(direction))]); % Update text
            update_cellfigure % Update figure
        end
        
        % Select direction
        function popup_axes2dview_cell_Callback(source,~)
            dim_ = size(cell_array);
            % Get direction
            direction=source.Value;
            % Set slider min, max
            minor_step = 1/(dim_(direction)-1);
            major_step = 0.1;
            set(Slider_axes_cell_2dview,'Min',1,'Max',dim_(direction),'SliderStep', [minor_step, major_step],'Value',1);
            % Update text
            set(Text_cell_slider,'String',['Position: ' num2str(1) '/' num2str(dim_(direction))]);
            % Update figure
            update_cellfigure
        end
        
        function update_cellfigure
            % Colors
            phase_ = unique(cell_array);
            n_ = length(phase_);
            dim_ = size(cell_array);
            for current_phase=1:1:n_
                RGB_phase.index(current_phase).rgb = [color_phase(current_phase,1) color_phase(current_phase,2) color_phase(current_phase,3)];
            end
            direction = Popup_cell_slider.Value;
            pos_ = round(Slider_axes_cell_2dview.Value);
            Position_slice(direction)=pos_;
            if direction==1
                % Initializaion
                slice_color = zeros(dim_(2),dim_(3),3); % RGB color map
                slice_r = zeros(dim_(2),dim_(3)); % Red color map
                slice_g = zeros(dim_(2),dim_(3)); % Green color map
                slice_b = zeros(dim_(2),dim_(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:n_
                    code_tmp =phase_(current_phase); % Current phase code
                    slice_r(cell_array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(cell_array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(cell_array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            elseif direction==2
                % Initializaion
                slice_color = zeros(dim_(1),dim_(3),3); % RGB color map
                slice_r = zeros(dim_(1),dim_(3)); % Red color map
                slice_g = zeros(dim_(1),dim_(3)); % Green color map
                slice_b = zeros(dim_(1),dim_(3)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:n_
                    code_tmp =phase_(current_phase); % Current phase code
                    slice_r(cell_array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(cell_array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(cell_array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            elseif direction==3
                % Initializaion
                slice_color = zeros(dim_(1),dim_(2),3); % RGB color map
                slice_r = zeros(dim_(1),dim_(2)); % Red color map
                slice_g = zeros(dim_(1),dim_(2)); % Green color map
                slice_b = zeros(dim_(1),dim_(2)); % Blue color map
                % Attribute RGB colors for each voxel
                for current_phase=1:1:n_
                    code_tmp =phase_(current_phase); % Current phase code
                    slice_r(cell_array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                    slice_g(cell_array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                    slice_b(cell_array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                end
            end
            slice_color(:,:,1)=slice_r; % Attribute RGB color
            slice_color(:,:,2)=slice_g;
            slice_color(:,:,3)=slice_b;
            % Display the slice
            image(slice_color,'parent',axes_cell_2dview);
            % Get axis position
            axe_position =axes_cell_2dview.Position;
            set(axes_cell_2dview ,'Position', axe_position);
            % Remove tick and label
            set(axes_cell_2dview,'xtick',[],'ytick',[]);
            set(axes_cell_2dview,'xticklabel',[],'yticklabel',[]);
            % Fit the axes box
            axis(axes_cell_2dview,'tight');
            % Aspect ratio is 1:1
            axis(axes_cell_2dview,'equal');
        end
        % Initialize
        update_cellfigure
    end

    function assemblecell_Callback(~,~)
        % In-plane dimensions
        n_2=min_inplane(1);
        n_3=min_inplane(2);
        % Initialize
        if celldata(1).heterogeneous
            dim_ = size(celldata(1).array);
            n_left = dim_(1);
        else
            n_left = celldata(1).length_voxel;
            celldata(1).array = zeros(n_left,n_2,n_3);
        end
        if celldata(2).heterogeneous
            dim_ = size(celldata(2).array);
            n_separator = dim_(1);
        else
            n_separator = celldata(2).length_voxel;
            celldata(2).array = zeros(n_separator,n_2,n_3);            
        end
        if celldata(3).heterogeneous
            dim_ = size(celldata(3).array);
            n_right = dim_(1);
        else
            n_right = celldata(3).length_voxel;
            celldata(3).array = zeros(n_right,n_2,n_3);            
        end        
        
        % Dimension of the half cell
        n_1=n_left+n_separator+n_right;
        domain_size_cell=[n_1, n_2, n_3];
        
        % Initialisation
        cell_array = zeros(domain_size_cell);
        
        % Left electrode
        left_ = celldata(1).array*0;
        phase_id = [double(cell2mat(table_phaseleft.Data(:,1))) double(cell2mat(table_phaseleft.Data(:,2))) double(cell2mat(table_phaseleft.Data(:,3)))];
        %phase_id = cell2mat(table_phaseleft.Data);
        [n_phase,~]=size(phase_id);
        for k=1:1:n_phase
            left_ (celldata(1).array ==  phase_id(k,1)) = phase_id(k,3);
        end
        cell_array(1:n_left,:,:) = left_;
        clear left_;
        
        % Separator
        Separator_ = celldata(2).array*0;
        phase_id = [double(cell2mat(table_phaseseparator.Data(:,1))) double(cell2mat(table_phaseseparator.Data(:,2))) double(cell2mat(table_phaseseparator.Data(:,3)))];
        %phase_id = cell2mat(table_phaseseparator.Data);
        [n_phase,~]=size(phase_id);
        for k=1:1:n_phase
            Separator_ (celldata(2).array ==  phase_id(k,1)) = phase_id(k,3);
        end
        cell_array(n_left+1:n_left+n_separator,:,:) = Separator_;
        clear Separator_;    
        
        % Right electrode
        Right_ = celldata(3).array*0;
        
        phase_id = [double(cell2mat(table_phaseright.Data(:,1))) double(cell2mat(table_phaseright.Data(:,2))) double(cell2mat(table_phaseright.Data(:,3)))];
        %phase_id = cell2mat(table_phaseright.Data);
        [n_phase,~]=size(phase_id);
        for k=1:1:n_phase
            Right_ (celldata(3).array ==  phase_id(k,1)) = phase_id(k,3);
        end
        cell_array(n_left+n_separator+1:end,:,:) = Right_;
        clear Right_;            
                
        % % FENICS INFORMATION
        % Dimension
        dim_ = size(cell_array);
        celldata(4).Microstructure_dimension=[voxel_size_um; dim_(1); dim_(2); dim_(3)]; % Dimension
        % Separator bounds
        celldata(4).Separator_bounds = [n_left; n_left+n_separator];
        % Min-max coordinates of each domain
        % cell_id = table_cell_assign.Data(:,1);
        cell_assignment = table_cell_assign.Data(:,2);
        
        % Left electrode: solid domain
        codes = find(strcmp(cell_assignment,'Left electrode: solid'));
        [celldata(4).leftelectrode_min_max, celldata(4).left_binary_array] = find_extremums_binaryarray(cell_array, codes);
        
        % All electrolyte domains
        codes=[-1 -1 -1];
        tmp = find(strcmp(cell_assignment,'Left electrode: electrolyte'));
        if ~isempty(tmp)
            codes(1) = tmp;
        end
        tmp = find(strcmp(cell_assignment,'Separator: electrolyte'));
        if ~isempty(tmp)
            codes(2) = tmp;
        end
        tmp = find(strcmp(cell_assignment,'Right electrode: electrolyte'));
        if ~isempty(tmp)
            codes(3) = tmp;
        end
        codes(codes==-1)=[]; % Remove -1
        [celldata(4).electrolyte_min_max, celldata(4).electrolyte_binary_array]= find_extremums_binaryarray(cell_array, codes);
        
        % Right electrode: solid domain
        codes = find(strcmp(cell_assignment,'Right electrode: solid'));
        [celldata(4).rightelectrode_min_max, celldata(4).right_binary_array] = find_extremums_binaryarray(cell_array, codes);

        % Update GUI
        set(pushbutton_visualizecell,'enable','on');
        set(pushbutton_savecell,'enable','on');     
    end

    function checkbox_crop_Callback(~,~)
        if celldata(1).heterogeneous
            [y0, y1] = center_cropped (celldata(1).array, 2, min_inplane(1));
            [z0, z1] = center_cropped (celldata(1).array, 3, min_inplane(2));
            celldata(1).array = celldata(1).array(:,y0:y1,z0:z1);
        end
        if celldata(2).heterogeneous
            [y0, y1] = center_cropped (celldata(2).array, 2, min_inplane(1));
            [z0, z1] = center_cropped (celldata(2).array, 3, min_inplane(2));
            celldata(2).array = celldata(2).array(:,y0:y1,z0:z1);
        end
        if celldata(3).heterogeneous
            [y0, y1] = center_cropped (celldata(3).array, 2, min_inplane(1));
            [z0, z1] = center_cropped (celldata(3).array, 3, min_inplane(2));
            celldata(3).array = celldata(3).array(:,y0:y1,z0:z1);
        end
        Checkdimension_cell % Update
        phaseid_cell % Update
        set(Checkbox_crop,'enable','off');
    end

    function applymorphologyopening_Callback(~,~)
        % Erosion + dilatation
        if celldata(1).heterogeneous && Checkbox_erosiondilatation_left.Value
            [celldata(1).array] = Function_erosion_dilatation(celldata(1).array,celldata(1).codesolid ,celldata(1).codepore );
        end
        if celldata(2).heterogeneous && Checkbox_erosiondilatation_separator.Value
            [celldata(2).array] = Function_erosion_dilatation(celldata(2).array,celldata(2).codesolid ,celldata(2).codepore );
        end        
        if celldata(3).heterogeneous && Checkbox_erosiondilatation_right.Value
            [celldata(3).array] = Function_erosion_dilatation(celldata(3).array,celldata(3).codesolid ,celldata(3).codepore );
        end
        % Convert solid and pore cluster(s) in a unique solid and unique pore cluster + clean microstructure
        clean_connection=1;
        if celldata(1).heterogeneous && Checkbox_uniquecluster_left.Value
            dim_ = size(celldata(1).array);
            array_tmp=Function_cleanmicrostructure(celldata(1).array,celldata(1).codesolid,celldata(1).codepore,clean_connection);
            unique_=unique(array_tmp);
            if length(unique_)==2
                % Create binary matrix
                celldata(1).array=zeros(dim_(1),dim_(2),dim_(3)) + celldata(1).codepore;
                celldata(1).array(array_tmp == celldata(1).codesolid + 1) = celldata(1).codesolid;
                clear array_tmp;
            else
                disp 'ERROR: THERE ARE STILL SOME NON ISOLATED CLUSTERS IN THE MICROSTRUCTURE';
            end
        end        
        if celldata(2).heterogeneous && Checkbox_uniquecluster_separator.Value
            dim_ = size(celldata(2).array);
            array_tmp=Function_cleanmicrostructure(celldata(2).array,celldata(2).codesolid,celldata(2).codepore,clean_connection);
            unique_=unique(array_tmp);
            if length(unique_)==2
                % Create binary matrix
                celldata(2).array=zeros(dim_(1),dim_(2),dim_(3)) + celldata(2).codepore;
                celldata(2).array(array_tmp == celldata(2).codesolid + 1) = celldata(2).codesolid;
                clear array_tmp;
            else
                disp 'ERROR: THERE ARE STILL SOME NON ISOLATED CLUSTERS IN THE MICROSTRUCTURE';
            end
        end     
        if celldata(3).heterogeneous && Checkbox_uniquecluster_right.Value
            dim_ = size(celldata(3).array);
            array_tmp=Function_cleanmicrostructure(celldata(3).array,celldata(3).codesolid,celldata(3).codepore,clean_connection);
            unique_=unique(array_tmp);
            if length(unique_)==2
                % Create binary matrix
                celldata(3).array=zeros(dim_(1),dim_(2),dim_(3)) + celldata(3).codepore;
                celldata(3).array(array_tmp == celldata(3).codesolid + 1) = celldata(3).codesolid;
                clear array_tmp;
            else
                disp 'ERROR: THERE ARE STILL SOME NON ISOLATED CLUSTERS IN THE MICROSTRUCTURE';
            end
        end           
        phaseid_cell % Update
    end

%% FUNCTIONS
    function [] = phaseid_cell()
        if celldata(1).heterogeneous
            Phase_=unique(celldata(1).array); % Get phase code
            number_=length(Phase_); % Get number of phase
            total_number_voxel=numel(celldata(1).array);
            volumefraction=zeros(number_,1);
            for n=1:1:number_
                volumefraction(n,1)=sum(sum(sum(celldata(1).array==Phase_(n))))/total_number_voxel;
            end
            Volumefraction_data = [num2cell(Phase_) num2cell(volumefraction) num2cell(Phase_)];
            table_phaseleft.Data=Volumefraction_data;
        else
            Volumefraction_data = [num2cell(0) num2cell(1) num2cell(0)];
            table_phaseleft.Data=Volumefraction_data;            
        end
        if celldata(2).heterogeneous
            Phase_=unique(celldata(2).array); % Get phase code
            number_=length(Phase_); % Get number of phase
            total_number_voxel=numel(celldata(2).array);
            volumefraction=zeros(number_,1);
            for n=1:1:number_
                volumefraction(n,1)=sum(sum(sum(celldata(2).array==Phase_(n))))/total_number_voxel;
            end
            Volumefraction_data = [num2cell(Phase_) num2cell(volumefraction) num2cell(Phase_)];
            table_phaseseparator.Data=Volumefraction_data;
        else
            Volumefraction_data = [num2cell(0) num2cell(1) num2cell(0)];
            table_phaseseparator.Data=Volumefraction_data;
        end
        if celldata(3).heterogeneous
            Phase_=unique(celldata(3).array); % Get phase code
            number_=length(Phase_); % Get number of phase
            total_number_voxel=numel(celldata(3).array);
            volumefraction=zeros(number_,1);
            for n=1:1:number_
                volumefraction(n,1)=sum(sum(sum(celldata(3).array==Phase_(n))))/total_number_voxel;
            end
            Volumefraction_data = [num2cell(Phase_) num2cell(volumefraction) num2cell(Phase_)];
            table_phaseright.Data=Volumefraction_data;
        else
            Volumefraction_data = [num2cell(0) num2cell(1) num2cell(0)];
            table_phaseright.Data=Volumefraction_data;             
        end
    end

    function [x0, x1] = center_cropped (array_, dimension, threshold)
        dim_ = size(array_);
        if dim_(dimension) == threshold
            x0=1; x1=dim_(dimension);
        else
            dx = dim_(dimension)-threshold;
            x0=round(1+dx/2); x1=threshold+x0-1;
        end
    end
                
    function Checkdimension_cell
        set(Checkbox_crop,'enable','on');
        inplane1 =[];
        inplane2 =[];
        if celldata(1).heterogeneous
            dim_ = size(celldata(1).array);
            set(Text_dimension_TP_left,'String',num2str(dim_(1)),'ForegroundColor','k');
            Text_dimension_IP1_left.String=num2str(dim_(2));
            Text_dimension_IP2_left.String=num2str(dim_(3));
            inplane1=[inplane1 dim_(2)];
            inplane2=[inplane2 dim_(3)];
        else
            set(Text_dimension_TP_left,'String',num2str(celldata(1).length_voxel),'ForegroundColor','k');
            set(Text_dimension_IP1_left,'String','n/a','ForegroundColor','k');
            set(Text_dimension_IP2_left,'String','n/a','ForegroundColor','k');                
        end
        if celldata(2).heterogeneous
            dim_ = size(celldata(2).array);
            set(Text_dimension_TP_separator,'String',num2str(dim_(1)),'ForegroundColor','k');
            Text_dimension_IP1_separator.String=num2str(dim_(2));
            Text_dimension_IP2_separator.String=num2str(dim_(3));
            inplane1=[inplane1 dim_(2)];
            inplane2=[inplane2 dim_(3)];
        else
            set(Text_dimension_TP_separator,'String',num2str(celldata(2).length_voxel),'ForegroundColor','k');
            set(Text_dimension_IP1_separator,'String','n/a','ForegroundColor','k');
            set(Text_dimension_IP2_separator,'String','n/a','ForegroundColor','k');     
        end
        if celldata(3).heterogeneous
            dim_ = size(celldata(3).array);
            set(Text_dimension_TP_right,'String',num2str(dim_(1)),'ForegroundColor','k');
            Text_dimension_IP1_right.String=num2str(dim_(2));
            Text_dimension_IP2_right.String=num2str(dim_(3));
            inplane1=[inplane1 dim_(2)];
            inplane2=[inplane2 dim_(3)];            
        else
            set(Text_dimension_TP_right,'String',num2str(celldata(3).length_voxel),'ForegroundColor','k');
            set(Text_dimension_IP2_right,'String','n/a','ForegroundColor','k');
            set(Text_dimension_IP1_right,'String','n/a','ForegroundColor','k');              
        end        
        
        % Check dimension identical
        unique_inplane1 = unique(inplane1);
        unique_inplane2 = unique(inplane2);
        incorrect_dimension = false;
        min_inplane = [min(unique_inplane1) min(unique_inplane2)];
        if length(unique_inplane1)~=1
            incorrect_dimension = true;
            if celldata(1).heterogeneous
                set(Text_dimension_IP1_left,'ForegroundColor','r');
            end
            if celldata(2).heterogeneous
                set(Text_dimension_IP1_separator,'ForegroundColor','r');
            end
            if celldata(3).heterogeneous
                set(Text_dimension_IP1_right,'ForegroundColor','r');
            end
        else
            if celldata(1).heterogeneous
                set(Text_dimension_IP1_left,'ForegroundColor','k');
            end
            if celldata(2).heterogeneous
                set(Text_dimension_IP1_separator,'ForegroundColor','k');
            end
            if celldata(3).heterogeneous
                set(Text_dimension_IP1_right,'ForegroundColor','k');
            end            
        end
        if length(unique_inplane2)~=1
            incorrect_dimension = true;
            if celldata(1).heterogeneous
                set(Text_dimension_IP2_left,'ForegroundColor','r');
            end
            if celldata(2).heterogeneous
                set(Text_dimension_IP2_separator,'ForegroundColor','r');
            end
            if celldata(3).heterogeneous
                set(Text_dimension_IP2_right,'ForegroundColor','r');
            end
        else
            if celldata(1).heterogeneous
                set(Text_dimension_IP2_left,'ForegroundColor','k');
            end
            if celldata(2).heterogeneous
                set(Text_dimension_IP2_separator,'ForegroundColor','k');
            end
            if celldata(3).heterogeneous
                set(Text_dimension_IP2_right,'ForegroundColor','k');
            end   
        end        

        if incorrect_dimension
            set(Checkbox_crop,'enable','on');
            enable_ = 'off';
        else
            set(Checkbox_crop,'enable','off');
            enable_ = 'on';
        end
        set(Checkbox_erosiondilatation_left,'enable',enable_);
        set(Checkbox_erosiondilatation_separator,'enable',enable_);        
        set(Checkbox_erosiondilatation_right,'enable',enable_);
        set(Checkbox_uniquecluster_left,'enable',enable_);
        set(Checkbox_uniquecluster_separator,'enable',enable_);
        set(Checkbox_uniquecluster_right,'enable',enable_);
        set(pushbutton_applymorphologyopening,'enable',enable_);
        if celldata(1).heterogeneous==false
            set(Checkbox_erosiondilatation_left,'enable','off','visible','off');        
            set(Checkbox_uniquecluster_left,'enable','off','visible','off');   
        end
        if celldata(2).heterogeneous==false
            set(Checkbox_erosiondilatation_separator,'enable','off','visible','off');        
            set(Checkbox_uniquecluster_separator,'enable','off','visible','off'); 
        end
        if celldata(3).heterogeneous==false
            set(Checkbox_erosiondilatation_right,'enable','off','visible','off');        
            set(Checkbox_uniquecluster_right,enable','off','visible','off');     
        end        
    end

    function [extremums, array_binary] = find_extremums_binaryarray (array_,codes)
        n_ = length(codes);
        all_index = []; % Initialize
        for iter=1:1:n_
            current_code = codes(iter);
            index_=find(array_==current_code); % Find all index
            all_index = [all_index; index_];
        end
        [ix,iy,iz]=ind2sub(size(array_),all_index); % Convert to matrix index
        x_min=min(ix); y_min=min(iy); z_min=min(iz); % Take min and max
        x_max=max(ix); y_max=max(iy); z_max=max(iz);
        % % Binary array
        array_binary = zeros(size(array_));
        array_binary(all_index)=1; 
        % % Extremums
        % Convert to coordinates
        x_min=x_min-1; y_min=y_min-1; z_min=z_min-1;
        extremums=[x_min, y_min, z_min; x_max, y_max, z_max];
    end

%%
%% "ISO2MESH OPTIONS" TAB
%%
% Smoothing parameter
% For actual microstructure: iteration_smoothing=2 with lowpass
% For generated spheres-based microstructure: iteration_smoothing=5 with laplacian

%% PRE-ASSIGNED PARAMETERS
choice_surfacemesh=[" " "cgalmesh"];
choice_smoothsurfacemesh=[" " "lowpass" "laplacian" "laplacianhc"];

%% GUI
% Title
Text_tab_iso2mesh = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Iso2mesh options');

% Iso2mesh version
Text_iso2mesh_version = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',['Iso2mesh version: ' iso2meshver],...
    'enable','off','visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.01 0.2 0.03]);

% Surface mesh
Text_surfacemesh = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','1) Surface mesh from 3D array',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.88 0.8 0.03]);
line_surfacemesh = annotation(tab_Iso2mesh,'line',[description_tab_fromleft 0.975],[0.88 0.88],'visible','on');
Text_surfacemeshmethod = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Choose method:',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.82 0.25 0.04]);
Popupmenu_choice_surfacemesh = uicontrol('Parent', tab_Iso2mesh, 'Style', 'popupmenu', 'String', choice_surfacemesh,...
    'visible','on','enable','on','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.80 0.25 0.04],...
    'Callback',{@choice_surfacemeshmethod_Callback});
table_surfamesh_options= uitable('Parent', tab_Iso2mesh,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.35 0.66 0.625 0.2]);
table_surfamesh_options.ColumnName = {'Parameter','Description                                                                                                          ','Value'}; % Column name
table_surfamesh_options.ColumnEditable = [false false true]; % Select column editable
table_surfamesh_options.ColumnWidth = {'auto', 'auto','auto'}; % Auto width
table_surfamesh_options.RowName = []; % Remove row name

% Smooth surface mesh
Text_smoothsurfacemesh = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','2) Smooth surface mesh',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.6 0.8 0.03]);
line_smoothsurfacemesh = annotation(tab_Iso2mesh,'line',[description_tab_fromleft 0.975],[0.6 0.6],'visible','on');
Text_smoothsurfacemeshmethod = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String','Choose method:',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.025 0.54 0.25 0.04]);
Popupmenu_choice_smoothsurfacemesh = uicontrol('Parent', tab_Iso2mesh, 'Style', 'popupmenu', 'String', choice_smoothsurfacemesh,...
    'visible','on','enable','on','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.025 0.52 0.25 0.04],...
    'Callback',{@choice_smoothsurfacemeshmethod_Callback});
table_smoothmesh_options= uitable('Parent', tab_Iso2mesh,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.35 0.38 0.625 0.2]);
table_smoothmesh_options.ColumnName = {'Parameter','Description                                                                                                          ','Value'}; % Column name
table_smoothmesh_options.ColumnEditable = [false false true]; % Select column editable
table_smoothmesh_options.ColumnWidth = {'auto', 'auto','auto'}; % Auto width
table_smoothmesh_options.RowName = []; % Remove row name

% Volume mesh
Text_volumemesh = uicontrol('Parent', tab_Iso2mesh, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','3) Volumetric mesh',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.32 0.8 0.03]);
line_volumemesh = annotation(tab_Iso2mesh,'line',[description_tab_fromleft 0.975],[0.32 0.32],'visible','on');
table_volumemesh_options= uitable('Parent', tab_Iso2mesh,'enable','off','Visible','off','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.35 0.1 0.625 0.2]);
table_volumemesh_options.ColumnName = {'Parameter','Description                                                                                                          ','Value'}; % Column name
table_volumemesh_options.ColumnEditable = [false false true]; % Select column editable
table_volumemesh_options.ColumnWidth = {'auto', 'auto','auto'}; % Auto width
table_volumemesh_options.RowName = []; % Remove row name

%% CALLBACKS
function choice_smoothsurfacemeshmethod_Callback(~,~)
    if Popupmenu_choice_smoothsurfacemesh.Value==1
        options.method_smoothsurfacemesh='none';
        turn_GUI_smoothsurfacemesh_lowpass('off');
        turn_GUI_smoothsurfacemesh_laplacian('off');
        turn_GUI_smoothsurfacemesh_laplacianhc('off');
    end    
    if Popupmenu_choice_smoothsurfacemesh.Value==2
        options.method_smoothsurfacemesh='lowpass';
        turn_GUI_smoothsurfacemesh_laplacian('off');
        turn_GUI_smoothsurfacemesh_laplacianhc('off');
        turn_GUI_smoothsurfacemesh_lowpass('on');        
    end
    if Popupmenu_choice_smoothsurfacemesh.Value==3
        options.method_smoothsurfacemesh='laplacian';
        turn_GUI_smoothsurfacemesh_lowpass('off');
        turn_GUI_smoothsurfacemesh_laplacianhc('off');
        turn_GUI_smoothsurfacemesh_laplacian('on');
    end
    if Popupmenu_choice_smoothsurfacemesh.Value==4
        options.method_smoothsurfacemesh='laplacianhc';
        turn_GUI_smoothsurfacemesh_lowpass('off');
        turn_GUI_smoothsurfacemesh_laplacian('off');
        turn_GUI_smoothsurfacemesh_laplacianhc('on');
    end
    if Popupmenu_choice_smoothsurfacemesh.Value>1 && Popupmenu_choice_surfacemesh.Value>1
        turn_GUI_volumetricmesh('on');
    else
        turn_GUI_volumetricmesh('off');
    end
    if Popupmenu_choice_surfacemesh.Value>1 && Popupmenu_choice_smoothsurfacemesh.Value>1
        celldata(4).iso2mesh_options_set = true;
        if celldata(4).exist
            set(pushbutton_generatemesh,'enable','on');
        else
            set(pushbutton_generatemesh,'enable','off');
        end
    else
        set(pushbutton_generatemesh,'enable','off');
    end
end

function choice_surfacemeshmethod_Callback(~,~)
    if Popupmenu_choice_surfacemesh.Value==1
        options.method_surfacemesh='none';
        turn_GUI_surfacemesh_cgalmesh('off');
    end    
    if Popupmenu_choice_surfacemesh.Value==2
        options.method_surfacemesh='cgalmesh';
        turn_GUI_surfacemesh_cgalmesh('on');
    end
    if Popupmenu_choice_smoothsurfacemesh.Value>1 && Popupmenu_choice_surfacemesh.Value>1
        turn_GUI_volumetricmesh('on');
    else
        turn_GUI_volumetricmesh('off');
    end
    if Popupmenu_choice_surfacemesh.Value>1 && Popupmenu_choice_smoothsurfacemesh.Value>1
        celldata(4).iso2mesh_options_set = true;
        if celldata(4).exist
            set(pushbutton_generatemesh,'enable','on');
        else
            set(pushbutton_generatemesh,'enable','off');
        end
    else
        set(pushbutton_generatemesh,'enable','off');
    end
end

%% FUNCTIONS
    function turn_GUI_volumetricmesh(value)
        set(table_volumemesh_options,'enable',value,'visible',value);
        str_keepratio = 'Percentage of elements being kept after the simplification';
        str_maxvol = 'Maximum tetrahedra element volume';
        data_ = [{'keepratio';'maxvol'} {str_keepratio;str_maxvol} num2cell([1;10])];
        table_volumemesh_options.Data=data_;
    end

    function turn_GUI_surfacemesh_cgalmesh(value)
        set(table_surfamesh_options,'enable',value,'visible',value);
        str_radbound = 'Max radius of the Delaunay sphere (control the facet size at the surface)';
        str_distbound = 'Maximum deviation from the specified isosurfaces';
        str_isovalues = 'Isovalues where the levelset is defined';
        data_ = [{'radbound';'distbound';'Isovalues'} {str_radbound;str_distbound;str_isovalues} num2cell([1;1;1])];
        table_surfamesh_options.Data=data_;
    end

    function turn_GUI_smoothsurfacemesh_lowpass(value)
        set(table_smoothmesh_options,'enable',value,'visible',value);
        str_iteration = 'Number of smoothing iteration';
        str_useralpha = 'Scaler, smoothing parameter, v(k+1)=(1-alpha)*v(k)+alpha*mean(neighbors)';
        data_ = [{'Iteration number';'useralpha'} {str_iteration;str_useralpha} num2cell([2;0.5])];
        table_smoothmesh_options.Data=data_;
    end
    function turn_GUI_smoothsurfacemesh_laplacian(value)
        set(table_smoothmesh_options,'enable',value,'visible',value);
        str_iteration = 'Number of smoothing iteration';
        str_useralpha = 'Scaler, smoothing parameter, v(k+1)=(1-alpha)*v(k)+alpha*mean(neighbors)';
        data_ = [{'Iteration number';'useralpha'} {str_iteration;str_useralpha} num2cell([2;0.5])];
        table_smoothmesh_options.Data=data_;
    end
    function turn_GUI_smoothsurfacemesh_laplacianhc(value)
        set(table_smoothmesh_options,'enable',value,'visible',value);
        str_iteration = 'Number of smoothing iteration';
        str_useralpha = 'Scaler, smoothing parameter, v(k+1)=(1-alpha)*v(k)+alpha*mean(neighbors)';
        str_userbeta = 'Scaler, smoothing parameter, for laplacianhc';
        data_ = [{'Iteration number';'useralpha';'userbeta'} {str_iteration;str_useralpha;str_userbeta} num2cell([2;0.5;1])];
        table_smoothmesh_options.Data=data_;
    end

%%
%% "CREATE MESH" TAB
%%

%% PRE-ASSIGNED PARAMETERS
str_region = {'Left electrode solid','Left electrode: electrolyte','Separator: solid','Separator: electrolyte','Right electrode: solid','Right electrode: electrolyte'};

%% GUI
% Title
Text_tab_generate = uicontrol('Parent', tab_generate, 'Style', 'text','Units','normalized','Position',[description_tab_fromleft description_tab_frombottom description_tab_xlenght description_tab_ylenght],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Create mesh');

% Create mesh
Text_generate = uicontrol('Parent', tab_generate, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','1) Create mesh from the 3D array',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.88 0.8 0.03]);
line_generate = annotation(tab_generate,'line',[description_tab_fromleft 0.975],[0.88 0.88],'visible','on');
pushbutton_generatemesh = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Click to create mesh','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.325 0.77 0.3 0.1],'Callback', @generatemesh_Callback);

% Check regions
Text_checkregions = uicontrol('Parent', tab_generate, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','2) Check regions',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.71 0.8 0.03]);
line_checkregions = annotation(tab_generate,'line',[description_tab_fromleft 0.975],[0.71 0.71],'visible','on');
pushbutton_inspectregions = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Inspect regions','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.05 0.6 0.2 0.1],'Callback', @inspectregions_Callback);
pushbutton_correctregions = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Correct regions','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.275 0.6 0.2 0.1],'Callback', @correctregions_Callback);
Text_incorrectregion_number = uicontrol('Parent', tab_generate, 'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',{'- If all regions look correct after visual inspection, do not click on "Correct regions".',...
                                                                                                                                                  '- If a region has incorrect isolated cells, fulfill the information in the above table. If the region shows only incorrect isolated cells, enter the first parameter only and keep the two other empty. If the region shows both correct connected cells and incorrect isolated cells, enter all the parameters. Once done, re-click on "Inspect regions" to actualize the regions.'},...
    'visible','on','BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.05 0.43 0.925 0.16]);

table_region_incorrect = uitable('Parent', tab_generate,'enable','on','Visible','on','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.5 0.6 0.475 0.1]);
table_region_incorrect.ColumnName = {'Parameter','Value'}; % Column name
table_region_incorrect.ColumnEditable = [false true]; % Select column editable
table_region_incorrect.ColumnWidth = {'auto', 'auto'}; % Auto width
table_region_incorrect.RowName = []; % Remove row name

% Assemble electrolyte and solid domains
Text_assembledomains = uicontrol('Parent', tab_generate, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','3) Assign regions',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.39 0.8 0.03]);
line_assembledomains = annotation(tab_generate,'line',[description_tab_fromleft 0.975],[0.39 0.39],'visible','on');
Popup_region_id = uicontrol('Parent', tab_generate,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', {'Region number '},'Units','normalized','Position', [0.025 0.335 0.1 0.04],'enable','on','Visible','on');
Popup_region_assign = uicontrol('Parent', tab_generate,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', str_region,'Units','normalized','Position', [0.135 0.335 0.215 0.04],'enable','on','Visible','on');
pushbutton_validate_regionassign = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Validate assignment','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.36 0.345 0.2 0.03],'Callback', @validateregionassign_Callback);
table_region_assign = uitable('Parent', tab_generate,'enable','on','Visible','on','FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',[0.57 0.2 0.405 0.18]);
table_region_assign.ColumnName = {'Region number','Region assignment'}; % Column name
table_region_assign.ColumnEditable = [false false]; % Select column editable
table_region_assign.ColumnWidth = {'auto', 'auto'}; % Auto width
table_region_assign.RowName = []; % Remove row name
pushbutton_assembledomain = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Assign solid and electrolyte','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.025 0.26 0.3250 0.075],'Callback', @assemble_solidelectrolyte_Callback);

% Save and analyse
Text_saveanalysis = uicontrol('Parent', tab_generate, 'Style', 'text','FontSize',font_size_large_GUI,'FontName',font_name_GUI,'String','4) Analyse and save meshes',...
    'visible','on','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.1 0.15 0.8 0.03]);
line_saveanalysis = annotation(tab_generate,'line',[description_tab_fromleft 0.975],[0.15 0.15],'visible','on');
pushbutton_analysemeshes = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Analyse mesh (optional)','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.15 0.04 0.3 0.1],'Callback', @analyse_mesh_Callback);
pushbutton_savemesh = uicontrol('Parent', tab_generate,'Style', 'pushbutton','Value',0,'String','Save meshes','FontSize',font_size_large_GUI,'FontName',font_name_GUI,...
    'enable','off','BackgroundColor','w','HorizontalAlignment','center','Units','normalized','Position', [0.55 0.04 0.3 0.1],'Callback', @savemesh_Callback);

%% CALLBACKS
    function generatemesh_Callback(~,~)
        % Parameter for v2s
        parameters_ = cell2mat(table_surfamesh_options.Data(:,3));
        if strcmp(options.method_surfacemesh,'cgalmesh')
            options.radbound = parameters_(1);
            options.distbound = parameters_(2);
            options.isovalues = parameters_(3);
        end
        
        % Parameters for smoothsurf
        parameters_ = cell2mat(table_smoothmesh_options.Data(:,3));
        if strcmp(options.method_smoothsurfacemesh,'lowpass')
            options.iteration = parameters_(1);
            options.useralpha = parameters_(2);
        elseif strcmp(options.method_smoothsurfacemesh,'laplacian')
            options.iteration = parameters_(1);
            options.useralpha = parameters_(2);
        elseif strcmp(options.method_smoothsurfacemesh,'laplacianhc')
            options.iteration = parameters_(1);
            options.useralpha = parameters_(2);
            options.userbeta = parameters_(3);
        end
        
        % Parameters for surf2mesh
        parameters_ = cell2mat(table_volumemesh_options.Data(:,3));
        options.keepratio = parameters_(1);
        options.maxvol = parameters_(2);
        [node,elem,face] = mesh_generation_withIso2mesh(cell_array, options); % Create volumetric mesh
        
        % Region numerotation and sorting
        [region_] = mesh_regionnumerotation(node,elem,face);        
        number_region_vol = length(region_);
      
        % Update GUI
        Popup_region_id.String = [{'Region number '}; num2cell([1:1:number_region_vol]')]; % Update GUI string
        set(Popup_region_assign,'String', ['Cell assignment' str_region]);
        set(pushbutton_validate_regionassign,'enable','on');
        set(pushbutton_inspectregions,'enable','on');
        set(pushbutton_assembledomain,'enable','on');
        set(table_region_incorrect,'enable','on');
        data_ = [{'Incorrect region number';'Incorrect isolated element for x >';'Incorrect isolated element for x <'} {' ';' ';' '}];
        table_region_incorrect.Data = data_;
        table_region_assign.Data = []; % Reset
    end

    function inspectregions_Callback(~,~)
        % Region numerotation and sorting
        [region_] = mesh_regionnumerotation(node,elem,face);      
        number_region_vol = length(region_);
        Popup_region_id.String = [{'Region number '}; num2cell([1:1:number_region_vol]')]; % Update GUI string
        mesh_inspectregions(region_,options)
        set(pushbutton_correctregions,'enable','on');
    end

    function correctregions_Callback(~,~)
        wrong_region = str2num(cell2mat(table_region_incorrect.Data(1,2)));
        region_numerotation=unique(elem(:,5));
        wrong_region = region_numerotation(wrong_region);
        threshold_right = str2num(cell2mat(table_region_incorrect.Data(2,2)));
        threshold_left = str2num(cell2mat(table_region_incorrect.Data(3,2)));
        if ~isempty(threshold_right) && ~isempty(threshold_left)
            [elem] = function_select_isolatedcell(wrong_region, elem, node, threshold_left, threshold_right);
            wrong_region = 9999; % Overwritte
        end
        [elem] = function_assign_isolatedcell(elem,wrong_region);
    end

    function validateregionassign_Callback(~,~)
        region_id_string = cell2mat(Popup_region_id.String(2:end));
        region_id_value = Popup_region_id.Value-1;
        region_assignment_string = Popup_region_assign.String(2:end);
        region_assignment_value = Popup_region_assign.Value-1;
        if region_id_value>=1 && region_assignment_value>=1
            val_ = table_region_assign.Data;
            if isempty(val_)
                val_ = [num2cell(region_id_string(region_id_value)) region_assignment_string(region_assignment_value)];
            else
                if sum(unique(str2num(cell2mat(val_(:,1)))) == str2num(region_id_string(region_id_value)))== 0 % Not yet assigned
                    new_line = [num2cell(region_id_string(region_id_value)) region_assignment_string(region_assignment_value)];
                    val_=[val_ ; new_line];
                end
            end
            table_region_assign.Data= val_; % Update table
            if length(table_region_assign.Data) == length(region_id_string)
                set(pushbutton_savemesh,'enable','on');
                set(pushbutton_analysemeshes,'enable','on');
            else
                set(pushbutton_savemesh,'enable','off');
                set(pushbutton_analysemeshes,'enable','off');
            end
        end
    end

    function assemble_solidelectrolyte_Callback(~,~)
        region_number = table_region_assign.Data(:,1);
        region_name = table_region_assign.Data(:,2);
                
        % Left electrode solid
        idx_codes = find(strcmp(region_name,'Left electrode: solid'));
        left_solid = str2num(cell2mat(region_number(idx_codes)));
        
        % Right electrode solid
        idx_codes = find(strcmp(region_name,'Right electrode: solid'));
        right_solid = str2num(cell2mat(region_number(idx_codes)));
        
        % Electrolyte
        idx_codes = find(strcmp(region_name,'Separator: electrolyte'));
        electrolyte = str2num(cell2mat(region_number(idx_codes)));
        idx_codes = find(strcmp(region_name,'Left electrode: electrolyte'));
        if ~isempty(idx_codes)
            electrolyte = [electrolyte str2num(cell2mat(region_number(idx_codes)))];
        end
        idx_codes = find(strcmp(region_name,'Right electrode: electrolyte'));
        if ~isempty(idx_codes)
            electrolyte = [electrolyte str2num(cell2mat(region_number(idx_codes)))];
        end
        
        region_numerotation=unique(elem(:,5));
        left_solid = region_numerotation(left_solid);
        right_solid = region_numerotation(right_solid);
        electrolyte = region_numerotation(electrolyte);

        % Left electrode
        [node_leftelectrode,elem_leftelectrode,face_leftelectrode,corresponding_index_leftelectrode,elem_leftelectrode_commonindex,face_leftelectrode_commonindex] = Function_NodeFaceElem_region(node,elem,face,left_solid);
        % Electrolyte: from separator and both electrodes
        [node_electrolyte,elem_electrolyte,face_electrolyte,corresponding_index_electrolyte,elem_electrolyte_commonindex,face_electrolyte_commonindex] = Function_NodeFaceElem_region(node,elem,face,electrolyte);
        % Right electrode: solid phase only
        [node_rightelectrode,elem_rightelectrode,face_rightelectrode,corresponding_index_rightelectrode,elem_rightelectrode_commonindex,face_rightelectrode_commonindex] = Function_NodeFaceElem_region(node,elem,face,right_solid);

        % Show domains
        if options.save.dommain_left
            function_show_mesh(node_leftelectrode,elem_leftelectrode,face_leftelectrode,'Left electrode','Left_electrode',options);
        end
        if options.save.dommain_electrolyte
            function_show_mesh(node_electrolyte,elem_electrolyte,face_electrolyte,'Electrolyte','Electrolyte',options);
        end
        if options.save.dommain_right
            function_show_mesh(node_rightelectrode,elem_rightelectrode,face_rightelectrode,'Right electrode','Right_electrode',options);
        end  

        % Node at the interfaces
        if options.save.nodeinterface_fig || options.save.nodeinterface_png || options.save.nodeinterface_mat
            function_get_commonnode_interface(node_electrolyte,node_leftelectrode,face_leftelectrode, node_rightelectrode,face_rightelectrode,options)
        end
        
        %% FROM JEFFEREY ALLEN: some hard-coded value to fix
        % % All facets (external and internal)
        % calculate facets for entire mesh
        facets = Function_calculate_facets(elem(:,1:4));
        % calculate facets for the sub meshes
        facets_leftelectrode = Function_calculate_facets(elem_leftelectrode(:,1:4));
        facets_electrolyte = Function_calculate_facets(elem_electrolyte(:,1:4));
        facets_rightelectorde = Function_calculate_facets(elem_rightelectrode(:,1:4));
        % calculate facets fot the sub meshes using the index for the entire mesh
        facets_leftelectrode_commonindex = Function_calculate_facets(elem_leftelectrode_commonindex(:,1:4));
        facets_electrolyte_commonindex = Function_calculate_facets(elem_electrolyte_commonindex(:,1:4));
        facets_rightelectorde_commonindex = Function_calculate_facets(elem_rightelectrode_commonindex(:,1:4));        
        
        % % Find interface facets
        % Find LeftAM and Electrolyte interface
        [~,leftAM_electrolyte_facets_index,~] = intersect(facets,intersect(facets_electrolyte_commonindex,facets_leftelectrode_commonindex,'rows'),'rows');
        % Find RightAM and Electrolyte interface
        [~,rightAM_electrolyte_facets_index,~] = intersect(facets,intersect(facets_electrolyte_commonindex,facets_rightelectorde_commonindex,'rows'),'rows');
        % Find Interior Facets for domains
        [~,LeftAM_facets_index,~] = intersect(facets, facets_leftelectrode_commonindex,'rows');
        [~,Electrolyte_facets_index,~] = intersect(facets, facets_electrolyte_commonindex,'rows');
        [~,RightAM_facets_index,~] = intersect(facets, facets_rightelectorde_commonindex,'rows');

        % % Find boundary facets
        % define the search tol
        bound_tol = 1.0;
        % Sort the faces columnwise
        face_leftelectrode_commonindex_sorted = sort(face_leftelectrode_commonindex(:,1:3),2);
        face_rightelectrode_commonindex_sorted = sort(face_rightelectrode_commonindex(:,1:3),2);
        % find the x value associated with each facets nodes
        left_boundary_x_values = node(face_leftelectrode_commonindex_sorted);
        right_boundary_x_values = node(face_rightelectrode_commonindex_sorted);
        % define the range of x values
        min_values = min(node);
        max_values = max(node);
        % find the boundary
        [left_boundary_rows,~] = find(left_boundary_x_values <= min_values(1) + bound_tol );
        [right_boundary_rows,~] = find(right_boundary_x_values >= max_values(1) - bound_tol );
        % filter the unique values
        left_boundary_rows = unique(left_boundary_rows);
        right_boundary_rows = unique(right_boundary_rows);
        % find where the boundary facets in the entire facets list
        [~,left_boundary_facets_index,~] = intersect(facets,face_leftelectrode_commonindex_sorted(left_boundary_rows,1:3),'rows');
        [~,right_boundary_facets_index,~] = intersect(facets,face_rightelectrode_commonindex_sorted(right_boundary_rows,1:3),'rows');
        % remove the facets that are apart of the interface
        [~,right_boundary_interface_index,~] = intersect(right_boundary_facets_index,rightAM_electrolyte_facets_index);
        right_boundary_facets_index(right_boundary_interface_index) = [];
        % define the boundary markers
        boundary_markers = zeros(length(facets),1);
        boundary_markers(left_boundary_facets_index)=1; % HARD CODED VALUE
        boundary_markers(right_boundary_facets_index)=4; % HARD CODED VALUE
        
        Fig_ = figure;
        Fig_.Name= 'Common facets';
        Fig_.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig_,'position',scrsz); % Full screen figure
        % - Create axes
        axes_ = axes('Parent',Fig_);
        hold(axes_,'on');
        % All facets
%         trisurf(face_leftelectrode_commonindex(:,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','blue','FaceAlpha',0.1) 
%         trisurf(face_rightelectrode_commonindex(:,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','blue','FaceAlpha',0.1)
        % Boundary Facets
        trisurf(facets(left_boundary_facets_index,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','green','FaceAlpha',0.5)
        trisurf(facets(right_boundary_facets_index,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','red','FaceAlpha',0.5)
        % Interface Facets
        trisurf(facets(leftAM_electrolyte_facets_index,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','magenta','FaceAlpha',0.5)
        trisurf(facets(rightAM_electrolyte_facets_index,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','cyan','FaceAlpha',0.5)
        % trisurf(face_leftelectrode_commonindex(left_boundary_rows,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','green','FaceAlpha',0.5)
        % trisurf(face_rightelectrode_commonindex(right_boundary_rows,1:3),node(:,1),node(:,2),node(:,3),'Facecolor','red','FaceAlpha',0.5)
        axis equal; axis tight; view(3);
        filename = 'Common_facets';
        hold(axes_,'off');
        fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
        % .fig
        savefig(Fig_,[fullpath filename])
        % .png
        saveas(Fig_,[fullpath filename],'png')
        
        % % Create the facet function containing boundaries, interfaces, internal domains
        % Create facets marker
        facet_markers = zeros(length(facets),1); % HARD CODED VALUE
        facet_markers(LeftAM_facets_index)=1; % HARD CODED VALUE
        facet_markers(Electrolyte_facets_index)=2; % HARD CODED VALUE
        facet_markers(RightAM_facets_index)=3; % HARD CODED VALUE
        facet_markers(left_boundary_facets_index)=4; % HARD CODED VALUE
        facet_markers(leftAM_electrolyte_facets_index)=5; % HARD CODED VALUE
        facet_markers(rightAM_electrolyte_facets_index)=6; % HARD CODED VALUE
        facet_markers(right_boundary_facets_index)=7; % HARD CODED VALUE
        
        % % Save Boundary Markers
        fullpath=[options.save.mainfolder options.save.meshdatafolder];
        save([fullpath,'Facets.mat'],'facet_markers')
        set(pushbutton_savemesh,'enable','on');
    end

    function analyse_mesh_Callback(~,~)
        fullpath=[options.save.mainfolder options.save.analysemeshfolder];
        if ~exist(fullpath,'dir')
            mkdir(fullpath);
        end
        % % Euler characteristic
        % The surface of a finite, convex, 3-dimensional object is topologically equivalent to the surface of a nice, smooth, round 3-dimensional ball,
        % and X=V-E+F is a topological invariant = 2
        % It holds true for any mesh that can be smoothly deformed to a sphere.
        % Thus, the number of holes G in the mesh controlled the Euler number H=2*(1-X)
        [Euler_characteristic_left,Number_of_verteces_left,Number_of_edges_left,Number_of_faces_left]=mesheuler(face_leftelectrode);
        Number_of_holes_left=1-(Euler_characteristic_left/2);
        [Euler_characteristic_electrolyte,Number_of_verteces_electrolyte,Number_of_edges_electrolyte,Number_of_faces_electrolyte]=mesheuler(face_electrolyte);
        Number_of_holes_electrolyte=1-(Euler_characteristic_electrolyte/2);
        [Euler_characteristic_right,Number_of_verteces_right,Number_of_edges_right,Number_of_faces_right]=mesheuler(face_rightelectrode);
        Number_of_holes_right=1-(Euler_characteristic_right/2);
        
        Table_Euler = table({'Number of verteces V';'Number of edges E';'Number of faces F';'Euler characteristic X=V-E+F';'Number of holes 1-X/2'},...
            [Number_of_verteces_left;Number_of_edges_left;Number_of_faces_left;Euler_characteristic_left;Number_of_holes_left],...
            [Number_of_verteces_electrolyte;Number_of_edges_electrolyte;Number_of_faces_electrolyte;Euler_characteristic_electrolyte;Number_of_holes_electrolyte],...
            [Number_of_verteces_right;Number_of_edges_right;Number_of_faces_right;Euler_characteristic_right;Number_of_holes_right],...
            'VariableNames',{'Property' 'Left_solid' 'Electrolyte' 'Right_solid'});%
        
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Euler_characteristics';
        DATA_writetable.sheet(1).table=Table_Euler;
        % Save function
        Function_Writetable(fullpath,'Euler_characteristics',DATA_writetable)
        if options.analysemesh.left
            tolerance_boundary = [1 1 1;-1 1 1]'*1;
            function_analyse_mesh(node_leftelectrode, elem_leftelectrode, face_leftelectrode, voxel_size_um, tolerance_boundary, celldata(4).left_binary_array, 'left electrode (solid)','left_electrode', options.domainfolder_left, options);
        end
        if options.analysemesh.electrolyte
            tolerance_boundary = [1 1 1;1 1 1]'*1;
            function_analyse_mesh(node_electrolyte, elem_electrolyte, face_electrolyte, voxel_size_um, tolerance_boundary, celldata(4).electrolyte_binary_array, 'electrolyte','electrolyte', options.domainfolder_electrolyte, options)
        end
        if options.analysemesh.right
            tolerance_boundary = [-1 1 1;1 1 1]'*1;
            function_analyse_mesh(node_rightelectrode, elem_rightelectrode, face_rightelectrode, voxel_size_um, tolerance_boundary, celldata(4).right_binary_array, 'right electrode (solid)','right_electrode', options.domainfolder_right, options)
        end        
    end

    function savemesh_Callback(~,~)
        % Save each domain
        if options.save.microstructurestl==1 || options.save.microstructurebinarystl==1 || options.save.microstructuremesh==1
            function_save_mesh(node_leftelectrode,elem_leftelectrode,'Left_electrode',options);
            function_save_mesh(node_electrolyte,elem_electrolyte,'Electrolyte',options);
            function_save_mesh(node_rightelectrode,elem_rightelectrode,'Right_electrode',options);
        end
        % Save all domain in one mesh
        if options.save.cellmesh
            fullpath=[options.save.mainfolder options.save.meshdatafolder];
            if ~exist(fullpath,'dir')
                mkdir(fullpath);
            end
            elem_=elem(:,1:4); % Remove useless last column
            elem_=elem_-1; % Start at 0 for python index convention
            save([fullpath,'Nodes.mat'],'node')
            save([fullpath,'Tetrahedron.mat'],'elem_')
            
            sub_markers=elem(:,5);
            %make sure to reorder the regions. 1:left_active, 2:electrolyte/seperator, % 3:right_active
            %sub_markers(elem(:,5)==region_numerotation(4))=1;
            %sub_markers(elem(:,5)==region_numerotation(3))=2;
            %sub_markers(elem(:,5)==region_numerotation(1))=2;
            %sub_markers(elem(:,5)==region_numerotation(2))=3;
            save([fullpath,'Subdomains.mat'],'sub_markers')
        end
    end

%% FUNCTIONS
    function [elem] = function_select_isolatedcell(wrong_region, elem, node, threshold_region_left, threshold_region_right)
        new_region = 9999; % This region contain all the isolated irrelevant meshed and is corrected later in the code
        % Loop over all elements
        for current_elem=1:1:length(elem)
            % Select region
            current_region=elem(current_elem,5);
            if current_region==wrong_region
                % select node
                node1=elem(current_elem,1);
                node2=elem(current_elem,2);
                node3=elem(current_elem,3);
                node4=elem(current_elem,4);
                % Select coordinates
                x1=node(node1,1);
                x2=node(node2,1);
                x3=node(node3,1);
                x4=node(node4,1);
                % Change region
                if ~isempty(threshold_region_right)
                    if x1>threshold_region_right || x2>threshold_region_right || x3>threshold_region_right || x4>threshold_region_right
                        elem(current_elem,5)=new_region;
                    end
                end
                if ~isempty(threshold_region_left)
                    if x1<threshold_region_left || x2<threshold_region_left || x3<threshold_region_left || x4<threshold_region_left
                        elem(current_elem,5)=new_region;
                    end
                end
            end
        end
    end

    function [elem] = function_assign_isolatedcell(elem,wrong_region)
        if wrong_region==9999 % hard-coded value
            tmp = unique(elem(:,5));
            tmp(tmp==9999)=[];
            other_region_max = max(tmp);
            other_region_min = min(tmp);
        else
            other_region_max = max(elem(:,5));
            other_region_min = min(elem(:,5));
        end
        for pass_=1:1:2
            incorrect_cell=1e9; % Initialize
            while incorrect_cell>0
                for other_region=other_region_max:-1:other_region_min
                    if other_region~=wrong_region
                        change_=1e9; % Initialization
                        previous_incorrect_cell = incorrect_cell;
                        while change_~=0
                            
                            % All node that belong to the wrong zone
                            idx_wrongregion=find(elem(:,5)==wrong_region);
                            all_node_wrongregion=[elem(idx_wrongregion,1); elem(idx_wrongregion,2); elem(idx_wrongregion,3); elem(idx_wrongregion,4)];
                            all_node_wrongregion=unique(all_node_wrongregion);
                            
                            % Check if they belong to other zone
                            idx_otherregion=find(elem(:,5)==other_region);
                            all_node_otherregion=[elem(idx_otherregion,1); elem(idx_otherregion,2); elem(idx_otherregion,3); elem(idx_otherregion,4)];
                            all_node_otherregion=unique(all_node_otherregion);
                            
                            % Node of wrong region that belong also to current other region
                            idx_tmp = ismember(all_node_wrongregion,all_node_otherregion);
                            idx_tmp = find(idx_tmp==1);
                            index_node_WrongInOther = all_node_wrongregion(idx_tmp);
                            
                            index_node_WrongInOther_a = ismember(elem(:,1),index_node_WrongInOther);
                            index_node_WrongInOther_b = ismember(elem(:,2),index_node_WrongInOther);
                            index_node_WrongInOther_c = ismember(elem(:,3),index_node_WrongInOther);
                            index_node_WrongInOther_d = ismember(elem(:,4),index_node_WrongInOther);
                            idx_a = find(index_node_WrongInOther_a==1);
                            idx_b = find(index_node_WrongInOther_b==1);
                            idx_c = find(index_node_WrongInOther_c==1);
                            idx_d = find(index_node_WrongInOther_d==1);
                            idx_all=[idx_a; idx_b; idx_c; idx_d];
                            idx_all=unique(idx_all);
                            
                            tmp = elem(idx_all,5);
                            tmp2 = find(tmp(:,1)==wrong_region);
                            idx_all_ = idx_all(tmp2);
                            
                            incorrect_cell=length(idx_all_);
                            fprintf ('- Pass %i. Check incorrect region with region %i, number of incorrect cell %i\n',pass_,other_region,incorrect_cell);
                            
                            change_ = incorrect_cell - previous_incorrect_cell;
                            previous_incorrect_cell = incorrect_cell;
                            
                            elem(idx_all_,5)=other_region;
                        end
                    end
                end
            end
        end
    end

end