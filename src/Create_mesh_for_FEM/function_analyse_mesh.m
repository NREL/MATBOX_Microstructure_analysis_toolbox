function [] = function_analyse_mesh(node_region, elem_region, face_region, voxel_size_um, tolerance_boundary_activeinterface, binary_array, str_region_name, str_region_name_sav, domainfolder, options)

% Save path
fullpath=[options.save.mainfolder options.save.analysemeshfolder domainfolder];
if ~exist(fullpath,'dir')
    mkdir(fullpath);
end

% Cumulative and distribution function shared parameters
parameters_distributionfigure.figureposition = [100 100 1500 800];
parameters_distributionfigure.fontname = options.fontname;
parameters_distributionfigure.grid = options.grid;
parameters_distributionfigure.minorgrid = options.minorgrid;
parameters_distributionfigure.fullpath = fullpath;
% Smoothing function (moving average) parameters (inherit from options)
parameters_smoothingmovingaverage.smooth_cumulative_fct = options.smooth_cumulative_fct;
parameters_smoothingmovingaverage.minimum_array_length = options.minimum_array_length;
parameters_smoothingmovingaverage.number_point = options.number_point;
parameters_smoothingmovingaverage.moving_rangeratio = options.moving_rangeratio;
parameters_smoothingmovingaverage.moving_range = options.moving_range;
parameters_smoothingmovingaverage.enforce_samelength = options.enforce_samelength;
parameters_smoothingmovingaverage.origin = options.origin;
parameters_smoothingmovingaverage.boundary_behavior = options.boundary_behavior;
% 3D patch figure shared parameters
parameters_3Dpatchfigure.fullpath = fullpath;
parameters_3Dpatchfigure.figureposition = [200 100 1000 800];
parameters_3Dpatchfigure.fontname = options.fontname;
parameters_3Dpatchfigure.xlabel = 'Axe 1 (voxels length)';
parameters_3Dpatchfigure.ylabel = 'Axe 2 (voxels length)';
parameters_3Dpatchfigure.zlabel = 'Axe 3 (voxels length)';
parameters_3Dpatchfigure.grid = options.grid;
parameters_3Dpatchfigure.minorgrid = options.minorgrid;
parameters_3Dpatchfigure.colormap = 'jet';
parameters_3Dpatchfigure.lighting = 'flat';
data_3Dpatchfigure.FaceColor = 'flat';
data_3Dpatchfigure.EdgeColor = 'none';
data_3Dpatchfigure.face = face_region(:,1:3);
data_3Dpatchfigure.vertices = node_region;

%%
%% CALCULATE PROPERTIES
%%

%% VOLUME AND SURFACE

% Cell volume
[cell_volume, cell_centroid] = Function_CalculateCellVolume(node_region,elem_region(:,1:4)); % Mesh volume information
volume_allmesh = sum(cell_volume); % Volume of the phase
volume_allmesh = volume_allmesh * (voxel_size_um^3); % Convert in um3

% Mesh surface information
[Angle_faces, ~, Area_faces,~,~] = Function_CalculateFacetsAngle(node_region,face_region(:,1:3),elem_region(:,1:4), cell_centroid);
% Beta in development: facet normal
% [Angle_faces, ~, Area_faces,Normal_faces,~] = Function_CalculateFacetsAngle_beta(node_region,face_region(:,1:3),elem_region(:,1:4), cell_centroid);

% Minimum and maximum angles
minimum_angle_faces = min(Angle_faces,[],2);
maximum_angle_faces = max(Angle_faces,[],2);
% Surface
surface_allmesh = sum(Area_faces);
surface_allmesh = surface_allmesh * (voxel_size_um^2); % Convert in um3

% Angle between adjacent facet normal
% [verteces_belong_to_face,Normal_adjacentfaces_maxdelta] = Funcion_CalculateFacetNormalVariation(node_region,face_region(:,1:3),Normal_faces);


%% SPECIFIC SURFACE AREA
% Specific surface area SSA: domain surface / domain volume
% Effective specific surface area ESSA: domain surface / effective volume

% With whole mesh
SSA_allmesh = surface_allmesh/volume_allmesh;
effectivevolume_allmesh = function_effectivevolume(node_region);
effectivevolume_allmesh = effectivevolume_allmesh * (voxel_size_um^3); % Convert in um3
selected_surface_withboundaries = zeros(length(face_region),1) + 1;
ESSA_allmesh = surface_allmesh/effectivevolume_allmesh;
str_wholemesh_a = ['Surface=' num2str(surface_allmesh,'%1.1f') '\mum^{2}' ', Volume=' num2str(volume_allmesh,'%1.1f') '\mum^{3}' ', Effective volume=' num2str(effectivevolume_allmesh,'%1.1f') '\mum^{3}'];
str_wholemesh_b = ['Specific surface area=' num2str(SSA_allmesh,'%1.3f') '\mum^{-1}' ', Effective specific surface area=' num2str(ESSA_allmesh,'%1.3f') '\mum^{-1}'];

% With mesh subregion (active interface)
[node_subregion,elem_subregion,face_subregion, index_subregion_node, index_subregion_elem, index_subregion_face] = Function_subregion(node_region,elem_region,face_region(:,1:3),tolerance_boundary_activeinterface);
selected_surface_activeinterface = zeros(length(face_region),1);
selected_surface_activeinterface(index_subregion_face)=1;
surfacesubregionactiveinterface = sum(Area_faces(index_subregion_face));
volumesubregionactiveinterface = sum(cell_volume(index_subregion_elem));
surfacesubregionactiveinterface = surfacesubregionactiveinterface * (voxel_size_um^2); % Convert in um2
volumesubregionactiveinterface = volumesubregionactiveinterface * (voxel_size_um^3); % Convert in um3
SSA_subregionactiveinterface=surfacesubregionactiveinterface/volumesubregionactiveinterface; % (domain surface / domain volume)
effectivevolumesubregionactiveinterface = function_effectivevolume(node_subregion); % Effective volume
effectivevolumesubregionactiveinterface = effectivevolumesubregionactiveinterface * (voxel_size_um^3); % Convert in um3
ESSA_subregionactiveinterface = surfacesubregionactiveinterface/effectivevolumesubregionactiveinterface;  % (domain surface / effective volume)
str_activeinterface_a = ['Surface=' num2str(surfacesubregionactiveinterface,'%1.1f') '\mum^{2}' ', Volume=' num2str(volumesubregionactiveinterface,'%1.1f') '\mum^{3}' ', Effective volume=' num2str(effectivevolumesubregionactiveinterface,'%1.1f') '\mum^{3}'];
str_activeinterface_b = ['Specific surface area=' num2str(SSA_subregionactiveinterface,'%1.3f') '\mum^{-1}' ', Effective specific surface area=' num2str(ESSA_subregionactiveinterface,'%1.3f') '\mum^{-1}'];

% With mesh subregion (without boundaries)
tolerance_boundary = [1 1 1;1 1 1]'*1;
[node_subregion,elem_subregion,face_subregion, index_subregion_node, index_subregion_elem, index_subregion_face] = Function_subregion(node_region,elem_region,face_region(:,1:3),tolerance_boundary);
selected_surface_withoutboundaries = zeros(length(face_region),1);
selected_surface_withoutboundaries(index_subregion_face)=1;
surfacesubregionwithoutboundaries = sum(Area_faces(index_subregion_face));
volumesubregionwithoutboundaries = sum(cell_volume(index_subregion_elem));
surfacesubregionwithoutboundaries = surfacesubregionwithoutboundaries * (voxel_size_um^2); % Convert in um2
volumesubregionwithoutboundaries = volumesubregionwithoutboundaries * (voxel_size_um^3); % Convert in um3
SSA_subregionwithoutboundaries=surfacesubregionwithoutboundaries/volumesubregionwithoutboundaries; % (domain surface / domain volume)
effectivevolumesubregionwithoutboundaries = function_effectivevolume(node_subregion); % Effective volume
effectivevolumesubregionwithoutboundaries = effectivevolumesubregionwithoutboundaries * (voxel_size_um^3); % Convert in um3
ESSA_subregionwithoutboundaries = surfacesubregionwithoutboundaries/effectivevolumesubregionwithoutboundaries;  % (domain surface / effective volume)
str_withoutboundary_a = ['Surface=' num2str(surfacesubregionwithoutboundaries,'%1.1f') '\mum^{2}' ', Volume=' num2str(volumesubregionwithoutboundaries,'%1.1f') '\mum^{3}' ', Effective volume=' num2str(effectivevolumesubregionwithoutboundaries,'%1.1f') '\mum^{3}'];
str_withoutboundary_b = ['Specific surface area=' num2str(SSA_subregionwithoutboundaries,'%1.3f') '\mum^{-1}' ', Effective specific surface area=' num2str(ESSA_subregionwithoutboundaries,'%1.3f') '\mum^{-1}'];

% With 3D array (to be compared with the mesh subregion without boundaries)
[~,surface_cube] = Function_Specificsurface_direct_Algorithm(binary_array);
index_=find(binary_array==1);
volume_cube = length(index_);
surface_cube = surface_cube * (voxel_size_um^2); % Convert in um2
volume_cube = volume_cube * (voxel_size_um^3); % Convert in um3
SSA_cube = surface_cube/volume_cube; % (domain surface / domain volume)
[ix,iy,iz]=ind2sub(size(binary_array),index_); % Convert to matrix index
x_min=min(ix)-1; y_min=min(iy)-1; z_min=min(iz)-1; % Take min and max
x_max=max(ix); y_max=max(iy); z_max=max(iz);
effectivevolumecube = (x_max-x_min) * (y_max-y_min) * (z_max-z_min); % Effective volume
effectivevolumecube = effectivevolumecube * (voxel_size_um^3); % Convert in um3
ESSA_cube = surface_cube/effectivevolumecube;  % (domain surface / effective volume)

% Comparison
surface_ratio = surface_cube/surfacesubregionwithoutboundaries;
volume_ratio = volume_cube/volumesubregionwithoutboundaries;
SSA_ratio = SSA_cube/SSA_subregionwithoutboundaries;
effectivevolume_ratio = effectivevolumecube/effectivevolumesubregionwithoutboundaries;
ESSA_ratio = ESSA_cube/ESSA_subregionwithoutboundaries;

Table_specificsurfacearea = table({'Domain surface S';'Domain volume V';'Specific surface area S/V';'Effective volume Veff';'Effective specific surface area S/Veff'},...
    {'um^{2}';'um^{3}';'um^{-1}';'um^{3}';'um^{-1}'},...
    [surface_allmesh;volume_allmesh;SSA_allmesh;effectivevolume_allmesh;ESSA_allmesh],...
    [surfacesubregionactiveinterface;volumesubregionactiveinterface;SSA_subregionactiveinterface;effectivevolumesubregionactiveinterface;ESSA_subregionactiveinterface],...
    [surfacesubregionwithoutboundaries;volumesubregionwithoutboundaries;SSA_subregionwithoutboundaries;effectivevolumesubregionwithoutboundaries;ESSA_subregionwithoutboundaries],...
    [surface_cube;volume_cube;SSA_cube;effectivevolumecube;ESSA_cube],...
    [surface_ratio;volume_ratio;SSA_ratio;effectivevolume_ratio;ESSA_ratio],...
    'VariableNames',{'Property' 'Unit' 'Full_mesh' 'Active_interface' 'Mesh_wo_boundaries' 'From_3D_array' 'Ratio_adimentionsal'});%

% Prepare the data
clear DATA_writetable
DATA_writetable.sheet(1).name='Specific_surface_area';
DATA_writetable.sheet(1).table=Table_specificsurfacearea;
% Save function
Function_Writetable(fullpath,'Specific_surface_area',DATA_writetable)


%% MESH QUALITY
%  compute the Joe-Liu mesh quality measure of an N-D mesh (N<=3)
%  input:
%     node:  node coordinates of the mesh (nn x 3)
%     elem:  element table of an N-D mesh (ne x (N+1))
%  output:
%     quality: a vector of the same length as size(elem,1), with
%            each element being the Joe-Liu mesh quality metric (0-1) of
%            the corresponding element. A value close to 1 represents
%            higher mesh quality (1 means equilateral tetrahedron);
%            a value close to 0 means nearly degenerated element.
%  reference:
%     A. Liu, B. Joe, Relationship between tetrahedron shape measures,
%                     BIT 34 (2) (1994) 268-287.
mesh_quality=meshquality(node_region,elem_region);

%%
%% FIGURES
%%

% 3D figure: minimum facet angles
parameters_3Dpatchfigure.figurename = ['Mimimum facets angle, ' str_region_name];
parameters_3Dpatchfigure.title = {'Minimum angle of the surface facets',str_region_name};
parameters_3Dpatchfigure.filename = ['Minimum_angle_surfacefacets_3D_' str_region_name_sav];
data_3Dpatchfigure.FaceVertexCData = minimum_angle_faces;
function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)

% 3D figure: facet area
parameters_3Dpatchfigure.figurename = ['Facets area, ' str_region_name];
parameters_3Dpatchfigure.title = {'Facets area',str_region_name};
parameters_3Dpatchfigure.filename = ['Area_surfacefacets_3D_' str_region_name_sav];
data_3Dpatchfigure.FaceVertexCData = Area_faces;
function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)

% % 3D figure: facet normal, direction 1
% parameters_3Dpatchfigure.figurename = ['Facets normal, ' str_region_name];
% parameters_3Dpatchfigure.title = {'Facets normal (direction 1)',str_region_name};
% parameters_3Dpatchfigure.filename = [ 'Nornal_surfacefacets_3D_dir1_' str_region_name_sav];
% data_3Dpatchfigure.FaceVertexCData = Normal_faces(:,1);
% function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)
% 
% % 3D figure: facet normal, direction 2
% parameters_3Dpatchfigure.figurename = ['Facets normal, ' str_region_name];
% parameters_3Dpatchfigure.title = {'Facets normal (direction 2)',str_region_name};
% parameters_3Dpatchfigure.filename = [ 'Nornal_surfacefacets_3D_dir2_' str_region_name_sav];
% data_3Dpatchfigure.FaceVertexCData = Normal_faces(:,2);
% function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)
% 
% % 3D figure: facet normal, direction 3
% parameters_3Dpatchfigure.figurename = ['Facets normal, ' str_region_name];
% parameters_3Dpatchfigure.title = {'Facets normal (direction 3)',str_region_name};
% parameters_3Dpatchfigure.filename = [ 'Nornal_surfacefacets_3D_dir3_' str_region_name_sav];
% data_3Dpatchfigure.FaceVertexCData = Normal_faces(:,3);
% function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)
% 
% % 3D figure: facet normal variation
% parameters_3Dpatchfigure.figurename = ['Facet normal variation, ' str_region_name];
% parameters_3Dpatchfigure.title = {'Facet normal variation',str_region_name};
% parameters_3Dpatchfigure.filename = ['Facet_normal_variation_3D_' str_region_name_sav];
% data_3Dpatchfigure.FaceVertexCData = Normal_adjacentfaces_maxdelta;
% function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)

% % Calculate cumulative and distribution function
% parameters_smoothingmovingaverage.round_value = 1;
% [results, ~] = Function_probability_density(Normal_adjacentfaces_maxdelta,[],parameters_smoothingmovingaverage);
% parameters_distributionfigure.figurename =  ['Normal variation, ' str_region_name];
% parameters_distributionfigure.filename = ['Facetnormalvar_distribution_' str_region_name_sav];
% parameters_distributionfigure.subaxe1_title = ['Facet normal variation: cumulative fct., ' str_region_name];
% parameters_distributionfigure.subaxe2_title = ['Facet normal variation: distribution fct., ' str_region_name];
% parameters_distributionfigure.xlabel = 'Facet normal variation (deg)';
% data_distributionfigure.cumulative_fct = results.cumulative_fct;
% data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
% data_distributionfigure.probability_density_fct = results.probability_density_fct;
% data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
% data_distributionfigure.x50 = results.x50;
% data_distributionfigure.smoothed_x50 = results.smoothed_x50;
% data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
% data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
% data_distributionfigure.unit = 'deg'; 
% function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);

Fig_= figure; % Create figure
Fig_.Name= ['Triangle angles at the mesh surface, ' str_region_name];
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', [460 428 660 520]);
axes_ = axes('Parent',Fig_); % Create axes
hold(axes_,'on');
t_=title (' ','FontName',options.fontname,'FontSize',16); % - Title
t_.String= {'Triangle angles at the mesh surface',str_region_name};
h_=plot(minimum_angle_faces,maximum_angle_faces); % Curves
set(h_,'LineStyle','none','Marker','x'); % Colors
xlabel('Minimum angle (deg)'); % - Axis label
ylabel('Maximum angle (deg)');
grid(axes_,options.grid); % Display grid
set(axes_,'XMinorGrid',options.minorgrid,'YMinorGrid',options.minorgrid); % Display grid for minor thicks also
set(axes_,'FontName',options.fontname,'FontSize',14); % - Fontname and fontsize
hold(axes_,'off');
% Save figures
filename = ['Minimum_and_maximum_angle_' str_region_name_sav];
savefig(Fig_,[fullpath filename]) % .fig
saveas(Fig_,[fullpath filename],'png') % .png
close(Fig_)

% Calculate cumulative and distribution function
parameters_smoothingmovingaverage.round_value = 1;
[results, ~] = Function_probability_density(minimum_angle_faces,[],parameters_smoothingmovingaverage);
parameters_distributionfigure.figurename =  ['Minimum triangle angles at the mesh surface, ' str_region_name];
parameters_distributionfigure.filename = ['Surfacetriangle_Minimumangles_distribution_' str_region_name_sav];
parameters_distributionfigure.subaxe1_title = ['Minimum surface triangle angles: cumulative fct., ' str_region_name];
parameters_distributionfigure.subaxe2_title = ['Minimum surface triangle angles: distribution fct., ' str_region_name];
parameters_distributionfigure.xlabel = 'Minimum surface angle (deg)';
data_distributionfigure.cumulative_fct = results.cumulative_fct;
data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
data_distributionfigure.probability_density_fct = results.probability_density_fct;
data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
data_distributionfigure.x50 = results.x50;
data_distributionfigure.smoothed_x50 = results.smoothed_x50;
data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
data_distributionfigure.unit = 'deg'; 
function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);

Fig_= figure; % Create figure
Fig_.Name= ['Minimum triangle angles and areas at the mesh surface, ' str_region_name];
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', [460 428 660 520]);
axes_ = axes('Parent',Fig_); % Create axes
hold(axes_,'on');
t_=title (' ','FontName',options.fontname,'FontSize',16); % Title
t_.String= {'Minimum triangle angles and areas at the mesh surface',str_region_name};
h_=plot(Area_faces,minimum_angle_faces); % Curves
set(h_,'LineStyle','none','Marker','x'); % Colors
xlabel('Facet area - normalized with the voxel length'); % - Axis label
ylabel('Minimum angle (deg)');
grid(axes_,options.grid); % Display grid
set(axes_,'XMinorGrid',options.minorgrid,'YMinorGrid',options.minorgrid); % Display grid for minor thicks also
set(axes_,'FontName',options.fontname,'FontSize',14); % Fontname and fontsize
hold(axes_,'off');
% Save figures
filename = ['Minimum_angle_and_facetarea_' str_region_name_sav];
savefig(Fig_,[fullpath filename]) % .fig
saveas(Fig_,[fullpath filename],'png') % .png
close(Fig_)

% Calculate cumulative and distribution function
parameters_smoothingmovingaverage.round_value = 2;
[results, ~] = Function_probability_density(Area_faces,[],parameters_smoothingmovingaverage);
parameters_distributionfigure.figurename =  ['Facet area at the mesh surface, ' str_region_name];
parameters_distributionfigure.filename = ['Surfacetriangle_area_distribution_' str_region_name_sav];
parameters_distributionfigure.subaxe1_title = ['Mesh surface area: cumulative fct., ' str_region_name];
parameters_distributionfigure.subaxe2_title = ['Mesh surface area: distribution fct., ' str_region_name];
parameters_distributionfigure.xlabel = 'Facet area (normalized with voxel surface)';
data_distributionfigure.cumulative_fct = results.cumulative_fct;
data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
data_distributionfigure.probability_density_fct = results.probability_density_fct;
data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
data_distributionfigure.x50 = results.x50;
data_distributionfigure.smoothed_x50 = results.smoothed_x50;
data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
data_distributionfigure.unit = ''; 
function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);

% Calculate cumulative and distribution function
parameters_smoothingmovingaverage.round_value = 2;
[results, ~] = Function_probability_density(cell_volume,[],parameters_smoothingmovingaverage);
parameters_distributionfigure.figurename = ['Cell volume of ' str_region_name];
parameters_distributionfigure.filename = ['Cell_volume_distribution_' str_region_name_sav];
parameters_distributionfigure.subaxe1_title = ['Cell volume: cumulative fct., ' str_region_name];
parameters_distributionfigure.subaxe2_title = ['Cell volume: distribution fct., ' str_region_name];
parameters_distributionfigure.xlabel = {'Cell volume','normalized with initial voxel cube size'};
data_distributionfigure.cumulative_fct = results.cumulative_fct;
data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
data_distributionfigure.probability_density_fct = results.probability_density_fct;
data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
data_distributionfigure.x50 = results.x50;
data_distributionfigure.smoothed_x50 = results.smoothed_x50;
data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
data_distributionfigure.unit = ''; 
function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);

% distance_transform=bwdist(~binary_array,'euclidean'); % Distance to the boundary (using the cubic mesh);
% % Match cell centroid with distance to the boundary
% [n_cell,~]=size(elem_region); % Number of cells
% cell_distancefromsurface=zeros(n_cell,1);
% position = round(cell_centroid+0.5);
% for k=1:1:n_cell
%     cell_distancefromsurface(k,1)=distance_transform(position(k,1),position(k,2),position(k,3)); % Deduce distance
% end

% Calculate cumulative and distribution function
parameters_smoothingmovingaverage.round_value = 2;
[results, ~] = Function_probability_density(mesh_quality,[],parameters_smoothingmovingaverage);
parameters_distributionfigure.figurename = ['Tetrahedron quality of ' str_region_name];
parameters_distributionfigure.filename = ['Tetrahedron_quality_distribution_' str_region_name_sav];
parameters_distributionfigure.subaxe1_title = ['Tetrahedron quality: cumulative fct., ' str_region_name];
parameters_distributionfigure.subaxe2_title = ['Tetrahedron quality: distribution fct., ' str_region_name];
parameters_distributionfigure.xlabel = {'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'};
data_distributionfigure.cumulative_fct = results.cumulative_fct;
data_distributionfigure.smoothed_cumulative_fct = results.smoothed_cumulative_fct;
data_distributionfigure.probability_density_fct = results.probability_density_fct;
data_distributionfigure.smoothed_probability_density_fct = results.smoothed_probability_density_fct;
data_distributionfigure.x50 = results.x50;
data_distributionfigure.smoothed_x50 = results.smoothed_x50;
data_distributionfigure.integral_probability_density_fct = results.integral_probability_density_fct;
data_distributionfigure.integral_smoothed_probability_density_fct = results.integral_smoothed_probability_density_fct;
data_distributionfigure.unit = ''; 
function_probability_distribution_figure(data_distributionfigure,parameters_distributionfigure);

Fig_= figure; % Create figure
Fig_.Name= ['Tetrahedron quality with cell volume of ' str_region_name];
Fig_.Color='white'; % Background colour
set(Fig_, 'Position', [460 428 660 520]);
axes_ = axes('Parent',Fig_); % Create axes
hold(axes_,'on');
t_=title (' ','FontName',options.fontname,'FontSize',16); % Title
t_.String= {'Tetrahedron quality with cell volume',str_region_name};
h_=plot(cell_volume,mesh_quality); % Curves
set(h_,'LineStyle','none','Marker','x'); % Colors
xlabel('Cell volume (normalized with initial voxel cube size)'); % Axis label
ylabel({'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'});
grid(axes_,options.grid); % Display grid
set(axes_,'XMinorGrid',options.minorgrid,'YMinorGrid',options.minorgrid); % Display grid for minor thicks also
set(axes_,'FontName',options.fontname,'FontSize',14); % - Fontname and fontsize
hold(axes_,'off');
% Save figures
filename = ['Tetrahedron_quality_withcellvolume_' str_region_name_sav];
savefig(Fig_,[fullpath filename]) % .fig
saveas(Fig_,[fullpath filename],'png') % .png
close(Fig_)

%         Fig_= figure; % Create figure
%         Fig_.Name= ['Tetrahedron quality with cell distance from the surface of ' str_region_name];
%         Fig_.Color='white'; % Background colour
%         set(Fig_, 'Position', [460 428 660 520]);
%         axes_ = axes('Parent',Fig_); % Create axes
%         hold(axes_,'on');
%         t_=title (' ','FontName',options.fontname,'FontSize',16); % Title
%         t_.String= {'Tetrahedron quality with cell distance from the surface',str_region_name};
%         h_=plot(cell_distancefromsurface,mesh_quality); % Curves
%         set(h_,'LineStyle','none','Marker','x'); % Colors
%         xlabel('Distance from the surface (normalized with voxel length)'); % - Axis label
%         ylabel({'Tetrahedron quality: ranges ','from 0 (needle-like) to 1 (equilateral tetrahedron)'});
%         grid(axes_,options.grid); % Display grid
%         set(axes_,'XMinorGrid',options.minorgrid,'YMinorGrid',options.minorgrid); % Display grid for minor thicks also
%         set(axes_,'FontName',options.fontname,'FontSize',14); % Fontname and fontsize
%         hold(axes_,'off');
%         % Save figures
%         filename = ['Tetrahedron_quality_withdistancefromsurface_' str_region_name_sav];
%         savefig(Fig_,[fullpath filename]) % .fig
%         saveas(Fig_,[fullpath filename],'png') % .png
%         close(Fig_)
%
%         Fig_= figure; % Create figure
%         Fig_.Name= ['Cell volume along thickness of ' str_region_name];
%         Fig_.Color='white'; % Background colour
%         set(Fig_, 'Position', [460 428 660 520]);
%         axes_ = axes('Parent',Fig_); % Create axes
%         hold(axes_,'on');
%         t_=title (' ','FontName',options.fontname,'FontSize',16); % Title
%         t_.String= {'Cell volume from surface to volume center',str_region_name};
%         h_=plot(cell_distancefromsurface(:,1),cell_volume(:,1)); % Curves
%         set(h_,'LineStyle','none','Marker','x'); % Colors
%         xlabel('Distance from the surface (normalized with voxel length)'); % Axis label
%         ylabel('Cell volume (normalized with initial voxel cube size)');
%         grid(axes_,options.grid); % Display grid
%         set(axes_,'XMinorGrid',options.minorgrid,'YMinorGrid',options.minorgrid); % Display grid for minor thicks also
%         set(axes_,'FontName',options.fontname,'FontSize',14); % Fontname and fontsize
%         hold(axes_,'off');
%         % Save figures
%         filename = ['Cell_volume_fromsrufacetobulk_' str_region_name_sav];
%         savefig(Fig_,[fullpath filename]) % .fig
%         saveas(Fig_,[fullpath filename],'png') % .png
%         close(Fig_)

% 3D figure: selected surface for the specific surface calculation
parameters_3Dpatchfigure.figurename = ['Facets used for the specific surface area, ' str_region_name];
parameters_3Dpatchfigure.title = {'Active interface, facets used are labeled 1',str_activeinterface_a,str_activeinterface_b,str_region_name};
parameters_3Dpatchfigure.filename = [ 'Specificsurfacearea_activeinterface' str_region_name_sav];
data_3Dpatchfigure.FaceVertexCData = selected_surface_activeinterface;
parameters_3Dpatchfigure.colormap = 'cool';
function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)

% 3D figure: selected surface for the specific surface calculation
parameters_3Dpatchfigure.figurename = ['Facets used for the specific surface area, ' str_region_name];
parameters_3Dpatchfigure.title = {'Interface w/o domain''s boundaries, facets used are labeled 1',str_withoutboundary_a,str_withoutboundary_b,str_region_name};
parameters_3Dpatchfigure.filename = [ 'Specificsurfacearea_withoutboundaries' str_region_name_sav];
data_3Dpatchfigure.FaceVertexCData = selected_surface_withoutboundaries;
parameters_3Dpatchfigure.colormap = 'cool';
function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)

% 3D figure: selected surface for the specific surface calculation
parameters_3Dpatchfigure.figurename = ['Facets used for the specific surface area, ' str_region_name];
parameters_3Dpatchfigure.title = {'Interface w/ domain''s boundaries, facets used are labeled 1',str_wholemesh_a,str_wholemesh_b,str_region_name};
parameters_3Dpatchfigure.filename = [ 'Specificsurfacearea_withboundaries' str_region_name_sav];
data_3Dpatchfigure.FaceVertexCData = selected_surface_withboundaries;
parameters_3Dpatchfigure.colormap = 'cool';
function_3Dfigure_withpatch(data_3Dpatchfigure,parameters_3Dpatchfigure)


%% FUNCTIONS

    function [effective_volume] = function_effectivevolume(node_)
        coor_x = node_(:,1); % Node coordinates
        coor_y = node_(:,2);
        coor_z = node_(:,3);
        x_min = min(coor_x); x_max = max(coor_x); % Extremums
        y_min = min(coor_y); y_max = max(coor_y);
        z_min = min(coor_z); z_max = max(coor_z);
        effective_volume = (x_max-x_min) * (y_max-y_min) * (z_max-z_min); % Effective volume
    end

end

