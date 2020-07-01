function [] = function_get_commonnode_interface(node_electrolyte,node_leftelectrode,face_leftelectrode, node_rightelectrode,face_rightelectrode,options)

% Tolerance
tolerance_interface = options.tolerance_interface;

% Electrolyte with left electrode
node_electrolyte_belong_Leftelectrode=ismembertol(node_electrolyte,node_leftelectrode,tolerance_interface,'ByRows',true); % find node that belong both to electrolyte and left electrode
% Find index
index= node_electrolyte_belong_Leftelectrode==1;
% Extract coordinates
Node_coordinate_interface_LeftAM_EL=node_electrolyte(index,:); % Express in voxel length. Array will be re-scaled in FEniCS
% Number of verteces at the interface
[number_common_vertece_interface_LeftAM_EL,~]=size(Node_coordinate_interface_LeftAM_EL);

% Electrolyte with right electrode
node_electrolyte_belong_RightElectrode=ismembertol(node_electrolyte,node_rightelectrode,tolerance_interface,'ByRows',true); % find node that belong both to electrolyte and left electrode
% Find index
index= node_electrolyte_belong_RightElectrode==1;
% Extract coordinates
Node_coordinate_interface_EL_RightAM=node_electrolyte(index,:); % Express in voxel length. Array will be re-scaled in FEniCS
% Number of verteces at the interface
[number_common_vertece_interface_EL_RightAM,~]=size(Node_coordinate_interface_EL_RightAM);

disp '>>> Verteces at the interface Left electrode - Electrolyte'
fprintf('Number of common verteces: %i\n',number_common_vertece_interface_LeftAM_EL)
disp '>>> Verteces at the interface Electrolyte - Right electrode'
fprintf('Number of common verteces: %i\n',number_common_vertece_interface_EL_RightAM)
disp ' '

% Save points
if options.save.nodeinterface_mat
    fullpath=[options.save.mainfolder options.save.meshdatafolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    save([fullpath,'Verteces_coordinate_interface_LeftAM_EL.mat'],'Node_coordinate_interface_LeftAM_EL')
    save([fullpath,'Verteces_coordinate_interface_EL_RightAM.mat'],'Node_coordinate_interface_EL_RightAM')
end

% Figure
Fig_ = figure;
Fig_.Name= 'Common vertces at the interfaces';
Fig_.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig_,'position',scrsz); % Full screen figure
% - Create axes
axes_ = axes('Parent',Fig_);
hold(axes_,'on');
% - Title
t_=title (' ','FontName','Times New Roman','FontSize',16);
t_.String= 'Common vertces at the interfaces';
% - Plot graphs
hold on;
h_leftelectrode_electrolyte = scatter3(Node_coordinate_interface_LeftAM_EL(:,1),Node_coordinate_interface_LeftAM_EL(:,2),Node_coordinate_interface_LeftAM_EL(:,3));
h_electrolyte_rightelectrode = scatter3(Node_coordinate_interface_EL_RightAM(:,1),Node_coordinate_interface_EL_RightAM(:,2),Node_coordinate_interface_EL_RightAM(:,3));
h_leftelectrode=plotmesh(node_leftelectrode,face_leftelectrode(:,1:3));
h_rightelectrode=plotmesh(node_rightelectrode,face_rightelectrode(:,1:3));
% Set color and transparency
h_leftelectrode.FaceAlpha=0.5;
h_rightelectrode.FaceAlpha=0.5;
h_leftelectrode_electrolyte.MarkerEdgeColor='b';
h_electrolyte_rightelectrode.MarkerEdgeColor='r';
axis equal; axis tight;
% Color map
colormap('gray');
hold off;
% - Axis label
xlabel('Axe 1 (voxels length)');
ylabel('Axe 2 (voxels length)');
zlabel('Axe 3 (voxels length)');
% Create legend
lgd = legend(axes_,'show',{'Left active material - Electrolyte','Electrolyte - Right active material','Left and right electrodes'},'Location','best');
title(lgd,'Common verteces at interfaces')
% - Grid
grid(axes_,'on'); % Display grid
set(axes_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
% - Fontname and fontsize
set(axes_,'FontName','Times New Roman','FontSize',14);
% - Figure has been done
hold(axes_,'off');
% Save figures
fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
if ~exist(fullpath,'dir')
    mkdir(fullpath);
end
filename = 'Verteces_at_the_interfaces';
% .fig
if options.save.nodeinterface_fig
    savefig(Fig_,[fullpath filename])
end
% .png
if options.save.nodeinterface_png
    saveas(Fig_,[fullpath filename],'png')
end
close(Fig_)

end

