function [node,elem,face] = mesh_generation_withIso2mesh(array_, options)

%% CREATE SURFACE MESH
% [no,el,regions,holes]=v2s(image_,isovalues,opt,method);

% image_   : a LOGICAL array
% isovalues: a list of isovalues where the levelset is defined

% opt      :  In v2s, method can be set to 'cgalmesh' in addition to those allowed by vol2surf.
%   if method is 'cgalsurf' or 'cgalpoly':
% 	     opt=a float number>1: max radius of the Delaunay sphere(element size)
% 	     opt.radbound: same as above, max radius of the Delaunay sphere
% 	     opt.distbound: maximum deviation from the specified isosurfaces
% 	     opt(1,2,...).radbound: same as above, for each levelset
% 	if method is 'simplify':
% 	     opt=a float number<1: compression rate for surf. simplification
% 	     opt.keepratio=a float less than 1: same as above, same for all surf.
% 	     opt(1,2,..).keepratio: setting compression rate for each levelset
% 	opt(1,2,..).maxsurf: 1 - only use the largest disjointed surface
% 				0 - use all surfaces for that levelset
%   opt(1,2,..).side: - 'upper': threshold at upper interface
%                       'lower': threshold at lower interface
% 	opt(1,2,..).maxnode: - the maximum number of surface node per levelset
% 	opt(1,2,..).holes: user specified holes interior pt list
% 	opt(1,2,..).regions: user specified regions interior pt list
% 	opt(1,2,..).surf.{node,elem}: add additional surfaces
% 	opt(1,2,..).{A,B}: linear transformation for each surface
% 	opt.autoregion: if set to 1, vol2surf will try to determine
%                   the interior points for each closed surface automatically

% In v2s, method can be set to 'cgalmesh' in addition to those allowed by vol2surf.
% method  : - if method is 'simplify', iso2mesh will first call
% 		      binsurface to generate a voxel-based surface mesh and then
%             use meshresample/meshcheckrepair to create a coarser mesh;
% 		    - if method is 'cgalsurf', iso2mesh will call the surface
% 		      extraction program from CGAL to make surface mesh
% 		    - if method is not specified, 'cgalsurf' is assumed by default

% Function parameter
method_surfacemesh = options.method_surfacemesh;
opt.radbound=options.radbound;
opt.distbound=options.distbound;
isovalues_=options.isovalues;

tmp=uint8(array_); % Convert in unsigned integer 8bits
[nodes_coordinates,triangular_faces,regions,holes]=v2s(tmp,isovalues_,opt,method_surfacemesh);


if options.save.initialsurfacemesh
    Fig_ = figure; % Figure
    Fig_.Name= 'Initial triangle-based surface mesh';
    Fig_.Color='white'; % Background colour
    set(Fig_, 'Position', [200 100 1000 800]);
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= 'Initial triangle-based surface mesh';
    % - Plot graphs
    hold on;
    plotmesh(nodes_coordinates,triangular_faces);
    axis equal; axis tight;
    hold off;
    % - Axis label
    xlabel('Axe 1 (voxels length)');
    ylabel('Axe 2 (voxels length)');
    zlabel('Axe 3 (voxels length)');
    % - Grid
    grid(axes_,'on'); % Display grid
    set(axes_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(axes_,'FontName','Times New Roman','FontSize',14);
    % - Figure has been done
    hold(axes_,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.viewalldomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = 'Initial_surface_mesh';
    % .fig
    if options.save.alldomain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.alldomain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end


%% REGION VERIFICATION

number_region=length(regions); % Number of region
number_holes=length(holes); % Number of region
fprintf ('- SURFACE MESH CREATED with v2s: number of region %i \n',(number_region));
fprintf ('                               : number of holes %i \n',(number_holes));
disp 'list of interior points for all sub-region, [x,y,z]:'
disp(regions)


%% SMOOTH SURFACE MESH

% cell structure of length size(node), conn{n} contains a list of all neighboring node ID for node n
conn=meshconn(triangular_faces(:,1:3),size(nodes_coordinates,1));
mask_=[];

% Parameters and call

if strcmp(options.method_smoothsurfacemesh,'lowpass')
    smooth_method = 'lowpass';
    iteration_smoothing = options.iteration;
    useralpha = options.useralpha;
    smoothed_nodes_coordinates=smoothsurf(nodes_coordinates,mask_,conn,iteration_smoothing,useralpha,smooth_method);
elseif strcmp(options.method_smoothsurfacemesh,'laplacian')
    smooth_method = 'laplacian';
    iteration_smoothing = options.iteration;
    useralpha = options.useralpha;
    smoothed_nodes_coordinates=smoothsurf(nodes_coordinates,mask_,conn,iteration_smoothing,useralpha,smooth_method);
elseif strcmp(options.method_smoothsurfacemesh,'laplacianhc')
    smooth_method = 'laplacianhc';
    iteration_smoothing = options.iteration;
    useralpha = options.useralpha;
    userbeta = options.userbeta;
    smoothed_nodes_coordinates=smoothsurf(nodes_coordinates,mask_,conn,iteration_smoothing,useralpha,smooth_method);
end

if options.save.smoothedsurfacemesh
    Fig_ = figure; % Figure
    Fig_.Name= 'Smoothed surface mesh';
    Fig_.Color='white'; % Background colour
    set(Fig_, 'Position', [200 100 1000 800]);
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= 'Smoothed surface mesh';
    % - Plot graphs
    hold on;
    plotmesh(smoothed_nodes_coordinates,triangular_faces);
    axis equal; axis tight;
    hold off;
    % - Axis label
    xlabel('Axe 1 (voxels length)');
    ylabel('Axe 2 (voxels length)');
    zlabel('Axe 3 (voxels length)');
    % - Grid
    grid(axes_,'on'); % Display grid
    set(axes_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(axes_,'FontName','Times New Roman','FontSize',14);
    % - Figure has been done
    hold(axes_,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.viewalldomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = 'Smoothed_surface_mesh';
    % .fig
    if options.save.alldomain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.alldomain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close 'Smoothed surface mesh'
end

%% CREATE VOLUMETRIC MESH

Domain_size = size(array_);
x0=0; y0=0; z0=0;
x1=Domain_size(1); y1=Domain_size(2); z1=Domain_size(3);

% Function parameter
keepratio = options.keepratio;
maxvol = options.maxvol;
[node,elem,face]=surf2mesh(smoothed_nodes_coordinates,triangular_faces,[x0-1 y0-1 z0-1],[x1+1 y1+1 z1+1],keepratio,maxvol,regions);


%% PLOT VOLUMETRIC MESH

if options.save.volumetricmesh
    Fig_ = figure;
    Fig_.Name= 'Volumetric mesh (with transparency)';
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure% - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= 'Volumetric mesh (with transparency)';
    % - Plot graphs
    hold on;
    %plotmesh(node,face(:,1:3));
    plottetra(node,elem,'facealpha',0.5);
    axis equal; axis tight;
    hold off;
    % - Axis label
    xlabel('Axe 1 (voxels length)');
    ylabel('Axe 2 (voxels length)');
    zlabel('Axe 3 (voxels length)');
    % - Grid
    grid(axes_,'on'); % Display grid
    set(axes_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(axes_,'FontName','Times New Roman','FontSize',14);
    % - Figure has been done
    hold(axes_,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.viewalldomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = 'Volumetric_mesh_transparent';
    % .fig
    if options.save.alldomain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.alldomain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_);
end

%% MESH VIEW SLICE PER SLICE

if options.save.volumetricmesh_slices
    % Position along x
    str_x=['x<' num2str(round(Domain_size(1)/2))];
    Fig_ = figure;
    Fig_.Name= 'Volumetric mesh - slice along axe 1';
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_x = axes('Parent',Fig_);
    hold(axes_x,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= 'Volumetric mesh - slice along axe 1';
    % - Plot graphs
    hold(axes_x,'on');
    h_x=plotmesh(node,elem,str_x,'Parent',axes_x);
    n_=length(h_x.FaceVertexCData); % Color scale
    h_x.FaceVertexCData=[1:1:n_]'; % Color scale
    axis(axes_x,'equal')
    axis(axes_x,'tight')
    hold off;
    xlim(axes_x,[0 Domain_size(1)])
    ylim(axes_x,[0 Domain_size(2)])
    zlim(axes_x,[0 Domain_size(3)])
    % - Axis label
    xlabel('Axe 1 (voxels length)');
    ylabel('Axe 2 (voxels length)');
    zlabel('Axe 3 (voxels length)');
    % - Grid
    grid(axes_x,'on'); % Display grid
    set(axes_x,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(axes_x,'FontName','Times New Roman','FontSize',14);
    % - Figure has been done
    view(axes_x,[90,0])
    hold(axes_x,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.viewalldomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end    
    filename = 'Volumetric_mesh_slice_axe1';
    % .fig
    if options.save.alldomain_slice_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.alldomain_slice_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_);
    
    % Position along y
    str_y=['y<' num2str(round(Domain_size(2)/2))];
    Fig_ = figure;
    Fig_.Name= 'Volumetric mesh - slice along axe 2';
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_y = axes('Parent',Fig_);
    hold(axes_y,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= 'Volumetric mesh - slice along axe 2';
    % - Plot graphs
    hold(axes_y,'on');
    h_y=plotmesh(node,elem,str_y,'Parent',axes_y);
    n_=length(h_y.FaceVertexCData);
    h_y.FaceVertexCData=[1:1:n_]';
    axis(axes_y,'equal')
    axis(axes_y,'tight')
    hold off;
    xlim(axes_y,[0 Domain_size(1)])
    ylim(axes_y,[0 Domain_size(2)])
    zlim(axes_y,[0 Domain_size(3)])
    % - Axis label
    xlabel('Axe 1 (voxels length)');
    ylabel('Axe 2 (voxels length)');
    zlabel('Axe 3 (voxels length)');
    % - Grid
    grid(axes_y,'on'); % Display grid
    set(axes_y,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(axes_y,'FontName','Times New Roman','FontSize',14);
    % - Figure has been done
    view(axes_y,[180,0])
    hold(axes_y,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.viewalldomainfolder];
    filename = 'Volumetric_mesh_slice_axe2';
    % .fig
    if options.save.alldomain_slice_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.alldomain_slice_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_);
   
    % Position along z
    str_z=['z<' num2str(round(Domain_size(3)/2))];
    Fig_ = figure;
    Fig_.Name= 'Volumetric mesh - slice along axe 3';
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_z = axes('Parent',Fig_);
    hold(axes_z,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= 'Volumetric mesh - slice along axe 3';
    % - Plot graphs
    hold on;
    h_z=plotmesh(node,elem,str_z,'Parent',axes_z);
    n_=length(h_z.FaceVertexCData);
    h_z.FaceVertexCData=[1:1:n_]';
    axis(axes_z,'equal')
    axis(axes_z,'tight')
    hold off;
    xlim(axes_z,[0 Domain_size(1)])
    ylim(axes_z,[0 Domain_size(2)])
    zlim(axes_z,[0 Domain_size(3)])
    % - Axis label
    xlabel('Axe 1 (voxels length)');
    ylabel('Axe 2 (voxels length)');
    zlabel('Axe 3 (voxels length)');
    % - Grid
    grid(axes_z,'on'); % Display grid
    set(axes_z,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    % - Fontname and fontsize
    set(axes_z,'FontName','Times New Roman','FontSize',14);
    % - Figure has been done
    view(axes_z,[0,90])
    hold(axes_z,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.viewalldomainfolder];
    filename = 'Volumetric_mesh_slice_axe3';
    % .fig
    if options.save.alldomain_slice_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.alldomain_slice_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_);
end

end

