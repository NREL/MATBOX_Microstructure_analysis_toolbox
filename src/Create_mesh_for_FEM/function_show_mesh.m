function [] = function_show_mesh(node_region,elem_region,face_region,str_region_name,str_region_name_sav,options)

if  options.save.dommain_plain
    % Figure
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh of ' str_region_name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh of ' str_region_name];
    % - Plot graphs
    hold on;
    plotmesh(node_region,face_region(:,1:3));
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
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end

if options.save.dommain_light
    % Figure
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh of ' str_region_name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh of ' str_region_name];
    % - Plot graphs
    hold on;
    plotmesh(node_region,face_region(:,1:3));
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
    % Light
    camlight left
    lighting flat
    view(3)
    % - Figure has been done
    hold(axes_,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_light_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end

if options.save.dommain_lightnomesh
    % Figure
    Fig_ = figure;
    Fig_.Name=  ['Volumetric mesh of ' str_region_name];
    Fig_.Color='white'; % Background colour
    %set(Fig_, 'Position', [200 100 1000 800]);
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh of ' str_region_name];
    % - Plot graphs
    hold on;
    plotmesh(node_region,face_region(:,1:3),'LineStyle','none');
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
    % Light
    camlight left
    lighting flat
    view(3)
    % Color map
    %colormap('jet');
    % - Figure has been done
    hold(axes_,'off');
    % Save figures
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_noedge_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end

if options.save.dommain_lightnomesh_jet
    % Figure
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh of ' str_region_name];
    Fig_.Color='white'; % Background colour
    %set(Fig_, 'Position', [200 100 1000 800]);
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh of ' str_region_name];
    % - Plot graphs
    hold on;
    plotmesh(node_region,face_region(:,1:3),'LineStyle','none');
    axis equal;
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
    % Light
    camlight left
    lighting flat
    view(3)
    % Color map
    colormap('jet');
    % - Figure has been done
    hold(axes_,'off'); axis tight;
    % Save figures
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_jet_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end

if options.save.dommain_transparent
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh with transparency of ' str_region_name];
    Fig_.Color='white'; % Background colour
    %set(Fig_, 'Position', [200 100 1000 800]);
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh with transparency of ' str_region_name];
    % - Plot graphs
    hold on;
    %plotmesh(node,face(:,1:3));
    plottetra(node_region,elem_region,'facealpha',0.5);
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
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_transparent_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end

if options.save.dommain_slice
    % Find middle position
    x_ = node_region(:,1);
    y_ = node_region(:,2);
    z_ = node_region(:,3);
    x_middle = (min(x_)+max(x_))/2;
    y_middle = (min(y_)+max(y_))/2;
    z_middle = (min(z_)+max(z_))/2;
      
    % Position along x
    str_x=['x<' num2str(x_middle)];
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh - slice along axe 1 of ' str_region_name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_x = axes('Parent',Fig_);
    hold(axes_x,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh - slice along axe 1 of ' str_region_name];
    % - Plot graphs
    hold(axes_x,'on');
    h_x=plotmesh(node_region,elem_region,str_x,'Parent',axes_x);
    n_=length(h_x.FaceVertexCData); % Color scale
    h_x.FaceVertexCData=[1:1:n_]'; % Color scale
    axis(axes_x,'equal')
    axis(axes_x,'tight')
    hold off;
    xlim(axes_x,[0 max(x_)])
    ylim(axes_x,[0 max(y_)])
    zlim(axes_x,[0 max(z_)])
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
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_slice_axe1_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
        
    % Position along y
    str_y=['y<' num2str(y_middle)];
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh - slice along axe 2 of ' str_region_name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_y = axes('Parent',Fig_);
    hold(axes_y,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh - slice along axe 2 of ' str_region_name];
    % - Plot graphs
    hold(axes_y,'on');
    h_y=plotmesh(node_region,elem_region,str_y,'Parent',axes_y);
    n_=length(h_y.FaceVertexCData);
    h_y.FaceVertexCData=[1:1:n_]';
    axis(axes_y,'equal')
    axis(axes_y,'tight')
    hold off;
    xlim(axes_y,[0 max(x_)])
    ylim(axes_y,[0 max(y_)])
    zlim(axes_y,[0 max(z_)])    
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
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_slice_axe2_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
    
    % Position along z
    str_z=['z<' num2str(z_middle)];
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh - slice along axe 3 of ' str_region_name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure     % - Create axes
    axes_z = axes('Parent',Fig_);
    hold(axes_z,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh - slice along axe 3 of ' str_region_name];
    % - Plot graphs
    hold on;
    h_z=plotmesh(node_region,elem_region,str_z,'Parent',axes_z);
    n_=length(h_z.FaceVertexCData);
    h_z.FaceVertexCData=[1:1:n_]';
    axis(axes_z,'equal')
    axis(axes_z,'tight')
    hold off;
    xlim(axes_z,[0 max(x_)])
    ylim(axes_z,[0 max(y_)])
    zlim(axes_z,[0 max(z_)])    
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
    fullpath=[options.save.mainfolder options.save.vieweachdomainfolder];
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    filename = ['Volumetric_mesh_slice_axe3_' str_region_name_sav];
    % .fig
    if options.save.dommain_fig
        savefig(Fig_,[fullpath filename])
    end
    % .png
    if options.save.domain_png
        saveas(Fig_,[fullpath filename],'png')
    end
    close(Fig_)
end

end

