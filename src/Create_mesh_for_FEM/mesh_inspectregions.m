function [] = mesh_inspectregions(region_,options)

number_region_vol = length(region_);
fullpath=[options.save.mainfolder options.save.viewregion];
if ~exist(fullpath,'dir')
    mkdir(fullpath);
else
    A=dir(fullpath);
    for k = 1:length(A)
        delete([ fullpath  A(k).name]);
    end
end
 
%% PLOT REGION
for k=1:1:number_region_vol
    
    node_tmp = region_(k).node;
    face_tmp = region_(k).face;
    
    % Figure
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh, region ' num2str(k)];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh, region ' num2str(k)];
    % - Plot graphs
    hold on;
    plotmesh(node_tmp,face_tmp(:,1:3));
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
    if options.save.regions_fig == 1 || options.save.regions_png == 1
        filename = ['Volumetric_mesh_region' num2str(k)];
        % .fig
        if options.save.regions_fig
            savefig(Fig_,[fullpath filename])
        end
        % .png
        if options.save.regions_png
            saveas(Fig_,[fullpath filename],'png')
        end
    end
end

end

