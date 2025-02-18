function [] = function_showmesh(node_,elem_,face_,subdomain_,str_mesh_name,options)
% plotmesh and plottetra functions are from Iso2mesh package

% Save options
opt.savefig_infig = options.savefig;
if options.savepng
    opt.savefig_informat = {'png'};
else
    opt.savefig_informat = [];
end
opt.fig_infig = options.fig_infig;
opt.overwritte = false;
opt.fig_format = {'png'};

if strcmp(options.coloris,'Phase Id')
    % Initialize
    [nface,~] = size(face_);
    face_phase_id = zeros(nface,1);
    
    % Identify face belongs to cell
    s_face = sort(face_,2);
    s_e_123 = sort([elem_(:,1) elem_(:,2) elem_(:,3)],2);
    s_e_124 = sort([elem_(:,1) elem_(:,2) elem_(:,4)],2);
    s_e_134 = sort([elem_(:,1) elem_(:,3) elem_(:,4)],2);
    s_e_234 = sort([elem_(:,2) elem_(:,3) elem_(:,4)],2);
    
    [Lia, Lib] = ismember(s_face, s_e_123, 'rows');
    id = Lia==1;
    face_phase_id(id) = subdomain_(Lib(id));

    [Lia, Lib] = ismember(s_face, s_e_124, 'rows');
    id = Lia==1;
    face_phase_id(id) = subdomain_(Lib(id));
    
    [Lia, Lib] = ismember(s_face, s_e_134, 'rows');
    id = Lia==1;
    face_phase_id(id) = subdomain_(Lib(id));
    
    [Lia, Lib] = ismember(s_face, s_e_234, 'rows');
    id = Lia==1;
    face_phase_id(id) = subdomain_(Lib(id));    
    
    % Assign face color based on subdomain id
    unique_id = unique(subdomain_);
    n_id = length(unique_id);
    if n_id>10 % Likely polychristalline
        c = rand(n_id,1);
    else
        c = linspace(0.1,0.9,n_id)';
    end
    face_color = zeros(1,nface);
    for k=1:1:n_id
        idx = find( face_phase_id==unique_id(k));
        face_color(1,idx) = c(k,1);
    end
    
    % Is this a bug? Annoying but face being only required for visualization and not meshing, it's not harmful.
    id_emptyface = find(face_color==0);
    face_color(:,id_emptyface)=[];
    face_(id_emptyface,:)=[];
end

if strcmp(options.show,'Volume')
    Fig_ = figure;
    Fig_.Name= ['Volumetric mesh of ' str_mesh_name];
    Fig_.Color='white'; % Background colour
    scrsz = get(0,'ScreenSize'); % Screen resolution
    set(Fig_,'position',scrsz); % Full screen figure
    % - Create axes
    axes_ = axes('Parent',Fig_);
    hold(axes_,'on');
    % - Title
    t_=title (' ','FontName','Times New Roman','FontSize',16);
    t_.String= ['Volumetric mesh of ' str_mesh_name];
    % - Plot graphs
    if options.transparency == 0
        if strcmp(options.edges,'Visible')
            h=plotmesh(node_,face_(:,1:3));
        elseif strcmp(options.edges,'Not visible')
            h=plotmesh(node_,face_(:,1:3),'LineStyle','none');
        end
    else
        h=plottetra(node_,elem_,'facealpha',options.transparency);
    end
    axis equal; axis tight;
    % - Axis label
    if options.axislabel
        xlabel('Axe 1 (voxels length)');
        ylabel('Axe 2 (voxels length)');
        zlabel('Axe 3 (voxels length)');
    end
    % - Grid
    if options.grid
        grid(axes_,'on'); % Display grid
        set(axes_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
    end
    % - Fontname and fontsize
    set(axes_,'FontName','Times New Roman','FontSize',14);
    % Light
    camlight left
    lighting(axes_, options.lighting)
    view(3)
    % Color map
    if strcmp(options.coloris,'Phase Id')
        h.CData = face_color(1,:);
    end
    colormap(options.colormap);
    % - Figure has been done
    hold(axes_,'off');
    % Save
    function_savefig(Fig_, options.folder, str_mesh_name, opt);
end

if strcmp(options.show,'Slice')
    % Find middle position
    x_ = node_(:,1);
    y_ = node_(:,2);
    z_ = node_(:,3);
    x_middle = (min(x_)+max(x_))/2;
    y_middle = (min(y_)+max(y_))/2;
    z_middle = (min(z_)+max(z_))/2;
    
    for k=1:1:3
        if k==1
            str_x=['x<' num2str(x_middle)]; %  Position along x
        elseif k==2
            str_x=['y<' num2str(y_middle)]; %  Position along y
        elseif k==3
            str_x=['z<' num2str(z_middle)]; %  Position along z
        end
        Fig_ = figure;
        Fig_.Name= ['Volumetric mesh of ' str_mesh_name ' - slice along axe ' num2str(k)];
        Fig_.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig_,'position',scrsz); % Full screen figure
        % - Create axes
        axes_ = axes('Parent',Fig_);
        hold(axes_,'on');
        % - Title
        t_=title (' ','FontName','Times New Roman','FontSize',16);
        t_.String= ['Volumetric mesh of ' str_mesh_name ' - slice along axe ' num2str(k)];
        % - Plot graphs
        h=plotmesh(node_,elem_,str_x,'Parent',axes_);
        n_=length(h.FaceVertexCData); % Color scale
        h.FaceVertexCData=[1:1:n_]'; % Color scale
        axis(axes_,'equal')
        axis(axes_,'tight')
        xlim(axes_,[0 max(x_)])
        ylim(axes_,[0 max(y_)])
        zlim(axes_,[0 max(z_)])
        % - Axis label
        if options.axislabel
            xlabel('Axe 1 (voxels length)');
            ylabel('Axe 2 (voxels length)');
            zlabel('Axe 3 (voxels length)');
        end
        % - Grid
        if options.grid
            grid(axes_,'on'); % Display grid
            set(axes_,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on'); % Display grid for minor thicks also
        end
        % - Fontname and fontsize
        set(axes_,'FontName','Times New Roman','FontSize',14);
        % - Figure has been done
        if k==1
            view(axes_,[90,0])
        elseif k==2
            view(axes_,[180,0])
        elseif k==3
            view(axes_,[0,90])
        end
        % Light
        camlight left
        lighting(axes_, options.lighting)
        % Color map
        colormap(options.colormap);
        hold(axes_,'off');
        % Save figures
        %function_savefig(Fig_, options.folder, [str_mesh_name '_axe_' num2str(k)], opt);
    end
end

end