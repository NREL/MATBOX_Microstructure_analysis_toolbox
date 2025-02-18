function [BW,newtype,foo] = fct_concavityBased_LocalizedErosion(BW,p)

foo=[];
sz = size(BW);
dimension = length(sz);
newtype = 'same';

if isempty(p.savefolder) || isempty(p.filename)
    p.video = false;
end

if p.video
    if ~exist(p.savefolder,'dir') % Folder existence is checked, and created if necessary
        mkdir(p.savefolder);
    end
    video_format = 'mpeg-4';
    fullpath = fullfile(p.savefolder,p.filename);
    video_handle = VideoWriter(fullpath,video_format);
    set(video_handle,'Quality',100);
    set(video_handle,'FrameRate',p.framerate);
    open(video_handle)
    
    if dimension==2
        Fig = figure;
        Fig.Name= 'Concavity-based localized erosion';
        Fig.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig,'position',scrsz);
        tiles = tiledlayout(Fig,2,2);
        tiles.TileSpacing = 'compact';
        tiles.Padding = 'compact';
        xplot = []; yplot_mean = []; yplot_max = [];
        yplot_vf = []; yplot_nparticle = [];
        if strcmp(p.metric,'solidity')
            str_metric = '1 - solidity';
        elseif strcmp(p.metric,'circularity')
            str_metric = '1 - circularity';
        end
    else
        Fig = uifigure;
        viewer = viewer3d('Parent',Fig);
        viewer.RenderingQuality = 'High';
        Vol = volshow(BW,Parent=viewer);
        viewer.CameraPosition = [66.5589  377.6661  183.0159];
    end
end

best_iter = [];
best_BW = BW;
best_m = 1e9;
best_vf = [];
vest_n = [];
col = colororder;
will_break_nextiter = false;
initial_vf = sum(sum(sum(BW)));
%next_save = p.save_tif_each;
Previous_Edges = zeros(sz);
Previous_Asolidity = zeros(sz);
Previous_Acircularity = zeros(sz);
for k_iter=1:1:p.max_iter

    %tic

    %if k_iter==next_save
    %    next_save = next_save + p.save_tif_each;
    %    function_save_tif(uint8(BW),[p.savefolder '\iter_' num2str(k_iter) '.tif']);
    %end

    % Reduce surface roughness
    if p.surface_simplification % Smooth circularity and solidity as a result
        % Erorison dilation
        distance_map = bwdist(~BW,p.surface_simplification_distancemethod);
        BW(distance_map <= p.surface_simplification_distance+0.01)=0;
        distance_map = bwdist(BW,p.surface_simplification_distancemethod);
        BW(distance_map <= p.surface_simplification_distance+0.01)=1;
    end

    % Find edges
    [Edges] = Find_edgesBW(BW); % 0: background, 1: BW, 2: BW and edges of BW
    %[Asolidity, Acircularity] = Calculate_concavitycircularity(BW,Edges,p);
    if p.video
        tmp2 = reshape(Edges,[1 numel(Edges)]);
    end

    if k_iter==1 % Calculate metric every where at the edges
        % Calculate 1-solidity or 1-circularity (measure of concavity)
        [Asolidity, Acircularity] = Calculate_concavitycircularity(BW,Edges,p);
    else % Only update metric near modified edges
        Overlapp_edges = Previous_Edges + Edges; % 4: Edges did not move
        Previous_Edges = Edges; % For next iteration
        % Modified edges
        cond1 = Edges==2;
        cond2 = Overlapp_edges~=4;
        Mod_Edges = zeros(sz);
        Mod_Edges(cond1.*cond2==1)=1;
        % Find all nearby edges
        dmap = bwdist(Mod_Edges,"chessboard");
        cond3 = dmap>p.wndrange+1;
        idx_edges_notmodified = cond1.*cond3==1;
        Edges(idx_edges_notmodified)=3; % 0: background, 1: BW, 2: BW and modified edges of BW, 3: BW and not modified edges of BW
        % Calculate 1-solidity or 1-circularity (measure of concavity)
        [Asolidity, Acircularity] = Calculate_concavitycircularity(BW,Edges,p);
        % Re-introduce values for unmodified edges
        if strcmp(p.metric,'solidity')
            Asolidity(idx_edges_notmodified)=Previous_Asolidity(idx_edges_notmodified);
        elseif strcmp(p.metric,'circularity')
            Acircularity(idx_edges_notmodified)=Previous_Acircularity(idx_edges_notmodified);
        end
    end

    % For next iteration
    Previous_Asolidity = Asolidity;
    Previous_Acircularity = Acircularity; 

    % Video
    if p.video
        if dimension==2
            ax = nexttile(tiles,1);
            if isempty(best_BW)
                imagesc(BW); colormap gray; axis equal; axis tight;
            else
                if sz(1)<2*sz(2)
                    imagesc([BW; best_BW]); colormap gray; axis equal; axis tight;
                else
                    imagesc([BW best_BW]); colormap gray; axis equal; axis tight;
                end
            end                

            ax = nexttile(tiles,2);
            if strcmp(p.metric,'solidity')
                n = length(unique(Asolidity));
                cmap = turbo(n); cmap(1,:)=[1 1 1];
                imagesc(Asolidity); colormap(ax,cmap); colorbar(ax); axis equal; axis tight;
                tmp = reshape(Asolidity,[1 numel(Asolidity)]);
            elseif strcmp(p.metric,'circularity')
                n = length(unique(Acircularity));
                cmap = turbo(n); cmap(1,:)=[1 1 1];
                imagesc(Acircularity); colormap(ax,cmap); colorbar(ax); axis equal; axis tight;
                tmp = reshape(Acircularity,[1 numel(Acircularity)]);
            end
            title(ax,str_metric,'Fontsize',12,'Fontname','Times new roman');
            %tmp2 = reshape(Edges,[1 numel(Edges)]);
            tmp(tmp2~=2)=[]; % Only values at edges

            ax = nexttile(tiles,3);
            cla(ax,"reset");
            xplot = [xplot k_iter];
            yplot_max = [yplot_max max(tmp)];
            yplot_mean = [yplot_mean mean(tmp)];   
            if k_iter == 1
                ylim_max = min([ceil(yplot_max*10)/10 + 0.2 1]);
                ylim_mean = min([ceil(yplot_mean*10)/10 + 0.2 1]);
            end
            hold(ax,'on');
            grid on;
            yyaxis left;
            plot(xplot,yplot_max,'LineWidth',2,'DisplayName','Maxmimum','LineStyle','-','Marker','none');
            ylabel([str_metric ', maximum']);
            ylim([0 ylim_max]);
            if ~isempty(best_iter)
                plot(best_iter,best_m,'LineWidth',2,'Marker','^','MarkerFaceColor','k','MarkerSize',12,'MarkerEdgeColor','k');
            end

            yyaxis right;
            plot(xplot,yplot_mean,'LineWidth',2,'DisplayName','Mean','LineStyle','--','Marker','none');
            ylabel([str_metric ', mean']);
            ylim([0 ylim_mean]);
            % if strcmp(p.metric,'solidity')
            %     plot(xplot,ones(1,numel(xplot))*Asolidity_remove_above,'LineWidth',2,'LineStyle','--','Color','k','DisplayName','Stopping condition');
            % elseif strcmp(p.metric,'circularity')
            %     plot(xplot,ones(1,numel(xplot))*Acircularity_remove_above,'LineWidth',2,'LineStyle','--','Color','k','DisplayName','Stopping condition');
            % end
            xlabel('Iteration');
            set(ax,'Fontsize',12,'Fontname','Times new roman');
            hold(ax,'off');

            ax = nexttile(tiles,4);
            cla(ax,"reset");
            hold(ax,'on');
            yyaxis left;
            yplot_vf = [yplot_vf sum(sum(sum(BW)))/numel(BW)];
            plot(xplot,yplot_vf,'LineWidth',2,'LineStyle','-','Marker','none');
            if ~isempty(best_iter)
                plot(best_iter,best_vf,'LineWidth',2,'Marker','^','MarkerFaceColor','k','MarkerSize',12,'MarkerEdgeColor','k');
            end

            grid on;
            xlabel('Iteration');
            ylabel('Volume fraction');
            set(ax,'Fontsize',12,'Fontname','Times new roman');
            yyaxis right
            C = bwlabeln(BW,4);
            yplot_nparticle = [yplot_nparticle length(unique(C))-1];
            plot(xplot,yplot_nparticle,'LineWidth',2,'LineStyle','--','Marker','none');
            if ~isempty(best_iter)
                plot(best_iter,best_n,'LineWidth',2,'Marker','^','MarkerFaceColor','k','MarkerSize',12,'MarkerEdgeColor','k');
            end            
            ylabel('Number of clusters (particles)');
            hold(ax,'off');

            sgtitle(['Iteration: ' num2str(k_iter,'%i')],'Fontsize',16,'Fontname','Times new roman')

            if k_iter==1
                saveas(Fig,[p.savefolder '\Initial'],'png');
            end

        else
            Fig = uifigure;
            viewer = viewer3d('Parent',Fig);
            viewer.RenderingQuality = 'High';
            Vol = volshow(BW,Parent=viewer);
            viewer.CameraPosition = [200  360  260];
        end

        stored_frame(k_iter) = getframe(Fig);
        writeVideo(video_handle,stored_frame(k_iter))
        drawnow limitrate nocallbacks
    end

    if will_break_nextiter
        break
    end

    % Erode locally based on asolidity or acircularity
    if strcmp(p.metric,'solidity')
        max_asolidity = max(max(max(Asolidity)));
        if max_asolidity < best_m
            best_iter = k_iter;
            best_m = max_asolidity;
            best_BW = BW;
            best_vf = sum(sum(sum(BW)))/numel(BW);
            if dimension==2
                C = bwlabel(BW,4);
            else
                C = bwlabeln(BW,6);
            end
            best_n = length(unique(C))-1;
        end
        if max_asolidity>p.metric_stoppingcondition
            %BW(Asolidity >= p.metric_stoppingcondition) = 0;
            BW(Asolidity >= (1-p.remove_ratiomax) * max_asolidity) = 0;
        else
            if p.video
                will_break_nextiter = true;
            else
                break
            end

        end
    elseif strcmp(p.metric,'circularity')
        max_acircularity = max(max(max(Acircularity)));
        if max_acircularity < best_m
            best_iter = k_iter;
            best_m = max_acircularity;
            best_BW = BW;
            best_vf = sum(sum(sum(BW)))/numel(BW);
            if dimension==2
                C = bwlabel(BW,4);
            else
                C = bwlabeln(BW,6);
            end
            best_n = length(unique(C))-1;            
        end        
        if max_acircularity>p.metric_stoppingcondition
            %BW( Acircularity >= p.metric_stoppingcondition) = 0;
            BW( Acircularity >= (1-p.remove_ratiomax) * max_acircularity) = 0;
        else
            if p.video
                will_break_nextiter = true;
            else
                break
            end
        end
    end

    if sum(sum(sum(BW))) < p.vf_stoppingcondition*initial_vf
        if p.video
            will_break_nextiter = true;
        else
            break
        end
    end

    %toc
end

if p.keepbest
    BW = best_BW;
end

if p.video && dimension==2
    saveas(Fig,[p.savefolder '\Final'],'png');
    savefig(Fig,[p.savefolder '\Final']);
end

end