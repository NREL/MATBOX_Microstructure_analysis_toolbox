function [BW] = ConcavityBased_LocalizedErosion(BW,p)

if p.video
    video_format = 'mpeg-4';
    fullpath = fullfile(p.savefolder,'Erosion video');
    video_handle = VideoWriter(fullpath,video_format);
    set(video_handle,'Quality',100);
    set(video_handle,'FrameRate',p.framerate);
    open(video_handle)

    dimension = length(size(BW));
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
        if p.check_asolidity
            str_metric = '1 - solidity';
        elseif p.check_acircularity
            str_metric = '1 - circularity';
        end
    else
        Fig = uifigure;
        viewer = viewer3d('Parent',Fig);
        viewer.RenderingQuality = 'High';
        Vol = volshow(BW,Parent=viewer);
        viewer.CameraPosition = [89.2395  -64.1046   99.4612];
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
next_save = p.save_tif_each;
for k_iter=1:1:p.max_iter

    if k_iter==next_save
        next_save = next_save + p.save_tif_each;
        function_save_tif(uint8(BW),[p.savefolder '\iter_' num2str(k_iter) '.tif']);
    end

    % Reduce surface roughness
    if p.surface_simplification % Smooth circularity and solidity as a result
        % Erorison dilation
        distance_map = bwdist(~BW,"chessboard");
        BW(distance_map<=1.01)=0;
        distance_map = bwdist(BW,"chessboard");
        BW(distance_map<=1.01)=1;
    end

    % Find edges
    [Edges] = Find_edgesBW(BW);

    % Calculate 1-solidity or 1-circularity (measure of concavity)
    [Asolidity, Acircularity] = Calculate_concavitycircularity(BW,Edges,p);

    % Video
    if p.video
        if dimension==2
            ax = nexttile(tiles,1);
            if isempty(best_BW)
                imagesc(BW); colormap gray; axis equal; axis tight;
            else
                imagesc([BW best_BW]); colormap gray; axis equal; axis tight;
            end                

            ax = nexttile(tiles,2);
            if p.check_asolidity
                n = length(unique(Asolidity));
                cmap = turbo(n); cmap(1,:)=[1 1 1];
                imagesc(Asolidity); colormap(ax,cmap); colorbar(ax); axis equal; axis tight;
                tmp = reshape(Asolidity,[1 numel(Asolidity)]);
            elseif p.check_acircularity
                n = length(unique(Acircularity));
                cmap = turbo(n); cmap(1,:)=[1 1 1];
                imagesc(Acircularity); colormap(ax,cmap); colorbar(ax); axis equal; axis tight;
                tmp = reshape(Acircularity,[1 numel(Acircularity)]);
            end
            title(ax,str_metric,'Fontsize',12,'Fontname','Times new roman');
            tmp2 = reshape(Edges,[1 numel(Edges)]);
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
            % if p.check_asolidity
            %     plot(xplot,ones(1,numel(xplot))*Asolidity_remove_above,'LineWidth',2,'LineStyle','--','Color','k','DisplayName','Stopping condition');
            % elseif p.check_acircularity
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
            viewer.CameraPosition = [89.2395  -64.1046   99.4612];
        end

        stored_frame(k_iter) = getframe(Fig);
        writeVideo(video_handle,stored_frame(k_iter))
    end

    if will_break_nextiter
        break
    end

    % Erode locally based on asolidity or acircularity
    if p.check_asolidity
        max_asolidity = max(max(max(Asolidity)));
        if max_asolidity < best_m
            best_iter = k_iter;
            best_m = max_asolidity;
            best_BW = BW;
            best_vf = sum(sum(sum(BW)))/numel(BW);
            C = bwlabeln(BW,4);
            best_n = length(unique(C))-1;
        end
        if max_asolidity>p.Asolidity_remove_above
            %BW(Asolidity >= p.Asolidity_remove_above) = 0;
            BW(Asolidity >= p.Asolidity_remove_ratiomax * max_asolidity) = 0;
        else
            if p.video
                will_break_nextiter = true;
            else
                break
            end

        end
    elseif p.check_acircularity
        max_acircularity = max(max(max(Acircularity)));
        if max_acircularity < best_m
            best_iter = k_iter;
            best_m = max_acircularity;
            best_BW = BW;
            best_vf = sum(sum(sum(BW)))/numel(BW);
            C = bwlabeln(BW,4);
            best_n = length(unique(C))-1;            
        end        
        if max_acircularity>p.Acircularity_remove_above
            %BW( Acircularity >= p.Acircularity_remove_above) = 0;
            BW( Acircularity >= p.Acircularity_remove_ratiomax * max_acircularity) = 0;
        else
            if p.video
                will_break_nextiter = true;
            else
                break
            end
        end
    end

    if sum(sum(sum(BW))) < p.min_vf*initial_vf
        if p.video
            will_break_nextiter = true;
        else
            break
        end
    end
end

BW = best_BW;

saveas(Fig,[p.savefolder '\Final'],'png');
savefig(Fig,[p.savefolder '\Final']);
end