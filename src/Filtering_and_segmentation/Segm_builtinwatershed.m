function [Instances] = Segm_builtinwatershed(BW,p)

%% REFERENCE
% Use MATLAB built-in function watershed from:
% Meyer, Fernand, "Topographic distance and watershed lines,” Signal Processing , Vol. 38, July 1994, pp. 113-125.

%% DOMAIN SIZE
sz=size(BW);
dimension=length(sz);
if dimension==2 % Case 2D
    visualize_2D = true;
else
    visualize_2D = false;
end

%% WATERSHED ALGORITHM
% We use the matlab built-in function bwdist to calculate the distance transform (or distance map) with an euclidean distance
% bwdist(BW) : For each pixel in BW, the distance transform assigns a number that is the distance between that pixel and the nearest nonzero pixel of BW
% then, bwdist(~BW,'euclidean') :
dmap=bwdist(~BW,'euclidean');

% Topography analogy:
dmap = -dmap; % Inverse the distance sign: 0 is the reference altitude, minimum is the lower altitude point

% Apply watershed
Instances = watershed(dmap);

Instances(~BW)=0;

% Fill watershed lines
tmp = zeros(size(Instances));
tmp(Instances==0)=1;
missingpoints = double(BW).*double(tmp);
id_missingpoints = find(missingpoints);
if ~isempty(id_missingpoints)>0 % Yes, assign to nearest label lake
    [~,idx] = bwdist(~tmp);
    Instances(id_missingpoints) = Instances(idx(id_missingpoints));
end

%% CORRECT IDENTIFICATION
if p.cpsd_refining
    chess_pattern_removal = false; % Not required for Watershed method
    P1 = []; P2 = []; P3 = []; % Not required for Watershed method
    details_convergence = false;
    [Instances] = Function_correct_DPSD_identification(BW, Instances, p.cpsd_refining, details_convergence, chess_pattern_removal, P1, P2, P3, false);

    if visualize_2D
        n_label = length(unique(Instances))-1;
        ax = nexttile(4);
        hold(ax,'on');
        imagesc(ax,Instances);
        cmap(1,:)=[1 1 1];
        colormap(ax,cmap);
        axis(ax,'tight');
        axis(ax,'equal');
        xlabel('Voxels');
        ylabel('Voxels');
        set(ax,'FontName','Times New Roman','FontSize',14); % Fontname and fontsize
        title({'After oversegmentation correction',['Number of bassins identified: ' num2str(n_label,'%i')]},'FontName','Times New Roman','FontSize',14,'Parent',ax);
    end
end

end