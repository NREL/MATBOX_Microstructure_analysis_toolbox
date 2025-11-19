function [stats,vf,ratio,nvoxel] = Charact_Volumefractions_algorithm(array,label,n_voxel,nanoporosity,wetting)

stats = [];
vf = [];
ratio = [];
nvoxel = [];

%% Number of voxel and background
% No need as we provide n_voxel instead
% n_voxel_background = 0;
% if ~isempty(backgroundlabel)
%     n_voxel_background = sum(sum(sum( array==backgroundlabel )));
% end
% n_voxel = numel(array) - n_voxel_background;

%% Binary image
BWlog = array==label;
BW = single(BWlog);

if isempty(nanoporosity)
    n_voxel_label = sum(sum(sum( BW )));
    nvoxel.n_voxel_label = n_voxel_label;
    vf.phase_label = n_voxel_label/n_voxel;

    stats.nanoporosity = [0 0 0 0 0 0];
    stats.wetting = [0 0 0 0 0 0];
    vf.phase_solid = 0;
    vf.phase_pore_idealwetting = 0;
    vf.phase_pore_partialwetting = 0;

    ratio.solid = 0;
    ratio.poreidealwetting = 0;
    ratio.porepartialwetting = 0;    

    nvoxel.onlypore = 0;
    nvoxel.onlysolid = 0;
    nvoxel.mixed = 0;
    nvoxel.solid = 0;
    nvoxel.pore_idealwetting = 0;
    nvoxel.pore_partialwetting = 0;   
    
else

    %% Bulk quantities
    vals = nanoporosity(BWlog);
    vals = round(double(vals),4);
    vals = reshape(vals,[1 numel(vals)]);
    if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
        stats.nanoporosity = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
    elseif numunique(vals)==0
        stats.nanoporosity = [NaN NaN NaN NaN NaN NaN];
    else
        stats.nanoporosity = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
    end

    vals = wetting(BWlog);
    vals = round(double(vals),4);
    vals = reshape(vals,[1 numel(vals)]);
    if numunique(vals)==1 % Std can provide non-zero (numerical error) if vals is a long array with the same repeating value
        stats.wetting = [vals(1) vals(1) vals(1) 0 vals(1) vals(1)];
    elseif numunique(vals)==0
        stats.wetting = [NaN NaN NaN NaN NaN NaN];
    else
        stats.wetting = [mean(vals) median(vals) mode(vals) std(vals) min(vals) max(vals)];
    end

    %% Volume fractions
    % Label
    n_voxel_label = sum(sum(sum( BW )));
    nvoxel.n_voxel_label = n_voxel_label;
    vf.phase_label = n_voxel_label/n_voxel;

    % Solid
    vf.phase_solid = sum(sum(sum( BW.*(1-nanoporosity) ))) / n_voxel;
    % Pore (ideal and partial wetting)
    vf.phase_pore_idealwetting = sum(sum(sum( BW.*nanoporosity ))) / n_voxel;
    vf.phase_pore_partialwetting = sum(sum(sum( BW.*nanoporosity.*wetting )))/n_voxel;

    %% Ratio
    ratio.solid = vf.phase_solid/vf.phase_label;
    ratio.poreidealwetting = vf.phase_pore_idealwetting/vf.phase_label;
    ratio.porepartialwetting = vf.phase_pore_partialwetting/vf.phase_label;

    %% Number of voxels for categories
    % To be combined with results from other phases
    nvoxel.onlypore = sum(sum(sum( BWlog.*(nanoporosity==1) )));
    nvoxel.onlysolid = sum(sum(sum( BWlog.*(nanoporosity==0) )));

    cond1 = nanoporosity>0;
    cond2 = nanoporosity<1;
    nvoxel.mixed = sum(sum(sum( BWlog.*cond1.*cond2 )));

    %% Number of voxels for true porosity
    % To be combined with results from other phases
    nvoxel.solid = vf.phase_solid * n_voxel;
    nvoxel.pore_idealwetting = vf.phase_pore_idealwetting * n_voxel;
    nvoxel.pore_partialwetting = vf.phase_pore_partialwetting * n_voxel;

end

end
