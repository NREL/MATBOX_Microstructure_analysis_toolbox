function [BW, BW_beforeerosiondownscaling] = Segm_erosion_and_downscaling(BW,p)

sz = size(BW);
dimension = length(sz);

%% CHECK DIMENSION COMPATIBILITY
% Crop last few slices so that downscaled image has integer dimension
if p.downscaling && p.downscaling_factor~=1 && round(p.downscaling_factor)==p.downscaling_factor
    tmp = [sz(1)-100:1:sz(1)];
    tmp2 = mod(tmp,p.downscaling_factor);
    idx = max(find(tmp2==0));
    %idx = max(find(tmp2==min(tmp2)));
    cropx = length(tmp)-idx;
    tmp = [sz(2)-100:1:sz(2)];
    tmp2 = mod(tmp,p.downscaling_factor);
    idx = max(find(tmp2==0));
    %idx = max(find(tmp2==min(tmp2)));
    
    cropy = length(tmp)-idx;
    if dimension==2
        BW(end-cropx+1:end, :) = [];
        BW(:, end-cropy+1:end) = [];
    else
        tmp = [sz(3)-100:1:sz(3)];
        tmp2 = mod(tmp,p.downscaling_factor);
        idx = max(find(tmp2==0));
        %idx = max(find(tmp2==min(tmp2)));
        cropz = length(tmp)-idx;
        BW(end-cropx+1:end,:,:) = [];
        BW(:,end-cropy+1:end,:) = [];
        BW(:,:,end-cropz+1:end) = [];
    end   
end


%% DOWNSCALE/EROSION
BW_beforeerosiondownscaling = BW;

for step=1:1:2
    if (step==1 && p.downscaling_before_erosion) || (step==2 && ~p.downscaling_before_erosion)

        if p.downscaling && p.downscaling_factor~=1
            parameters_scaling.scaling_factor = p.downscaling_factor;
            parameters_scaling.label_or_greylevel = 'Label';
            parameters_scaling.background = 0;
            % Scale
            BW = function_scaling(BW,parameters_scaling);
            if p.downscaling_show
                if dimension==2
                    zdisplay=1;
                else
                    sz = size(BW);
                    zdisplay = round(sz(3)/2);
                end
                figure; imagesc(BW(:,:,zdisplay)); axis equal; colormap copper; title('Downscaling'); pause(0.1);
            end
        end

    else

        if p.erosion && p.erosion_distance > 0
            tolerance = 0.1;
            dmap = bwdist(~BW,'Euclidean');
            BW(dmap<=p.erosion_distance+tolerance)=0;
            if p.erosion_display
                if dimension==2
                    zdisplay=1;
                else
                    sz = size(BW);
                    zdisplay = round(sz(3)/2);
                end
                figure; imagesc(BW(:,:,zdisplay)); axis equal; colormap copper; title('Erosion'); pause(0.1);
            end
        end

    end
end


end