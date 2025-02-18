function [FOV2ndsize, FOV2ndsize_lengthstr] = RVE_2ndsize(pRVE,sz,voxel_size)

if length(sz)==2 || sz(3)==1
    dimension = 2;
    sz(3)=1;
else
    dimension = 3;
end

if strcmp(pRVE.RVEconvergence_Crop,'1 direction')
    if strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 1, from x min to x max') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 1, from x max to x min') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 1, from both ends')
        FOV2ndsize = sz(1);
    elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 2, from x min to x max') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 2, from x max to x min') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 2, from both ends')
        FOV2ndsize = sz(2);
    elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from x min to x max') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from x max to x min') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from both ends')
        FOV2ndsize = sz(3);
    end
    FOV2ndsize_lengthstr = 'length';
    FOV2ndsize = FOV2ndsize * voxel_size;

elseif strcmp(pRVE.RVEconvergence_Crop,'2 directions')
    if dimension == 3
        if strcmp(pRVE.RVEconvergence_thatis,'And keep direction 1 uncropped')
            FOV2ndsize = (sz(2)*sz(3))^(1/2);
        elseif strcmp(pRVE.RVEconvergence_thatis,'And keep direction 2 uncropped')
            FOV2ndsize = (sz(1)*sz(3))^(1/2);
        elseif strcmp(pRVE.RVEconvergence_thatis,'And keep direction 3 uncropped')
            FOV2ndsize = (sz(1)*sz(2))^(1/2);
        end
    else
        if strcmp(pRVE.RVEconvergence_thatis,'And keep direction 1 uncropped')
            error('Should not reach this point')
        elseif strcmp(pRVE.RVEconvergence_thatis,'And keep direction 2 uncropped')
            error('Should not reach this point')
        elseif strcmp(pRVE.RVEconvergence_thatis,'And keep direction 3 uncropped')
            FOV2ndsize = (sz(1)*sz(2))^(1/2);
        end
    end
    FOV2ndsize_lengthstr = 'square root';
    FOV2ndsize = FOV2ndsize * voxel_size;

elseif strcmp(pRVE.RVEconvergence_Crop,'3 directions')
    if dimension == 3
        FOV2ndsize = [];
        FOV2ndsize_lengthstr = [];
    else
        error('Should not reach this point')
    end
end

end