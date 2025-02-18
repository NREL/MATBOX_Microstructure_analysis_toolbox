function [FOVbounds] = RVE_convergence_FOVbounds(pRVE,sz)

initial_volume = prod(sz);
minimum_volume = initial_volume*pRVE.RVEconvergence_CropUntil/100;

if length(sz)==2 || sz(3)==1
    dimension = 2;
    sz(3)=1;
else
    dimension = 3;
end

% Initialize
FOVbounds = [];
nFOV = 0;
if strcmp(pRVE.RVEconvergence_Crop,'1 direction')
    volume_toremove = initial_volume*pRVE.RVEconvergence_Ateachstep_x/100;
    nextvolume = initial_volume-volume_toremove;
    while nextvolume>minimum_volume
        nFOV = nFOV+1;
        ratio = (initial_volume - nFOV*volume_toremove)/initial_volume;

        FOVbounds(nFOV,1) = 1;
        FOVbounds(nFOV,2) = sz(1);
        FOVbounds(nFOV,3) = 1;
        FOVbounds(nFOV,4) = sz(2);
        FOVbounds(nFOV,5) = 1;
        FOVbounds(nFOV,6) = sz(3);

        if strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 1, from x min to x max')
            FOVbounds(nFOV,1) = round(sz(1)-sz(1)*ratio);
        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 1, from x max to x min')
            FOVbounds(nFOV,2) = round(sz(1)*ratio);
        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 1, from both ends')
            FOVbounds(nFOV,1) = round(sz(1)/2 - sz(1)*ratio/2);
            FOVbounds(nFOV,2) = round(sz(1)/2 + sz(1)*ratio/2);

        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 2, from y min to y max')
            FOVbounds(nFOV,3) = round(sz(2)-sz(2)*ratio);
        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 2, from y max to y min')
            FOVbounds(nFOV,4) = round(sz(2)*ratio);
        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 2, from both ends')
            FOVbounds(nFOV,3) = round(sz(2)/2 - sz(2)*ratio/2);
            FOVbounds(nFOV,4) = round(sz(2)/2 + sz(2)*ratio/2);

        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from z min to z max') && dimension == 3
            FOVbounds(nFOV,5) = round(sz(3)-sz(3)*ratio);
        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from z max to z min') && dimension == 3
            FOVbounds(nFOV,6) = round(sz(3)*ratio);
        elseif strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from both ends') && dimension == 3
            FOVbounds(nFOV,5) = round(sz(3)/2 - sz(3)*ratio/2);
            FOVbounds(nFOV,6) = round(sz(3)/2 + sz(3)*ratio/2);
        end

        if isempty(FOVbounds)
            warning('RVE convergence options incompatible with a 2D image');
            return
        else
            if dimension == 2
                FOVbounds(nFOV,5) = 1;
                FOVbounds(nFOV,6) = 1;
            end
        end

        nextvolume = initial_volume - (nFOV+1)*volume_toremove;
    end

elseif strcmp(pRVE.RVEconvergence_Crop,'2 directions')
    volume_toremove = initial_volume*pRVE.RVEconvergence_Ateachstep_x/100;
    nextvolume = initial_volume-volume_toremove;
    while nextvolume>minimum_volume
        nFOV = nFOV+1;
        ratio = ((initial_volume - nFOV*volume_toremove)/initial_volume)^(1/2);

        if strcmp(pRVE.RVEconvergence_thatis,'And keep direction 1 uncropped') && dimension == 3
            FOVbounds(nFOV,1) = 1;
            FOVbounds(nFOV,2) = sz(1);
            FOVbounds(nFOV,3) = round(sz(2)/2 - sz(2)*ratio/2);
            FOVbounds(nFOV,4) = round(sz(2)/2 + sz(2)*ratio/2);
            FOVbounds(nFOV,5) = round(sz(3)/2 - sz(3)*ratio/2);
            FOVbounds(nFOV,6) = round(sz(3)/2 + sz(3)*ratio/2);
        elseif strcmp(pRVE.RVEconvergence_thatis,'And keep direction 2 uncropped') && dimension == 3
            FOVbounds(nFOV,1) = round(sz(1)/2 - sz(1)*ratio/2);
            FOVbounds(nFOV,2) = round(sz(1)/2 + sz(1)*ratio/2);
            FOVbounds(nFOV,3) = 1;
            FOVbounds(nFOV,4) = sz(2);
            FOVbounds(nFOV,5) = round(sz(3)/2 - sz(3)*ratio/2);
            FOVbounds(nFOV,6) = round(sz(3)/2 + sz(3)*ratio/2);
        elseif strcmp(pRVE.RVEconvergence_thatis,'And keep direction 3 uncropped')
            FOVbounds(nFOV,1) = round(sz(1)/2 - sz(1)*ratio/2);
            FOVbounds(nFOV,2) = round(sz(1)/2 + sz(1)*ratio/2);
            FOVbounds(nFOV,3) = round(sz(2)/2 - sz(2)*ratio/2);
            FOVbounds(nFOV,4) = round(sz(2)/2 + sz(2)*ratio/2);
            FOVbounds(nFOV,5) = 1;
            FOVbounds(nFOV,6) = sz(3);
        end

        if isempty(FOVbounds)
            warning('RVE convergence options incompatible with a 2D image');
            return
        end

        nextvolume = initial_volume - (nFOV+1)*volume_toremove;
    end

elseif strcmp(pRVE.RVEconvergence_Crop,'3 directions')
    if dimension == 3
        volume_toremove = initial_volume*pRVE.RVEconvergence_Ateachstep_x/100;
        nextvolume = initial_volume-volume_toremove;
        while nextvolume>minimum_volume
            nFOV = nFOV+1;
            ratio = ((initial_volume - nFOV*volume_toremove)/initial_volume)^(1/3);
            FOVbounds(nFOV,1) = round(sz(1)/2 - sz(1)*ratio/2);
            FOVbounds(nFOV,2) = round(sz(1)/2 + sz(1)*ratio/2);
            FOVbounds(nFOV,3) = round(sz(2)/2 - sz(2)*ratio/2);
            FOVbounds(nFOV,4) = round(sz(2)/2 + sz(2)*ratio/2);
            FOVbounds(nFOV,5) = round(sz(3)/2 - sz(3)*ratio/2);
            FOVbounds(nFOV,6) = round(sz(3)/2 + sz(3)*ratio/2);

            next_ratio = ((initial_volume - (nFOV+1)*volume_toremove)/initial_volume)^(1/3);
            next_x0 = round(sz(1)/2 - sz(1)*next_ratio/2);
            next_x1 = round(sz(1)/2 + sz(1)*next_ratio/2);
            next_y0 = round(sz(2)/2 - sz(2)*next_ratio/2);
            next_y1 = round(sz(2)/2 + sz(2)*next_ratio/2);
            next_z0 = round(sz(3)/2 - sz(3)*next_ratio/2);
            next_z1 = round(sz(3)/2 + sz(3)*next_ratio/2);
            nextvolume = (next_x1-next_x0+1)*(next_y1-next_y0+1)*(next_z1-next_z0+1);
        end
    else
        warning('RVE convergence options incompatible with a 2D image');
        return
    end
end