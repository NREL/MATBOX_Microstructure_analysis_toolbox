function [n_nestedRVE,xcrop] = Function_NestedRVE2(RVEparameters,Domain_size)

n_nestedRVE = 0;
initial_volume = prod(Domain_size);
minimum_volume = initial_volume*RVEparameters.RVEconvergence_CropUntil/100;

if strcmp(RVEparameters.RVEconvergence_Crop,'1 direction')
   
        volume_toremove = initial_volume*RVEparameters.RVEconvergence_Ateachstep_x/100;
        nextvolume = initial_volume-volume_toremove;
        while nextvolume>minimum_volume
            n_nestedRVE = n_nestedRVE+1;
            ratio = (initial_volume - n_nestedRVE*volume_toremove)/initial_volume;
            xcrop(n_nestedRVE,1) = 1;
            xcrop(n_nestedRVE,2) = Domain_size(1);
            xcrop(n_nestedRVE,3) = 1;
            xcrop(n_nestedRVE,4) = Domain_size(2);
            xcrop(n_nestedRVE,5) = 1;
            xcrop(n_nestedRVE,6) = Domain_size(3);

            if strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 1, from x min to x max')
                xcrop(n_nestedRVE,1) = round(Domain_size(1)-Domain_size(1)*ratio);
            elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 1, from x max to x min')
                xcrop(n_nestedRVE,2) = round(Domain_size(1)*ratio);
            elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 1, from both ends')
                xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);

            elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 2, from y min to y max')
                xcrop(n_nestedRVE,3) = round(Domain_size(2)-Domain_size(2)*ratio);
            elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 2, from y max to y min')
                xcrop(n_nestedRVE,4) = round(Domain_size(2)*ratio);
            elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 2, from both ends')
                xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
            end

            if Domain_size(3)>1
                if strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 3, from z min to z max')
                    xcrop(n_nestedRVE,5) = round(Domain_size(3)-Domain_size(3)*ratio);
                elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 3, from z max to z min')
                    xcrop(n_nestedRVE,6) = round(Domain_size(3)*ratio);
                elseif strcmp(RVEparameters.RVEconvergence_thatis,'Shrink along direction 3, from both ends')
                    xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
                    xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);
                end
            else
                n_nestedRVE = -1;
                warning('RVE convergence options incompatible with a 2D image');
                return
            end

            nextvolume = initial_volume - (n_nestedRVE+1)*volume_toremove;
        end    

elseif strcmp(RVEparameters.RVEconvergence_Crop,'2 directions')

    if strcmp(RVEparameters.RVEconvergence_Ateachstep,'remove x %vol of the initial volume ')
        volume_toremove = initial_volume*RVEparameters.RVEconvergence_Ateachstep_x/100;
        nextvolume = initial_volume-volume_toremove;
        while nextvolume>minimum_volume
            n_nestedRVE = n_nestedRVE+1;
            ratio = ((initial_volume - n_nestedRVE*volume_toremove)/initial_volume)^(1/2);
            if strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 1 uncropped')
                xcrop(n_nestedRVE,1) = 1;
                xcrop(n_nestedRVE,2) = Domain_size(1);
                xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
                xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);
            elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 2 uncropped')
                xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,3) = 1;
                xcrop(n_nestedRVE,4) = Domain_size(2);
                xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
                xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);
            elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 3 uncropped')
                xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,5) = 1;
                xcrop(n_nestedRVE,6) = Domain_size(3);
            end
            nextvolume = initial_volume - (n_nestedRVE+1)*volume_toremove;

            if Domain_size(3)==1 && ~strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 3 uncropped')
                n_nestedRVE = -1;
                warning('RVE convergence options incompatible with a 2D image');
                return
            end
        end
    else
        % if strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 1 uncropped')
        %     nextvolume = Domain_size(1) * Domain_size(2)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(3)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100);
        % elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 2 uncropped')
        %     nextvolume = Domain_size(2) * Domain_size(1)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(3)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100);
        % elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 3 uncropped')
        %     nextvolume = Domain_size(3) * Domain_size(1)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(2)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100);
        % end
        % while nextvolume>minimum_volume
        %     n_nestedRVE = n_nestedRVE+1;
        %     volume_toremove = initial_volume-nextvolume;
        %     ratio = ((initial_volume - volume_toremove)/initial_volume)^(1/2);
        %     if strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 1 uncropped')
        %         xcrop(n_nestedRVE,1) = 1;
        %         xcrop(n_nestedRVE,2) = Domain_size(1);
        %         xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
        %         xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
        %         xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
        %         xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);
        %     elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 2 uncropped')
        %         xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
        %         xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);
        %         xcrop(n_nestedRVE,3) = 1;
        %         xcrop(n_nestedRVE,4) = Domain_size(2);
        %         xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
        %         xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);
        %     elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 3 uncropped')
        %         xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
        %         xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);
        %         xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
        %         xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
        %         xcrop(n_nestedRVE,5) = 1;
        %         xcrop(n_nestedRVE,6) = Domain_size(3);
        %     end
        %     if strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 1 uncropped')
        %         nextvolume = Domain_size(1) * Domain_size(2)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(3)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100);
        %     elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 2 uncropped')
        %         nextvolume = Domain_size(2) * Domain_size(1)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(3)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100);
        %     elseif strcmp(RVEparameters.RVEconvergence_thatis,'And keep direction 3 uncropped')
        %         nextvolume = Domain_size(3) * Domain_size(1)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(2)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100);
        %     end
        % end
    end

elseif strcmp(RVEparameters.RVEconvergence_Crop,'3 directions')
    if Domain_size(3)>1
        if strcmp(RVEparameters.RVEconvergence_Ateachstep,'remove x %vol of the initial volume ')
            volume_toremove = initial_volume*RVEparameters.RVEconvergence_Ateachstep_x/100;
            nextvolume = initial_volume-volume_toremove;
            while nextvolume>minimum_volume
                n_nestedRVE = n_nestedRVE+1;
                ratio = ((initial_volume - n_nestedRVE*volume_toremove)/initial_volume)^(1/3);
                xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);
                xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
                xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
                xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);

                next_ratio = ((initial_volume - (n_nestedRVE+1)*volume_toremove)/initial_volume)^(1/3);
                next_x0 = round(Domain_size(1)/2 - Domain_size(1)*next_ratio/2);
                next_x1 = round(Domain_size(1)/2 + Domain_size(1)*next_ratio/2);
                next_y0 = round(Domain_size(2)/2 - Domain_size(2)*next_ratio/2);
                next_y1 = round(Domain_size(2)/2 + Domain_size(2)*next_ratio/2);
                next_z0 = round(Domain_size(3)/2 - Domain_size(3)*next_ratio/2);
                next_z1 = round(Domain_size(3)/2 + Domain_size(3)*next_ratio/2);
                nextvolume = (next_x1-next_x0+1)*(next_y1-next_y0+1)*(next_z1-next_z0+1);
            end
        else
            % nextvolume = Domain_size(1)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(2)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(3)*(1-RVEparameters.RVEconvergence_Ateachstep_x/100);
            % while nextvolume>minimum_volume
            %     n_nestedRVE = n_nestedRVE+1;
            %     volume_toremove = initial_volume-nextvolume;
            %     ratio = ((initial_volume - volume_toremove)/initial_volume)^(1/3);
            %     xcrop(n_nestedRVE,1) = round(Domain_size(1)/2 - Domain_size(1)*ratio/2);
            %     xcrop(n_nestedRVE,2) = round(Domain_size(1)/2 + Domain_size(1)*ratio/2);
            %     xcrop(n_nestedRVE,3) = round(Domain_size(2)/2 - Domain_size(2)*ratio/2);
            %     xcrop(n_nestedRVE,4) = round(Domain_size(2)/2 + Domain_size(2)*ratio/2);
            %     xcrop(n_nestedRVE,5) = round(Domain_size(3)/2 - Domain_size(3)*ratio/2);
            %     xcrop(n_nestedRVE,6) = round(Domain_size(3)/2 + Domain_size(3)*ratio/2);
            %     nextvolume = Domain_size(1)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(2)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100) * Domain_size(3)*(1 - (n_nestedRVE+1)*RVEparameters.RVEconvergence_Ateachstep_x/100);
            % end
        end
    else
        n_nestedRVE = -1;
        warning('RVE convergence options incompatible with a 2D image');
        return
    end

end