function [M] = Function_assign_voxels_based_on_contact(M,idx_assign,background_code)

sz = size(M);
[I1,I2,I3] = ind2sub(sz,idx_assign);

for k=1:1:length(idx_assign) % Loop over all voxels to re-assign
    x=I1(k); y=I2(k); z=I3(k); % Coordinates of current voxel
    neighbors_voxels = [];
    if x>1
        neighbors_voxels = [neighbors_voxels; M(x-1,y,z)];
    end
    if x<sz(1)
        neighbors_voxels = [neighbors_voxels; M(x+1,y,z)];
    end
    if y>1
        neighbors_voxels = [neighbors_voxels; M(x,y-1,z)];
    end
    if y<sz(2)
        neighbors_voxels = [neighbors_voxels; M(x,y+1,z)];
    end
    if z>1
        neighbors_voxels = [neighbors_voxels; M(x,y,z-1)];
    end
    if z<sz(3)
        neighbors_voxels = [neighbors_voxels; M(x,y,z+1)];
    end    
    
    % Id not to be reassigned
    no_assignment_id = abs(M(idx_assign(k)));
    neighbors_voxels(neighbors_voxels==no_assignment_id)=[]; % Remove them
    neighbors_voxels(neighbors_voxels<0)=[]; % Remove them
    
    % Prioritize assign to a solid phase, not to the background
    neighbors_voxels(neighbors_voxels==background_code)=[]; % Remove them
    
    % Sort to find phase that shares the most facet with the voxel
    [C,~,ic] = unique(neighbors_voxels);
    % Reassign
    if ~isempty(C)
        counts = accumarray(ic,1);
        value_counts = [C, counts];
        value_counts = sortrows(value_counts,-2);
        M(idx_assign(k)) = value_counts(1,1);
    else
        M(idx_assign(k)) = background_code;
    end
    
end


end

