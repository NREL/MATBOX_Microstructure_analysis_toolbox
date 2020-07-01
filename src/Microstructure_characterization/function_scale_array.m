function [rescaled_array] = function_scale_array(array, initial_step, new_step, phaseinfo)
% Custom rescale of a 3D array based upon imresize3 and bwdist MATLAB bult-in function
% Works with n-phase array: no artifical intermediate phase is generated
% Preserve shape: no staircase visible when upscaling

% Using direcly imresize3(array,scaling_factor,'linear')
%    may create artifacts (artifical intermediate phase) when the number of phase is >2.
% Using direcly imresize3(array,scaling_factor,'nearest')
%    will not refine the interface at high scaling factor, making apparent the voxel initial size

scaling_factor = initial_step/new_step;
domain_size = size(array);
[number_phase,~] = size(phaseinfo);

if scaling_factor~=1
    % Create a binary volume
    binary_volume = zeros(domain_size);
    binary_volume(array~=0)=1; % 0 is the background
    % Upscale or downscale it
    allphases=imresize3(binary_volume,scaling_factor,'linear');
    idx_void = allphases<=0.5;
    allphases(idx_void)=0;
    allphases(allphases>0.5)=-1; % The shape we want to reach, with the phase information
   
    for current_phase=2:1:number_phase % Do not pass over background
        % Create a binary volume
        binary_volume = zeros(domain_size);
        binary_volume(array==phaseinfo(current_phase,1))=1;
        % Upscale or downscale it
        phase_resized=imresize3(binary_volume,scaling_factor,'linear');
        idx_phase = phase_resized>0.5;
        phase_resized(idx_phase)=1;
        phase_resized(phase_resized~=1)=0;
        allphases(logical(phase_resized))=phaseinfo(current_phase,1);
    end
    
    % -1 values are the voxels that refine the interface, but they still miss the phase information
    assigned_resizedvolume = zeros(size(allphases));
    assigned_resizedvolume(allphases>0)=1;
    [~,Idx_distance_map] = bwdist(assigned_resizedvolume);
    unassigned_idx = find(allphases==-1);
    Idx_distance_map(unassigned_idx);
    allphases(unassigned_idx) = allphases(Idx_distance_map(unassigned_idx));
    
    % And finally for the pore
    allphases(idx_void)=0;
    rescaled_array = allphases;
    clear allphases;
end

end

