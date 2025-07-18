function [Microstructure_resized] = function_scaling(Microstructure,p)
%function_scaling scales a 3D array
% [Microstructure_resized] = function_scaling(Microstructure,p)
% Inputs:
% - Microstructure: a 3D array
% - p: a structure
%   p.scaling_factor: new voxel size = current voxel size * p.scaling_factor
%                     A scaling factor value above 1 is for downscaling, while a scaling factor value below 1 is upscaling.
%   p.label_or_greylevel: 'Label' for segmented image, or 'Grey level' for  unsegmented image
%   p.background: value of the backgroun (or pore phase), only required for 'Label'   

if p.scaling_factor~=1
    p.scaling_factor=1/p.scaling_factor; 
    domain_size = size(Microstructure);
    dimension = length(domain_size);
    if strcmp(p.label_or_greylevel,'Label')

        if p.scaling_factor > 1
            % Using direcly imresize3(Microstructure,scaling_factor,'linear')
            % may create artifacts (artifical intermediate phase) when the number of phase is >2.
            % Using direcly imresize3(Microstructure,scaling_factor,'nearest')
            % will not refine the interface at high scaling factor, making apparent the voxel initial size

            % Create a binary volume
            binary_volume = zeros(domain_size);
            binary_volume(Microstructure~=p.background)=1;
            % Upscale or downscale it
            if dimension==2
                allphases=imresize(binary_volume,p.scaling_factor,'bilinear');
            else
                allphases=imresize3(binary_volume,p.scaling_factor,'linear');
            end
            idx_void = allphases<=0.5;
            allphases(idx_void)=0;
            allphases(allphases>0.5)=-1; % The shape we want to reach, with the phase information
            %allphases(allphases>0.5)=1;

            phases_id = unique(Microstructure);
            phases_id(phases_id == p.background)=[];
            number_phase = length(phases_id); % Minus background
            for current_phase=1:1:number_phase
                % Create a binary volume
                binary_volume = zeros(domain_size);
                binary_volume(Microstructure==phases_id(current_phase))=1;
                % Upscale or downscale it
                if dimension==2
                    phase_resized=imresize(binary_volume,p.scaling_factor,'bilinear');
                else
                    phase_resized=imresize3(binary_volume,p.scaling_factor,'linear');
                end
                idx_phase = phase_resized>0.5;
                phase_resized(idx_phase)=1;
                phase_resized(phase_resized~=1)=0;
                allphases(logical(phase_resized))=phases_id(current_phase);
            end

            % -1 values are the voxels that refine the interface, but they still miss the phase information
            unassigned_idx = find(allphases==-1);
            if ~isempty(unassigned_idx)
                assigned_resizedvolume = zeros(size(allphases));
                assigned_resizedvolume(allphases>0)=1;
                [~,Idx_distance_map] = bwdist(assigned_resizedvolume);
                %Idx_distance_map(unassigned_idx);
                allphases(unassigned_idx) = allphases(Idx_distance_map(unassigned_idx));
            end

            % And finally for the pore
            allphases(idx_void)=p.background;
            Microstructure_resized = allphases;
            clear allphases;

        else
            if dimension==2
                Microstructure_resized = imresize(Microstructure,p.scaling_factor,'nearest');
            else
                Microstructure_resized = imresize3(Microstructure,p.scaling_factor,'nearest');
            end
        end

    elseif strcmp(p.label_or_greylevel,'Grey level')
        if dimension==2
            Microstructure_resized = imresize(Microstructure,p.scaling_factor,'bilinear');
        else
            Microstructure_resized = imresize3(Microstructure,p.scaling_factor,'linear');
            %Microstructure_resized = imresize3(Microstructure,p.scaling_factor,'nearest');
        end
        
    end
    
    % Data type
    if strcmp(class(Microstructure),'uint8')
        Microstructure_resized = uint8(Microstructure_resized);
    elseif strcmp(class(Microstructure),'uint16')
        Microstructure_resized = uint16(Microstructure_resized);
    end
    
else
    Microstructure_resized=Microstructure;
end

end

