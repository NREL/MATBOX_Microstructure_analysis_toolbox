function [Mb] = Function_clean_voxelconnection_multiphase(M,background_code,onlycheck,verification)

sz = size(M);
Ma = zeros(sz);
phases = unique(M); % all phases id
phases(phases==background_code)=[]; % Remove background id from phases id
n_phase = length(phases); % Number of phase

%% IDENTIFY VERTECE-VERTECE AND EDGE-EDGE CONNECTIONS
for k=1:1:n_phase % Loop over all phase
    % % Select subdomain
    phase = phases(k);
    idx = find(M==phase);
    [I1,I2,I3] = ind2sub(sz,idx);
    x_min = min(I1); x_max = max(I1);
    y_min = min(I2); y_max = max(I2);
    z_min = min(I3); z_max = max(I3);
    subdomain = M(x_min:x_max,y_min:y_max,z_min:z_max);
    BW=zeros(size(subdomain));
    BW(subdomain==phase)=1;
          
    % % Identify vertece-vertece and edge-edge connections
    [BW,index_removed_voxels] = Function_clean_voxelconnection(BW);
    
    subdomain(index_removed_voxels)=-phase; % Mark the removed voxels
    Ma(x_min:x_max,y_min:y_max,z_min:z_max)=subdomain;
        
end
idx_illconnected = find(Ma<0);
number_illconnected_voxels = length(idx_illconnected)
percent_illconnected_voxels = number_illconnected_voxels/sum(sum(sum(M~=background_code))) * 100;


%% CORRECT CONNECTIONS
if onlycheck
    Mb=M;
else
    if isempty(idx_illconnected)
        Mb=Ma;
    else
        [Mb] = Function_assign_voxels_based_on_contact(Ma,idx_illconnected,background_code);
    end
end

if verification
    for k=1:1:n_phase % Loop over all phase
        % % Select subdomain
        phase = phases(k);
        idx = find(Mb==phase);
        [I1,I2,I3] = ind2sub(sz,idx);
        x_min = min(I1); x_max = max(I1);
        y_min = min(I2); y_max = max(I2);
        z_min = min(I3); z_max = max(I3);
        subdomain = Mb(x_min:x_max,y_min:y_max,z_min:z_max);
        BW=zeros(size(subdomain));
        BW(subdomain==phase)=1;
        % % Identify vertece-vertece and edge-edge connections
        [~,index_removed_voxels] = Function_clean_voxelconnection(BW);
        if ~isempty(index_removed_voxels)
            length(index_removed_voxels)
            warning('Function_clean_voxelconnection_multiphase: some voxel are still incorrect')
        end
    end
end

end

