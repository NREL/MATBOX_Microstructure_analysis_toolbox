function [M] = Function_remove_isolated_clusters(M, background_code, sizefilter)
% Isolated voxels are assigned to adjacent phase

phases = double(unique(M)); % all phases id
n_phase = length(phases); % Number of phase
conn=6;
for k_phase=1:1:n_phase
    idx_assign=[];
    BW=zeros(size(M));
    BW(M==phases(k_phase))=1;
    [L] = bwlabeln(BW,conn);
    [C,~,ic] = unique(L);
    counts = accumarray(ic,1);
    value_counts = [C, counts];
    % Remove background
    idx=find(value_counts(:,1)==0);
    value_counts(idx,:)=[];
    % Find id for which number of voxel is below or equal with sizefilter
    id1 = find( value_counts(:,2) <=sizefilter );
    if ~isempty(id1)
        id2 = value_counts(id1,1);
        for k=1:1:length(id2)
            idx_assign = [idx_assign; find( L==id2(k) )];
        end
    end   
    [M] = Function_assign_voxels_based_on_contact(M,idx_assign,background_code);
end

end

