function [M] = Function_convert_to_unique_cluster(M,background_code, choice_pore_uniquecluster,choice_unionsolid_uniquecluster,verification)

phases = double(unique(M)); % all phases id
phases(phases==background_code)=[]; % Remove background id from phases id
n_phase = length(phases); % Number of solid phase
conn=6; % Connectivity face-to-face only

%% PORE CONVERTED IN A UNIQUE CONNECTED CLUSTER
if choice_pore_uniquecluster
    % Select largest pore cluster
    BW=zeros(size(M));
    BW(M==background_code)=1;
    [L] = bwlabeln(BW,conn);
    [C,~,ic] = unique(L);
    counts = accumarray(ic,1);
    value_counts = [C, counts];
    % Remove background
    idx=find(value_counts(:,1)==0);
    value_counts(idx,:)=[];
    value_counts = sortrows(value_counts,-2);
    %Largest_porecluster = zeros(size(M));
    %Largest_porecluster(L==value_counts(1,1))=1; 
    
    % Remaining cluster will be assigned to nearest solid phase
    cond_1 = M==background_code;
    cond_2 = L~=value_counts(1,1);
    idx_assign = find(cond_1+cond_2==2);
    [M] = Function_assign_voxels_based_on_contact(M,idx_assign,background_code);
    
    if verification
        BW=zeros(size(M));
        BW(M==background_code)=1;
        [L] = bwlabeln(BW,conn);
        [C] = unique(L);
        if length(C)>2
            warning('Function_convert_to_unique_cluster: ''pore converted in a unique connected cluster'' operation is incorrect')
        end
    end
    
end

%% UNION OF SOLID PHASE CONVERTED IN A UNIQUE CONNECTED CLUSTER
if choice_unionsolid_uniquecluster
    % Select largest pore cluster
    BW=zeros(size(M));
    BW(M~=background_code)=1;
    [L] = bwlabeln(BW,conn);
    [C,~,ic] = unique(L);
    counts = accumarray(ic,1);
    value_counts = [C, counts];
    % Remove background
    idx=find(value_counts(:,1)==0);
    value_counts(idx,:)=[];
    value_counts = sortrows(value_counts,-2);
    Largest_solidcluster = zeros(size(M));
    Largest_solidcluster(L==value_counts(1,1))=1;
    
    % Outside from the largest solid cluster, assign to pore
    tmp = zeros(size(M))-1;
    tmp(Largest_solidcluster==0)=background_code;
    
    % Within this largest solid cluster, assign the different solid phases
    for k_phase=1:1:n_phase
        cond_1 = Largest_solidcluster==1;
        cond_2 = M==phases(k_phase);
        idx = find(cond_1+cond_2==2);
        tmp(idx)=phases(k_phase);
    end
    M = tmp; clear tmp;
    idx = find(M==-1);
    if ~isempty(idx)
        warning('Function_convert_to_unique_cluster: some voxels could not be fixed. They are assigned to the background') % Not sure if case possible
        M(idx)=opt.background_code;
    end
    
    if verification
        BW=zeros(size(M));
        BW(M~=background_code)=1;
        [L] = bwlabeln(BW,conn);
        [C] = unique(L);
        if length(C)>2
            warning('Function_convert_to_unique_cluster: ''pore converted in a unique connected cluster'' operation is incorrect')
        end
    end    
end

end
