function [M] = Function_correct_microstructure(M,opt)
% Inputs:
% opt.background_code
%
% opt.todo_clean_connection
% opt.choice_cleansolidconnection
%
% opt.todo_convert_to_uniquecluster
% opt.choice_pore_uniquecluster
% opt.choice_unionsolid_uniquecluster

phases = double(unique(M)); % all phases id
phases(phases==opt.background_code)=[]; % Remove background id from phases id
n_phase = length(phases); % Number of phase

verification = true; % Mostly for debugging
number_change = 1;
old_M = M;
while number_change>0
    % Clean connection
    if opt.todo_clean_connection
        if strcmp(opt.choice_cleansolidconnection, 'All solid phases in one step')
            BW = zeros(size(M));
            BW(M~=opt.background_code)=1;
            [BW] = Function_clean_voxelconnection_multiphase(BW,opt.background_code,false,verification); % Apply to the solid phase union
            % Reassign over solid phases
            tmp = zeros(size(M))-1;
            tmp(BW==0)=opt.background_code;
            for k_phase=1:1:n_phase
                cond_1 = BW==1;
                cond_2 = M==phases(k_phase);
                idx = find(cond_1+cond_2==2);
                tmp(idx)=phases(k_phase);
            end
            M = tmp; clear tmp;
            idx = find(M==-1);
            if ~isempty(idx)
                warning('Function_correct_microstructure/clean connection: some voxels could not be fixed. They are assigned to the background') % Not sure if case possible
                M(idx)=opt.background_code;
            end
        elseif strcmp(opt.choice_cleansolidconnection, 'Solid phase per solid phase')
            [M] = Function_clean_voxelconnection_multiphase(M,opt.background_code,false,verification);
        end
    end
    
    % Convert to unique cluster
    if opt.todo_convert_to_uniquecluster
        [M] = Function_convert_to_unique_cluster(M,opt.background_code,opt.choice_pore_uniquecluster,opt.choice_unionsolid_uniquecluster,verification);
    end
    
    % Check while condition
    if opt.todo_clean_connection && opt.todo_convert_to_uniquecluster
        number_change = sum(sum(sum( old_M~= M )))
        old_M=M;
    else
        number_change = 0; % No need to continue
    end
    
    
end




% 
% dd
% 
% 
% 
% % Domain size
% Domain_size=size(M);
% 
% % New code (it can't be zero, as zero if for the non-connected clusters)
% new_code=code+1;
% new_complementary_code=complementary_code+1;
% 
% if ispc
%     % Calculate Checkcum
%     % Allow to exit the loop if the cleaning process goes back and forth
%     Checksum_state0=[];
%     for current_=1:1:Domain_size(3)
%         current_checksum = DataHash(M(:,:,current_));
%         Checksum_state0=[Checksum_state0 current_checksum];
%     end
%     Checksum_list.state(1).hash=Checksum_state0;
% end
% 
% % Initialisation and loop until every phase has only one unique connected cluster
% unallocated_voxel=1e9;
% iteration=1;
% while unallocated_voxel>0
%     
%     %% Connectivity: take largest cluster for both phase
%     
%     % Complementary phase cluster(s)
%     binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     if iteration==1
%         binary_phase(M == complementary_code) = 1;
%     else
%         binary_phase(M == new_complementary_code) = 1;
%     end    
%     L_complementary = bwlabeln(binary_phase,6);
%     % Number of cluster
%     n_cluster =length(unique(L_complementary))-1;
%     % Number of voxel for each cluster
%     allcluster_complementaryphase=zeros(n_cluster,2);
%     for k=1:1:n_cluster
%         allcluster_complementaryphase(k,1)=k;
%         allcluster_complementaryphase(k,2)=sum(sum(sum(L_complementary==k)));
%     end 
%     % Sort cluster (from larger to smaller)
%      allcluster_complementaryphase = sortrows(allcluster_complementaryphase,-2);
%        
%     % Phase cluster(s)
%     binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     if iteration==1
%         binary_phase(M == code) = 1;
%     else
%         binary_phase(M == new_code) = 1;
%     end        
%     L_phase = bwlabeln(binary_phase,6);
%     % Number of cluster
%     n_cluster =length(unique(L_phase))-1;
%     % Number of voxel for each cluster
%     allcluster_phase=zeros(n_cluster,2);
%     for k=1:1:n_cluster
%         allcluster_phase(k,1)=k;
%         allcluster_phase(k,2)=sum(sum(sum(L_phase==k)));
%     end
%     % Sort cluster (from larger to smaller)
%     allcluster_phase = sortrows(allcluster_phase,-2);
%     
%     % New code attributed to the larger cluster
%     M=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     M(L_phase==allcluster_phase(1))=new_code;
%     M(L_complementary==allcluster_complementaryphase(1))=new_complementary_code;
%     
%     
%     %% Isolated clusters are converted to the phase that surrounds them
%     
%     % Update number of Unallocated voxels
%     unallocated_voxel = sum(sum(sum(M==0)));
%     
%     if unallocated_voxel>0
%         % Take larger cluster
%         binary_complementaryphase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%         binary_complementaryphase(M==new_complementary_code)=1;
%         binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%         binary_phase(M==new_code)=1;
%         % Calculate distance
%         distance_complementaryphase = bwdist(binary_complementaryphase);
%         distance_phase = bwdist(binary_phase);
%         % Closer to complementary phase
%         delta_distance=distance_complementaryphase-distance_phase;
%         % Allocate the non connected voxels
%         binary_1=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%         binary_1(M==0)=1;
%         binary_2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%         binary_2(delta_distance>=0)=1;
%         binary_3=binary_1.*binary_2;
%         M(binary_3==1)=new_code;
%         
%         binary_1=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%         binary_1(M==0)=1;
%         binary_2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%         binary_2(delta_distance<0)=1;
%         binary_3=binary_1.*binary_2;
%         M(binary_3==1)=new_complementary_code;
%     end
%     
%     % Updated cluster
%     binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_phase(M == new_complementary_code) = 1;
%     L_complementary = bwlabeln(binary_phase,6);
%     binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_phase(M == new_code) = 1;
%     L_phase = bwlabeln(binary_phase,6);
% 
%         
%     %% Vertece-vertece and line-line connections are removed
%     
%     if cleanconnection==1
%         M = Function_clean_microstructure_onephase(M,new_code,new_complementary_code);
%     end
%     
%     %% Check exit
%     
%     % Update number of Unallocated voxels
%     unallocated_voxel = sum(sum(sum(M==0)));
% 
%     if ispc
%         % Checksum
%         Checksum_currentstate=[];
%         for current_=1:1:Domain_size(3)
%             current_checksum = DataHash(M(:,:,current_));
%             Checksum_currentstate=[Checksum_currentstate current_checksum];
%         end
%         
%         % Compare with previous checksum
%         exit_loop_checksum=0;
%         for previous_state=1:1:iteration
%             Checksum_previousstate = Checksum_list.state(previous_state).hash;
%             check_checksum = isequal(Checksum_previousstate,Checksum_currentstate);
%             if check_checksum==1
%                 exit_loop_checksum=1;
%             end
%         end
%         % Save current checksum in the list
%         Checksum_list.state(iteration+1).hash=Checksum_currentstate;
%         % Exit loop due to checksum equality
%         if (exit_loop_checksum==1 && unallocated_voxel>0)
%             break
%         end
%     end
%     
%     % Next iteration number
%     iteration=iteration+1;
%     
% end
% 
% %% Treat last isolated voxels
% 
% % Complementary phase cluster(s)
% binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% binary_phase(M == new_complementary_code) = 1;
% L_complementary = bwlabeln(binary_phase,6);
% % Number of cluster
% n_cluster =length(unique(L_complementary))-1;
% % Number of voxel for each cluster
% allcluster_complementaryphase=zeros(n_cluster,2);
% for k=1:1:n_cluster
%     allcluster_complementaryphase(k,1)=k;
%     allcluster_complementaryphase(k,2)=sum(sum(sum(L_complementary==k)));
% end
% % Sort cluster (from larger to smaller)
% allcluster_complementaryphase = sortrows(allcluster_complementaryphase,-2);
% 
% % Phase cluster(s)
% binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% binary_phase(M == new_code) = 1;
% L_phase = bwlabeln(binary_phase,6);
% % Number of cluster
% n_cluster =length(unique(L_phase))-1;
% % Number of voxel for each cluster
% allcluster_phase=zeros(n_cluster,2);
% for k=1:1:n_cluster
%     allcluster_phase(k,1)=k;
%     allcluster_phase(k,2)=sum(sum(sum(L_phase==k)));
% end
% % Sort cluster (from larger to smaller)
% allcluster_phase = sortrows(allcluster_phase,-2);
% 
% % New code attributed to the larger cluster
% M=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% M(L_phase==allcluster_phase(1))=new_code;
% M(L_complementary==allcluster_complementaryphase(1))=new_complementary_code;
% 
% 
% % Isolated clusters are converted to the phase that surrounds them
% % Update number of Unallocated voxels
% unallocated_voxel = sum(sum(sum(M==0)));
% if unallocated_voxel>0
%     % Take larger cluster
%     binary_complementaryphase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_complementaryphase(M==new_complementary_code)=1;
%     binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_phase(M==new_code)=1;
%     % Calculate distance
%     distance_complementaryphase = bwdist(binary_complementaryphase);
%     distance_phase = bwdist(binary_phase);
%     % Closer to complementary
%     delta_distance=distance_complementaryphase-distance_phase;
%     % Allocate the non connected voxels
%     binary_1=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_1(M==0)=1;
%     binary_2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_2(delta_distance>=0)=1;
%     binary_3=binary_1.*binary_2;
%     M(binary_3==1)=new_code;
%     binary_1=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_1(M==0)=1;
%     binary_2=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
%     binary_2(delta_distance<0)=1;
%     binary_3=binary_1.*binary_2;
%     M(binary_3==1)=new_complementary_code;
% end
% 
% unallocated_voxel = sum(sum(sum(M==0)))
% 
% 
% % Updated cluster
% binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% binary_phase(M == new_complementary_code) = 1;
% L_complementary = bwlabeln(binary_phase,6);
% binary_phase=zeros(Domain_size(1),Domain_size(2),Domain_size(3));
% binary_phase(M == new_code) = 1;
% L_phase = bwlabeln(binary_phase,6);


end

