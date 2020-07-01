function [verteces_belong_to_face,Normal_adjacentfaces_maxdelta] = Funcion_CalculateFacetNormalVariation(node,face,Normal_faces)

% Number of faces
[number_face,~]=size(face);

% Number of verteces
[number_verteces,~]=size(node);

% Initialisation
Normal_adjacentfaces_maxdelta=zeros(number_face,1);
verteces_belong_to_face=zeros(number_verteces,100);
count_verteces=zeros(number_verteces+1,1);

%% STEP 1: identify belonging of each verteces

% Loop over all faces
for k=1:1:number_face
    % Verteces index
    verteceA=face(k,1);
    verteceB=face(k,2);
    verteceC=face(k,3);
    % Update count
    count_verteces(verteceA+1,1)=count_verteces(verteceA+1,1)+1;
    count_verteces(verteceB+1,1)=count_verteces(verteceB+1,1)+1;
    count_verteces(verteceC+1,1)=count_verteces(verteceC+1,1)+1;
    % Update vertece belonging
    verteces_belong_to_face(verteceA+1,count_verteces(verteceA+1,1))=k;
    verteces_belong_to_face(verteceB+1,count_verteces(verteceB+1,1))=k;
    verteces_belong_to_face(verteceC+1,count_verteces(verteceC+1,1))=k;
end

%% STEP 2: Calculate normal difference

% Loop over all faces
for k=1:1:number_face
    % Verteces index
    verteceA=face(k,1);
    verteceB=face(k,2);
    verteceC=face(k,3);
    % Vertces belonging
    vertce_belonging=[verteces_belong_to_face(verteceA+1,:) verteces_belong_to_face(verteceB+1,:) verteces_belong_to_face(verteceC+1,:)];
    % Remove 0
    vertce_belonging(vertce_belonging==0)=[];
    % Remove current face
    vertce_belonging(vertce_belonging==k)=[];
    
    % Current normal
    normal_currentface=Normal_faces(k,:);
    
    %Loop over all adjacent facets
    n_=length(vertce_belonging);
    Normal_adjacentfaces_maxdelta(k,1)=0;
    for kk=1:1:n_
        % Adjacent normal
        normal_adjacent=Normal_faces(vertce_belonging(kk),:);
        
        % Calculation of the angle
        % https://www.mathworks.com/matlabcentral/answers/90668-how-to-calculate-a-angle-between-two-vectors-in-3d
        %angle_originA = rad2deg(atan2(norm(cross(AB,AC)),dot(AB,AC)));
        % https://www.cs.berkeley.edu/~wkahan/Mindless.pdf), section 12 "Mangled Angles."
        angle_= rad2deg( 2*atan2(  norm( normal_currentface*norm(normal_adjacent) - norm(normal_currentface)*normal_adjacent) , norm( normal_currentface*norm(normal_adjacent) + norm(normal_currentface)*normal_adjacent) ) );
        % Keep only the maximum angle
        Normal_adjacentfaces_maxdelta(k,1)=max(Normal_adjacentfaces_maxdelta(k,1),angle_);
        
    end
    
    
    
end

% % OLD
% % Loop over all faces
% for k=1:1:number_face
%     
%     % Current normal
%     normal_currentface=Normal_faces(k,:);
%     
%     % First point
%     current_face_pointA=ones(number_face,1).*face(k,1);
%     % Difference - take only the min
%     difference=min(abs(current_face_pointA-face),[],2);
%     % Find the 0
%     idx_A=find(difference==0);
%     
%     % Second point
%     current_face_pointB=ones(number_face,1).*face(k,2);
%     % Difference - take only the min
%     difference=min(abs(current_face_pointB-face),[],2);
%     % Find the 0
%     idx_B=find(difference==0);  
% 
%     % Third point
%     current_face_pointC=ones(number_face,1).*face(k,3);
%     % Difference - take only the min
%     difference=min(abs(current_face_pointC-face),[],2);
%     % Find the 0
%     idx_C=find(difference==0);
%     
%     % Total index
%     idx_ABC=[idx_A;idx_B;idx_C];
%     % Remove k
%     idx_ABC(idx_ABC==k)=[];
%     
%     % Loop over all adjacent facets
%     n_=length(idx_ABC);
%     Normal_adjacentfaces_maxdelta(k,1)=0;
%     for kk=1:1:n_
%      % Adjacent normal
%      normal_adjacent=Normal_faces(idx_ABC(kk),:);
%  
%      % Calculation of the angle
%      % https://www.mathworks.com/matlabcentral/answers/90668-how-to-calculate-a-angle-between-two-vectors-in-3d
%      %angle_originA = rad2deg(atan2(norm(cross(AB,AC)),dot(AB,AC)));
%      % https://www.cs.berkeley.edu/~wkahan/Mindless.pdf), section 12 "Mangled Angles."
%      angle_= rad2deg( 2*atan2(  norm( normal_currentface*norm(normal_adjacent) - norm(normal_currentface)*normal_adjacent) , norm( normal_currentface*norm(normal_adjacent) + norm(normal_currentface)*normal_adjacent) ) );
%      % Keep only the maximum angle
%      Normal_adjacentfaces_maxdelta(k,1)=max(Normal_adjacentfaces_maxdelta(k,1),angle_);
%      
%     end
%     
% end

end