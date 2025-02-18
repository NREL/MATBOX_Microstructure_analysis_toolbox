function [Position_line, TPBcounter, LOC_tpbl] = Function_TPBL_direct_Algorithm(binary_phase, labels)

label_A = labels(1);
label_B = labels(2);
label_C = labels(3);

Domain_size=size(binary_phase); % Domain size of the microstructure

% % Position of the line, Initialization
% % One step is equal to half voxel length
% Position_line.direction(1).value=zeros(Domain_size(1),2);
% Position_line.direction(2).value=zeros(Domain_size(2),2);
% Position_line.direction(3).value=zeros(Domain_size(3),2);
% % Position (axis)
% Position_line.direction(1).value(:,1)=(1:1:Domain_size(1))' - 0.5;
% Position_line.direction(2).value(:,1)=(1:1:Domain_size(2))' - 0.5;
% Position_line.direction(3).value(:,1)=(1:1:Domain_size(3))' - 0.5;

% To keep indexing coherent during reshape
Allposition_1 = zeros(size(binary_phase));
for k=1:1:Domain_size(1)
    Allposition_1(k,:,:) = k;
end
Allposition_2 = zeros(size(binary_phase));
for k=1:1:Domain_size(2)
    Allposition_2(:,k,:) = k;
end
Allposition_3 = zeros(size(binary_phase));
for k=1:1:Domain_size(3)
    Allposition_3(:,:,k) = k;
end

% Initialize position
pos_1 = [];
pos_2 = [];
pos_3 = [];
%  Initialize number of TPB
TPBcounter = 0;

% Find all contact line parallel with axe 3
Squarelist = zeros( (Domain_size(1)-1) * (Domain_size(2)-1), 4);  
Squarelist_1 = zeros( (Domain_size(1)-1) * (Domain_size(2)-1), 4);  
Squarelist_2 = zeros( (Domain_size(1)-1) * (Domain_size(2)-1), 4);  
for k=1:1:Domain_size(3)
    
    slice_tmp = binary_phase(:,:,k);
    Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    tmp = logical(sum(Squarelist==label_A,2)) + logical(sum(Squarelist==label_B,2)) + logical(sum(Squarelist==label_C,2));
    idx = find(tmp==3);
    n = length(idx);
    TPBcounter = TPBcounter + n;

    pos_3 = [pos_3; ones(n,1)*k];

    slice_1 = Allposition_1(:,:,k);
    Squarelist_1(:,1) = reshape(slice_1(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist_1(:,2) = reshape(slice_1(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist_1(:,3) = reshape(slice_1(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist_1(:,4) = reshape(slice_1(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    pos_1 = [pos_1; round(mean(Squarelist_1(idx,:),2))];  

    slice_2 = Allposition_2(:,:,k);
    Squarelist_2(:,1) = reshape(slice_2(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist_2(:,2) = reshape(slice_2(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist_2(:,3) = reshape(slice_2(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    Squarelist_2(:,4) = reshape(slice_2(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    pos_2 = [pos_2; round(mean(Squarelist_2(idx,:),2))];  

    % % Alternative, but do not get position of ALL contact length (only 1/3)
    % slice_tmp = binary_phase(:,:,k);
    % Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % tmp = logical(sum(Squarelist==label_A,2)) + logical(sum(Squarelist==label_B,2)) + logical(sum(Squarelist==label_C,2));
    % TPBcounter = TPBcounter + sum(tmp==3); 

    % % Alternative, but assume binary phase has only 3 different labels (while there can be more)
    % slice_tmp = binary_phase(:,:,k);
    % Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);
    % Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(2)-1),1]);    
    % numberunique = sum(diff(sort(Squarelist,2),1,2)~=0,2)+1; % Number of unique value per rows (https://www.mathworks.com/matlabcentral/answers/427942-unique-values-per-row)
    % TPBcounter = TPBcounter + sum(numberunique==3);
end

% Find all contact line parallel with axe 1
Squarelist = zeros( (Domain_size(2)-1) * (Domain_size(3)-1), 4);  
Squarelist_2 = zeros( (Domain_size(2)-1) * (Domain_size(3)-1), 4);  
Squarelist_3 = zeros( (Domain_size(2)-1) * (Domain_size(3)-1), 4);  
for k=1:1:Domain_size(1)

    slice_tmp = squeeze(binary_phase(k,:,:));
    Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    tmp = logical(sum(Squarelist==label_A,2)) + logical(sum(Squarelist==label_B,2)) + logical(sum(Squarelist==label_C,2));
    idx = find(tmp==3);
    n = length(idx);
    TPBcounter = TPBcounter + n;

    pos_1 = [pos_1; ones(n,1)*k];

    slice_2 = squeeze(Allposition_2(k,:,:));
    Squarelist_2(:,1) = reshape(slice_2(2:end,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist_2(:,2) = reshape(slice_2(2:end,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist_2(:,3) = reshape(slice_2(1:end-1,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist_2(:,4) = reshape(slice_2(1:end-1,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    pos_2 = [pos_2; round(mean(Squarelist_2(idx,:),2))];  

    slice_3 = squeeze(Allposition_3(k,:,:));
    Squarelist_3(:,1) = reshape(slice_3(2:end,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist_3(:,2) = reshape(slice_3(2:end,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist_3(:,3) = reshape(slice_3(1:end-1,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    Squarelist_3(:,4) = reshape(slice_3(1:end-1,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    pos_3 = [pos_3; round(mean(Squarelist_3(idx,:),2))];      
    
    %slice_tmp = squeeze(binary_phase(k,:,:));
    %Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    %Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    %Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    %Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(2)-1) * (Domain_size(3)-1),1]);
    %tmp = logical(sum(Squarelist==label_A,2)) + logical(sum(Squarelist==label_B,2)) + logical(sum(Squarelist==label_C,2));
    %TPBcounter = TPBcounter + sum(tmp==3); 
    %numberunique = sum(diff(sort(Squarelist,2),1,2)~=0,2)+1;
    %TPBcounter = TPBcounter + sum(numberunique==3); 
end

% Find all contact line parallel with axe 2
Squarelist = zeros( (Domain_size(1)-1) * (Domain_size(3)-1), 4);  
Squarelist_1 = zeros( (Domain_size(1)-1) * (Domain_size(3)-1), 4);  
Squarelist_3 = zeros( (Domain_size(1)-1) * (Domain_size(3)-1), 4);  
for k=1:1:Domain_size(2)

    slice_tmp = squeeze(binary_phase(:,k,:));
    Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    tmp = logical(sum(Squarelist==label_A,2)) + logical(sum(Squarelist==label_B,2)) + logical(sum(Squarelist==label_C,2));
    idx = find(tmp==3);
    n = length(idx);
    TPBcounter = TPBcounter + n;

    pos_2 = [pos_2; ones(n,1)*k];

    slice_1 = squeeze(Allposition_1(:,k,:));
    Squarelist_1(:,1) = reshape(slice_1(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist_1(:,2) = reshape(slice_1(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist_1(:,3) = reshape(slice_1(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist_1(:,4) = reshape(slice_1(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);    
    pos_1 = [pos_1; round(mean(Squarelist_1(idx,:),2))]; 

    slice_3 = squeeze(Allposition_3(:,k,:));
    Squarelist_3(:,1) = reshape(slice_3(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist_3(:,2) = reshape(slice_3(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist_3(:,3) = reshape(slice_3(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    Squarelist_3(:,4) = reshape(slice_3(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);    
    pos_3 = [pos_3; round(mean(Squarelist_3(idx,:),2))];     

    %slice_tmp = squeeze(binary_phase(:,k,:));
    %Squarelist(:,1) = reshape(slice_tmp(2:end,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    %Squarelist(:,2) = reshape(slice_tmp(2:end,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    %Squarelist(:,3) = reshape(slice_tmp(1:end-1,1:end-1),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    %Squarelist(:,4) = reshape(slice_tmp(1:end-1,2:end),[(Domain_size(1)-1) * (Domain_size(3)-1),1]);
    %tmp = logical(sum(Squarelist==label_A,2)) + logical(sum(Squarelist==label_B,2)) + logical(sum(Squarelist==label_C,2));
    %TPBcounter = TPBcounter + sum(tmp==3); 
    %numberunique = sum(diff(sort(Squarelist,2),1,2)~=0,2)+1;
    %TPBcounter = TPBcounter + sum(numberunique==3); 
end

id = sub2ind(size(binary_phase),pos_1,pos_2,pos_3);
LOC_tpbl=zeros(size(binary_phase));
LOC_tpbl(id)=1;

[GC,GR] = groupcounts(pos_1);
tmp = zeros(Domain_size(1),1);
for k=1:1:Domain_size(1)
    id = find(GR==k);
    if ~isempty(id)
        tmp(k)=GC(id);
    else
        tmp(k)=0;
    end
end
Position_line.direction(1).value(:,1) = [1:1:Domain_size(1)]'; % GR
Position_line.direction(1).value(:,2) = tmp; % GC

[GC,GR] = groupcounts(pos_2);
tmp = zeros(Domain_size(2),1);
for k=1:1:Domain_size(2)
    id = find(GR==k);
    if ~isempty(id)
        tmp(k)=GC(id);
    else
        tmp(k)=0;
    end
end
Position_line.direction(2).value(:,1) = [1:1:Domain_size(2)]'; % GR
Position_line.direction(2).value(:,2) = tmp; % GC

[GC,GR] = groupcounts(pos_3);
tmp = zeros(Domain_size(3),1);
for k=1:1:Domain_size(3)
    id = find(GR==k);
    if ~isempty(id)
        tmp(k)=GC(id);
    else
        tmp(k)=0;
    end
end
Position_line.direction(3).value(:,1) = [1:1:Domain_size(3)]'; % GR
Position_line.direction(3).value(:,2) = tmp; % GC

end