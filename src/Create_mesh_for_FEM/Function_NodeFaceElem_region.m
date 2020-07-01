function [node_region_coordinate,elem_region_newindex,face_region_newindex,corresponding_index,elem_region,face_region] = Function_NodeFaceElem_region(node,elem,face,region_number)

% % Step1: find index
% Find all cell for the selected region
number_different_region=length(region_number);
index_region=[];
for k=1:1:number_different_region
   index_region_tmp=find(elem(:,5)==region_number(k));
   index_region=[index_region; index_region_tmp];
end
%index_region=find(elem(:,5)==region_number);
elem_region=elem(index_region,:);
% Find node number
node_region_index=unique(elem_region(:,1:4));
% Find face region
face_region_tmp=ismember(face(:,1:3), node_region_index); % find face that have any region node
sum_= sum(face_region_tmp,2); 
index_=find(sum_==3); % All three node of the face belong to the region
face_region=face(index_,:);
% Node coordinates
node_region_coordinate=node(node_region_index,:);

% % Step2: Build new index, starting from 1
n_=max(node_region_index);
max_=length(node_region_index);
corresponding_index=zeros(n_,2);
tmp=zeros(n_,1);
corresponding_index(:,1)=1:1:n_;
tmp(node_region_index)=1;
for k=n_:-1:1
    if tmp(k)==1
        tmp(k)=max_;
        max_=max_-1;
    end
end
corresponding_index(:,2)=tmp;

% % Step3: re-ordering
[n_cellregion,~]=size(elem_region);
elem_region_newindex=elem_region;
for k=1:1:n_cellregion
    for kk=1:1:4
        elem_region_newindex(k,kk)=corresponding_index(elem_region_newindex(k,kk),2);
    end
end
[n_faceregion,~]=size(face_region);
face_region_newindex=face_region;
for k=1:1:n_faceregion
    for kk=1:1:3
        face_region_newindex(k,kk)=corresponding_index(face_region_newindex(k,kk),2);
    end
end

end

