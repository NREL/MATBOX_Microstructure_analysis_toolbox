function [region_] = mesh_regionnumerotation(node,elem,face)

%% REGION NUMEROTATION
region_numerotation=unique(elem(:,5));
number_region_vol=length(region_numerotation); % Number of region
number_element_perregion=zeros(number_region_vol,1); % Number of element per region
for k=1:1:number_region_vol
    number_element_perregion(k,1)=sum(elem(:,5)==region_numerotation(k));
end
fprintf ('Number of region %i \n',(number_region_vol));
disp('  region numerotation:');
disp(region_numerotation)
disp('  number of element per region:');
disp(number_element_perregion)

%% SORT BY REGION 

for k=1:1:number_region_vol
    [node_,elem_,face_,corresponding_index_,elem_commonindex,face_commonindex] = Function_NodeFaceElem_region(node,elem,face,region_numerotation(k));
    % Mean position
    mean_x=mean(node_(:,1));
    % Save in a tructure
    region_(k).node = node_;
    region_(k).elem = elem_;
    region_(k).face = face_;
    region_(k).mean_x = mean_x;
    region_(k).corresponding_index = corresponding_index_;
    region_(k).elem_commonindex = elem_commonindex;
    region_(k).face_commonindex = face_commonindex;
end

end

