function [node_subregion,elem_subregion,face_subregion, index_subregion_node, index_subregion_elem, index_subregion_face] = Function_subregion(node_region,elem_region,face_region,tolerance)

% Node coordinates
coor_x = node_region(:,1);
coor_y = node_region(:,2);
coor_z = node_region(:,3);

% Extremums
x_min = min(coor_x); x_max = max(coor_x); 
y_min = min(coor_y); y_max = max(coor_y); 
z_min = min(coor_z); z_max = max(coor_z); 

% Find nodes 
index_x_left = find(coor_x-x_min >= tolerance(1,1));
index_x_right = find(x_max-coor_x >= tolerance(1,2));
index_y_left = find(coor_y-y_min >= tolerance(2,1));
index_y_right = find(y_max-coor_y >= tolerance(2,2));
index_z_left = find(coor_z-z_min >= tolerance(3,1));
index_z_right = find(z_max-coor_z >= tolerance(3,2));
index_subregion_x = intersect(index_x_left,index_x_right);
index_subregion_y = intersect(index_y_left,index_y_right);
index_subregion_z = intersect(index_z_left,index_z_right);
index_subregion_node = intersect(index_subregion_z, intersect(index_subregion_x,index_subregion_y));
node_subregion = node_region(index_subregion_node,:);

% Find faces
tmp = ismember(face_region,index_subregion_node);
tmp = sum(tmp,2);
index_subregion_face = find(tmp==3);
face_subregion = face_region(index_subregion_face,:);

% Find elem
tmp = ismember(elem_region(:,1:4),index_subregion_node);
tmp = sum(tmp,2);
index_subregion_elem = find(tmp==4);
elem_subregion = elem_region(index_subregion_elem,:);

end

