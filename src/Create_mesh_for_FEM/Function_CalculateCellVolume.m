function [ cell_volume,cell_centroid ] = Function_CalculateCellVolume(node,elem)

% Number of cells
[n_cell,~]=size(elem);

% Initialisation cell volume
cell_volume=zeros(n_cell,1);
% Cell center of mass
cell_centroid=zeros(n_cell,3);

% Loop over cell
for k=1:1:n_cell
    % Verteces index
    index_point_A=elem(k,1);
    index_point_B=elem(k,2);
    index_point_C=elem(k,3);
    index_point_D=elem(k,4);
    % Verteces coordinates
    A = node(index_point_A,:);
    B = node(index_point_B,:);
    C = node(index_point_C,:);
    D = node(index_point_D,:);

    % Distance
    AD=A-D;
    BD=B-D;
    CD=C-D;
    
    % Cross product
    cross_product=cross(BD,CD);
    
    % Dot product
    dot_product = dot(AD,cross_product);
    
    % Volume
    cell_volume(k,1)=norm(dot_product)/6;
    
    % Location
    cell_centroid(k,:)=(A+B+C+D)/4;
end



end

