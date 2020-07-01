function [index_border_label,B1,B2,B3] = Function_identify_labelsedges(label_map, background)


labels = unique(label_map);
labels(labels==background)=[];
Domain_size = size(label_map);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end
edge_map = zeros(Domain_size); % Initialization
for k=1:1:length(labels)
    n = sum(sum(sum(label_map==labels(k)))); % Number of voxel
    idx = find(label_map==labels(k)); % Find index for the label
    [C1,C2,C3] = ind2sub(Domain_size,idx); % Get all coordinates
    for current_voxel=1:1:n
        % Get back voxel coordinate
        x_=C1(current_voxel); y_=C2(current_voxel); z_=C3(current_voxel);
        % Check x-
        if x_>1
            if label_map(x_-1,y_,z_)~=label_map(x_,y_,z_) && label_map(x_-1,y_,z_)~=background
                edge_map(x_,y_,z_)=1;
            end
        end
        % Check y-
        if y_>1
            if label_map(x_,y_-1,z_)~=label_map(x_,y_,z_) && label_map(x_,y_-1,z_)~=background
                edge_map(x_,y_,z_)=1;
            end
        end
        % Check z-
        if z_>1
            if label_map(x_,y_,z_-1)~=label_map(x_,y_,z_) && label_map(x_,y_,z_-1)~=background
                edge_map(x_,y_,z_)=1;
            end
        end
        % Check x+
        if x_<Domain_size(1)
            if label_map(x_+1,y_,z_)~=label_map(x_,y_,z_) && label_map(x_+1,y_,z_)~=background
                edge_map(x_,y_,z_)=1;
            end
        end
        % Check y+
        if y_<Domain_size(2)
            if label_map(x_,y_+1,z_)~=label_map(x_,y_,z_) && label_map(x_,y_+1,z_)~=background
                edge_map(x_,y_,z_)=1;
            end
        end
        % Check z+
        if z_<Domain_size(3)
            if label_map(x_,y_,z_+1)~=label_map(x_,y_,z_) && label_map(x_,y_,z_+1)~=background
                edge_map(x_,y_,z_)=1;
            end
        end
    end
end

% Find index for the border voxel
index_border_label = find(edge_map==1);
% Get all coordinates
[B1,B2,B3] = ind2sub(Domain_size,index_border_label);

end

