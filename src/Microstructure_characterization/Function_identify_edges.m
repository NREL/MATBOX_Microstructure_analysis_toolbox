function [index_border_phase,B1,B2,B3] = Function_identify_edges(binary_phase)

% % For the complementary phase
% Number of voxel
number_voxel_complementary = sum(sum(sum(binary_phase==0)));
% Find index for the phase
index_complementary_phase = find(binary_phase==0);
% Get all coordinates
Domain_size = size(binary_phase);
[~,number_of_dimension]=size(Domain_size);
if number_of_dimension==2 % Case 2D
    Domain_size=[Domain_size(1) Domain_size(2) 1];
end
[C1,C2,C3] = ind2sub(Domain_size,index_complementary_phase);

for current_voxel=1:1:number_voxel_complementary
    % Get back voxel coordinate
    x_=C1(current_voxel); y_=C2(current_voxel); z_=C3(current_voxel);
    
    % Check x-
    if x_>1
        if binary_phase(x_-1,y_,z_)==1
            binary_phase(x_,y_,z_)=2;
        end
    end
    % Check y-
    if y_>1
        if binary_phase(x_,y_-1,z_)==1
            binary_phase(x_,y_,z_)=2;
        end
    end
    % Check z-
    if z_>1
        if binary_phase(x_,y_,z_-1)==1
            binary_phase(x_,y_,z_)=2;
        end
    end
    % Check x+
    if x_<Domain_size(1)
        if binary_phase(x_+1,y_,z_)==1
            binary_phase(x_,y_,z_)=2;
        end
    end
    % Check y+
    if y_<Domain_size(2)
        if binary_phase(x_,y_+1,z_)==1
            binary_phase(x_,y_,z_)=2;
        end
    end
    % Check z+
    if z_<Domain_size(3)
        if binary_phase(x_,y_,z_+1)==1
            binary_phase(x_,y_,z_)=2;
        end
    end
    
end
% Find index for the border voxel
index_border_phase = find(binary_phase==2);
% Get all coordinates
[B1,B2,B3] = ind2sub(Domain_size,index_border_phase);

end

