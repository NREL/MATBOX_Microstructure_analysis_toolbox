function [microstructure] = generate_menger_sponge(number_of_iteration,Cube_subdivision)
% [geometry] = generate_menger_sponge(number_of_iteration,Cube_subdivision)
% Inputs:
% - number_of_iteration: number of iteration of the Menger sponge
% - Cube_subdivision: subdivide each cube (up-scaling)

% Initialisation
microstructure = ones(3^(number_of_iteration)*Cube_subdivision,3^(number_of_iteration)*Cube_subdivision,3^(number_of_iteration)*Cube_subdivision);
Domain_size = size(microstructure);

% Iteration process
for current_iteration = 1:1:number_of_iteration
    % Set number of subcube along one domain dimension
    subcube_number = 3^(current_iteration-1);
    % Set size of the subcube
    subcube_size = Domain_size(1)/subcube_number;
    % Set size of the subsubcube
    subsubcube_size = subcube_size/3;
    % Loop over all subcubes
    for current_subcube_dimension_x=1:1:subcube_number
        for current_subcube_dimension_y=1:1:subcube_number
            for current_subcube_dimension_z=1:1:subcube_number
                % Coordinates of the subcube
                x_min = (current_subcube_dimension_x-1)*subcube_size+1;
                x_max = x_min+subcube_size-1;
                y_min = (current_subcube_dimension_y-1)*subcube_size+1;
                y_max = y_min+subcube_size-1;                
                z_min = (current_subcube_dimension_z-1)*subcube_size+1;
                z_max = z_min+subcube_size-1;
                % Remove voxels according to the Menger sponge geometry definition
                % 7 sub cubes muse be removed
                x_remove_start = x_min+subsubcube_size;
                x_remove_end = x_remove_start+subsubcube_size-1;
                y_remove_start = y_min+subsubcube_size;
                y_remove_end = y_remove_start+subsubcube_size-1;
                z_remove_start = z_min+subsubcube_size;
                z_remove_end = z_remove_start+subsubcube_size-1;                
                % Cube central 1
                microstructure(x_remove_start:x_remove_end,y_remove_start:y_remove_end,z_remove_start:z_remove_end)=0;
                
                % Cube boundary 1 (x+)
                microstructure(x_remove_end+1:x_max,y_remove_start:y_remove_end,z_remove_start:z_remove_end)=0;
                % Cube boundary 2 (x-)
                microstructure(x_min:x_remove_start-1,y_remove_start:y_remove_end,z_remove_start:z_remove_end)=0;            

                % Cube boundary 3 (y+)
                microstructure(x_remove_start:x_remove_end,y_remove_end+1:y_max,z_remove_start:z_remove_end)=0;
                % Cube boundary 4 (y-)
                microstructure(x_remove_start:x_remove_end,y_min:y_remove_start-1,z_remove_start:z_remove_end)=0;
                
                % Cube boundary 5 (z+)
                microstructure(x_remove_start:x_remove_end,y_remove_start:y_remove_end,z_remove_end+1:z_max)=0;
                % Cube boundary 6 (z-)
                microstructure(x_remove_start:x_remove_end,y_remove_start:y_remove_end,z_min:z_remove_start-1)=0;

            end
        end
    end
end
microstructure=uint8(microstructure); % Convert in 8 bits

end
