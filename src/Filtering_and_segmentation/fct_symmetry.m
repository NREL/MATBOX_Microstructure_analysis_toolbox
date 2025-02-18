function [M] = fct_symmetry(M,p)

%% Check dimension compatibility first
sz = size(M);
dimension = length(sz);
sz_background = size(p.Background);
dimension_background = length(sz_background);
if dimension~=dimension_background
    warning('Dimension incompatiblity in fct_symmetry. Exit function.')
    return
end
if sum(sz==sz_background)~=dimension
    warning('Dimension incompatiblity in fct_symmetry. Exit function.')
    return
end

%% Algorithm
if dimension == 2
    sz = [sz 1];
end

keep_iterate = true;
while keep_iterate
    keep_iterate = false;

    if dimension==2 || (dimension == 3 && p.perslice)

        for z=1:1:sz(3) % Loop over slices
            slice_M = M(:,:,z);
            slice_background = p.Background(:,:,z);

            [dmap_background, idx_back] = bwdist(~slice_background);
            % figure; imagesc(dmap_background); colormap turbo; axis equal; axis tight;


            idx = find(slice_background);
            n = length(idx);
            [IX, IY] = ind2sub(sz,idx);
            for k=1:1:n % Loop over all pixel within background
                background_x = IX(k);
                background_y = IY(k);
                %background_dist = dmap_background(background_x,background_y);

                % Find nearest interface pixel from the current background pixel
                interface_idx = idx_back(background_x,background_y);
                [interface_x, interface_y] = ind2sub(sz,interface_idx);

                % Find symmetric point and replace
                material_x = interface_x + interface_x-background_x;
                material_y = interface_y + interface_y-background_y;

                if material_x ==0
                    material_x=1;
                elseif material_x == sz(1)+1
                    material_x = sz(1);
                end
                if material_y ==0
                    material_y=1;
                elseif material_y == sz(2)+1
                    material_y = sz(2);
                end
                if material_x<1 || material_y<1 || material_x>sz(1) || material_y>sz(2)
                    keep_iterate = true; % We will apply symmetry on the result of the previous symmetry
                else
                    if slice_background(material_x,material_y)~=1
                        % Overwritte background value
                        slice_M(background_x,background_y) = slice_M(material_x, material_y);
                        slice_background(background_x,background_y) = 0;
                    else
                        keep_iterate = true; % We will apply symmetry on the result of the previous symmetry
                    end
                end

                % Alternative way to find the symmetric point (slower)
                % % Find all pixels within material domain at same distance from interface
                % [dmap_material, ~] = bwdist(slice_background);
                % % figure; imagesc(dmap_material); colormap turbo; axis equal; axis tight;
                % tmp = find(abs(dmap_material-background_dist)<0.5);
                % if ~isempty(tmp)
                %     %dmap_material_currentdist = zeros(sz);
                %     %dmap_material_currentdist(tmp)=1;
                %     % figure; imagesc(dmap_material_currentdist); colormap grey; axis equal; axis tight;
                %     [material_x, material_y] = ind2sub(sz,tmp);
                %
                %     % Calculate all distances from these material points to the interface point
                %     all_dist = ((material_x-interface_x).^2 + (material_y-interface_y).^2).^0.5;
                %
                %     % Select pixel from these material points nearest to the interface point
                %     delta_dist = abs(all_dist-background_dist);
                %     id_material = find( delta_dist==min(delta_dist) );
                %     id_material = id_material(1);
                %
                %     % Overwritte background value
                %     slice_M(background_x,background_y) = slice_M(material_x(id_material), material_y(id_material));
                %     slice_background(background_x,background_y) = 0;
                % else
                %     keep_iterate = true; % We will apply symmetry on the result of the previous symmetry
                % end

            end
            % Re-insert
            M(:,:,z) = slice_M;
            p.Background(:,:,z) = slice_background;
        end


    else

        slice_M = M;
        slice_background = p.Background;

        [dmap_background, idx_back] = bwdist(~slice_background);
        % figure; imagesc(dmap_background); colormap turbo; axis equal; axis tight;


        idx = find(slice_background);
        n = length(idx);
        [IX, IY, IZ] = ind2sub(sz,idx);
        for k=1:1:n % Loop over all pixel within background
            background_x = IX(k);
            background_y = IY(k);
            background_z = IZ(k);
            %background_dist = dmap_background(background_x,background_y,background_z);

            % Find nearest interface pixel from the current background pixel
            interface_idx = idx_back(background_x,background_y,background_z);
            [interface_x, interface_y, interface_z] = ind2sub(sz,interface_idx);

            % Find symmetric point and replace
            material_x = interface_x + interface_x-background_x;
            material_y = interface_y + interface_y-background_y;
            material_z = interface_z + interface_z-background_z;

            if material_x ==0
                material_x=1;
            elseif material_x == sz(1)+1
                material_x = sz(1);
            end
            if material_y ==0
                material_y=1;
            elseif material_y == sz(2)+1
                material_y = sz(2);
            end
            if material_z ==0
                material_z=1;
            elseif material_z == sz(3)+1
                material_z = sz(3);
            end            
            if material_x<1 || material_y<1 || material_z<1 || material_x>sz(1) || material_y>sz(2) || material_z>sz(3)
                keep_iterate = true; % We will apply symmetry on the result of the previous symmetry
            else
                if slice_background(material_x,material_y,material_z)~=1
                    % Overwritte background value
                    slice_M(background_x,background_y,background_z) = slice_M(material_x, material_y,material_z);
                    slice_background(background_x,background_y,background_z) = 0;
                else
                    keep_iterate = true; % We will apply symmetry on the result of the previous symmetry
                end
            end

            % Alternative way to find the symmetric point (slower)
            % % Find all pixels within material domain at same distance from interface
            % [dmap_material, ~] = bwdist(slice_background);
            % % figure; imagesc(dmap_material); colormap turbo; axis equal; axis tight;
            % tmp = find(abs(dmap_material-background_dist)<0.5);
            % if ~isempty(tmp)
            %     %dmap_material_currentdist = zeros(sz);
            %     %dmap_material_currentdist(tmp)=1;
            %     % figure; imagesc(dmap_material_currentdist); colormap grey; axis equal; axis tight;
            %     [material_x, material_y] = ind2sub(sz,tmp);
            %
            %     % Calculate all distances from these material points to the interface point
            %     all_dist = ((material_x-interface_x).^2 + (material_y-interface_y).^2).^0.5;
            %
            %     % Select pixel from these material points nearest to the interface point
            %     delta_dist = abs(all_dist-background_dist);
            %     id_material = find( delta_dist==min(delta_dist) );
            %     id_material = id_material(1);
            %
            %     % Overwritte background value
            %     slice_M(background_x,background_y) = slice_M(material_x(id_material), material_y(id_material));
            %     slice_background(background_x,background_y) = 0;
            % else
            %     keep_iterate = true; % We will apply symmetry on the result of the previous symmetry
            % end

        end
        % Re-insert
        M = slice_M;
        p.Background = slice_background;



    end

    if sum(sum(sum(p.Background)))==0
        keep_iterate = false;
    end

end

% Convert to uint8 or uint16
[M] = fct_intconvert(M);

end
