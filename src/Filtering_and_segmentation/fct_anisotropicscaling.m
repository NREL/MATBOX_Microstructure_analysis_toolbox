function [M,newtype,foo] = fct_anisotropicscaling(M,p)

sz = size(M);
foo=[];
newtype = 'same';
dimension = length(sz);

if any(p.scaling_factors~=1)
    %p.scaling_factors=1/p.scaling_factors; % I have to harmonize this will all the different scaling in this toolbox...

    if dimension == 2
        A22 = p.scaling_factors(1); % not an error
        A11= p.scaling_factors(2); % not an error
        % Set parameter
        A = [A11 0 0;
            0 A22 0;
            0 0 1];
        tform = affinetform2d(A);
        if strcmp(p.datatype,'Grey level') || strcmp(p.datatype,'Channel')
            M = imwarp(M,tform,'linear');
        elseif app.Anisotropy_GreylevelButton.Value
            M = imwarp(M,tform,'nearest');
        end

        % Padding
        if A11~=0
            M(1,:) = M(2,:);
            M(end,:) = M(end-1,:);
        end
        if A22~=0
            M(:,1) = M(:,2);
            M(:,end) = M(:,end-1);
        end           

    else
        A22 = p.scaling_factors(1); % not an error
        A11= p.scaling_factors(2); % not an error
        A33 = p.scaling_factors(3);
        % Set parameter
        A = [A11 0 0 0;
            0 A22 0 0;
            0 0 A33 0;
            0 0 0 1];
        tform = affinetform3d(A);
        if strcmp(p.datatype,'Grey level') || strcmp(p.datatype,'Channel')
            M = imwarp(M,tform,'linear');
        elseif app.Anisotropy_GreylevelButton.Value
            M = imwarp(M,tform,'nearest');
        end
    
        % Padding
        if A11~=0
            M(1,:,:) = M(2,:,:);
            M(end,:,:) = M(end-1,:,:);
        end
        if A22~=0
            M(:,1,:) = M(:,2,:);
            M(:,end,:) = M(:,end-1,:);
        end 
        if A33~=0
            M(:,:,1) = M(:,:,2);
            M(:,:,end) = M(:,:,end-1);
        end
    end

    % Convert to uint8 or uint16
    if min(min(min(M)))>=0
        if max(max(max(M)))<=255
            M=uint8(M);
        elseif max(max(max(M)))<=65535
            M=uint16(M);
        end
    else
        if max(max(max(abs(M))))<=255
            M=int8(M);
        elseif max(max(max(abs(M))))<=65535
            M=int16(M);
        end
    end

end
