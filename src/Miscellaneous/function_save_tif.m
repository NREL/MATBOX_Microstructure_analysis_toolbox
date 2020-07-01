function [] = function_save_tif(array3D, pathname)

Domain_size=size(array3D);

datatype = whos('array3D');
if strcmp(datatype.class,'double') || strcmp(datatype.class,'single')
    array3D=im2uint16(array3D); % Prevent imwrite to convert in uint8
end

% Initialise
img = array3D(:,:,1);
imwrite(img, pathname) % if double, convert to uint8
if length(Domain_size)==3
    % Append next image
    for ii = 2:1:Domain_size(3)
        img = array3D(:,:,ii);
        imwrite(img, pathname, 'writemode', 'append');
    end
end

end