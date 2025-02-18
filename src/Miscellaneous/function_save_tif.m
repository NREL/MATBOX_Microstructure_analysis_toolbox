function [] = function_save_tif(array3D, pathname)

sz=size(array3D);

array3D(array3D<0)=0;
[array3D] = fct_intconvert(array3D);

datatype = whos('array3D');

if strcmp(datatype.class,'uint32') || strcmp(datatype.class,'int32')
    % imwrite does not support 32 bits. We will have to use Tiff

    % Create a Tiff object and set up tags
    t = Tiff(pathname, 'w');
    tagstruct.ImageLength = size(array3D, 1);
    tagstruct.ImageWidth = size(array3D, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);

    % Write the data and close the Tiff object
    t.write(array3D(:,:,1));
    t.close();

    if length(sz)==3
        % Append next image
        for ii = 2:1:sz(3)
            t = Tiff(pathname,'a');
            tagstruct.ImageLength = size(array3D, 1);
            tagstruct.ImageWidth = size(array3D, 2);
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample = 32;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.RowsPerStrip = 16;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.Software = 'MATLAB';
            t.setTag(tagstruct);
            t.write(array3D(:,:,ii));
            t.close();
        end
    end

else

    % Initialise
    img = array3D(:,:,1);
    imwrite(img, pathname) % if double, convert to uint8
    if length(sz)==3
        % Append next image
        for ii = 2:1:sz(3)
            img = array3D(:,:,ii);
            imwrite(img, pathname, 'writemode', 'append');
        end
    end

end

end