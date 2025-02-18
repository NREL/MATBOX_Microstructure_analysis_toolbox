function [sz, outcome] = function_checksize_tif(FileTif)

% Initialize
sz = [];
outcome = false;

if isfile(FileTif) % Check existence
    [~,~,ext] = fileparts(FileTif);
    
    if strcmp(ext,'.tif') || strcmp(ext,'.tiff')  % Check file extension
        InfoImage=imfinfo(FileTif);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        NumberImages=length(InfoImage);
        if NumberImages==1
            sz = [nImage mImage];
        else
            sz = [nImage mImage NumberImages];
        end
        outcome = true;
        return % Exit function
    else
        warning('function_checksize_tif: File imported is not a .tif or a .tiff file');
        return % Exit function
    end
   
else
    warning('function_checksize_tif: File does not exist');
    return % Exit function
end

end