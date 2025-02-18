function [array,outcome] = function_load_tif(FileTif,datatype_array)
%function_load_tif loads a tif file in an array
% [array,outcome] = function_load_tif(FileTif,datatype_array)
% Inputs: 
% - FileTif: full path of the tif file
% - (optional) datatype_array: numerical format of array: 'double', 'single', 'logical','int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'  
%   By default: 'uint16'
% Outputs:
% - array: MATLAB array with the user-defined datatype
% - outcome: true if algorithm has been done successfully, false otherwise.

outcome = false; % Initialize

% Set default argument
switch nargin
    case 1
        datatype_array = 'uint16';
    case 2
        foo=1;
    otherwise
        warning('function_load_tif: wrong number of arguments (1 or 2).')
        help function_load_tif
        return % Exit function
end

if isfile(FileTif) % Check existence
    [~,~,ext] = fileparts(FileTif);
    
    if strcmp(ext,'.tif') || strcmp(ext,'.tiff')  % Check file extension
        InfoImage=imfinfo(FileTif);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        NumberImages=length(InfoImage);
        array=zeros(nImage,mImage,NumberImages,datatype_array);
        TifLink = Tiff(FileTif, 'r');
        for i=1:NumberImages
            TifLink.setDirectory(i);
            array(:,:,i)=TifLink.read();
        end
        TifLink.close();
        outcome = true; % Initialize
        return % Exit function
    else
        warning('function_load_tif: File imported is not a .tif or a .tiff file');
        help function_load_tif
        return % Exit function
    end
   
else
    warning('function_load_tif: File to load does not exist');
    help function_load_tif
    return % Exit function
end

end

