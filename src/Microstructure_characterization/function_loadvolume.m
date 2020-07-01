function [ Phase_microstructure, outcome] = function_loadvolume( filefullpath, format, structurefield )
%function_loadvolume load a volume file
%   Phase_microstructure, outcome = function_loadvolume(filefullpath, format)
%   filefullpath must be the absolute path of the file to import. Accepted extention are .mat and .tif
%   format can be 'double', 'uint16' or 'uint8'
%   structurefield field name of the structure, in case a .mat structure is loaded
%   outcome informs about the success or failure of the loading.    

% Get file path, file name, and file extentsion
[~,~,ext] = fileparts(filefullpath);

Phase_microstructure = []; % Initialize
outcome.success = true; % Initialize

% Check file exist
if exist(filefullpath, 'file') == 2
    % Check extention
    if strcmp(ext, '.mat')
        dataloaded = load (filefullpath); % Load file
        if isnumeric(dataloaded)
            Phase_microstructure = dataloaded; % Already an array
        elseif isstruct(dataloaded) % A structure
            names = fieldnames(dataloaded); % Get all fields
            found_array = false; % Initialise
            if strcmp(structurefield, 'none') % No field given
                for k=1:1:length(names) % Loop over all fields
                    val = getfield(dataloaded, char(names(k)));
                    if isnumeric(val)
                        Phase_microstructure = val;
                        found_array=true;
                        break % Exit for loop
                    end
                end
                if found_array==false
                    outcome.success = false;
                    outcome.reason = 'Failure: mat file exists, but is a structure with no array inside';
                    return
                end
            else
                if isfield(dataloaded,structurefield)
                    val = dataloaded.structurefield;
                    if isnumeric(val)
                        Phase_microstructure = val;
                    else
                        outcome.success = false;
                        outcome.reason = 'Failure: mat file exists, with the expected field structure, but the field is not an array';
                        return
                    end
                else
                    outcome.success = false;
                    outcome.reason = 'Failure: mat file exists, but is a structure with not the expected field';
                    return                    
                end
            end
        end
        % Convert to desire format
        if strcmp(format, 'uint8')
            Phase_microstructure = uint8(Phase_microstructure);
        elseif strcmp(format, 'uint16')
            Phase_microstructure = uint16(Phase_microstructure);
        elseif strcmp(format, 'double')
            Phase_microstructure = double(Phase_microstructure);
        else
            outcome.success = false;
            outcome.reason = 'Failure: mat file exists, but format required is unexpected';
            return
        end

    elseif strcmp(ext, '.tif')
        InfoImage=imfinfo(filefullpath);
        mImage=InfoImage(1).Width; % In-plane dimension 1
        nImage=InfoImage(1).Height; % In-plane dimension 2
        NumberImages=length(InfoImage); % Third dimension 
        Phase_microstructure=zeros(nImage,mImage,NumberImages, format); % Initialize
        TifLink = Tiff(filefullpath, 'r'); % Read file
        for i=1:NumberImages % Loop over slice images
            TifLink.setDirectory(i);
            Phase_microstructure(:,:,i)=TifLink.read();
        end
        TifLink.close();
    else
        outcome.success = false;
        outcome.reason = 'Failure: file exist but extension is neither .mat, neither .tif';   
        return
    end
    
else
    outcome.success = false;
    outcome.reason = 'Failure: file does not exist';
    return
end

end

