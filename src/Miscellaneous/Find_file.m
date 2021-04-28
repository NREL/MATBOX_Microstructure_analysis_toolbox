function [] = Find_file(filename,str_end,errormessage)
if ispc
    separation_folder = '\';
else
    separation_folder = '/';
end
path_app = mfilename('fullpath');
higherlevelfolder = extractBetween(path_app,path_app(1:5),str_end,'Boundaries','inclusive');
filepath_path = [char(higherlevelfolder) separation_folder 'Documentation' separation_folder filename];
    if exist(filepath_path,'file')
    open(filepath_path);
    else
        warning(['MATLAB did not find the file ' filename]);
        if ~isempty(errormessage)
            warning(errormessage);
        end
    end
end

