function [] = Function_Writetable(Current_folder,filename,DATA_writetable)

%% THIRD-PARTY FILE

% xlscol
% Downloaded from the MATLAB FILE EXCHANGE
%
% License
% Copyright (c) 2010, Kevin Crosby
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%% CLOSE EXCEL GHOST PROCESSES
% See https://www.mathworks.com/matlabcentral/answers/98261-why-do-some-excel-processes-fail-to-terminate-after-using-xlsread-xlswrite-or-xlsfinfo-in-matlab

% If you have ghost excel processes that clog the memory, you can kill all excel processes with
% system('taskkill /F /IM EXCEL.EXE')
% But it will aslo kill active excel processes, so consider it as a last resort choice.
% Current option is to use the command taskkill (see %% GOST PROCESS REMOVAL: before and after sections)

if ispc

    %% GHOST PROCESS REMOVAL: before
    [~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
    tasks=sscanf(tasks,'%s');
    exe=strfind(tasks,'.EXE');
    console=strfind(tasks,'Console');
    pid_before=zeros(1,size(exe,2));
    for i=1:1:size(exe,2)
        pid_before(1,i)=convertCharsToStrings(tasks((exe(i)+4):(console(i)-1)));
    end

    %% CREATE THE XLS FILE
    warning('off','MATLAB:xlswrite:AddSheet'); % Suppress warning when the requested worksheet does not exist and is created
    full_path_xls=[Current_folder filename '.xls']; % Path of the xlsx file
    if exist(full_path_xls, 'file') % If the Excel file exits, it is deleted
        delete(full_path_xls);
    end
    Excel = actxserver('Excel.Application'); % Open Excel Automation server
    ExcelWorkbook = Excel.workbooks.Add; % Create an xls file
    ExcelWorkbook.SaveAs(full_path_xls,1); % Save it
    ExcelWorkbook.Close(false); % Close the workbook
    
    %% PREPARE SHEETS (NUMBER AND HEADER FONT)
    Number_of_sheet = numel(DATA_writetable.sheet); % Number of sheet
    invoke(Excel.Workbooks,'Open',full_path_xls); % Open the Excel file
    Sheets = Excel.ActiveWorkBook.Sheets; % Attribute a handle to the Sheets
    if isempty(DATA_writetable.sheet(1).name)
        DATA_writetable.sheet(1).name='1';
    end
    name_ = DATA_writetable.sheet(1).name;
    if length(name_)>31
        name_=function_remove_emptyandspecialcharacter_string(name_);
    end
    if length(name_)>31
        name_=name_(1:31);
    end
    Sheets.Item(Sheets.Count).Name = name_; % Name the first sheet
    for current_sheet = 2:1:Number_of_sheet % Add extra sheet, with their name
        Sheets = Excel.ActiveWorkBook.Sheets;
        Sheets.Add([], Sheets.Item(Sheets.Count));
        if isempty(DATA_writetable.sheet(current_sheet).name)
            DATA_writetable.sheet(current_sheet).name=num2str(current_sheet);
        end
        name_ = DATA_writetable.sheet(current_sheet).name;
        if length(name_)>31
            name_=function_remove_emptyandspecialcharacter_string(name_);
        end
        if length(name_)>31
            name_=name_(1:31);
        end
        Sheets.Item(Sheets.Count).Name = name_;
    end
    Sheets = Excel.ActiveWorkBook.Sheets; % Update the Sheets handle
    for current_sheet = 1:1:Number_of_sheet
        % Header location (i.e Excel range of the header)
        dimension = size(DATA_writetable.sheet(current_sheet).table);
        if dimension(2)~=0
            letter = xlscol(dimension(2));
            corner_1='A1';
            corner_2= [letter num2str(1)];
            range = [corner_1 ':' corner_2];
            eSheet = Sheets.get('Item',current_sheet); eSheet.Activate; % Make the current sheet active.
            HeaderRange = get(Excel.Activesheet, 'Range', range);  % Select the header
            HeaderRange.Font.Bold = true; % Make the header bold
        end
    end
    
    %% SAVE AND CLOSE BOTH FILE AND EXCEL SERVER
    invoke(Excel.ActiveWorkbook,'Save'); % Save the Excel file
    Excel.Quit; Excel.delete; clear Excel; % Close the Excel Automation server. Will fail is other app is using the file
    
    %% GHOST PROCESS REMOVAL: after
    [~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
    tasks=sscanf(tasks,'%s');
    exe=strfind(tasks,'.EXE');
    console=strfind(tasks,'Console');
    pid_after=zeros(1,size(exe,2));
    for i=1:1:size(exe,2)
        pid_after(1,i)=convertCharsToStrings(tasks((exe(i)+4):(console(i)-1)));
    end
    if size(pid_after,2)>size(pid_before,2)
        command=['taskkill /f /PID ',mat2str(pid_after(1,end))];
        [status,cmdout] = system(command); % no print
    end

    %% GHOST PROCESS REMOVAL: before
    [~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
    tasks=sscanf(tasks,'%s');
    exe=strfind(tasks,'.EXE');
    console=strfind(tasks,'Console');
    pid_before=zeros(1,size(exe,2));
    for i=1:1:size(exe,2)
        pid_before(1,i)=convertCharsToStrings(tasks((exe(i)+4):(console(i)-1)));
    end

    %% WRITE DATA FOR EACH SHEET
    for current_sheet = 1:1:Number_of_sheet % Loop over all sheet
        dimension = size(DATA_writetable.sheet(current_sheet).table); % Get the data range for excel (Do not specify the range may bug the function, especially if the number of row is greater than 1000)
        if dimension(2)~=0
            dimension(1)=dimension(1)+1; % % Number of row. The 1st row which contains the names is not counted with the size function
            letter = xlscol(dimension(2)); % Excel letter which corresponds to the number of column (xlscol)
            % Data range for the headline
            corner_1='A1';
            corner_2= [letter num2str(1)];
            range_headline = [corner_1 ':' corner_2];
            % Data range for the values
            corner_1='A2';
            corner_2= [letter num2str(dimension(1))];
            range_values = [corner_1 ':' corner_2];
            % Extract the headline from the table
            Headline = DATA_writetable.sheet(current_sheet).table.Properties.VariableNames;
            % Extract the numeric values from the table
            Values = table2cell(DATA_writetable.sheet(current_sheet).table);
            
            %Values = DATA_writetable.sheet(current_sheet).table{:,Headline};
            % Write headline to the xls file
            [status1,message1]   = xlswrite(full_path_xls,Headline,current_sheet,range_headline);
            % Write Numeric values to the xls file
            [status2,message2]   = xlswrite(full_path_xls,Values,current_sheet,range_values);
            if status1~=1 || status2~=1
                disp 'WARNING: the  xlswrite function did not work. Results are not saved in a xls file';
                disp(message1);
                disp(message2);
            end
        end
    end

    %% GHOST PROCESS REMOVAL: after
    [~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
    tasks=sscanf(tasks,'%s');
    exe=strfind(tasks,'.EXE');
    console=strfind(tasks,'Console');
    pid_after=zeros(1,size(exe,2));
    for i=1:1:size(exe,2)
        pid_after(1,i)=convertCharsToStrings(tasks((exe(i)+4):(console(i)-1)));
    end
    if size(pid_after,2)>size(pid_before,2)
        command=['taskkill /f /PID ',mat2str(pid_after(1,end))];
        [status,cmdout] = system(command); % no print
    end

end

end