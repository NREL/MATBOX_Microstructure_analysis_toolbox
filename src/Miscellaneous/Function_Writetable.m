function [] = Function_Writetable(Current_folder,filename,DATA_writetable)

Number_of_sheet = numel(DATA_writetable.sheet); % Number of sheet
if ispc % One table per sheet
    f = fullfile(Current_folder,[filename '.xlsx']); % Valid full name for all platforms
    if exist(f, 'file') % If the file exits, it is deleted
        delete(f);
    end

    for current_sheet = 1:1:Number_of_sheet % Add extra sheet, with their name
        if isempty(DATA_writetable.sheet(current_sheet).name)
            DATA_writetable.sheet(current_sheet).name=num2str(current_sheet);
        end
        name_ = DATA_writetable.sheet(current_sheet).name;

        lists = {':' '\' '/' '?' '*' '[' ']'}; % These characters can't be used in a excel sheet name
        if length(name_)>31
            new = '';
        else
            new = '-';
        end
        for k=1:1:length(lists)
            name_ = strrep(name_,char(lists(k)),new);
        end

        if length(name_)>31
            name_=name_(1:31);
        end

        if ~isempty(DATA_writetable.sheet(current_sheet).table)
            Headline = DATA_writetable.sheet(current_sheet).table.Properties.VariableNames;
            Values = table2cell(DATA_writetable.sheet(current_sheet).table);
            T=[Headline;Values];
            writecell(T,f,'Sheet',name_);
            %writecell(T,f,'Sheet',current_sheet);
        end
    end

else % One csv per sheet
    for current_sheet = 1:1:Number_of_sheet % Add extra sheet, with their name
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
        f = fullfile(Current_folder,[filename '_' name_ '.csv']); % Valid full name for all platforms
        if exist(f, 'file') % If the file exits, it is deleted
            delete(f);
        end        
        if ~isempty(DATA_writetable.sheet(current_sheet).table)
            Headline = DATA_writetable.sheet(current_sheet).table.Properties.VariableNames;
            Values = table2cell(DATA_writetable.sheet(current_sheet).table);
            T=[Headline;Values];
            writecell(T,f);
        end
    end
end

end