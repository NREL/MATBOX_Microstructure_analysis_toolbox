function [ str_new ] = function_remove_emptyandspecialcharacter_string(str_old)
%function_remove_emptyandspecialcharacter_string is used to change
%   string 'str1 str2-str3' to 'str1_str2str3'

str_new=str_old; % Initialization
number_of_suppression = 0; % Initialise suppression counter
for current_iteration = 1:1:length(str_old) % Loop over all letter
    current_position_new_string = current_iteration - number_of_suppression;
    current_letter = str_old(current_iteration); % Select letter of old string
    if strcmp(current_letter,' ') % Check
        str_new(current_position_new_string) = '_';
    end
    if strcmp(current_letter,',') % Check
        str_new(current_position_new_string) = '_';
    end    
    if strcmp(current_letter,':') % Check
        str_new(current_position_new_string) = '_';
    end        
    if strcmp(current_letter,'-')
        str_new(current_position_new_string) = '';
        number_of_suppression = number_of_suppression+1; % Increment suppression counter
    end
    if strcmp(current_letter,'.')
        str_new(current_position_new_string) = '';
        number_of_suppression = number_of_suppression+1; % Increment suppression counter
    end    
    
    if strcmp(current_letter,'/') || strcmp(current_letter,'(') || strcmp(current_letter,')')% Check
        str_new(current_position_new_string) = '_';
    end
    

    
end

str_old=str_new;
clear str_new
previous_letter='none';
k=0;
for current_iteration = 1:1:length(str_old) % Loop over all letter
    current_letter = str_old(current_iteration); % Select letter of old string
    if strcmp(current_letter,'_') &&  strcmp(previous_letter,'_')% Check
        foo=1;
    else
        k=k+1;
        str_new(k)=str_old(current_iteration);
    end
    previous_letter = current_letter;
end

str_old=str_new;
clear str_new
if strcmp(str_old(end) ,'_')
    str_new=str_old(1:end-1);
else
    str_new=str_old;
end

end

