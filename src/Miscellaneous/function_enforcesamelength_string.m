function [ string_withspecificlength ] = function_enforcesamelength_string( string_withanylength, specified_length)
%function_enforcesamelength_strings add spaces to a string to match the specified lenght

string_length = length(string_withanylength);
if string_length<specified_length
    string_withspecificlength = string_withanylength; % Initialisation
    for p=1:1:specified_length-string_length
        string_withspecificlength = [string_withspecificlength ' '];
    end

elseif string_length==specified_length
      string_withspecificlength = string_withanylength; % Keep as it is
else
    foo=1;
end


end

