function [] = function_savefig(Fig, Current_folder, filename, OPTIONS)
%function_savefig is a custom file to save figure

nargin; % Number of input variable when the function is call
if nargin == 3 
    OPTIONS.savefig_infig = true;
    OPTIONS.savefig_informat = {'png'};
end

[ filename ] = function_remove_emptyandspecialcharacter_string(filename);

fullpath=[Current_folder filename]; % Path without extension

if ~isfield(OPTIONS,'overwritte')
    OPTIONS.overwritte=true;
end

if OPTIONS.savefig_infig == true
    if ~isfile([fullpath '.fig']) || OPTIONS.overwritte
        savefig(Fig,fullpath) % Save in fig format
    else
        tmp=fullpath;
        while isfile([tmp '.fig'])
            tmp=[tmp '_bis'];
        end
        savefig(Fig,tmp)
    end
end

number_format = length(OPTIONS.savefig_informat); % Save the figure in different format
if number_format>0
    for k=1:1:number_format % Loop over the different format
        current_format = char(OPTIONS.savefig_informat(k)); % Select format
        switch current_format
            case 'png' % Standard png save
                if ~isfile([fullpath '.png']) || OPTIONS.overwritte
                    saveas(Fig,fullpath,'png');
                else
                    tmp=fullpath;
                    while isfile([tmp '.png'])
                        tmp=[tmp '_bis'];
                    end
                    saveas(Fig,tmp,'png');
                end
            otherwise
                warning('Unexpected format type. No save created.')
        end
    end
end

end

