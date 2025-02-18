function [] = function_savefig(Fig, Current_folder, filename, opt)
%function_savefig is a custom file to save figure

nargin; % Number of input variable when the function is call
if nargin == 3 
    opt.fig_infig = true;
    opt.fig_format = {'png'};
end

[ filename ] = function_remove_emptyandspecialcharacter_string(filename);

fullpath=fullfile(Current_folder, filename); % Path without extension

if ~isfield(opt,'overwritte')
    opt.overwritte=true;
end

% To be removed
% opt.fig_infig = true;
% opt.fig_format = {'png'};

if opt.fig_infig
    if ~isfile([fullpath '.fig']) || opt.overwritte
        savefig(Fig,fullpath) % Save in fig format
    else
        tmp=fullpath;
        while isfile([tmp '.fig'])
            tmp=[tmp '_bis'];
        end
        savefig(Fig,tmp)
    end
end


number_format = length(opt.fig_format); % Save the figure in different format
if number_format>0
    for k=1:1:number_format % Loop over the different format
        current_format = char(opt.fig_format(k)); % Select format
        switch current_format
            case 'png' % Standard png save
                if ~isfile([fullpath '.png']) || opt.overwritte
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

