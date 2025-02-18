function [] = Function_subdomains_display_and_save2(RVE,k_RVE,RVEparameters,number_phase,number_phase_todo,propertyname,Sub_folder_RVE,opt,infovol,dimension,p)

propertyname(1) = upper(propertyname(1)); % Uppercase for first letter

if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
    fprintf('> %s Representative Volume Element (RVE) analysis: #%i\n\n',propertyname,k_RVE);
else
    fprintf('> %s Convergence analysis: #%i\n\n',propertyname,k_RVE);
end
if RVEparameters.disp_parRVE
    disp(RVEparameters)
end

current_phase_todo = 0;
phasename_todo = cell(number_phase_todo,1);
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
        fprintf('       For the phase: %s:\n\n',char(infovol.phasename(current_phase,1)));
        disp(RVE(k_RVE).phase(current_phase_todo).statistics)
    end
end

if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
    if dimension == 3
        disp '    Cubic root length of the Representative Volume Element (RVE):';
    else
        disp '    Square root length of the Representative Section Area (RSA):';
    end
    str_name = '_RVE';
    disp(RVE(k_RVE).Equal_size1)
    if strcmp(RVEparameters.type,'C')
        if dimension == 3
            disp '    Square root length of the Representative Section Area (RSA):';
        else
            disp '    Representative length (RL):';
        end
        disp(RVE(k_RVE).Equal_size2);
    elseif strcmp(RVEparameters.type,'D')
        disp '    Representative length (RL):';
        disp(RVE(k_RVE).Equal_size2);
    end
else
    if dimension == 3
        disp '    Cubic root length to reach convergence:';
    else
        disp '    Square root length to reach convergence:';
    end
    str_name = '_CONV';
    disp(RVE(k_RVE).Equal_size1)
    if strcmp(RVEparameters.type,'G')
        disp '    Length to reach convergence:';
        disp(RVE(k_RVE).Equal_size2);
    elseif strcmp(RVEparameters.type,'H')
        if dimension == 3
            disp '    Square root length to reach convergence:';
        else
            disp '    Length to reach convergence:';
        end
        disp(RVE(k_RVE).Equal_size2);
    end    
end

propertyname = function_remove_emptyandspecialcharacter_string(propertyname);

if opt.save.xls
    filename = [propertyname str_name]; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    % Data : RVE parameters
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
        DATA_writetable.sheet(1).name='RVE parameters';
    else
        DATA_writetable.sheet(1).name='Convergence parameters';
    end
    DATA_writetable.sheet(1).table=RVE(k_RVE).parameters;
    % Data : Subdomains info
    DATA_writetable.sheet(2).name='Subdomains info';
    DATA_writetable.sheet(2).table=RVE(k_RVE).info;
    % Data : Subdomains results
    DATA_writetable.sheet(3).name='Subdomains results';
    DATA_writetable.sheet(3).table=RVE(k_RVE).results;
    % Data : Subdomains statistics
    for current_phase_todo=1:1:number_phase_todo
        sheet=3+current_phase_todo;
        DATA_writetable.sheet(sheet).name=['Statistics_' char(phasename_todo(current_phase_todo))];
        DATA_writetable.sheet(sheet).table=RVE(k_RVE).phase(current_phase_todo).statistics;
    end
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
        % Data : Size of the Representative Volume elements
        if dimension == 3
            DATA_writetable.sheet(sheet+1).name='RVE cubic root (<)';
            DATA_writetable.sheet(sheet+1).table=RVE(k_RVE).Inferior_size1;
            DATA_writetable.sheet(sheet+2).name='RVE cubic root (=)';
            DATA_writetable.sheet(sheet+2).table=RVE(k_RVE).Equal_size1;
            DATA_writetable.sheet(sheet+3).name='RVE cubic root (>)';
            DATA_writetable.sheet(sheet+3).table=RVE(k_RVE).Superior_size1;
            if strcmp(RVEparameters.type,'C')
                DATA_writetable.sheet(sheet+4).name='RSA square root (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='RSA square root (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='RSA square root (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            elseif strcmp(RVEparameters.type,'D')
                DATA_writetable.sheet(sheet+4).name='Representative length (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Representative length (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Representative length (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            end
        else
            DATA_writetable.sheet(sheet+1).name='RSA square root (<)';
            DATA_writetable.sheet(sheet+1).table=RVE(k_RVE).Inferior_size1;
            DATA_writetable.sheet(sheet+2).name='RSA square root (=)';
            DATA_writetable.sheet(sheet+2).table=RVE(k_RVE).Equal_size1;
            DATA_writetable.sheet(sheet+3).name='RSA square root (>)';
            DATA_writetable.sheet(sheet+3).table=RVE(k_RVE).Superior_size1;
            if strcmp(RVEparameters.type,'C')
                DATA_writetable.sheet(sheet+4).name='Representative length (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Representative length (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Representative length (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            elseif strcmp(RVEparameters.type,'D')
                DATA_writetable.sheet(sheet+4).name='Representative length (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Representative length (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Representative length (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            end
        end

    else
        % Data : Convergence
        if dimension == 3
            DATA_writetable.sheet(sheet+1).name='Cubic root to convergence (<)';
            DATA_writetable.sheet(sheet+1).table=RVE(k_RVE).Inferior_size1;
            DATA_writetable.sheet(sheet+2).name='Cubic root to convergence (=)';
            DATA_writetable.sheet(sheet+2).table=RVE(k_RVE).Equal_size1;
            DATA_writetable.sheet(sheet+3).name='Cubic root to convergence (>)';
            DATA_writetable.sheet(sheet+3).table=RVE(k_RVE).Superior_size1;
            if strcmp(RVEparameters.type,'G')
                DATA_writetable.sheet(sheet+4).name='Length to convergence (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Length to convergence (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Length to convergence (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            elseif strcmp(RVEparameters.type,'H')
                DATA_writetable.sheet(sheet+4).name='Square root to convergence (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Square root to convergence (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Square root to convergence (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            end
        else
            DATA_writetable.sheet(sheet+1).name='Square root to convergence (<)';
            DATA_writetable.sheet(sheet+1).table=RVE(k_RVE).Inferior_size1;
            DATA_writetable.sheet(sheet+2).name='Square root to convergence (=)';
            DATA_writetable.sheet(sheet+2).table=RVE(k_RVE).Equal_size1;
            DATA_writetable.sheet(sheet+3).name='Square root to convergence (>)';
            DATA_writetable.sheet(sheet+3).table=RVE(k_RVE).Superior_size1;
            if strcmp(RVEparameters.type,'G')
                DATA_writetable.sheet(sheet+4).name='Length to convergence (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Length to convergence (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Length to convergence (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            elseif strcmp(RVEparameters.type,'H')
                DATA_writetable.sheet(sheet+4).name='Length to convergence (<)';
                DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Inferior_size2;
                DATA_writetable.sheet(sheet+5).name='Length to convergence (=)';
                DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Equal_size2;
                DATA_writetable.sheet(sheet+6).name='Length to convergence (>)';
                DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Superior_size2;
            end
        end
    end
    % Save function
    Function_Writetable(Sub_folder_RVE,filename,DATA_writetable)
end
    
end