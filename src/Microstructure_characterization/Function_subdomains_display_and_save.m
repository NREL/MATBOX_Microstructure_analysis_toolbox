function [] = Function_subdomains_display_and_save(OPTIONS,INFO,RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)

propertyname(1) = upper(propertyname(1)); % Uppercase for first letter

if (OPTIONS.displaytext==1)
    fprintf('> %s Representative Volume Element analysis: RVE #%i\n\n',propertyname,n_RVE);
    if isfield(INFO,'showrveparameters')
        showrveparameters=INFO.showrveparameters;
    else
       showrveparameters=true; 
    end
    if showrveparameters
        disp(RVE(n_RVE).parameters)
    end

    for current_phase=1:1:number_phase
        fprintf('       For the phase: %s:\n\n',INFO.phase(current_phase).name);
        disp(RVE(n_RVE).phase(current_phase).statistics)
    end
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C')
        disp '    Criterion to reach RVE size';
        disp(RVE(n_RVE).criterion)
        disp '    Equivalent cubic length (micrometer) of the Representative Volume Element (RVE):';
        disp(RVE(n_RVE).size1)
        if strcmp(RVEparameters.type,'C')
            disp '    Equivalent section length (micrometer) of the Representative Volume Element (RVE):';
            disp(RVE(n_RVE).size2)
        end
    end
end

propertyname = function_remove_emptyandspecialcharacter_string(propertyname);

if OPTIONS.save_xls==true
    filename = [propertyname '_RVE']; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    % Data : RVE parameters
    DATA_writetable.sheet(1).name='RVE parameters';
    DATA_writetable.sheet(1).table=RVE(n_RVE).parameters;
    % Data : Subdomains info
    DATA_writetable.sheet(2).name='Subdomains info';
    DATA_writetable.sheet(2).table=RVE(n_RVE).info;
    % Data : Subdomains results
    DATA_writetable.sheet(3).name='Subdomains results';
    DATA_writetable.sheet(3).table=RVE(n_RVE).results;
    % Data : Subdomains statistics
    for current_phase=1:1:number_phase
        sheet=3+current_phase;
        DATA_writetable.sheet(sheet).name=['Statistics_' INFO.phase(current_phase).name];
        DATA_writetable.sheet(sheet).table=RVE(n_RVE).phase(current_phase).statistics;
    end
    % Data : Size of the Representative Volume elements
    DATA_writetable.sheet(sheet+1).name='Cubic length of the RVE';
    DATA_writetable.sheet(sheet+1).table=RVE(n_RVE).size1;
    if strcmp(RVEparameters.type,'C')
        DATA_writetable.sheet(sheet+2).name='Section length of the RVE';
        DATA_writetable.sheet(sheet+2).table=RVE(n_RVE).size2;
    end
    % Save function
    Function_Writetable(Sub_folder_RVE,filename,DATA_writetable)
end
    
end