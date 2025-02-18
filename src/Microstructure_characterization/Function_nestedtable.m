function [RVE] = Function_nestedtable(RVE,k_RVE,RVEparameters,number_phase,propertyname,Result_nestedRVE,Sub_folder_RVE,opt,infovol,p)

propertyname(1) = upper(propertyname(1)); % Uppercase for first letter

%% CREATE TABLES

str_table={['Equivalent length: FOV cubic root ' infovol.unit], []};
length_table = {'cubic root',[]};
if strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'H')
    str_table(2)={['Equivalent length: FOV square root ' infovol.unit]};
    length_table(2) = {'square root'};
    ky = [1 2]; nsize=2;
elseif strcmp(RVEparameters.type,'D')
    str_table(2)={['FOV length ' infovol.unit]};
    length_table(2) = {'length'};
    ky = [1 3]; nsize=2;    
elseif strcmp(RVEparameters.type,'G')
    str_table(2)={['FOV length ' infovol.unit]};
    length_table(2) = {'length'};
    ky = [1 3]; nsize=2;
else
    ky =1; nsize=1;
end

if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
    thresholds = RVEparameters.threshold_std;
else
    thresholds = RVEparameters.threshold_reldiff;
end
n_threshold = length(thresholds);

Variable_name_table=cell(1,n_threshold+1);
for k_threshold=1:1:n_threshold % Loop over all phases
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
        Variable_name_table(k_threshold+1) = {['RelStdDev ' num2str(thresholds(k_threshold),'%1.2f')]};
    else
        Variable_name_table(k_threshold+1) = {['RelDiff ' num2str(thresholds(k_threshold),'%1.2f')]};
    end
end

sz_table = size(Result_nestedRVE(:,:,1,1,1));
varTypes = cell(1,n_threshold+1);
for k=1:1:n_threshold+1
    varTypes(k)={'double'};
end

for k=1:1:3 % <, = , >
    current_phase_todo = 0;
    for current_phase=1:1:number_phase % Loop over all phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            for ksize=1:1:nsize
                Variable_name_table(1) = {char(str_table(ksize))};
                T = table('Size',sz_table,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
                T.Variables = Result_nestedRVE(:,:,current_phase_todo,k,ky(ksize));
                RVE(k_RVE).nestedanalysis(current_phase_todo,ksize,k).table = T;
            end
        end
    end
end

%% SAVE TABLES
if opt.save.xls
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D') 
        filename = [propertyname '_RVEnested']; % Filename without extension
    else
        filename = [propertyname '_CONVnested']; % Filename without extension
    end        
    % Prepare the data
    clear DATA_writetable
    sheet = 0;
    for kk=1:1:3
        if kk==1
            k=2; str_ = ' =';
        elseif kk==2
            k=1; str_ = ' <';
        else
            k=3; str_ = ' >';
        end
        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over all phases
            if p.todo(current_phase)
                current_phase_todo=current_phase_todo+1;
                for ksize=1:1:nsize
                    sheet=sheet+1;
                    DATA_writetable.sheet(sheet).name=[char(infovol.phasename(current_phase)) ' ' char(length_table(ksize)) str_];
                    DATA_writetable.sheet(sheet).table=RVE(k_RVE).nestedanalysis(current_phase_todo,ksize,k).table;
                end
            end
        end
    end
    % Save function
    Function_Writetable(Sub_folder_RVE,filename,DATA_writetable)
end

end