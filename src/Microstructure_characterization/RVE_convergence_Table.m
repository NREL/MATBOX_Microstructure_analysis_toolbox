function [RVE] = RVE_convergence_Table(p, RVE)

p.propertyname(1) = upper(p.propertyname(1)); % Uppercase for first letter

%% CREATE TABLES
str_table = {['FOV ' char(p.Length_name(1)) ' ' p.infovol.unit], []};
if p.nsize == 2
    str_table(2) = {['FOV ' char(p.Length_name(2)) ' ' p.infovol.unit]};
end

n_threshold = length(p.thresholds);
Variable_name_table=cell(1,n_threshold+1);
for k_threshold=1:1:n_threshold % Loop over all phases
    Variable_name_table(k_threshold+1) = {['RelStdDev ' num2str(p.thresholds(k_threshold),'%1.2f')]};
end

sz_table = size(p.Result_RVEconv(:,:,1,1,1));
varTypes = cell(1,n_threshold+1);
for k=1:1:n_threshold+1
    varTypes(k)={'double'};
end

for k=1:1:3 % <, = , >
    current_domain_todo = 0;
    for current_domain=1:1:p.number_domain % Loop over all phases
        if p.todo(current_domain)
            current_domain_todo=current_domain_todo+1;
            for ksize=1:1:p.nsize
                Variable_name_table(1) = {char(str_table(ksize))};
                T = table('Size',sz_table,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
                T.Variables = p.Result_RVEconv(:,:,current_domain_todo,k,ksize);              
                RVE.RVEconvergenceanalysis(current_domain_todo,ksize,k).table = T;
            end
        end
    end
end

%% SAVE TABLES
if p.opt.save.xls
    filename = [p.propertyname '_RVE_FOVconv']; % Filename without extension
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
        current_domain_todo = 0;
        for current_domain=1:1:p.number_domain % Loop over all phases
            if p.todo(current_domain)
                current_domain_todo=current_domain_todo+1;
                for ksize=1:1:p.nsize
                    sheet=sheet+1;
                    DATA_writetable.sheet(sheet).name=[char(p.infovol.phasename(current_domain)) ' ' char(p.Length_name(ksize)) str_];
                    DATA_writetable.sheet(sheet).table=RVE.RVEconvergenceanalysis(current_domain_todo,ksize,k).table;
                end
            end
        end
    end
    % Save function
    Function_Writetable(p.Sub_folder_RVE,filename,DATA_writetable)
end

end