function [] = RVE_displayandsave(RVE,k_RVE,kMetric,pRVE,number_domain,number_domain_todo,propertyname,Sub_folder_RVE,opt,infovol,p)

propertyname(1) = upper(propertyname(1)); % Uppercase for first letter

if strcmp(pRVE.analysis,'Independent subvolumes')
    fprintf('> %s representativity analysis: #%i\n\n',propertyname,k_RVE);
else
    fprintf('> %s convergence analysis: #%i\n\n',propertyname,k_RVE);
end
if pRVE.disp_parRVE
    disp(pRVE)
end

current_domain_todo = 0;
domainname_todo = cell(number_domain_todo,1);
for current_domain=1:1:number_domain % Loop over all phases
    if p.todo(current_domain)
        current_domain_todo=current_domain_todo+1;
        domainname_todo(current_domain_todo,1) = infovol.phasename(current_domain,1);
        fprintf('       For the domain: %s:\n\n',char(infovol.phasename(current_domain,1)));
        disp(RVE(k_RVE).Metric(kMetric).res.phase(current_domain_todo).statistics);
    end
end

if strcmp(pRVE.analysis,'Independent subvolumes')
    str_name = '_RVE';
else
    str_name = '_CONV';
end

% First size
str1 = ['    ' RVE(k_RVE).Representativity_analysis ' analysis: result is FOV ' RVE(k_RVE).FOVlength_str ]; disp(str1);
disp(RVE(k_RVE).Metric(kMetric).res.Equal_size1)

% Second size
if ischar(RVE(k_RVE).Representativity_analysis_size2)
    str2 = ['    ' RVE(k_RVE).Representativity_analysis_size2 ' analysis: result is FOV ' RVE(k_RVE).Representativity_2ndlength_str ]; disp(str2);
    disp(RVE(k_RVE).Metric(kMetric).res.Equal_size2)
end

propertyname = function_remove_emptyandspecialcharacter_string(propertyname);


if opt.save.xls
    filename = [propertyname str_name]; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    % Data : RVE parameters
    if strcmp(pRVE.analysis,'Independent subvolumes')
        DATA_writetable.sheet(1).name='Representativity parameters';
    else
        DATA_writetable.sheet(1).name='Convergence parameters';
    end
    DATA_writetable.sheet(1).table=RVE(k_RVE).Metric(kMetric).res.parameters;
    % Data : Subdomains info
    DATA_writetable.sheet(2).name='Subdomains info';
    DATA_writetable.sheet(2).table=RVE(k_RVE).info;
    % Data : Subdomains results
    DATA_writetable.sheet(3).name='Subdomains results';
    DATA_writetable.sheet(3).table=RVE(k_RVE).Metric(kMetric).res.results;
    % Data : Subdomains statistics
    for current_domain_todo=1:1:number_domain_todo
        sheet=3+current_domain_todo;
        DATA_writetable.sheet(sheet).name=['Statistics_' char(domainname_todo(current_domain_todo))];
        DATA_writetable.sheet(sheet).table=RVE(k_RVE).Metric(kMetric).res.phase(current_domain_todo).statistics;
    end

    if strcmp(pRVE.analysis,'Independent subvolumes')
        % First size
        idx = find(RVE(k_RVE).Representativity_analysis == '(');
        short = RVE(k_RVE).Representativity_analysis(idx+1:end-1);
        DATA_writetable.sheet(sheet+1).name=[short ' ' RVE(k_RVE).FOVlength_str ' (<)'];
        DATA_writetable.sheet(sheet+1).table=RVE(k_RVE).Metric(kMetric).res.Inferior_size1;
        DATA_writetable.sheet(sheet+2).name=[short ' ' RVE(k_RVE).FOVlength_str ' (=)'];
        DATA_writetable.sheet(sheet+2).table=RVE(k_RVE).Metric(kMetric).res.Equal_size1;
        DATA_writetable.sheet(sheet+3).name=[short ' ' RVE(k_RVE).FOVlength_str ' (>)'];
        DATA_writetable.sheet(sheet+3).table=RVE(k_RVE).Metric(kMetric).res.Superior_size1;
        % Second size
        if ischar(RVE(k_RVE).Representativity_analysis_size2)
            idx = find(RVE(k_RVE).Representativity_analysis_size2 == '(');
            short = RVE(k_RVE).Representativity_analysis_size2(idx+1:end-1);
            DATA_writetable.sheet(sheet+4).name=[short ' ' RVE(k_RVE).Representativity_2ndlength_str ' (<)'];
            DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Metric(kMetric).res.Inferior_size2;
            DATA_writetable.sheet(sheet+5).name=[short ' ' RVE(k_RVE).Representativity_2ndlength_str ' (=)'];
            DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Metric(kMetric).res.Equal_size2;
            DATA_writetable.sheet(sheet+6).name=[short ' ' RVE(k_RVE).Representativity_2ndlength_str ' (>)'];
            DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Metric(kMetric).res.Superior_size2;
        end
    else
        % First size
        DATA_writetable.sheet(sheet+1).name=[RVE(k_RVE).FOVlength_str 'to convergence (<)'];
        DATA_writetable.sheet(sheet+1).table=RVE(k_RVE).Metric(kMetric).res.Inferior_size1;
        DATA_writetable.sheet(sheet+2).name=[RVE(k_RVE).FOVlength_str 'to convergence (=)'];
        DATA_writetable.sheet(sheet+2).table=RVE(k_RVE).Metric(kMetric).res.Equal_size1;
        DATA_writetable.sheet(sheet+3).name=[RVE(k_RVE).FOVlength_str 'to convergence (>)'];
        DATA_writetable.sheet(sheet+3).table=RVE(k_RVE).Metric(kMetric).res.Superior_size1;

        % Second size
        if ischar(RVE(k_RVE).Representativity_analysis_size2)
            % Second size
            DATA_writetable.sheet(sheet+4).name=[RVE(k_RVE).Representativity_2ndlength_str 'to convergence (<)'];
            DATA_writetable.sheet(sheet+4).table=RVE(k_RVE).Metric(kMetric).res.Inferior_size2;
            DATA_writetable.sheet(sheet+5).name=[RVE(k_RVE).Representativity_2ndlength_str 'to convergence (=)'];
            DATA_writetable.sheet(sheet+5).table=RVE(k_RVE).Metric(kMetric).res.Equal_size2;
            DATA_writetable.sheet(sheet+6).name=[RVE(k_RVE).Representativity_2ndlength_str 'to convergence (>)'];
            DATA_writetable.sheet(sheet+6).table=RVE(k_RVE).Metric(kMetric).res.Superior_size2;
        end
    end

    % Save function
    Function_Writetable(Sub_folder_RVE,filename,DATA_writetable)
end
    
end