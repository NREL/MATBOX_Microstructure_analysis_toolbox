function [Res] = RVE_tableresult(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, Difference_convergence, Size_convergence, pRVE,FOVlength_str,Representativity_2ndlength_str,infovol,pMET)

% Table - result for each subvolume
FOVlength_str = [upper(FOVlength_str(1)) FOVlength_str(2:end)];
if ischar(Representativity_2ndlength_str)
    Representativity_2ndlength_str = [upper(Representativity_2ndlength_str(1)) Representativity_2ndlength_str(2:end)];
    Variable_name_table = {'Subdomain_Id',[FOVlength_str ' ' infovol.unit],[Representativity_2ndlength_str ' ' infovol.unit],'n/ a'};
else
    Variable_name_table = {'Subdomain_Id',[FOVlength_str ' ' infovol.unit],'n/a','n/ a'};
end
Variable_name_table_tmp = Variable_name_table;

number_domain = length(pMET.domain_name);
for current_domain=1:number_domain % Loop over all phases
    Variable_name_table(4+current_domain) = pMET.domain_name(current_domain);
end
Res.results = array2table(Property_eachsubdomain,'VariableNames',Variable_name_table);

% Table - Statistics calculated from these results (i.e. the final output of the subdomains analysis)
%Variable_name_table = [Variable_name_table 'Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum','Standard error of the mean','95% confidence interval','N required for 1% error','N required for 5% error','N required for 10% error'];
Variable_name_table = [Variable_name_table_tmp(2:end) 'Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum','Maximum-Minimum','Standard error of the mean','95% confidence interval'];
for current_domain=1:1:number_domain
    %Res.phase(current_domain).statistics = table(Property_subdomains_statistics(:,1,current_domain),Property_subdomains_statistics(:,2,current_domain),Property_subdomains_statistics(:,3,current_domain),Property_subdomains_statistics(:,4,current_domain),Property_subdomains_statistics(:,5,current_domain),Property_subdomains_statistics(:,6,current_domain),Property_subdomains_statistics(:,7,current_domain),Property_subdomains_statistics(:,8,current_domain),Property_subdomains_statistics(:,9,current_domain),Property_subdomains_statistics(:,13,current_domain),2*Property_subdomains_statistics(:,14,current_domain),Property_subdomains_statistics(:,15,current_domain),Property_subdomains_statistics(:,16,current_domain),Property_subdomains_statistics(:,17,current_domain),...
    %    'VariableNames',Variable_name_table);
    Res.phase(current_domain).statistics = table(Property_subdomains_statistics(:,1,current_domain),Property_subdomains_statistics(:,2,current_domain),Property_subdomains_statistics(:,3,current_domain),Property_subdomains_statistics(:,4,current_domain),Property_subdomains_statistics(:,5,current_domain),Property_subdomains_statistics(:,6,current_domain),Property_subdomains_statistics(:,7,current_domain),Property_subdomains_statistics(:,8,current_domain),Property_subdomains_statistics(:,9,current_domain),Property_subdomains_statistics(:,10,current_domain),Property_subdomains_statistics(:,14,current_domain),2*Property_subdomains_statistics(:,15,current_domain),...
        'VariableNames',Variable_name_table);    
end

% Table - Size of the Representative Volume element
% Size_RVE=zeros(n_threshold,number_domain,3,2);
% *,*,1 RVE size inferior to
% *,*,2 RVE size equal to
% *,*,3 RVE size is superior to
% *,*,*,1: with size1
% *,*,*,2: with size2
if strcmp(pRVE.analysis,'Independent subvolumes')
    Variable_name_table=cell(1,number_domain+2);

    if strcmp(pRVE.threshold_subs_choice,'Relative standard deviation')
        Variable_name_table(1) = {'Relative standard deviation threshold (%)'};
    elseif strcmp(pRVE.threshold_subs_choice,'Standard deviation')
        if ~isempty(pMET.unit)
            Variable_name_table(1) = {['Standard deviation threshold ' pMET.unit]};
        else
            Variable_name_table(1) = {'Standard deviation threshold'};
        end
    elseif strcmp(pRVE.threshold_subs_choice,'Maximum - minimum')
        if ~isempty(pMET.unit)
            Variable_name_table(1) = {['Maximum - minimum ' pMET.unit]};
        else
            Variable_name_table(1) = {'Maximum - minimum'};
        end        
    end

    for current_domain=1:number_domain % Loop over all phases
        Variable_name_table(current_domain+1) = pMET.domain_name(current_domain);
    end
    
    thresholds = pRVE.threshold_subs_val;
    unit_RVE = cell(length(thresholds),1);
    for k_threshold=1:1:length(thresholds)
        unit_RVE(k_threshold) = {infovol.unit};
    end

    Variable_name_table(end) = {[FOVlength_str ', unit']};

    sz = [length(thresholds) number_domain+2];
    varTypes = cell(1,number_domain+2);
    for k=1:number_domain+1
        varTypes(k)={'double'};
    end
    varTypes(end)={'string'};
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);

    T(:,1).Variables = thresholds';
    T(:,2:number_domain+1).Variables = Size_RVE(:,:,1,1);
    T(:,end) = unit_RVE;
    Res.Inferior_size1 = T;

    T(:,2:number_domain+1).Variables = Size_RVE(:,:,2,1);
    Res.Equal_size1 = T;

    T(:,2:number_domain+1).Variables = Size_RVE(:,:,3,1);
    Res.Superior_size1 = T;    

    if ischar(Representativity_2ndlength_str)
        Variable_name_table(end) = {[Representativity_2ndlength_str ', unit']};
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
        T(:,1).Variables = thresholds';
        T(:,2:number_domain+1).Variables = Size_RVE(:,:,1,2);
        T(:,end) = unit_RVE;
        Res.Inferior_size2 = T;
        T(:,2:number_domain+1).Variables = Size_RVE(:,:,2,2);
        Res.Equal_size2 = T;
        T(:,2:number_domain+1).Variables = Size_RVE(:,:,3,2);
        Res.Superior_size2 = T;
    end
    
else
    % Table - relative diffference
    Variable_name_table = Variable_name_table_tmp(2:3);
    for current_domain=1:current_domain % Loop over all phases
        Variable_name_table(2+current_domain) = pMET.domain_name(current_domain);
    end
    Res.Difference_convergence = array2table(Difference_convergence,'VariableNames',Variable_name_table);

    % Table - Size of the relative difference convergence
    % Size_convergence=zeros(n_threshold,number_domain,3,2);
    % *,*,1 convergence size inferior to
    % *,*,2 convergence size equal to
    % *,*,3 convergence size is superior to
    % *,*,*,1: with size1
    % *,*,*,2: with size2
    Variable_name_table=cell(1,number_domain+2);
    if strcmp(pRVE.threshold_onesub_choice,'Relative difference')
        Variable_name_table(1) = {'Relative difference threshold (%)'};
    elseif strcmp(pRVE.threshold_onesub_choice,'Difference')
        if ~isempty(pMET.unit)
            Variable_name_table(1) = {['Difference threshold ' pMET.unit]};
        else
            Variable_name_table(1) = {'Difference threshold'};
        end     
    end

    for current_domain=1:1:number_domain % Loop over all phases
        Variable_name_table(current_domain+1) = pMET.domain_name(current_domain);
    end
    
    thresholds = pRVE.threshold_onesub_val;
    unit_RVE = cell(length(thresholds),1);
    for k_threshold=1:1:length(thresholds)
        unit_RVE(k_threshold) = {infovol.unit};
    end

    Variable_name_table(end) = {[FOVlength_str ', unit']};

    sz = [length(thresholds) number_domain+2];
    varTypes = cell(1,number_domain+2);
    for k=1:number_domain+1
        varTypes(k)={'double'};
    end
    varTypes(end)={'string'};
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);

    T(:,1).Variables = thresholds';
    T(:,2:number_domain+1).Variables = Size_convergence(:,:,1,1);
    T(:,end) = unit_RVE;
    Res.Inferior_size1 = T;

    T(:,2:number_domain+1).Variables = Size_convergence(:,:,2,1);
    Res.Equal_size1 = T;

    T(:,2:number_domain+1).Variables = Size_convergence(:,:,3,1);
    Res.Superior_size1 = T;   

    if ischar(Representativity_2ndlength_str)
        Variable_name_table(end) = {[Representativity_2ndlength_str ', unit']};
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
        T(:,1).Variables = thresholds';
        T(:,2:number_domain+1).Variables = Size_convergence(:,:,1,2);
        T(:,end) = unit_RVE;
        Res.Inferior_size2 = T;
        T(:,2:number_domain+1).Variables = Size_convergence(:,:,2,2);
        Res.Equal_size2 = T;
        T(:,2:number_domain+1).Variables = Size_convergence(:,:,3,2);
        Res.Superior_size2 = T;
    end      
end

% Table - RVE parameters
MyFieldNames = fieldnames(pRVE);
cell_str=[];
for k=1:length(MyFieldNames)
    cell_str =[cell_str; {num2str(getfield(pRVE,MyFieldNames{k}))}];
end
Res.parameters = table(MyFieldNames,cell_str,'VariableNames',{'RVE_parameters' 'Values'});

end

