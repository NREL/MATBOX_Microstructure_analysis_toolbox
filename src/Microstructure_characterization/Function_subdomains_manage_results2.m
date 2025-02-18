function [RVE] = Function_subdomains_manage_results2(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,RVE,k_RVE,number_phase,number_phase_todo,infovol,dimension,p)

% Table - result for each subvolume
if dimension == 3
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'E') || strcmp(RVEparameters.type,'F')
        Variable_name_table={'Subdomain_Id',['Equivalent_cubicroot_length_' infovol.unit],'n/a','n/ a'};
    elseif strcmp(RVEparameters.type,'G')
        Variable_name_table={'Subdomain_Id',['Equivalent_cubicroot_length_' infovol.unit],'n/a',['Length_' infovol.unit]};
    elseif strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'H')
        Variable_name_table={'Subdomain_Id',['Equivalent_cubicroot_length_' infovol.unit],['Equivalent_squareroot_length_' infovol.unit],'n/ a'};
    elseif strcmp(RVEparameters.type,'D')
        Variable_name_table={'Subdomain_Id',['Equivalent_cubicroot_length_' infovol.unit],'n/a',['Length_' infovol.unit]};
    end
else
    if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'E') || strcmp(RVEparameters.type,'F')
        Variable_name_table={'Subdomain_Id',['Equivalent_squareroot_length_' infovol.unit],'n/a','n/ a'};    
    elseif strcmp(RVEparameters.type,'G')
        Variable_name_table={'Subdomain_Id',['Equivalent_squareroot_length_' infovol.unit],'n/a',['Length_' infovol.unit]};
    elseif strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'H')
        Variable_name_table={'Subdomain_Id',['Equivalent_squareroot_length_' infovol.unit],['Length_' infovol.unit],'n/ a'};
    elseif strcmp(RVEparameters.type,'D')
        Variable_name_table={'Subdomain_Id',['Equivalent_squareroot_length_' infovol.unit],'n/a',['Length_' infovol.unit]};
    end
end
Variable_name_table_tmp = Variable_name_table;

current_phase_todo = 0;
phasename_todo = cell(number_phase_todo,1);
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        Variable_name_table(4+current_phase_todo) = infovol.phasename(current_phase,1);
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
    end
end
RVE(k_RVE).results = array2table(Property_eachsubdomain,'VariableNames',Variable_name_table);

% Table - Statistics calculated from these results (i.e. the final output of the subdomains analysis)
%Variable_name_table = [Variable_name_table 'Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum','Standard error of the mean','95% confidence interval','N required for 1% error','N required for 5% error','N required for 10% error'];
Variable_name_table = [Variable_name_table_tmp(2:end) 'Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum','Standard error of the mean','95% confidence interval'];

for current_phase_todo=1:1:number_phase_todo
    %RVE(k_RVE).phase(current_phase_todo).statistics = table(Property_subdomains_statistics(:,1,current_phase_todo),Property_subdomains_statistics(:,2,current_phase_todo),Property_subdomains_statistics(:,3,current_phase_todo),Property_subdomains_statistics(:,4,current_phase_todo),Property_subdomains_statistics(:,5,current_phase_todo),Property_subdomains_statistics(:,6,current_phase_todo),Property_subdomains_statistics(:,7,current_phase_todo),Property_subdomains_statistics(:,8,current_phase_todo),Property_subdomains_statistics(:,9,current_phase_todo),Property_subdomains_statistics(:,13,current_phase_todo),2*Property_subdomains_statistics(:,14,current_phase_todo),Property_subdomains_statistics(:,15,current_phase_todo),Property_subdomains_statistics(:,16,current_phase_todo),Property_subdomains_statistics(:,17,current_phase_todo),...
    %    'VariableNames',Variable_name_table);
    RVE(k_RVE).phase(current_phase_todo).statistics = table(Property_subdomains_statistics(:,1,current_phase_todo),Property_subdomains_statistics(:,2,current_phase_todo),Property_subdomains_statistics(:,3,current_phase_todo),Property_subdomains_statistics(:,4,current_phase_todo),Property_subdomains_statistics(:,5,current_phase_todo),Property_subdomains_statistics(:,6,current_phase_todo),Property_subdomains_statistics(:,7,current_phase_todo),Property_subdomains_statistics(:,8,current_phase_todo),Property_subdomains_statistics(:,9,current_phase_todo),Property_subdomains_statistics(:,13,current_phase_todo),2*Property_subdomains_statistics(:,14,current_phase_todo),...
        'VariableNames',Variable_name_table);    
end

% Table - Size of the Representative Volume element
% Size_RVE=zeros(n_threshold,number_phase_todo,3,2);
% *,*,1 RVE size inferior to
% *,*,2 RVE size equal to
% *,*,3 RVE size is superior to
% *,*,*,1: with size1
% *,*,*,2: with size2
if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
    Variable_name_table=cell(1,number_phase_todo+2);
    Variable_name_table(1) = {'Relative standard deviation threshold (%)'};
    for current_phase_todo=1:1:number_phase_todo % Loop over all phases
        Variable_name_table(current_phase_todo+1) = phasename_todo(current_phase_todo);
    end
    
    std_thresholds = RVEparameters.threshold_std;
    unit_RVE = cell(length(std_thresholds),1);
    for k_threshold=1:1:length(std_thresholds)
        unit_RVE(k_threshold) = {infovol.unit};
    end

    if dimension == 3
        Variable_name_table(end) = {'Cubic root, unit'};
    else
        Variable_name_table(end) = {'Square root, unit'};
    end

    sz = [length(std_thresholds) number_phase_todo+2];
    varTypes = cell(1,number_phase_todo+2);
    for k=1:1:number_phase_todo+1
        varTypes(k)={'double'};
    end
    varTypes(end)={'string'};
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);

    T(:,1).Variables = std_thresholds';
    T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,1,1);
    T(:,end) = unit_RVE;

    RVE(k_RVE).Inferior_size1 = T;
    T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,2,1);
    RVE(k_RVE).Equal_size1 = T;
    T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,3,1);
    RVE(k_RVE).Superior_size1 = T;    
    if strcmp(RVEparameters.type,'C')
        if dimension == 3
            Variable_name_table(end) = {'Square root, unit'};
        else
            Variable_name_table(end) = {'length, unit'};
        end
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
        T(:,1).Variables = std_thresholds';
        T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,1,2);
        T(:,end) = unit_RVE;
        RVE(k_RVE).Inferior_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,2,2);
        RVE(k_RVE).Equal_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,3,2);
        RVE(k_RVE).Superior_size2 = T;
    elseif strcmp(RVEparameters.type,'D')
        Variable_name_table(end) = {'Length, unit'};
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
        T(:,1).Variables = std_thresholds';
        T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,1,2);
        T(:,end) = unit_RVE;
        RVE(k_RVE).Inferior_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,2,2);
        RVE(k_RVE).Equal_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_RVE(:,:,3,2);
        RVE(k_RVE).Superior_size2 = T;
    end
end

% Table - Derivative
% if strcmp(RVEparameters.type,'E') || strcmp(RVEparameters.type,'F') || strcmp(RVEparameters.type,'G') || strcmp(RVEparameters.type,'H')
%     Variable_name_table={['Equivalent cubicroot length ' infovol.unit]};
%     for current_phase_todo=1:1:current_phase_todo % Loop over all phases
%         Variable_name_table(1+current_phase_todo) = phasename_todo(current_phase_todo);
%     end
%     RVE(k_RVE).derivative_convergence_cubicroot = array2table(derivative_convergence(:,:,1),'VariableNames',Variable_name_table);
%     if strcmp(RVEparameters.type,'G')
%         Variable_name_table={['Length ' infovol.unit]};
%         for current_phase_todo=1:1:current_phase_todo % Loop over all phases
%             Variable_name_table(1+current_phase_todo) = phasename_todo(current_phase_todo);
%         end
%         RVE(k_RVE).derivative_convergence_length = array2table(derivative_convergence(:,:,2),'VariableNames',Variable_name_table);
%     elseif strcmp(RVEparameters.type,'H')
%         Variable_name_table={['Square root ' infovol.unit]};
%         for current_phase_todo=1:1:current_phase_todo % Loop over all phases
%             Variable_name_table(1+current_phase_todo) = phasename_todo(current_phase_todo);
%         end
%         RVE(k_RVE).derivative_convergence_length = array2table(derivative_convergence(:,:,2),'VariableNames',Variable_name_table);
%     end    
% end

% Table - relative diffference
if strcmp(RVEparameters.type,'E') || strcmp(RVEparameters.type,'F') || strcmp(RVEparameters.type,'G') || strcmp(RVEparameters.type,'H')
    if dimension == 3
        if strcmp(RVEparameters.type,'G')
            Variable_name_table={['Equivalent cubicroot length ' infovol.unit], ['Length_' infovol.unit]};
        elseif strcmp(RVEparameters.type,'H')
            Variable_name_table={['Equivalent cubicroot length ' infovol.unit], ['Equivalent squareroot length ' infovol.unit]};
        else
            Variable_name_table={['Equivalent cubicroot length ' infovol.unit], 'n/a'};
        end
    else
        if strcmp(RVEparameters.type,'G')
            Variable_name_table={['Equivalent squareroot length ' infovol.unit], ['Length_' infovol.unit]};
        elseif strcmp(RVEparameters.type,'H')
            Variable_name_table={['Equivalent squareroot length ' infovol.unit], ['Length ' infovol.unit]};
        else
            Variable_name_table={['Equivalent squareroot length ' infovol.unit], 'n/a'};
        end
    end
    for current_phase_todo=1:1:current_phase_todo % Loop over all phases
        Variable_name_table(2+current_phase_todo) = phasename_todo(current_phase_todo);
    end
    RVE(k_RVE).relativedifference_convergence = array2table(relativedifference_convergence,'VariableNames',Variable_name_table);
end

% Table - Size of the relative difference convergence
% Size_convergence=zeros(n_threshold,number_phase_todo,3,2);
% *,*,1 convergence size inferior to
% *,*,2 convergence size equal to
% *,*,3 convergence size is superior to
% *,*,*,1: with size1
% *,*,*,2: with size2
if strcmp(RVEparameters.type,'E') || strcmp(RVEparameters.type,'F') || strcmp(RVEparameters.type,'G') || strcmp(RVEparameters.type,'H')
    Variable_name_table=cell(1,number_phase_todo+2);
    Variable_name_table(1) = {'Relative difference threshold (%)'};
    for current_phase_todo=1:1:number_phase_todo % Loop over all phases
        Variable_name_table(current_phase_todo+1) = phasename_todo(current_phase_todo);
    end
    
    reldiff_thresholds = RVEparameters.threshold_reldiff;
    unit_RVE = cell(length(reldiff_thresholds),1);
    for k_threshold=1:1:length(reldiff_thresholds)
        unit_RVE(k_threshold) = {infovol.unit};
    end

    if dimension == 3
        Variable_name_table(end) = {'Cubic root, unit'};
    else
        Variable_name_table(end) = {'Square root, unit'};
    end  

    sz = [length(reldiff_thresholds) number_phase_todo+2];
    varTypes = cell(1,number_phase_todo+2);
    for k=1:1:number_phase_todo+1
        varTypes(k)={'double'};
    end
    varTypes(end)={'string'};
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);

    T(:,1).Variables = reldiff_thresholds';
    T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,1,1);
    T(:,end) = unit_RVE;

    RVE(k_RVE).Inferior_size1 = T;
    T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,2,1);
    RVE(k_RVE).Equal_size1 = T;
    T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,3,1);
    RVE(k_RVE).Superior_size1 = T;    
    if strcmp(RVEparameters.type,'G')
        Variable_name_table(end) = {'Length, unit'};
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
        T(:,1).Variables = reldiff_thresholds';
        T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,1,2);
        T(:,end) = unit_RVE;
        RVE(k_RVE).Inferior_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,2,2);
        RVE(k_RVE).Equal_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,3,2);
        RVE(k_RVE).Superior_size2 = T;
    elseif strcmp(RVEparameters.type,'H')
        if dimension == 3
            Variable_name_table(end) = {'Equivalent squareroot length, unit'};
        else
            Variable_name_table(end) = {'Length, unit'};
        end
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',Variable_name_table);
        T(:,1).Variables = reldiff_thresholds';
        T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,1,2);
        T(:,end) = unit_RVE;
        RVE(k_RVE).Inferior_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,2,2);
        RVE(k_RVE).Equal_size2 = T;
        T(:,2:number_phase_todo+1).Variables = Size_convergence(:,:,3,2);
        RVE(k_RVE).Superior_size2 = T;
    end    
end

% Table - RVE parameters
MyFieldNames = fieldnames(RVEparameters);
cell_str=[];
for k=1:1:length(MyFieldNames)
    cell_str =[cell_str; {num2str(getfield(RVEparameters,MyFieldNames{k}))}];
end
RVE(k_RVE).parameters = table(MyFieldNames,cell_str,'VariableNames',{'RVE_parameters' 'Values'});

end

