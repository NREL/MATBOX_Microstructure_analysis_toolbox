function [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, RVEparameters,RVE,Size_RVE,n_RVE,number_phase,INFO)

%% MANAGING RESULTS

% Table - result for each subvolume
if strcmp(RVEparameters.type,'E')
    Variable_name_table={'Subdomain_Id','Equivalent_cubic_length_um','Length_um'};
else
    Variable_name_table={'Subdomain_Id','Equivalent_cubic_length_um','Section_length_um'};
end
for current_phase=1:1:number_phase
    Variable_name_table(3+current_phase)={INFO.phase(current_phase).filename};
end
RVE(n_RVE).results = array2table(Property_eachsubdomain,'VariableNames',Variable_name_table);

% Table - Statistics calculated from these results (i.e. the final output of the subdomains analysis)
clear Variable_name_table;
if strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'E')
    if strcmp(RVEparameters.type,'C')
        Variable_name_table={'Equivalent_cubic_length_um','Section_length_um','Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum'};
    elseif strcmp(RVEparameters.type,'E')
        Variable_name_table={'Equivalent_cubic_length_um','Length_um','Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum'};
    end
    for current_phase=1:1:number_phase
        RVE(n_RVE).phase(current_phase).statistics = table(Property_subdomains_statistics(:,1,current_phase),Property_subdomains_statistics(:,2,current_phase),Property_subdomains_statistics(:,3,current_phase),Property_subdomains_statistics(:,4,current_phase),Property_subdomains_statistics(:,5,current_phase),Property_subdomains_statistics(:,6,current_phase),Property_subdomains_statistics(:,7,current_phase),Property_subdomains_statistics(:,8,current_phase),...
            'VariableNames',Variable_name_table);
    end
else
    Variable_name_table={'Equivalent_cubic_length_um','Number_of_subdomain','Mean','Standard_deviation','Relative_standard_deviation_percents','Minimum','Maximum'};
    for current_phase=1:1:number_phase
        RVE(n_RVE).phase(current_phase).statistics = table(Property_subdomains_statistics(:,1,current_phase),Property_subdomains_statistics(:,3,current_phase),Property_subdomains_statistics(:,4,current_phase),Property_subdomains_statistics(:,5,current_phase),Property_subdomains_statistics(:,6,current_phase),Property_subdomains_statistics(:,7,current_phase),Property_subdomains_statistics(:,8,current_phase),...
            'VariableNames',Variable_name_table);
    end
end

% Table - Size of the Representative Volume element
RVE(n_RVE).criterion =  table({RVEparameters.threshold_std},{RVEparameters.threshold_numbersubvolumes},'VariableNames',{'Standard_deviation_criterion_percents' 'Minimal_number_subdomains' });
if strcmp(RVEparameters.type,'C')
    RVE(n_RVE).size2 = table({INFO.phase.name}',Size_RVE(1,:,2)',Size_RVE(2,:,2)',Size_RVE(3,:,2)',...
        'VariableNames',{'Phase' 'RVE_size_is_inferior_to' 'RVE_size_is_equal_to' 'RVE_size_is_superior_to'});
end
RVE(n_RVE).size1 = table({INFO.phase.name}',Size_RVE(1,:,1)',Size_RVE(2,:,1)',Size_RVE(3,:,1)',...
    'VariableNames',{'Phase' 'RVE_size_is_inferior_to' 'RVE_size_is_equal_to' 'RVE_size_is_superior_to'});

% Table - RVE parameters
MyFieldNames = fieldnames(RVEparameters);
cell_str=[];
for k=1:1:length(MyFieldNames)
    cell_str =[cell_str; {num2str(getfield(RVEparameters,MyFieldNames{k}))}];
end
RVE(n_RVE).parameters = table(fieldnames(RVEparameters),cell_str,'VariableNames',{'RVE_parameters' 'Values'});





end

