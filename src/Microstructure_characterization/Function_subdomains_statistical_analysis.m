function [Property_subdomains_statistics, Size_RVE] = Function_subdomains_statistical_analysis(number_group_size,number_phase,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters)

Property_subdomains_statistics = zeros(number_group_size,8,number_phase);
% Each row corresponds to a different ideal size
% column 1 : size 1 (always the equivalent cubic length)
% column 2 : size 2
% column 3 : number of subdomains with the same size
% column 4 : mean value
% column 5 : standard deviation, only calculated if number of subdomains is >2
% column 6 : standard deviation expressed in % of the mean value
% column 7 : minimun
% column 8 : maximum
% column 9-11: subdomain aspect ratio (AR)
for current_phase=1:1:number_phase
    for group_subdomain = 1:1:number_group_size
        size1 = GROUP_SUBDOMAIN.id(group_subdomain).equivalent_cubic_length;
        size2 = GROUP_SUBDOMAIN.id(group_subdomain).equivalent_square_length;
        AR = GROUP_SUBDOMAIN.id(group_subdomain).Aspectratio;
        % Number of subdomain sharing this size
        number_subdomain_withingroup = length(GROUP_SUBDOMAIN.id(group_subdomain).domain_id);
        value_tmp =[];
        for current_subdomain=1:1:number_subdomain_withingroup
            subdomain_id=GROUP_SUBDOMAIN.id(group_subdomain).domain_id(current_subdomain);
            value_tmp=[value_tmp Property_eachsubdomain(subdomain_id,current_phase+3)];
        end
        Property_subdomains_statistics(group_subdomain,1,current_phase)=size1;
        Property_subdomains_statistics(group_subdomain,2,current_phase)=size2;
        Property_subdomains_statistics(group_subdomain,3,current_phase)=number_subdomain_withingroup;
        Property_subdomains_statistics(group_subdomain,4,current_phase)=nanmean(value_tmp);
        Property_subdomains_statistics(group_subdomain,5,current_phase)=nanstd(value_tmp);
        Property_subdomains_statistics(group_subdomain,6,current_phase)=Property_subdomains_statistics(group_subdomain,5,current_phase)*100/Property_subdomains_statistics(group_subdomain,4,current_phase);
        Property_subdomains_statistics(group_subdomain,7,current_phase)=min(value_tmp);
        Property_subdomains_statistics(group_subdomain,8,current_phase)=max(value_tmp);
        Property_subdomains_statistics(group_subdomain,9,current_phase)=AR(1);
        Property_subdomains_statistics(group_subdomain,10,current_phase)=AR(2);
        Property_subdomains_statistics(group_subdomain,11,current_phase)=AR(3);
        
    end
    % Subdomains are sorted with size
    Property_subdomains_statistics (:,:,current_phase) = sortrows(Property_subdomains_statistics(:,:,current_phase),1);
end
% Size are set in micrometer
Property_subdomains_statistics(:,1,:)=Property_subdomains_statistics(:,1,:)*voxel_size/1000;
Property_subdomains_statistics(:,2,:)=Property_subdomains_statistics(:,2,:)*voxel_size/1000;

% The representative volume element is reached when the standard deviation
% is lower or equal to a RVEparameters.threshold_std percent of the mean value
% and calculated with a number of subdomains >= RVEparameters.threshold_numbersubvolumes
Size_RVE=zeros(3,number_phase,2);
% Row 1 : RVE size inferior to
% Row 2 : RVE size equal to
% Row 3 : RVE size is superior to
% *,*,1 with size1
% *,*,2 with size2
% If RVE size (row 2) can't be defined, than inferior or superior limit are defined
for k_size=1:1:2
    for current_phase=1:1:number_phase
        RVE_found=0;
        for k=1:1:(number_group_size-1)
            if Property_subdomains_statistics(k,3,current_phase)>1 % 1 domain case (type  'D', 'E')
                val1 = Property_subdomains_statistics(k,6,current_phase); % Relative std
                val2 = Property_subdomains_statistics(k+1,6,current_phase);
                n1 = Property_subdomains_statistics(k,3,current_phase); % Number of subdomins
                n2 = Property_subdomains_statistics(k+1,3,current_phase);
                s1 = Property_subdomains_statistics(k,k_size,current_phase); % size
                s2 = Property_subdomains_statistics(k+1,k_size,current_phase);
                % RVE is inferior to all the subdomains size analysed
                if val1 <RVEparameters.threshold_std && n1 >=RVEparameters.threshold_numbersubvolumes && RVE_found==0
                    Size_RVE(1,current_phase,k_size)=s1;
                    RVE_found=1;
                end
                % RVE size is founded
                if val1 ==RVEparameters.threshold_std && n1 >=RVEparameters.threshold_numbersubvolumes && RVE_found==0
                    Size_RVE(2,current_phase,k_size)=s1;
                    RVE_found=1;
                end
                % RVE size is founded
                if val1 >=RVEparameters.threshold_std && val2 <=RVEparameters.threshold_std && n1 >=RVEparameters.threshold_numbersubvolumes && n2 >=RVEparameters.threshold_numbersubvolumes && RVE_found==0
                    Size_RVE(2,current_phase,k_size)= interp1([val1 val2],[s1 s2],RVEparameters.threshold_std,'linear');
                    RVE_found=1;
                end
            end
        end
        % RVE size is superior to all subdomains analysed
        if RVE_found==0
            for k=1:1:number_group_size
                if Property_subdomains_statistics(k,3,current_phase)>=RVEparameters.threshold_numbersubvolumes
                    Size_RVE(3,current_phase,k_size)=Property_subdomains_statistics(k,k_size,current_phase);
                end
            end
        end
    end
end

end

