function [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis2(number_group_size,number_phase_todo,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters,dimension)

%% You have to have the Statistics Toolbox to use the tinv function
% To avoid this requirement, here is a replacement
% Source: https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
% Variables: 
% t: t-statistic
% v: degrees of freedom
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
% This calculates the inverse t-distribution (parameters given the
%   probability ‘alpha’ and degrees of freedom ‘v’: 
t_inv = @(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - tdist1T(tval,v)), 5);  % T-Statistic Given Probability ‘alpha’ & Degrees-Of-Freedom ‘v’

%% STATS
Property_subdomains_statistics = zeros(number_group_size,17,number_phase_todo);
% Each row corresponds to a different ideal size
% column 1 : size 1: equivalent cubic length
% column 2 : size 2: equivalent square root length
% column 3 : size 3: length
% column 4 : number of subdomains with the same size
% column 5 : mean value
% column 6 : standard deviation, only calculated if number of subdomains is >2
% column 7 : standard deviation expressed in % of the mean value
% column 8 : minimun
% column 9 : maximum
% column 10-12: subdomain aspect ratio (AR)
% column 13: Standard Error Of The Mean std(y)/sqrt(N);
% column 14: 95% error interval
% column 15: number of volume required to have 1% error
% column 16: number of volume required to have 5% error
% column 17: number of volume required to have 10% error

for current_phase_todo=1:1:number_phase_todo
    for group_subdomain = 1:1:number_group_size
        AR = GROUP_SUBDOMAIN.id(group_subdomain).Aspectratio;
        % Number of subdomain sharing this size
        number_subdomain_withingroup = length(GROUP_SUBDOMAIN.id(group_subdomain).domain_id);
        value_tmp =[];
        for current_subdomain=1:1:number_subdomain_withingroup
            subdomain_id=GROUP_SUBDOMAIN.id(group_subdomain).domain_id(current_subdomain);
            value_tmp=[value_tmp Property_eachsubdomain(subdomain_id,current_phase_todo+4)];
        end
        Property_subdomains_statistics(group_subdomain,1,current_phase_todo)= GROUP_SUBDOMAIN.id(group_subdomain).equivalent_cubic_length;
        Property_subdomains_statistics(group_subdomain,2,current_phase_todo)= GROUP_SUBDOMAIN.id(group_subdomain).equivalent_square_length;
        Property_subdomains_statistics(group_subdomain,3,current_phase_todo)= GROUP_SUBDOMAIN.id(group_subdomain).length;

        Property_subdomains_statistics(group_subdomain,4,current_phase_todo)=number_subdomain_withingroup;
        Property_subdomains_statistics(group_subdomain,5,current_phase_todo)=nanmean(value_tmp);
        Property_subdomains_statistics(group_subdomain,6,current_phase_todo)=nanstd(value_tmp);
        Property_subdomains_statistics(group_subdomain,7,current_phase_todo)=Property_subdomains_statistics(group_subdomain,6,current_phase_todo)*100/Property_subdomains_statistics(group_subdomain,5,current_phase_todo);
        Property_subdomains_statistics(group_subdomain,8,current_phase_todo)=min(value_tmp);
        Property_subdomains_statistics(group_subdomain,9,current_phase_todo)=max(value_tmp);
        Property_subdomains_statistics(group_subdomain,10,current_phase_todo)=AR(1);
        Property_subdomains_statistics(group_subdomain,11,current_phase_todo)=AR(2);
        if dimension == 3
            Property_subdomains_statistics(group_subdomain,12,current_phase_todo)=AR(3);
        else
            Property_subdomains_statistics(group_subdomain,12,current_phase_todo)=0;
        end
        Property_subdomains_statistics(group_subdomain,13,current_phase_todo)=Property_subdomains_statistics(group_subdomain,6,current_phase_todo)/sqrt(number_subdomain_withingroup);
        % ts = tinv([0.025  0.975],number_subdomain_withingroup-1);      % T-Score. Require Statistics Toolbox
        % Property_subdomains_statistics(group_subdomain,14,current_phase_todo) = abs(ts(1))*Property_subdomains_statistics(group_subdomain,13,current_phase_todo);   
        % Instead:
        if number_subdomain_withingroup>1
            ts = t_inv(0.025,number_subdomain_withingroup-1);
            Property_subdomains_statistics(group_subdomain,14,current_phase_todo) = ts*Property_subdomains_statistics(group_subdomain,13,current_phase_todo); % 2.5% confidence interval
        end
        %Property_subdomains_statistics(group_subdomain,15,current_phase_todo) = Property_subdomains_statistics(group_subdomain,6,current_phase_todo)^2 / (0.01^2 * Property_subdomains_statistics(group_subdomain,5,current_phase_todo)^2);
        %Property_subdomains_statistics(group_subdomain,16,current_phase_todo) = Property_subdomains_statistics(group_subdomain,6,current_phase_todo)^2 / (0.05^2 * Property_subdomains_statistics(group_subdomain,5,current_phase_todo)^2);
        %Property_subdomains_statistics(group_subdomain,17,current_phase_todo) = Property_subdomains_statistics(group_subdomain,6,current_phase_todo)^2 / (0.10^2 * Property_subdomains_statistics(group_subdomain,5,current_phase_todo)^2);

    end
    % Subdomains are sorted with size
    Property_subdomains_statistics (:,:,current_phase_todo) = sortrows(Property_subdomains_statistics(:,:,current_phase_todo),1);
end
% Size are set in pysical unit
Property_subdomains_statistics(:,1,:)=Property_subdomains_statistics(:,1,:)*voxel_size;
Property_subdomains_statistics(:,2,:)=Property_subdomains_statistics(:,2,:)*voxel_size;
Property_subdomains_statistics(:,3,:)=Property_subdomains_statistics(:,3,:)*voxel_size;

% The representative volume element is reached when the standard deviation
% is lower or equal to a RVEparameters.threshold_std percent of the mean value
% and calculated with a number of subdomains >= RVEparameters.threshold_numbersubvolumes
n_threshold = length(RVEparameters.threshold_std);
Size_RVE=zeros(n_threshold,number_phase_todo,3,2);
% *,*,1 RVE size inferior to
% *,*,2 RVE size equal to
% *,*,3 RVE size is superior to
% *,*,*,1: with size1
% *,*,*,2: with size2
% By definition, case D,E,F are not representativity analysis
% If RVE size can't be defined, than inferior or superior limit are defined
if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
    derivative_convergence = [];
    relativedifference_convergence = [];
    Size_convergence = [];
    if dimension == 3
        ktmpsize = [1 1]; % Volume
        if strcmp(RVEparameters.type,'C')
            ktmpsize = [1 2]; % Square root
        elseif strcmp(RVEparameters.type,'D')
            ktmpsize = [1 3]; % Length
        end
    else
        ktmpsize = [2 2]; % Square root
        if strcmp(RVEparameters.type,'C')
            ktmpsize = [2 3]; % Length
        elseif strcmp(RVEparameters.type,'D')
            ktmpsize = [2 3]; % Length
        end
    end

    for k_size=1:1:2
        for k_threshold=1:1:n_threshold
            for current_phase_todo=1:1:number_phase_todo

                % Select value for current phase, and remove analys performed on a too low number of subdomains
                tmp = Property_subdomains_statistics(:,:,current_phase_todo);
                [n_group_size,~]=size(tmp);

                if ~isempty(tmp)
                    idx=find(tmp(:,7) > RVEparameters.threshold_std(k_threshold));
                    if ~isempty(idx)
                        if idx(end)==n_group_size % RVE >
                            Size_RVE(k_threshold,current_phase_todo,1,k_size) = 0;
                            Size_RVE(k_threshold,current_phase_todo,2,k_size) = 0;
                            Size_RVE(k_threshold,current_phase_todo,3,k_size) = tmp(end,ktmpsize(k_size));
                        else % RVE =
                            x1 = tmp(idx(end),ktmpsize(k_size)); x2 = tmp(idx(end)+1,ktmpsize(k_size)); 
                            y1 = tmp(idx(end),7); y2 = tmp(idx(end)+1,7); 
                            x = interp1([y1 y2],[x1 x2],RVEparameters.threshold_std(k_threshold),'linear');
                            Size_RVE(k_threshold,current_phase_todo,1,k_size) = 0;
                            Size_RVE(k_threshold,current_phase_todo,2,k_size) = x;
                            Size_RVE(k_threshold,current_phase_todo,3,k_size) = 0;                            
                        end
                    else % RVE <
                        Size_RVE(k_threshold,current_phase_todo,1,k_size) = tmp(1,ktmpsize(k_size));
                        Size_RVE(k_threshold,current_phase_todo,2,k_size) = 0;
                        Size_RVE(k_threshold,current_phase_todo,3,k_size) = 0;
                    end
                end

            end
        end
    end
else
    derivative_convergence = zeros(number_group_size-2,number_phase_todo+1,2);
    % derivative_convergence(:,1,1) = Property_subdomains_statistics(2:end-1,1,1);
    % derivative_convergence(:,1,2) = Property_subdomains_statistics(2:end-1,3,1);    
    % 
    % if dimension == 3
    %     kefgh = 1; % Volume
    %     kh = 2;
    % else
    %     kefgh = 2; % Square
    %     kh = 3;
    % end
    % for current_phase_todo=1:1:number_phase_todo
    %     for k=2:1:(number_group_size-1)
    %         x1 = Property_subdomains_statistics(k-1,kefgh,current_phase_todo);
    %         x2 = Property_subdomains_statistics(k+1,kefgh,current_phase_todo);            
    %         y1 = Property_subdomains_statistics(k-1,5,current_phase_todo);
    %         y2 = Property_subdomains_statistics(k+1,5,current_phase_todo);
    %         derivative_convergence(k-1,current_phase_todo+1,1) = abs((y2-y1)/(x2-x1));
    %     end
    % end
    % if strcmp(RVEparameters.type,'G')
    %     for current_phase_todo=1:1:number_phase_todo
    %         for k=2:1:(number_group_size-1)
    %             x1 = Property_subdomains_statistics(k-1,3,current_phase_todo);
    %             x2 = Property_subdomains_statistics(k+1,3,current_phase_todo);
    %             y1 = Property_subdomains_statistics(k-1,5,current_phase_todo);
    %             y2 = Property_subdomains_statistics(k+1,5,current_phase_todo);
    %             derivative_convergence(k-1,current_phase_todo+1,2) = abs((y2-y1)/(x2-x1));
    %         end
    %     end
    % end
    % if strcmp(RVEparameters.type,'H')
    %     for current_phase_todo=1:1:number_phase_todo
    %         for k=2:1:(number_group_size-1)
    %             x1 = Property_subdomains_statistics(k-1,kh,current_phase_todo);
    %             x2 = Property_subdomains_statistics(k+1,kh,current_phase_todo);
    %             y1 = Property_subdomains_statistics(k-1,5,current_phase_todo);
    %             y2 = Property_subdomains_statistics(k+1,5,current_phase_todo);
    %             derivative_convergence(k-1,current_phase_todo+1,2) = abs((y2-y1)/(x2-x1));
    %         end
    %     end
    % end    


    relativedifference_convergence = zeros(number_group_size,number_phase_todo+2);
    if dimension == 3
        k1 = 1; % Cubic root
        k2g = 3; % Length
        k2h = 2; % % Square root
    else
        k1 = 2; % Square root
        k2g = 3; % Length
        k2h = 3; % Length
    end

    relativedifference_convergence(:,1) = Property_subdomains_statistics(:,k1,1); % Cubic root
    if strcmp(RVEparameters.type,'G')
        relativedifference_convergence(:,2) = Property_subdomains_statistics(:,k2g,1); % Length
    end
    if strcmp(RVEparameters.type,'H')
        relativedifference_convergence(:,2) = Property_subdomains_statistics(:,k2h,1); % Square root
    end

    for current_phase_todo=1:1:number_phase_todo
        lastvalue = Property_subdomains_statistics(end,5,current_phase_todo);
        for k=1:1:number_group_size
            currentvalue = Property_subdomains_statistics(k,5,current_phase_todo);
            relativedifference_convergence(k,current_phase_todo+2) = abs(100*(currentvalue-lastvalue)/lastvalue);
        end
    end

    % Find convergence threshold
    n_threshold = length(RVEparameters.threshold_reldiff);
    Size_convergence=zeros(n_threshold,number_phase_todo,3,2); % *,*, < = >, size 1 size 2
    for k_threshold=1:1:n_threshold
        for current_phase_todo=1:1:number_phase_todo
            idx=find(relativedifference_convergence(:,current_phase_todo+2) > RVEparameters.threshold_reldiff(k_threshold));
            if ~isempty(idx)
                idx=[idx(end) idx(end)+1];
                if idx(end)+1 == number_group_size
                    Size_convergence(k_threshold,current_phase_todo,3,1)= Property_subdomains_statistics(end,k1,current_phase_todo); % Convergence size is > x
                    if strcmp(RVEparameters.type,'G')
                        Size_convergence(k_threshold,current_phase_todo,3,2)=Property_subdomains_statistics(end,k2g,current_phase_todo); % Convergence size is > x
                    elseif strcmp(RVEparameters.type,'H')
                        Size_convergence(k_threshold,current_phase_todo,3,2)=Property_subdomains_statistics(end,k2h,current_phase_todo); % Convergence size is > x
                    end
                else
                    x1 = relativedifference_convergence(idx(1),1); x2 = relativedifference_convergence(idx(2),1);
                    y1 = relativedifference_convergence(idx(1),current_phase_todo+2); y2 = relativedifference_convergence(idx(2),current_phase_todo+2);
                    x = interp1([y1 y2],[x1 x2],RVEparameters.threshold_reldiff(k_threshold),'linear');
                    Size_convergence(k_threshold,current_phase_todo,2,1)=x; % Convergence size is = x
                    if strcmp(RVEparameters.type,'G') || strcmp(RVEparameters.type,'H')
                        x1 = relativedifference_convergence(idx(1),2); x2 = relativedifference_convergence(idx(2),2);
                        x = interp1([y1 y2],[x1 x2],RVEparameters.threshold_reldiff(k_threshold),'linear');
                        Size_convergence(k_threshold,current_phase_todo,2,2)=x; % Convergence size is = x
                    end
                end
            else
                Size_convergence(k_threshold,current_phase_todo,1,1)=Property_subdomains_statistics(1,k1,current_phase_todo); % Convergence size is < x
                if strcmp(RVEparameters.type,'G')
                    Size_convergence(k_threshold,current_phase_todo,1,2)=Property_subdomains_statistics(1,k2g,current_phase_todo); % Convergence size is < x
                elseif strcmp(RVEparameters.type,'H')
                    Size_convergence(k_threshold,current_phase_todo,1,2)=Property_subdomains_statistics(1,k2h,current_phase_todo);
                end
            end
        end
    end
end

% We can also determine a RVE from the error analysis
% N = 4*std^2 / (relative error^2 * mean^2)
% For N=1 and a given relative error, what should be the std ?
% Knowing std=f(subdomain volume), deduce the RVE size.
% This method is proposed in : T. Kanit et al. / International Journal of Solids and Structures 40 (2003) 3647–3679, page 3672
% To do



end

