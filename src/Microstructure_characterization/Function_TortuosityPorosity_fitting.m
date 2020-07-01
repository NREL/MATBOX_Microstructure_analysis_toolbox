function [results] = Function_TortuosityPorosity_fitting(epsilon,tau, form)
% [results] = Function_TortuosityPorosity_fitting(epsilon,tau, form)
% or
% [results] = Function_TortuosityPorosity_fitting(epsilon,tau)
% - form  = 'Generalized_ArchieLaw' (by default if form is not specified)
% find alpha and gamma so that tau = alpha * epsilon^(-gamma)
% results.alpha and results.gamma
% - form  = 'Empirical_Bruggeman' (by default if function argument epsilon,tau a unique doublet)
% find gamma and k so that tau = epsilon^(-gamma)
% results.gamma

results.success = false; % Initialize

nargin; % Number of input variable when the function is call
if nargin == 2 
    form = 'Generalized_ArchieLaw';
end

% Make sure data are one-column array
epsilon = reshape(epsilon,[numel(epsilon),1]);
tau = reshape(tau,[numel(tau),1]);
    
% Check input have correct dimension
n =length(epsilon);
n_tmp = length(tau);
if n ~= n_tmp
    warning('error in [results] = Function_TortuosityPorosity_fitting(epsilon,tau, form) : epsilon and tau do not have the same length')
    return
end
clear n_tmp

% Check for NaN
tmp1 = find(isnan(epsilon)); % Should never happen
tmp2 = find(isnan(tau)); % Can happen if there is no percolation path
tmp = unique([tmp1; tmp2]);
if ~isempty(tmp)
    epsilon(tmp)=[];
    tau(tmp)=[];
end
clear tmp1 tmp2 tmp

% Check for 0 (to prevent log(0) error)
tmp1 = find(epsilon==0); % Should never happen
tmp2 = find(tau==0); % Should never happen
tmp = unique([tmp1; tmp2]);
if ~isempty(tmp)
    epsilon(tmp)=[];
    tau(tmp)=[];
end
clear tmp1 tmp2 tmp

% Calculate
n = length(epsilon);
if n==0 % Empty data
    warning('error in [results] = Function_TortuosityPorosity_fitting(epsilon,tau, form) : epsilon and tau are empty, are NaN, or 0')
    return
elseif n==1 % One unique doublet (epsilon, tau): enforce tau = epsilon^(-gamma)
    results.gamma = -log(tau)/log(epsilon);
    results.form = 'Empirical_Bruggeman';
    results.equation = 'tau = epsilon^(-gamma)';
    results.success = true;
else
    results.form = form;
    if strcmp(form, 'Empirical_Bruggeman')
        p = polyfit(log(epsilon),log(tau),1);
        results.gamma = -p(1);
        results.equation = 'tau = epsilon^(-gamma)';
        results.epsilon_fit = linspace(min(epsilon), max(epsilon));
        results.tau_fit = results.epsilon_fit.^(-results.gamma);        
        results.success = true;
    elseif strcmp(form, 'Generalized_ArchieLaw')
        p = polyfit(log(epsilon),log(tau),1);
        results.alpha = exp(p(2));
        results.gamma = -p(1);
        results.equation = 'tau = alpha * epsilon^(-gamma)';
        results.epsilon_fit = linspace(min(epsilon), max(epsilon));
        results.tau_fit = results.alpha * results.epsilon_fit.^(-results.gamma);
        results.success = true;
    end
end

end

