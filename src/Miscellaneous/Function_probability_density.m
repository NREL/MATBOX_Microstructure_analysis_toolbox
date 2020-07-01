function [results, outcome] = Function_probability_density(x,w,parameters)
%Function_probability_density Calculate cumulative function and probability density function
%   [results, outcome] = Function_probability_density(x,w,parameters)
%   or
%   [results, outcome] = Function_probability_density(x,w)
%   inputs
%      x: 1D array of length n (value)
%      w: 1D array of length n (weight), or [].
%   outputs
%      results: structure
%      results.cumulative_fct: 2D array
%      results.probability_density_fct: 2D array, calculated from cumulative_fct
%      results.integral_probability_density_fct: scalar equal to trapz(probability_density_fct(:,1),probability_density_fct(:,2)). It should be theoritically equal to 1.
%      results.smoothed_cumulative_fct: 2D array
%      results.smoothed_probability_density_fct: 2D array, calculated from smoothed_cumulative_fct
%      results.integral_smoothed_probability_density_fct: scalar equal to trapz(smoothed_probability_density_fct(:,1),smoothed_probability_density_fct(:,2)). It should be theoritically equal to 1.
%      outcome: true if algorithm has been done successfully, false otherwise.
%   optional input
%      parameters: structure
%      parameters.round_value: integer >=0. If ~=0, x =round(x,round_value) in order to reduce the number of unique values of x. If ==0, x is unchanged.
%      parameters.smooth_cumulative_fct: true or false. If true, a moving average filter is applied to the cumulative function and
%                                        output results.smoothed_cumulative_fct, results.smoothed_probability_density_fct, and
%                                        results.integral_smoothed_probability_density_fct are calculated with the function function_smooth_movingaveragefilter
%                                        for which additional parameters are required.
%                                        If false, these three outputs are empty.
% Author: Francois L.E. Usseglio-Viretta, National Renewable Energy Laboratory (NREL)

%% PARAMETERS
s_warning = warning; % Save the current warning settings.
warning('on') % Enable all warnings.
switch nargin
    case 2 % argument "parameters" is missing
        % Use default parameters
        round_value = 0;
        smooth_cumulative_fct = false;
    case 3
        if isstruct(parameters) % Check parameters is a structure, then check all fields exist within it.
            if isfield(parameters,'round_value') &&...
                    isfield(parameters,'smooth_cumulative_fct')
                round_value = parameters.round_value;
                smooth_cumulative_fct = parameters.smooth_cumulative_fct;
            else
                error('Function_probability_density: 3rd argument ''parameters'' must be a structure. Type help Function_probability_density for details.')
            end
        else
            error('Function_probability_density: 3rd argument ''parameters'' must be a structure. Type help Function_probability_density for details.')
        end
    otherwise
        error('Function_probability_density: wrong number of arguments (2 or 3). Type help Function_probability_density for details.')
end

%% CHECK PARAMETERS TYPE, BOUNDS and DIMENSION
% x and w
if ~isfloat(x) || ~isfloat(w)
    error('Function_probability_density: 1st argument ''x'' must be a 1D array (value). 2nd argument ''w'' must be a 1D array (weight) or [] (no weight)')
else
    if (isempty(w) && isvector(x) || (numel(x) == numel(w) && isvector(x) && isvector(w)))==false
        error('Function_probability_density: 1st argument ''x'' must be a 1D array (value). 2nd argument ''w'' must be a 1D array (weight) or [] (no weight)')
    end
end

% round_value
if ~isfloat(round_value) || numel(round_value)~=1 || round(round_value)~=round_value || round_value<0
    error('Function_probability_density: 1st argument ''round_value'' must be a positive integer or 0.')
end

% smooth_cumulative_fct
if ~islogical(smooth_cumulative_fct)
    error('Function_probability_density: 2nd argument ''parameters.smooth_cumulative_fct'' must be true or false.');
end

%% CALCULATE THE SUM FUNCTION
results.smoothed_cumulative_fct = []; % Initialize
results.smoothed_x50 = [];
results.smoothed_probability_density_fct = [];
results.integral_smoothed_probability_density_fct = [];

if round_value~=-1
    x = round(x,round_value); % Round the array, to reduce the number of unique values
end
if isempty(w)
    w = x*0+1; % Weight 1. Ensure x and w have same shape.
end
[results.cumulative_fct, results.x50] = function_calculate_cumulative_fct (x, w);
if smooth_cumulative_fct
    % Step 3: smooth cummulative function
    if isfield(parameters,'minimum_array_length') &&...
            isfield(parameters,'number_point')  &&...
            isfield(parameters,'moving_rangeratio')  &&...
            isfield(parameters,'moving_range')  &&...
            isfield(parameters,'enforce_samelength')  &&...
            isfield(parameters,'origin') &&...
            isfield(parameters,'boundary_behavior')
        [smoothed_x, smoothed_y, outcome_smoothing] = function_smooth_movingaveragefilter(results.cumulative_fct(:,1),results.cumulative_fct(:,2),parameters);
    else
        [smoothed_x, smoothed_y, outcome_smoothing] = function_smooth_movingaveragefilter(results.cumulative_fct(:,1),results.cumulative_fct(:,2)); % Use default parameters of function_smooth_movingaveragefilter
    end
    if outcome_smoothing
        results.smoothed_cumulative_fct = [smoothed_x smoothed_y];
        results.smoothed_x50 = interp1(results.smoothed_cumulative_fct(:,2),results.smoothed_cumulative_fct(:,1),0.5);
    end
end

%% CALCULATE THE PROBABILITY DENSITY FUNCTION

[results.probability_density_fct, results.integral_probability_density_fct] = function_calculate_probability_density_fct(results.cumulative_fct(:,1),results.cumulative_fct(:,2));
if ~isempty(results.smoothed_cumulative_fct)
    [results.smoothed_probability_density_fct, results.integral_smoothed_probability_density_fct] = function_calculate_probability_density_fct(results.smoothed_cumulative_fct(:,1),results.smoothed_cumulative_fct(:,2));
else
    results.smoothed_probability_density_fct = [];
    results.integral_smoothed_probability_density_fct = [];
end

outcome = true; % Success
warning(s_warning) % Restore the saved warning state structure


%% FUNCTIONS
    function [cumulative_fct, x50] = function_calculate_cumulative_fct (x,w)
        
        unique_ = unique(x); % Unique values
        cumulative_fct = zeros(length(unique_),2); % Initialisation (x, cumulative fct, smoothed cumulative fct)
        cumulative_fct(:,1)=unique_;        
        % Step 1: probability randomly chosen point p from the array a of lenght n is equal with the unique values v
        %         p(p==v)=1/n*sum(a==v),
        pdi_=zeros(length(unique_),2);
        pdi_(:,1)=unique_;
        all_weight = unique(w);
        if length(all_weight)==1 && all_weight==1
            number_element=numel(x);
            for current_value=1:1:length(unique_)
                pdi_(current_value,2)= sum(sum(sum( x==unique_(current_value) ))) / number_element;
            end
        else
            total_weight = sum(w);
            for current_value=1:1:length(unique_)
                idx = find(x==unique_(current_value));
                pdi_(current_value,2)= sum(sum(sum( w(idx) ))) / total_weight;
            end            
        end
        % Step 2: The sum function P(D)= 1-(P(0<=d<=D)=sum(p(di)), di<=D)
        %for current_value=1:1:length(unique_)
        %    cumulative_fct(current_value,2)=sum(pdi_(current_value:end,2));
        %end
        % Much more faster: start form the last value
        cumulative_fct(end,2)=pdi_(end,2);
        for current_value=length(unique_)-1:-1:1
            cumulative_fct(current_value,2)=cumulative_fct(current_value+1,2)+pdi_(current_value,2);
        end
        % Find x50.
        x50 = interp1(cumulative_fct(:,2),cumulative_fct(:,1),0.5);
    end

    function [probability_density_fct, integral_pdf] = function_calculate_probability_density_fct (x, cumulative_fct)
        n_=length(x);
        probability_density_fct = zeros(n_,2); % Initialize
        integral_pdf=0;
        if n_>1
            f_=cumulative_fct; % For visibility sake
            for current_=2:1:n_-1
                probability_density_fct(current_,2)= (f_(current_+1)-f_(current_-1))/(x(current_+1)-x(current_-1));
            end
            probability_density_fct(1,2)= (f_(2)-f_(1))/(x(2)-x(1));
            probability_density_fct(n_,2)= (f_(n_)-f_(n_-1))/(x(n_)-x(n_-1));
            probability_density_fct(:,2)=-probability_density_fct(:,2);
            probability_density_fct(:,1) = x;
            integral_pdf = trapz(probability_density_fct(:,1),probability_density_fct(:,2)); % Probability density function integral (should be equal to 1)
        end
    end
end