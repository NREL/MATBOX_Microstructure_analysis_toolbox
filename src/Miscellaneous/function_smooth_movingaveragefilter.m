function [smoothed_x, smoothed_y, outcome] = function_smooth_movingaveragefilter(x, y, parameters)
%function_smooth_movingaveragefilter Moving average filter
%   [smoothed_x, smoothed_y, outcome] = function_smooth_movingaveragefilter(x,y,parameters)
%   or
%   [smoothed_x, smoothed_y, outcome] = function_smooth_movingaveragefilter(x,y)
%   inputs
%      x: 1D array of length n
%      y: 1D array of length n
%   outputs
%      smoothed_x: 1D array of length m
%      smoothed_y: 1D array of length m
%      outcome: true if algorithm has been done successfully, false otherwise.
%   optional input
%      parameters: structure
%      parameters.minimum_array_length: do not apply the smoothing algorithm is n <= to this value
%      parameters.number_point: Input array are interpolated in a uniformly spaced array (i.e., linspace(x(1),x(end),number_point)
%                               and then the moving average algorithm is applied.
%      parameters.enforce_samelength: true or false
%                                     if true m=n (i.e., smoothed_y=interp1(smoothed_x,smoothed_y,x) and smoothed_x=x;
%                                     if false m=number_point
%      r point moving average filter: r is defined using one of the two arguments below
%      parameters.moving_range: r=parameters.moving_range. Set parameters.moving_range=0 if you want to use the other argument below.
%      parameters.moving_rangeratio: r=round(parameters.moving_rangeratio*number_point), with 0 < parameters.moving_rangeratio < 1.
%                                    Set parameters.moving_rangeratio=0 if you want to use the other argument above.
%      parameters.origin: 'start', 'end', or 'symmetrical'  
%                         'start': for a r point moving average filter: y(k) = (y(k)+y(k+1)+...+y(k+r-1))/r
%                         'end': for a r point moving average filter: y(k) = (y(k-r+1)+...+y(k-1)+y(k))/r
%                         'symmetrical': for a r point moving average filter: y(k) = (y(k-(r-1)/2)+...+y(k-1)+y(k)+y(k+1)+...+y(k+(r-1)/2))/r
%                                        If r is an even number, then r=r+1 to be odd (only for 'symmetrical')
%      parameters.boundary_behavior: 'asymmetrical' or 'keep origin'
%                                    'asymmetrical': For points closed to the boundaries, an asymetrical moving average filter is used to keep constant the number of points used.
%                                                    E.g. for a 5 point moving average filter: y(2) = (y(1)+y(2)+y(3)+y(4)+y(5))/5
%                                    'keep origin' : For points closed to the boundaries, the range r is reduced to keep the same behaviror (start, end or symmetrical).
% Author: Francois L.E. Usseglio-Viretta, National Renewable Energy Laboratory (NREL)

%% PARAMETERS
s_warning = warning; % Save the current warning settings.
warning('on') % Enable all warnings.
switch nargin
    case 2 % argument "parameters" is missing
        % Use default parameters
        minimum_array_length = 5;
        number_point = length(x);
        moving_range = 0;
        moving_rangeratio = 0.05;
        enforce_samelength = false;
        origin = 'symmetrical';
        boundary_behavior = 'keep origin';
    case 3
        if isstruct(parameters) % Check parameters is a structure, then check all fields exist within it.
            if isfield(parameters,'minimum_array_length') &&...
                    isfield(parameters,'number_point')  &&...
                    isfield(parameters,'moving_rangeratio')  &&...
                    isfield(parameters,'moving_range')  &&...
                    isfield(parameters,'enforce_samelength')  &&...
                    isfield(parameters,'origin') &&...
                    isfield(parameters,'boundary_behavior')
                minimum_array_length = parameters.minimum_array_length;
                number_point = parameters.number_point;
                moving_rangeratio = parameters.moving_rangeratio;
                moving_range = parameters.moving_range;
                enforce_samelength = parameters.enforce_samelength;
                origin = parameters.origin;
                boundary_behavior = parameters.boundary_behavior;
            else
                error('function_smooth_movingaveragefilter: 3rd argument ''parameters'' must be a structure. Type help function_smooth_movingaveragefilter for details.')
            end
        else
            error('function_smooth_movingaveragefilter: 3rd argument ''parameters'' must be a structure. Type help function_smooth_movingaveragefilter for details.')
        end
    otherwise
        error('function_smooth_movingaveragefilter: wrong number of arguments (2 or 3). Type help function_smooth_movingaveragefilter for details.')
end

%% CHECK PARAMETERS TYPE, BOUNDS and DIMENSION
% x and y
if ~isfloat(x) || ~isfloat(y)
    error('function_smooth_movingaveragefilter: 1st and 2nd arguments ''x'' and ''y'' must be a 1D array of same length.')
else
    if (isvector(x) && isvector(y) && numel(x) == numel(y))==false
        error('function_smooth_movingaveragefilter: 1st and 2nd arguments ''x'' and ''y'' must be a 1D array of same length.')
    end
end

% minimum_array_length
if ~isfloat(minimum_array_length) || numel(minimum_array_length)~=1 || round(minimum_array_length)~=minimum_array_length || minimum_array_length<1 
    error('function_smooth_movingaveragefilter: 3rd argument ''parameters.minimum_array_length'' must be a positive integer.')
end
if length(x)<=minimum_array_length
    warning('function_smooth_movingaveragefilter: no smoothing performed as the array is too short.')
    smoothed_x = x; smoothed_y = y; outcome = false; warning(s_warning); return % Failure: no smoothing performed (output=input), but program continues.
end

% number_point
if ~isfloat(number_point) || numel(number_point)~=1 || round(number_point)~=number_point || number_point<1 
    error('function_smooth_movingaveragefilter: 3rd argument ''parameters.number_point'' must be a positive integer.');
end

% origin
if (ischar(origin) && (strcmp(origin,'start') || strcmp(origin,'end') || strcmp(origin,'symmetrical')))==false
    error('function_smooth_movingaveragefilter: 3rd argument ''parameters.origin'' must be ''start'' or ''end'' or ''symmetrical''.');
end

% boundary_behavior
if (ischar(boundary_behavior) && (strcmp(boundary_behavior,'asymmetrical') || strcmp(boundary_behavior,'keep origin')))==false
    error('function_smooth_movingaveragefilter: 3rd argument ''parameters.boundary_behavior'' must be ''asymmetrical'' or ''keep origin''.');
end

% enforce_samelength
if ~islogical(enforce_samelength)
    error('function_smooth_movingaveragefilter: 3rd argument ''parameters.enforce_samelength'' must be true or false.');
end

% moving_rangeratio and  moving_range
if ~isfloat(moving_rangeratio) || ~isfloat(moving_range) || numel(moving_rangeratio)~=1 || numel(moving_range)~=1
    error('function_smooth_movingaveragefilter: 3rd argument ''parameters.moving_rangeratio'' and ''parameters.moving_range'' must be numbers.')
else
    if moving_rangeratio~=0 && moving_range~=0
        error('function_smooth_movingaveragefilter: 3rd argument ''parameters.moving_rangeratio'' and ''parameters.moving_range'' cannot be both ~=0. Choose only one ~=0.')
    else
        if moving_range~=0
            if moving_range>=1 && round(moving_range)==moving_range
                r = moving_range;
            else
                error('function_smooth_movingaveragefilter: 3rd argument ''parameters.moving_range'' must be an integer >0. (or =0 if user chooses ''parameters.moving_rangeratio'' instead)')
            end
        elseif moving_rangeratio~=0
            if moving_rangeratio>0 && moving_rangeratio<1
                r = round(moving_rangeratio*number_point); % Ratio of the array length
            else
                error 'function_smooth_movingaveragefilter: 3rd argument ''parameters.moving_rangeratio'' must be a real >0, <1. (or =0 if user chooses ''parameters.moving_range'' instead)')
            end
        end
    end
end
r=max(r,2); % Enforce a minimum value
if strcmp(origin,'symmetrical') && mod(r,2)==0
    r=r+1; % even -> odd
end
if r>=number_point
    warning('function_smooth_movingaveragefilter: no smoothing performed as the range is too high.')
    smoothed_x = x; smoothed_y = y; outcome = false; warning(s_warning); return % Failure: no smoothing performed (output=input), but program continues.
end    

%% ALGORITHM
% x,y: initial array, length=n
% xx,yy: interpolated arrays, uniform spacing, length=number_point
% xx,yyy: averaged arrays, length=number_point
xx = linspace(x(1),x(end),number_point)'; % Uniformly spaced array (each point has the same weight for the subsequent average caluclation)
yy=interp1(x,y,xx); % Interpolate
yyy=zeros(number_point,1); % Initialize

% Moving average filter
switch origin
    case  'start'
        sum_=0; % Initialize
        for k=0:1:(r-1)
            sum_ = sum_ + yy(1+k : number_point-r+1+k);
        end
        yyy(1:number_point-r+1) = sum_/r; % Average value
        if strcmp(boundary_behavior,'keep origin')
            for k=number_point-r+2:number_point
                yyy(k) = sum(yy(k:end)) / (number_point-k+1);
            end
        elseif strcmp(boundary_behavior,'asymmetrical')
            yyy(number_point-r+2:number_point) = sum(yy(number_point-r+1:end)) / r;
        end
        
    case 'end'
        sum_=0; % Initialize
        for k=0:1:(r-1)
            sum_ = sum_ + yy(r-k : number_point-k);
        end
        yyy(r:number_point) = sum_/r; % Average value
        if strcmp(boundary_behavior,'keep origin')
            for k=r-1:-1:1
                yyy(k) = sum(yy(1:k)) / k;
            end
        elseif strcmp(boundary_behavior,'asymmetrical')
            yyy(1:r-1) = sum(yy(1:r)) / r;
        end

    case 'symmetrical'
        r_ = (r-1)/2;
        sum_= yy(r_+1:end-r_); % Initialize
        for k=0:1:(r_-1)
            sum_ = sum_ +  yy(1+k : number_point-r+1+k);
        end
        for k=0:1:(r_-1)
            sum_ = sum_ +  yy(r-k : number_point-k);
        end        
        yyy(r_+1:end-r_) = sum_/r; % Average value
        
        
        if strcmp(boundary_behavior,'keep origin')
            for k=1:1:r_
                r__ = k-1;
                yyy(k) = sum(yy(k-r__:1:k+r__))/(2*r__+1);
            end
            for k=number_point-r_+1:1:number_point
                r__ = number_point-k;
                yyy(k) = sum(yy(k-r__:1:k+r__))/(2*r__+1);                
            end
        elseif strcmp(boundary_behavior,'asymmetrical')
            yyy(number_point-r+2:number_point) = sum(yy(number_point-r+1:end)) / r;
            yyy(1:r-1) = sum(yy(1:r)) / r;
        end
end

if enforce_samelength
    smoothed_x = x;
    smoothed_y = interp1(xx,yyy,x); % Interpolate back to the initial x array
else
    smoothed_x = xx;
    smoothed_y = yyy;
end
 
outcome = true; % Success
warning(s_warning) % Restore the saved warning state structure

end

