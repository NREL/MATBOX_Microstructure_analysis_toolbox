function [histogram_array] = function_calculate_histogram(array)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% colunm 1: values
% colunm 2: sum(array)=value
% colunm 3: normalized (i.e. colunm 2/max(colunm 2))

[tmp, ~, ic]=unique(array); % Get all values of the array
n=length(tmp); % Number of unique value
histogram_array=zeros(n,3); % Initialise histogram array
histogram_array(:,1) = tmp; % Column 1

% colunm 2: histogram sum(array)=value
%for k=1:1:n
%    histogram_array(k,2)=sum(sum(sum( array==unique_array(k,1) )));
%end
histogram_array(:,2) = accumarray(ic,1); % Faster

% colunm 3: normalized (i.e. colunm 2/max(colunm 2))
histogram_array(:,3)=histogram_array(:,2)/max(histogram_array(:,2));

end

