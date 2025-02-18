function [GL] = function_saturate(GL,Threshold_H,Threshold_L)
% % Threshold will be identified using a Dichotomy approach.
% Faster than compute histogram.
max_=max(max(max(GL)));
min_=min(min(min(GL)));
number_voxel = numel(GL);
max_iteration=100;
target_error=0.1; % in percent
if Threshold_H~=0
    threshold_highervalue=max_;
    error_=1e9; iteration_no=0;
    while error_>target_error && iteration_no<max_iteration % Dichotomy loop
        iteration_no=iteration_no+1; % Update iteration number
        current_threshold=min_+round((max_-min_)/2); % Current threshold
        % Calculate volume
        volume_ratio = 100*sum(sum(sum(GL>current_threshold)))/number_voxel;
        % Check error
        rel_error=volume_ratio-Threshold_H;
        error_=abs(rel_error);
        % Keep going if error is still too high
        if error_<target_error
            % Exit the while loop
            break
        else
            % Update the interval
            if rel_error>0
                min_= current_threshold;
            else
                max_=current_threshold;
            end
        end
    end
    threshold_highervalue=current_threshold;
end
if Threshold_L~=0
    max_=max(max(max(GL)));
    min_=min(min(min(GL)));
    threshold_lowervalue=min_;
    error_=1e9; iteration_no=0;
    while error_>target_error && iteration_no<max_iteration % Dichotomy loop
        iteration_no=iteration_no+1; % Update iteration number
        current_threshold=min_+round((max_-min_)/2); % Current threshold
        % Calculate volume
        volume_ratio = 100*sum(sum(sum(GL<current_threshold)))/number_voxel;
        % Check error
        rel_error=volume_ratio-Threshold_L;
        error_=abs(rel_error);
        % Keep going if error is still too high
        if error_<target_error
            % Exit the while loop
            break
        else
            % Update the interval
            if rel_error>0
                max_= current_threshold;
            else
                min_=current_threshold;
            end
        end
    end
    threshold_lowervalue=current_threshold;
end
% Apply changes
if Threshold_H~=0
    GL(GL>=threshold_highervalue)=threshold_highervalue;
end
if Threshold_L~=0
    GL(GL<=threshold_lowervalue)=threshold_lowervalue;
end

end