function [max_height,max_times,max_indices] = max_heights(ankle_height,threshold,sr)
%max_heights Detects max heights in the walking cycle (for ankle height)
%   takes :
%   - ankle_height: raw ankle height from mocap
%   - threshold : threshold under which we do not consider the point a max
%   - sr: sample rate
%   returns :
%   - maxs_height a vector containing the max height for each max
%   - max_times a vector containing time of each maximum
%   - max_indices a vector containing the indices of each maximum
    
    % take derivative
    ankle_speed = diff(ankle_height);
    
    %remove the noise from the derivative
    Low = 15/(sr/2);
    [b1, a1] = butter(2,Low,'low'); 
    ankle_speed_flt = filter(b1, a1, ankle_speed);
    
    % find zero-crossings of the filtered derivative
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    raw_zero_indices_speed = zci(ankle_speed_flt);
    height_when_diff_0 = ankle_height(raw_zero_indices_speed);
    
    % only consider high enough values
    max_indices = raw_zero_indices_speed(height_when_diff_0 > threshold);
    max_times = max_indices * (1/sr);
    max_height = ankle_height(max_indices);
end

