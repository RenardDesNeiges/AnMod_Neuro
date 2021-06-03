function [cycle,cycle_index,cycle_time] = get_cycle(ankle,hip,sr)
%GET_CYCLE Returns the time of each walking cycle for a single leg's ankle
%and hip y coordinates
%   Detailed explanation goes here
    normalized_y_ankle = ankle-hip;
    normalized_y_ankle_speed = diff(normalized_y_ankle);
    
    %filter the  finite differences vector to remove noise
    Low = 7/(sr/2);
    [b1, a1] = butter(2,Low,'low'); 
    normalized_y_ankle_speed_flt = filter(b1, a1, normalized_y_ankle_speed);
    
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

    zeros_ankle = (zci(normalized_y_ankle));
    zeros_ankle = zeros_ankle(1:end-1);
    
    derivative_when_zero = normalized_y_ankle_speed_flt(zeros_ankle);
    cycle_index = zeros_ankle( derivative_when_zero < 0);
    cycle = (1/sr) * cycle_index;
    
    cycle_time = zeros(size(cycle)-1);
    
    for i = 2:1:(size(cycle)-1)
        cycle_time(i) = cycle(i) - cycle(i-1);
    end
end

