function [onset,offset] = onset_offset_extraction(EMG,delta_t,start_swing,start_stance)
%ONSET_OFFSET EXTRACTION Summary of this function goes here
%   Detailed explanation goes here
    mean_swing = mean(EMG(start_swing:start_stance,1));
    sd_swing = std(EMG(start_swing:start_stance,1));
    threshold = mean_swing+2*sd_swing;
    offset = [];
    onset = [];
    for i=1:length(EMG)-1
        if (EMG(i,1)<=threshold & EMG(i+1,1) >= threshold)
            if(isempty(onset))
                onset=[onset i+1];
            elseif (i-onset(end)>delta_t)
                onset=[onset i+1];
            end
        end
        if (EMG(i,1)>= threshold & EMG(i+1,1) <= threshold )
            if(isempty(offset))
                offset=[offset i];
            elseif (i-offset(end)>delta_t)
                offset=[offset i];
            end
        end
    end
end