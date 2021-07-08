function [onset_id,offset_id,onset_t,offset_t,onset_dt,offset_dt] = ...
    onset_offset_extraction(EMG,delta_t,start_swing,start_stance,sr)
%ONSET_OFFSET EXTRACTION Computes onsets and ofset indices for an EMG
%signal
%   Parameters
%       - EMG : EMG time series
%       - delta_t : expected burst length
    
    mean_swing = mean(EMG(start_swing:start_stance,1));
    sd_swing = std(EMG(start_swing:start_stance,1));
    threshold = mean_swing+2*sd_swing;
    onset_id = [];
    offset_id = [];
    
    for i=1:length(EMG)-1
        if (EMG(i,1)<=threshold & EMG(i+1,1) >= threshold)
            if(isempty(onset_id))
                onset_id=[onset_id i+1];
            elseif (i-onset_id(end)>delta_t)
                onset_id=[onset_id i+1];
            end
        end
        if (EMG(i,1)>= threshold & EMG(i+1,1) <= threshold )
            if(isempty(offset_id))
                offset_id=[offset_id i];
            elseif (i-offset_id(end)>delta_t)
                offset_id=[offset_id i];
            end
        end
    end
    
    onset_t = onset_id * (1/sr);
    offset_t = offset_id * (1/sr);
    
    onset_dt = diff(onset_t);
    offset_dt = diff(offset_t);
end