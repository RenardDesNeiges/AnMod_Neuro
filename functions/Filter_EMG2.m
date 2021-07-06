function [filtered_EMG] = Filter_EMG2(EMG)
%FILTER_EMG2 Summary of this function goes here
%   Detailed explanation goes here

% Step 1 - Slightly filter (band pass 10-2000 Hz) --> avec notre sampling
% frequency on peut aller jusqu'à 1000
Band  = [0.050, 0.450];
[b1, a1] = butter(2, Band, 'Bandpass');   
signal_step1 = filter(b1, a1, EMG);



% Step 2 - Signal rectification
signal_step2 = abs(signal_step1);



% Step 3 - Low pass filter (10 Hz)
Low = 0.010;
[b3,a3]=butter(2,Low,'low'); 
filtered_EMG = filter(b3,a3,signal_step2);

end

