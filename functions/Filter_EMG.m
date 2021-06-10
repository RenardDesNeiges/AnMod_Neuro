function [filtered_EMG] = Filter_EMG(EMG,sr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Step 1 - Slightly filter (band pass 10-2000 Hz) --> avec notre sampling
% frequency on peut aller jusqu'à 1000
Band  = [10/(sr/2), 999/(sr/2)];
[b1, a1] = butter(2, Band, 'Bandpass');   
signal_step1 = filter(b1, a1, EMG);



% Step 2 - High pass filter (30 Hz)
High = 10/(sr/2);
[b2,a2]=butter(2,High,'high'); 
signal_step2 = filter(b2,a2,signal_step1);



% Step 3 - Signal rectification
signal_step3 = abs(signal_step2);


% Step 4 - Band stop filter (around 50 Hz)
BandStop = [45/(sr/2), 55/(sr/2)];
[b3,a3]=butter(2,BandStop,'stop');
signal_step4 = filter(b3,a3,signal_step3);


% Step 5 - Low pass filter (10 Hz)
Low = 10/(sr/2);
[b4,a4]=butter(2,Low,'low'); 
filtered_EMG = filter(b4,a4,signal_step4);


end

