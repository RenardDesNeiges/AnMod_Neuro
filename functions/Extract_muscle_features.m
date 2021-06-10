function [mean_amp,RMS,integ] = Extract_muscle_features(EMG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mean_amp = mean(EMG);
RMS = rms(EMG);
integ = cumtrapz(EMG);
end

