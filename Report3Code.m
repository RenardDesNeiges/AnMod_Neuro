%%
addpath(genpath('..'))
clear
% Loading the data

load('H01_TDM_2kmh.mat')


%%
%EMG signals preparation
% extensor : sol, ST , VLAT,MG
% flexor : TA ,II, RF
sr = data.EMG_sr;
tmin = 0;
tmax = 242900/sr;
t = linspace(tmin,tmax,242900);

Flexors_left = [data.LTA, data.LIl, data.LRF];
Flexors_right = [data.RTA, data.RIl, data.RRF];
Extensors_left = [data.LSol, data.LST, data.LVLat, data.LMG];
Extensors_right = [data.RSol, data.RST, data.RVLat, data.RMG];

[r_F,c_F]=size(Flexors_left);
Flexors_left_filtered = [];
Flexors_right_filtered = [];

for i=1:c_F
    inter_l = Filter_EMG(Flexors_left(:,i),sr);
    Flexors_left_filtered= [Flexors_left_filtered, inter_l];
    inter_r = Filter_EMG(Flexors_right(:,i),sr);
    Flexors_right_filtered= [Flexors_right_filtered, inter_r];
end


[r_E,c_E]=size(Extensors_left);
Extensors_left_filtered = [];
Extensors_right_filtered = [];

for i=1:c_E
    inter_l = Filter_EMG(Extensors_left(:,i),sr);
    Extensors_left_filtered= [Extensors_left_filtered, inter_l];
    inter_r = Filter_EMG(Extensors_right(:,i),sr);
    Extensors_right_filtered= [Extensors_right_filtered, inter_r];
end


%% Filtered muscle EMG plots
plot (t,Flexors_left_filtered(:,:));
% IL et RF inutiles (activity ~ 0)
figure;
plot (t,Extensors_left_filtered(:,:));
% ST inutile (activity ~ 0)
%% muscle activity params extraction
%Flexors
[r_F,c_F]=size(Flexors_left);
Flexors_MEAN_left = [];
Flexors_RMS_left = [];
Flexors_integral_left = [];
Flexors_MEAN_right = [];
Flexors_RMS_right = [];
Flexors_integral_right = [];

for i=1:c_F
    [mean_l,RMS_l,integ_l] = Extract_muscle_features(Flexors_left_filtered(:,i));
    Flexors_MEAN_left= [Flexors_MEAN_left, mean_l];
    Flexors_RMS_left= [Flexors_RMS_left, RMS_l];
    Flexors_integral_left= [Flexors_integral_left, integ_l];
    [mean_r,RMS_r,integ_r] = Extract_muscle_features(Flexors_right_filtered(:,i));
    Flexors_MEAN_right= [Flexors_MEAN_right, mean_r];
    Flexors_RMS_right= [Flexors_RMS_right, RMS_r];
    Flexors_integral_right= [Flexors_integral_right, integ_r];
    
end

%Extensors
[r_F,c_F]=size(Flexors_left);
Extensors_MEAN_left = [];
Extensors_RMS_left = [];
Extensors_integral_left = [];
Extensors_MEAN_right = [];
Extensors_RMS_right = [];
Extensors_integral_right = [];

for i=1:c_F
    [mean_l,RMS_l,integ_l] = Extract_muscle_features(Extensors_left_filtered(:,i));
    Extensors_MEAN_left= [Extensors_MEAN_left, mean_l];
    Extensors_RMS_left= [Extensors_RMS_left, RMS_l];
    Extensors_integral_left= [Extensors_integral_left, integ_l];
    [mean_r,RMS_r,integ_r] = Extract_muscle_features(Extensors_right_filtered(:,i));
    Extensors_MEAN_right= [Extensors_MEAN_right, mean_r];
    Extensors_RMS_right= [Extensors_RMS_right, RMS_r];
    Extensors_integral_right= [Extensors_integral_right, integ_r];
    
end
%%
%One of the most common transformations used is the integration of 
%the absolute values of the amplitudes of the EMG spikes.
%Through this transformation, it has been found that the area 
%under the graph of the absolute integral of the EMG is linearly
%proportional to the strength of the muscle contraction.
plot (t,Flexors_integral_left(:,:));


