%% Script for the analysis of H01_TDM_2kmh.mat

%clear workspace, load the data
addpath(genpath('..'))
clear 
load('H01_TDM_2kmh.mat') % load the dataset
s = struct; %feature structure

%% Gait events detection
marker_sr = data.marker_sr;

% Right:
hip_r = data.RHIP;
knee_r = data.RKNE;
ankle_r = data.RANK;
toe_r = data.RTOE;

% Left:
hip_l = data.LHIP;
knee_l = data.LKNE;
ankle_l = data.LANK;
toe_l = data.LTOE;

[cycle_r,cycle_index_r,cycle_time_r] = get_cycle(ankle_r(:,2),hip_r(:,2),marker_sr);
[cycle_l,cycle_index_l,cycle_time_l] = get_cycle(ankle_l(:,2),hip_l(:,2),marker_sr);

avg_cycle_time = mean([cycle_time_r,cycle_time_l]);
var_cycle_time = var([cycle_time_r,cycle_time_l]);
%%
clf
clc

[pitch_foot_angle,pitch_angular_velocity] = ...
        foot_pitch_vel(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));
plot(pitch_foot_angle()*0.04)
hold on 


[stance_starts_r,swing_starts_r] = ...
        swing_stance(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));
    
[stance_starts_l,swing_starts_l] = ...
        swing_stance(toe_l(:,2),ankle_l(:,2),toe_l(:,3),ankle_l(:,3));

plot(stance_starts_r,0,'or')
plot(swing_starts_r,0,'ok')

%% setting structure parameters 
s.avg_cycle_time = avg_cycle_time;  % in seconds
s.var_cycle_time = var_cycle_time;  % in seconds
s.velocity = 2/3.6;                 % in meter/second
