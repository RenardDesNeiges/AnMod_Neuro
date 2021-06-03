%% Script for the analysis of H01_TDM_2kmh.mat

%clear workspace, load the data
clear 
load('H01_TDM_2kmh.mat') % load the dataset
s = struct; %feature structure

%% Gait events detection

marker_sr = data.marker_sr;

% Right:
hip_r = data.RHIP;
knee_r = data.RKNE;
ankle_r = data.RANK;

% Left:
hip_l = data.LHIP;
knee_l = data.LKNE;
ankle_l = data.LANK;

[cycle_r,cycle_index_r,cycle_time_r] = get_cycle(ankle_r(:,2),hip_r(:,2),marker_sr);
[cycle_l,cycle_index_l,cycle_time_l] = get_cycle(ankle_l(:,2),hip_l(:,2),marker_sr);

avg_cycle_time = mean([cycle_time_r,cycle_time_l]);
var_cycle_time = var([cycle_time_r,cycle_time_l]);

%% setting structure parameters 
s.avg_cycle_time = avg_cycle_time;
s.var_cycle_time = var_cycle_time;
s.velocity = 2/3.6;
