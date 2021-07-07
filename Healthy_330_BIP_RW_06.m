%% Script for the analysis of H01_TDM_2kmh.mat

%clear workspace, load the data
addpath(genpath('..'))
clear 
load('Healthy_330_BIP_RW_06.mat') % load the dataset 
velocity = 2; %velocity in km/h
s = struct; %feature structure
%% Gait events detection
marker_sr = data.marker_sr;
emg_sr = data.EMG_sr;

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

% average cycle time
[cycle_r,cycle_index_r,cycle_time_r] = get_cycle(ankle_r(:,2),hip_r(:,2),marker_sr);
[cycle_l,cycle_index_l,cycle_time_l] = get_cycle(ankle_l(:,2),hip_l(:,2),marker_sr);

avg_cycle_time = mean([cycle_time_r,cycle_time_l]);

var_cycle_time = var([cycle_time_r,cycle_time_l]);

% swing - stance phase segmentation

[stance_starts_indices_r,swing_starts_indices_r,swing_stange_seg_r] = ...
        swing_stance(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));
    
[stance_starts_indices_l,swing_starts_indices_l,swing_stange_seg_l] = ...
        swing_stance(toe_l(:,2),ankle_l(:,2),toe_l(:,3),ankle_l(:,3));

% getting stance times
stance_starts_r = stance_starts_indices_r * (1/marker_sr);
stance_starts_l = stance_starts_indices_l * (1/marker_sr);
swing_starts_r = swing_starts_indices_r * (1/marker_sr);
swing_starts_l = swing_starts_indices_l * (1/marker_sr);

[pitch_foot_angle,pitch_angular_velocity] = ...
       foot_pitch_vel(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));
   
avg_stance_proportion = mean([swing_stange_seg_r,swing_stange_seg_l]);
var_stance_proportion = var([swing_stange_seg_r,swing_stange_seg_l]);

%%
%EMG signals preparation
% extensor : EDL, VL, ST, MG 
% flexor : FHL, RF, BF, TA
sr = data.EMG_sr;
%sr= EMG.sampFq;
tmin = 0;
tmax = 242900/sr;
t = linspace(tmin,tmax,242900);
delta_time = 1200.0 ;

Flexors = [data.FHL, data.RF, data.BF, data.TA];
Extensors = [data.EDL, data.VL, data.ST, data.MG];


[r_F,c_F]=size(Flexors);
Flexors_filtered = [];

for i=1:c_F
    inter_l = Filter_EMG(Flexors(:,i),sr);
    Flexors_filtered= [Flexors_filtered, inter_l];

end


[r_E,c_E]=size(Extensors);
Extensors_filtered = [];

for i=1:c_E
    inter_l = Filter_EMG(Extensors(:,i),sr);
    Extensors_filtered= [Extensors_filtered, inter_l];

end


%% Filtered muscle EMG plots
plot (t,Flexors_filtered(:,:));
% IL et RF inutiles (activity ~ 0)
figure;
plot (t,Extensors_filtered(:,:));
% ST inutile (activity ~ 0)
