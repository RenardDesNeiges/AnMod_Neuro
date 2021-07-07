%% Script for the analysis of H01_TDM_2kmh.mat

%clear workspace, load the data
addpath(genpath('..'))
clear 
load('H01_TDM_2kmh_20_incl.mat') % load the dataset
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
% extensor : sol, ST , VLAT,MG
% flexor : TA ,II, RF
sr = data.EMG_sr;
%sr= EMG.sampFq;
tmin = 0;
tmax = 244070/sr;
t = linspace(tmin,tmax,244070);
delta_time = 1200.0 ;

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
[r_E,c_E]=size(Extensors_left);
Extensors_MEAN_left = [];
Extensors_RMS_left = [];
Extensors_integral_left = [];
Extensors_MEAN_right = [];
Extensors_RMS_right = [];
Extensors_integral_right = [];

for i=1:c_E
    [mean_l,RMS_l,integ_l] = Extract_muscle_features(Extensors_left_filtered(:,i));
    Extensors_MEAN_left= [Extensors_MEAN_left, mean_l];
    Extensors_RMS_left= [Extensors_RMS_left, RMS_l];
    Extensors_integral_left= [Extensors_integral_left, integ_l];
    [mean_r,RMS_r,integ_r] = Extract_muscle_features(Extensors_right_filtered(:,i));
    Extensors_MEAN_right= [Extensors_MEAN_right, mean_r];
    Extensors_RMS_right= [Extensors_RMS_right, RMS_r];
    Extensors_integral_right= [Extensors_integral_right, integ_r];
    
end

%% Muscle activity parameters : onset and offset of extensors and flexors
%Extensors:
[Extensors_LSol_Onset,Extensors_LSol_Offset]= onset_offset_extraction(Extensors_left_filtered(:,1), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Extensors_RSol_Onset,Extensors_RSol_Offset]= onset_offset_extraction(Extensors_right_filtered(:,1), delta_time, swing_starts_indices_r, stance_starts_indices_r);

[Extensors_LST_Onset,Extensors_LST_Offset]= onset_offset_extraction(Extensors_left_filtered(:,2), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Extensors_RST_Onset,Extensors_RST_Offset]= onset_offset_extraction(Extensors_right_filtered(:,2), delta_time, swing_starts_indices_r, stance_starts_indices_r);

[Extensors_LVlat_Onset,Extensors_LVlat_Offset]= onset_offset_extraction(Extensors_left_filtered(:,3), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Extensors_RVlat_Onset,Extensors_RVlat_Offset]= onset_offset_extraction(Extensors_right_filtered(:,3), delta_time, swing_starts_indices_r, stance_starts_indices_r);

[Extensors_LMG_Onset,Extensors_LMG_Offset]= onset_offset_extraction(Extensors_left_filtered(:,4), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Extensors_RMG_Onset,Extensors_RMG_Offset]= onset_offset_extraction(Extensors_right_filtered(:,4), delta_time, swing_starts_indices_r, stance_starts_indices_r);

%Flexors:
[Flexors_LTA_Onset,Flexors_LTA_Offset]= onset_offset_extraction(Flexors_left_filtered(:,1), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Flexors_RTA_Onset,Flexors_RTA_Offset]= onset_offset_extraction(Flexors_right_filtered(:,1), delta_time, swing_starts_indices_r, stance_starts_indices_r);

[Flexors_LIl_Onset,Flexors_LIl_Offset]= onset_offset_extraction(Flexors_left_filtered(:,2), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Flexors_RIl_Onset,Flexors_RIl_Offset]= onset_offset_extraction(Flexors_right_filtered(:,2), delta_time, swing_starts_indices_r, stance_starts_indices_r);

[Flexors_LRF_Onset,Flexors_LRF_Offset]= onset_offset_extraction(Flexors_left_filtered(:,3), delta_time, stance_starts_indices_l, swing_starts_indices_l);
[Flexors_RRF_Onset,Flexors_RRF_Offset]= onset_offset_extraction(Flexors_right_filtered(:,3), delta_time, swing_starts_indices_r, stance_starts_indices_r);

%% setting structure parameters 
s.avg_cycle_time = avg_cycle_time;                  % in seconds        1
s.var_cycle_time = var_cycle_time;                  % in seconds        2
s.velocity = velocity/3.6;                          % in meter/second   3
s.avg_stance_proportion = avg_stance_proportion;    % unitless          4
s.var_stance_proportion = var_stance_proportion;    % unitless          5
s.avg_step_height = avg_step_height;                % in mm             6
s.var_step_height = var_step_height;                % in mm             7
s.height_disymmetry = height_disymmetry;            % unitless          8
s.knee_amplitude_l = knee_amplitude_l;              % in radients       9
s.knee_amplitude_r = knee_amplitude_r;              % in radients       10
s.knee_amp_asymetry = knee_amp_asymetry;            % in radients       11

% EMG features
%Extensors:
s.Extensors_LSol_MEAN = Extensors_MEAN_left(1);
s.Extensors_LST_MEAN = Extensors_MEAN_left(2);
s.Extensors_LVlat_MEAN = Extensors_MEAN_left(3);
s.Extensors_LMG_MEAN = Extensors_MEAN_left(4);

s.Extensors_LSol_RMS = Extensors_RMS_left(1);
s.Extensors_LST_RMS = Extensors_RMS_left(2);
s.Extensors_LVlat_RMS = Extensors_RMS_left(3);
s.Extensors_LMG_RMS = Extensors_RMS_left(4);

s.Extensors_LSol_integral = Extensors_integral_left(1);
s.Extensors_LST_integral = Extensors_integral_left(2);
s.Extensors_LVlat_integral = Extensors_integral_left(3);
s.Extensors_LMG_integral = Extensors_integral_left(4);

s.Extensors_RSol_MEAN = Extensors_MEAN_right(1);
s.Extensors_RST_MEAN = Extensors_MEAN_right(2);
s.Extensors_RVlat_MEAN = Extensors_MEAN_right(3);
s.Extensors_RMG_MEAN = Extensors_MEAN_right(4);

s.Extensors_RSol_RMS = Extensors_RMS_right(1);
s.Extensors_RST_RMS = Extensors_RMS_right(2);
s.Extensors_RVlat_RMS = Extensors_RMS_right(3);
s.Extensors_RMG_RMS = Extensors_RMS_right(4);

s.Extensors_RSol_integral = Extensors_integral_right(1);
s.Extensors_RST_integral = Extensors_integral_right(2);
s.Extensors_RVlat_integral = Extensors_integral_right(3);
s.Extensors_RMG_integral = Extensors_integral_right(4);

s.Extensors_LSol_Onset = Extensors_LSol_Onset;
s.Extensors_LST_Onset = Extensors_LST_Onset;
s.Extensors_LVlat_Onset = Extensors_LVlat_Onset;
s.Extensors_LMG_Onset = Extensors_LMG_Onset;

s.Extensors_RSol_Onset = Extensors_RSol_Onset;
s.Extensors_RST_Onset = Extensors_RST_Onset;
s.Extensors_RVlat_Onset = Extensors_RVlat_Onset;
s.Extensors_RMG_Onset = Extensors_RMG_Onset;

s.Extensors_LSol_Offset = Extensors_LSol_Offset;
s.Extensors_LST_Offset = Extensors_LST_Offset;
s.Extensors_LVlat_Offset = Extensors_LVlat_Offset;
s.Extensors_LMG_Offset = Extensors_LMG_Offset;

s.Extensors_RSol_Offset = Extensors_RSol_Offset;
s.Extensors_RST_Offset = Extensors_RST_Offset;
s.Extensors_RVlat_Offset = Extensors_RVlat_Offset;
s.Extensors_RMG_Offset = Extensors_RMG_Offset;

%Flexors:
s.Flexors_LTA_MEAN = Flexors_MEAN_left(1);
s.Flexors_LIl_MEAN = Flexors_MEAN_left(2);
s.Flexors_LRF_MEAN = Flexors_MEAN_left(3);

s.Flexors_LTA_RMS = Flexors_RMS_left(1);
s.Flexors_LIl_RMS = Flexors_RMS_left(2);
s.Flexors_LRF_RMS = Flexors_RMS_left(3);

s.Flexors_LTA_integral = Flexors_integral_left(1);
s.Flexors_LIl_integral = Flexors_integral_left(2);
s.Flexors_LRF_integral = Flexors_integral_left(3);

s.Flexors_RTA_MEAN = Flexors_MEAN_right(1);
s.Flexors_RIl_MEAN = Flexors_MEAN_right(2);
s.Flexors_RRF_MEAN = Flexors_MEAN_right(3);

s.Flexors_RTA_RMS = Flexors_RMS_right(1);
s.Flexors_RIl_RMS = Flexors_RMS_right(2);
s.Flexors_RRF_RMS = Flexors_RMS_right(3);

s.Flexors_RTA_integral = Flexors_integral_right(1);
s.Flexors_RIl_integral = Flexors_integral_right(2);
s.Flexors_RRF_integral = Flexors_integral_right(3);

s.Flexors_LTA_Onset = Flexors_LTA_Onset;
s.Flexors_LIl_Onset = Flexors_LIl_Onset;
s.Flexors_LRF_Onset = Flexors_LRF_Onset;

s.Flexors_RTA_Onset = Flexors_RTA_Onset;
s.Flexors_RIl_Onset = Flexors_RIl_Onset;
s.Flexors_RRF_Onset = Flexors_RRF_Onset;

s.Flexors_LTA_Offset = Flexors_LTA_Offset;
s.Flexors_LIl_Offset = Flexors_LIl_Offset;
s.Flexors_LRF_Offset = Flexors_LRF_Offset;

s.Flexors_RTA_Offset = Flexors_RTA_Offset;
s.Flexors_RIl_Offset = Flexors_RIl_Offset;
s.Flexors_RRF_Offset = Flexors_RRF_Offset;

%% exporting the data to a file
save('./features/H01_TDM_2kmh_20_incl_features.mat','s')
