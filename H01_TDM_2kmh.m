%% Script for the analysis of H01_TDM_2kmh.mat

%clear workspace, load the data
addpath(genpath('..'))

clear 
name = 'H01_TDM_2kmh';
load(strcat(name,'.mat')) % load the dataset 
velocity = 2; %velocity in km/h
s = struct; %feature structure
show_plots = true; %set to true to display results
%dock the figures by default to prevent mental breakdown
set(0,'DefaultFigureWindowStyle','docked') 

%% Gait events detection
clc
close all

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical mean and variance estimates of cycle time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cycle_r,cycle_index_r,cycle_time_r] = ...
    get_cycle(ankle_r(:,2),hip_r(:,2),marker_sr);
[cycle_l,cycle_index_l,cycle_time_l] = ...
    get_cycle(ankle_l(:,2),hip_l(:,2),marker_sr);

avg_cycle_time = mean([cycle_time_r,cycle_time_l]);
var_cycle_time = var([cycle_time_r,cycle_time_l]-avg_cycle_time);

if show_plots == true
    figure
    % plotting cycle length segmentation
    normalized_y_ankle = ankle_r(:,2)-hip_r(:,2);
    normalized_y_ankle_speed = diff(normalized_y_ankle);
    
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
    hold on
	plot((0:1:1000)*(1/marker_sr),normalized_y_ankle(500:1500));
    plot(cycle_r( (cycle_r < 15) & (cycle_r > 5) ) - 5,-0.1,'or');
    xlabel("time [s]")
    ylabel('normalized ankle position (ankle-hip dist) [m]')
    t = title(strcat("cycle start detection for dataset : ", ...
        name)); % avoids interpreting "_" as latex for indice
    set(t,'Interpreter','none')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% swing - stance phase segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if show_plots == true
    figure
    % plotting stance - swing segmentation
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
	plot(times,pitch_foot_angle(501:1500)');
    hold on
	plot(stance_starts_r(( ...
        stance_starts_r<1500*(1/marker_sr)) & ...
        (stance_starts_r>500*(1/marker_sr))) - ...
        500/marker_sr,-0.1,'or');
    plot(swing_starts_r(( ...
        swing_starts_r<1500*(1/marker_sr)) & ...
        (swing_starts_r>500*(1/marker_sr))) - ...
        500/marker_sr,-1,'ob');
    xlabel("time [s]")
    ylabel('pitch angle [rad]')
    t = title(strcat(...
        "foot angle with detection of toe-off/heel-contact",...
        " events for time series : ",...
        name)); % avoids interpreting _ as latex
    set(t,'Interpreter','none')
end

avg_stance_proportion = mean([swing_stange_seg_r,swing_stange_seg_l]);
var_stance_proportion = var([swing_stange_seg_r,swing_stange_seg_l]-...
                                                    avg_stance_proportion);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[max_height_r,max_height_times_r,max_height_indices_r] = ...
    max_heights(ankle_r(:,3),350,marker_sr);
[max_height_l,max_height_times_l,max_height_indices_l] = ...
    max_heights(ankle_l(:,3),350,marker_sr);

% getting stance times
avg_step_height = mean([max_height_r;max_height_l]);
var_step_height = var([max_height_r;max_height_l]-avg_step_height);
height_disymmetry = abs((mean(max_height_r)-mean(max_height_l))/...
                (mean(max_height_r)+mean(max_height_l)));

if show_plots == true
    figure
    ankle_height = ankle_r(501:1500,3);
    [max_height,max_times,max_indices] = ...
        max_heights(ankle_height,350,marker_sr);
    % plotting stance - swing segmentation
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
	plot(times,ankle_height')
    hold on
    plot(max_indices * (1/marker_sr),max_height,'or')
    xlabel("time [s]")
    %ylabel('pitch angle [rad]')
    t = title(strcat(...
        "step (ankle) height with event detection",...
        " events for time series : ",...
        name)); % avoids interpreting _ as latex
    set(t,'Interpreter','none')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude of oscillation for knee joint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[knee_angular_velocity_l,knee_angle_l] = knee_pitch_vel(hip_l,knee_l,ankle_l);
[knee_angular_velocity_r,knee_angle_r] = knee_pitch_vel(hip_r,knee_r,ankle_r);

knee_amplitude_l = 2*sqrt(mean((abs(knee_angle_l-mean(knee_angle_l)))));
knee_amplitude_r = 2*sqrt(mean((abs(knee_angle_r-mean(knee_angle_r)))));
knee_amp_asymetry = abs((knee_amplitude_l-knee_amplitude_r)/...
    (knee_amplitude_l+knee_amplitude_r));

if show_plots == true
    figure
    [angular_velocity,angle] = knee_pitch_vel(hip_l(501:1500,:),...
        knee_l(501:1500,:),ankle_l(501:1500,:));
    % plotting stance - swing segmentation
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
	plot(times,angle)
    hold on
    xlabel("time [s]")
    ylabel('knee angle [rad]')
    t = title(strcat(...
        "knee angle over time, time series : ",...
        name)); % avoids interpreting _ as latex
    set(t,'Interpreter','none')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude of oscillation for hip joint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hip_angular_velocity_l,hip_angle_l] = hip_angle_vel(hip_l,knee_l);
[hip_angular_velocity_r,hip_angle_r] = hip_angle_vel(hip_r,knee_r);

hip_amplitude_l = 2*sqrt(mean((abs(hip_angle_l-mean(hip_angle_l)))));
hip_amplitude_r = 2*sqrt(mean((abs(hip_angle_r-mean(hip_angle_r)))));
hip_amp_asymetry = abs((hip_amplitude_l-hip_amplitude_r)/...
                                (hip_amplitude_l+hip_amplitude_r));
                            
if show_plots == true
    figure
    [~,angle] = hip_angle_vel(hip_l(501:1500,:),knee_l(501:1500,:));
    % plotting stance - swing segmentation
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
	plot(times,angle)
    hold on
    xlabel("time [s]")
    ylabel('hip angle [rad]')
    t = title(strcat(...
        "hip angle over time, time series : ",...
        name)); % avoids interpreting _ as latex
    set(t,'Interpreter','none')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ankle angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% knee_right = data.RKNE;
% ankle_right = data.RANK;
% toe_right = data.RTOE;
% 
% knee_left = data.LKNE;
% ankle_left = data.LANK;
% toe_left = data.LTOE;
% 
% % sagital plane :
% vec1_r = knee_right - ankle_right;
% vec2_r = ankle_right - toe_right;
% prod_r = dot(vec1_r(:,2,3), vec2_r(:,2,3)) / dot(vecnorm(vec1_r(:,2,3)), vecnorm(vec2_r(:,2,3)));
% 
% vec1_l = knee_left - ankle_left;
% vec2_l = ankle_left - toe_left;
% prod_l = dot(vec1_l(:,2,3), vec2_l(:,2,3)) / dot(vecnorm(vec1_l(:,2,3)), vecnorm(vec2_l(:,2,3)));
% 
% angle_r = acos(prod_r);
% angle_l = acos(prod_l);
% 
% if show_plot == true
%     figure 
%     set(gcf,'color','w');
%     times = (1/marker_sr) * (1:1:1000);
%     plot(times,ankle_r(501:1500))
%     hold on
%     plot(times, ankle_l(501:1500))
%     xlabel("time [s]")
%     ylabel("Ankle angle [rad]")
% end

%%
%EMG signals preparation
% extensor : sol, ST , VLAT,MG
% flexor : TA ,II, RF
sr = data.EMG_sr;
%sr= EMG.sampFq;
tmin = 0;
tmax = 242900/sr;
t = linspace(tmin,tmax,242900);
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

save(strcat('./features/',name,'_features.mat'),'s')
