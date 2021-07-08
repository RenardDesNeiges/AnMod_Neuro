%% Script for the feature extraction from non-human primate time series

%clear workspace, load the data
addpath(genpath('..'))

clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% | Dataset name .                    | Condition                         | 
% | --------------------------------- | --------------------------------- |
% | Elektra_20190425_TM20_004         | Healthy, 2kmh walk                |
% | Elektra_20190425_TM30_002         | Healthy, 3kmh walk                |
% | Elektra_20190425_TM40_005         | Healthy, 4kmh walk                |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%change the time series name here to get features from another dataset
name = 'Elektra_20190425_TM20_004'; 

load(strcat(name,'.mat')) % load the dataset 
velocity = 2; %velocity in km/h
s = struct; %feature structure
show_plots = true; %set to true to display results
%dock the figures by default to prevent mental breakdown
set(0,'DefaultFigureWindowStyle','docked') 

% Gait events detection (from motion capture data)
clc
close all

marker_sr = Kinematic.sampFq;
emg_sr = EMG.sampFq;

% Right:
% to check that this is the right part of the dataset just run Kinematic.KINnames(20:22);
hip_r = Kinematic.data(1000:6000,20:22); 
hip_r = fillmissing(hip_r,'previous');          %get rid of NaN values
knee_r = Kinematic.data(1000:6000,23:25);
knee_r = fillmissing(knee_r,'previous');
ankle_r = Kinematic.data(1000:6000,26:28);
ankle_r = fillmissing(ankle_r,'previous');
toe_r = Kinematic.data(1000:6000,29:31);
toe_r = fillmissing(toe_r,'previous');

% Left:
hip_l = Kinematic.data(1000:6000,5:7);
hip_l = fillmissing(hip_l,'previous');
knee_l = Kinematic.data(1000:6000,8:10);
knee_l = fillmissing(knee_l,'previous');
ankle_l = Kinematic.data(1000:6000,11:13);
ankle_l = fillmissing(ankle_l,'previous');
toe_l = Kinematic.data(1000:6000,14:16);
toe_l = fillmissing(toe_l,'previous');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical mean and variance estimates of cycle time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cycle_r,cycle_index_r,cycle_time_r] = ...
    get_cycle(ankle_r(:,2),hip_r(:,2),marker_sr,7);
[cycle_l,cycle_index_l,cycle_time_l] = ...
    get_cycle(ankle_l(:,2),hip_l(:,2),marker_sr,7);

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
    if size(cycle_r( (cycle_r < 15) & (cycle_r > 5) )) > 0
        plot(cycle_r( (cycle_r < 15) & (cycle_r > 5) ) - 5,0,'or');
    end
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

[pitch_foot_angle_l,pitch_angular_velocity_l] = ...
       foot_pitch_vel(toe_l(:,2),ankle_l(:,2),toe_l(:,3),ankle_l(:,3));
[pitch_foot_angle_r,pitch_angular_velocity_r] = ...
       foot_pitch_vel(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));

   
foot_amplitude_l = 2*sqrt(mean((abs(pitch_foot_angle_l-...
                                    mean(pitch_foot_angle_l)))));
foot_amplitude_r = 2*sqrt(mean((abs(pitch_foot_angle_r-...
                                    mean(pitch_foot_angle_r)))));
foot_amp_asymetry = abs((foot_amplitude_l-foot_amplitude_r)/...
                                (foot_amplitude_l+foot_amplitude_r));
   
if show_plots == true
    figure
    % plotting stance - swing segmentation
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
	plot(times,pitch_foot_angle_r(501:1500)');
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


h_threshold = 0;

[max_height_r,max_height_times_r,max_height_indices_r] = ...
    max_heights(ankle_r(:,3),h_threshold,marker_sr);
[max_height_l,max_height_times_l,max_height_indices_l] = ...
    max_heights(ankle_l(:,3),h_threshold,marker_sr);

% getting stance times
avg_step_height = mean([max_height_r;max_height_l]);
var_step_height = var([max_height_r;max_height_l]-avg_step_height);
height_disymmetry = abs((mean(max_height_r)-mean(max_height_l))/...
                (mean(max_height_r)+mean(max_height_l)));

if show_plots == true
    figure
    ankle_height = ankle_r(501:1500,3);
    [max_height,max_times,max_indices] = ...
        max_heights(ankle_height,h_threshold,marker_sr);
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

[knee_angular_velocity_l,knee_angle_l] = ...
    knee_pitch_vel(hip_l,knee_l,ankle_l);
[knee_angular_velocity_r,knee_angle_r] = ...
    knee_pitch_vel(hip_r,knee_r,ankle_r);

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
% amplitude of oscillation for ankle joint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ankle_angular_velocity_l,ankle_angle_l] = ...
       ankle_pitch_vel(knee_l,ankle_l,toe_l);
[ankle_angular_velocity_r,ankle_angle_r] = ...
    ankle_pitch_vel(knee_r,ankle_r,toe_r);

ankle_amplitude_l = 2*sqrt(mean((abs(ankle_angle_l-mean(ankle_angle_l)))));
ankle_amplitude_r = 2*sqrt(mean((abs(ankle_angle_r-mean(ankle_angle_r)))));
ankle_amp_asymetry = abs((ankle_amplitude_l-ankle_amplitude_r)/...
                                (ankle_amplitude_l+ankle_amplitude_r));

if show_plots == true
    figure
    [~,angle] = ...
        ankle_pitch_vel(knee_l(501:1500,:),ankle_l(501:1500,:),...
        toe_l(501:1500,:));
    % plotting stance - swing segmentation
    set(gcf,'color','w');
    times = (1/marker_sr) * (1:1:1000);
	plot(times,angle)
    hold on
    xlabel("time [s]")
    ylabel('ankle angle [rad]')
    t = title(strcat(...
        "ankle angle over time, time series : ",...
        name)); % avoids interpreting _ as latex
    set(t,'Interpreter','none')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlation between different angle signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ankle_knee_angle_corr,ankle_knee_angle_asymetry] = ...
    correlate_leg_signals(ankle_angle_l,ankle_angle_r, ...
    knee_angle_l,knee_angle_r);

[ankle_foot_angle_corr,ankle_foot_angle_asymetry] = ...
    correlate_leg_signals(ankle_angle_l,ankle_angle_r, ...
    pitch_foot_angle_l',pitch_foot_angle_r');

[ankle_hip_angle_corr,ankle_hip_angle_asymetry] = ...
    correlate_leg_signals(ankle_angle_l,ankle_angle_r, ...
    hip_angle_l,hip_angle_r);

[foot_knee_angle_corr,foot_knee_angle_asymetry] = ...
    correlate_leg_signals(pitch_foot_angle_l',pitch_foot_angle_r', ...
    knee_angle_l,knee_angle_r);

[hip_knee_angle_corr,hip_knee_angle_asymetry] = ...
    correlate_leg_signals(hip_angle_l,hip_angle_r, ...
    knee_angle_l,knee_angle_r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EMG signals preparation (filtering)
% extensor : sol, ST , VLAT,MG
% flexor : TA ,II, RF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
emg_sr = EMG.sampFq;


R_Glut = EMG.data(1,:);
R_Il = EMG.data(2,:);
R_RF = EMG.data(3,:);
R_ST = EMG.data(4,:);
R_MG = EMG.data(5,:);
R_TA = EMG.data(6,:);
R_EDL = EMG.data(7,:);
R_FHL = EMG.data(8,:);


tmin = 0;
tmax = size(R_TA,1)/emg_sr;
t = linspace(tmin,tmax,242900);
delta_time = 1200.0 ;

%%
Flexors_left = [data.LTA, data.LRF];
Flexors_right = [data.RTA, data.RRF];
Extensors_left = [data.LSol, data.LST, data.LVLat, data.LMG];
Extensors_right = [data.RSol, data.RST, data.RVLat, data.RMG];

[~,c_F]=size(Flexors_left);
Flexors_left_filtered = [];
Flexors_right_filtered = [];

for i=1:c_F
    inter_l = Filter_EMG(Flexors_left(:,i),emg_sr);
    Flexors_left_filtered= [Flexors_left_filtered, inter_l];
    inter_r = Filter_EMG(Flexors_right(:,i),emg_sr);
    Flexors_right_filtered= [Flexors_right_filtered, inter_r];
end


[~,c_E]=size(Extensors_left);
Extensors_left_filtered = [];
Extensors_right_filtered = [];

for i=1:c_E
    inter_l = Filter_EMG(Extensors_left(:,i),emg_sr);
    Extensors_left_filtered= [Extensors_left_filtered, inter_l];
    inter_r = Filter_EMG(Extensors_right(:,i),emg_sr);
    
    Extensors_right_filtered= [Extensors_right_filtered, inter_r];
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% muscle activity params extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flexors

[r_F,c_F]=size(Flexors_left);
Flexors_MEAN_left = [];
Flexors_RMS_left = [];
Flexors_integral_left = [];
Flexors_MEAN_right = [];
Flexors_RMS_right = [];
Flexors_integral_right = [];

for i=1:c_F
    [mean_l,RMS_l,integ_l] = ...
        Extract_muscle_features(Flexors_left_filtered(:,i));
    Flexors_MEAN_left= [Flexors_MEAN_left, mean_l];
    Flexors_RMS_left= [Flexors_RMS_left, RMS_l];
    Flexors_integral_left= [Flexors_integral_left, integ_l];
    [mean_r,RMS_r,integ_r] = ...
        Extract_muscle_features(Flexors_right_filtered(:,i));
    Flexors_MEAN_right= [Flexors_MEAN_right, mean_r];
    Flexors_RMS_right= [Flexors_RMS_right, RMS_r];
    Flexors_integral_right= [Flexors_integral_right, integ_r];
    
end

% Extensors

[r_E,c_E]=size(Extensors_left);
Extensors_MEAN_left = [];
Extensors_RMS_left = [];
Extensors_integral_left = [];
Extensors_MEAN_right = [];
Extensors_RMS_right = [];
Extensors_integral_right = [];

for i=1:c_E
    [mean_l,RMS_l,integ_l] = ...
        Extract_muscle_features(Extensors_left_filtered(:,i));
    Extensors_MEAN_left= [Extensors_MEAN_left, mean_l];
    Extensors_RMS_left= [Extensors_RMS_left, RMS_l];
    Extensors_integral_left= [Extensors_integral_left, integ_l];
    [mean_r,RMS_r,integ_r] = ...
        Extract_muscle_features(Extensors_right_filtered(:,i));
    Extensors_MEAN_right= [Extensors_MEAN_right, mean_r];
    Extensors_RMS_right= [Extensors_RMS_right, RMS_r];
    Extensors_integral_right= [Extensors_integral_right, integ_r];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Muscle activity parameters : onset and offset 
% of extensors and flexors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extensors:

[Extensors_LSol_Onset_id,Extensors_LSol_Offset_id,...
    Extensors_LSol_Onset_t,Extensors_LSol_Offset_t,...
    Extensors_LSol_Onset_dt,Extensors_LSol_Offset_dt] = ...
    onset_offset_extraction(Extensors_left_filtered(:,1), delta_time, ...
    stance_starts_indices_l, swing_starts_indices_l,emg_sr);

avg_LSol_Onset_dt = mean(Extensors_LSol_Onset_dt);
avg_LSol_Offset_dt = mean(Extensors_LSol_Offset_dt);
var_LSol_Onset_dt = var(Extensors_LSol_Onset_dt);
var_LSol_Offset_dt = var(Extensors_LSol_Offset_dt);

[Extensors_RSol_Onset_id,Extensors_RSol_Offset_id,...
    Extensors_RSol_Onset_t,Extensors_RSol_Offset_t,...
    Extensors_RSol_Onset_dt,Extensors_RSol_Offset_dt] = ...
    onset_offset_extraction(Extensors_right_filtered(:,1), delta_time, ...
    swing_starts_indices_r, stance_starts_indices_r,emg_sr);

avg_RSol_Onset_dt = mean(Extensors_RSol_Onset_dt);
avg_RSol_Offset_dt = mean(Extensors_RSol_Offset_dt);
var_RSol_Onset_dt = var(Extensors_RSol_Onset_dt);
var_RSol_Offset_dt = var(Extensors_RSol_Offset_dt);

[Extensors_LST_Onset_id,Extensors_LST_Offset_id,...
    Extensors_LST_Onset_t,Extensors_LST_Offset_t,...
    Extensors_LST_Onset_dt,Extensors_LST_Offset_dt]= ...
    onset_offset_extraction(Extensors_left_filtered(:,2), delta_time, ...
    stance_starts_indices_l, swing_starts_indices_l,emg_sr);

avg_LST_Onset_dt = mean(Extensors_LST_Onset_dt);
avg_LST_Offset_dt = mean(Extensors_LST_Offset_dt);
var_LST_Onset_dt = var(Extensors_LST_Onset_dt);
var_LST_Offset_dt = var(Extensors_LST_Offset_dt);

[Extensors_RST_Onset_id,Extensors_RST_Offset_id,...
    Extensors_RST_Onset_t,Extensors_RST_Offset_t,...
    Extensors_RST_Onset_dt,Extensors_RST_Offset_dt]= ...
    onset_offset_extraction(Extensors_right_filtered(:,2), delta_time, ...
    swing_starts_indices_r, stance_starts_indices_r,emg_sr);

avg_RST_Onset_dt = mean(Extensors_RST_Onset_dt);
avg_RST_Offset_dt = mean(Extensors_RST_Offset_dt);
var_RST_Onset_dt = var(Extensors_RST_Onset_dt);
var_RST_Offset_dt = var(Extensors_RST_Offset_dt);

[Extensors_LVlat_Onset_id,Extensors_LVlat_Offset_id,...
    Extensors_LVlat_Onset_t,Extensors_LVlat_Offset_t,...
    Extensors_LVlat_Onset_dt,Extensors_LVlat_Offset_dt]= ...
    onset_offset_extraction(Extensors_left_filtered(:,3), delta_time, ...
    stance_starts_indices_l, swing_starts_indices_l,emg_sr);

avg_LVlat_Onset_dt = mean(Extensors_LVlat_Onset_dt);
avg_LVlat_Offset_dt = mean(Extensors_LVlat_Offset_dt);
var_LVlat_Onset_dt = var(Extensors_LVlat_Onset_dt);
var_LVlat_Offset_dt = var(Extensors_LVlat_Offset_dt);

[Extensors_RVlat_Onset_id,Extensors_RVlat_Offset_id,...
    Extensors_RVlat_Onset_t,Extensors_RVlat_Offset_t,...
    Extensors_RVlat_Onset_dt,Extensors_RVlat_Offset_dt]= ...
    onset_offset_extraction(Extensors_right_filtered(:,3), delta_time, ...
    swing_starts_indices_r, stance_starts_indices_r,emg_sr);

avg_RVlat_Onset_dt = mean(Extensors_RVlat_Onset_dt);
avg_RVlat_Offset_dt = mean(Extensors_RVlat_Offset_dt);
var_RVlat_Onset_dt = var(Extensors_RVlat_Onset_dt);
var_RVlat_Offset_dt = var(Extensors_RVlat_Offset_dt);

[Extensors_LMG_Onset_id,Extensors_LMG_Offset_id,...
    Extensors_LMG_Onset_t,Extensors_LMG_Offset_t,...
    Extensors_LMG_Onset_dt,Extensors_LMG_Offset_dt]= ...
    onset_offset_extraction(Extensors_left_filtered(:,4), delta_time, ...
    stance_starts_indices_l, swing_starts_indices_l,emg_sr);

avg_LMG_Onset_dt = mean(Extensors_LMG_Onset_dt);
avg_LMG_Offset_dt = mean(Extensors_LMG_Offset_dt);
var_LMG_Onset_dt = var(Extensors_LMG_Onset_dt);
var_LMG_Offset_dt = var(Extensors_LMG_Offset_dt);

[Extensors_RMG_Onset_id,Extensors_RMG_Offset_id,...
    Extensors_RMG_Onset_t,Extensors_RMG_Offset_t,...
    Extensors_RMG_Onset_dt,Extensors_RMG_Offset_dt]= ...
    onset_offset_extraction(Extensors_right_filtered(:,4), delta_time, ...
    swing_starts_indices_r, stance_starts_indices_r,emg_sr);

avg_RMG_Onset_dt = mean(Extensors_RMG_Onset_dt);
avg_RMG_Offset_dt = mean(Extensors_RMG_Offset_dt);
var_RMG_Onset_dt = var(Extensors_RMG_Onset_dt);
var_RMG_Offset_dt = var(Extensors_RMG_Offset_dt);

%Flexors:

[Flexors_LTA_Onset_id,Flexors_LTA_Offset_id,...
    Flexors_LTA_Onset_t,Flexors_LTA_Offset_t,...
    Flexors_LTA_Onset_dt,Flexors_LTA_Offset_dt]= ...
    onset_offset_extraction(Flexors_left_filtered(:,1), delta_time, ...
    stance_starts_indices_l, swing_starts_indices_l,emg_sr);

avg_LTA_Onset_dt = mean(Flexors_LTA_Onset_dt);
avg_LTA_Offset_dt = mean(Flexors_LTA_Offset_dt);
var_LTA_Onset_dt = var(Flexors_LTA_Onset_dt);
var_LTA_Offset_dt = var(Flexors_LTA_Offset_dt);

[Flexors_RTA_Onset_id,Flexors_RTA_Offset_id,...
    Flexors_RTA_Onset_t,Flexors_RTA_Offset_t,...
    Flexors_RTA_Onset_dt,Flexors_RTA_Offset_dt]= ...
    onset_offset_extraction(Flexors_right_filtered(:,1), delta_time, ...
    swing_starts_indices_r, stance_starts_indices_r,emg_sr);

avg_RTA_Onset_dt = mean(Flexors_RTA_Onset_dt);
avg_RTA_Offset_dt = mean(Flexors_RTA_Offset_dt);
var_RTA_Onset_dt = var(Flexors_RTA_Onset_dt);
var_RTA_Offset_dt = var(Flexors_RTA_Offset_dt);

[Flexors_LRF_Onset_id,Flexors_LRF_Offset_id,...
    Flexors_LRF_Onset_t,Flexors_LRF_Offset_t,...
    Flexors_LRF_Onset_dt,Flexors_LRF_Offset_dt]= ...
    onset_offset_extraction(Flexors_left_filtered(:,2), delta_time, ...
    stance_starts_indices_l, swing_starts_indices_l,emg_sr);

avg_LRF_Onset_dt = mean(Flexors_LRF_Onset_dt);
avg_LRF_Offset_dt = mean(Flexors_LRF_Offset_dt);
var_LRF_Onset_dt = var(Flexors_LRF_Onset_dt);
var_LRF_Offset_dt = var(Flexors_LRF_Offset_dt);

[Flexors_RRF_Onset_id,Flexors_RRF_Offset_id,...
    Flexors_RRF_Onset_t,Flexors_RRF_Offset_t,...
    Flexors_RRF_Onset_dt,Flexors_RRF_Offset_dt]= ...
    onset_offset_extraction(Flexors_right_filtered(:,2), delta_time, ...
    swing_starts_indices_r, stance_starts_indices_r,emg_sr);

avg_RRF_Onset_dt = mean(Flexors_LRF_Onset_dt);
avg_RRF_Offset_dt = mean(Flexors_LRF_Offset_dt);
var_RRF_Onset_dt = var(Flexors_LRF_Onset_dt);
var_RRF_Offset_dt = var(Flexors_LRF_Offset_dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtered muscle EMG plots with offset/onset detections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if show_plots
    figure
    set(gcf,'color','w');
    
    subplot(2,1,1)
    plot(t(10000:30000),Flexors_left_filtered(10000:30000,1));
    hold on
    %plot(Flexors_LTA_Onset_times, ...
    %    Flexors_left_filtered(Flexors_LTA_Onset_indices) ,'or');
    xlabel("time [s]")
    ylabel('ankle angle [rad]')
    tx = title(strcat(...
        "TA EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    
    subplot(2,1,2)
    plot(t(10000:30000),Flexors_left_filtered(10000:30000,2));
    xlabel("time [s]")
    ylabel('ankle angle [rad]')
    tx = title(strcat(...
        "II EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    
    figure
    
    subplot(4,1,1)
    plot (t(10000:30000),Extensors_left_filtered(10000:30000,1));
    xlabel("time [s]")
    ylabel('activity [mV]')
    tx = title(strcat(...
        "sol EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    
    subplot(4,1,2)
    plot (t(10000:30000),Extensors_left_filtered(10000:30000,2));
    xlabel("time [s]")
    ylabel('activity [mV]')
    tx = title(strcat(...
        "ST EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    
    subplot(4,1,3)
    plot (t(10000:30000),Extensors_left_filtered(10000:30000,3));
    xlabel("time [s]")
    ylabel('activity [mV]')
    tx = title(strcat(...
        "VLAT EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    
    subplot(4,1,4)
    plot (t(10000:30000),Extensors_left_filtered(10000:30000,4));
    
    xlabel("time [s]")
    ylabel('activity [mV]')
    tx = title(strcat(...
        "MG EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    set(gcf,'color','w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting structure parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.avg_cycle_time = avg_cycle_time;                  % in seconds        1
s.var_cycle_time = var_cycle_time;                  % in seconds        2
s.velocity = velocity/3.6;                          % in meter/second   3
s.avg_stance_proportion = avg_stance_proportion;    % unitless          4
s.var_stance_proportion = var_stance_proportion;    % unitless          5
s.avg_step_t = avg_step_height;                % in mm             6
s.var_step_height = var_step_height;                % in mm             7
s.height_disymmetry = height_disymmetry;            % unitless          8
s.foot_amplitude_l = foot_amplitude_l;              % in radients       9
s.foot_amplitude_r = foot_amplitude_r;              % in radients       10
s.foot_amp_asymetry = foot_amp_asymetry;            % in radients       11
s.knee_amplitude_l = knee_amplitude_l;              % in radients       12
s.knee_amplitude_r = knee_amplitude_r;              % in radients       13
s.knee_amp_asymetry = knee_amp_asymetry;            % in radients       14
s.ankle_amplitude_l = ankle_amplitude_l;            % in radients       15
s.ankle_amplitude_r = ankle_amplitude_r;            % in radients       16
s.ankle_amp_asymetry = ankle_amp_asymetry;          % in radients       17
s.ankle_knee_angle_corr = ankle_knee_angle_corr;    % unitless          18
s.ankle_knee_angle_asymetry = ...                   % unitless          19
        ankle_knee_angle_asymetry;
s.ankle_foot_angle_corr = ankle_foot_angle_corr;    % unitless          20
s.ankle_foot_angle_asymetry = ...                   % unitless          21
    ankle_foot_angle_asymetry;
s.ankle_hip_angle_corr = ankle_hip_angle_corr;      % unitless          22
s.ankle_hip_angle_asymetry = ...
    ankle_hip_angle_asymetry;                       % unitless          23
s.foot_knee_angle_corr = foot_knee_angle_corr;      % unitless          24
s.foot_knee_angle_asymetry = ...
    foot_knee_angle_asymetry;                       % unitless          25
s.hip_knee_angle_corr = hip_knee_angle_corr;        % unitless          26
s.hip_knee_angle_asymetry = hip_knee_angle_asymetry;% unitless          27

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

%Flexors:
s.Flexors_LTA_MEAN = Flexors_MEAN_left(1);
s.Flexors_LRF_MEAN = Flexors_MEAN_left(2);

s.Flexors_LTA_RMS = Flexors_RMS_left(1);
s.Flexors_LRF_RMS = Flexors_RMS_left(2);

s.Flexors_LTA_integral = Flexors_integral_left(1);
s.Flexors_LRF_integral = Flexors_integral_left(2);

s.Flexors_RTA_MEAN = Flexors_MEAN_right(1);
s.Flexors_RRF_MEAN = Flexors_MEAN_right(2);

s.Flexors_RTA_RMS = Flexors_RMS_right(1);
s.Flexors_RRF_RMS = Flexors_RMS_right(2);

s.Flexors_RTA_integral = Flexors_integral_right(1);
s.Flexors_RRF_integral = Flexors_integral_right(2);

% Onset and Offset EMG features

% Flexors:

s.avg_LSol_Onset_dt = avg_LSol_Onset_dt;
s.avg_LSol_Offset_dt = avg_LSol_Offset_dt;
s.var_LSol_Onset_dt = var_LSol_Onset_dt;
s.var_LSol_Offset_dt = var_LSol_Offset_dt;

s.avg_RSol_Onset_dt = avg_RSol_Onset_dt;
s.avg_RSol_Offset_dt = avg_RSol_Offset_dt;
s.var_RSol_Onset_dt = var_RSol_Onset_dt;
s.var_RSol_Offset_dt = var_RSol_Offset_dt;

s.avg_LST_Onset_dt = avg_LST_Onset_dt;
s.avg_LST_Offset_dt = avg_LST_Offset_dt;
s.var_LST_Onset_dt = var_LST_Onset_dt;
s.var_LST_Offset_dt = var_LST_Offset_dt;

s.avg_RST_Onset_dt = avg_RST_Onset_dt;
s.avg_RST_Offset_dt = avg_RST_Offset_dt;
s.var_RST_Onset_dt = var_RST_Onset_dt;
s.var_RST_Offset_dt = var_RST_Offset_dt;

s.avg_LVlat_Onset_dt = avg_LVlat_Onset_dt;
s.avg_LVlat_Offset_dt = avg_LVlat_Offset_dt;
s.var_LVlat_Onset_dt = var_LVlat_Onset_dt;
s.var_LVlat_Offset_dt = var_LVlat_Offset_dt;

s.avg_RVlat_Onset_dt = avg_RVlat_Onset_dt;
s.avg_RVlat_Offset_dt = avg_RVlat_Offset_dt;
s.var_RVlat_Onset_dt = var_RVlat_Onset_dt;
s.var_RVlat_Offset_dt = var_RVlat_Offset_dt;

s.avg_LMG_Onset_dt = avg_LMG_Onset_dt;
s.avg_LMG_Offset_dt = avg_LMG_Offset_dt;
s.var_LMG_Onset_dt = var_LMG_Onset_dt;
s.var_LMG_Offset_dt = var_LMG_Offset_dt;

s.avg_RMG_Onset_dt = avg_RMG_Onset_dt;
s.avg_RMG_Offset_dt = avg_RMG_Offset_dt;
s.var_RMG_Onset_dt = var_RMG_Onset_dt;
s.var_RMG_Offset_dt = var_RMG_Offset_dt;

%Flexors:

s.avg_LTA_Onset_dt = avg_LTA_Onset_dt;
s.avg_LTA_Offset_dt = avg_LTA_Offset_dt;
s.var_LTA_Onset_dt = var_LTA_Onset_dt;
s.var_LTA_Offset_dt = var_LTA_Offset_dt;


s.avg_RTA_Onset_dt = avg_RTA_Onset_dt;
s.avg_RTA_Offset_dt = avg_RTA_Offset_dt;
s.var_RTA_Onset_dt = var_RTA_Onset_dt;
s.var_RTA_Offset_dt = var_RTA_Offset_dt;

s.avg_LRF_Onset_dt = avg_LRF_Onset_dt;
s.avg_LRF_Offset_dt = avg_LRF_Offset_dt;
s.var_LRF_Onset_dt = var_LRF_Onset_dt;
s.var_LRF_Offset_dt = var_LRF_Offset_dt;


s.avg_LRF_Onset_dt = avg_LRF_Onset_dt;
s.avg_LRF_Offset_dt = avg_LRF_Offset_dt;
s.var_LRF_Onset_dt = var_LRF_Onset_dt;
s.var_LRF_Offset_dt = var_LRF_Offset_dt;

if show_plots
    figure
    feature_vec = cell2mat(struct2cell(s));
    bar(log(feature_vec))
    tx = title(strcat(...
        "log scale feature vector, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    set(gcf,'color','w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exporting the data to a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(strcat('./features/',name,'_features.mat'),'s')

