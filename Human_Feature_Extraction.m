%% Script for the feature extraction from human patient time series

%clear workspace, load the data
addpath(genpath('..'))

clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% | Dataset name         | Species | Condition                            |
% | -------------------- | ------- | -------------------------------------|
% | H01_TDM_2kmh         | Human   | Healthy, 2kmh walk                   |
% | H01_TDM_35kmh        | Human   | Healthy, 3.5kmh walk                 |
% | H01_TDM_2kmh_20_incl | Human   | Healthy, 2kmh walk on a slope (20°)  |
% | DM002_TDM_08_2kmh    | Human   | SCI, EES treatement, 2kmh walk       |
% | DM002_TDM_08_1kmh    | Human   | SCI, EES treatement, 1khm walk       |
% | DM002_TDM_1kmh_NoEES | Human   | SCI, no EES treatement, 1kmh walk    |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%change the time series name here to get features from another dataset
name = 'H01_TDM_2kmh'; 

load(strcat(name,'.mat')) % load the dataset 
velocity = 2; %velocity in km/h
s = struct; %feature structure
show_plots = true; %set to true to display results
%dock the figures by default to prevent mental breakdown
set(0,'DefaultFigureWindowStyle','docked') 

% Gait events detection (from motion capture data)
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

knee_angle_l(isinf(knee_angle_l)|isnan(knee_angle_l))=0;
knee_angle_r(isinf(knee_angle_r)|isnan(knee_angle_r))=0;

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

hip_angle_l(isinf(hip_angle_l)|isnan(hip_angle_l))=0;
hip_angle_r(isinf(hip_angle_r)|isnan(hip_angle_r))=0;

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

ankle_angle_l(isinf(ankle_angle_l)|isnan(ankle_angle_l))=0;
ankle_angle_r(isinf(ankle_angle_r)|isnan(ankle_angle_r))=0;

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
    knee_angle_l,knee_angle_r)

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

emg_sr = data.EMG_sr;

tmin = 0;
tmax = size(data.RTA,1)/emg_sr;
t = linspace(tmin,tmax,242900);
delta_time = 1200.0 ;

Flexors_left = [data.LTA];
Extensors_left = [data.LMG];

[~,c_F]=size(Flexors_left);
Flexors_left_filtered = [];
Flexors_right_filtered = [];

for i=1:c_F
    inter_l = Filter_EMG(Flexors_left(:,i),emg_sr);
    Flexors_left_filtered= [Flexors_left_filtered, inter_l];
end


[~,c_E]=size(Extensors_left);
Extensors_left_filtered = [];


for i=1:c_E
    inter_l = Filter_EMG(Extensors_left(:,i),emg_sr);
    Extensors_left_filtered= [Extensors_left_filtered, inter_l];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% muscle activity params extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flexors

[r_F,c_F]=size(Flexors_left);
Flexors_MEAN_left = [];
Flexors_RMS_left = [];
Flexors_integral_left = [];

for i=1:c_F
    [mean_l,RMS_l,integ_l] = ...
        Extract_muscle_features(Flexors_left_filtered(:,i));
    Flexors_MEAN_left= [Flexors_MEAN_left, mean_l];
    Flexors_RMS_left= [Flexors_RMS_left, RMS_l];
    Flexors_integral_left= [Flexors_integral_left, integ_l];    
end

% Extensors

[r_E,c_E]=size(Extensors_left);
Extensors_MEAN_left = [];
Extensors_RMS_left = [];
Extensors_integral_left = [];

for i=1:c_E
    [mean_l,RMS_l,integ_l] = ...
        Extract_muscle_features(Extensors_left_filtered(:,i));
    Extensors_MEAN_left= [Extensors_MEAN_left, mean_l];
    Extensors_RMS_left= [Extensors_RMS_left, RMS_l];
    Extensors_integral_left= [Extensors_integral_left, integ_l];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Muscle activity parameters : onset and offset 
% of extensors and flexors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extensors:

[Extensors_LMG_Onset_id,Extensors_LMG_Offset_id,...
    Extensors_LMG_Onset_t,Extensors_LMG_Offset_t,...
    Extensors_LMG_Onset_dt,Extensors_LMG_Offset_dt]= ...
    onset_offset_extraction(Extensors_left_filtered(:,1), delta_time, ...
    swing_starts_indices_l, stance_starts_indices_l, emg_sr);

avg_LMG_Onset_dt = mean(Extensors_LMG_Onset_dt);
avg_LMG_Offset_dt = mean(Extensors_LMG_Offset_dt);
var_LMG_Onset_dt = var(Extensors_LMG_Onset_dt);
var_LMG_Offset_dt = var(Extensors_LMG_Offset_dt);

%Flexors:

[Flexors_LTA_Onset_id,Flexors_LTA_Offset_id,...
    Flexors_LTA_Onset_t,Flexors_LTA_Offset_t,...
    Flexors_LTA_Onset_dt,Flexors_LTA_Offset_dt]= ...
    onset_offset_extraction(Flexors_left_filtered(:,1), delta_time, ...
    swing_starts_indices_l, stance_starts_indices_l, emg_sr);

avg_LTA_Onset_dt = mean(Flexors_LTA_Onset_dt);
avg_LTA_Offset_dt = mean(Flexors_LTA_Offset_dt);
var_LTA_Onset_dt = var(Flexors_LTA_Onset_dt);
var_LTA_Offset_dt = var(Flexors_LTA_Offset_dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtered muscle EMG plots with offset/onset detections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if show_plots
    figure
    set(gcf,'color','w');
    
    subplot(1,1,1)
    plot(t(10000:30000),Flexors_left_filtered(10000:30000,1));
    hold on
    plot(Flexors_LTA_Onset_t((Flexors_LTA_Onset_t > t(10000)) &  (Flexors_LTA_Onset_t < t(30000))),.02,'or')
    plot(Flexors_LTA_Offset_t((Flexors_LTA_Offset_t > t(10000)) &  (Flexors_LTA_Offset_t < t(30000))),.02,'ob')
    xlabel("time [s]")
    ylabel('ankle angle [rad]')
    tx = title(strcat(...
        "TA EMG, time series : ",...
        name)); % avoids interpreting _ as latex
    set(tx,'Interpreter','none')
    
    figure
   
    subplot(1,1,1)
    plot (t(10000:30000),Extensors_left_filtered(10000:30000,1));
    hold on
    plot(Extensors_LMG_Onset_t((Extensors_LMG_Onset_t > t(10000)) &  (Extensors_LMG_Onset_t < t(30000))),.02,'or')
    plot(Extensors_LMG_Offset_t((Extensors_LMG_Offset_t > t(10000)) &  (Extensors_LMG_Offset_t< t(30000))),.02,'ob')
    
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
s.avg_step_t = avg_step_height;                     % in mm             6
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
s.Extensors_LMG_MEAN = Extensors_MEAN_left(1);
s.Extensors_LMG_RMS = Extensors_RMS_left(1);
s.Extensors_LMG_integral = Extensors_integral_left(1);

%Flexors:
s.Flexors_LTA_MEAN = Flexors_MEAN_left(1);
s.Flexors_LTA_RMS = Flexors_RMS_left(1);
s.Flexors_LTA_integral = Flexors_integral_left(1);

% Onset and Offset EMG features

% Flexors:

s.avg_LMG_Onset_dt = avg_LMG_Onset_dt;
s.avg_LMG_Offset_dt = avg_LMG_Offset_dt;
s.var_LMG_Onset_dt = var_LMG_Onset_dt;
s.var_LMG_Offset_dt = var_LMG_Offset_dt;

%Flexors:

s.avg_LTA_Onset_dt = avg_LTA_Onset_dt;
s.avg_LTA_Offset_dt = avg_LTA_Offset_dt;
s.var_LTA_Onset_dt = var_LTA_Onset_dt;
s.var_LTA_Offset_dt = var_LTA_Offset_dt;

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

