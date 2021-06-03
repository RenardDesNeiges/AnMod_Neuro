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
%hold off
foot_y = toe_r(:,2)-ankle_r(:,2);
foot_z = toe_r(:,3)-ankle_r(:,3);
pitch_foot_angle = pi-mod(atan2(foot_z,foot_y),2*pi);
pitch_angular_velocity = diff(pitch_foot_angle);
plot(pitch_foot_angle(:)*0.04)
hold on 
%plot(pitch_angular_velocity(:))



Low = 5/(marker_sr/2);
[b4,a4]=butter(4,Low,'low'); 
flt_abs = filter(b4,a4,(abs(pitch_angular_velocity(:))));
flt = filter(b4,a4,(pitch_angular_velocity(:)));
%plot(flt());

zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

raw_zero_indices_pitch = zci(pitch_angular_velocity)
pitch_when_diff_0 = pitch_foot_angle(raw_zero_indices_pitch)
zero_indices_pitch = zero_indices_pitch(raw_zero_indices_pitch > -2e-3)

plot(zero_indices_pitch,0,'o')

logic = (flt < 0.006);
%plot(logic()*0.02)
%hold on
%plot(foot_y)
%plot(foot_z)

%% setting structure parameters 
s.avg_cycle_time = avg_cycle_time;  % in seconds
s.var_cycle_time = var_cycle_time;  % in seconds
s.velocity = 2/3.6;                 % in meter/second
