%% Script for the plotting of H01_TDM_2kmh.mat data

%clear workspace, load the data
addpath(genpath('..'))
clear 
load('H01_TDM_2kmh.mat') % load the dataset
s = struct; %feature structure

%% Getting the data

l_ankle = data.LANK;
r_ankle = data.RANK;
l_hip = data.LHIP;
r_hip = data.RHIP;
l_kne = data.LKNE;
r_kne = data.RKNE;
l_toe = data.LTOE;
r_toe = data.RTOE;

marker_sr = data.marker_sr;


%% Displacement



displace = [zeros(1,size(r_ankle,1)); ...
            1:1:size(r_ankle,1); ...
            zeros(1,size(r_ankle,1));]' * 100 * (2/3.6)*(1/marker_sr);
        
c_l_ankle = l_ankle + displace;
c_r_ankle = r_ankle + displace;
c_l_hip = l_hip + displace;
c_r_hip = r_hip + displace;
c_l_kne = l_kne + displace;
c_r_kne = r_kne + displace;
c_l_toe = l_toe + displace;
c_r_toe = r_toe + displace;
        
color_map = [1:1:size(r_ankle,1)]';

%% swing - stance phase segmentation

[stance_starts_indices_r,swing_starts_indices_r] = ...
        swing_stance(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));
    
[stance_starts_indices_l,swing_starts_indices_l] = ...
        swing_stance(toe_l(:,2),ankle_l(:,2),toe_l(:,3),ankle_l(:,3));

    
color_left = []
color_right = []

%% Displaying the motion capture data
clf
hold on 
pcshow(c_l_ankle,color_map)
pcshow(c_r_ankle,color_map)
pcshow(c_l_hip,color_map)
pcshow(c_r_hip,color_map)
pcshow(c_l_kne,color_map)
pcshow(c_r_kne,color_map)
pcshow(c_l_toe,color_map)
pcshow(c_r_toe,color_map)
hold off