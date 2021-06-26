%% Script to test the efficiency of FFT binning

%clear workspace, load the data
addpath(genpath('..'))
clear 
load('H01_TDM_2kmh.mat') % load the dataset 
velocity = 2; %velocity in km/h
s = struct; %feature structure
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


%% Testing DFT binning

[pitch_foot_angle,pitch_angular_velocity] = ...
    foot_pitch_vel(toe_r(:,2),ankle_r(:,2),toe_r(:,3),ankle_r(:,3));

foot_pitch_fft = fftshift(fft(pitch_foot_angle-mean(pitch_foot_angle)));
fft_real = real(foot_pitch_fft);
fft_img = imag(foot_pitch_fft);

% Displaying the fft 
subplot(2,3,1)
plot(pitch_foot_angle(1:1000))
title('signal')
subplot(2,3,2)
plot(pitch_angular_velocity(1:1000))
title('derivative')
subplot(2,3,3)
plot(real(foot_pitch_fft))
title('Re(fft)')
subplot(2,3,4)
plot(imag(foot_pitch_fft))
title('Im(fft)')

re = [];
im = [];
for i = 1:500:(size(foot_pitch_fft,1)-500)
    re = [re, sum(fft_real(i:i+500))];
    im = [im, sum(fft_img(i:i+500))];
end


subplot(2,3,5)
imagesc(im)
title("Im(fft) vector")
subplot(2,3,6)
imagesc(re)
title("Re(fft) vector")