%%

clear
% Loading the data

load('H01_TDM_2kmh.mat')


%%
% Parameters
EMG = data.LTA;
sr = data.EMG_sr;
tmin = 0;
tmax = 242900/sr;
t = linspace(tmin,tmax,242900);

%% Preparation of the EMG signal

close all

% Step 1 - Slightly filter (band pass 10-2000 Hz) --> avec notre sampling
% frequency on peut aller jusqu'à 1000
Band  = [10/(sr/2), 999/(sr/2)];
[b1, a1] = butter(2, Band, 'Bandpass');   
signal_step1 = filter(b1, a1, EMG);
figure
plot(t', signal_step1)
xlim([tmin,tmax])


% Step 3 - High pass filter (30 Hz)
High = 10/(sr/2);
[b2,a2]=butter(2,High,'high'); 
signal_step2 = filter(b2,a2,signal_step1);
figure
plot(t',signal_step2)
xlim([tmin,tmax])



% Step 2 - Signal rectification
signal_step3 = abs(signal_step2);
figure
plot(t', signal_step3)
xlim([tmin,tmax])


% Step 4 - Band stop filter (around 50 Hz)
BandStop = [45/(sr/2), 55/(sr/2)];
[b3,a3]=butter(2,BandStop,'stop');
signal_step4 = filter(b3,a3,signal_step3);
figure
plot(t',signal_step4)
xlim([tmin,tmax])

% Step 5 - Low pass filter (10 Hz)
Low = 10/(sr/2);
[b4,a4]=butter(2,Low,'low'); 
signal_step5 = filter(b4,a4,signal_step4);
figure
plot(t',signal_step5)
xlim([tmin,tmax])

%% Gait events detection

%%
% Plotting the markers

% Parameters
marker_sr = data.marker_sr;
% Right:
hip = data.RHIP;
knee = data.RKNE;
ankle = data.RANK;

close all
%plot(hip(:,2))
normalized_y_ankle = ankle(:,2)-hip(:,2);
normalized_y_ankle_speed = diff(normalized_y_ankle);

Low = 7/(marker_sr/2);
[b1, a1] = butter(2,Low,'low'); 
normalized_y_ankle_speed_flt = filter(b1, a1, normalized_y_ankle_speed);




zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

zeros_ankle = (zci(normalized_y_ankle));
zeros_ankle = zeros_ankle(1:end-1);

derivative_when_zero = normalized_y_ankle_speed_flt(zeros_ankle);
cycle = (1/marker_sr) * zeros_ankle( derivative_when_zero < 0);

plot(normalized_y_ankle)
hold on
plot(normalized_y_ankle_speed_flt*15.0)

plot(cycle,0,'o')

%knee_low = filter(b2,a2,hip);
%plot(knee(:,2))
%figure
%plot(knee_low(:,2))


