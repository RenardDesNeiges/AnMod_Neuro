function [stance_starts_indices,swing_starts_indices,swing_stance_seg] = ...
                            swing_stance(toe_y,ankle_y,toe_z,ankle_z)
%swing_stance Detects swing and stance from mocap data
%   takes :
%   - toe y coordinate
%   - toe z coordinate
%   - ankle y coordinate
%   - ankle z coordinate
%   returns :
%   - stance_starts a vector containing time of the begining of
%   stance events
%   - swing_starts a vector containing time of the begining of
%   swing events
%   - swing_stance_seg a vector containing stance segmentation 

    [pitch_foot_angle,pitch_angular_velocity] = ...
        foot_pitch_vel(toe_y,ankle_y,toe_z,ankle_z);


    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

    raw_zero_indices_pitch = zci(pitch_angular_velocity);
    pitch_when_diff_0 = pitch_foot_angle(raw_zero_indices_pitch);
    maxs = raw_zero_indices_pitch(pitch_when_diff_0 > -0.2);
    mins = raw_zero_indices_pitch(pitch_when_diff_0 < -0.4);

    stance_starts_indices = [];
    swing_starts_indices = [];

    prev_indice = -1000;


     for i = 1:1:size(maxs,1)
        if maxs(i)-prev_indice > 100
           stance_starts_indices = [stance_starts_indices,maxs(i)];
        end
        prev_indice = maxs(i);
     end

     prev_indice = -1000;
     for i = 1:1:size(mins,1)
        if mins(i)-prev_indice > 100
           swing_starts_indices = [swing_starts_indices,mins(i)];
        end
        prev_indice = mins(i);
     end
     
     stance_labelled_indices = [stance_starts_indices; ...
         -ones(size(stance_starts_indices))];
     swing_labelled_indices = [swing_starts_indices; ...
         ones(size(swing_starts_indices))];
     transition_indices = [stance_labelled_indices,...
         swing_labelled_indices];
     [~,idx] = sort(transition_indices(1,:));
     transition_indices = transition_indices(:,idx);
     
     swing_stance_seg = zeros(size(toe_y,2));
     for i = transition_indices
         index = i(1);
         class = i(2);
         update = swing_stance_seg + class * ...
             [ones(1,index),zeros(1,size(toe_y,1)-index)];
         if max(update)-min(update) < 2
             swing_stance_seg = update;
         end
     end
     
     if max(swing_stance_seg) < 1
        swing_stance_seg = swing_stance_seg + 1;
     end
end

