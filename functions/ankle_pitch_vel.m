function [angular_velocity,angle] = ankle_pitch_vel(knee,ankle,toe)
%knee_pitch_vel Computes angle and angular velocity of the ankle joint
    shank = knee - ankle;
    foot = ankle - toe;

    shank = shank(:,2:3)';
    foot = foot(:,2:3)';

    angle = acos(dot(foot,shank)./...
        (vecnorm(foot,2,1).*vecnorm(shank,2,1)));
    angular_velocity = diff(angle);
end

