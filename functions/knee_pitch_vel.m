function [angular_velocity,angle] = knee_pitch_vel(hip,knee,ankle)
%knee_pitch_vel Computes angle and angular velocity of the knee joint
    thigh = hip - knee;
    shank = knee - ankle;

    thigh = thigh(:,2:3)';
    shank = shank(:,2:3)';
    
    angle = acos(dot(thigh,shank)./...
        (vecnorm(thigh,2,1).*vecnorm(shank,2,1)));
    
    angular_velocity = diff(angle);
end

