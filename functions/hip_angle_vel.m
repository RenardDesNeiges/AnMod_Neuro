function [angular_velocity,angle] = hip_angle_vel(hip,knee)
%knee_pitch_vel Computes angle and angular velocity of the hip joint
    thigh = hip - knee;
    thigh = thigh(:,:)';
    vertical = ones(size(thigh)).*[0,0,1]';

    angle = acos(dot(thigh,vertical)./...
        (vecnorm(thigh,2,1).*vecnorm(vertical,2,1)));
    angle = angle + ((angle> pi/2)*-pi);
    angular_velocity = diff(angle);
    

end

