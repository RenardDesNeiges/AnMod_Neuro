function [pitch_foot_angle,pitch_angular_velocity] = foot_pitch_vel(toe_y,ankle_y,toe_z,ankle_z)
%foot_angular_velocity Compute foot pitch angle and velocity from motion
%capture data

    foot_y = toe_y-ankle_y;
    foot_z = toe_z-ankle_z;

    pitch_foot_angle = pi-mod(atan2(foot_z,foot_y),2*pi);
    pitch_angular_velocity = diff(pitch_foot_angle);

end

