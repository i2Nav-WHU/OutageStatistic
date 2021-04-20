function yaw_smooth = smooth_angle(yaw)
% this function is used for yaw interpolation
for i = 1:length(yaw)-1
    while yaw(i+1)-yaw(i) > 180
        yaw(i+1) = yaw(i+1)-360;
    end
    
    while yaw(i+1)-yaw(i) < -180
        yaw(i+1) = yaw(i+1)+360;
    end
end
yaw_smooth = yaw;