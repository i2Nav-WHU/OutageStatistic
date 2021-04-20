function yaw_wrap = wrap_angle(yaw)
% this function is used to set angle in range [-180, 180] deg
for i = 1:length(yaw)
    while yaw(i) > 180
        yaw(i) = yaw(i)-360;
    end
    
    while yaw(i) < -180
        yaw(i) = yaw(i)+360;
    end
end
yaw_wrap = yaw;